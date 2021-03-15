using Pkg
Pkg.activate(joinpath(@__DIR__, "../src/mgs/"))

using Symbolics
using LinearAlgebra
using Test
using SparseArrays
using Random
using BenchmarkTools
using Plots
using StaticArrays

function mgs!(A::AbstractMatrix{T}, R::UpperTriangular{T}) where {T}
    n = size(A)[1]
    @assert n == size(A)[2] == size(R)[1] == size(R)[2]
    for j = 1:n
        for k = 1:j-1
            R[k,j] = A[:,k]'*A[:,j]
            A[:,j] = A[:,j] - R[k,j]*A[:,k]
        end
        R[j,j] = norm(A[:,j], 2)
        A[:,j] ./= R[j,j]
    end
    return A,R
end

function mgs_flat!(A::AbstractMatrix{T}, R::AbstractVector{T}) where {T}
    for j = 1:n
        for k = 1:j-1
            R[triu_ind(k,j,n)] = A[:,k]'*A[:,j]
            A[:,j] = A[:,j] - R[triu_ind(k,j,n)]*A[:,k]
        end
        R[triu_ind(j,j,n)] = sqrt(transpose(A[:,j]) * A[:,j])
        A[:,j] ./= R[triu_ind(j,j,n)]
    end
    return A,R
end

function mgs_sparse!(A::AbstractMatrix{T}, R::AbstractVector{T}) where {T}
    Q = copy(A)
    for j = 1:n
        Q[:,j] = A[:,j]
        for k = 1:j-1
            R[triu_ind(k,j,n)] = A[:,k]'*q
            A[:,j] .= A[:,j] - R[triu_ind(k,j,n)]*A[:,k]
            A[:,j] = simplify(A[:,j])
        end
        R[triu_ind(j,j,n)] = sqrt(transpose(A[:,j]) * A[:,j])
        A[:,j] ./= R[triu_ind(j,j,n)]
        A[:,j] = simplify(A[:,j])
    end
    return A,R
end

function triu_ind(i::Int, j::Int, n::Int)
    @assert i <= j
    return Int((n+n-i+1)*i/2) -n + j
end

function triangularize(Rf::AbstractVector{T}, n::Int) where {T}
    R = zeros(n,n)
    for i = 1:n
        for j = i:n
            R[i,j] = Rf[triu_ind(i,j,n)]
        end
    end
    return R
end

function triangularize_vect(rs)
    n = length(rs)
    R = zeros(n,n)
    for j = 1:n
        R[1:j,j] .= rs[j]
    end
    return UpperTriangular(R)
end

function triangularize_atom(rs,n::Int)
    R = zeros(n,n)
    off = 1
    for j = 1:n
        for k = 1:j
            R[k,j] = rs[off][1]
            off += 1
        end
    end
    return UpperTriangular(R)
end

function test_triu(n::Int)
    R = UpperTriangular(zeros(n,n))
    for i = 1:n
        for j = i:n
            R[i,j] = triu_ind(i,j,n)
        end
    end
    return R
end


function cgs_best!(A::AbstractMatrix{T}, n::Int) where {T}
    qi_expr = []
    r_expr  = []
    qu_expr = []
    re_expr = []
    q_expr  = []
    for j = 1:n
        # Initialize qu
        @variables Af[1:nnz(A)]
        Ae = similar(A, eltype(Af))
        Ae.nzval .= Af
        a = Ae[:,j]
        qu = zeros(n) + Vector(a)
        push!(qi_expr, eval(build_function(simplify.(qu), Af)[2]))

        # Compute rk,j
        @variables Af[1:nnz(A)]
        Ae = similar(A, eltype(Af))
        Ae.nzval .= Af
        a = Ae[:,j]
        @variables qk[1:n]
        r = [transpose(qk)*a]
        push!(r_expr, eval(build_function(simplify.(r), qk, Af)[2]))

        # Compute q_unnormalized
        @variables qu[1:n]
        @variables qk[1:n]
        @variables r[1:1]
        qu_ = qu - r[1]*qk
        push!(qu_expr, eval(build_function(simplify.(qu_), qu, qk, r)[2]))
    end

    # Compute rj,j
    @variables qu[1:n]
    re = [sqrt(transpose(qu)*qu)]
    re_expr = eval(build_function(simplify.(re), qu)[2])

    # Compute q_normalized
    @variables qu[1:n]
    @variables re[1:1]
    q = qu ./ re[1]
    q_expr = eval(build_function(simplify.(q), qu, re)[2])

    return qi_expr, r_expr, qu_expr, re_expr, q_expr
end


function add_expr!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A::SparseMatrixCSC{T,Int}, j::Int, n::Int) where {T}
    # Initialize qu
    @variables Af[1:nnz(A)]
    Ae = similar(A, eltype(Af))
    Ae.nzval .= Af
    a = Ae[:,j]
    qu = zeros(n) + Vector(a)
    push!(qi_expr, eval(build_function(simplify.(qu), Af)[2]))

    for k = 1:j-1
        # Compute rk,j
        @variables Af[1:nnz(A)]
        Ae = similar(A, eltype(Af))
        Ae.nzval .= Af
        a = Ae[:,j]
        @variables qk[1:n]
        r = [transpose(qk)*a]
        push!(r_expr[j], eval(build_function(simplify.(r), qk, Af)[2]))
        # push!(r_expr[j], simplify.(r))

        # Compute q_unnormalized
        @variables qu[1:n]
        @variables qk[1:n]
        @variables r[1:1]
        qu_ = qu - r[1]*qk
        push!(qu_expr[j], eval(build_function(simplify.(qu_), qu, qk, r)[2]))
    end

    # Compute rj,j
    @variables qu[1:n]
    re = [sqrt(transpose(qu)*qu)]
    push!(re_expr, eval(build_function(simplify.(re), qu)[2]))

    # Compute q_normalized
    @variables qu[1:n]
    @variables re[1:1]
    q = qu ./ re[1]
    push!(q_expr, eval(build_function(simplify.(q), qu, re)[2]))
    return nothing
end

function qr_cgs2!(
    qi_expr::Vector{Any},
    r_expr::Vector{Any},
    qu_expr::Vector{Any},
    re_expr::Any,
    q_expr::Any,
    Af::AbstractVector{T},
    qs::SVector{n,Vq},
    rs::Vector{Vr},
    ) where {n, T, Vq <: AbstractVector{T}, Vr <: Vector{T}}
    off = 1
    for j in eachindex(qi_expr)
        # qi
        qi_expr[j](qs[j], Af)

        for k in eachindex((1:j-1))
            # rk
            r_expr[j](rs[off], qs[k], Af)
            # qu
            qu_expr[j](qs[j], qs[j], qs[k], rs[off])
            off += 1
        end
        # re
        re_expr(rs[off], qs[j])
        # q
        q_expr(qs[j], qs[j], rs[off])
        off += 1
    end
    return nothing
end



# @benchmark qr_fast2!($q_expr, $qu_expr, $r_expr, $re_expr, $r_gen, $qu_gen, $re_gen, $A, $qs, $rs)
# @code_warntype qr_fast!(q_expr, qu_expr, r_expr, re_expr, A, qs, rs)

T = Float64
n = 43
A = sprand(n,n,0.16)
while rank(A) < n
    A = sprand(n,n,0.16)
end
Af = A.nzval
Q = deepcopy(A)
b = rand(n)
A\b
qi_expr, r_expr, qu_expr, re_expr, q_expr = cgs_best!(A,n)

qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
rs = [zeros(T,1) for i=1:Int((n+1)*n/2)]

@time qr_cgs2!(qi_expr, r_expr, qu_expr, re_expr, q_expr, Af, qs, rs)
# @code_warntype qr_cgs2!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A, qs, rs)
@benchmark qr_cgs2!(qi_expr, r_expr, qu_expr, re_expr, q_expr, Af, qs, rs)

qs
triangularize_atom(rs,n)

norm(Matrix(A) - reshape(vcat(qs...), (n,n))*triangularize_atom(rs,n), Inf)
Matrix(A)
reshape(vcat(qs...), (n,n))*triangularize_atom(rs,n)

A1 = deepcopy(A)
A1.nzval .+= rand(nnz(A))
A1f = A1.nzval
qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
rs = [zeros(T,1) for i=1:Int((n+1)*n/2)]
@time qr_cgs2!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A1f, qs, rs)
norm(Matrix(A1) - hcat(qs...)*triangularize_atom(rs,n), Inf)


A2 = deepcopy(A)
A2.nzval .+= rand(nnz(A))
A2f = A2.nzval
@benchmark A2\b # 160 μs
qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
rs = [zeros(T,1) for i=1:Int((n+1)*n/2)]
@benchmark qr_cgs2!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A2f, qs, rs)
Q = hcat(qs...)
R = triangularize_atom(rs,n)
@benchmark R\(Q'b) # 1.5μs



function fff(M)
    for k = 1:M
        qr_cgs2!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A2f, qs, rs)
    end
    return nothing
end
@profiler fff(10000)

a = 11111
a = 11111
a = 11111




# Basic test
n = 43
R = UpperTriangular(Float64.(rand(-n:n,n,n)))
A_ = rand(-n:n, n,n)
A = Float64.(A_'*A_)
b = Float64.(rand(-n:n,n))

Q = deepcopy(A)
mgs!(Q,R)
@test norm(A - Q*R, Inf) < 1e-8

Q = deepcopy(A)
Rf = zeros(Int((n+1)*n/2))
mgs_flat!(Q,Rf)
R = triangularize(Rf,n)
@test norm(A - Q*R, Inf) < 1e-8


@benchmark A\b # 160 μs
Q = deepcopy(A)
R = UpperTriangular(zeros(n,n))
@benchmark mgs!(Q,R) # 600μs
@benchmark R\(Q'b) # 1.5μs
