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
    q_expr = []
    qu_expr = []
    r_expr = []
    re_expr = []
    for j = 1:n
        # # j = 1
        # @variables a[1:n]
        # q = a
        # r = zeros(eltype(q), 1)
        # r[end] = sqrt(transpose(q)*q)
        # q = q ./ r[end]
        #
        # q_expr = build_function(q, a)[2]
        # r_expr = build_function(r, a)[2]
        #
        # # j = 2
        # @variables a[1:n]
        # @variables q1[1:n]
        # q = a
        # r = zeros(eltype(a),2)
        # r[1] = transpose(q1)*q
        # q = q - r[1]*q1
        # r[end] = sqrt(transpose(q)*q)
        # q = q ./ r[end]

        # j = j
        add_expr!(q_expr, qu_expr, r_expr, re_expr, A, j, n)
    end
    return SVector{n,Any}(q_expr), SVector{n,Any}(qu_expr), SVector{n,Any}(r_expr), SVector{n,Any}(re_expr)
end

function add_expr!(q_expr, qu_expr, r_expr, re_expr, A::SparseMatrixCSC{T,Int}, j::Int, n::Int) where {T}
    # Compute r1,j ... rj-1,j
    @variables Af[1:nnz(A)]
    Ae = similar(A, eltype(Af))
    Ae.nzval .= Af
    a = Ae[:,j]
    @variables qp[1:n*(j-1)]
    qs = [qp[(i-1)*n .+ (1:n)] for i=1:j-1]
    q = a
    r = zeros(eltype(a), j-1)
    for i = 1:j-1
        r[i] = transpose(qs[i])*q
    end
    push!(r_expr, eval(build_function(simplify.(r), Af, qs...)[2]))

    # Compute rjj and q_unormalized
    @variables Af[1:nnz(A)]
    Ae = similar(A, eltype(Af))
    Ae.nzval .= Af
    a = Ae[:,j]
    @variables qp[1:n*(j-1)]
    @variables r[1:j-1]
    qs = [qp[(i-1)*n .+ (1:n)] for i=1:j-1]
    qu = zeros(n) + Vector(a)
    for i = 1:j-1
        qu = qu - r[i]*qs[i]
    end
    re = [sqrt(transpose(qu)*qu)]
    push!(qu_expr, eval(build_function(simplify.(qu), Af, qs..., r)[2]))
    push!(re_expr, eval(build_function(simplify.(re), Af, qs..., r)[2]))

    # Compute q_normalized
    @variables qu[1:n]
    @variables re[1:1]
    q = qu ./ re[1]
    push!(q_expr, eval(build_function(simplify.(q), qu, re)[2]))
    return nothing
end

function qr_fast2!(q_expr::SVector{n,Any}, qu_expr::SVector{n,Any}, r_expr::SVector{n,Any},
    re_expr::SVector{n,Any}, r_gen::SVector{n,Any}, qu_gen::SVector{n,Any}, re_gen::SVector{n,Any},
    A::AbstractMatrix{T}, qs::SVector{n,Vq},
    rs::SVector{n,Vr}) where {n, Vq <: AbstractVector{T}, Vr <: AbstractVector{T}}
    for j in eachindex(q_expr)
        # r
        fr! = r_expr[j]
        # fr!(rs[j], A.nzval, qs[1:j-1]...)
        eval(r_gen[j])
        # rs[j] .= rsj

        # qu, re
        # fqu! = qu_expr[j]
        # fqu!(qs[j], A.nzval, qs[1:j-1]..., rs[j][1:end-1])
        eval(qu_gen[j])

        re = rs[j][end:end]
        # fre! = re_expr[j]
        # fre!(re, A.nzval, qs[1:j-1]..., rs[j][1:end-1])
        # eval(re_gen[j])
        rs[j][end:end] .= re

        # q
        fq! = q_expr[j]
        fq!(qs[j], qs[j], rs[j][end:end])
    end
    return qs, rs
end

function get_r_gen(n::Int)
    gen = []
    for j = 1:n
        e = "r_expr[$j](rs[$j], A.nzval"
        for i = 1:j-1
            e *= ", qs[$(i)]"
        end
        e *= ")"
        # push!(gen, e)
        push!(gen, Meta.parse(e))
    end
    return SVector{n,Any}(gen)
end

function get_qu_gen(n::Int)
    gen = []
    for j = 1:n
        # fqu!(qs[j], A.nzval, qs[1:j-1]..., rs[j][1:end-1])

        e = "qu_expr[$j](qs[$j], A.nzval"
        for i = 1:j-1
            e *= ", qs[$(i)]"
        end
        e *= ", rs[$j][1:end-1])"
        # push!(gen, e)
        push!(gen, Meta.parse(e))
    end
    return SVector{n,Any}(gen)
end

function get_re_gen(n::Int)
    gen = []
    for j = 1:n
        # fre!(re, A.nzval, qs[1:j-1]..., rs[j][1:end-1])

        e = "re_expr[$j](re, A.nzval"
        for i = 1:j-1
            e *= ", qs[$(i)]"
        end
        e *= ", rs[$j][1:end-1])"
        # push!(gen, e)
        push!(gen, Meta.parse(e))
    end
    return SVector{n,Any}(gen)
end

r_gen = get_r_gen(n)
qu_gen = get_qu_gen(n)
re_gen = get_re_gen(n)

@benchmark qr_fast2!($q_expr, $qu_expr, $r_expr, $re_expr, $r_gen, $qu_gen, $re_gen, $A, $qs, $rs)
# @code_warntype qr_fast!(q_expr, qu_expr, r_expr, re_expr, A, qs, rs)

T = Float64
n = 8
A = sprand(n,n,0.16)
while rank(A) < n
    A = sprand(n,n,0.16)
end
Q = deepcopy(A)
b = rand(n)
A\b
q_expr, qu_expr, r_expr, re_expr = cgs_best!(A,n)

qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
rs = SVector{n,Vector{T}}([zeros(i) for i=1:n])
@time qr_fast2!(q_expr, qu_expr, r_expr, re_expr, r_gen, qu_gen, re_gen, A, qs, rs)
@benchmark qr_fast2!($q_expr, $qu_expr, $r_expr, $re_expr, $r_gen, $qu_gen, $re_gen, $A, $qs, $rs)

qs
rs
triangularize_vect(rs)


norm(Matrix(A) - reshape(vcat(qs...), (n,n))*triangularize_vect(rs), Inf)
Matrix(A)
reshape(vcat(qs...), (n,n))*triangularize_vect(rs)

SVector{n,Any}(q_expr)



qq = ones(SizedVector{n})
rr = 10*ones(SizedVector{1})
@benchmark q_expr[3]($qq, $qq, $rr)
qq

plot([nnz(Q[:,i]) for i=1:n])
plot!([nnz(A[:,i]) for i=1:n])
sparse(R)





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
