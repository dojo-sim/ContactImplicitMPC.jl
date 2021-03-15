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

function cgs_code_gen!(A::AbstractMatrix{T}, n::Int) where {T}
    qi_expr = []
    r_expr  = []
    qu_expr = []
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

function cgs!(
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

function triu_perm(k::Int, j::Int)
    ind = Int((j-1)*j/2)+k
    return ind::Int
end

function qr_solve!(qs::SVector{n,V}, rs::Vector{V}, x::SizedVector{n,T}, b::SizedVector{n,T}) where {n,T,V<:Vector{T}}
    for j = 1:n
        x[j] = transpose(qs[j])*b
    end
    for j = n:-1:1
        for k = j+1:n
            x[j] -= rs[triu_perm(j,k)][1]*x[k]
        end
        x[j] /= rs[triu_perm(j,j)][1]
    end
    return nothing
end




T = Float64
n = 20
A = sprand(n,n,0.16)
while rank(A) < n
    A = sprand(n,n,0.16)
end
Af = A.nzval
Q = deepcopy(A)
b = rand(n)
A\b
qi_expr, r_expr, qu_expr, re_expr, q_expr = cgs_code_gen!(A,n)

qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
rs = [zeros(T,1) for i=1:Int((n+1)*n/2)]

@time cgs!(qi_expr, r_expr, qu_expr, re_expr, q_expr, Af, qs, rs)
# @code_warntype cgs!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A, qs, rs)
@benchmark cgs!(qi_expr, r_expr, qu_expr, re_expr, q_expr, Af, qs, rs)

qs
triangularize_atom(rs,n)

norm(A - reshape(vcat(qs...), (n,n))*triangularize_atom(rs,n), Inf)




A1 = deepcopy(A)
A1.nzval .+= rand(nnz(A))
A1f = A1.nzval
qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
rs = [zeros(T,1) for i=1:Int((n+1)*n/2)]
@time cgs!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A1f, qs, rs)
norm(A1 - hcat(qs...)*triangularize_atom(rs,n), Inf)



@benchmark qr_solve!(qs, rs, c, b)
norm(c - x, Inf)

function fff(M)
    for k = 1:M
        # cgs!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A2f, qs, rs)
        qr_solve!(qs, rs, c, b)
    end
    return nothing
end
@profiler fff(10000)
