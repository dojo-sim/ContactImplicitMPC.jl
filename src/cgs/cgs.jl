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

mutable struct CGSData13{n,T}
    n::Int
    A::SparseMatrixCSC{T,Int}
    qs::SVector{n,Vector{T}}
    rs::Vector{Array{T,1}}
    qi_expr::Vector{Any}
    r_expr::Vector{Any}
    qu_expr::Vector{Any}
    re_expr::Any
    q_expr::Any
end

function triangularize(rs,n::Int)
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

function triu_perm(k::Int, j::Int)
    ind = Int((j-1)*j/2)+k
    return ind::Int
end

function CGSData13!(A::AbstractMatrix{T}, n::Int) where {T}
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

    # return qi_expr, r_expr, qu_expr, re_expr, q_expr
    qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
    rs = [zeros(T,1) for i=1:Int((n+1)*n/2)]
    return CGSData13{n,T}(n, A, qs, rs, qi_expr, r_expr, qu_expr, re_expr, q_expr)
end

function cgs!(cgs_data::CGSData13{n,T}, A::SparseMatrixCSC{T,Int}) where {n,T}
    # Unpack
    qs = cgs_data.qs
    rs = cgs_data.rs
    qi_expr = cgs_data.qi_expr
    r_expr = cgs_data.r_expr
    qu_expr = cgs_data.qu_expr
    re_expr = cgs_data.re_expr
    q_expr = cgs_data.q_expr

    off = 1
    for j = 1:n
        # qi
        qi_expr[j](qs[j], A.nzval)
        for k in eachindex((1:j-1))
            # rk
            r_expr[j](rs[off], qs[k], A.nzval)
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

function qr_solve!(cgs_data::CGSData13{n,T}, x::SizedVector{n,T}, b::SizedVector{n,T}) where {n,T}
    qr_solve!(cgs_data.qs, cgs_data.rs, x, b)
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


# # Generate a sparse random non-singular matrix
# T = Float64
# n = 20
# Random.seed!(100)
# A = sprand(n,n,0.16)
# while rank(A) < n
#     A = sprand(n,n,0.16)
# end
#
# # Generate the code to handle matrices with A's sparsity pattern
# cgs_data = CGSData13!(A,n)
#
# # Compute the Q and R matrices
# @time cgs!(cgs_data, A)
# @benchmark cgs!(cgs_data, A)
# @test norm(A - hcat(cgs_data.qs...)*triangularize(cgs_data.rs,n), Inf) < 1e-10
#
# # Test it on a different matrix with the same sparsity pattern
# A1 = deepcopy(A)
# A1.nzval .+= rand(nnz(A))
# @time cgs!(cgs_data, A1)
# @test norm(A1 - hcat(cgs_data.qs...)*triangularize(cgs_data.rs,n), Inf) < 1e-10
#
# # Execute the simple backslash solve
# b = rand(SizedVector{n,T})
# c = zeros(SizedVector{n,T})
# @benchmark qr_solve!(cgs_data, c, b)
# x = A1\b
# @test norm(c - x, Inf) < 1e-10


# function profile(M)
#     for k = 1:M
#         # cgs!(qi_expr, r_expr, qu_expr, re_expr, q_expr, A2f, qs, rs)
#         qr_solve!(qs, rs, c, b)
#     end
#     return nothing
# end
# @profiler profile(10000)
