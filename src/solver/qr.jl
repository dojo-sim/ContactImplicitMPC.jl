"""
    QR factorization and solve methods
"""
abstract type QRSolver{n,T} <: LinearSolver end

function triangularize(rs::AbstractVector, n::Int)
    R = zeros(n, n)
    off = 1
    for j = 1:n
        for k = 1:j
            R[k, j] = rs[off][1]
            off += 1
        end
    end
    return UpperTriangular(R)
end

function triu_perm(k::Int, j::Int)
    ind = Int((j - 1) * j / 2) + k
    return ind::Int
end

function qr_solve!(gs_solver::QRSolver, x, b)
    qr_solve!(gs_solver.qs, gs_solver.rs, x, b)
end

function qr_solve!(qs, rs, x, b)
    n = length(b)
    for j = 1:n
        x[j] = transpose(qs[j]) * b
    end
    for j = n:-1:1
        for k = j+1:n
            x[j] -= rs[triu_perm(j, k)][1] * x[k]
        end
        x[j] /= rs[triu_perm(j, j)][1]
    end
    return nothing
end

function qr_matrix_solve!(gs_solver::QRSolver, X::Matrix{T}, B::Matrix{T}) where {T}
    qr_matrix_solve!(gs_solver.qs, gs_solver.rs, X, B)
end

function qr_matrix_solve!(qs, rs, X::Matrix{T}, B::Matrix{T}) where {T}
    n,m = size(B)
    for j in eachindex((1:n))
        for k in eachindex((1:m))
            @inbounds @views X[j,k] = transpose(qs[j]) * B[:,k]
        end
    end
    for j = n:-1:1
        for k = j+1:n
            @inbounds @views @. X[j,:] -= rs[triu_perm(j, k)][1] * X[k,:]
        end
        @inbounds @views @. X[j,:] /= rs[triu_perm(j, j)][1]
    end
    return nothing
end

"""
    Static Dense Modified Gram-Schmidt
"""
mutable struct SDMGSSolver{n,T} <: QRSolver{n,T}
    as::Vector{SVector{n,T}}
    qs::Vector{SVector{n,T}}
    rs::Vector{T}
    xv::Vector{T}
    xs::SVector{n,T}
end

function SDMGSSolver(n::Int; T::DataType=Float64)
    as = Vector{SVector{n,T}}([zeros(SVector{n,T}) for i=1:n])
    qs = Vector{SVector{n,T}}([zeros(SVector{n,T}) for i=1:n])
    rs = zeros(T,Int((n + 1) * n / 2))
    xv = zeros(T,n)
    xs = zeros(SVector{n,T})
    return SDMGSSolver{n,T}(as, qs, rs, xv, xs)
end

function SDMGSSolver(A::AbstractMatrix{T}) where T
    n, m = size(A)
    @assert n == m
    gs_solver = SDMGSSolver(n; T=T)
    factorize!(gs_solver, A)
    return gs_solver
end

"""
    Gram-Schmidt algorithm perform on A.
"""
function factorize!(gs_solver::SDMGSSolver{n,T}, A::AbstractMatrix{T}) where {n,T}
    for j in eachindex(1:n)
        @inbounds @views gs_solver.as[j] = A[:,j]
    end
    factorize!(gs_solver)
    return nothing
end

"""
    Gram-Schmidt algorithm perform on a.
"""
function factorize!(gs_solver::SDMGSSolver{n,T}, a::Vector{SVector{n,T}}) where {n,T}
    gs_solver.as .= a
    factorize!(gs_solver)
    return nothing
end

"""
    Gram-Schmidt algorithm perform on gs_solver.a.
"""
function factorize!(gs_solver::SDMGSSolver{n,T}) where {n,T}
    # Unpack
    as = gs_solver.as
    qs = gs_solver.qs
    rs = gs_solver.rs

    off = 1
    for j in eachindex(1:n)
        # qi
        qs[j] = as[j]
        for k in eachindex(1:j-1)
            # rk
            @inbounds rs[off] = dot(qs[j], qs[k])
            # qu
            @inbounds qs[j] -= qs[k] .* rs[off]
            off += 1
        end
        # re
        rs[off] = norm(qs[j],2)
        # q
        qs[j] /= rs[off]
        off += 1
    end
    return nothing
end

"""
    QR Back substitution after Gram-Schmidt.
"""
function qr_solve!(gs_solver::SDMGSSolver{n,T}, b) where {n,T}
    qs = gs_solver.qs
    rs = gs_solver.rs
    xv = gs_solver.xv

    for j in eachindex(1:n)
        xv[j] = transpose(qs[j]) * b
    end
    for j = n:-1:1
        for k = j+1:n
            xv[j] -= rs[triu_perm(j, k)] * xv[k]
        end
        xv[j] /= rs[triu_perm(j, j)]
    end
    gs_solver.xs = xv
    return nothing
end


"""
    QR solver
"""

function linear_solve!(solver::QRSolver, x::Vector{T}, A::SparseMatrixCSC{T,Int}, b::Vector{T}) where T
    factorize!(solver, A)
    qr_solve!(solver, x, b)
end

function linear_solve!(solver::QRSolver, X::Matrix{T}, A::AbstractMatrix{T}, B::Matrix{T}) where T
    factorize!(solver, A)
    qr_matrix_solve!(solver, X, B)
end

function sdmgs_solver(A::Array{T, 2}) where T
    SDMGSSolver(A)
end
