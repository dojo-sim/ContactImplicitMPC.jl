"""triu_perm
    QR factorization and solve methods
"""
abstract type GSData{n,T} end

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

function qr_solve!(gs_data::GSData, x, b)
    qr_solve!(gs_data.qs, gs_data.rs, x, b)
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

function qr_matrix_solve!(gs_data::GSData, X::Matrix{T}, B::Matrix{T}) where {T}
    qr_matrix_solve!(gs_data.qs, gs_data.rs, X, B)
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
    Modified Gram Schmidt
"""
mutable struct MGSData{n,T} <: GSData{n,T}
    n::Int
    A::SparseMatrixCSC{T,Int}
    qs::SVector{n,Vector{T}}
    rs::Vector{Array{T,1}}
    qi_expr::Vector{Any}
    r_expr::Any
    qu_expr::Any
    re_expr::Any
    q_expr::Any
end

function MGSData(A::AbstractMatrix{T}, n::Int) where {T}
    qi_expr = []
    for j = 1:n
        # Initialize qu
        @variables Af[1:nnz(A)]
        Ae = similar(A, eltype(Af))
        Ae.nzval .= Af
        a = Ae[:,j]
        qu = zeros(n) + Vector(a)
        push!(qi_expr, eval(build_function(simplify.(qu), Af)[2]))
    end

    # Compute rk,j
    @variables qu[1:n]
    @variables qk[1:n]
    r = [transpose(qk) * qu]
    r_expr = eval(build_function(simplify.(r), qu, qk)[2])

    # Compute q_unnormalized
    @variables qu[1:n]
    @variables qk[1:n]
    @variables r[1:1]
    qu_ = qu - r[1] * qk
    qu_expr = eval(build_function(simplify.(qu_), qu, qk, r)[2])

    # Compute rj,j
    @variables qu[1:n]
    re = [sqrt(transpose(qu) * qu)]
    re_expr = eval(build_function(simplify.(re), qu)[2])

    # Compute q_normalized
    @variables qu[1:n]
    @variables re[1:1]
    q = qu ./ re[1]
    q_expr = eval(build_function(simplify.(q), qu, re)[2])

    qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
    rs = [zeros(T, 1) for i=1:Int((n + 1) * n / 2)]
    return MGSData{n,T}(n, A, qs, rs, qi_expr, r_expr, qu_expr, re_expr, q_expr)
end

function factorize!(mgs_data::MGSData{n,T}, A::SparseMatrixCSC{T,Int}) where {n,T}
    # Unpack
    qs = mgs_data.qs
    rs = mgs_data.rs
    qi_expr = mgs_data.qi_expr
    r_expr = mgs_data.r_expr
    qu_expr = mgs_data.qu_expr
    re_expr = mgs_data.re_expr
    q_expr = mgs_data.q_expr

    off = 1
    for j = 1:n
        # qi
        qi_expr[j](qs[j], A.nzval)
        for k in eachindex((1:j-1))
            # rk
            r_expr(rs[off], qs[j], qs[k])
            # qu
            qu_expr(qs[j], qs[j], qs[k], rs[off])
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

"""
    Classic Gram Schmidt
"""
mutable struct CGSData{n,T} <: GSData{n,T}
    n::Int
    A::SparseMatrixCSC{T,Int}
    qs::SVector{n,Vector{T}}
    rs::Vector{Array{T,1}}
    qi_expr::Vector{Any}
    r_expr::Vector{Any}
    qu_expr::Any
    re_expr::Any
    q_expr::Any
end

function CGSData(A::AbstractMatrix{T}, n::Int) where {T}
    qi_expr = []
    r_expr  = []
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
        r = [transpose(qk) * a]
        push!(r_expr, eval(build_function(simplify.(r), qk, Af)[2]))
    end

    # Compute q_unnormalized
    @variables qu[1:n]
    @variables qk[1:n]
    @variables r[1:1]
    qu_ = qu - r[1] * qk
    qu_expr = eval(build_function(simplify.(qu_), qu, qk, r)[2])

    # Compute rj,j
    @variables qu[1:n]
    re = [sqrt(transpose(qu) * qu)]
    re_expr = eval(build_function(simplify.(re), qu)[2])

    # Compute q_normalized
    @variables qu[1:n]
    @variables re[1:1]
    q = qu ./ re[1]
    q_expr = eval(build_function(simplify.(q), qu, re)[2])

    # return qi_expr, r_expr, qu_expr, re_expr, q_expr
    qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
    rs = [zeros(T,1) for i=1:Int((n + 1) * n / 2)]
    return CGSData{n,T}(n, A, qs, rs, qi_expr, r_expr, qu_expr, re_expr, q_expr)
end

function factorize!(cgs_data::CGSData{n,T}, A::SparseMatrixCSC{T,Int}) where {n,T}
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
            qu_expr(qs[j], qs[j], qs[k], rs[off])
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


"""
    Dense Modified Gram-Schmidt
"""
mutable struct DMGSData{n,T} <: GSData{n,T}
    n::Int
    A::Matrix{T}
    qs::SVector{n,Vector{T}}
    rs::Vector{Array{T,1}}
end

function DMGSData(A::Array{T, 2}, n::Int) where {T}
    qs = SVector{n,Vector{T}}([zeros(n) for i=1:n])
    rs = [zeros(T,1) for i=1:Int((n + 1) * n / 2)]
    return DMGSData{n,T}(n, A, qs, rs)
end

function factorize!(mgs_data::DMGSData{n,T}, A::Array{T, 2}) where {n,T}
    # Unpack
    qs = mgs_data.qs
    rs = mgs_data.rs

    off = 1
    for j in eachindex((1:n))
        # qi
        @inbounds @. @views qs[j] .= A[:,j]
        for k in eachindex((1:j-1))
            # rk
            # @inbounds @views rs[off] .= transpose(qs[j])*qs[k]
            @inbounds @views rs[off] .= LinearAlgebra.BLAS.dot(qs[j], qs[k])
            # qu
            @inbounds @. @views qs[j] .-= qs[k] .* rs[off]
            off += 1
        end
        # re
        @inbounds @views rs[off] .= LinearAlgebra.BLAS.nrm2(qs[j])
        # q
        @inbounds @. @views qs[j] ./= rs[off]
        off += 1
    end
    return nothing
end


"""
    Static Dense Modified Gram-Schmidt
"""
mutable struct SDMGSData{n,T} <: GSData{n,T}
    as::Vector{SVector{n,T}}
    qs::Vector{SVector{n,T}}
    rs::Vector{T}
    xv::Vector{T}
    xs::SVector{n,T}
end

function SDMGSData(n::Int; T::DataType=Float64)
    as = Vector{SVector{n,T}}([zeros(SVector{n,T}) for i=1:n])
    qs = Vector{SVector{n,T}}([zeros(SVector{n,T}) for i=1:n])
    rs = zeros(T,Int((n + 1) * n / 2))
    xv = zeros(T,n)
    xs = zeros(SVector{n,T})
    return SDMGSData{n,T}(as, qs, rs, xv, xs)
end

function SDMGSData(A::AbstractMatrix{T}) where T
    n, m = size(A)
    @assert n == m
    gs_data = SDMGSData(n; T=T)
    factorize!(gs_data, A)
    return gs_data
end

"""
    Gram-Schmidt algorithm perform on A.
"""
function factorize!(gs_data::SDMGSData{n,T}, A::AbstractMatrix{T}) where {n,T}
    for j in eachindex(1:n)
        @inbounds @views gs_data.as[j] = A[:,j]
    end
    factorize!(gs_data)
    return nothing
end

"""
    Gram-Schmidt algorithm perform on a.
"""
function factorize!(gs_data::SDMGSData{n,T}, a::Vector{SVector{n,T}}) where {n,T}
    gs_data.as .= a
    factorize!(gs_data)
    return nothing
end

"""
    Gram-Schmidt algorithm perform on gs_data.a.
"""
function factorize!(gs_data::SDMGSData{n,T}) where {n,T}
    # Unpack
    as = gs_data.as
    qs = gs_data.qs
    rs = gs_data.rs

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
function qr_solve!(gs_data::SDMGSData{n,T}, b) where {n,T}
    qs = gs_data.qs
    rs = gs_data.rs
    xv = gs_data.xv

    for j in eachindex(1:n)
        xv[j] = transpose(qs[j]) * b
    end
    for j = n:-1:1
        for k = j+1:n
            xv[j] -= rs[triu_perm(j, k)] * xv[k]
        end
        xv[j] /= rs[triu_perm(j, j)]
    end
    gs_data.xs = xv
    return nothing
end


"""
    QR solver
"""
abstract type QRSolver <: LinearSolver end

mutable struct EmptySolver <: LinearSolver
    F::Any
end

mutable struct CGSSolver <: QRSolver
    F::CGSData
end

mutable struct MGSSolver <: QRSolver
    F::MGSData
end

mutable struct DMGSSolver <: QRSolver
    F::DMGSData
end

mutable struct SDMGSSolver <: QRSolver
    F::SDMGSData
end

function linear_solve!(solver::QRSolver, x::Vector{T}, A::SparseMatrixCSC{T,Int}, b::Vector{T}) where T
    factorize!(solver.F, A)
    qr_solve!(solver.F, x, b)
end

function linear_solve!(solver::DMGSSolver, x::Vector{T}, A::Array{T, 2}, b::Vector{T}) where T
    factorize!(solver.F, A)
    qr_solve!(solver.F, x, b)
end

function linear_solve!(solver::QRSolver, X::Matrix{T}, A::AbstractMatrix{T}, B::Matrix{T}) where T
    factorize!(solver.F, A)
    qr_matrix_solve!(solver.F, X, B)
end

function empty_solver(A::Any)
    EmptySolver(A)
end

function mgs_solver(A::SparseMatrixCSC{T,Int}) where T
    MGSSolver(MGSData(A, size(A, 1)))
end

function cgs_solver(A::SparseMatrixCSC{T,Int}) where T
    CGSSolver(CGSData(A, size(A, 1)))
end

function dmgs_solver(A::Array{T, 2}) where T
    DMGSSolver(DMGSData(A, size(A, 1)))
end

function sdmgs_solver(A::Array{T, 2}) where T
    SDMGSSolver(SDMGSData(A))
end

mgs_solver(A::Array{T, 2}) where T = mgs_solver(sparse(A))
cgs_solver(A::Array{T, 2}) where T = cgs_solver(sparse(A))
