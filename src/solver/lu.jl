"""
    LU solver
"""
mutable struct LUSolver{T} <: LinearSolver
    lu::LU{T, <:AbstractMatrix{T}}
    A::Array{T, 2}
end

function lu_solver(A)
    LUSolver(
        lu!(zero(Array(A)), check = false),
        zero(Array(A)))
end

# function factorize!(solver::LUSolver{T}, A::SparseMatrixCSC{T,Int}) where T
#     solver.A .= A # sparse -> dense
#     solver.lu = lu!(solver.A, check = false)
# end

function factorize!(solver::LUSolver{T}, A::AbstractMatrix{T}) where T
    solver.A .= A # sparse -> dense
    solver.lu = lu!(solver.A, check = false)
end

function linear_solve!(solver::LUSolver{T}, x::Vector{T}, A::Matrix{T},
        b::Vector{T}; reg::T = 0.0, fact::Bool = true) where T
    #reg is useless but used for Mehrotra
    fact && factorize!(solver, A)
    ldiv!(x, solver.lu, b) # solve
end

function linear_solve!(solver::LUSolver{T}, x::Vector{T}, A::SparseMatrixCSC{T,Int},
        b::Vector{T}; reg::T = 0.0, fact::Bool = true) where T
    #reg is useless but used for Mehrotra
    fact && factorize!(solver, A)
    ldiv!(x, solver.lu, b) # solve
end

function linear_solve!(solver::LUSolver{T}, x::Matrix{T}, A::AbstractMatrix{T},
        B::AbstractMatrix{T}; reg::T = 0.0, fact::Bool = true) where T
    #reg is useless but used for Mehrotra
    fact && factorize!(solver, A)
    ldiv!(x, solver.lu, B) # solve
end
