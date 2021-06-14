"""
    LU solver
"""
mutable struct LUSolver14{T} <: LinearSolver
    lu::LU{T, <:AbstractMatrix{T}}
    A::Array{T, 2}
end

function lu_solver(A)
    LUSolver14(
        lu!(zero(Array(A)), check = false),
        zero(Array(A)))
end

function factorize!(solver::LUSolver14{T}, A::SparseMatrixCSC{T,Int}) where T
    solver.A .= A # sparse -> dense
    solver.lu = lu!(solver.A, check = false)
end

function factorize!(solver::LUSolver14{T}, A::AbstractMatrix{T}) where T
    solver.A .= A # sparse -> dense
    solver.lu = lu!(solver.A, check = false)
end

function linear_solve!(solver::LUSolver14{T}, x::Vector{T}, A::SparseMatrixCSC{T,Int}, b::Vector{T}) where T
    factorize!(solver, A)
    ldiv!(x, solver.lu, b) # solve
end

function linear_solve!(solver::LUSolver14{T}, x::Matrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    factorize!(solver, A)
    ldiv!(x, solver.lu, B) # solve
end
