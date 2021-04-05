"""
    LU solver
"""
struct LUSolver{T} <: LinearSolver
    A::Array{T, 2}
end

function lu_solver(A)
    LUSolver(zero(Array(A)))
end

function linear_solve!(solver::LUSolver{T}, x::Vector{T}, A::SparseMatrixCSC{T,Int}, b::Vector{T}) where T
    solver.A .= A # sparse -> dense
    ldiv!(x, lu!(solver.A, check = false), b) # solve
end

function linear_solve!(solver::LUSolver{T}, x::Matrix{T}, A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    solver.A .= A # sparse -> dense
    ldiv!(x, lu!(solver.A, check = false), B) # solve
end
