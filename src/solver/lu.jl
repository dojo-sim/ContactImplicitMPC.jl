"""
    LU solver
"""
function lu_solver(A::SparseMatrixCSC{T,Int}) where T
    A = Array(A) 
    lu_solver(A)
end

function linear_solve!(solver::LUSolver{T}, x::Vector{T}, A::SparseMatrixCSC{T,Int},
        b::Vector{T}; reg::T = 0.0, fact::Bool = true) where T
    linear_solve!(solver, x, Array(A), b, reg=reg, fact=fact)
end