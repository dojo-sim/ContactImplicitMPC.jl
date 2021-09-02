abstract type LinearSolver end

mutable struct EmptySolver <: LinearSolver
    F::Any
end

function empty_solver(A::Any)
    EmptySolver(A)
end
