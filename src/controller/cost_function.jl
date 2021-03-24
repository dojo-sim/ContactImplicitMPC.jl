mutable struct CostFunction{Q,U,C,B}
    q::Vector{Q}
    u::Vector{U}
    γ::Vector{C}
    b::Vector{B}
end

function CostFunction(H::Int, dim::Dimensions;
    q = [Diagonal(zeros(SizedVector{dim.q})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{dim.u})) for t = 1:H],
    γ = [Diagonal(zeros(SizedVector{dim.c})) for t = 1:H],
    b = [Diagonal(zeros(SizedVector{dim.b})) for t = 1:H])
    return CostFunction(q, u, γ, b)
end
