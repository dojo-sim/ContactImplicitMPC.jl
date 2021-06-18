abstract type Objective end

mutable struct TrackingObjective{Q,U,C,B} <: Objective
    q::Vector{Q}
    u::Vector{U}
    γ::Vector{C}
    b::Vector{B}
end

function TrackingObjective(model, env, H::Int;
    q = [Diagonal(zeros(SizedVector{model.dim.q})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{model.dim.u})) for t = 1:H],
    γ = [Diagonal(zeros(SizedVector{model.dim.c})) for t = 1:H],
    b = [Diagonal(zeros(SizedVector{model.dim.c * friction_dim(env)})) for t = 1:H])
    return TrackingObjective(q, u, γ, b)
end

mutable struct TrackingVelocityObjective{Q,V,U,C,B} <: Objective
    q::Vector{Q}
    v::Vector{V}
    u::Vector{U}
    γ::Vector{C}
    b::Vector{B}
end

function TrackingVelocityObjective(model, env, H::Int;
    q = [Diagonal(zeros(SizedVector{model.dim.q})) for t = 1:H],
    v = [Diagonal(zeros(SizedVector{model.dim.q})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{model.dim.u})) for t = 1:H],
    γ = [Diagonal(zeros(SizedVector{model.dim.c})) for t = 1:H],
    b = [Diagonal(zeros(SizedVector{model.dim.c * friction_dim(env)})) for t = 1:H])
    return TrackingVelocityObjective(q, v, u, γ, b)
end
