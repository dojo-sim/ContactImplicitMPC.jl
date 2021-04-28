abstract type Objective end

mutable struct TrackingObjective{Q,U,C,B} <: Objective
    q::Vector{Q}
    u::Vector{U}
    γ::Vector{C}
    b::Vector{B}
end

function TrackingObjective(H::Int, dim::Dimensions;
    q = [Diagonal(zeros(SizedVector{dim.q})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{dim.u})) for t = 1:H],
    γ = [Diagonal(zeros(SizedVector{dim.c})) for t = 1:H],
    b = [Diagonal(zeros(SizedVector{dim.b})) for t = 1:H])
    return TrackingObjective(q, u, γ, b)
end

function gradient!(res, obj::TrackingObjective, core, traj, ref_traj)
    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])
        delta!(core.Δγ[t], traj.γ[t], ref_traj.γ[t])
        delta!(core.Δb[t], traj.b[t], ref_traj.b[t])

        res.q2[t] .+= obj.q[t] * core.Δq[t]
        res.u1[t] .+= obj.u[t] * core.Δu[t]
        res.γ1[t] .+= obj.γ[t] * core.Δγ[t]
        res.b1[t] .+= obj.b[t] * core.Δb[t]
    end
end

function hessian!(hess, obj::TrackingObjective)
    for t = 1:length(obj.u)
        # Cost function
        hess.obj_q2[t] .+= obj.q[t]
        hess.obj_u1[t] .+= obj.u[t]
        hess.obj_γ1[t] .+= obj.γ[t]
        hess.obj_b1[t] .+= obj.b[t]
    end
end

mutable struct TrackingVelocityObjective{Q,V,U,C,B} <: Objective
    q::Vector{Q}
    v::Vector{V}
    u::Vector{U}
    γ::Vector{C}
    b::Vector{B}
end

function TrackingVelocityObjective(H::Int, dim::Dimensions;
    q = [Diagonal(zeros(SizedVector{dim.q})) for t = 1:H],
    v = [Diagonal(zeros(SizedVector{dim.q})) for t = 1:H],
    u = [Diagonal(zeros(SizedVector{dim.u})) for t = 1:H],
    γ = [Diagonal(zeros(SizedVector{dim.c})) for t = 1:H],
    b = [Diagonal(zeros(SizedVector{dim.b})) for t = 1:H])
    return TrackingVelocityObjective(q, v, u, γ, b)
end

function gradient!(res, obj::TrackingVelocityObjective, core, traj, ref_traj)
    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])
        delta!(core.Δγ[t], traj.γ[t], ref_traj.γ[t])
        delta!(core.Δb[t], traj.b[t], ref_traj.b[t])

        res.q2[t] .+= obj.q[t] * core.Δq[t]
        res.u1[t] .+= obj.u[t] * core.Δu[t]
        res.γ1[t] .+= obj.γ[t] * core.Δγ[t]
        res.b1[t] .+= obj.b[t] * core.Δb[t]

        # velocity
        res.q2[t] .+= obj.v[t] * (traj.q[t+2] - traj.q[t+1])
        t == 1 && continue
        res.q2[t-1] .-= obj.v[t] * (traj.q[t+2] - traj.q[t+1])
    end
end

function hessian!(hess, obj::TrackingVelocityObjective)
    for t = 1:length(obj.u)
        # Cost function
        hess.obj_q2[t] .+= obj.q[t]
        hess.obj_u1[t] .+= obj.u[t]
        hess.obj_γ1[t] .+= obj.γ[t]
        hess.obj_b1[t] .+= obj.b[t]

        # velocity
        hess.obj_q2[t] .+= obj.v[t]
        t == 1 && continue
        hess.obj_q2[t-1] .+= obj.v[t]
        hess.obj_q1q2[t-1] .-= obj.v[t]
        hess.obj_q2q1[t-1] .-= obj.v[t]
    end
end
