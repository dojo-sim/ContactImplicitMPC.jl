abstract type Policy end

function policy(p::Policy, traj::ContactTraj, t::T) where T
    @warn "policy not defined"
    return nothing
end

"""
    no policy
"""
struct NoPolicy{T} <: Policy
    u::Vector{T}
end

function no_policy(model::ContactDynamicsModel)
    NoPolicy(zeros(model.dim.u))
end

function policy(p::NoPolicy, x, traj::ContactTraj, t)
    return p.u
end

"""
    open-loop policy
"""
struct OpenLoop{U, T} <: Policy
    u::Vector{U} # nominal controls
    t::Vector{T} # time trajectory
    N_sample::Int
end

function open_loop_policy(u, h; N_sample = 1)
    OpenLoop(u, [(t - 1) * h for t = 1:length(u)], N_sample)
end

function policy(p::OpenLoop, x, traj, t)
    k = searchsortedlast(p.t, t)
    return p.u[k] ./ p.N_sample
end
