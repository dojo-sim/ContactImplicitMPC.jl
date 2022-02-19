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

function no_policy(model::ContactModel)
    NoPolicy(zeros(model.dim.u))
end

function policy(p::NoPolicy, x, traj::ContactTraj, t)
    return p.u
end

"""
    open-loop policy
"""
mutable struct OpenLoop{U} <: Policy
    u::Vector{U} # nominal controls
    idx::Int
    cnt::Int
    N_sample::Int
end

function open_loop_policy(u; N_sample = 1)
    OpenLoop(u, 0, N_sample, N_sample)
end

function policy(p::OpenLoop, x, traj, t)
    # reset
    if t == 1
        p.idx = 0
        p.cnt = p.N_sample
    end

    if p.cnt == p.N_sample
        p.idx += 1
        p.cnt = 0
    end

    p.cnt += 1

    return p.u[p.idx] ./ p.N_sample
end

"""
    control saturation
"""

control_saturation(u, uL, uU) = min.(max.(uL, u), uU)
