abstract type Disturbances end

function disturbances(d::Disturbances, x, t)
    @warn "disturbances not defined"
    return nothing
end

"""
    no disturbances
"""
struct NoDisturbances{T} <: Disturbances
    w::Vector{T}
end

function no_disturbances(model::ContactDynamicsModel)
    NoDisturbances(zeros(model.dim.w))
end

function disturbances(d::NoDisturbances, x, t)
    return d.w
end

"""
    open-loop disturbances
"""
struct OpenLoopDisturbance{W, T} <: Disturbances
    w::Vector{W} # nominal disturbances
    t::Vector{T} # time trajectory
end

function open_loop_disturbances(w, h)
    OpenLoopDisturbance(w, [(t - 1) * h for t = 1:length(w)])
end

function disturbances(d::OpenLoopDisturbance, x, t)
    k = searchsortedlast(d.t, t)
    return p.w[k]
end
