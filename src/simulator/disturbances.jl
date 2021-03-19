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
