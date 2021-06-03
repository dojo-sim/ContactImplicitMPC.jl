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

function no_disturbances(model::ContactModel)
    NoDisturbances(zeros(model.dim.w))
end

function disturbances(d::NoDisturbances, x, t)
    return d.w
end

"""
    open-loop disturbances
"""
mutable struct OpenLoopDisturbance{W} <: Disturbances
    w::Vector{W} # nominal disturbances
    idx::Int
    cnt::Int
    N_sample::Int
end

function open_loop_disturbances(w, N_sample=Int)
    OpenLoopDisturbance(w, 0, N_sample, N_sample)
end

function disturbances(d::OpenLoopDisturbance, x, t)
    # With N_sample = 3
    # t   = 1   2   3   4   5   6   7   8   9   10
    # idx = 1   1   1  12   2   2  23   3   3   34
    # cnt = 0   1   2  30   1   2  30   1   2   30

    # reset
    if t == 1
        d.idx = 0
        d.cnt = d.N_sample
    end

    if d.cnt == d.N_sample
        d.idx += 1
        d.cnt = 0
    end

    d.cnt += 1

    return d.w[d.idx] ./ d.N_sample
end

"""
    impulse disturbances
"""
mutable struct ImpulseDisturbance{W} <: Disturbances
    w::Vector{W} # nominal disturbances
    idx::Vector{Int}
end

function impulse_disturbances(w, idx)
    ImpulseDisturbance(w, idx)
end

function disturbances(d::ImpulseDisturbance, x, t)
    for (i, w) in enumerate(d.w)
        if d.idx[i] == t
            return w
        end
    end
    return zero(d.w[1])
end

"""
    random disturbances
    Sample random disturbance forces uniformly between 0 and -w_amplitude (one-sided).
    Apply these forces along the x (and y) axes depending on the model considered.
"""
struct RandomDisturbance{W, T} <: Disturbances
    w_amp::Vector{T} # disturbance amplitude
    w::Vector{W}     # disturbances
    t::Vector{T}     # time trajectory
end

# function random_disturbances(model::ContactModel, w_amp::Vector{T}, H::Int, h::T) where {T}
#     @assert length(w_amp) == model.dim.w
#     w = [-rand(model.dim.w) * w_amp for i=1:H]
#     RandomDisturbance(w_amp, w, [(t - 1) * h for t = 1:H])
# end

function random_disturbances(model::ContactModel, w_amp::Vector{T}, H::Int, h::T) where {T}
    if length(w_amp) == model.dim.w
        w = [rand(model.dim.w) .* w_amp for i=1:H]
        # wr = rand(model.dim.w) .* w_amp
        # w = [wr for i=1:H]
    elseif length(w_amp) == 1
        w = [rand(model.dim.w) .* w_amp[1] for i=1:H]
        # wr = rand(model.dim.w) .* w_amp[1]
        # w = [wr for i=1:H]
    else
        @warn "w_amp is not of the correct size."
    end
    RandomDisturbance(w_amp, [rand(model.dim.w) .* w_amp[1] for i=1:H], [(t - 1) * h for t = 1:H])
end

function disturbances(d::RandomDisturbance, x, t)
    k = searchsortedlast(d.t, t)
    return d.w[k]
end
