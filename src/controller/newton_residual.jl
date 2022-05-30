
# Configurations and contact forces
struct NewtonResidualConfigurationForce{T,vq2,vu1,vγ1,vb1,vd,vI,vq0,vq1} <: NewtonResidual
    r::Vector{T}                           # residual

    q2::Vector{vq2}                    # rsd objective views
    u1::Vector{vu1}                    # rsd objective views
    γ1::Vector{vγ1}                    # rsd objective views
    b1::Vector{vb1}                    # rsd objective views

    rd::Vector{vd}                         # rsd dynamics lagrange multiplier views
    rI::Vector{vI}                         # rsd dynamics -I views [q2, γ1, b1]

    q0::Vector{vq0}                        # rsd dynamics q0 views
    q1::Vector{vq1}                        # rsd dynamics q1 views
end

function NewtonResidualConfigurationForce(model::Model, env::Environment, H::Int)

    nq = model.nq # configuration
    nu = model.nu # control
    nc = model.nc # contact
    nb = nc * friction_dim(env) # linear friction
    nd = nq + nc + nb # implicit dynamics constraint
    nr = nq + nu + nc + nb# + nd # size of a one-time-step block

    off = 0
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iγ = SizedVector{nc}(off .+ (1:nc)); off += nc # index of the impact γ1
    ib = SizedVector{nb}(off .+ (1:nb)); off += nb # index of the linear friction b1
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2

    iz = vcat(iq, iγ, ib) # index of the IP solver solution [q2, γ1, b1]

    iν = SizedVector{nd}(1:nd); #off += nd # index of the dynamics lagrange multiplier ν1

    r = zeros(H * (nr + nd))

    u1  = [view(r, (t - 1) * nr .+ iu) for t = 1:H]
    γ1  = [view(r, (t - 1) * nr .+ iγ) for t = 1:H]
    b1  = [view(r, (t - 1) * nr .+ ib) for t = 1:H]
    q2  = [view(r, (t - 1) * nr .+ iq) for t = 1:H]
    rI  = [view(r, (t - 1) * nr .+ iz) for t = 1:H]

    rd  = [view(r, H * nr + (t - 1) * nd .+ iν) for t = 1:H]

    q0  = [view(r, (t - 3) * nr .+ iq) for t = 3:H]
    q1  = [view(r, (t - 2) * nr .+ iq) for t = 2:H]

    T = eltype(r)

    return NewtonResidualConfigurationForce{T, eltype.((q2, u1, γ1, b1))...,eltype.((rd, rI, q0, q1))...}(
        r, q2, u1, γ1, b1, rd, rI, q0, q1)
end

# Configurations
struct NewtonResidualConfiguration{T,vq2,vu1,vd,vI,vq0,vq1} <: NewtonResidual
    r::Vector{T}                           # residual

    q2::Vector{vq2}                    # rsd objective views
    u1::Vector{vu1}                    # rsd objective views

    rd::Vector{vd}                         # rsd dynamics lagrange multiplier views
    rI::Vector{vI}                         # rsd dynamics -I views [q2]
    q0::Vector{vq0}                        # rsd dynamics q0 views
    q1::Vector{vq1}                        # rsd dynamics q1 views
end

function NewtonResidualConfiguration(model::Model, env::Environment, H::Int)
    nq = model.nq # configuration
    nu = model.nu # control
    nc = model.nc # contact
    nd = nq    # implicit dynamics constraint
    nr = nq + nu # size of a one-time-step block

    off = 0
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iz = copy(iq) # index of the IP solver solution [q2]

    iν = SizedVector{nd}(1:nd) # index of the dynamics lagrange multiplier ν1

    r = zeros(H * (nr + nd))

    u1  = [view(r, (t - 1) * nr .+ iu) for t = 1:H]
    q2  = [view(r, (t - 1) * nr .+ iq) for t = 1:H]
    rI  = [view(r, (t - 1) * nr .+ iz) for t = 1:H]

    rd  = [view(r, H * nr + (t - 1) * nd .+ iν) for t = 1:H]

    q0  = [view(r, (t - 3) * nr .+ iq) for t = 3:H]
    q1  = [view(r, (t - 2) * nr .+ iq) for t = 2:H]

    T = eltype(r)

    return NewtonResidualConfiguration{T, eltype.((q2, u1))...,eltype.((rd, rI, q0, q1))...}(
        r, q2, u1, rd, rI, q0, q1)
end

function NewtonResidual(model::Model, env::Environment, H::Int;
    mode = :configurationforce)

    if mode == :configurationforce
        return NewtonResidualConfigurationForce(model, env, H)
    elseif mode == :configuration
        NewtonResidualConfiguration(model, env, H)
    else
        @error "mode not implemented"
    end
end


function residual!(res::NewtonResidual, core::Newton,
    ν::Vector{Nν}, im_traj::ImplicitTrajectory, traj::ContactTraj, ref_traj::ContactTraj, window::Vector{Int}) where Nν

    # unpack
    opts = core.opts
    obj = core.obj
    res.r .= 0.0

    # Objective
    gradient!(res, obj, core, traj, ref_traj)

    for (i, t) in enumerate(window[1:end-2])
        # Lagrangian
        i >= 3 && mul!(res.q2[i-2], transpose(im_traj.δq0[t]), ν[i], 1.0, 1.0)
        i >= 2 && mul!(res.q2[i-1], transpose(im_traj.δq1[t]), ν[i], 1.0, 1.0)
        mul!(res.u1[i], transpose(im_traj.δu1[t]), ν[i], 1.0, 1.0)

        # Implicit dynamics
        res.rd[i] .+= im_traj.d[t]

        # Minus Identity term #∇qk1, ∇γk, ∇bk
        res.rI[i] .-= ν[i]
    end

    return nothing
end

function update_traj!(traj_cand::ContactTraj, traj::ContactTraj,
        ν_cand::Vector, ν::Vector, Δ::NewtonResidualConfigurationForce, α::T) where T

    H = traj_cand.H

    for t = 1:H
        traj_cand.q[t+2] .= traj.q[t+2] .- α .* Δ.q2[t]
        traj_cand.u[t] .= traj.u[t] .- α .* Δ.u1[t]
        traj_cand.γ[t] .= traj.γ[t] .- α .* Δ.γ1[t]
        traj_cand.b[t] .= traj.b[t] .- α .* Δ.b1[t]

        ν_cand[t] .= ν[t] .- α .* Δ.rd[t]
    end

    update_z!(traj_cand)
    update_θ!(traj_cand)

    return nothing
end

function update_traj!(traj_cand::ContactTraj, traj::ContactTraj,
        ν_cand::Vector, ν::Vector, Δ::NewtonResidualConfiguration, α::T) where T

    H = traj_cand.H

    for t = 1:H
        traj_cand.q[t+2] .= traj.q[t+2] .- α .* Δ.q2[t]
        traj_cand.u[t] .= traj.u[t] .- α .* Δ.u1[t]

        ν_cand[t] .= ν[t] .- α .* Δ.rd[t]
    end

    update_z!(traj_cand)
    update_θ!(traj_cand)

    return nothing
end

function gradient!(res::NewtonResidualConfigurationForce{T,vq2,vu1,vγ1,vb1,vd,vI,vq0,vq1},
        obj::TrackingObjective{Q,U,C,B},
        core::Newton{T,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS},
        traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
        ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}) where {T,vq2,vu1,vγ1,vb1,vd,vI,vq0,vq1,Q,U,C,B,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS}

    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])
        delta!(core.Δγ[t], traj.γ[t], ref_traj.γ[t])
        delta!(core.Δb[t], traj.b[t], ref_traj.b[t])

        # res.q2[t] .+= obj.q[t] * core.Δq[t]
        mul!(res.q2[t], obj.q[t], core.Δq[t], 1.0, 1.0)
        # res.u1[t] .+= obj.u[t] * core.Δu[t]
        mul!(res.u1[t], obj.u[t], core.Δu[t], 1.0, 1.0)
        # res.γ1[t] .+= obj.γ[t] * core.Δγ[t]
        mul!(res.γ1[t], obj.γ[t], core.Δγ[t], 1.0, 1.0)
        # res.b1[t] .+= obj.b[t] * core.Δb[t]
        mul!(res.b1[t], obj.b[t], core.Δb[t], 1.0, 1.0)
    end
end


function gradient!(res::NewtonResidualConfiguration{T,vq2,vu1,vd,vI,vq0,vq1},
        obj::TrackingObjective{Q,U,C,B},
        core::Newton{T,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS},
        traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
        ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}) where {T,vq2,vu1,vd,vI,vq0,vq1,Q,U,C,B,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS}

    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])

        # res.q2[t] .+= obj.q[t] * core.Δq[t]
        # res.u1[t] .+= obj.u[t] * core.Δu[t]
        mul!(res.q2[t], obj.q[t], core.Δq[t], 1.0, 1.0)
        mul!(res.u1[t], obj.u[t], core.Δu[t], 1.0, 1.0)
    end
end

function gradient!(res::NewtonResidualConfigurationForce{T,vq2,vu1,vγ1,vb1,vd,vI,vq0,vq1},
        obj::TrackingVelocityObjective{Q,V,U,C,B},
        core::Newton{T,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS},
        traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
        ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}) where {T,vq2,vu1,vγ1,vb1,vd,vI,vq0,vq1,Q,V,U,C,B,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS}

    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2] + obj.q_target[t])
        # delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])
        delta!(core.Δγ[t], traj.γ[t], ref_traj.γ[t])
        delta!(core.Δb[t], traj.b[t], ref_traj.b[t])

        # res.q2[t] .+= obj.q[t] * core.Δq[t]
        mul!(res.q2[t], obj.q[t], core.Δq[t], 1.0, 1.0)
        # res.u1[t] .+= obj.u[t] * core.Δu[t]
        mul!(res.u1[t], obj.u[t], core.Δu[t], 1.0, 1.0)
        # res.γ1[t] .+= obj.γ[t] * core.Δγ[t]
        mul!(res.γ1[t], obj.γ[t], core.Δγ[t], 1.0, 1.0)
        # res.b1[t] .+= obj.b[t] * core.Δb[t]
        mul!(res.b1[t], obj.b[t], core.Δb[t], 1.0, 1.0)

        # velocity
        # res.q2[t] .+= obj.v[t] * (traj.q[t+2] - traj.q[t+1])
        mul!(res.q2[t], obj.v[t], traj.q[t+2], 1.0, 1.0)
        mul!(res.q2[t], obj.v[t], traj.q[t+1], -1.0, 1.0)

        t == 1 && continue
        # res.q2[t-1] .-= obj.v[t] * (traj.q[t+2] - traj.q[t+1])
        mul!(res.q2[t-1], obj.v[t], traj.q[t+2], -1.0, 1.0)
        mul!(res.q2[t-1], obj.v[t], traj.q[t+1], 1.0, 1.0)
    end
end

function gradient!(res::NewtonResidualConfiguration{T,vq2,vu1,vd,vI,vq0,vq1},
    obj::TrackingVelocityObjective{Q,V,U,C,B},
    core::Newton{T,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS},
    traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
    ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}) where {T,vq2,vu1,vd,vI,vq0,vq1,Q,V,U,C,B,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS}

    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2] + obj.q_target[t])
        # delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])

        mul!(res.q2[t], obj.q[t], core.Δq[t], 1.0, 1.0)
        mul!(res.u1[t], obj.u[t], core.Δu[t], 1.0, 1.0)

        # velocity
        mul!(res.q2[t], obj.v[t], traj.q[t+2], 1.0, 1.0)
        mul!(res.q2[t], obj.v[t], traj.q[t+1], -1.0, 1.0)
        mul!(res.q2[t], obj.v[t], obj.v_target[t], -1.0, 1.0)

        t == 1 && continue
        mul!(res.q2[t-1], obj.v[t], traj.q[t+2], -1.0, 1.0)
        mul!(res.q2[t-1], obj.v[t], traj.q[t+1], 1.0, 1.0)
        mul!(res.q2[t-1], obj.v[t], obj.v_target[t], 1.0, 1.0)
    end
end
