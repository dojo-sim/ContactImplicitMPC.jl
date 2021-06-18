
# Configurations and forces
struct NewtonJacobianConfigurationForce{T,Vq,Vu,Vγ,Vb,VI,VIT,Vq0,Vq0T,Vq1,Vq1T,Vu1,Vu1T,Vreg} <: NewtonJacobian
    R#::SparseMatrixCSC{T,Int}                 # jacobian

    obj_q2::Vector{Vq}                          # obj views
    obj_u1::Vector{Vu}                          # obj views
    obj_γ1::Vector{Vγ}                          # obj views
    obj_b1::Vector{Vb}                          # obj views

    obj_q1q2::Vector{Vq}
    obj_q2q1::Vector{Vq}

    IV::Vector{VI}                          # dynamics -I views [q2, γ1, b1]
    ITV::Vector{VIT}                        # dynamics -I views [q2, γ1, b1] transposed
    q0::Vector{Vq0}                         # dynamics q0 views
    q0T::Vector{Vq0T}                       # dynamics q0 views transposed
    q1::Vector{Vq1}                         # dynamics q1 views
    q1T::Vector{Vq1T}                       # dynamics q1 views transposed
    u1::Vector{Vu1}                         # dynamics u1 views
    u1T::Vector{Vu1T}                       # dynamics u1 views transposed
    reg::Vector{Vreg}                       # dual regularization views
end

function NewtonJacobianConfigurationForce(model::ContactModel, env::Environment, H::Int)
    dim = model.dim

    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
    nc = dim.c # contact
    nb = nc * friction_dim(env) # linear friction
    nd = nq + nc + nb # implicit dynamics constraint
    nr = nq + nu + nc + nb + nd # size of a one-time-step block

    off = 0
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iγ = SizedVector{nc}(off .+ (1:nc)); off += nc # index of the impact γ1
    ib = SizedVector{nb}(off .+ (1:nb)); off += nb # index of the linear friction b1
    iν = SizedVector{nd}(off .+ (1:nd)); off += nd # index of the dynamics lagrange multiplier ν1
    iz = vcat(iq, iγ, ib) # index of the IP solver solution [q2, γ1, b1]
    iθ = vcat(iq .- 2nr, iq .- nr, iu) # index of the IP solver data [q0, q1, u1]

    R = spzeros(H * nr, H * nr)

    obj_q2  = [view(R, (t - 1) * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H]
    obj_u1  = [view(R, (t - 1) * nr .+ iu, (t - 1) * nr .+ iu) for t = 1:H]
    obj_γ1  = [view(R, (t - 1) * nr .+ iγ, (t - 1) * nr .+ iγ) for t = 1:H]
    obj_b1  = [view(R, (t - 1) * nr .+ ib, (t - 1) * nr .+ ib) for t = 1:H]
    obj_q1q2  = [view(R, (t - 1) * nr .+ iq, t * nr .+ iq) for t = 1:H-1]
    obj_q2q1  = [view(R, t * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H-1]

    IV  = [view(R, CartesianIndex.((t - 1) * nr .+ iz, (t - 1) * nr .+ iν)) for t = 1:H]
    ITV = [view(R, CartesianIndex.((t - 1) * nr .+ iν, (t - 1) * nr .+ iz)) for t = 1:H]
    q0  = [view(R, (t - 1) * nr .+ iν, (t - 3) * nr .+ iq) for t = 3:H]
    q0T = [view(R, (t - 3) * nr .+ iq, (t - 1) * nr .+ iν) for t = 3:H]
    q1  = [view(R, (t - 1) * nr .+ iν, (t - 2) * nr .+ iq) for t = 2:H]
    q1T = [view(R, (t - 2) * nr .+ iq, (t - 1) * nr .+ iν) for t = 2:H]
    u1  = [view(R, (t - 1) * nr .+ iν, (t - 1) * nr .+ iu) for t = 1:H]
    u1T = [view(R, (t - 1) * nr .+ iu, (t - 1) * nr .+ iν) for t = 1:H]
    reg = [view(R, CartesianIndex.((t - 1) * nr .+ iν, (t - 1) * nr .+ iν)) for t = 1:H] # TODO: Cartesian indices to only grab diagonals

    return NewtonJacobianConfigurationForce{eltype(R),
        eltype.((obj_q2, obj_u1, obj_γ1, obj_b1))...,
        eltype.((IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg))...}(
        R, obj_q2, obj_u1, obj_γ1, obj_b1,
        obj_q1q2, obj_q2q1,
        IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg)
end

# Configurations
struct NewtonJacobianConfiguration{T,Vq,Vu,VI,VIT,Vq0,Vq0T,Vq1,Vq1T,Vu1,Vu1T,Vreg} <: NewtonJacobian
    R#::SparseMatrixCSC{T,Int}                 # jacobian

    obj_q2::Vector{Vq}                          # obj views
    obj_u1::Vector{Vu}                          # obj views

    obj_q1q2::Vector{Vq}
    obj_q2q1::Vector{Vq}

    IV::Vector{VI}                          # dynamics -I views [q2]
    ITV::Vector{VIT}                        # dynamics -I views [q2] transposed
    q0::Vector{Vq0}                         # dynamics q0 views
    q0T::Vector{Vq0T}                       # dynamics q0 views transposed
    q1::Vector{Vq1}                         # dynamics q1 views
    q1T::Vector{Vq1T}                       # dynamics q1 views transposed
    u1::Vector{Vu1}                         # dynamics u1 views
    u1T::Vector{Vu1T}                       # dynamics u1 views transposed
    reg::Vector{Vreg}                       # dual regularization views
end

function NewtonJacobianConfiguration(model::ContactModel, env::Environment, H::Int)
    dim = model.dim

    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
    nc = dim.c # contact
    nd = nq # implicit dynamics constraint
    nr = nq + nu + nd # size of a one-time-step block

    off = 0
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iν = SizedVector{nd}(off .+ (1:nd)); off += nd # index of the dynamics lagrange multiplier ν1
    iz = copy(iq)                                  # index of the IP solver solution [q2]
    iθ = vcat(iq .- 2nr, iq .- nr, iu) # index of the IP solver data [q0, q1, u1]

    R = spzeros(H * nr, H * nr)

    obj_q2  = [view(R, (t - 1) * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H]
    obj_u1  = [view(R, (t - 1) * nr .+ iu, (t - 1) * nr .+ iu) for t = 1:H]

    obj_q1q2  = [view(R, (t - 1) * nr .+ iq, t * nr .+ iq) for t = 1:H-1]
    obj_q2q1  = [view(R, t * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H-1]

    IV  = [view(R, CartesianIndex.((t - 1) * nr .+ iz, (t - 1) * nr .+ iν)) for t = 1:H]
    ITV = [view(R, CartesianIndex.((t - 1) * nr .+ iν, (t - 1) * nr .+ iz)) for t = 1:H]
    q0  = [view(R, (t - 1) * nr .+ iν, (t - 3) * nr .+ iq) for t = 3:H]
    q0T = [view(R, (t - 3) * nr .+ iq, (t - 1) * nr .+ iν) for t = 3:H]
    q1  = [view(R, (t - 1) * nr .+ iν, (t - 2) * nr .+ iq) for t = 2:H]
    q1T = [view(R, (t - 2) * nr .+ iq, (t - 1) * nr .+ iν) for t = 2:H]
    u1  = [view(R, (t - 1) * nr .+ iν, (t - 1) * nr .+ iu) for t = 1:H]
    u1T = [view(R, (t - 1) * nr .+ iu, (t - 1) * nr .+ iν) for t = 1:H]
    reg = [view(R, CartesianIndex.((t - 1) * nr .+ iν, (t - 1) * nr .+ iν)) for t = 1:H] # TODO: Cartesian indices to only grab diagonals

    return NewtonJacobianConfiguration{eltype(R),
        eltype.((obj_q2, obj_u1))...,
        eltype.((IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg))...}(
        R, obj_q2, obj_u1,
        obj_q1q2, obj_q2q1,
        IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg)
end

function NewtonJacobian(model::ContactModel, env::Environment, H::Int;
    mode = :configurationforce)

    if mode == :configurationforce
        return NewtonJacobianConfigurationForce(model, env, H)
    elseif mode == :configuration
        NewtonJacobianConfiguration(model, env, H)
    else
        @error "mode not implemented"
    end
end

function initialize_jacobian!(jac::NewtonJacobian, obj::Objective, H::Int)

    fill!(jac.R, 0.0)

    # Objective
    hessian!(jac, obj)

    for t = 1:H
        # Implicit dynamics
        jac.IV[t] .-= 1.0
        jac.ITV[t] .-= 1.0
    end

    return nothing
end

function update_jacobian!(jac::NewtonJacobian, im_traj::ImplicitTraj, obj::Objective,
    H::Int, β::T) where T

    # reset
    for t = 1:H
        if t >= 3
            fill!(jac.q0[t-2], 0.0)
            fill!(jac.q0T[t-2], 0.0)
        end

        if t >= 2
            fill!(jac.q1[t-1], 0.0)
            fill!(jac.q1T[t-1], 0.0)
        end

        fill!(jac.u1[t], 0.0)
        fill!(jac.u1T[t], 0.0)

        # Dual regularization
        fill!(jac.reg[t], 0.0)
    end

    for t = 1:H
        if t >= 3
            jac.q0[t-2]  .+= im_traj.δq0[t]
            jac.q0T[t-2] .+= im_traj.δq0[t]'
        end

        if t >= 2
            jac.q1[t-1]  .+= im_traj.δq1[t]
            jac.q1T[t-1] .+= im_traj.δq1[t]'
        end

        jac.u1[t]  .+= im_traj.δu1[t]
        jac.u1T[t] .+= im_traj.δu1[t]'

        # Dual regularization
        jac.reg[t] .-= 1.0 * β * im_traj.ip[t].κ # TODO sort the κ stuff, maybe make it a prameter of this function
    end

    return nothing
end

function jacobian!(jac::NewtonJacobian, im_traj::ImplicitTraj, obj::Objective,
    H::Int, β::T) where T

    initialize_jacobian!(jac, obj, H)
    update_jacobian!(jac, im_traj, obj, H, β)

    return nothing
end

function hessian!(hess::NewtonJacobianConfigurationForce, obj::TrackingObjective)
    for t = 1:length(obj.u)
        # Cost function
        hess.obj_q2[t] .+= obj.q[t]
        hess.obj_u1[t] .+= obj.u[t]
        hess.obj_γ1[t] .+= obj.γ[t]
        hess.obj_b1[t] .+= obj.b[t]
    end
end

function hessian!(hess::NewtonJacobianConfiguration, obj::TrackingObjective)
    for t = 1:length(obj.u)
        # Cost function
        hess.obj_q2[t] .+= obj.q[t]
        hess.obj_u1[t] .+= obj.u[t]
    end
end

function hessian!(hess::NewtonJacobianConfigurationForce, obj::TrackingVelocityObjective)
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

function hessian!(hess::NewtonJacobianConfiguration, obj::TrackingVelocityObjective)
    for t = 1:length(obj.u)
        # Cost function
        hess.obj_q2[t] .+= obj.q[t]
        hess.obj_u1[t] .+= obj.u[t]

        # velocity
        hess.obj_q2[t] .+= obj.v[t]
        t == 1 && continue
        hess.obj_q2[t-1] .+= obj.v[t]
        hess.obj_q1q2[t-1] .-= obj.v[t]
        hess.obj_q2q1[t-1] .-= obj.v[t]
    end
end
