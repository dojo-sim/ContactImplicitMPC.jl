# Configurations and forces
struct NewtonJacobianConfigurationForce{T,Vq,Vu,Vγ,Vb,VI,VIT,Vq0,Vq0T,Vq1,Vq1T,Vu1,Vu1T,Vpr,Vdu} <: NewtonJacobian
    R::SparseMatrixCSC{T,Int}                 # jacobian
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
    reg_pr::Vpr                             # primal regularization views
    reg_du::Vdu                             # dual regularization views
end

function NewtonJacobianConfigurationForce(model::Model, env::Environment, H::Int)
    nq = model.nq # configuration
    nu = model.nu # control
    nw = model.nw # disturbance
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
    iθ = vcat(iq .- 2nr, iq .- nr, iu) # index of the IP solver data [q0, q1, u1]

    # iν = SizedVector{nd}(off .+ (1:nd)); off += nd # index of the dynamics lagrange multiplier ν1
    iν = SizedVector{nd}(1:nd) # index of the dynamics lagrange multiplier ν1

    R = spzeros(H * (nr + nd), H * (nr + nd))
    # R = zeros(H * (nr + nd), H * (nr + nd))

    obj_u1  = [view(R, (t - 1) * nr .+ iu, (t - 1) * nr .+ iu) for t = 1:H]
    obj_γ1  = [view(R, (t - 1) * nr .+ iγ, (t - 1) * nr .+ iγ) for t = 1:H]
    obj_b1  = [view(R, (t - 1) * nr .+ ib, (t - 1) * nr .+ ib) for t = 1:H]
    obj_q2  = [view(R, (t - 1) * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H]
    obj_q1q2  = [view(R, (t - 1) * nr .+ iq, t * nr .+ iq) for t = 1:H-1]
    obj_q2q1  = [view(R, t * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H-1]

    IV  = [view(R, CartesianIndex.((t - 1) * nr .+ iz, H * nr + (t - 1) * nd .+ iν)) for t = 1:H]
    ITV = [view(R, CartesianIndex.(H * nr + (t - 1) * nd .+ iν, (t - 1) * nr .+ iz)) for t = 1:H]
    q0  = [view(R, H * nr + (t - 1) * nd .+ iν, (t - 3) * nr .+ iq) for t = 3:H]
    q0T = [view(R, (t - 3) * nr .+ iq, H * nr + (t - 1) * nd .+ iν) for t = 3:H]
    q1  = [view(R, H * nr + (t - 1) * nd .+ iν, (t - 2) * nr .+ iq) for t = 2:H]
    q1T = [view(R, (t - 2) * nr .+ iq, H * nr + (t - 1) * nd .+ iν) for t = 2:H]
    u1  = [view(R, collect(H * nr + (t - 1) * nd .+ iν), collect((t - 1) * nr .+ iu)) for t = 1:H]
    u1T = [view(R, collect((t - 1) * nr .+ iu), collect(H * nr + (t - 1) * nd .+ iν)) for t = 1:H]
    reg_pr = view(R, CartesianIndex.(1:H*nr, 1:H*nr))
    reg_du = view(R, CartesianIndex.(H * nr .+ (1:H * nd), H * nr .+ (1:H * nd)))

    return NewtonJacobianConfigurationForce(
        R, obj_q2, obj_u1, obj_γ1, obj_b1,
        obj_q1q2, obj_q2q1,
        IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg_pr, reg_du)
end

# Configurations
struct NewtonJacobianConfiguration{T,Vq,Vu,VI,VIT,Vq0,Vq0T,Vq1,Vq1T,Vu1,Vu1T,Vpr,Vdu} <: NewtonJacobian
    R::SparseMatrixCSC{T,Int}                 # jacobian
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
    reg_pr::Vpr                             # primal regularization views
    reg_du::Vdu                             # dual regularization views
end

function NewtonJacobianConfiguration(model::Model, env::Environment, H::Int)

    nq = model.nq # configuration
    nu = model.nu # control
    nw = model.nw # disturbance
    nc = model.nc # contact
    nd = nq # implicit dynamics constraint
    nr = nq + nu # size of a one-time-step block

    off = 0
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iz = copy(iq)                                  # index of the IP solver solution [q2]
    iθ = vcat(iq .- 2nr, iq .- nr, iu) # index of the IP solver data [q0, q1, u1]

    iν = SizedVector{nd}(1:nd) # index of the dynamics lagrange multiplier ν1

    R = spzeros(H * (nr + nd), H * (nr + nd))
    # R = zeros(H * (nr + nd), H * (nr + nd))

    obj_u1  = [view(R, (t - 1) * nr .+ iu, (t - 1) * nr .+ iu) for t = 1:H]
    obj_q2  = [view(R, (t - 1) * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H]

    obj_q1q2  = [view(R, (t - 1) * nr .+ iq, t * nr .+ iq) for t = 1:H-1]
    obj_q2q1  = [view(R, t * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H-1]

    IV  = [view(R, CartesianIndex.((t - 1) * nr .+ iz, H * nr + (t - 1) * nd .+ iν)) for t = 1:H]
    ITV = [view(R, CartesianIndex.(H * nr + (t - 1) * nd .+ iν, (t - 1) * nr .+ iz)) for t = 1:H]
    q0  = [view(R, H * nr + (t - 1) * nd .+ iν, (t - 3) * nr .+ iq) for t = 3:H]
    q0T = [view(R, (t - 3) * nr .+ iq, H * nr + (t - 1) * nd .+ iν) for t = 3:H]
    q1  = [view(R, H * nr + (t - 1) * nd .+ iν, (t - 2) * nr .+ iq) for t = 2:H]
    q1T = [view(R, (t - 2) * nr .+ iq, H * nr + (t - 1) * nd .+ iν) for t = 2:H]
    u1  = [view(R, collect(H * nr + (t - 1) * nd .+ iν), collect((t - 1) * nr .+ iu)) for t = 1:H]
    u1T = [view(R, collect((t - 1) * nr .+ iu), collect(H * nr + (t - 1) * nd .+ iν)) for t = 1:H]
    reg_pr = view(R, CartesianIndex.(1:H*nr, 1:H*nr))
    reg_du = view(R, CartesianIndex.(H * nr .+ (1:H * nd), H * nr .+ (1:H * nd)))

    return NewtonJacobianConfiguration(
        R, obj_q2, obj_u1,
        obj_q1q2, obj_q2q1,
        IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg_pr, reg_du)
end

function NewtonJacobian(model::Model, env::Environment, H::Int;
    mode = :configurationforce)

    if mode == :configurationforce
        return NewtonJacobianConfigurationForce(model, env, H)
    elseif mode == :configuration
        NewtonJacobianConfiguration(model, env, H)
    else
        @error "mode not implemented"
    end
end

function initialize_jacobian!(jac::NewtonJacobian, obj::Objective, H::Int; update_hessian::Bool=true)

    fill!(jac.R, 0.0)
    # fill!(jac.R[301:480,:], 0.0)
    # fill!(jac.R[:,301:480], 0.0)

    # Objective
    update_hessian && hessian!(jac, obj)

    for t = 1:H
        # Implicit dynamics
        jac.IV[t] .-= 1.0
        jac.ITV[t] .-= 1.0
    end

    return nothing
end

function update_jacobian!(jac::NewtonJacobian, im_traj::ImplicitTrajectory, obj::Objective,
    H::Int, β::T, window::Vector{Int}) where T

    for (i, t) in enumerate(window[1:end-2])
        if i >= 3
            jac.q0[i-2]  .+= im_traj.δq0[t]
            jac.q0T[i-2] .+= transpose(im_traj.δq0[t])
        end

        if i >= 2
            jac.q1[i-1]  .+= im_traj.δq1[t]
            jac.q1T[i-1] .+= transpose(im_traj.δq1[t])
        end

        jac.u1[i]  .+= im_traj.δu1[t]
        jac.u1T[i] .+= transpose(im_traj.δu1[t])

        # Dual regularization
        # jac.reg_pr .+= 1.0 * β * im_traj.ip[t].κ[1]
        jac.reg_du .-= β * im_traj.ip[t].κ[1]
    end

    return nothing
end

function jacobian!(jac::NewtonJacobian, im_traj::ImplicitTrajectory, obj::Objective,
        H::Int, β::T, window::Vector{Int}; update_hessian::Bool=true) where T

    initialize_jacobian!(jac, obj, H; update_hessian=update_hessian)
    update_jacobian!(jac, im_traj, obj, H, β, window)

    return nothing
end

function hessian!(jac::NewtonJacobianConfigurationForce, obj::TrackingObjective)
    for t = 1:length(obj.u)
        # Cost function
        jac.obj_q2[t] .+= obj.q[t]
        jac.obj_u1[t] .+= obj.u[t]
        jac.obj_γ1[t] .+= obj.γ[t]
        jac.obj_b1[t] .+= obj.b[t]
    end
end

function hessian!(jac::NewtonJacobianConfiguration, obj::TrackingObjective)
    for t = 1:length(obj.u)
        # Cost function
        jac.obj_q2[t] .+= obj.q[t]
        jac.obj_u1[t] .+= obj.u[t]
    end
end

function hessian!(jac::NewtonJacobianConfigurationForce, obj::TrackingVelocityObjective)
    for t = 1:length(obj.u)
        # Cost function
        jac.obj_q2[t] .+= obj.q[t]
        jac.obj_u1[t] .+= obj.u[t]
        jac.obj_γ1[t] .+= obj.γ[t]
        jac.obj_b1[t] .+= obj.b[t]

        # velocity
        jac.obj_q2[t] .+= obj.v[t]
        t == 1 && continue
        jac.obj_q2[t-1] .+= obj.v[t]
        jac.obj_q1q2[t-1] .-= 1.0*obj.v[t]
        jac.obj_q2q1[t-1] .-= 1.0*obj.v[t]
    end
end

function hessian!(jac::NewtonJacobianConfiguration, obj::TrackingVelocityObjective)
    for t = 1:length(obj.u)
        # Cost function
        jac.obj_q2[t] .+= obj.q[t]
        jac.obj_u1[t] .+= obj.u[t]

        # velocity
        jac.obj_q2[t] .+= obj.v[t]
        t == 1 && continue
        jac.obj_q2[t-1] .+= obj.v[t]
        jac.obj_q1q2[t-1] .-= 1.0*obj.v[t]
        jac.obj_q2q1[t-1] .-= 1.0*obj.v[t]
    end
end
