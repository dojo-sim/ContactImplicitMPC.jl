struct LMPResidual{vq2,vu1,vγ1,vb1,vz,vν}
    r                                 # residual

    q2::Vector{vq2}                    # rsd objective views
    u1::Vector{vu1}                    # rsd objective views
    γ1::Vector{vγ1}                    # rsd objective views
    b1::Vector{vb1}                    # rsd objective views

    z::Vector{vz}                         # rsd dynamics -I views [q2, γ1, b1]
    ν::Vector{vν}                         # rsd dynamics lagrange multiplier views
end

function LMPResidual(r, H::Int, dim::Dimensions)
    nq = dim.q # configuration
    nu = dim.u # control
    nc = dim.c # contact
    nb = dim.b # linear friction
    nd = nq + nc + nb # implicit dynamics constraint
    nr = nq + nu + nc + nb # size of a one-time-step block

    vr = view(r, :)

    off = 0
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iγ = SizedVector{nc}(off .+ (1:nc)); off += nc # index of the impact γ1
    ib = SizedVector{nb}(off .+ (1:nb)); off += nb # index of the linear friction b1
    iν = SizedVector{nd}(1:nd);
    iz = vcat(iq, iγ, ib) # index of the IP solver solution [q2, γ1, b1]

    q2  = [view(vr, (t - 1) * nr .+ iq) for t = 1:H]
    u1  = [view(vr, (t - 1) * nr .+ iu) for t = 1:H]
    γ1  = [view(vr, (t - 1) * nr .+ iγ) for t = 1:H]
    b1  = [view(vr, (t - 1) * nr .+ ib) for t = 1:H]

    z  = [view(vr, (t - 1) * nr .+ iz) for t = 1:H]

    ν  = [view(vr, nr * H + (t - 1) * nd .+ iν) for t = 1:H]

    return LMPResidual{eltype.((q2, u1, γ1, b1))...,
        eltype.((z, ν))...}(vr, q2, u1, γ1, b1, z, ν)
end

function residual!(res::LMPResidual, cost::CostFunction,
    im_traj::ImplicitTraj, traj, ref_traj::ContactTraj)

    fill!(res.r, 0.0) # TODO: maybe not necessary

    for t = 1:traj.H
        # Lagrangian
        # Cost function
        res.q2[t] .+= cost.q[t] * (traj.q[t+2] .- ref_traj.q[t+2])
        res.u1[t] .+= cost.u[t] * (traj.u[t] .- ref_traj.u[t])
        res.γ1[t] .+= cost.γ[t] * (traj.γ[t] .- ref_traj.γ[t])
        res.b1[t] .+= cost.b[t] * (traj.b[t] .- ref_traj.b[t])

        # Minus Identity term # ∇qk1, ∇γk, ∇bk
        res.z[t] .-= traj.ν[t]

        t >= 3 ? res.q2[t-2] .+= im_traj.δq0[t]' * traj.ν[t] : nothing
        t >= 2 ? res.q2[t-1] .+= im_traj.δq1[t]' * traj.ν[t] : nothing
                   res.u1[t] .+= im_traj.δu1[t]' * traj.ν[t]

        # Constraints
        res.ν[t] .+= im_traj.d[t] #TODO set: .=
    end

    return nothing
end

struct LMPJacobian{Vq,Vu,Vγ,Vb,VI,VIT,Vq0,Vq0T,Vq1,Vq1T,Vu1,Vu1T}
    R                 # jacobian

    obj_q2::Vector{Vq}                          # obj views
    obj_u1::Vector{Vu}                          # obj views
    obj_γ1::Vector{Vγ}                          # obj views
    obj_b1::Vector{Vb}                          # obj views

    IV::Vector{VI}                          # dynamics -I views [q2, γ1, b1]
    ITV::Vector{VIT}                        # dynamics -I views [q2, γ1, b1] transposed
    q0::Vector{Vq0}                         # dynamics q1 views
    q0T::Vector{Vq0T}                       # dynamics q1 views transposed
    q1::Vector{Vq1}                         # dynamics q1 views
    q1T::Vector{Vq1T}                       # dynamics q1 views transposed
    u1::Vector{Vu1}                         # dynamics u1 views
    u1T::Vector{Vu1T}                       # dynamics u1 views transposed
end

function LMPJacobian(R, H::Int, dim::Dimensions)
    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
    nc = dim.c # contact
    nb = dim.b # linear friction
    nd = nq + nc + nb # implicit dynamics constraint
    nr = nq + nu + nc + nb # size of a one-time-step block

    vR = view(R, :, :)

    off = 0
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iγ = SizedVector{nc}(off .+ (1:nc)); off += nc # index of the impact γ1
    ib = SizedVector{nb}(off .+ (1:nb)); off += nb # index of the linear friction b1
    iν = SizedVector{nd}(1:nd);
    iz = vcat(iq, iγ, ib) # index of the IP solver solution [q2, γ1, b1]

    obj_q2  = [view(vR, (t - 1) * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H]
    obj_u1  = [view(vR, (t - 1) * nr .+ iu, (t - 1) * nr .+ iu) for t = 1:H]
    obj_γ1  = [view(vR, (t - 1) * nr .+ iγ, (t - 1) * nr .+ iγ) for t = 1:H]
    obj_b1  = [view(vR, (t - 1) * nr .+ ib, (t - 1) * nr .+ ib) for t = 1:H]

    IV  = [view(vR, (t - 1) * nr .+ iz, nr * H + (t - 1) * nd .+ iν) for t = 1:H] #TODO: CartesianIndex
    ITV = [view(vR, nr * H + (t - 1) * nd .+ iν, (t - 1) * nr .+ iz) for t = 1:H]
    q0  = [view(vR, nr * H + (t - 1) * nd .+ iν, (t - 3) * nr .+ iq) for t = 3:H] #TODO check
    q0T = [view(vR, (t - 3) * nr .+ iq, nr * H + (t - 1) * nd .+ iν) for t = 3:H]
    q1  = [view(vR, nr * H + (t - 1) * nd .+ iν, (t - 2) * nr .+ iq) for t = 2:H] #TODO check
    q1T = [view(vR, (t - 2) * nr .+ iq, nr * H + (t - 1) * nd .+ iν) for t = 2:H]
    u1  = [view(vR, nr * H + (t - 1) * nd .+ iν, (t - 1) * nr .+ iu) for t = 1:H]
    u1T = [view(vR, (t - 1) * nr .+ iu, nr * H + (t - 1) * nd .+ iν) for t = 1:H]

    return LMPJacobian{
        eltype.((obj_q2, obj_u1, obj_γ1, obj_b1))...,
        eltype.((IV, ITV, q0, q0T, q1, q1T, u1, u1T))...}(
        vR, obj_q2, obj_u1, obj_γ1, obj_b1, IV, ITV, q0, q0T, q1, q1T, u1, u1T)
end

function jacobian!(jac::LMPJacobian, cost::CostFunction, im_traj::ImplicitTraj)
    # reset
    fill!(jac.R, 0.0) #TODO: maybe remove

    # unpack
    H = length(im_traj.ip)

    for t = 1:H
        # Cost function
        jac.obj_q2[t] .+= cost.q[t]
        jac.obj_u1[t] .+= cost.u[t]
        jac.obj_γ1[t] .+= cost.γ[t]
        jac.obj_b1[t] .+= cost.b[t]

        # Implicit dynamics
        jac.IV[t][diagind(jac.IV[t])] .-= 1.0 #TODO: replace with correct view
        jac.ITV[t][diagind(jac.ITV[t])] .-= 1.0

        # TODO: ^ perform only once
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
    end

    return nothing
end

struct LMPResidualCache <: Cache
    res::LMPResidual
    cost::CostFunction
    traj::Trajectory
    ref_traj::ContactTraj
    im_traj::ImplicitTraj
    model::ContactDynamicsModel
end

struct LMPJacobianCache <: Cache
    jac::LMPJacobian
    cost::CostFunction
    im_traj::ImplicitTraj
end

function r_lmp!(r, z, θ, κ, cache)
    # update motion planning trajectory
    cache.traj.x .= z
    update_z!(cache.traj) #TODO fix views
    update_θ!(cache.traj)

    # compute implicit dynamics
    implicit_dynamics!(cache.im_traj, cache.model, cache.traj, κ = cache.traj.κ)
    residual!(cache.res, cache.cost, cache.im_traj, cache.traj, cache.ref_traj)
end

function rz_lmp!(rz, z, θ, cache)
    # dynamics are already updated
    jacobian!(cache.jac, cache.cost, cache.im_traj)
end

struct LMPSolver
    num_var::Int
    num_pr::Int
    num_du::Int
    num_data::Int
    traj
    ip::InteriorPoint
end

function linear_motion_planning_solver(model, ref_traj, cost;
    κ = 1.0e-4,
    opts = InteriorPointOptions(
            κ_init = 1.0,
            κ_tol = 1.0,
            r_tol = 1.0e-3,
            max_iter_outer = 1,
            res_norm = 2,
            reg = false,
            diff_sol = false))

    # time
    H = ref_traj.H
    h = ref_traj.h

    # model dimensions
    nq = model.dim.q
    nu = model.dim.u
    nc = model.dim.c
    nb = model.dim.b

    # time step sizes
    nd = nq + nc + nb
    nr = nq + nu + nc + nb + nd

    # planning sizes
    num_pr = H * (nq + nu + nc + nb)
    num_du = H * nd
    num_var = num_pr_mp + num_du_mp
    num_data = 0

    # implicit dynamics
    im_traj = ImplicitTraj(ref_traj, model)

    # initial guess
    τ0 = vcat([[ref_traj.q[t+2]; ref_traj.u[t]; ref_traj.γ[t]; ref_traj.b[t]] for t = 1:H]...)
    ν0 = zeros(H * nd)
    x0 = [τ0; ν0]

    # interior-point solver
    ip = interior_point(copy(x0), zeros(num_data),
             idx_ineq = collect(1:0),
             idx_pr = collect(1:num_pr),
             idx_du = collect(num_pr .+ (1:num_du)),
             opts = opts)

    # motion-planning trajectory
    mp_traj = trajectory_x(model, ip.z, ref_traj.q[1], ref_traj.q[2], H, h)

    # residual
    res = LMPResidual(ip.r, H, model.dim)
    res_cand = LMPResidual(ip.r̄, H, model.dim)

    # Jacobian
    jac = LMPJacobian(ip.rz, H, model.dim)

    # caches
    r_cache = LMPResidualCache(res, cost, mp_traj, ref_traj, im_traj, model)
    r̄_cache = LMPResidualCache(res_cand, cost, mp_traj, ref_traj, im_traj, model)
    rz_cache = LMPJacobianCache(jac, cost, im_traj)

    ip.r_cache = r_cache
    ip.r̄_cache = r̄_cache
    ip.rz_cache = rz_cache

    # methods
    ip.methods.r! = r_lmp!
    ip.methods.rz! = rz_lmp!

    LMPSolver(num_var, num_pr, num_du, num_data, mp_traj, ip)
end

function lmp!(lmp_solver::LMPSolver)
    interior_point!(lmp_solver.ip)
end
