model = get_model("quadruped")
κ = 1.0e-4
ref_traj = get_trajectory("quadruped", "gait1")
ref_traj.κ .= κ
H = ref_traj.H
h = 0.1
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# Test Jacobian!
cost = CostFunction(H, model.dim,
    q = [Diagonal(1.0e-2 *
        ([0.02, 0.02, 1.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15])) for t = 1:H],
    u = [Diagonal(3.0e-2 * ones(nu)) for t = 1:H],
    γ = [Diagonal(1.0e-6 * ones(nc)) for t = 1:H],
    b = [Diagonal(1.0e-6 * ones(nb)) for t = 1:H])

im_traj = ImplicitTraj(ref_traj, model)
num_var_mp = H * (nq + nu + nc + nb) + H * nd
num_data_mp = 0

τ0 = vcat([[ref_traj.q[t+2]; ref_traj.u[t]; ref_traj.γ[t]; ref_traj.b[t]] for t = 1:H]...)
ν0 = zeros(H * nd)
x0 = [τ0; ν0]

opts = InteriorPointOptions(
        κ_init = 1.0,
        κ_tol = 1.0,
        r_tol = 1.0e-8,
        diff_sol = false)

ip = interior_point(x0, zeros(num_data_mp),
         idx_ineq = collect(1:0),
         opts = opts)

mp_traj = trajectory_x(model, ip.z, ref_traj.q[1], ref_traj.q[2], H, h)

mp_traj.x .= [τ0; ν0] + 0.1 * randn(num_var_mp)
update_z!(mp_traj)
update_θ!(mp_traj)

implicit_dynamics!(im_traj, model, mp_traj)

function NewtonResidual(r, H::Int, dim::Dimensions)
    nq = dim.q # configuration
    nu = dim.u # control
    nc = dim.c # contact
    nb = dim.b # linear friction
    nd = nq + nc + nb # implicit dynamics constraint
    nr = nq + nu + nc + nb + nd # size of a one-time-step block

    vr = view(r, :)
    off = 0
    iq = SizedVector{nq}(off .+ (1:nq)); off += nq # index of the configuration q2
    iu = SizedVector{nu}(off .+ (1:nu)); off += nu # index of the control u1
    iγ = SizedVector{nc}(off .+ (1:nc)); off += nc # index of the impact γ1
    ib = SizedVector{nb}(off .+ (1:nb)); off += nb # index of the linear friction b1
    iν = SizedVector{nd}(off .+ (1:nd)); off += nd # index of the dynamics lagrange multiplier ν1
    iz = vcat(iq, iγ, ib) # index of the IP solver solution [q2, γ1, b1]

    q2  = [view(vr, (t - 1) * nr .+ iq) for t = 1:H]
    u1  = [view(vr, (t - 1) * nr .+ iu) for t = 1:H]
    γ1  = [view(vr, (t - 1) * nr .+ iγ) for t = 1:H]
    b1  = [view(vr, (t - 1) * nr .+ ib) for t = 1:H]

    rd  = [view(vr, (t - 1) * nr .+ iν) for t = 1:H]
    rI  = [view(vr, (t - 1) * nr .+ iz) for t = 1:H]
    q0  = [view(vr, (t - 3) * nr .+ iq) for t = 3:H]
    q1  = [view(vr, (t - 2) * nr .+ iq) for t = 2:H]

    T = eltype(r)

    return NewtonResidual{T, eltype.((q2, u1, γ1, b1))...,eltype.((rd, rI, q0, q1))...}(
        1.0, q2, u1, γ1, b1, rd, rI, q0, q1)
end

function NewtonJacobian(R, H::Int, dim::Dimensions)
    nq = dim.q # configuration
    nu = dim.u # control
    nw = dim.w # disturbance
    nc = dim.c # contact
    nb = dim.b # linear friction
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

    obj_q2  = [view(R, (t - 1) * nr .+ iq, (t - 1) * nr .+ iq) for t = 1:H]
    obj_u1  = [view(R, (t - 1) * nr .+ iu, (t - 1) * nr .+ iu) for t = 1:H]
    obj_γ1  = [view(R, (t - 1) * nr .+ iγ, (t - 1) * nr .+ iγ) for t = 1:H]
    obj_b1  = [view(R, (t - 1) * nr .+ ib, (t - 1) * nr .+ ib) for t = 1:H]

    IV  = [view(R, (t - 1) * nr .+ iz, (t - 1) * nr .+ iν) for t = 1:H]
    ITV = [view(R, (t - 1) * nr .+ iν, (t - 1) * nr .+ iz) for t = 1:H]
    q0  = [view(R, (t - 1) * nr .+ iν, (t - 3) * nr .+ iq) for t = 3:H]
    q0T = [view(R, (t - 3) * nr .+ iq, (t - 1) * nr .+ iν) for t = 3:H]
    q1  = [view(R, (t - 1) * nr .+ iν, (t - 2) * nr .+ iq) for t = 2:H]
    q1T = [view(R, (t - 2) * nr .+ iq, (t - 1) * nr .+ iν) for t = 2:H]
    u1  = [view(R, (t - 1) * nr .+ iν, (t - 1) * nr .+ iu) for t = 1:H]
    u1T = [view(R, (t - 1) * nr .+ iu, (t - 1) * nr .+ iν) for t = 1:H]
    reg = [view(R, (t - 1) * nr .+ iν, (t - 1) * nr .+ iν) for t = 1:H] # TODO: Cartesian indices to only grab diagonals

    return NewtonJacobian{eltype(R),
        eltype.((obj_q2, obj_u1, obj_γ1, obj_b1))...,
        eltype.((IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg))...}(
        R, obj_q2, obj_u1, obj_γ1, obj_b1, IV, ITV, q0, q0T, q1, q1T, u1, u1T, reg)
end

function residual!(res::NewtonResidual, model::ContactDynamicsModel, core::Newton,
    im_traj::ImplicitTraj, traj, ref_traj::ContactTraj)

    # unpack
    opts = core.opts
    cost = core.cost
    # res.r .= 0.0

    for t = 1:traj.H
        # Cost function
        delta!(core.Δq[t], traj.q[t+2], ref_traj.q[t+2])
        delta!(core.Δu[t], traj.u[t], ref_traj.u[t])
        delta!(core.Δγ[t], traj.γ[t], ref_traj.γ[t])
        delta!(core.Δb[t], traj.b[t], ref_traj.b[t])

        res.q2[t] .+= cost.q[t] * core.Δq[t]
        res.u1[t] .+= cost.u[t] * core.Δu[t]
        res.γ1[t] .+= cost.γ[t] * core.Δγ[t]
        res.b1[t] .+= cost.b[t] * core.Δb[t]

        # Implicit dynamics
        res.rd[t] .+= im_traj.d[t]

        # Minus Identity term #∇qk1, ∇γk, ∇bk
        res.rI[t] .-= traj.ν[t]
        # Implicit function theorem part #∇qk_1, ∇qk, ∇uk
        # t >= 3 ? mul!(res.q0[t-2], impl.δq0[t]', ν[t]) : nothing
        # t >= 2 ? mul!(res.q1[t-1], impl.δq1[t]', ν[t]) : nothing
        # mul!(res.u1[t], impl.δu1[t]', ν[t])
        t >= 3 ? res.q0[t-2] .+= im_traj.δq0[t]' * traj.ν[t] : nothing
        t >= 2 ? res.q1[t-1] .+= im_traj.δq1[t]' * traj.ν[t] : nothing
        res.u1[t] .+= im_traj.δu1[t]' * traj.ν[t]
    end

    return nothing
end

function jacobian!(jac::NewtonJacobian, model::ContactDynamicsModel,
    core::Newton, im_traj::ImplicitTraj)

    # unpack
    H = length(im_traj.ip)
    cost = core.cost
    opts = core.opts

    fill!(jac.R, 0.0) # TODO: remove

    for t = 1:H
        # Cost function
        jac.obj_q2[t] .+= cost.q[t]
        jac.obj_u1[t] .+= cost.u[t]
        jac.obj_γ1[t] .+= cost.γ[t]
        jac.obj_b1[t] .+= cost.b[t]

        # Implicit dynamics
        jac.IV[t][diagind(jac.IV[t])]   .-= 1.0
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

        # Dual regularization
        # jac.reg[t][diagind(jac.reg[t])] .-= opts.β * im_traj.lin[t].κ # TODO sort the κ stuff, maybe make it a prameter of this function
    end

    return nothing
end

res = NewtonResidual(ip.r, H, model.dim)
res_cand = NewtonResidual(ip.r̄, H, model.dim)
jac = NewtonJacobian(ip.rz, H, model.dim)
core = Newton(H, h, model, cost = cost)

function r2!(r, z, θ, κ)
    # update motion planning trajectory
    mp_traj.x .= z
    update_z!(mp_traj)
    update_θ!(mp_traj)

    # compute implicit dynamics
    implicit_dynamics!(im_traj, model, core.traj, κ = core.traj.κ)
    residual!(res, model, core, im_traj, mp_traj, ref_traj)
end
τ0
r2!(nothing, copy([τ0; ν0]), nothing, nothing)
ip.r
res.q2
im_traj.d
res.r
norm(mp_traj.x - [τ0; ν0])
function r̄!(r, z, θ, κ)
    # update motion planning trajectory
    mp_traj.x .= z
    update_z!(mp_traj)
    update_θ!(mp_traj)

    # compute implicit dynamics
    implicit_dynamics!(im_traj, model, core.traj, κ = core.traj.κ)
    residual!(res_cand, model, core, im_traj, mp_traj, ref_traj)
end

r̄!(nothing, ip.z, nothing, nothing)
ip.r̄

function rz!(rz, z, θ)
    # dynamics are already updated
    jacobian!(jac, model, core, im_traj)
end

# nd = nq + nc + nb
# vz = [view(τ, (t - 1) * nd .+ (1:nd)) for t = 1:H]
# q0 = copy(ref_traj.q[1])
# q1 = copy(ref_traj.q[2])
# vq = [t == 1 ? view(q0, 1:nq) : (t == 2 ? view(q1, 1:nq) : view(τ, (t - 3) * nd .+ (1:nq))) for t = 1:H+2]
# vγ = [view(τ, (t - 1) * nd + nq .+ (1:nc)) for t = 1:H]
# vb = [view(τ, (t - 1) * nd + nq + nc .+ (1:nb)) for t = 1:H]
#
# norm(vcat(vz...) - τ)
# norm(vcat(vq...) - vcat(ref_traj.q...))
# norm(vcat(vγ...) - vcat(ref_traj.γ...))
# norm(vcat(vb...) - vcat(ref_traj.b...))


a = zeros(5)

function update_a()
    a .= 1.0
end
update_a()
a
