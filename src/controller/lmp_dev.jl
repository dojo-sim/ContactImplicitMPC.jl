include("motion_planning_trajectory.jl")
include("linearized_motion_planning.jl")

model = get_model("quadruped")
κ = 1.0e-4
ref_traj = get_trajectory("quadruped", "gait1")
ref_traj.κ .= κ
H = ref_traj.H
h = ref_traj.h
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
num_pr_mp = H * (nq + nu + nc + nb)
num_du_mp = H * nd
num_var_mp = num_pr_mp + num_du_mp
num_data_mp = 0

τ0 = vcat([[ref_traj.q[t+2]; ref_traj.u[t]; ref_traj.γ[t]; ref_traj.b[t]] for t = 1:H]...)
ν0 = zeros(H * nd)
x0 = [τ0; ν0]

opts = InteriorPointOptions(
        κ_init = 1.0,
        κ_tol = 1.0,
        r_tol = 1.0e-3,
        max_iter_outer = 1,
        res_norm = 2,
        reg = false,
        diff_sol = false)

ip = interior_point(copy(x0), zeros(num_data_mp),
         idx_ineq = collect(1:0),
         idx_pr = collect(1:num_pr_mp),
         idx_du = collect(num_pr_mp .+ (1:num_du_mp)),
         opts = opts)

mp_traj = trajectory_x(model, ip.z, ref_traj.q[1], ref_traj.q[2], H, h)
# x0 = [τ0; ν0] + 0.01 * randn(num_var_mp)
maximum([norm(mp_traj.q[t] - ref_traj.q[t]) for t = 1:length(ref_traj.q)])
maximum([norm(mp_traj.u[t] - ref_traj.u[t]) for t = 1:length(ref_traj.u)])
maximum([norm(mp_traj.γ[t] - ref_traj.γ[t]) for t = 1:length(ref_traj.γ)])
maximum([norm(mp_traj.b[t] - ref_traj.b[t]) for t = 1:length(ref_traj.b)])

res = LMPResidual(ip.r, H, model.dim)
res_cand = LMPResidual(ip.r̄, H, model.dim)
jac = LMPJacobian(ip.rz, H, model.dim)

r_cache = LMPResidualCache(res, cost, mp_traj, ref_traj, im_traj, model)
r̄_cache = LMPResidualCache(res_cand, cost, mp_traj, ref_traj, im_traj, model)
rz_cache = LMPJacobianCache(jac, cost, im_traj)

ip.methods.r! = r_lmp!
ip.methods.rz! = rz_lmp!
ip.r_cache = r_cache
ip.r̄_cache = r̄_cache
ip.rz_cache = rz_cache
#
# ip.z .= copy([τ0; ν0])
# r_lmp!(ip.r, ip.z, ip.θ, ip.κ[1], ip.r_cache)
# r_lmp!(ip.r̄, ip.z, ip.θ, ip.κ[1], ip.r̄_cache)
#
# ip.r_cache.traj.x
# ip.r_cache.traj.q
# ip.r_cache.im_traj.d
# implicit_dynamics!(im_traj, model, mp_traj)
# norm(mp_traj.θ[1][1:2nq + nu + nw + 1] - ref_traj.θ[1][1:2nq + nu + nw + 1])
# im_traj.d[1]
# ip.r_cache.traj
# norm(ip.r, Inf)
# norm(r_cache.res.r - ip.r)
# im_traj.ip[1].θ
# r_lmp!(ip.r̄, ip.z, ip.θ, ip.κ[1], ip.r̄_cache)
# norm(r̄_cache.res.r - ip.r̄)
#
# rz_lmp!(ip.rz, ip.z, ip.θ, ip.rz_cache)
# norm(ip.rz_cache.jac.R - ip.rz)
# rank(Array(ip.rz_cache.jac.R))
# rank(ip.rz)
# size(ip.rz)
# norm(ip.r)
# norm(ip.r̄)

ip.z .= copy([τ0; ν0])
ip.z[1:2] .+= 1.0e-2 * rand(2)
# re - ip.r
# ja - ip.rz
r_lmp!(ip.r, ip.z, ip.θ, ip.κ[1], ip.r_cache)
r_lmp!(ip.r̄, ip.z, ip.θ, ip.κ[1], ip.r̄_cache)
rz_lmp!(ip.rz, ip.z, ip.θ, ip.rz_cache)
norm(ip.r)
norm(ip.r̄)
ip.Δ .= ip.rz \ ip.r
α = 1.0
ip.z̄ .= ip.z - α * ip.Δ
r_lmp!(ip.r̄, ip.z̄, ip.θ, ip.κ[1], ip.r̄_cache)
norm(ip.r̄)
ip.z .= ip.z̄

status = interior_point!(ip)

solver = linear_motion_planning_solver(model, ref_traj, cost)
status = lmp!(solver)
