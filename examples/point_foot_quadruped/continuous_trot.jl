# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC
using LinearAlgebra
using Quaternions
using BenchmarkTools

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.open(vis)

@show Threads.nthreads()

# include("continuous_policy.jl")
include("continuous_policy_v2.jl")

# ## Simulation
s = get_simulation("point_foot_quadruped", "flat_3D_lc", "flat")
model = s.model
env = s.env

# ref_traj.z
# ref_traj.θ
# ref_traj.q
# ref_traj.u
# ref_traj.H
# κ = 2e-4
# N = 200
# z = [[deepcopy(ref_traj.z) for i = 1:N]...;]
# θ = [[deepcopy(ref_traj.θ) for i = 1:N]...;]
# rz = [LinearizedStep(s, z[t], θ[t], κ).rz for t = 1:N*H]
# rθ = [LinearizedStep(s, z[t], θ[t], κ).rθ for t = 1:N*H]
# idx_z = indices_z(s)
# idx_θ = indices_θ(model, nf=1)
# ∂z∂θ = [-rz[t] \ rθ[t] for t = 1:N*H]
# ∂q3∂q1 = [∂z∂θ[t][idx_z.q, idx_θ.q1] for t = 1:N*H]
# ∂q3∂q2 = [∂z∂θ[t][idx_z.q, idx_θ.q2] for t = 1:N*H]
# ∂q2∂q1 = [zeros(model.nq, model.nq) for t = 1:N*H]
# ∂q2∂q2 = [I(model.nq) for t = 1:N*H]
# ∂q2∂u = [zeros(model.nq, model.nu) for t = 1:N*H]
# ∂q3∂u = [∂z∂θ[t][idx_z.q, idx_θ.u] for t = 1:N*H]
#
# A = [[∂q2∂q1[t] ∂q3∂q1[t]; ∂q2∂q2[t] ∂q3∂q2[t]] for t = 1:N*H]
# B = [[∂q2∂u[t]; ∂q3∂u[t]] for t = 1:N*H]
# Q = [0.01 * I(2model.nq) for t = 1:N*H+1]
# R = [1000 * I(model.nu) for t = 1:N*H]
#
# A = [[1 1; 0 sin(2π*t/1000)] for t = 1:N*H]
# B = [[0 1.]' for t = 1:N*H]
# Q = [1 * I(2) for t = 1:N*H+1]
# R = [1 * I(1) for t = 1:N*H]
#
# K, P = tvlqr(A, B, Q, R)


# plot(hcat([reshape(K[t], 2) for t = 1:N*H]...)')
# function tvlqr(A, B, Q, R)
#     T = length(Q)
#
#     P = [zero(A[1]) for t = 1:T]
#     K = [zero(B[1]') for t = 1:T-1]
#     P[T] = Q[T]
#
#     for t = T-1:-1:1
#         K[t] = (R[t] + B[t]' * P[t+1] *  B[t]) \ (B[t]' * P[t+1] * A[t])
#         P[t] = (Q[t] + K[t]' * R[t] * K[t]
#                 + (A[t] - B[t] * K[t])' * P[t+1] * (A[t] - B[t] * K[t]))
#     end
#
#     return K, P
# end


# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
# joinpath(module_dir(), "src/dynamics/point_foot_quadruped/gaits/inplace_trot_v0.jld2"),
	joinpath(module_dir(), "src/dynamics/point_foot_quadruped/gaits/inplace_trot_v1.jld2"),
    load_type = :split_traj_alt));


H = ref_traj.H
h = ref_traj.h

# ## MPC setup
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 4000
κ_mpc = 2.0e-4

v0 = 0.2
obj = TrackingVelocityObjective(model, env, H_mpc,
    v = [Diagonal(1e-3 * [[1,1,1]; 1e+3*[1,1,1]; fill([1,1,1], 4)...]) for t = 1:H_mpc],
	q = [relative_state_cost(1e-0*[1e-2,1e-2,1], 3e-1*[1,1,1], 1e-0*[0.2,0.2,1]) for t = 1:H_mpc],
	u = [Diagonal(3e-3 * vcat(fill([1,1,1], 4)...)) for t = 1:H_mpc],
	v_target = [1/ref_traj.h * [v0;0;0; 0;0;0; v0;0;0; v0;0;0; v0;0;0; v0;0;0] for t = 1:H_mpc],)

p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
	ip_opts = InteriorPointOptions(
					undercut = 5.0,
					κ_tol = κ_mpc,
					r_tol = 1.0e-4, # TODO maybe relax this parameter
					diff_sol = true,
					solver = :empty_solver,
					max_time = 1e5),
    n_opts = NewtonOptions(
        r_tol = 3e-5,
        max_time=10.0e-1,
		solver=:ldl_solver,
        threads=false,
        verbose=false,
        max_iter = 5),
    mpc_opts = CIMPCOptions(
		# live_plotting=true
		));

# ## Disturbances
w = [[0.0,0.0,0.0] for i=1:H_sim/N_sample]
d = open_loop_disturbances(w, N_sample)

# ## Initial conditions
q1_sim, v1_sim = initial_conditions(ref_traj);

# ## Simulator
sim = simulator(s, H_sim, h=h_sim, policy=p, dist=d);

# ## Simulate
q1_sim0 = deepcopy(q1_sim)
RoboDojo.simulate!(sim, q1_sim0, v1_sim)

# ## Visualize
set_light!(vis)
set_floor!(vis, grid=true)
set_background!(vis)
anim = visualize!(vis, model, sim.traj.q; Δt=h_sim)
# anim = visualize!(vis, model, ref_traj.q; Δt=h)

# # ## Timing result
# # Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio
plot(sim.stats.policy_time, xlabel="timestep", ylabel="mpc time (s)",
	ylims=[-0.001, 0.1],
	label="", linetype=:steppost)

sim
sim.traj

plt = plot()
plot!(plt, hcat(Vector.([(sim.traj.q[i+1][1:1] - sim.traj.q[i][1:1]) / sim.h for i=1:H_sim])...)')
plot!(plt, hcat(Vector.([(sim.traj.q[i+1][7:7] - sim.traj.q[i][7:7]) / sim.h for i=1:H_sim])...)')
plot!(plt, hcat(Vector.([(sim.traj.q[i+1][10:10] - sim.traj.q[i][10:10]) / sim.h for i=1:H_sim])...)')
plot!(plt, hcat(Vector.([(sim.traj.q[i+1][13:13] - sim.traj.q[i][13:13]) / sim.h for i=1:H_sim])...)')
plot!(plt, hcat(Vector.([(sim.traj.q[i+1][16:16] - sim.traj.q[i][16:16]) / sim.h for i=1:H_sim])...)')
