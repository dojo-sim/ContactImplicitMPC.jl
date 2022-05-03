# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC
using LinearAlgebra
using Quaternions

@show Threads.nthreads()

include("continuous_policy_v2.jl")

# ## Simulation
s = get_simulation("centroidal_quadruped", "flat_3D_lc", "flat")
model = s.model
env = s.env

# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
	joinpath(module_dir(), "src/dynamics/centroidal_quadruped/gaits/inplace_trot_v4.jld2"),
    # joinpath(module_dir(), "src/dynamics/centroidal_quadruped/gaits/stand_euler_v0.jld2"),
    load_type = :split_traj_alt));


H = ref_traj.H
h = ref_traj.h

# ## MPC setup
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 3000
κ_mpc = 2.0e-4

v0 = 0.0
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
        max_time=1.0e-2,
		solver=:ldl_solver,
        threads=false,
        verbose=false,
        max_iter = 5),
    mpc_opts = CIMPCOptions(
		# live_plotting=true
		));

# ## Disturbances
idx_d1 = 20
idx_d2 = idx_d1 + 350
idx_d3 = idx_d2 + 350
idx_d4 = idx_d3 + 350
idx_d5 = idx_d4 + 350
idx_d6 = idx_d5 + 350
idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5, idx_d6]
# impulses = [[0,0,10.0], [0,0,-10.0], [0,5,0.0], [0,-5,0.0], [5,0,0.0], [-5,0,0.0]]
impulses = [1.0 * [15.0,25.0,30.0], [0,0,0.0], [0,0,0.0], [0,0,0.0], [0,0,0.0], [0,0,0.0]]
d = impulse_disturbances(impulses, idx)

w = [[0.0,0.0,0.0] for i=1:H_sim/N_sample]
# w = [[0,0.07,0.0] for i=1:H_sim/N_sample]
# w = [[0,0.0,0.0] for i=1:H_sim/N_sample]
# d = open_loop_disturbances(w, N_sample)

# ## Initial conditions
q1_sim, v1_sim = initial_conditions(ref_traj);

# ## Simulator
# u_ref = deepcopy(ref_traj.u[1])
# p = open_loop_policy([u_ref for i=1:1000], N_sample=N_sample)
# p = open_loop_policy([(fill(ref_traj.u, 100)...)...], N_sample=N_sample)
# sim = simulator(s, H_sim, h=h_sim, dist=d);
sim = simulator(s, H_sim, h=h_sim, policy=p, dist=d);


using BenchmarkTools
# ## Simulate
q1_sim0 = deepcopy(q1_sim)
# q1_sim0[4:6] += [0.02, 0.04, 0.02]
Δx = 0.0*[0.1, 0.2, 0.1]
q1_sim0[1:3] += Δx
q1_sim0[7:9] += Δx
q1_sim0[10:12] += Δx
q1_sim0[13:15] += Δx
q1_sim0[16:18] += Δx

RoboDojo.simulate!(sim, q1_sim0, v1_sim)

# # ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.open(vis)

# ## Visualize
set_light!(vis)
set_floor!(vis, grid=true)
set_background!(vis)
anim = visualize!(vis, model, sim.traj.q; Δt=h_sim)

# u_force = impulses[1]
# u_force ./= norm(u_force)
# u_force .*= 0.35
# force_vis = ArrowVisualizer(vis[:force])
# setobject!(force_vis, MeshPhongMaterial(color=RGBA(1.0, 0.0, 0.0, 1.0)))
# settransform!(force_vis,
#             Point(sim.traj.q[1][1], sim.traj.q[1][2], sim.traj.q[1][3]),
#             Vec(u_force[1], u_force[2], u_force[3]),
#             shaft_radius=0.035,
#             max_head_radius=0.075)
# setvisible!(vis[:force], true)
# for t = 1:length(sim.traj.q)
#     MeshCat.atframe(anim, t) do
#         t < 75 ? setvisible!(vis[:force], true) : setvisible!(vis[:force], false)
#     end
# end
# MeshCat.setanimation!(vis, anim)


# 12 / 6000 == 0.01 / 5

# thr = 1.0e-3
# N = 18
# H_ref = length(ref_traj.q)
# foot1_ref = [([[q[6 .+ 3] > thr ? 0 : 1 for q in ref_traj.q] for t = 1:N]...)...]
# foot2_ref = [([[q[9 .+ 3] > thr ? 0 : 1 for q in ref_traj.q] for t = 1:N]...)...]
# foot3_ref = [([[q[12 .+ 3] > thr ? 0 : 1 for q in ref_traj.q] for t = 1:N]...)...]
# foot4_ref = [([[q[15 .+ 3] > thr ? 0 : 1 for q in ref_traj.q] for t = 1:N]...)...]
# time_ref = range(0, stop=h * (N * H_ref - 1), length=N * H_ref)

plt = plot(layout=(4, 1));

plt = plot!(plt, time_ref, foot1_ref,
    linetype=:steppost, color=:orange, width=3.0, label="", yaxis=nothing, subplot=1)
plt = plot!(plt, time_ref, foot2_ref,
    linetype=:steppost, color=:lightgreen, width=3.0, label="", yaxis=nothing, subplot=2)
plt = plot!(plt, time_ref, foot3_ref,
    linetype=:steppost, color=:cyan, width=3.0, label="", yaxis=nothing, subplot=3)
plt = plot!(plt, time_ref, foot4_ref,
    linetype=:steppost, color=:magenta, width=3.0, label="", yaxis=nothing, subplot=4)


# foot1_sim = [q[6 .+ 3] > thr ? 0.01 : 0.99 for q in sim.traj.q]
# foot2_sim = [q[9 .+ 3] > thr ? 0.01 : 0.99 for q in sim.traj.q]
# foot3_sim = [q[12 .+ 3] > thr ? 0.01 : 0.99 for q in sim.traj.q]
# foot4_sim = [q[15 .+ 3] > thr ? 0.01 : 0.99 for q in sim.traj.q]
# time_sim = range(0, stop=h_sim * (H_sim - 1), length=H_sim)

# plt = plot!(plt, time_sim, foot1_sim[2:end-1],
#     linetype=:steppost, color=:black, width=2.0, label="", yaxis=nothing, subplot=1)
# plt = plot!(time_sim, foot2_sim[2:end-1],
#     linetype=:steppost, color=:black, width=2.0, label="", yaxis=nothing, subplot=2)
# plt = plot(plt, time_sim, foot3_sim[2:end-1],
#     linetype=:steppost, color=:black, width=2.0, label="", yaxis=nothing, subplot=3)
# plt = plot!(plt, time_sim, foot4_sim[2:end-1],
#     linetype=:steppost, color=:black, width=2.0, label="", yaxis=nothing, subplot=4, xlabel="time (s)")

# display(plt)
# savefig(plt, "/home/taylor/Downloads/centroidal_quadruped_push_recovery.png")

# # # ## Timing result
# # # Julia is [JIT-ed](https://en.wikipedia.org/wiki/Just-in-time_compilation) so re-run the MPC setup through Simulate for correct timing results.
process!(sim.stats, N_sample) # Time budget
H_sim * h_sim / sum(sim.stats.policy_time) # Speed ratio
plot(sim.stats.policy_time, xlabel="timestep", ylabel="mpc time (s)", label="", linetype=:steppost)

# # convert_frames_to_video_and_gif("centroidal_push_recovery")













# sim.policy.ref_traj
# for i = 1:100
# 	@show 1 + Int(floor((i-1)/N_sample)) % 10
# end


# sim.traj
# sim.policy.ref_traj
