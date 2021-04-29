# Reference trajectory
model = deepcopy(ContactControl.get_model("quadruped", surf = "flat"))
model.μ_world = 0.1

ref_traj = deepcopy(ContactControl.get_trajectory("quadruped", "jump_v1", load_type = :split_traj_alt))
ContactControl.update_friction_coefficient!(ref_traj, model)

H = ref_traj.H
h = ref_traj.h

for t = 1:H
	r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)
	@test norm(r) < 1.0e-4
end

ϕ_func(model, q1)
model.μ_world = 1.0

# initial conditions
N_sample = 2
h_sim = h / N_sample
q1 = SVector{model.dim.q}(ref_traj.q[2])
q0 = q1 - (ref_traj.q[2] - ref_traj.q[1]) / N_sample
# simulator
sim = ContactControl.simulator(model, q0, q1, h_sim, 2 * H,
	p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u], N_sample = N_sample),
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-8, κ_init = 1.0e-7, solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
status = ContactControl.simulate!(sim, verbose = false)

include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
vis = Visualizer()
render(vis)
anim = visualize_robot!(vis, model, ref_traj)
anim = visualize_robot!(vis, model, sim.traj)

anim = visualize_meshrobot!(vis, model, ref_traj)
anim = visualize_force!(vis, model, ref_traj, anim=anim, h=h_sim)

plot(hcat(ref_traj.γ...)', linetype = :steppost)
plot(hcat(ref_traj.b...)', linetype = :steppost)
plot(hcat(ref_traj.u...)', linetype = :steppost)

plot(hcat(ref_traj.q...)', labels = "")
