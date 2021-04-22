model = ContactControl.get_model("inverted_pendulum")

# time
h = 0.01
T = 100

# initial conditions
q0 = @SVector [-0.1]
q1 = @SVector [0.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-8),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
status = ContactControl.simulate!(sim)
@test status

# visualize
include(joinpath(@__DIR__, "..", "dynamics", "inverted_pendulum", "visuals.jl"))
vis = Visualizer()
render(vis)
add_walls!(vis, model)
anim = visualize_robot!(vis, model, sim.traj)
# anim = visualize_force!(vis, model, sim.traj, anim=anim, h=h_sim)
