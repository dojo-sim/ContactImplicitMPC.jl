model = ContactControl.get_model("particle")

# time
h = 0.01
T = 100

## DROP
# initial conditions
q0 = @SVector [0.0, 0.0, 0.0]
q1 = @SVector [0.0, 0.0, 0.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-8),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
status = ContactControl.simulate!(sim)
@test status

include(joinpath(pwd(), "src/dynamics/particle/visuals.jl"))

vis = Visualizer()
render(vis)
visualize!(vis, model, sim.traj.q,
	Δt = h, r = 0.1)
