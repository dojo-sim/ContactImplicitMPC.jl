# simulation
s = get_simulation("particle", "quadratic_bowl_3D_lc", "quadratic")
s.model.μ_world = 0.1

# time
h = 0.01
T = 1000

# initial conditions
q1 = @SVector [1.0, 0.5, 2.0]
q0 = @SVector [1.1, 0.5, 2.0]

# simulator
sim = ContactControl.simulator(s, q0, q1, h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-6,
		diff_sol = true),
	sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
@time status = ContactControl.simulate!(sim)
@test status

plot(hcat(sim.traj.q[1:3:T]...)', label = ["x" "y" "z"], legend = :bottomleft)
@show ϕ_func(s.model, s.env, sim.traj.q[end])
@show s.env.surf(sim.traj.q[end])
@show sim.traj.q[end]

include(joinpath(module_dir(), "src/dynamics/particle/visuals.jl"))

vis = Visualizer()
render(vis)
# open(vis)
visualize!(vis, s.model, sim.traj.q,
	Δt = h, r = 0.1)
plot_surface!(vis, s.env)

@test abs(sim.traj.q[end][1]) < 0.05
@test abs(sim.traj.q[end][2]) < 0.05
@test abs(sim.traj.q[end][3]) < 0.001
