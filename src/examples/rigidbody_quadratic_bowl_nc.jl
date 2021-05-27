s = get_simulation("rigidbody", "quadratic_bowl_3D_nc", "quadratic_bowl_nc")
s.model.μ_world = 0.1

# time
h = 0.01
T = 500

# Rn + quaternion space
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
			collect([(1:3)..., (8:num_var(s.model, s.env))...]),
			collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
			collect((4:7)),
			collect((4:6)))

# initial conditions
r0 = [0.5; 0.75; 2.0]
v0 = [7.5; 5.0; 0.0]

quat0 = [1.0; 0.0; 0.0; 0.0]
ω0 = [0.0; 0.0; 0.0]

q0 = SVector{s.model.dim.q}([r0; quat0])
q1 = SVector{s.model.dim.q}([r0 + v0 * h; 0.5 * h * L_multiply(quat0) * [sqrt((2.0 / h)^2.0 - ω0' * ω0); ω0]])

@assert norm(q0[4:7]) ≈ 1.0
@assert norm(q1[4:7]) ≈ 1.0
@assert ϕ_func(s.model, s.env, q1)[1] > 0.0

# simulator
sim = ContactControl.simulator(s, q0, q1, h, T,
	space = rq_space,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-6,
		diff_sol = true,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)
@test status

include(joinpath(module_dir(), "src/dynamics/rigidbody/visuals.jl"))
vis = Visualizer()
render(vis)
# open(vis)
visualize!(vis, s.model, sim.traj.q, Δt = h)
plot_surface!(vis, s.env, ylims = (-2, 2), xlims = (-2, 2))


@assert all([norm(q[4:7]) ≈ 1.0 for q in sim.traj.q])

plot(hcat(sim.traj.q...)[4:7, :]')
