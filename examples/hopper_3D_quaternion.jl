s = get_simulation("hopper_3D_quaternion", "flat_3D_nc", "flat_nc")

# time
h = 0.1
T = 50

# Rn + quaternion space
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
			collect([(1:3)..., (8:num_var(s.model, s.env))...]),
			collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
			collect((4:7)),
			collect((4:6)))

# initial conditions
p0 = [0.0; 0.0; 0.5]
v0 = [0.0; 0.0; 0.0]
quat0 = [1.0; 0.0; 0.0; 0.0]
ω0 = [0.0; 1.0; 0.0]

q0 = SVector{s.model.dim.q}([p0; quat0; 0.5])
q1 = SVector{s.model.dim.q}([p0 + v0 * h; 0.5 * h * L_multiply(quat0) * [sqrt((2.0 / h)^2.0 - ω0' * ω0); ω0]; 0.5])

@assert norm(q0[4:7]) ≈ 1.0
@assert norm(q1[4:7]) ≈ 1.0

# simulator
sim = ContactImplicitMPC.simulator(s, q0, q1, h, T,
	space = rq_space,
	ip_opts = ContactImplicitMPC.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-6,
		diff_sol = true,
		solver = :lu_solver),
	sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactImplicitMPC.simulate!(sim)
@test status

include(joinpath(pwd(), "src/dynamics/hopper_3D_quaternion/visuals.jl"))
vis = Visualizer()
render(vis)
# open(vis)
visualize_robot!(vis, s.model, sim.traj.q, h = h)

@assert all([norm(q[4:7]) ≈ 1.0 for q in sim.traj.q])
