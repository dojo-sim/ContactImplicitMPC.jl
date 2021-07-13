s = get_simulation("box_mrp", "flat_3D_nc", "flat_nc")
s.model.μ_world = 1.0

# time
h = 0.01
T = 200

# initial conditions
r0 = [0.0; 0.0; 1.0]
v0 = [5.0; 0.0; 0.0]

q0 = SVector{s.model.dim.q}([r0; zeros(3)])
q1 = SVector{s.model.dim.q}([r0 + v0 * h; zeros(3)])

# simulator
sim = ContactControl.simulator(s, q0, q1, h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-6,
		diff_sol = false,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)
@test status

include(joinpath(module_dir(), "src/dynamics/box_mrp/visuals.jl"))
vis = Visualizer()
render(vis)
# open(vis)
visualize!(vis, s.model, sim.traj.q, Δt = h)

@assert all([norm(q[4:7]) ≈ 1.0 for q in sim.traj.q])

r0 = [0.0; 0.0; 0.5]
quat0 = [1.0; 0.0; 0.0; 0.0]
q0 = SVector{s.model.dim.q}([r0; quat0])

r1 = [0.0; 0.0; 0.5 * sqrt(2.0)]
_quat1 = UnitQuaternion(RotX(0.25 * π))
# quat1 = [1.0; 0.0; 0.0; 0.0]
quat1 = [_quat1.w; _quat1.x; _quat1.y; _quat1.z]
q1 = SVector{s.model.dim.q}([r1; quat1])

visualize!(vis, s.model, [q1], Δt = h)
