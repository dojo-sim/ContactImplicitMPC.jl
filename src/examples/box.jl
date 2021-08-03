s = get_simulation("box", "flat_3D_nc", "flat_nc")
s.model.μ_world = 1.0

# time
h = 0.1
T = 10

# Rn + quaternion space
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
			collect([(1:3)..., (8:num_var(s.model, s.env))...]),
			collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
			[collect((4:7))],
			[collect((4:6))])

# initial conditions
r0 = [-1.5; 0.0; 1.5]
v0 = [8.0; 0.0; 0.0]

quat0 = [1.0; 0.0; 0.0; 0.0]
_quat0 = UnitQuaternion(RotX(0.0 * π))
quat0 = [_quat0.w; _quat0.x; _quat0.y; _quat0.z]
ω0 = [0.0; 0.0; 0.0]

q0 = SVector{s.model.dim.q}([r0; quat0])
q1 = SVector{s.model.dim.q}([r0 + v0 * h; 0.5 * h * L_multiply(quat0) * [sqrt((2.0 / h)^2.0 - ω0' * ω0); ω0]])

@assert norm(q0[4:7]) ≈ 1.0
@assert norm(q1[4:7]) ≈ 1.0

# simulator
sim = ContactControl.simulator(s, q0, q1, h, T,
	space = rq_space,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-6,
		diff_sol = true,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time ContactControl.simulate!(sim)
@test all([norm(q[4:7]) ≈ 1.0 for q in sim.traj.q])
include(joinpath(module_dir(), "src/dynamics/box/visuals.jl"))
vis = Visualizer()
render(vis)
open(vis)
visualize!(vis, s.model, sim.traj.q, Δt = h)
settransform!(vis["/Cameras/default"],
		compose(Translation(0.0, -95.0, -1.0), LinearMap(RotY(0.0 * π) * RotZ(-π / 2.0))))
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 20)
x_range = range(-3.0, stop = 3.0, length = 100)
surface_points = [Point(x, 0.0, s.env.surf(x)) for x in x_range]
surface_mat = LineBasicMaterial(color=RGBA(0.0, 0.0, 0.0, 1.0), linewidth=2.5)
setobject!(vis[:lines][:surface], MeshCat.Line(surface_points, surface_mat))
