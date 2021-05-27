s_lc = ContactControl.get_simulation("particle", "flat_3D_lc", "flat")
s_nc = ContactControl.get_simulation("particle", "flat_3D_nc", "flat_nc")

# time
h = 0.01
T = 250

# initial conditions
q0 = @SVector [0.0, 0.0, 1.0]
# q1 = @SVector [0.075, 0.05, 1.0]
q1 = @SVector [0.0, 0.085, 1.0]

ip_opts = ContactControl.InteriorPointOptions(
	r_tol = 1.0e-8, κ_tol = 1.0e-8)

sim_opts = ContactControl.SimulatorOptions(warmstart = false)

sim_lc = ContactControl.simulator(s_lc, copy(q0), copy(q1), h, T,
	ip_opts = ip_opts,
	sim_opts = sim_opts)

sim_nc = ContactControl.simulator(s_nc, copy(q0), copy(q1), h, T,
	ip_opts = ip_opts,
	sim_opts = sim_opts)

# simulate
@test status_ls = ContactControl.simulate!(sim_lc)
@test status_nc = ContactControl.simulate!(sim_nc)

# visualize
include(joinpath(module_dir(), "src/dynamics/particle/visuals.jl"))
vis = Visualizer()
open(vis)

default_background!(vis)

r = 0.1
setobject!(vis["particle_lc"],
	GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
	convert(Float32, r)),
	MeshPhongMaterial(color = RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, 1.0)))

setobject!(vis["particle_nc"],
	GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
	convert(Float32, r)),
	MeshPhongMaterial(color = RGBA(51.0 / 255.0, 1.0, 1.0, 1.0)))

anim = MeshCat.Animation(convert(Int, floor(1.0 / h)))

for t = 1:T
	MeshCat.atframe(anim, t) do
		settransform!(vis["particle_lc"], MeshCat.Translation(sim_lc.traj.q[t][1:3]...))
		settransform!(vis["particle_nc"], MeshCat.Translation(sim_nc.traj.q[t][1:3]...))
	end
end

MeshCat.setanimation!(vis, anim)

setobject!(vis[:lc_line],
	MeshCat.Line([Point{3}(q) for q in sim_lc.traj.q],
	LineBasicMaterial(color=RGBA(1.0, 153.0 / 255.0, 51.0 / 255.0, 1.0), linewidth=5)))
setobject!(vis[:nc_line],
	MeshCat.Line([Point{3}(q) for q in sim_nc.traj.q],
	LineBasicMaterial(color=RGBA(51.0 / 255.0, 1.0, 1.0, 1.0), linewidth=5)))
setobject!(vis[:ref_line],
	MeshCat.Line([Point{3}([q0[1]; q0[2]; -0.025]), Point{3}([q1[1] * 75.0; q1[2] * 75.0; -0.025])],
	LineBasicMaterial(color=RGBA(0.0, 0.0, 0.0, 1.0), linewidth=5)))
