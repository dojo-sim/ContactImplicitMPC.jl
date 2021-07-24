const ContactControl = Main
include(joinpath(module_dir(), "src/dynamics/box_mrp/visuals.jl"))
vis = Visualizer()
# render(vis)
open(vis)

# s = get_simulation("box_mrp", "flat_3D_nc", "flat_nc")
s = get_simulation("box_mrp", "flat_3D_lc", "flat_lc")
s.model.μ_world = 1.0

# time
h = 0.01
T = 175

# initial conditions
r0 = [0.0; 0.0; 1.0]
v0 = [5.0; 0.0; 0.0]

q0 = SVector{s.model.dim.q}([r0; zeros(3)])
q1 = SVector{s.model.dim.q}([r0 + v0 * h; zeros(3)])

q0 = SVector{s.model.dim.q}([0.0, 0.50, 1.0, 0.1, 0.1, 0.0])
q1 = SVector{s.model.dim.q}([0.0, 0.55, 1.0, 0.1, 0.1, 0.0])

# simulator
soft_tol = 1e-2
hard_tol = 1e-8
sim = ContactControl.simulator(s, q0, q1, h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = soft_tol, κ_tol = soft_tol,
		diff_sol = false,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)
@test status
l = 10
plot_surface!(vis, s.env, xlims=[-1,1], ylims=[-l,l])
plot_lines!(vis, s.model, sim.traj.q)
anim = visualize_robot!(vis, s.model, sim.traj)
visualize_force!(vis, model, env, sim.traj, anim=anim, h=0.0)

build_robot!(vis, model, name=:ghost0, α=0.2)
build_robot!(vis, model, name=:ghost1, α=0.4)
build_robot!(vis, model, name=:ghost2, α=1.0)

set_robot!(vis, model, sim.traj.q[30], name=:ghost0)
set_robot!(vis, model, sim.traj.q[34], name=:ghost1)
set_robot!(vis, model, sim.traj.q[38], name=:ghost2)


filename = "box_soft_visuals"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)





################################################################################
# Just visuals
################################################################################

# # time
# h = 0.01
# T = 155
#
# # initial conditions
# r0 = [0.0; 0.0; 1.0]
# v0 = [5.0; 0.0; 0.0]
#
# q0 = SVector{s.model.dim.q}([r0; zeros(3)])
# q1 = SVector{s.model.dim.q}([r0 + v0 * h; zeros(3)])
#
# q0 = SVector{s.model.dim.q}([0.0, -8.50, 1.50, 0.10, 0.10, 0.00])
# q1 = SVector{s.model.dim.q}([0.0, -8.35, 1.47, 0.10, 0.10, 0.00])
#
# # simulator
# soft_tol = 1e-2
# hard_tol = 1e-8
# sim = ContactControl.simulator(s, q0, q1, h, T,
# 	ip_opts = ContactControl.InteriorPointOptions(
# 		r_tol = soft_tol, κ_tol = soft_tol,
# 		diff_sol = false,
# 		solver = :lu_solver),
# 	sim_opts = ContactControl.SimulatorOptions(warmstart = false))
#
# # simulate
# @time status = ContactControl.simulate!(sim)
# @test status
# plot_surface!(vis, s.env, xlims=[-2,3], ylims=[-10,8])
# plot_lines!(vis, s.model, sim.traj.q)
# anim = visualize_robot!(vis, s.model, sim.traj)
