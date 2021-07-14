const ContactControl = Main
include(joinpath(module_dir(), "src/dynamics/particle_2D/visuals.jl"))
vis = Visualizer()
open(vis)

################################################################################
# Parameters
################################################################################
# time
h = 0.01
T = 175

# initial conditions
q0 = SVector{s.model.dim.q}([0.00, 0.50])
q1 = SVector{s.model.dim.q}([0.05, 0.50])

# Simulation tolerance
tol = 1e-8

################################################################################
# Linearized Cone
################################################################################
s = get_simulation("particle_2D", "flat_2D_lc", "flat_lc")
s.model.μ_world = 0.5

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = tol, κ_tol = tol,
		diff_sol = false,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot_surface!(vis, s.env, xlims=[-10,10], ylims=[-0.3,0.3])
plot_lines!(vis, s.model, sim.traj.q, name = :LClines)
anim = visualize_robot!(vis, s.model, sim.traj, name = :LC)
visualize_force!(vis, s.model, s.env, sim.traj, anim=anim, shift=-0.25, name = :LC)

plt = plot(legend=false)
plot!(plt, hcat(Vector.([q[1:end] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u]*N_sample)...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ]*N_sample)...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b]*N_sample)...)', color=:black, linewidth=4.0)


################################################################################
# Non Linear Cone
################################################################################
s = get_simulation("particle_2D", "flat_2D_nc", "flat_nc")
s.model.μ_world = 0.5

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = tol, κ_tol = tol,
		diff_sol = false,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot_surface!(vis, s.env, xlims=[-10,10], ylims=[-0.3,0.3])
plot_lines!(vis, s.model, sim.traj.q, name = :NClines)
visualize_robot!(vis, s.model, sim.traj, name = :NC, anim = anim)
visualize_force!(vis, s.model, s.env, sim.traj, anim=anim, shift=-0.25, name = :NC)

plt = plot(legend=false)
plot!(plt, hcat(Vector.([q[1:end] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u]*N_sample)...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ]*N_sample)...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b]*N_sample)...)', color=:black, linewidth=4.0)







# build_robot!(vis, s.model, name=:ghost0, α=0.2)
# build_robot!(vis, s.model, name=:ghost1, α=0.4)
# build_robot!(vis, s.model, name=:ghost2, α=1.0)
#
# set_robot!(vis, s.model, sim.traj.q[30], name=:ghost0)
# set_robot!(vis, s.model, sim.traj.q[34], name=:ghost1)
# set_robot!(vis, s.model, sim.traj.q[38], name=:ghost2)

# filename = "particle_2D_mdp"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
