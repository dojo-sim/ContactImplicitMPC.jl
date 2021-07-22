const ContactControl = Main
include(joinpath(module_dir(), "src/dynamics/particle_2D/visuals.jl"))
vis = Visualizer()
open(vis)

include(joinpath(module_dir(), "src", "solver", "interior_point.jl"))
include(joinpath(module_dir(), "src", "solver", "interior_point_latest.jl"))

################################################################################
# Parameters
################################################################################
# time
h = 0.01
T = 175

# initial conditions
q0 = SVector{2}([0.00, 0.50])
q1 = SVector{2}([0.05, 0.50])

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
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)

plot(hcat(Vector.(sim.traj.z)...)')
plot(hcat(Vector.(sim.traj.θ)...)')

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
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b])...)', color=:black, linewidth=4.0)

plot(hcat(Vector.(sim.traj.z)...)')
plot(hcat(Vector.(sim.traj.θ)...)')

# ################################################################################
# # Save Trajectory
# ################################################################################
# traj = deepcopy(sim.traj)
# gait_path = joinpath(@__DIR__, "..", "dynamics", "particle_2D", "gaits", "gait_NC.jld2")
# @save gait_path traj
#
# # Reload trajectory
# res = JLD2.jldopen(gait_path)
# loaded_traj = res["traj"]
#
# traj = deepcopy(get_trajectory(s.model, s.env,
# 	joinpath(module_dir(), "src/dynamics/particle_2D/gaits/gait_NC.jld2"),
# 	load_type = :joint_traj))
# plot(hcat(Vector.(traj.q)...)')


################################################################################
# Non Linear Cone (Mehrotra)
################################################################################
s = get_simulation("particle_2D", "flat_2D_nc", "flat_nc")
s.model.μ_world = 0.5

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	ip_type = :mehrotra_latest,
	ip_opts = Mehrotra113Options(
		max_iter_inner = 100,
		r_tol = tol,
		κ_tol = tol,
		κ_reg = 1e-3,
		γ_reg = 1.0,
		ls_relaxation = Inf,
		diff_sol = false,
		solver = :lu_solver,
		verbose = true,
		),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot_surface!(vis, s.env, xlims=[-10,10], ylims=[-0.3,0.3])
plot_lines!(vis, s.model, sim.traj.q, name = :NC_latestlines)
visualize_robot!(vis, s.model, sim.traj, name = :NC_latest, anim = anim)
visualize_force!(vis, s.model, s.env, sim.traj, anim=anim, shift=-0.25, name = :NC_latest)

plt = plot(legend=false)
plot!(plt, hcat(Vector.([q[1:end] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b])...)', color=:black, linewidth=4.0)

plot(hcat(Vector.(sim.traj.z)...)')
plot(hcat(Vector.(sim.traj.θ)...)')

################################################################################
# Dev
################################################################################
s = get_simulation("particle_2D", "flat_2D_nc", "flat_nc")
s.model.μ_world = 0.5
model = s.model
env = s.env

traj = deepcopy(get_trajectory(s.model, s.env,
	joinpath(module_dir(), "src/dynamics/particle_2D/gaits/gait_NC.jld2"),
	load_type = :joint_traj))

plot(hcat(Vector.(traj.z)...)')
plot(hcat(Vector.(traj.θ)...)')

t = 122
z = deepcopy(traj.z[t])
θ = deepcopy(traj.θ[t])
r_tol = 1e-8
κ_tol = 1e-8

ip = mehrotra_latest(z, θ,
	ix = linearization_var_index(model, env)[1],
	iy1 = linearization_var_index(model, env)[2],
	iy2 = linearization_var_index(model, env)[3],
	idyn = linearization_term_index(model, env)[1],
	irst = linearization_term_index(model, env)[2],
	ibil = linearization_term_index(model, env)[3],
	idx_ineq = inequality_indices(model, env),
	idx_soc = soc_indices(model, env),
	r! = s.res.r!,
	rm! = s.res.rm!,
	rz! = s.res.rz!,
	rθ! = s.res.rθ!,
	rz = s.rz,
	rθ = s.rθ,
	opts = Mehrotra113Options(
		max_iter_inner = 20,
		r_tol = r_tol,
		κ_tol = κ_tol,
		κ_reg = 1e-3,
		γ_reg = 1.0,
		ls_relaxation = Inf,
		verbose = true,
		))

# ip = interior_point(z, θ,
# 	ix = linearization_var_index(model, env)[1],
# 	iy1 = linearization_var_index(model, env)[2],
# 	iy2 = linearization_var_index(model, env)[3],
# 	idyn = linearization_term_index(model, env)[1],
# 	irst = linearization_term_index(model, env)[2],
# 	ibil = linearization_term_index(model, env)[3],
# 	idx_ineq = inequality_indices(model, env),
# 	idx_soc = soc_indices(model, env),
# 	r! = s.res.r!,
# 	rm! = s.res.rm!,
# 	rz! = s.res.rz!,
# 	rθ! = s.res.rθ!,
# 	rz = s.rz,
# 	rθ = s.rθ,
# 	opts = InteriorPointOptions(
# 		# max_iter_inner = 20,
# 		r_tol = r_tol,
# 		κ_tol = κ_tol,
# 		solver = :lu_solver,
# 		diff_sol = false,
# 		verbose = true,
# 		))

z0 = 1e-1rand(length(z))
z_initialize!(z0, model, env, deepcopy(traj.q[t+1]))
interior_point_solve!(ip, z0, θ)
norm(ip.z - deepcopy(traj.z[t]))
norm(ip.z)
norm(traj.z[t])

ip.r[ip.idyn]
ip.r[ip.irst]
ip.r[ip.ibil]
ip.Δ[ip.idyn]
ip.Δ[ip.irst]
ip.Δ[ip.ibil]

nq = s.model.dim.q
ip.θ[1:nq]
ip.θ[nq .+ (1:nq)]
ip.z[1:nq]

a = 10
a = 10
a = 10
a = 10













################################################################################
# Non Linear Cone (InteriorPoint Latest)
################################################################################
s = get_simulation("particle_2D", "flat_2D_nc", "flat_nc")
s.model.μ_world = 0.5

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	ip_type = :interior_point_latest,
	ip_opts = ContactControl.InteriorPoint113Options(
		r_tol = tol, κ_tol = tol,
		diff_sol = false,
		verbose = false,
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
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b])...)', color=:black, linewidth=4.0)

plot(hcat(Vector.(sim.traj.z)...)')
plot(hcat(Vector.(sim.traj.θ)...)')


################################################################################
# Dev
################################################################################
s = get_simulation("particle_2D", "flat_2D_nc", "flat_nc")
s.model.μ_world = 0.5
model = s.model
env = s.env

traj = deepcopy(get_trajectory(s.model, s.env,
	joinpath(module_dir(), "src/dynamics/particle_2D/gaits/gait_NC.jld2"),
	load_type = :joint_traj))

plot(hcat(Vector.(traj.z)...)')
plot(hcat(Vector.(traj.θ)...)')

t = 62
z = deepcopy(traj.z[t])
θ = deepcopy(traj.θ[t])
r_tol = 1e-8
κ_tol = 1e-8

ip = interior_point_latest(z, θ,
	ix = linearization_var_index(model, env)[1],
	iy1 = linearization_var_index(model, env)[2],
	iy2 = linearization_var_index(model, env)[3],
	idyn = linearization_term_index(model, env)[1],
	irst = linearization_term_index(model, env)[2],
	ibil = linearization_term_index(model, env)[3],
	idx_ineq = inequality_indices(model, env),
	idx_soc = soc_indices(model, env),
	r! = s.res.r!,
	rm! = s.res.rm!,
	rz! = s.res.rz!,
	rθ! = s.res.rθ!,
	rz = s.rz,
	rθ = s.rθ,
	opts = InteriorPoint113Options(
		# max_iter_inner = 20,
		max_ls = 20,
		r_tol = r_tol,
		κ_tol = κ_tol,
		κ_scale = 1e-1,
		# κ_init = 1e-3,
		solver = :lu_solver,
		diff_sol = false,
		verbose = true,
		))

z0 = ones(length(z))
z_initialize!(z0, model, env, deepcopy(traj.q[t+1]))
interior_point_solve!(ip, z0, θ)
norm(ip.z - deepcopy(traj.z[t]))
norm(ip.z)
norm(traj.z[t])

ip.iterations






################################################################################
# Original IP
################################################################################
s = get_simulation("particle_2D", "flat_2D_nc", "flat_nc")
s.model.μ_world = 0.5
model = s.model
env = s.env

traj = deepcopy(get_trajectory(s.model, s.env,
	joinpath(module_dir(), "src/dynamics/particle_2D/gaits/gait_NC.jld2"),
	load_type = :joint_traj))

plot(hcat(Vector.(traj.z)...)')
plot(hcat(Vector.(traj.θ)...)')

t = 62
z = deepcopy(traj.z[t])
θ = deepcopy(traj.θ[t])
r_tol = 1e-8
κ_tol = 1e-8

ip = interior_point(z, θ,
	ix = linearization_var_index(model, env)[1],
	iy1 = linearization_var_index(model, env)[2],
	iy2 = linearization_var_index(model, env)[3],
	idyn = linearization_term_index(model, env)[1],
	irst = linearization_term_index(model, env)[2],
	ibil = linearization_term_index(model, env)[3],
	idx_ineq = inequality_indices(model, env),
	idx_soc = soc_indices(model, env),
	r! = s.res.r!,
	rm! = s.res.rm!,
	rz! = s.res.rz!,
	rθ! = s.res.rθ!,
	rz = s.rz,
	rθ = s.rθ,
	opts = InteriorPointOptions(
		# max_iter_inner = 20,
		max_ls = 20,
		r_tol = r_tol,
		κ_tol = κ_tol,
		κ_scale = 1e-1,
		# κ_init = 1e-3,
		solver = :lu_solver,
		diff_sol = false,
		verbose = true,
		))

z0 = ones(length(z))
z_initialize!(z0, model, env, deepcopy(traj.q[t+1]))
interior_point_solve!(ip, z0, θ)
norm(ip.z - deepcopy(traj.z[t]))
norm(ip.z)
norm(traj.z[t])

ip.iterations




################################################################################
# Visualize Robots
################################################################################
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
