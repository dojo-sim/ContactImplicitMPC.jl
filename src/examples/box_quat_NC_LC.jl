const ContactControl = Main
include(joinpath(module_dir(), "src/dynamics/box/visuals.jl"))
vis = Visualizer()
open(vis)

################################################################################
# Parameters
################################################################################
# time
h = 0.01
T = 80

# initial conditions
r0 = [-1.5; 0.0; 1.5]
v0 = [8.0; 0.0; 0.0]
quat0 = [1.0; 0.0; 0.0; 0.0]
_quat0 = UnitQuaternion(RotX(0.0 * π))
quat0 = [_quat0.w; _quat0.x; _quat0.y; _quat0.z]
ω0 = [0.8; 0.5; 0.2]

# @assert norm(q0[4:7]) ≈ 1.0
# @assert norm(q1[4:7]) ≈ 1.0

# Simulation tolerance
tol = 1e-8

################################################################################
# Linearized Cone
################################################################################
s = get_simulation("box", "flat_3D_lc", "flat_lc")
s.model.μ_world = 1.0
# Rn + quaternion space
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
	collect([(1:3)..., (8:num_var(s.model, s.env))...]),
	collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
	[collect((4:7))],
	[collect((4:6))])
q0 = SVector{s.model.dim.q}([r0; quat0])
q1 = SVector{s.model.dim.q}([r0 + v0 * h; 0.5 * h * L_multiply(quat0) * [sqrt((2.0 / h)^2.0 - ω0' * ω0); ω0]])

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	space = rq_space,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = tol, κ_tol = tol,
		diff_sol = false,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)
@test all([norm(q[4:7]) ≈ 1.0 for q in sim.traj.q])

plot_surface!(vis, s.env, xlims=[-10,10], ylims=[-10,10])
plot_lines!(vis, s.model, sim.traj.q, name = :LClines)
anim = visualize_robot!(vis, s.model, sim.traj, name = :LC)
visualize_force!(vis, s.model, s.env, sim.traj, anim=anim, shift=-0.25, name = :LC)

plt = plot(legend=false)
plot!(plt, hcat(Vector.([q[1:end] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b])...)', color=:black, linewidth=4.0)

# ################################################################################
# # Save Trajectory
# ################################################################################
# traj = deepcopy(sim.traj)
# gait_path = joinpath(@__DIR__, "..", "dynamics", "box", "gaits", "gait_LC.jld2")
# @save gait_path traj
#
# # Reload trajectory
# res = JLD2.jldopen(gait_path)
# loaded_traj = res["traj"]
#
# traj = deepcopy(get_trajectory(s.model, s.env,
# 	joinpath(module_dir(), "src/dynamics/box/gaits/gait_LC.jld2"),
# 	load_type = :joint_traj))
# plot(hcat(Vector.(traj.q)...)')


################################################################################
# Non Linear Cone
################################################################################
s = get_simulation("box", "flat_3D_nc", "flat_nc")
s.model.μ_world = 1.0
# Rn + quaternion space
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
	collect([(1:3)..., (8:num_var(s.model, s.env))...]),
	collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
	[collect((4:7))],
	[collect((4:6))])
q0 = SVector{s.model.dim.q}([r0; quat0])
q1 = SVector{s.model.dim.q}([r0 + v0 * h; 0.5 * h * L_multiply(quat0) * [sqrt((2.0 / h)^2.0 - ω0' * ω0); ω0]])

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	space = rq_space,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = tol, κ_tol = tol,
		diff_sol = false,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot_surface!(vis, s.env, xlims=[-10,10], ylims=[-10,10])
plot_lines!(vis, s.model, sim.traj.q, name = :NClines)
visualize_robot!(vis, s.model, sim.traj, name = :NC, anim = anim)
visualize_force!(vis, s.model, s.env, sim.traj, anim=anim, shift=-0.25, name = :NC)

plt = plot(legend=false)
plot!(plt, hcat(Vector.([q[1:end] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b])...)', color=:black, linewidth=4.0)

# ################################################################################
# # Save Trajectory
# ################################################################################
# traj = deepcopy(sim.traj)
# gait_path = joinpath(@__DIR__, "..", "dynamics", "box", "gaits", "gait_NC.jld2")
# @save gait_path traj
#
# # Reload trajectory
# res = JLD2.jldopen(gait_path)
# loaded_traj = res["traj"]
#
# traj = deepcopy(get_trajectory(s.model, s.env,
# 	joinpath(module_dir(), "src/dynamics/box/gaits/gait_NC.jld2"),
# 	load_type = :joint_traj))
# plot(hcat(Vector.(traj.q)...)')

################################################################################
# Linearized Cone Latest
################################################################################
s = get_simulation("box", "flat_3D_lc", "flat_lc")
s.model.μ_world = 1.0

# Rn + quaternion space
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
	collect([(1:3)..., (8:num_var(s.model, s.env))...]),
	collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
	[collect((4:7))],
	[collect((4:6))])
q0 = SVector{s.model.dim.q}([r0; quat0])
q1 = SVector{s.model.dim.q}([r0 + v0 * h; 0.5 * h * L_multiply(quat0) * [sqrt((2.0 / h)^2.0 - ω0' * ω0); ω0]])

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	space = rq_space,
	ip_type = :interior_point_latest,
	ip_opts = ContactControl.InteriorPoint116Options(
		r_tol = tol, κ_tol = tol,
		diff_sol = false,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)
@test all([norm(q[4:7]) ≈ 1.0 for q in sim.traj.q])

plot_surface!(vis, s.env, xlims=[-10,10], ylims=[-10,10])
plot_lines!(vis, s.model, sim.traj.q, name = :LClines)
anim = visualize_robot!(vis, s.model, sim.traj, name = :LC)
visualize_force!(vis, s.model, s.env, sim.traj, anim=anim, shift=-0.25, name = :LC)

plt = plot(legend=false)
plot!(plt, hcat(Vector.([q[1:end] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b])...)', color=:black, linewidth=4.0)



################################################################################
# Non Linear Cone Latest
################################################################################
s = get_simulation("box", "flat_3D_nc", "flat_nc")
s.model.μ_world = 1.0

# Rn + quaternion space
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
	collect([(1:3)..., (8:num_var(s.model, s.env))...]),
	collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
	[collect((4:7))],
	[collect((4:6))])
q0 = SVector{s.model.dim.q}([r0; quat0])
q1 = SVector{s.model.dim.q}([r0 + v0 * h; 0.5 * h * L_multiply(quat0) * [sqrt((2.0 / h)^2.0 - ω0' * ω0); ω0]])

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	space = rq_space,
	ip_type = :interior_point_latest,
	ip_opts = ContactControl.InteriorPoint116Options(
		r_tol = tol, κ_tol = tol,
		diff_sol = false,
		max_iter = 100,
		max_ls = 3,
		solver = :lu_solver),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)
@test all([norm(q[4:7]) ≈ 1.0 for q in sim.traj.q])

plot_surface!(vis, s.env, xlims=[-10,10], ylims=[-10,10])
plot_lines!(vis, s.model, sim.traj.q, name = :LClines)
anim = visualize_robot!(vis, s.model, sim.traj, name = :LC)
visualize_force!(vis, s.model, s.env, sim.traj, anim=anim, shift=-0.25, name = :LC)

plt = plot(legend=false)
plot!(plt, hcat(Vector.([q[1:end] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)
plot!(plt, hcat(Vector.([b[1:end] for b in sim.traj.b])...)', color=:black, linewidth=4.0)






################################################################################
# Dev
################################################################################
s = get_simulation("box", "flat_3D_nc", "flat_nc")
s.model.μ_world = 1.0
model = s.model
env = s.env
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
	collect([(1:3)..., (8:num_var(s.model, s.env))...]),
	collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
	[collect((4:7))],
	[collect((4:6))])

traj = deepcopy(get_trajectory(s.model, s.env,
	joinpath(module_dir(), "src/dynamics/box/gaits/gait_NC.jld2"),
	load_type = :joint_traj))

r_tol = 1e-8
κ_tol = 1e-8

iter = 0
for t = 1:T
	z = deepcopy(traj.z[t])
	θ = deepcopy(traj.θ[t])

	ip = interior_point_latest(z, θ,
		s = rq_space,
		oss = OptimizationSpace13(model, env),
		r! = s.res.r!,
		rz! = s.res.rz!,
		rθ! = s.res.rθ!,
		rz = s.rz,
		rθ = s.rθ,
		opts = InteriorPoint116Options(
			max_iter = 20,
			max_ls = 3,
			r_tol = r_tol,
			κ_tol = κ_tol,
			solver = :lu_solver,
			diff_sol = false,
			# verbose = true,
			))

	z0 = ones(length(z))
	z_initialize!(z0, model, env, deepcopy(traj.q[t+1]))
	interior_point_solve!(ip, z0, θ)
	norm(ip.z - deepcopy(traj.z[t]))

	iter += ip.iterations
	end
mean_iter = iter / T





################################################################################
# Original IP
################################################################################
s = get_simulation("box", "flat_3D_nc", "flat_nc")
s.model.μ_world = 1.0
model = s.model
env = s.env
rq_space = rn_quaternion_space(num_var(s.model, s.env) - 1, x -> Gz_func(s.model, s.env, x),
	collect([(1:3)..., (8:num_var(s.model, s.env))...]),
	collect([(1:3)..., (7:num_var(s.model, s.env)-1)...]),
	[collect((4:7))],
	[collect((4:6))])

traj = deepcopy(get_trajectory(s.model, s.env,
	joinpath(module_dir(), "src/dynamics/box/gaits/gait_NC.jld2"),
	load_type = :joint_traj))

r_tol = 1e-8
κ_tol = 1e-8

iter = 0
for t = 1:T
	z = deepcopy(traj.z[t])
	θ = deepcopy(traj.θ[t])

	ip = interior_point(z, θ,
		s = rq_space,
		oss = OptimizationSpace13(model, env),
		r! = s.res.r!,
		rz! = s.res.rz!,
		rθ! = s.res.rθ!,
		rz = s.rz,
		rθ = s.rθ,
		opts = InteriorPointOptions(
			# max_iter = 20,
			# max_ls = 3,
			r_tol = r_tol,
			κ_tol = κ_tol,
			κ_scale = 1e-1,
			# κ_init = 1e-3,
			solver = :lu_solver,
			diff_sol = false,
			# verbose = true,
			))

	z0 = ones(length(z))
	z_initialize!(z0, model, env, deepcopy(traj.q[t+1]))
	interior_point_solve!(ip, z0, θ)
	norm(ip.z - deepcopy(traj.z[t]))

	iter += ip.iterations
end
mean_iter = iter / T
36.0 / 12.6




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

# filename = "box_quat_soc"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
