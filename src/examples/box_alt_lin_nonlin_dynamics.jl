# const ContactControl = Main
# include(joinpath(module_dir(), "src/dynamics/box_alt/visuals.jl"))
# vis = Visualizer()
# open(vis)
#
# # time
# h = 0.01
# T = 175
#
# ################################################################################
# # Nonlinear simulation
# ################################################################################
#
# s = get_simulation("box_alt", "flat_3D_lc", "flat")
# s.model.μ_world = 1.0
#
#
# # initial conditions
# q0 = SVector{s.model.dim.q}([0.0, 0.50, 1.0, 0.1, 0.1, 0.0])
# q1 = SVector{s.model.dim.q}([0.0, 0.55, 1.0, 0.1, 0.1, 0.0])
#
# # simulator
# soft_tol = 1e-8
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
# l = 10
# plot_surface!(vis, env, xlims=[-1,1], ylims=[-l,l])
# plot_lines!(vis, model, sim.traj.q)
# anim = visualize_robot!(vis, model, sim.traj.q)
# visualize_force!(vis, model, env, sim.traj, anim=anim, h=0.0)
#
#
# ################################################################################
# # Reference trajectory
# ################################################################################
#
# s = get_simulation("box_alt", "flat_3D_lc", "flat")
# s.model.μ_world = 1.0
#
#
# # initial conditions
# q0 = SVector{s.model.dim.q}([0.0, 0.00, 0.50, 0.0, 0.0, 0.0])
# q1 = SVector{s.model.dim.q}([0.0, 0.00, 0.50, 0.0, 0.0, 0.0])
#
# # simulator
# tol = 1e-8
# sim = ContactControl.simulator(s, q0, q1, h, T,
# 	ip_opts = ContactControl.InteriorPointOptions(
# 		r_tol = tol, κ_tol = tol,
# 		diff_sol = false,
# 		solver = :lu_solver),
# 	sim_opts = ContactControl.SimulatorOptions(warmstart = false))
#
# # simulate
# @time status = ContactControl.simulate!(sim)
# @test status
# l = 10
# plot_surface!(vis, env, xlims=[-1,1], ylims=[-l,l])
# plot_lines!(vis, model, sim.traj.q)
# anim = visualize_robot!(vis, model, sim.traj.q)
# visualize_force!(vis, model, env, sim.traj, anim=anim, h=0.0)
#
# ref_traj = deepcopy(sim.traj)
#
#
#
#
# ################################################################################
# # Linear simulation
# ################################################################################
#
# s = get_simulation("box_alt", "flat_3D_lc", "flat")
# s.model.μ_world = 1.0
#
# # test reference
# for t = 1:H
# 	r = ContactControl.residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)
# 	@test norm(r) < 1.0e-4
# end
#
# im_traj = ImplicitTraj(ref_traj, s, κ=1e-6)
# ips = im_traj.ip
#
# sim_traj = deepcopy(ref_traj)
#
# q0 = SVector{s.model.dim.q}([0.0, 0.00, 0.50, 0.0, 0.0, 0.0])
# q1 = SVector{s.model.dim.q}([0.0, 0.00, 0.51, 0.0, 0.0, 0.0])
# sim_traj.q[1] = q0
# sim_traj.q[2] = q1
#
# for t = 1:1
# 	@show t
# 	ips[t].opts.verbose = true
# 	# z_initialize!(ips[t].z, model, env, sim_traj.q[t+1])
# 	θ_initialize!(sim_traj.θ[t], s.model, sim_traj.q[t], sim_traj.q[t+1],
# 		sim_traj.u[t], sim_traj.w[t], s.model.μ_world, h)
# 	interior_point_solve!(ips[t], sim_traj.z[t], sim_traj.θ[t])
# 	@show scn.(sim_traj.q[t])
# 	@show scn.(sim_traj.q[t+1])
# 	@show scn.(sim_traj.q[t+2])
# 	@show scn.(ips[t].z[1:s.model.dim.q])
# 	sim_traj.q[t+2] = ips[t].z[1:s.model.dim.q]
# 	@show scn.(sim_traj.q[t+2])
# end
#
#
# #
# # # initial conditions
# # q0 = SVector{s.model.dim.q}([0.0, 0.50, 1.0, 0.1, 0.1, 0.0])
# # q1 = SVector{s.model.dim.q}([0.0, 0.55, 1.0, 0.1, 0.1, 0.0])
# #
# # # simulator
# # soft_tol = 1e-8
# # hard_tol = 1e-8
# # sim = ContactControl.simulator(s, q0, q1, h, T,
# # 	ip_opts = ContactControl.InteriorPointOptions(
# # 		r_tol = soft_tol, κ_tol = soft_tol,
# # 		diff_sol = false,
# # 		solver = :lu_solver),
# # 	sim_opts = ContactControl.SimulatorOptions(warmstart = false))
# #
# # # simulate
# # @time status = ContactControl.simulate!(sim)
# # @test status
# # l = 10
# plot_surface!(vis, env, xlims=[-1,1], ylims=[-l,l])
# plot_lines!(vis, model, sim_traj.q)
# anim = visualize_robot!(vis, model, sim_traj.q)
# visualize_force!(vis, model, env, sim_traj, anim=anim, h=0.0)
#
#
#
#
#
#
#
#
# tt = 10
# sim.traj.q[tt]
# sim.traj.u[tt]
# sim.traj.γ[tt]
#
#
#
#
# s = get_simulation("box_alt", "flat_3D_lc", "flat")
# s.model.μ_world = 1.0
#
# # time
# h = 0.01
# T = 175
#
# # initial conditions
# q0 = SVector{s.model.dim.q}([0.0, 0.50, 1.0, 0.1, 0.1, 0.0])
# q1 = SVector{s.model.dim.q}([0.0, 0.55, 1.0, 0.1, 0.1, 0.0])
#
# q0 = SVector{s.model.dim.q}([0.0, 0.00, 0.0, 0.0, 0.0, 0.0])
# q1 = SVector{s.model.dim.q}([0.0, 0.00, 0.0, 0.0, 0.0, 0.0])
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
#
#
# # reference trajectory
# ref_traj = contact_trajectory(s.model, s.env, H, h)
# ref_traj.h
# qref = [0.0; 0.0; 0.5; 0.0; 0.0; 0.0]
# ur = zeros(s.model.dim.u)
# γr = s.model.g * s.model.m * h * ones(s.model.dim.c)
# br = zeros(s.model.dim.c * friction_dim(s.env))
# ψr = zeros(s.model.dim.c)
# ηr = zeros(s.model.dim.c * friction_dim(s.env))
# wr = zeros(s.model.dim.w)
#
# # set reference
# for t = 1:H
# 	ref_traj.z[t] = pack_z(s.model, s.env, qref, γr, br, ψr, ηr)
# 	ref_traj.θ[t] = pack_θ(s.model, qref, qref, ur, wr, s.model.μ_world, ref_traj.h)
# end
#
# # test reference
# for t = 1:H
# 	r = ContactControl.residual(model, env, ref_traj.z[t], ref_traj.θ[t], 0.0)
# 	@test norm(r) < 1.0e-4
# end
#
#
#
# ref_traj =
# im_traj = ImplicitTraj(ref_traj, s, κ = 1e-4)
#
# # simulate
# @time status = ContactControl.simulate!(sim)
# @test status
# l = 10
# plot_surface!(vis, env, xlims=[-1,1], ylims=[-l,l])
# plot_lines!(vis, model, sim.traj.q)
# anim = visualize_robot!(vis, model, sim.traj)
# visualize_force!(vis, model, env, sim.traj, anim=anim, h=0.0)
#
#
#
# # filename = "box_shadow"
# # MeshCat.convert_frames_to_video(
# #     "/home/simon/Downloads/$filename.tar",
# #     "/home/simon/Documents/$filename.mp4", overwrite=true)
# #
# # convert_video_to_gif(
# #     "/home/simon/Documents/$filename.mp4",
# #     "/home/simon/Documents/$filename.gif", overwrite=true)
