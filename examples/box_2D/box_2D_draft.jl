#
#
# # PREAMBLE
#
# # PKG_SETUP
#
# # ## Setup
#
# using ContactImplicitMPC
# using LinearAlgebra
#
# # ## Simulation
# s = get_simulation("box_2D", "flat_2D_lc", "flat_lc")
# model = s.model
# env = s.env
# h_sim = 0.05
# H_sim = 200
#
# # ## Initial conditions
# q1_sim = ContactImplicitMPC.SVector{model.dim.q}(8*h_sim, 1.5, 2*h_sim);
# q0_sim = ContactImplicitMPC.SVector{model.dim.q}(0, 1.5, 0);
#
# # ## Hard contact: relaxation = 1e-8
# κ_relax = 1.0e-8
# sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
#     ip_opts = InteriorPointOptions(
# 		γ_reg = 0.0,
# 		undercut = 4.0,
#         r_tol = 1.0e-8,
#         κ_tol = κ_relax,),
#     sim_opts = SimulatorOptions(warmstart = true));
#
# # ### Simulate
# status = simulate!(sim, verbose = true);
#
# # ### Visualizer
# vis = ContactImplicitMPC.Visualizer()
# ContactImplicitMPC.render(vis)
#
# # ### Visualize
# plot_surface!(vis, env)
# plot_lines!(vis, model, sim.traj.q)
# anim = visualize_robot!(vis, model, sim.traj);
#
#
# # ## Soft contact: relaxation = 2e-2
# κ_relax = 4.0e-2
# sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
#     ip_opts = InteriorPointOptions(
# 		γ_reg = 0.0,
# 		undercut = 4.0,
#         r_tol = 1.0e-8,
#         κ_tol = κ_relax,),
#     sim_opts = SimulatorOptions(warmstart = true));
#
# # ### Simulate
# status = simulate!(sim, verbose = true);
#
# # ### Visualizer
# vis = ContactImplicitMPC.Visualizer()
# ContactImplicitMPC.render(vis)
#
# # ### Visualize
# plot_surface!(vis, env)
# plot_lines!(vis, model, sim.traj.q)
# anim = visualize_robot!(vis, model, sim.traj);
#
#
# # ## Slightly soft contact: relaxation = 1e-4
# κ_relax = 1.0e-4
# sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
#     ip_opts = InteriorPointOptions(
# 		γ_reg = 0.0,
# 		undercut = 4.0,
#         r_tol = 1.0e-8,
#         κ_tol = κ_relax,),
#     sim_opts = SimulatorOptions(warmstart = true));
#
# # ### Simulate
# status = simulate!(sim, verbose = true);
#
# # ### Visualizer
# vis = ContactImplicitMPC.Visualizer()
# ContactImplicitMPC.render(vis)
#
# # ### Visualize
# plot_surface!(vis, env)
# plot_lines!(vis, model, sim.traj.q)
# anim = visualize_robot!(vis, model, sim.traj);
#
#
#
#
#
# ################################################################################
# # Box 2D
# ################################################################################
# dir = joinpath(@__DIR__, "box_2D")
# dir = joinpath(@__DIR__, "..", "..", "src", "dynamics", "box_2D")
# include(joinpath(dir, "model.jl"))
# model = deepcopy(box_2D)
#
# path_base = joinpath(dir, "dynamics/base.jld2")
# path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
#
# expr_base = generate_base_expressions(model)
# save_expressions(expr_base, path_base, overwrite=true)
# instantiate_base!(model, path_base)
#
# expr_dyn = generate_dynamics_expressions(model)
# save_expressions(expr_dyn, path_dyn, overwrite=true)
# instantiate_dynamics!(model, path_dyn)
#
# ################################################################################
# # Box 2D (flat + LC)
# ################################################################################
# dir_model = joinpath(module_dir(), "src/dynamics/box_2D")
# dir_sim   = joinpath(module_dir(), "src/simulation/box_2D")
# model = deepcopy(box_2D)
# env = deepcopy(flat_2D_lc)
# sim = Simulation(model, env)
#
# path_base = joinpath(dir_model, "dynamics/base.jld2")
# path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
# path_res = joinpath(dir_sim, "flat_lc/residual.jld2")
# path_jac = joinpath(dir_sim, "flat_lc/jacobians.jld2")
#
# instantiate_base!(sim.model, path_base)
# instantiate_dynamics!(sim.model, path_dyn)
#
# expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
#
# save_expressions(expr_res, path_res, overwrite=true)
# isfile(path_jac) && rm(path_jac)
# ContactImplicitMPC.JLD2.save(path_jac, "rz_sp", rz_sp, "rθ_sp", rθ_sp)
# instantiate_residual!(sim, path_res, path_jac)
#
#
#
#
#
#
#
#
#
#
#
#
# using Pkg
# using LinearAlgebra
# Pkg.activate(joinpath(@__DIR__, "..", ".."))
# const ContactImplicitMPC = Main
#
#
# vis = Visualizer()
# open(vis)
#
#
# # ## Simulation
# s = get_simulation("box_2D", "flat_2D_lc", "flat_lc")
# model = s.model
# env = s.env
#
# h_sim = 0.01
# H_sim = 300
#
# # ## Initial conditions
# q1_sim = ContactImplicitMPC.SVector{model.dim.q}(5*h_sim, 1.5, 2*h_sim);
# q0_sim = ContactImplicitMPC.SVector{model.dim.q}(0, 1.5, 0.00);
#
# # ## Simulator
# sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
#     ip_opts = InteriorPointOptions(
# 		γ_reg = 0.0,
# 		undercut = 5.0,
#         r_tol = 1.0e-8,
#         κ_tol = 1.0e-8,),
#     sim_opts = SimulatorOptions(warmstart = true));
#
# # ## Simulate
# status = simulate!(sim, verbose = true);
#
#
# # ## Visualize
# plot_surface!(vis, env, n=200)
# plot_lines!(vis, model, sim.traj.q)
# visualize_robot!(vis, model, sim.traj)
#
# # impact of central path parameter on simulation
# hard contact: what we use for sim
# soft contact
# slightly softened contact: what we use for planning
#
# # impact of central path parameter on gradients




# # ## Implicit Function Theorem
# function get_gradient(κ_relax, model, env)
# 	q1_sim = ContactImplicitMPC.SVector{model.dim.q}(0, 1.0, 0);
# 	q0_sim = ContactImplicitMPC.SVector{model.dim.q}(0, 1.0, 0);
#
# 	# ## Slightly soft contact: relaxation = 1e-4
# 	sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
# 	    ip_opts = InteriorPointOptions(
# 			γ_reg = 0.0,
# 			undercut = 4.,
# 	        r_tol = 1.0e-8,
# 	        κ_tol = κ_relax,),
# 	    sim_opts = SimulatorOptions(warmstart = true));
#
# 	# ### Simulate
# 	status = simulate!(sim, verbose = false);
#
# 	nz = num_var(model, env)
# 	nθ = num_data(model)
# 	rz = zeros(nz, nz)
# 	rθ = zeros(nz, nθ)
#
# 	z = sim.traj.z[end]
# 	θ = sim.traj.θ[end]
# 	rz!(sim.ip, rz, z, θ)
# 	s.res.rθ!(rθ, z, θ)
# 	jac = -rz \ rθ
#
# 	iq2 = index_q2(model, env)[1] # state x
# 	iu1 = index_u1(model)[1] # control x
# 	∂x∂ux = jac[iq2, iu1]
# 	return ∂x∂ux
# end
#
# # ## Contact smoothness -> gradient information
# κs = [exp(log(10)* (-8 + 0.1i)) for i = 1:60]
# grads = [get_gradient(exp(log(10)* (-8 + 0.1i)), model, env) for i = 1:60]
# plot(κs, grads, yaxis = :log, xaxis = :log)
# plot(κs, grads ./ κs)
#
# grad_mpc = get_gradient(sim_mpc, model, env)
# grad_soft = get_gradient(sim_soft, model, env)
#
# ### Visualizer
# vis = ContactImplicitMPC.Visualizer()
# ContactImplicitMPC.render(vis)
#
# # ### Visualize
# ContactImplicitMPC.plot_surface!(vis, env, ylims = (-1,1.), xlims = (-1,13.))
# anim = visualize_robot!(vis, model, sim_mpc.traj);
#
# function visualize_gradient
