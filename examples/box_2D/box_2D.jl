

# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC
using LinearAlgebra

# ## Simulation
s = get_simulation("box_2D", "flat_2D_lc", "flat_lc")
model = s.model
env = s.env
h_sim = 0.05
H_sim = 200

# ## Initial conditions
q1_sim = ContactImplicitMPC.SVector{model.dim.q}(8*h_sim, 1.5, 2*h_sim);
q0_sim = ContactImplicitMPC.SVector{model.dim.q}(0, 1.5, 0);

# ## Hard contact: relaxation = 1e-8
κ_relax = 1.0e-8
sim_hard = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    ip_opts = InteriorPointOptions(
		γ_reg = 0.0,
		undercut = 4.0,
        r_tol = 1.0e-8,
        κ_tol = κ_relax,),
    sim_opts = SimulatorOptions(warmstart = true));

# ### Simulate
status = simulate!(sim_hard, verbose = false);

# ### Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ### Visualize
ContactImplicitMPC.plot_surface!(vis, env, ylims = (-1,1.), xlims = (-1,13.))
ContactImplicitMPC.plot_lines!(vis, model, sim_hard.traj.q)
anim = visualize_robot!(vis, model, sim_hard.traj);


# ## Soft contact: relaxation = 2e-2
κ_relax = 4.0e-2
sim_soft = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    ip_opts = InteriorPointOptions(
		γ_reg = 0.0,
		undercut = 4.0,
        r_tol = 1.0e-8,
        κ_tol = κ_relax,),
    sim_opts = SimulatorOptions(warmstart = true));

# ### Simulate
status = simulate!(sim_soft, verbose = false);

### Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ### Visualize
ContactImplicitMPC.plot_surface!(vis, env, ylims = (-1,1.), xlims = (-1,13.))
ContactImplicitMPC.plot_lines!(vis, model, sim_soft.traj.q)
anim = visualize_robot!(vis, model, sim_soft.traj);


# ## Slightly soft contact: relaxation = 1e-4
κ_relax = 1.0e-4
sim_mpc = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    ip_opts = InteriorPointOptions(
		γ_reg = 0.0,
		undercut = 4.0,
        r_tol = 1.0e-8,
        κ_tol = κ_relax,),
    sim_opts = SimulatorOptions(warmstart = true));

# ### Simulate
status = simulate!(sim_mpc, verbose = false);

# ### Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ### Visualize
ContactImplicitMPC.plot_surface!(vis, env, ylims = (-1,1.), xlims = (-1,13.))
ContactImplicitMPC.plot_lines!(vis, model, sim_mpc.traj.q)
anim = visualize_robot!(vis, model, sim_mpc.traj);
