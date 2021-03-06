# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra
# ## Raibert Policy 
include("policy/3D.jl") 

# ## Simulation
s = get_simulation("hopper_3D", "sine2_3D_lc", "sinusoidal")
model_sim = s.model 
env_sim = s.env
nq = model_sim.nq
nu = model_sim.nu
nc = model_sim.nc
nw = model_sim.nw

s_model = get_simulation("hopper_3D", "flat_3D_lc", "flat")

# ## Setup
 
using ContactImplicitMPCH = 92
h = 0.01
N_sample = 5
h_sim = h / N_sample
H_sim = 5000

# ## Policy
v0 = [0.0; 0.2]
Tstance = 0.13 # measure using hop-in-place gait
Tflight = 0.62 # measure using hop-in-place gait
p = raibert_policy(s_model.model, v0=v0, Tstance=Tstance, Tflight=Tflight, h=h);

# ## Initial conditions
off0 = ContactImplicitMPC.SVector{nq}([0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5])
off1 = ContactImplicitMPC.SVector{nq}([v0[1] * h_sim, v0[2] * h_sim, 0.5, 0.0, 0.0, 0.0, 0.5])
q_ref = ContactImplicitMPC.SVector{nq}([0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0])
q0_sim = copy(q_ref) + off0
q1_sim = copy(q_ref) + off1

# ## Simulator
sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
      γ_reg = 0.0,
      undercut = Inf,
      r_tol = 1.0e-8,
      κ_tol = 1.0e-8,),
    sim_opts = SimulatorOptions(warmstart = true)
    );

# ## Simulate
@time status = simulate!(sim);

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
ContactImplicitMPC.plot_surface!(vis, s.env, n=200, xlims = [-1, 40]);
visualize_robot!(vis, model_sim, sim.traj, sample=20);
