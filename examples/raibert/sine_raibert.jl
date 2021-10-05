# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC 
using LinearAlgebra 
using StaticArrays

# ## Raibert Policy 
include("policy/2D.jl")

# ## Simulation
s_sim = get_simulation("hopper_2D", "sine2_2D_lc", "sinusoidal")
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model
env = s.env

# ## Reference Trajectory
ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
    load_type = :joint_traj))
H = ref_traj.H
h = ref_traj.h
N_sample = 5
h_sim = h / N_sample
H_sim = 1000# 100*H*N_sample

# ## Raibert policy
v0 = 0.2
Tstance = 0.13 # measure using hop-in-place gait
Tflight = 0.62 # measure using hop-in-place gait
p = raibert_policy(s_sim.model, v0=v0, Tstance=Tstance, Tflight=Tflight, h=h)

# ## Initial conditions
off0 = SVector{model.dim.q}([0.0, 0.5, 0.0, 0.0])
off1 = SVector{model.dim.q}([0*v0*h_sim, 0.5, 0.0, 0.0])
q_ref = SVector{model.dim.q}([0.0, 0.5, 0.0, 0.5])
q0_sim = copy(q_ref) + off0
q1_sim = copy(q_ref) + off1

# ## Simulator
sim = ContactImplicitMPC.simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactImplicitMPC.InteriorPointOptions(
		γ_reg = 0.0,
		undercut = Inf,
		r_tol = 1.0e-8,
		κ_tol = 1.0e-8,),
    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true));

# ## Simulate
status = ContactImplicitMPC.simulate!(sim, verbose = true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
ContactImplicitMPC.plot_surface!(vis, s.env, n=200, xlims = [-1, 40])
anim = visualize_robot!(vis, model, sim.traj, sample=5)