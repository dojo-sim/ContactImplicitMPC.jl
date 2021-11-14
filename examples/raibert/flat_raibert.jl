# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# ## Raibert Policy 
include("policy/2D.jl")

# ## Simulation
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model
env = s.env

# ## Reference Trajectory
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
    load_type = :joint_traj))
H = ref_traj.H
h = ref_traj.h
N_sample = 5
h_sim = h / N_sample
H_sim = 1000# 100*H*N_sample 

# ## Policy
v0 = 0.2
Tstance = 0.13 # measure using hop-in-place gait
Tflight = 0.62 # measure using hop-in-place gait
p = raibert_policy(s.model, v0=v0, Tstance=Tstance, Tflight=Tflight, h=h);

# ## Initial conditions
off0 = ContactImplicitMPC.SVector{model.nq}([0.0, 0.5, 0.0, 0.0])
off1 = ContactImplicitMPC.SVector{model.nq}([0.0 * v0 * h_sim, 0.5, 0.0, 0.0])
q_ref = ContactImplicitMPC.SVector{model.nq}([0.0, 0.5, 0.0, 0.5])
q0_sim = copy(q_ref) + off0
q1_sim = copy(q_ref) + off1

# ## Simulator
sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
      γ_reg = 0.0,
      undercut = Inf,
      r_tol = 1.0e-8,
      κ_tol = 1.0e-8,),
    sim_opts = SimulatorOptions(warmstart = true));

# ## Simulation
@time status = simulate!(sim, verbose = true)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize
anim = visualize_robot!(vis, model, sim.traj, sample=5);
