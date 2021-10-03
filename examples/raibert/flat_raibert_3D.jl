# PREAMBLE

# PKG_SETUP

# ## Setup

using ContactImplicitMPC 
using LinearAlgebra 
using StaticArrays

# ## Simulation
model_sim = get_model("hopper_3D", surf="flat");
nq = model_sim.dim.q
nu = model_sim.dim.u
nc = model_sim.dim.c
nb = model_sim.dim.b
nw = model_sim.dim.w

# ## Setup
H = 92
h = 0.01
N_sample = 5
h_sim = h / N_sample
H_sim = 5000

# ## Raibert policy
v0 = [0.0; 0.2]
Tstance = 0.13 # measure using hop-in-place gait
Tflight = 0.62 # measure using hop-in-place gait
include(joinpath(module_dir(), "src/controller/raibert_3D_policy.jl"));
p = raibert_policy(model_sim, v0=v0, Tstance=Tstance, Tflight=Tflight, h=h);

# ## Initial conditions
off0 = SVector{nq,T}([0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5])
off1 = SVector{nq,T}([v0[1]*h_sim, v0[2]*h_sim, 0.5, 0.0, 0.0, 0.0, 0.5])
q_ref = SVector{nq,T}([0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0])
q0_sim = copy(q_ref) + off0
q1_sim = copy(q_ref) + off1

# ## Disturbances
w = [zeros(nw) for t=1:Int(ceil(H_sim/N_sample))]
d = open_loop_disturbances(w);

# ## Simulator
sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    d = d,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    );

# ## Simulate
@time status = simulate!(sim);

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
open(vis)

# ## Visualize
anim = visualize_robot!(vis, model_sim, sim.traj, sample=20)
