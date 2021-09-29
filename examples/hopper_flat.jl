using ContactImplicitMPC 
using LinearAlgebra 
using StaticArrays

# Visualizer
vis = ContactImplicitMPC.Visualizer()
open(vis)

# Simulation
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model
env = s.env

# Reference Trajectory
ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
    load_type = :joint_traj))

H = ref_traj.H
h = ref_traj.h

# MPC setup 
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 1000# 10*H*N_sample #500*H*N_sample

# barrier parameter
κ_mpc = 2.0e-4

obj = TrackingObjective(model, env, H_mpc,
    q = [Diagonal(1.0e-1 * [0.1,3,1,3])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-3, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
	ip_opts = InteriorPointOptions(
					undercut = 5.0,
					κ_tol = κ_mpc,
					r_tol = 1.0e-8,
					diff_sol = true,
					solver = :empty_solver,
					max_time = 1e5,),
    n_opts = NewtonOptions(
        r_tol = 3e-4,
		# verbose = true,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        # altitude_impact_threshold = 0.05,
        # altitude_verbose = true,
        )
    )

# Initial conditions
q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))

# Simulator
sim = ContactImplicitMPC.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactImplicitMPC.InteriorPointOptions(
		γ_reg = 0.0,
		undercut = Inf,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8,),
    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true))


# Simulate
status = ContactImplicitMPC.simulate!(sim, verbose = true)

# Visualize 
anim = visualize_robot!(vis, model, sim.traj, sample=5)
anim = visualize_force!(vis, model, env, sim.traj, anim=anim, h=h_sim, sample = 5)

################################################################################
# Timing result
################################################################################
process!(sim)
# Time budget
ref_traj.h
# Time used on average
sim.stats.μ_dt
# Speed ratio
H_sim * h_sim / sum(sim.stats.dt)


