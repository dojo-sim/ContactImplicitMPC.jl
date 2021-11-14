# PREAMBLE

# PKG_SETUP

# ## Setup
 
using ContactImplicitMPC
using LinearAlgebra

# Define a special stride where x and z are updated.
function get_stride(model, traj)
    stride = zeros(SizedVector{model.nq})
    stride[1:2] = traj.q[end-1][1:2] - traj.q[1][1:2]
    return stride
end

# ## Simulation
s = get_simulation("hopper_2D", "flat_2D_lc", "flat");
s_sim = get_simulation("hopper_2D", "stairs3_2D_lc", "stairs");

# ## Reference Trajectory (stairs) 
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/hopper_stair_ref.jld2"),
    load_type=:split_traj_alt));

H = ref_traj.H
h = ref_traj.h
h_sim = h / N_sample
H_sim = 240*N_sample

# ## MPC parameters (stairs)
N_sample = 10
H_mpc = 10
κ_mpc = 2.0e-4
n_opts = NewtonOptions(
    r_tol = 3e-4,
    max_iter = 5)
mpc_opts = CIMPCOptions(
    altitude_update = true,
    altitude_impact_threshold = 0.1,
    altitude_verbose = true)

obj = TrackingVelocityObjective(s.model, s.env, H_mpc,
    v = [Diagonal(1e-3 * [1e-2,1,1,10]) for t = 1:H_mpc],
    q = [[Diagonal(1e-0 * [1e1,1e-1,1,1])   for t = 1:H_mpc-5]; [Diagonal(1e-1 * [1,1e-1,1e1,0.1])   for t = 1:5]],
    u = [Diagonal(1e-0 * [1e0, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1e-100 * ones(s.model.nc)) for t = 1:H_mpc],
    b = [Diagonal(1e-100 * ones(s.model.nc * friction_dim(s.env))) for t = 1:H_mpc])

p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = n_opts,
    mpc_opts = mpc_opts);

# ## Initial conditions (stairs)
q1_sim = ContactImplicitMPC.SVector{model.nq}(copy(ref_traj.q[2]))
q0_sim = ContactImplicitMPC.SVector{model.nq}(copy(q1_sim - (copy(ref_traj.q[2]) - copy(ref_traj.q[1])) / N_sample));

# ## Simulator (stairs)
sim_stair = ContactImplicitMPC.simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactImplicitMPC.InteriorPointOptions(
        γ_reg = 0.0,
        undercut = Inf,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8),
    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true));

# ## Simulate (stairs)
@time status = simulate!(sim_stair)

# ## Visualizer
vis = ContactImplicitMPC.Visualizer()
ContactImplicitMPC.render(vis)

# ## Visualize (stairs)
anim = visualize_robot!(vis, model, sim_stair.traj, sample=10, name=:Sim, α=1.0)
stairs!(vis)

# ## Reference trajectory (flip)
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/hopper_tall_flip_ref.jld2"),
    load_type=:split_traj_alt));

for t = 1:ref_traj.H+2
    ref_traj.q[t][1] += sim_stair.traj.q[end][1] # shift flip to top of stairs
end

H = ref_traj.H
h = ref_traj.h
h_sim = h / N_sample
H_sim = 64*N_sample

# ## MPC setup (flip)
obj = TrackingVelocityObjective(s.model, s.env, H_mpc,
    v = [Diagonal(1e-13 * [1e-2,1,1,10]) for t = 1:H_mpc],
    q = [[Diagonal(1e-10 * [1e1,1e1,1,1])   for t = 1:H_mpc-5]; [Diagonal(1e-11 * [1,1e1,1e1,0.1])   for t = 1:5]],
    u = [Diagonal(1e-0 * [1e0, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1e-100 * ones(s.model.nc)) for t = 1:H_mpc],
    b = [Diagonal(1e-100 * ones(s.model.nc * friction_dim(s.env))) for t = 1:H_mpc])

p = ci_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = n_opts,
    mpc_opts = mpc_opts,
    );
 
# ## Set initial configurations to simulation result
q0_sim = deepcopy(ContactImplicitMPC.SVector{model.nq}(sim_stair.traj.q[end-1]))
q1_sim = deepcopy(ContactImplicitMPC.SVector{model.nq}(sim_stair.traj.q[end]))

# ## Simulator (flip)
sim_flip = ContactImplicitMPC.simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactImplicitMPC.InteriorPointOptions(
        γ_reg = 0.0,
        undercut = Inf,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8),
    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true));

# ## Simulate (flip)
@time status = ContactImplicitMPC.simulate!(sim_flip)

# ## Visualize (flip)
anim = visualize_robot!(vis, model, sim_flip.traj, sample=10, name=:Sim, α=1.0)
anim = visualize_robot!(vis, model, ref_traj, anim=anim, name=:Ref, α=0.3)

# ## Visualize (parkour)
ref_traj_full = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/hopper_stairs_flip_ref.jld2"),
    load_type=:split_traj_alt));

N_sample = 10
sim_traj_full = [sim_stair.traj.q[1:end-2]; sim_flip.traj.q]

anim = visualize_robot!(vis, model,
    [sim_traj_full[1:N_sample:end]..., [sim_traj_full[end] for i = 1:50]...],
    name = :Sim, α = 1.0, h = h_sim*N_sample);

anim = visualize_robot!(vis, model,
    [ref_traj_full.q..., [ref_traj_full.q[end] for i = 1:50]...],
    anim = anim, name = :Ref, α = 0.15, h = h_sim*N_sample);

ContactImplicitMPC.stairs!(vis)

