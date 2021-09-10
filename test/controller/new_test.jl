


T = Float64
# vis = Visualizer()
# open(vis)
include(joinpath(module_dir(), "src", "dynamics", "quadruped", "visuals.jl"))

# s = get_simulation("quadruped", "flat_2D_lc", "flat")

ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))

im_traj = ImplicitTraj(ref_traj, s, ip_type=:mehrotra)

H = 10
H = ref_traj.H
h = ref_traj.h
q0_sim = SVector{s.model.dim.q}(deepcopy(ref_traj.q[1]))
q1_sim = SVector{s.model.dim.q}(deepcopy(ref_traj.q[2]))

p = open_loop_policy(ref_traj.u; N_sample = 1)

sim_b = simulator(s, q0_sim, q1_sim, h, H, p=p, ip_type=:interior_point,
    sim_opts = SimulatorOptions(warmstart=true))
sim_m = simulator(s, q0_sim, q1_sim, h, H, p=p, ip_type=:mehrotra,
    ip_opts = MehrotraOptions(max_iter=100),
    sim_opts = SimulatorOptions(warmstart=false))
