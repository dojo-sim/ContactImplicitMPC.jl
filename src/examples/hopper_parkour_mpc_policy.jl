# const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# Define a special stride where x and z are updated.
function get_stride(model::Hopper2D, traj::ContactTraj)
    stride = zeros(SizedVector{model.dim.q})
    stride[1:2] = traj.q[end-1][1:2] - traj.q[1][1:2]
    return stride
end

# get hopper model
model_sim = get_model("hopper_2D", surf="stairs")
model = get_model("hopper_2D")

# MPC parameters
N_sample = 10
H_mpc = 10
κ_mpc = 1.0e-4
n_opts = NewtonOptions(
    r_tol = 3e-4,
    max_iter = 5)
mpc_opts = LinearizedMPCOptions(
    altitude_update = true,
    altitude_impact_threshold = 0.1,
    altitude_verbose = true)


################################################################################
# Stair climbing
################################################################################

# get stair trajectory
ref_traj_ = get_trajectory(model,
    joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "parkour", "hopper_stair1.jld2"),
    load_type=:split_traj_alt)
ref_traj = deepcopy(ref_traj_)

# Horizon
H = ref_traj.H
h = ref_traj.h
h_sim = h / N_sample
H_sim = 254*N_sample

# Objective
obj = TrackingVelocityObjective(H_mpc, model.dim,
    v = [Diagonal(1e-3 * [0.1,1,1,10]) for t = 1:H_mpc],
    q = [[Diagonal(1.0e-0 * [0.3,0.3,0.3,1])   for t = 1:H_mpc]; [Diagonal(1.0e+2 * [5,0.3,0.3,0.1])   for t = 1:H_mpc]],
    u = [Diagonal(1.0e-1 * [1e0, 1e-0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-3 * ones(model.dim.b)) for t = 1:H_mpc])

# Policy
p = linearized_mpc_policy(ref_traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = n_opts,
    mpc_opts = mpc_opts,
    )

# Get initial configurations
q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

# Simulate the stair climbing
sim_stair = ContactControl.simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim_stair)


# plot_surface!(vis, model_sim.env, n=400)
anim = visualize_robot!(vis, model, sim_stair.traj, sample=10, name=:Sim, α=1.0)

# ghost
ref_traj_full = get_trajectory(model,
    joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "parkour", "hopper_stairs_3_flip_v3.jld2"),
    load_type=:split_traj_alt)
ref_traj_full = deepcopy(ref_traj_full)

idx = [1, 300, 500, 1100, 1300, 1900, 2100]
α = [0.2, 0.2, 0.4, 0.4, 0.6, 0.6, 0.6]
hopper_parkour_ghost!(vis, sim_stair, sim_stair.traj, ref_traj_full, idx = idx, α = α)

################################################################################
# Front flip
################################################################################

# get trajectory
ref_traj_ = get_trajectory(model,
    joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "parkour", "hopper_tall_flip5.jld2"),
    load_type=:split_traj_alt)
ref_traj = deepcopy(ref_traj_)
# offset the trajectory
for t = 1:ref_traj.H+2
    ref_traj.q[t][1] += sim_stair.traj.q[end][1]
end

# time
H = ref_traj.H
h = ref_traj.h
h_sim = h / N_sample
H_sim = 64*N_sample

# Objective
obj = TrackingVelocityObjective(H_mpc, model.dim,
    v = [Diagonal(1e-3 * [0.1,1,1,10]) for t = 1:H_mpc],
    q = [[Diagonal(1.0e-1 * [0.1,10,3,3]) for t = 1:H_mpc]; [Diagonal(1.0e+2 * [5,10,3,10])   for t = 1:H_mpc]],
    u = [Diagonal(1.0e-1 * [1e0, 5e-1]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-3 * ones(model.dim.b)) for t = 1:H_mpc])

# Policy
p = linearized_mpc_policy(ref_traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = n_opts,
    mpc_opts = mpc_opts,
    )

# Initial configurations
q0_sim = deepcopy(SVector{model.dim.q}(sim_stair.traj.q[end-1]))
q1_sim = deepcopy(SVector{model.dim.q}(sim_stair.traj.q[end]))

sim_flip = ContactControl.simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim_flip)


# plot_surface!(vis, model_sim.env, n=400)
anim = visualize_robot!(vis, model, sim_flip.traj, sample=10, name=:Sim, α=1.0)
# anim = visualize_robot!(vis, model, ref_traj, anim=anim, name=:Ref, α=0.3)
idx = [150, 190, 325, 450, 500, 550, 640]
α   = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 1.0]
hopper_parkour_ghost!(vis, sim_flip, sim_flip.traj, ref_traj_full, idx = idx, α = α, _name = "_flip")

# name = Symbol("Hopper" * "$t"  * "flip")
# build_robot!(vis, sim_flip.model, name=name, α = 1.0)
#
# set_robot!(vis, sim_flip.model, sim_flip.traj.q[t], name = name)

################################################################################
# Full trajectory
################################################################################

ref_traj_full = get_trajectory(model,
    joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "parkour", "hopper_stairs_3_flip_v3.jld2"),
    load_type=:split_traj_alt)

sim_traj_full = [sim_stair.traj.q[1:end-2]; sim_flip.traj.q]
anim = visualize_robot!(vis, model, sim_traj_full[1:10:end], name=:Sim, α=1.0)
anim = visualize_robot!(vis, model, ref_traj_full.q, anim=anim, name=:Ref, α=0.3)







# filename = "hopper_full_traj"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)

function hopper_parkour_ghost!(vis, sim, traj, ref_traj;
    idx = [1, 35, 50, 110, 130, 190, 210, 265, 270, 284, 295, 300, 305, 320],
    α = [convert(Float64, i / length(idx)) for i = 1:length(idx)],
    _name = "_stairs")

    plot_lines!(vis, sim.model, ref_traj.q[1:1:end], size = 5, offset = -0.5)
    stairs!(vis)
    settransform!(vis["/Cameras/default"],
            compose(Translation(0.0, -95.0, -1.0), LinearMap(RotY(0.0 * π) * RotZ(-π / 2.0))))
    setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 20)

    for (i, t) in enumerate(idx)
        name = Symbol("Hopper" * "$i" * _name)
        build_robot!(vis, sim.model, name=name, α = α[i])
        set_robot!(vis, sim.model, traj.q[t], name = name)
    end
end

# hopper_parkour_ghost!(vis, sim, ref_traj)

"Hopper" * "$t" * "_stairs"
