const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# Define a special stride
function get_stride(model::Hopper2D, traj::ContactTraj)
    stride = zeros(SizedVector{model.dim.q})
    stride[1:2] = traj.q[end-1][1:2] - traj.q[1][1:2]
    return stride
end

# get hopper model
model_sim = get_model("hopper_2D", surf="stairs")
model = get_model("hopper_2D", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj_ = get_trajectory(model,
    joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "parkour", "hopper_stairs_3_flip_v3.jld2"),
    # joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "parkour", "hopper_stair_v3.jld2"),
    # joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "parkour", "hopper_tall_flip5.jld2"),
    load_type=:split_traj_alt)
ref_traj = deepcopy(ref_traj_)
plot(hcat([γ[1:1] for γ in ref_traj.γ[1:end]]...)')

# Solve altitude issue
l1 = 41
l2 = 80
l3 = 80
l4 = 62
l5 = (ref_traj.H+2) - (l1+l2+l3+l4)
offset = [zeros(l1); 0.25*ones(l2), 0.50*ones(l3), 0.25*ones(l4)]
# for
plot(hcat([q[1:1] for q in ref_traj.q]...)')
plot!(hcat([q[4:4] for q in ref_traj.q]...)')

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 10
H_mpc = 10
h_sim = h / N_sample
H_sim = 1400

# barrier parameter
κ_mpc = 1.0e-4

# obj = TrackingObjective(H_mpc, model.dim,
#     q = [Diagonal(1.0e-1 * [10,30,1,3])   for t = 1:H_mpc],
#     u = [Diagonal(1.0e-0 * [1e-1, 1e0]) for t = 1:H_mpc],
#     γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
#     b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

# obj = TrackingObjective(H_mpc, model.dim,
#     q = [[Diagonal(1.0e-0 * [1,1,1,1])   for t = 1:H_mpc]; [Diagonal(1.0e+1 * [1,1,1,1])   for t = 1:H_mpc]],
#     u = [Diagonal(1.0e-1 * [1e-0, 1e1]) for t = 1:H_mpc],
#     γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
#     b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])
obj = TrackingVelocityObjective(H_mpc, model.dim,
    v = [Diagonal(1e-3 * [1,1,1,10]) for t = 1:H_mpc],
    # q = [[Diagonal(1.0e-1 * [0.1,10,3,3]) for t = 1:H_mpc]; [Diagonal(1.0e+2 * [0.1,10,3,10])   for t = 1:H_mpc]],
    q = [[Diagonal(1.0e-0 * [0.3,1,1,1])   for t = 1:H_mpc]; [Diagonal(1.0e+2 * [0.3,1,1,1])   for t = 1:H_mpc]],
    u = [Diagonal(1.0e-1 * [1e0, 1e-0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-3 * ones(model.dim.b)) for t = 1:H_mpc])

obj = TrackingVelocityObjective(H_mpc, model.dim,
    v = [Diagonal(1e-3 * [1,1,1,10]) for t = 1:H_mpc],
    # q = [[Diagonal(1.0e-1 * [0.1,10,3,3]) for t = 1:H_mpc]; [Diagonal(1.0e+2 * [0.1,10,3,10])   for t = 1:H_mpc]],
    q = [[Diagonal(1.0e-0 * [0.3,1,1,1])   for t = 1:H_mpc]; [Diagonal(1.0e+2 * [0.3,1,1,1])   for t = 1:H_mpc]],
    u = [Diagonal(1.0e-1 * [1e0, 1e-0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-3 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        altitude_update = true,
        altitude_impact_threshold = 0.1,
        altitude_verbose = true,
        )
    )


q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = ContactControl.simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim)


# plot_surface!(vis, model_sim.env, n=400)
anim = visualize_robot!(vis, model, sim.traj, sample=10, name=:Sim, α=1.0)
anim = visualize_robot!(vis, model, ref_traj, anim=anim, name=:Ref, α=0.3)

# anim = visualize_robot!(vis, model, sim.traj.q[1:20], name=:Sim, α=1.0)
# anim = visualize_robot!(vis, model, ref_traj.q[1:20], anim=anim, name=:Ref, α=0.3)

set_robot!(vis, model_sim, sim.traj.q[1], name=:Sim)
set_robot!(vis, model_sim, sim.traj.q[2], name=:Sim)
set_robot!(vis, model_sim, sim.traj.q[3], name=:Sim)

set_robot!(vis, model_sim, ref_traj.q[1], name=:Ref)
set_robot!(vis, model_sim, ref_traj.q[2], name=:Ref)
set_robot!(vis, model_sim, ref_traj.q[3], name=:Ref)


plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)


plot(hcat([γ[1:1] for γ in sim.traj.γ[400:430]]...)')
plot(hcat([q[1:4] for q in ref_traj.q[1:200]]...)')
plot(hcat([q[1:4] for q in sim.traj.q[1:1900]]...)')
# plot!(hcat([q[1:4] for q in sim.traj.q[1:10]]...)')
# filename = "hopper_flip4"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
