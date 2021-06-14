include(joinpath(@__DIR__, "..", "dynamics", "hopper_3D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
model_sim = get_model("hopper_3D", surf="sinusoidal")
model = get_model("hopper_3D")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
# ref_traj = get_trajectory("hopper_3D", "gait_in_place", load_type=:joint_traj)
ref_traj = get_trajectory("hopper_3D", "gait_forward", load_type=:joint_traj)
ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_3D/gaits/gait_forward.jld2"),
    load_type = :joint_traj))
# time
H = ref_traj.H
h = ref_traj.h
N_sample = 10
H_mpc = 20
h_sim = h / N_sample
H_sim = 12000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingObjective(H_mpc, model.dim,
    q = [Diagonal(1.0e-1 * [10,10,0.1,5e+1,5e+1,5e+1,10])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-1, 1e-1, 1e1]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

obj = TrackingObjective(H_mpc, model.dim,
    q = [Diagonal(1.0e-1 * [1,1,0.1,5e+1,5e+1,5e+1,10])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-1, 1e-1, 1e1]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # altitude_update = false,
        altitude_update = true,
        altitude_impact_threshold = 0.15,
        altitude_verbose = true,
        # live_plotting = true,
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

plot_surface!(vis, model_sim.env, n=200, ylims=[-0.3, 0.9])
visualize_robot!(vis, model, sim.traj)



plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

visualize_robot!(vis, model, sim.traj, sample=10)



filename = "hopper_3d_sine"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)

const ContactControl = Main
