const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "racecar", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
s = get_simulation("racecar", "flat_3D_lc", "flat")
s_load = get_simulation("racecar", "flat_3D_lc", "payload",
    model_variable_name="racecar_payload",
    dynamics_name="dynamics_payload")
model = s.model
env = s.env
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c

s_load.model.mb
s.model.mb

s_load.model.μ_world
s.model.μ_world

# get trajectory
ref_traj_ = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/racecar/gaits/drift.jld2"),
    load_type = :joint_traj))
ref_traj = deepcopy(ref_traj_)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = (H - H_mpc)*N_sample

# barrier parameter
κ_mpc = 2.0e-4

obj = TrackingObjective(model, env, H_mpc,
    q = [Diagonal(1.0e-1 * [100; 100; 1;  10; 10; 10; 10; ones(nq-7)])   for t = 1:H_mpc],
    u = [Diagonal(3.0e-1 * [20; 1; 1; 1; 1;]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(nc * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        # altitude_impact_threshold = 0.05,
        # altitude_verbose = true,
        )
    )

# p = open_loop_policy(ref_traj.u, N_sample = N_sample)

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = ContactControl.simulator(s_load, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8,
        diff_sol = true),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))


status = ContactControl.simulate!(sim, verbose = true)


plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

short_traj = sub_traj(ref_traj, Vector(1:Int(floor(H_sim/N_sample))))

plot_surface!(vis, env, n=200, xlims=[-1.0, 16.0], ylims=[-1.0, 18.0])
plot_lines!(vis, model, short_traj.q, col = false, name = :ref)
plot_lines!(vis, model, uncontrolled_traj.q, name = :uncontrolled)
plot_lines!(vis, model, controlled_traj.q, col = false, name = :controlled)
# plot_lines!(vis, model, sim.traj.q, name = :mpc)
anim = visualize_robot!(vis, model, short_traj, sample=1, α = 0.4, name = :ref)
visualize_robot!(vis, model, uncontrolled_traj, sample=N_sample, α = 0.4, anim = anim, name = :uncontrolled)
visualize_robot!(vis, model, controlled_traj, sample=N_sample, α = 1.0, anim = anim, name = :controlled)
# visualize_robot!(vis, model, sim.traj, sample=N_sample, α = 1.0, anim = anim, name = :mpc)
anim = visualize_force!(vis, model, env, sim.traj, anim=anim, h=h_sim, sample = 1)

plot(hcat([u[[1,2]] for u in Vector.(sim.traj.u)[1:end]]...)')


controlled_traj = deepcopy(sim.traj)
# uncontrolled_traj = deepcopy(sim.traj)

filename = "racecar_topdown"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/video/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/video/$filename.mp4",
    "/home/simon/Documents/video/$filename.gif", overwrite=true)

# const ContactControl = Main
