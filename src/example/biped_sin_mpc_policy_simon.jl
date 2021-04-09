include(joinpath(@__DIR__, "..", "dynamics", "biped", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
# model = get_model("quadruped")
model_sim = get_model("biped", surf="flat")
model = get_model("biped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
# ref_traj = get_trajectory("biped", "gait5", load_type=:split_traj, model=model)
# ref_traj = get_trajectory("biped", "biped_gait (1)", load_type=:split_traj, model=model)
ref_traj = get_trajectory("biped", "biped_gait (2)", load_type=:split_traj, model=model)
visualize!(vis, model, ref_traj.q, Δt=20*h/N_sample, name=:mpc)
visualize!(vis, model, ref_traj.q, Δt=h, name=:mpc)

ref_traj_copy = deepcopy(ref_traj)

H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 15
h_sim = h / N_sample
H_sim = 1000

# barrier parameter
κ_mpc = 1.0e-4

cost = CostFunction(H_mpc, model.dim,
    q = [Diagonal(1e-1 * [1.0, 0.01, 0.5, .15, .15, .15, .15, .05, .05]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [10; ones(nu-3); 1; 1]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, cost,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        altitude_update = true,
        altitude_impact_threshold = 0.10,
        # altitude_verbose = true,
        )
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

w_amp = [+0.02, -0.20]
sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    # d = random_disturbances(model, w_amp, H, h)
    # d = open_loop_disturbances([rand(model.dim.w) .* w_amp for i=1:H_sim]),
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    )

# @profiler status = simulate!(sim)
@time status = simulate!(sim)
# 4.88*2300/8500/(400*h_sim)



plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([b[1:nb] for b in sim.traj.b]*N_sample)...)', color=:red, linewidth=1.0)

visualize!(vis, model, sim.traj.q[1:N_sample:end], Δt=10*h/N_sample, name=:mpc)
visualize!(vis, model, sim.traj.q[1:N_sample:end], Δt=h, name=:mpc)
# draw_lines!(vis, model, sim.traj.q[1:N_sample:end])
plot_surface!(vis, model_sim.env)


a = 10
a = 10
a = 10

filename = "biped_best_so_far"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)

# const ContactControl = Main

# qq = []
# for q in ref_traj_copy.q
#     for i = 1:N_sample
#         push!(qq, q)
#     end
# end
# plot(hcat(qq...)[1:model.dim.q, 1:end]',
#     label = "", color = :black, width = 3.0)
# plot!(hcat(sim.traj.q...)[1:model.dim.q, 1:100]',
#     label = "", color = :cyan, width = 1.0, legend = :topleft)
#


# q_test = [ref_traj.q[1] + [0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.6, 0.6, ]]
# visualize!(vis, model, q_test, Δt=10*h/N_sample, name=:mpc)
