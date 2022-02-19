include(joinpath(@__DIR__, "..", "dynamics", "biped", "visuals.jl"))
T = Float64
vis = Visualizer()
render(vis)
open(vis)

# get hopper model
model_sim = get_model("biped", surf="flat")
model = get_model("biped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("biped", "gait1", load_type=:split_traj, model=model)

ref_traj_copy = deepcopy(ref_traj)

H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 1000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingObjective(H_mpc, model.dim,
    # q = [Diagonal(1e-1 * [1.0, 0.01, 0.5, .15, .15, .15, .15, .01, .01]) for t = 1:H_mpc],
    q = [Diagonal(1e-1 * [1.0, 0.01, 0.05, 1.5, 1.5, .15, .15, .0005, .0005]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [10; 1; 10; ones(nu-5); 10; 10]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

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
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
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
    # d = open_loop_disturbances([rand(model.dim.w) .* w_amp for i=1:H_sim]),
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    )

@time status = simulate!(sim)

l = 9
plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i][l:l], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.([q[l:l] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([b[1:nb] for b in sim.traj.b]*N_sample)...)', color=:red, linewidth=1.0)

plot_lines!(vis, model, sim.traj.q[1:N_sample:end])
plot_surface!(vis, model_sim.env)
anim = visualize_robot!(vis, model_sim, sim.traj)
anim = visualize_force!(vis, model_sim, sim.traj, anim=anim, h=h_sim)
