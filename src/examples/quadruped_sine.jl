const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
render(vis)

s_sim = get_simulation("quadruped", "sine1_2D_lc", "sinusoidal")
s = get_simulation("quadruped", "flat_2D_lc", "flat")
model = s.model
env = s.env

s_sim.model.μ_world
s.model.μ_world

ref_traj_ = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))
ref_traj = deepcopy(ref_traj_)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 1500 #4000 #5000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
    q = [Diagonal(1e-2 * [1.00; 0.02; 0.1; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    # v = [Diagonal(1e-4 * [0.01; 0.01; 0.01; 1.00 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    # q = [Diagonal(1e-2 * [1.0; 0.02; 0.01; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    # q = [Diagonal(1e-2 * [1.0; 0.02; 0.01; 0.5 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    # q = [Diagonal(1e-2 * [2.0; 0.02; 0.01; 0.5 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    # q = [Diagonal(1e-2 * [10.0; 0.02; 0.25; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    # q = [Diagonal(1e-2 * [10.0; 0.02; 0.25; 0.5 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    u = [Diagonal(1e-2 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5,
        verbose = true,
        ),
    mpc_opts = LinearizedMPCOptions(
        live_plotting=true,
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
sim = simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    # d = open_loop_disturbances([rand(model.dim.w) .* w_amp for i=1:H_sim]),
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = SimulatorOptions(warmstart = true)
    )

@time status = simulate!(sim)

# plt = plot(layout=(3,1), legend=false)
# plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
# plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

plot_surface!(vis, s_sim.env, ylims=[0.3, -0.05])
plot_lines!(vis, model, sim.traj.q[1:5:end])
anim = visualize_meshrobot!(vis, model, sim.traj, sample=5)
# anim = visualize_robot!(vis, model, sim.traj, anim=anim)
anim = visualize_force!(vis, model, sim.traj, anim=anim, h=h_sim)

# Display ghosts
t_ghosts = [1, 1333, 2666]
mvis_ghosts = []
for (i,t) in enumerate(t_ghosts)
    α = i/(length(t_ghosts)+1)
    name = Symbol("ghost$i")
    mvis = build_meshrobot!(vis, model, name=name, α=α)
    push!(mvis_ghosts, mvis)
end

for (i,t) in enumerate(t_ghosts)
    name = Symbol("ghost$i")
    set_meshrobot!(vis, mvis_ghosts[i], model, sim.traj.q[t], name=name)
end
