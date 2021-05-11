const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
# render(vis)

# get hopper model
model_payload = get_model("quadruped", surf="payload", dynamics="dynamics_payload")
model_no_payload = get_model("quadruped", surf="flat")

model = get_model("quadruped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd
M_fast(model_payload, zeros(nq))
M_fast(model_no_payload, zeros(nq))

# get trajectory
ref_traj_ = get_trajectory("quadruped", "gait2", load_type=:split_traj_alt, model=model)
ref_traj = deepcopy(ref_traj_)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 2100 #220 #5000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingObjective(H_mpc, model.dim,
    q = [Diagonal(1e-2 * [10; 0.02; 0.25; 0.25 * ones(nq-3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
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
        # live_plotting=true,
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
        )
    )

# Test open loop policy to see the impact of the payload
# p = open_loop_policy(deepcopy(ref_traj.u), N_sample=N_sample)

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim_payload = simulator(model_payload, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = SimulatorOptions(warmstart = true)
    )
@time status = simulate!(sim_payload)


sim_no_payload = simulator(model_no_payload, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = SimulatorOptions(warmstart = true)
    )
@time status = simulate!(sim_no_payload)


# plt = plot(layout=(3,1), legend=false)
# plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
# plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

plot_surface!(vis, model.env, ylims=[0.3, -0.05])
plot_lines!(vis, model, sim_no_payload.traj.q[1:1:end], name=:NoPayload, offset=-0.15)
plot_lines!(vis, model, sim_payload.traj.q[1:1:end], name=:Payload, offset=-0.15)
ext_ref_traj = repeat_ref_traj(ref_traj, model, 7; idx_shift = (1:1))
plot_lines!(vis, model, ext_ref_traj.q, offset=-0.17, name=:Ref, col=false,)

anim = visualize_meshrobot!(vis, model, sim_no_payload.traj, sample=5, name=:NoPayload)
anim = visualize_meshrobot!(vis, model, sim_payload.traj, anim=anim, sample=5, name=:Payload)
anim = visualize_payload!(vis, model, sim_payload.traj, anim=anim, sample=5, name=:Payload, object=:mesh)
# anim = visualize_force!(vis, model, sim_no_payload.traj, anim=anim, sample=5, h=h_sim, name=:NoPayload)
# anim = visualize_force!(vis, model, sim_payload.traj, anim=anim, sample=5, h=h_sim, name=:Payload)

anim = visualize_meshrobot!(vis, model, sim_payload.traj, sample=5, name=:Payload)
anim = visualize_payload!(vis, model, sim_payload.traj, anim=anim, sample=5, name=:Payload, object=:mesh)
settransform!(vis["/Cameras/default"],
        compose(Translation(0.0, -25.0, -1.0), LinearMap(RotY(0.0 * π) * RotZ(-π / 2.0))))
setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 20)


# Display ghosts
t_ghosts = [1]
mvis_ghosts = []
for (i,t) in enumerate(t_ghosts)
    α = i/(length(t_ghosts)+1)
    name = Symbol("ghost$i")
    mvis = build_meshrobot!(vis, model_payload, name=name, α=α)
    build_payload!(vis, model_payload, name=name, α=α)
    push!(mvis_ghosts, mvis)
end

for (i,t) in enumerate(t_ghosts)
    name = Symbol("ghost$i")
    set_meshrobot!(vis, mvis_ghosts[i], model, sim_payload.traj.q[t], name=name)
    set_payload!(vis, model, sim_payload.traj.q[t], name=name)
end

anim = visualize_payload!(vis, model, sim_payload.traj, anim=anim, sample=5, name=:Payload)


filename = "quadruped_cardboard"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)
