const ContactImplicitMPC = Main
include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
# render(vis)

# get simulation
s_no_load = get_simulation("quadruped", "flat_2D_lc", "flat")
s_load = get_simulation("quadruped", "flat_2D_lc", "payload",
    model_variable_name="quadruped_payload",
    dynamics_name="dynamics_payload")
model_no_load = s_no_load.model
model_load = s_load.model
env_no_load = s_no_load.env
env_load = s_load.env


nq = model.dim.q
nu = model.dim.u
nc = model.dim.c

@testset "Check Quadruped Load" begin
    @test (model_load.m_torso - model_no_load.m_torso) == 3.0
    @test (model_load.J_torso - model_no_load.J_torso) == 0.03
    q_ = rand(nq)
    @test (M_func(model_load, q_)[1,1] - M_func(model_no_load, q_)[1,1]) == 3.0
    @test (M_fast(model_load, q_)[1,1] - M_fast(model_no_load, q_)[1,1]) == 3.0
end

# get trajectory
ref_traj_ = deepcopy(ContactImplicitMPC.get_trajectory(s_no_load.model, s_no_load.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))
ref_traj = deepcopy(ref_traj_)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 2100

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingObjective(model_no_load, env_no_load, H_mpc,
    # q = [Diagonal(1e-2 * [10.0; 0.02; 0.25; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    q = [Diagonal(1e-2 * [10.0; 0.02; 0.25; 0.5 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s_no_load, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        solver = :ldl_solver,
        verbose = false,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
        )
    )

# Test open loop policy to see the impact of the load
# p = open_loop_policy(deepcopy(ref_traj.u), N_sample=N_sample)

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim_load = simulator(s_load, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6,
        diff_sol = false,
        verbose = false),
    sim_opts = SimulatorOptions(warmstart = true))
time = @elapsed status = ContactImplicitMPC.simulate!(sim_load)


sim_no_load = simulator(s_no_load, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6,
        diff_sol = false,
        verbose = false),
    sim_opts = SimulatorOptions(warmstart = true))
time = @elapsed status = ContactImplicitMPC.simulate!(sim_no_load)



# plt = plot(layout=(3,1), legend=false)
# plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
# plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)
# render(vis)
plot_surface!(vis, s_load.env, ylims=[0.3, -0.05])
plot_lines!(vis, model, sim_no_load.traj.q[1:1:end], name=:NoPayload, offset=-0.15)
plot_lines!(vis, model, sim_load.traj.q[1:1:end], name=:Payload, offset=-0.15)
ext_ref_traj = repeat_ref_traj(ref_traj, 7; idx_shift = (1:1))
plot_lines!(vis, model, ext_ref_traj.q, offset=-0.17, name=:Ref, col=false,)

anim = visualize_meshrobot!(vis, s_load.model, sim_load.traj, anim=anim, sample=1, name=:Payload)
anim = visualize_meshrobot!(vis, s_no_load.model, sim_no_load.traj, sample=1, name=:NoPayload)
anim = visualize_payload!(vis, model, sim_load.traj, anim=anim, sample=1, name=:Payload, object=:mesh)
# anim = visualize_force!(vis, model, sim_no_load.traj, anim=anim, sample=5, h=h_sim, name=:NoPayload)
# anim = visualize_force!(vis, model, sim_load.traj, anim=anim, sample=5, h=h_sim, name=:Payload)

anim = visualize_meshrobot!(vis, model, sim_load.traj, sample=1, name=:Payload)
anim = visualize_payload!(vis, model, sim_load.traj, anim=anim, sample=1, name=:Payload, object=:mesh)
# settransform!(vis["/Cameras/default"],
#         compose(Translation(0.0, -25.0, -1.0), LinearMap(RotY(0.0 * π) * RotZ(-π / 2.0))))
# setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 20)


# Display ghosts
t_ghosts = [1]
mvis_ghosts = []
for (i,t) in enumerate(t_ghosts)
    α = i/(length(t_ghosts)+1)
    name = Symbol("ghost$i")
    mvis = build_meshrobot!(vis, model_load, name=name, α=α)
    build_load!(vis, model_load, name=name, α=α)
    push!(mvis_ghosts, mvis)
end

for (i,t) in enumerate(t_ghosts)
    name = Symbol("ghost$i")
    set_meshrobot!(vis, mvis_ghosts[i], model, sim_load.traj.q[t], name=name)
    set_load!(vis, model, sim_load.traj.q[t], name=name)
end

anim = visualize_load!(vis, model, sim_load.traj, anim=anim, sample=5, name=:Payload)


filename = "quadruped_payload_slow"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)
