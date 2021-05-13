include(joinpath(@__DIR__, "..", "dynamics", "flamingo", "visuals.jl"))
T = Float64
vis = Visualizer()
render(vis)
open(vis)

# get model
include(joinpath(pwd(), "src/simulator/terrain/piecewise.jl"))
include(joinpath(pwd(), "src/dynamics/flamingo/model.jl"))
flamingo_piecewise = Flamingo(Dimensions(nq, nu, nw, nc, nb),
			  g, μ_world, μ_joint,
			  l_torso, d_torso, m_torso, J_torso,
			  l_thigh, d_thigh, m_thigh, J_thigh,
			  l_calf, d_calf, m_calf, J_calf,
			  l_foot, d_foot, m_foot, J_foot,
			  l_thigh, d_thigh, m_thigh, J_thigh,
			  l_calf, d_calf, m_calf, J_calf,
			  l_foot, d_foot, m_foot, J_foot,
			  zeros(nc),
			  BaseMethods(), DynamicsMethods(), ContactMethods(),
			  ResidualMethods(), ResidualMethods(),
			  SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			  SVector{nq}([zeros(3); 0.0 * μ_joint * ones(nq - 3)]),
			  Environment{R2}(piecewise_smoothed, d_piecewise_smoothed),
			  )

model_sim = deepcopy(flamingo_piecewise)

plot(x, model_sim.env.surf.(x))
dir = joinpath(pwd(), "src/dynamics/flamingo")

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")
path_res = joinpath(dir, "piecewise/residual.jld2")
path_jac = joinpath(dir, "piecewise/sparse_jacobians.jld2")
path_linearized = joinpath(dir, "piecewise/linearized.jld2")

instantiate_base!(model_sim, path_base)

expr_dyn = generate_dynamics_expressions(model_sim, derivs = true)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_dynamics!(model_sim, path_dyn, derivs = true)

expr_res, rz_sp, rθ_sp = generate_residual_expressions(model_sim, jacobians = :approx)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(model_sim, path_res, jacobians = :approx)
model_sim.spa.rz_sp = copy(rz_sp)
model_sim.spa.rθ_sp = copy(rθ_sp)

model = get_model("flamingo", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd
nz = num_var(model)
nθ = num_data(model)

# get trajectory
# ref_traj = get_trajectory("flamingo", "gait1", load_type=:split_traj_alt, model=model)
# ref_traj = get_trajectory("flamingo", "gait_forward_36_3", load_type=:split_traj_alt, model=model)
ref_traj = get_trajectory("flamingo", "gait_forward_36_4", load_type=:split_traj_alt, model=model)


H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 15
h_sim = h / N_sample
H_sim = 2500#15000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingVelocityObjective(H_mpc, model.dim,
    v = [Diagonal(1e-3 * [1e0,1,1e4,1,1,1,1,1e4,1e4]) for t = 1:H_mpc],
    q = [Diagonal(1e-1 * [3e2, 1e-6, 3e2, 1, 1, 1, 1, 0.1, 0.1]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [0.1; 0.1; 0.3; 0.3; ones(nu-6); 2; 2]) for t = 1:H_mpc],
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
        altitude_impact_threshold = 0.02,
        altitude_verbose = true,
        )
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    )

@time status = simulate!(sim)

# save trajectory
@save joinpath(pwd(), "src/dynamics/flamingo/simulations/piecewise.jld2") sim
@load joinpath(pwd(), "src/dynamics/flamingo/simulations/piecewise.jld2") sim

l = 9
lu = 1
plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.([q[l:l] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][lu:lu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[3,1], hcat(Vector.(vcat([fill(ref_traj.γ[i][1:nc], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[lu:lu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([b[1:nb] for b in sim.traj.b]*N_sample)...)', color=:red, linewidth=1.0)

plot_lines!(vis, model, sim.traj.q[1:N_sample:end], offset=-0.01)
plot_surface!(vis, model_sim.env, xlims=[-1, 9])
anim = visualize_meshrobot!(vis, model_sim, sim.traj, sample=10)
anim = visualize_force!(vis, model_sim, sim.traj, anim=anim, h=h_sim, sample=10)
open(vis)

filename = "flamingo_40_steps_sine"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)

# Ghost
flamingo_ghost!(vis, sim, piecewise_smoothed)

# Animation
anim, shift_traj = flamingo_animation!(vis, sim, piecewise_smoothed)
anim = visualize_force!(vis, sim.model, shift_traj, anim=anim, h=h_sim, sample=10)
