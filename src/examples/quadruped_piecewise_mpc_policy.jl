include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
render(vis)

# get hopper model
model_sim = get_model("quadruped", surf="piecewise")
model = get_model("quadruped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("quadruped", "gait0", load_type=:split_traj, model=model)
ref_traj_copy = deepcopy(ref_traj)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 100

# barrier parameter
κ_mpc = 1.0e-4

cost = CostFunction(H_mpc, model.dim,
    q = [Diagonal(1e-2 * [0.02; 0.02; 1.0; 0.25 * ones(nq-3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
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
        altitude_impact_threshold = 0.05,
        # altitude_verbose = true,
        )
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))

# policy(p, ref_traj.q[2], ref_traj, 1)
#
# z0 = rand(num_var(model_sim))
# θ0 = rand(num_data(model_sim))
# r0 = zeros(num_var(model_sim))
# rz0 = zeros(num_var(model_sim), num_var(model_sim))
# rθ0 = zeros(num_var(model_sim), num_data(model_sim))
#
#
# nq = model_sim.dim.q
# nu = model_sim.dim.u
# nc = model_sim.dim.c
# nb = model_sim.dim.b
# nz = num_var(model_sim)
# nθ = num_data(model_sim)
#
# # Declare variables
# @variables z[1:nz]
# @variables θ[1:nθ]
# @variables κ
#
# # Residual
# r = residual(model_sim, z, θ, κ)
# r = Symbolics.simplify.(r)
# rz = Symbolics.jacobian(r, z, simplify = true)
# # rθ = Symbolics.jacobian(r, θ, simplify = true) # TODO: sparse version
#
# rz_sp = similar(rz, T)
# # rθ_sp = similar(rθ, T)
#
# # Build function
# expr = Dict{Symbol, Expr}()
# # expr[:r]  = build_function(r, z, θ, κ)[2]
# expr[:rz] = build_function(rz, z, θ)[2]
# # expr[:rθ] = build_function(rθ, z, θ)[2]
#
# eval(expr[:rz])(rz_sp, z0, θ0)
#
# # model_sim.res.r!(r0, z0, θ0, 0.1)
# # model_sim.res.rz!(rz0, z0, θ0)
# # model_sim.res.rθ!(rθ0, z0, θ0)
#
# p = no_policy(model_sim)
w_amp = [+0.02, -0.20]
sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    # d = open_loop_disturbances([rand(model.dim.w) .* w_amp for i=1:H_sim]),
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = SimulatorOptions(warmstart = true)
    )

@time status = simulate!(sim)

plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

# visualize!(vis, model, sim.traj.q[1:N_sample:end], Δt=10*h/N_sample, name=:mpc)
plot_lines!(vis, model, sim.traj.q[1:N_sample:end])
plot_surface!(vis, model_sim.env)
anim = visualize_robot!(vis, model_sim, sim.traj)
anim = visualize_force!(vis, model_sim, sim.traj, anim=anim, h=h_sim)

settransform!(vis["/Cameras/default"],
	    compose(Translation(0.0, 0.0, -1.0), LinearMap(RotZ(-pi / 2.0))))

# filename = "quadruped_forces"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)

# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)

# const ContactControl = Main
