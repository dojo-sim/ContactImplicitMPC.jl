include(joinpath(@__DIR__, "..", "dynamics", "flamingo", "visuals.jl"))
T = Float64
vis = Visualizer()
# render(vis)
open(vis)

s = get_simulation("flamingo", "flat_2D_lc", "flat")
model = s.model
env = s.env
const ContactControl = Main
ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 15
h_sim = h / N_sample
H_sim = 3000 #35000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
    v = [Diagonal(1e-3 * [1e0,1,1e4,1,1,1,1,1e4,1e4]) for t = 1:H_mpc],
    q = [Diagonal(1e-1 * [3e2, 1e-6, 3e2, 1, 1, 1, 1, 0.1, 0.1]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [0.1; 0.1; 0.3; 0.3; ones(model.dim.u-6); 2; 2]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	ip_type = :interior_point,
	mode = :configurationforce,
	# mode = :configuration,
	# ip_type = :mehrotra,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        altitude_impact_threshold = 0.02,
        altitude_verbose = true,
        )
    )

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	# mode = :configurationforce,
	mode = :configuration,
	ip_type = :mehrotra,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5,
		max_time = ref_traj.h, # HARD REAL TIME
		),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        altitude_impact_threshold = 0.02,
        altitude_verbose = true,
        ),
	ip_opts = MehrotraOptions(
		max_iter_inner = 100,
		# verbose = true,
		r_tol = 1.0e-4,
		κ_tol = 1.0e-4,
		diff_sol = true,
		# κ_reg = 1e-3,
		# γ_reg = 1e-1,
		solver = :empty_solver,
		),
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

# sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
#     p = p,
#     ip_opts = MehrotraOptions(
#         r_tol = 1.0e-8,
#         κ_init = 1.0e-8,
#         κ_tol = 2.0e-8),
#     sim_opts = SimulatorOptions(warmstart = true),
# 	ip_type = :mehrotra,
#     )

sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true),
	ip_type = :interior_point,
    )

telap = @elapsed status = simulate!(sim, verbose = true)
# @profiler status = simulate!(sim, verbose = true)
H_sim * h_sim / (telap * 0.5)




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

plot_surface!(vis, env, xlims=[-0.5, 1.5], ylims = [-0.5, 0.5])
plot_lines!(vis, model, sim.traj.q)
anim = visualize_robot!(vis, model, sim.traj, sample=10)
anim = visualize_meshrobot!(vis, model, sim.traj, sample=10)
anim = visualize_force!(vis, model, env, sim.traj, anim=anim, h=h_sim, sample=10)





















# Test robustness
s = get_simulation("flamingo", "flat_2D_lc", "flat")
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

im_traj = ImplicitTraj(ref_traj, s,
	ip_type = :mehrotra,
	opts = MehrotraOptions(
		max_iter_inner = 100,
		κ_init = 1e-8,
		κ_tol = 2.0 * 1e-8,
		r_tol = 1.0e-8,
		diff_sol = true,
		verbose = true,
		solver = :empty_solver))

im_traj.ip[1].z .= ref_traj.z[1]
im_traj.ip[1].θ .= ref_traj.θ[1] .+ 1.0*[ones(length(ref_traj.θ[1])-2); zeros(2)]
interior_point_solve!(im_traj.ip[1])
im_traj.ip[1].iterations

cnt = 0
itl = []
sul = 0
for tt = 1:1000
	t = (tt % ref_traj.H) + 1
    Random.seed!(t)
	im_traj.ip[t].z .= ref_traj.z[t]
	im_traj.ip[t].θ .= ref_traj.θ[t] .+ 1e-3*[0.5 .- rand(length(ref_traj.θ[t])-2); zeros(2)]
	interior_point_solve!(im_traj.ip[t])
	it = im_traj.ip[t].iterations
	su = it < 100
    sul += Int(su)
    push!(itl, it)
end
sul
mean(itl)

filename = "flamingo_mehrotra"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)
