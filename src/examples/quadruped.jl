include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
vis = Visualizer()
open(vis)
const ContactControl = Main



s = get_simulation("quadruped", "flat_2D_lc", "flat")
model = s.model
env = s.env

ref_traj_ = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))
ref_traj = deepcopy(ref_traj_)
update_friction_coefficient!(ref_traj, model, env)

R = Dict{Symbol,Vector}(
	:dyn => [],
	:imp => [],
	:mdp => [],
	:fri => [],
	:bimp => [],
	:bmdp => [],
	:bfri => [],
	)
nz = num_var(model, env)
dyn = index_dyn(model, env, quat = true)
imp = index_imp(model, env, quat = true)
mdp = index_mdp(model, env, quat = true)
fri = index_fri(model, env, quat = true)
bimp = index_bimp(model, env, quat = true)
bmdp = index_bmdp(model, env, quat = true)
bfri = index_bfri(model, env, quat = true)
for t = 1:ref_traj.H

	z = deepcopy(ref_traj.z[t])
	θ = deepcopy(ref_traj.θ[t])
	r = zeros(nz)
	r = residual(model, env, z, θ, 0.0)
	s.res.r!(r, z, θ, 1e-6)
	push!(R[:dyn], norm(r[dyn]))
	push!(R[:imp], norm(r[imp]))
	push!(R[:mdp], norm(r[mdp]))
	push!(R[:fri], norm(r[fri]))
	push!(R[:bimp], norm(r[bimp]))
	push!(R[:bmdp], norm(r[bmdp]))
	push!(R[:bfri], norm(r[bfri]))
end
plt = plot()
plot!(plt, R[:dyn], label = "dyn")
plot!(plt, R[:imp], label = "imp")
plot!(plt, R[:mdp], label = "mdp")
plot!(plt, R[:fri], label = "fri")
plot!(plt, R[:bimp], label = "bimp")
plot!(plt, R[:bmdp], label = "bmdp")
plot!(plt, R[:bfri], label = "bfri")

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 10000 #4000 #3000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingObjective(model, env, H_mpc,
    q = [Diagonal(1e-2 * [1.0; 0.02; 0.25; 0.25 * ones(model.dim.q-3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	mode = :configuration,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
		# β_init = 1.0e-5,
        # solver = :ldl_solver,
        # verbose=true,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions())


# p = linearized_mpc_policy(ref_traj, s, obj,
#     H_mpc = H_mpc,
#     N_sample = N_sample,
#     κ_mpc = κ_mpc,
# 	# mode = :configurationforce,
# 	mode = :configuration,
# 	ip_type = :mehrotra,
#     n_opts = NewtonOptions(
# 		solver = :lu_solver,
# 		r_tol = 3e-4,
# 		max_iter = 5,
# 		max_time = ref_traj.h, # HARD REAL TIME
# 		),
#     mpc_opts = LinearizedMPCOptions(
#         # live_plotting=true,
#         # altitude_update = true,
#         # altitude_impact_threshold = 0.05,
#         # altitude_verbose = true,
#         ),
# 	ip_opts = MehrotraOptions(
# 		max_iter = 100,
# 		verbose = false,
# 		r_tol = 1.0e-4,
# 		κ_tol = 1.0e-4,
# 		diff_sol = true,
# 		# κ_reg = 1e-3,
# 		# γ_reg = 1e-1,
# 		solver = :empty_solver,
# 		),
#     )


p.im_traj.ip[10].opts.κ_tol


q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = ContactControl.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
		undercut = Inf,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-7,
        diff_sol = true),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim, verbose = true)
p.im_traj.ip[10].opts.r_tol
p.im_traj.ip[10].opts.κ_tol

################################################################################
# Timing result
################################################################################
process!(sim)
# Time budget
ref_traj.h
# Time used on average
sim.stats.μ_dt
# Speed ratio
H_sim * h_sim / sum(sim.stats.dt)


# @elapsed status = ContactControl.simulate!(sim)
# @profiler status = ContactControl.simulate!(sim)

plot_lines!(vis, model, sim.traj.q[1:1:end])
plot_surface!(vis, env, ylims=[0.3, -0.05])
anim = visualize_meshrobot!(vis, model, sim.traj, sample=5)
# anim = visualize_robot!(vis, model, sim.traj, anim=anim)
anim = visualize_force!(vis, model, env, sim.traj, anim=anim, h=h_sim, sample = 5)

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
