const ContactControl = Main
vis = Visualizer()
open(vis)
# render(vis)






dir = joinpath(module_dir(), "src", "dynamics", "planarpush")
model = deepcopy(planarpush)

path_base = joinpath(dir, "dynamics/base.jld2")
path_dyn = joinpath(dir, "dynamics/dynamics.jld2")

expr_base = generate_base_expressions(model, M_analytical = false)
save_expressions(expr_base, path_base, overwrite=true)
instantiate_model_base!(model, path_base)

expr_dyn = generate_dynamics_expressions(model)
save_expressions(expr_dyn, path_dyn, overwrite=true)
instantiate_model_dynamics!(model, path_dyn)




dir_model = joinpath(module_dir(), "src/dynamics/planarpush")
dir_sim   = joinpath(module_dir(), "src/simulation/planarpush")
model = deepcopy(planarpush)
env = deepcopy(flat_3D_lc)
sim = Simulation(model, env)

path_base = joinpath(dir_model, "dynamics/base.jld2")
path_dyn = joinpath(dir_model, "dynamics/dynamics.jld2")
path_res = joinpath(dir_sim, "flat/residual.jld2")
path_jac = joinpath(dir_sim, "flat/jacobians.jld2")

expr_contact = generate_contact_expressions(model, env, jacobians=false)
instantiate_contact_methods!(sim.con, expr_contact, jacobians=:full)
expr_contact = generate_contact_expressions(model, env, jacobians=true)
instantiate_contact_methods!(sim.con, expr_contact, jacobians=:approx)


# expr_base = generate_sim_base_expressions(model, env,)
# save_expressions(expr_base, path_sim_base, overwrite=true)
# instantiate_sim_base!(model, path_sim_base)
#
# expr_dyn = generate_sim_dynamics_expressions(model, env)
# save_expressions(expr_dyn, path_dyn, overwrite=true)
# instantiate_dynamics!(model, path_dyn)


expr_res, rz_sp, rθ_sp = generate_residual_expressions(sim.model, sim.env)
save_expressions(expr_res, path_res, overwrite=true)
@save path_jac rz_sp rθ_sp
instantiate_residual!(sim, path_res, path_jac)



include(joinpath(@__DIR__, "..", "dynamics", "planarpush", "visuals.jl"))

s = get_simulation("planarpush", "flat_3D_lc", "flat")
model = s.model
nq = model.dim.q
nu = model.dim.u
nw = model.dim.w
nc = model.dim.c
nf = friction_dim(s.env)
nb = nc * nf
nz = num_var(model, s.env)
nθ = num_data(model)


# time
h = 0.01
H = 50 #550
N_sample = 1

# initial conditions
q0 = @SVector [0.00, 0.00, 0.0, 0.0, -0.25, 5e-3, 0.00]
q1 = @SVector [0.00, 0.00, 0.0, 0.0, -0.25, 5e-3, 0.00]

# p = open_loop_policy(fill(SVector{nu}([10.0*h, -0.0]), H*2), N_sample=N_sample)
p = open_loop_policy(fill(SVector{nu}([40*h, -0.0]), H), N_sample=N_sample)

# simulator
sim0 = simulator(s, q0, q1, h, H,
				p = p,
				ip_opts = ContactControl.InteriorPointOptions(
					r_tol = 1.0e-8, κ_init=1e-6, κ_tol = 2.0e-6),
				sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
status = simulate!(sim0)
@test status

visualize_robot!(vis, model, sim0.traj)

plot(hcat([x[1:nq] for x in sim0.traj.q]...)')
plot(hcat([x[1:nu] for x in sim0.traj.u]...)')
plot(hcat([x[1:nw] for x in sim0.traj.w]...)')
plot(hcat([x[1:nc] for x in sim0.traj.γ]...)')
plot(hcat([x[1:nb] for x in sim0.traj.b]...)')
# plot(hcat([x[1:8] for x in sim0.traj.b]...)')


ref_traj = deepcopy(sim0.traj)
ref_traj.H

# MPC
N_sample = 1
H_mpc = 40
h_sim = h / N_sample
H_sim = 500

# barrier parameter
κ_mpc = 1.0e-4

# Aggressive
obj = TrackingVelocityObjective(model, s.env, H_mpc,
	q = [Diagonal(1.0e-4 * [1e2, 1e2, 1e-1, 1e-1, 1e-0, 1e-0, 1e-0,]) for t = 1:H_mpc],
	v = [Diagonal(1.0e-3 * ones(model.dim.q)) for t = 1:H_mpc],
	u = [Diagonal(1.0e-6 * ones(model.dim.u)) for t = 1:H_mpc],
	γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
	b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:H_mpc])

#
obj = TrackingVelocityObjective(model, s.env, H_mpc,
	q = [Diagonal(1.0e-2 * [1e2, 1e2, 1e-1, 1e-1, 1e-0, 1e-0, 1e-0,]) for t = 1:H_mpc],
	v = [Diagonal(1.0e-1 * ones(model.dim.q)) for t = 1:H_mpc],
	u = [Diagonal(1.0e-4 * ones(model.dim.u)) for t = 1:H_mpc],
	γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
	b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        solver = :ldl_solver,
		# max_iter = 5),
		max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
		# live_plotting=true
		))

# p = open_loop_policy(fill(SVector{nu}([40*h, -0.0]), H*2), N_sample=N_sample)
using Random
Random.seed!(100)
d = open_loop_disturbances([[0.05*rand(), 0.0, 0.4*rand()] for t=1:H_sim])

q0_sim = @SVector [0.00, 0.00, 0.0, 0.0, -0.25, 5e-3, 0.00]
q1_sim = @SVector [0.00, 0.00, 0.0, 0.0, -0.25, 5e-3, 0.00]

sim = ContactControl.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
	uL=-0.9*ones(model.dim.u),
	uU=+2.0*ones(model.dim.u),
    p = p,
	d = d,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6,
		diff_sol=true),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

sim.traj.u

@time status = ContactControl.simulate!(sim)
sim0
anim = visualize_robot!(vis, model, sim.traj, name=:sim, sample = N_sample, α=1.0)
anim = visualize_robot!(vis, model, ref_traj, name=:ref, anim=anim, sample = 1, α=0.5)

plot(hcat([x[1:2] ./ N_sample for x in ref_traj.u]...)')
scatter!(hcat([x[1:2] for x in sim.traj.u]...)')
sim.traj.u
# filename = "planarpush_precise"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
