model = ContactControl.get_model("inverted_pendulum")

# time
h = 0.01
H = 1000

# reference trajectory
ref_traj = contact_trajectory(H, h, model)
ref_traj.h
qref = zeros(model.dim.q)
ur = zeros(model.dim.u)
γr = zeros(model.dim.c)
br = zeros(model.dim.b)
ψr = zeros(model.dim.c)
ηr = zeros(model.dim.b)
wr = zeros(model.dim.w)

# set reference
for t = 1:H
	ref_traj.z[t] = pack_z(model, qref, γr, br, ψr, ηr)
	ref_traj.θ[t] = pack_θ(model, qref, qref, ur, wr, model.μ_world, ref_traj.h)
end

# test reference
for t = 1:H
	r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)#[model.dim.q .+ (1:model.dim.c)]
	@test norm(r) < 1.0e-4
end


# initial conditions
q0 = @SVector [-0.01]
q1 = @SVector [0.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, H,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-8),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
status = ContactControl.simulate!(sim)
@test status

# visualize
include(joinpath(@__DIR__, "..", "dynamics", "inverted_pendulum", "visuals.jl"))
vis = Visualizer()
render(vis)
add_walls!(vis, model)
anim = visualize_robot!(vis, model, sim.traj)
# anim = visualize_force!(vis, model, sim.traj, anim=anim, h=h_sim)

# MPC
N_sample = 5
H_mpc = 5
h_sim = h / N_sample
H_sim = 1000

# barrier parameter
κ_mpc = 1.0e-4

cost = CostFunction(H_mpc, model.dim,
    q = [Diagonal(100.0 * ones(model.dim.q)) for t = 1:H_mpc],
    u = [Diagonal(1.0e-2 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, cost,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        solver = :ldl_solver,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions())

d = impulse_disturbances([0.001 * ones(model.dim.w)], [100])

q1_sim = SVector{model.dim.q}([0.0])
q0_sim = SVector{model.dim.q}([0.0])

sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	d = d,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim)
anim = visualize_robot!(vis, model, sim.traj, sample = 1)

control_switch(model, sim.traj.q[end])
B_func(model, q0_sim)
sim.traj.u[end]
# plot(hcat(sim.traj.u...)', linetype = :steppost)
