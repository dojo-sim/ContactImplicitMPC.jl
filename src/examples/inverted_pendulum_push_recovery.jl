model = ContactControl.get_model("inverted_pendulum")

# time
h = 0.001
H = 1000

# reference trajectory
ref_traj = contact_trajectory(H, h, model)
ref_traj.h
qref = [0.0; 0.0]
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


#
# ϕ_func(model, q0)
# q0 = rand(nq)
# θ, d1 = q0
# mode = :d1
# k(z) = _kinematics(model, z, mode = mode)
# norm(ForwardDiff.jacobian(k, q0) - _jacobian(model, q0, mode = mode))
#
# ForwardDiff.jacobian(k, q0)
# _jacobian(model, q0, mode = mode)

# initial conditions
q0 = @SVector [0.0 * π, 0.0]
q1 = @SVector [0.0 * π, 0.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, H,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-5),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
status = ContactControl.simulate!(sim)
@test status

# plot(hcat(sim.traj.γ...)')
# plot(hcat(sim.traj.b...)')
# plot(hcat(sim.traj.u...)')

# visualize
include(joinpath(@__DIR__, "..", "dynamics", "inverted_pendulum", "visuals.jl"))
vis = Visualizer()
render(vis)
add_walls!(vis, model)
anim = visualize_robot!(vis, model, sim.traj, sample = 1)
# anim = visualize_force!(vis, model, sim.traj, anim=anim, h=h_sim)

# MPC
N_sample = 2
H_mpc = 10
h_sim = h / N_sample
H_sim = 5000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingVelocityObjective(H_mpc, model.dim,
    q = [Diagonal(1.0 * [10.0; 1.0]) for t = 1:H_mpc],
	v = [Diagonal(0.1 * [0.0001; 0.001] ./ (h^2.0)) for t = 1:H_mpc],
    u = [Diagonal(1.0 * [1.0; 1.0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        solver = :ldl_solver,
        max_iter = 10),
    mpc_opts = LinearizedMPCOptions())

d = impulse_disturbances([[-1.0; 0.0]], [10])

q1_sim = SVector{model.dim.q}([0.0, 0.0])
q0_sim = SVector{model.dim.q}([0.0, 0.0])

sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	d = d,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim)
# vis = Visualizer()
# render(vis)
# open(vis)
# add_walls!(vis, model)
anim = visualize_robot!(vis, model, sim.traj, sample = 1)

γ_max = maximum(hcat(sim.traj.γ...))
u_max = maximum(hcat(sim.traj.u...))

plot((hcat(sim.traj.γ...)./γ_max)', linetype = :steppost)
plot!((hcat(sim.traj.u...)./u_max)[2:2, :]', linetype = :steppost)
plot((hcat(sim.traj.u...)./u_max)[1:1, :]', linetype = :steppost)
plot(hcat(sim.traj.q...)')
