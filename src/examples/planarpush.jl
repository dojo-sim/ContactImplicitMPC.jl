const ContactControl = Main
vis = Visualizer()
open(vis)
# render(vis)



include(joinpath(@__DIR__, "..", "dynamics", "planarpush", "visuals.jl"))
model = get_model("planarpush")
nq = model.dim.q
nu = model.dim.u
nw = model.dim.w
nc = model.dim.c
nb = model.dim.b
nf = 4

# time
# h = 0.04
h = 0.002
H = 600
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 3000

# # reference trajectory
# ref_traj = contact_trajectory(H, h, model)
# ref_traj.h
# qref = [0.0; 0.0]
# ur = zeros(model.dim.u)
# γr = zeros(model.dim.c)
# br = zeros(model.dim.b)
# ψr = zeros(model.dim.c)
# ηr = zeros(model.dim.b)
# wr = zeros(model.dim.w)
#
# # set reference
# for t = 1:H
# 	ref_traj.z[t] = pack_z(model, qref, γr, br, ψr, ηr)
# 	ref_traj.θ[t] = pack_θ(model, qref, qref, ur, wr, model.μ_world, ref_traj.h)
# end
#
# # test reference
# for t = 1:H
# 	r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)#[model.dim.q .+ (1:model.dim.c)]
# 	@test norm(r) < 1.0e-4
# end

# initial conditions
q0 = @SVector [0.201, 0.20, 0.03, 0.0, -0.150, 0.21, 0.00]
q1 = @SVector [0.200, 0.20, 0.03, 0.0, -0.149, 0.21, 0.00]
q0 = @SVector [0.20, 0.20, 0.0, 0.0, -0.15, 0.28, 0.00]
q1 = @SVector [0.20, 0.20, 0.0, 0.0, -0.15, 0.28, 0.00]
# q2 = @SVector [0.02, 0.02, 0.3, 0.0, 0.11, 0.21, 0.6]


# C_fast(model, q0, (q1-q0)/h_sim)
#
# #
# h0 = [h_sim]
# d0 = zeros(nq)
# # q0 = zeros(nq)
# # q1 = zeros(nq)
# u1 = zeros(nu)
# w1 = zeros(nw)
# λ1 = zeros(nc+nc*nf)
# # q2 = zeros(nq)
#
# q0 = [0.0, 0.0, 0.3, 0.0, 0.1, 0.2, 0.6]
# q1 = [0.02, 0.02, 0.3, 0.0, 0.1, 0.2, 0.6]
# q2 = [0.0, 0.0, 0.3, 0.0, 0.1, 0.2, 0.6]
#
# d0 = model.dyn.d(h0, q0, q1, u1, w1, λ1, q2)
#
#
# a = 10
# a = 10
# a = 10
# model.
#
# nc = model.dim.c
# nb = model.dim.b
# nf = Int(nb / nc)
# np = dim(model.env)
#
# q0, q1, u1, w1, μ, h = unpack_θ(model, θ)
# q2, γ1, b1, ψ1, η1, s1, s2 = unpack_z(model, z)
#
# ϕ = ϕ_func(model, q2)
#
# k = kinematics(model, q2)
# λ1 = contact_forces(model, γ1, b1, q2, k)
# vT_stack = velocity_stack(model, q1, q2, k, h)
# ψ_stack = transpose(E_func(model)) * ψ1
#
# [model.dyn.d(h, q0, q1, u1, w1, λ1, q2);
#  s1 - ϕ;
#  vT_stack + ψ_stack - η1;
#  s2 .- (μ[1] * γ1 .- E_func(model) * b1);
#  γ1 .* s1 .- κ;
#  b1 .* η1 .- κ;
#  ψ1 .* s2 .- κ]


#
# model.con.ϕ(q0)
# model.con.ϕ(q0)
#
# r0 = rand(num_var(model))
# z0 = rand(num_var(model))
# θ0 = rand(num_data(model))
# κ0 = [1e-4]
# model.res.r!(r0, z0, θ0, κ0)
# model.res.r!(r0, z0, θ0, κ0)

# p = open_loop_policy(fill(SVector{nu}([0.0, -0.125]), H_sim), N_sample=N_sample)
# p = open_loop_policy(fill(SVector{nu}([0.01, -0.0]), H_sim*10), N_sample=N_sample)
p = open_loop_policy(fill(SVector{nu}([0.033, -0.0]), H_sim*10), N_sample=N_sample)

# simulator
sim = ContactControl.simulator(model, q0, q1, h, H,
	p = p,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_init=1e-6, κ_tol = 2.0e-6),
	sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
status = ContactControl.simulate!(sim)
@test status

visualize_robot!(vis, model, sim.traj.q[1:2])
visualize_robot!(vis, model, sim.traj)

plot(hcat([q[1:7] for q in sim.traj.q]...)')
sim.traj.q
[q[1] for q in sim.traj.q]



# MPC
N_sample = 2
H_mpc = 40
h_sim = h / N_sample
H_sim = 1000

# barrier parameter
κ_mpc = 1.0e-4

# SLOW RECOVERY
obj = TrackingVelocityObjective(H_mpc, model.dim,
	q = [Diagonal([12*(t/H_mpc)^2; 2.0*(t/H_mpc)^4]) for t = 1:H_mpc-0],
	v = [Diagonal([1; 0.01] ./ (h^2.0)) for t = 1:H_mpc-0],
	u = [Diagonal([100; 1]) for t = 1:H_mpc-0],
	γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc-0],
	b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])
# FAST RECOVERY
obj = TrackingVelocityObjective(H_mpc, model.dim,
    q = [Diagonal([12*(t/H_mpc)^2; 12*(t/H_mpc)^2]) for t = 1:H_mpc-0],
	v = [Diagonal([1; 0.01] ./ (h^2.0)) for t = 1:H_mpc-0],
    u = [Diagonal([100; 1]) for t = 1:H_mpc-0],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc-0],
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

idx_d1 = 20
idx_d2 = idx_d1 + 200
idx_d3 = idx_d2 + 80
idx_d4 = idx_d3 + 200
idx_d5 = idx_d4 + 30
idx = [idx_d1, idx_d2, idx_d3, idx_d4, idx_d5]
impulses = [[-5.5; 0.0], [+5.5; 0.0], [+5.5; 0.0], [-1.5; 0.0], [-4.5; 0.0]]
d = impulse_disturbances(impulses, idx)


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

sample = 1
anim = visualize_robot!(vis, model, sim.traj, sample = sample, α=0.5)
pθ_right = generate_pusher_traj(d, sim.traj, side=:right)
pθ_left  = generate_pusher_traj(d, sim.traj, side=:left)
visualize_disturbance!(vis, model, pθ_right, anim=anim, sample=sample, offset=0.05, name=:PusherRight)
visualize_disturbance!(vis, model, pθ_left,  anim=anim, sample=sample, offset=0.05, name=:PusherLeft)


γ_max = maximum(hcat(sim.traj.γ...))
u_max = maximum(hcat(sim.traj.u...))

plot((hcat(sim.traj.γ...) ./ γ_max)', linetype = :steppost)
plot!((hcat(sim.traj.u...) ./ u_max)[2:2, :]', linetype = :steppost)
plot((hcat(sim.traj.u...) ./ u_max)[1:1, :]', linetype = :steppost)
plot(hcat(sim.traj.q...)')

plot([pθ[1] for pθ in pθ_right])
plot!([pθ[1] for pθ in pθ_left])
plot!([q[1] for q in sim.traj.q[3:end]])



# Display highlights
t_highlights = [20,23,35,216]
α_highlights = ones(4)
α_pusher = [1.0, 0.0, 0.0, 0.0]
