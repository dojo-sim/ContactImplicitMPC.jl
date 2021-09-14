const ContactControl = Main
vis = Visualizer()
open(vis)
include(joinpath(@__DIR__, "..", "dynamics", "planarpush_2D", "visuals.jl"))

################################################################################
# Simulate reference trajectory
################################################################################
s = get_simulation("planarpush_2D", "flat_2D_lc", "flat")
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
H = 550
N_sample = 1

# initial conditions
q0 = @SVector [0.0, 0.0, -0.30, 5e-3]
q1 = @SVector [0.0, 0.0, -0.30, 5e-3]

# p = open_loop_policy(fill(SVector{nu}([10.0*h, -0.0]), H*2), N_sample=N_sample)
p = open_loop_policy([SVector{nu}([15*h*(0.5+sin(t/50))]) for t=1:H], N_sample=N_sample)
# p = open_loop_policy(fill(SVector{nu}([40*h,]), H), N_sample=N_sample)
# p = open_loop_policy(fill(SVector{nu}([32*h, -0.0]), H), N_sample=N_sample)
# p = open_loop_policy(fill(SVector{nu}([0*h, -0.0]), H), N_sample=N_sample)

# simulator
sim0 = simulator(s, q0, q1, h, H,
				p = p,
				ip_opts = ContactControl.InteriorPointOptions(
					γ_reg = 0.0,
					undercut = Inf,
					r_tol = 1.0e-8,
					κ_tol = 1.0e-8),
				sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
status = simulate!(sim0)
@test status

anim = visualize_robot!(vis, model, sim0.traj, α=0.7)
anim = visualize_force!(vis, model, env, sim0.traj, anim=anim, )
plot(hcat([x[3:3] for x in sim0.traj.q]...)')
plot(hcat([x[1:nu] for x in sim0.traj.u]...)')
plot(hcat([x[1:nw] for x in sim0.traj.w]...)')
plot(hcat([x[1:nc] for x in sim0.traj.γ]...)')
plot(hcat([x[1:nb] for x in sim0.traj.b]...)')
# plot(hcat([x[1:8] for x in sim0.traj.b]...)')

differentiate_solution!(sim0.ip)
sim0.ip.δz
sim0.ip.δz[iq2, iu]

################################################################################
# MPC Control
################################################################################

ref_traj = deepcopy(sim0.traj)
ref_traj.H

# MPC
N_sample = 5
H_mpc = 5
h_sim = h / N_sample
H_sim = min((H-H_mpc-1)*N_sample, 5500)
# H_sim = 1

# barrier parameter
κ_mpc = 2.0e-4
# κ_mpc = 1.0e-2

# Aggressive
obj = TrackingVelocityObjective(model, s.env, H_mpc,
	q = [Diagonal(1.0e-4 * [1e2, 1e2, 1e-0, 1e-0]) for t = 1:H_mpc],
	v = [Diagonal(1.0e-3 * ones(model.dim.q)) for t = 1:H_mpc],
	u = [Diagonal(1.0e-6 * ones(model.dim.u)) for t = 1:H_mpc],
	γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
	b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:H_mpc])

# #
# obj = TrackingVelocityObjective(model, s.env, H_mpc,
# 	q = [Diagonal(1.0e-0 * [1e0, 1e0, 1e-4, 1e-4]) for t = 1:H_mpc],
# 	v = [Diagonal(1.0e-1 * ones(model.dim.q)) for t = 1:H_mpc],
# 	u = [Diagonal(1.0e-5 * ones(model.dim.u)) for t = 1:H_mpc],
# 	γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
# 	b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
		r_tol = 3e-4,
		# β_init = 1e-8,
		β_init = 1e-5,
		solver = :lu_solver,
		# solver = :ldl_solver,
		verbose=true,
		max_iter = 5),
		# max_iter = 50),
    mpc_opts = LinearizedMPCOptions(
		# live_plotting=true
		))

# p.κ
# p.traj.κ[1]
# p.ref_traj.κ[1]
# p.im_traj.lin[1].κ[1]
# p.newton.traj.κ[1]
# p.newton.traj_cand.κ[1]
# p.im_traj.ip[1].κ[1]
# p.im_traj.ip[1].opts.κ_init[1]
# p.im_traj.ip[1].opts.κ_tol[1]


# p = open_loop_policy(fill(SVector{nu}([40*h, -0.0]), H*2), N_sample=N_sample)
# p = open_loop_policy([SVector{nu}([15*h*(0.5+sin(t/50))]) for t=1:H], N_sample=N_sample)
using Random
Random.seed!(100)
d = open_loop_disturbances([[0.0, -1.0*(0.25*rand()+0.5)] for t=1:H_sim], N_sample)

q0_sim = @SVector [0.0, 0.0, -0.30, 5e-3]
q1_sim = @SVector [0.0, 0.0, -0.30, 5e-3]

sim = ContactControl.simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	d = d,
    ip_opts = ContactControl.InteriorPointOptions(
		γ_reg = 0.0,
		undercut = Inf,
		r_tol = 1.0e-8,
        κ_tol = 1.0e-8,),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

sim.traj.H
sim.traj.u

@time status = ContactControl.simulate!(sim)
anim = visualize_robot!(vis, model, sim.traj, name=:sim, sample = N_sample, α=1.0)
anim = visualize_robot!(vis, model, ref_traj, name=:ref, anim=anim, sample = 1, α=0.5)

plot([t*h for t=0:H-1], hcat([x[1:1] ./ N_sample for x in ref_traj.u[1:H]]...)')
scatter!([t*h/N_sample for t=0:H_sim-1], hcat([x[1:1] for x in sim.traj.u[1:H_sim]]...)')

iu = Vector(2nq .+ (1:nu))
iq2 = Vector(1:nq)
differentiate_solution!(p.im_traj.ip[end])
plot(log.(10, [p.im_traj.ip[t].δz[iq2, iu][1] for t = 1:H]))
plot(log.(10, [p.im_traj.ip[t].δz[iq2, iu][3] for t = 1:H]))
p.im_traj.ip[H].δzs[iq2, iu]


num_var(model, env)
num_data(model)

r = p.im_traj.ip[1].r
rz = p.im_traj.ip[1].rz
rθ = p.im_traj.ip[1].rθ

size(rz.Dx)
size(rz.Dy1)
size(rz.Rx)
size(rz.Ry1)
size(Diagonal(rz.Ry2))
size(Diagonal(rz.y2))
size(Diagonal(rz.y1))

rz = [rz.Dx rz.Dy1 zeros(nq,8);
	  rz.Rx rz.Ry1 Diagonal(rz.Ry2);
	  zeros(8,nq) Diagonal(rz.y2) Diagonal(rz.y1)]
rθ = rθ.rθ0

vidyn = Vector(p.im_traj.ip[1].rz.idyn)
virst = Vector(p.im_traj.ip[1].rz.irst)
vibil = Vector(p.im_traj.ip[1].rz.ibil)
vix   = Vector(p.im_traj.ip[1].rz.ix)
viy1  = Vector(p.im_traj.ip[1].rz.iy1)
viy2  = Vector(p.im_traj.ip[1].rz.iy2)

cond(rz)
cond(rθ)
rank(rz)
rank(rθ)

plot(Gray.(1e10*abs.(rz[[vidyn; virst; vibil;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz[[vidyn;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz[[virst;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz[[vibil;], [vix; viy1; viy2]])))

plot(Gray.(abs.(rz[[vidyn; virst; vibil;], [vix; viy1; viy2]])))
plot(Gray.(abs.(rz[[virst;], [viy1; ]])))



# filename = "planarpush_precise"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
# why kappa is not 1e-4 in the ip solver ?
