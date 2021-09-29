# run hopper_in_place_hop.jl

include(joinpath(@__DIR__, "..", "dynamics", "hopper_3D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
# model_sim = get_model("hopper_2D", surf="sinusoidal")
model_sim = get_model("hopper_3D", surf="flat")
nq = model_sim.dim.q
nu = model_sim.dim.u
nc = model_sim.dim.c
nb = model_sim.dim.b
nw = model_sim.dim.w

# time
H = 92#traj.H
h = 0.01#traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 5000

# Select the initial speed
v0 = [0.0; 0.2]
Tstance = 0.13 # measure using hop-in-place gait
Tflight = 0.62 # measure using hop-in-place gait
include(joinpath(module_dir(), "src/controller/raibert_3D_policy.jl"))
p = raibert_policy(model_sim, v0=v0, Tstance=Tstance, Tflight=Tflight, h=h)

off0 = SVector{nq,T}([0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5])
off1 = SVector{nq,T}([v0[1]*h_sim, v0[2]*h_sim, 0.5, 0.0, 0.0, 0.0, 0.5])
q_ref = SVector{nq,T}([0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0])
q0_sim = copy(q_ref) + off0
q1_sim = copy(q_ref) + off1

# Disturbances
w = [zeros(nw) for t=1:Int(ceil(H_sim/N_sample))]
# w[75] += [5.0, 5.0, 0.0]
d = open_loop_disturbances(w)

sim = simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    d = d,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    )

@time status = simulate!(sim)


# plt = plot(layout=(3,1), legend=false)
# plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
# plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)
#
# plot_lines!(vis, model_sim, sim.traj.q[1:10:end])
# plot_surface!(vis, model_sim.env, n=200)
anim = visualize_robot!(vis, model_sim, sim.traj, sample=20)

# hopper_3D_euler = Hopper3D(Dimensions(nq, nu, nw, nc, nb),
# 		mb, ml, Jb, Jl,
# 		μ_world, μ_joint, g,
# 		:RotXYX, # Euler angles
# 		BaseMethods(), DynamicsMethods(), ContactMethods(),
# 		ResidualMethods(), ResidualMethods(),
# 		SparseStructure(spzeros(0, 0), spzeros(0, 0)),
# 		SVector{7}(zeros(7)),
# 		environment_3D(x -> 0.075 * sin(2π * x[1])),
# 	    )
#
# eul = RotXYX(0.5 * π, 0.0, 0.0)
# qq = Array(copy(q_ref))
# qq[4:6] = [eul.theta1; eul.theta2; eul.theta3]
#
# build_robot!(vis, hopper_3D_euler, name = model_name(hopper_3D_euler))
# set_robot!(vis, model_sim, q1_sim, name = model_name(model_sim))
#
# anim = visualize_force!(vis, model_sim, sim.traj, anim=anim, h=h_sim, sample=20)
#
# plot(hcat([q[[1,2,3,4]] for q in Vector.(sim.traj.q)[1:1200]]...)')
# plot(hcat([q[[3]] for q in Vector.(sim.traj.q)[1:end]]...)')
# plot(hcat([γ[[1]] for γ in Vector.(sim.traj.γ)[1:end]]...)')
# plot(hcat([u[[1,2]] for u in Vector.(sim.traj.u)[1:end]]...)')

# filename = "raibert_hopper_push"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
