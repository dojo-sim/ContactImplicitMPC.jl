const ContactImplicitMPC = Main
include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
s_sim = get_simulation("hopper_2D", "sine2_2D_lc", "sinusoidal")
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model
env = s.env

# get trajectory
ref_traj = deepcopy(ContactImplicitMPC.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_forward.jld2"),
    load_type = :joint_traj))
# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 100*H*N_sample #500*H*N_sample


# Select the initial speed
v0 = 0.2
Tstance = 0.13 # measure using hop-in-place gait
Tflight = 0.62 # measure using hop-in-place gait
p = raibert_policy(s_sim.model, v0=v0, Tstance=Tstance, Tflight=Tflight, h=h)

off0 = SVector{model.dim.q,T}([0.0, 0.5, 0.0, 0.0])
off1 = SVector{model.dim.q,T}([0*v0*h_sim, 0.5, 0.0, 0.0])
q_ref = SVector{model.dim.q,T}([0.0, 0.5, 0.0, 0.5])
q0_sim = copy(q_ref) + off0
q1_sim = copy(q_ref) + off1

sim = ContactImplicitMPC.simulator(s_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactImplicitMPC.InteriorPointOptions(
		γ_reg = 0.0,
		undercut = Inf,
		r_tol = 1.0e-8,
		κ_tol = 1.0e-8,),
    sim_opts = ContactImplicitMPC.SimulatorOptions(warmstart = true))


status = ContactImplicitMPC.simulate!(sim, verbose = true)


plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:model.dim.u] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:model.dim.c] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

plot_lines!(vis, model, sim.traj.q[1:10:end])
plot_surface!(vis, s_sim.env, n=200, xlims = [-1, 40])
anim = visualize_robot!(vis, model, sim.traj, sample=5)
anim = visualize_force!(vis, model, s_sim.env, sim.traj, anim=anim, h=h_sim, sample = 5)

plot(hcat([u[[1,2]] for u in Vector.(sim.traj.u)[1:end]]...)')


# filename = "hopper_ipopt"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)

# const ContactImplicitMPC = Main
