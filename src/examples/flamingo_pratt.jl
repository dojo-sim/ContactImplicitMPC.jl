include(joinpath(@__DIR__, "..", "dynamics", "flamingo", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

s = get_simulation("flamingo", "flat_2D_lc", "flat")
model = s.model
env = s.env
const ContactControl = Main
ref_traj_ = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))
ref_traj = deepcopy(ref_traj_)

H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 15
h_sim = h / N_sample
H_sim = 3000 #35000

p = pratt_policy(model, h_sim)

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	ip_opts = InteriorPointOptions(
		undercut = Inf,
		γ_reg = 0.0,
		# verbose = true,
        r_tol = 1.0e-8,
        κ_tol = 1.0e-8),
    sim_opts = SimulatorOptions(warmstart = true),
    )

@time status = simulate!(sim, verbose = true)
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
# anim = visualize_robot!(vis, model, sim.traj, sample=10)
anim = visualize_meshrobot!(vis, model, sim.traj, sample=10)
anim = visualize_force!(vis, model, env, sim.traj, anim=anim, h=h_sim, sample=10)

# filename = "flamingo_100_steps"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
#
#
# settransform!(vis["/Cameras/default"],
# 		compose(Translation(0.0, 0.5, -1.0),LinearMap(RotZ(-pi / 2.0))))
#
# # Ghost
# flamingo_ghost!(vis, sim, x -> 0.0)
#
# # Animation
# flamingo_animation!(vis, sim, x -> 0.0)
