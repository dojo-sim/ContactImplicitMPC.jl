include(joinpath(@__DIR__, "..", "dynamics", "unicycle", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# qet hopper simulation
s = get_simulation("unicycle", "flat_2D_lc", "flat")
model = s.model
env = s.env
nq = s.model.dim.q
nu = s.model.dim.u
nc = s.model.dim.c

H = 92
h = 0.008
κ = 1.0e-8

# Design open-loop control trajectory
q0_ref = SVector{nq,T}([0.10, 0.75, 0.50, 0.00])
q1_ref = SVector{nq,T}([0.10, 0.75, 0.50, 0.30])

u_ref = [[0.0, 10] for k = 1:25]
push!(u_ref,  [[0.0, -25.0] for k = 1:H-25]...);

# Simulate
sim = simulator(s, q0_ref, q1_ref, h, H;
    p = open_loop_policy(u_ref),
    d = no_disturbances(s.model),
    ip_opts = InteriorPointOptions(r_tol=1e-9, κ_tol=2κ),
    sim_opts = SimulatorOptions())

simulate!(sim)


plot_surface!(vis, s.env, xlims=[-3,3.])
anim = visualize_robot!(vis, s.model, sim.traj)
anim = visualize_force!(vis, s.model, s.env, sim.traj, anim = anim)

plot(hcat(Vector.([q[4:4] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot!(hcat(Vector.([u[1:2] for u in sim.traj.u])...)', color=:red, linewidth=1.0)



filename = "unicycle_slide"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)
