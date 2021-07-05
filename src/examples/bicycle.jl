include(joinpath(@__DIR__, "..", "dynamics", "bicycle", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# qet hopper simulation
s = get_simulation("bicycle", "flat_2D_lc", "flat")
model = s.model
env = s.env
nq = s.model.dim.q
nu = s.model.dim.u
nc = s.model.dim.c

H = 400
h = 0.01
κ = 1.0e-8

# Design open-loop control trajectory
q0_ref = SVector{nq,T}([0.10, 0.75, 0.10, 0.50, 0.50, 0.00, 0.00])
q1_ref = SVector{nq,T}([0.10, 0.75, 0.10, 0.50, 0.50, 0.00, 0.00])

u_ref = [[0.0, 0.0, 1.0, 0.0] for k = 1:100]
push!(u_ref,  [[0.0, 0.0, 0.0, 0.0] for k = 1:H-100]...);

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

plot(hcat(Vector.([q[6:6] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot(hcat(Vector.([q[7:7] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot!(hcat(Vector.([u[1:2] for u in sim.traj.u])...)', color=:red, linewidth=1.0)

# kinematics_1(model, q0_ref, body = :front)
# kinematics_1(model, q0_ref, body = :rear)
#
# kinematics_2(model, q0_ref, body = :front_hub)
# kinematics_2(model, q0_ref, body = :rear_hub)
#
# kinematics_3(model, q0_ref, body = :front_contact)
# kinematics_3(model, q0_ref, body = :rear_contact)


filename = "bicycle_slide"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)




vis = Visualizer()
open(vis)

filename = joinpath(module_dir(), "src", "dynamics", "bicycle", "mesh", "untitled.obj")
supra_obj = MeshFileObject(filename)
# supra_obj = MeshFileGeometry(filename)

i = 4
setobject!(vis["supra$i"], supra_obj)
settransform!(vis["supra$i"], compose(Translation(0,0,0.1), LinearMap(RotX(π/2))))
default_background!(vis)
