include(joinpath(@__DIR__, "..", "dynamics", "racecar", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# qet hopper simulation
s = get_simulation("racecar", "flat_3D_lc", "flat")
model = s.model
env = s.env
nq = s.model.dim.q
nu = s.model.dim.u
nc = s.model.dim.c

H = 200
h = 0.01
κ = 1.0e-8

# Design open-loop control trajectory
v0 = 3.0
ω0 = v0 / model.rw
#                       X     Y     Z     M1    M2    M3    α     L1    L2    L3    L4    θ1    θ2    θ3    θ4
q0_ref = SVector{nq,T}([0.00, 0.00, 0.35, 0.02, 0.03, 0.02, 0.40, 0.15, 0.15, 0.15, 0.15, 0.00, 0.00, 0.00, 0.00])
q1_ref = SVector{nq,T}([v0*h, 0.00, 0.35, 0.02, 0.03, 0.02, 0.40, 0.15, 0.15, 0.15, 0.15, ω0*h, ω0*h, ω0*h, ω0*h])

u_ref = [[0.0, 0.00, 0.00, 2.50, 2.50] for k = 1:100]
push!(u_ref,  [[0.0, 0.0, 0.0, 0.0, 0.0] for k = 1:H-100]...);

# Simulate
sim = simulator(s, q0_ref, q1_ref, h, H;
    p = open_loop_policy(u_ref),
    d = no_disturbances(s.model),
    ip_opts = InteriorPointOptions(r_tol=1e-9, κ_tol=2κ),
    sim_opts = SimulatorOptions())

simulate!(sim)


plot_surface!(vis, s.env, xlims=[-10,10.], ylims=[-10, 10.])
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


filename = "racecar_drift"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)

q0_ref = SVector{nq,T}([0.10, 0.00, 0.75, 0.00, 0.00, 0.10, 0.10, 0.50, 0.50, 0.50, 0.50, 0.00, 0.00, 0.00, 0.00])

build_robot!(vis, model, name = :test)
set_robot!(vis, model, q0_ref, name = :test)



vis = Visualizer()
open(vis)

filename = joinpath(module_dir(), "src", "dynamics", "racecar", "mesh", "untitled.obj")
supra_obj = MeshFileObject(filename)
supra_obj = MeshFileGeometry(filename)

setobject!(vis["supra$i"], supra_obj)
settransform!(vis["supra$i"], compose(Translation(0,0,0.1), LinearMap(RotX(π/2))))
default_background!(vis)
