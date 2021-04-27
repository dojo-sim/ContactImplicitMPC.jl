include(joinpath(@__DIR__, "..", "dynamics", "hopper_3D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
render(vis)

# qet hopper model
model = get_model("hopper_3D")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

H = 92
# h = 0.04
# h = 0.02
h = 0.01
κ = 1.0e-8

# Design open-loop control trajectory
q0_ref = SVector{nq,T}([0.0, 0.0, 0.50, 0.0, 0.0, 0.0, 0.5])
q1_ref = SVector{nq,T}([0.0, 0.0, 0.50, 0.0, 0.0, 0.0, 0.5])
# α = 0.02937 # N=45 h=0.04
# α = 0.01468 # N=45 h=0.02
α = 0.0077 # N=46 h=0.01
# u_ref = []
# push!(u_ref,  [[0.0,  5.0*α*model.g*(model.mb+model.ml)/2] for k=1:6]...);
# push!(u_ref,  [[0.0, -0.9*α*model.g*(model.mb+model.ml)/2] for k=1:8]...);
# push!(u_ref,  [[0.0,  0.2*α*model.g*(model.mb+model.ml)/2] for k=1:14]...);
# push!(u_ref,  [[0.0,  2.1*α*model.g*(model.mb+model.ml)/2] for k=1:H-1-length(u_ref)]...);

u_ref = []
push!(u_ref,  [[0.0, 0.0,   5.0*α*model.g*(model.mb+model.ml)/2] for k=1:2*6]...);
push!(u_ref,  [[0.0, 0.0,  -0.60*α*model.g*(model.mb+model.ml)/2] for k=1:2*10]...);
push!(u_ref,  [[0.0, 0.0,   0.14*α*model.g*(model.mb+model.ml)/2] for k=1:2*15]...);
push!(u_ref,  [[0.0, 0.0,   2.19*α*model.g*(model.mb+model.ml)/2] for k=1:H-length(u_ref)]...);

contact_trajectory(H, h, model)

# Simulate
sim = simulator(model, q0_ref, q1_ref, h, H;
    p = open_loop_policy(u_ref),
    d = no_disturbances(model),
    r! = model.res.r!,
    rz! = model.res.rz!,
    rθ! = model.res.rθ!,
    rz = model.spa.rz_sp,
    rθ = model.spa.rθ_sp,
    ip_opts = InteriorPointOptions(r_tol=1e-9, κ_tol=2κ),
    sim_opts = SimulatorOptions())
simulate!(sim)

# Correct trajectory: set the orientation to 0,0,0
# for t = 1:H+2
#     sim.traj.q[t][[1,2,4,5,6]] .= 0.0
# end
# for t = 1:H
#     sim.traj.b[t] .= 0.0
# end

plot(hcat(Vector.(sim.traj.u)...)')
plot(hcat(Vector.(sim.traj.q)...)'[:,4:6])
anim = visualize_robot!(vis, model, sim.traj)
anim = visualize_force!(vis, model, sim.traj, anim=anim)

# Check ~perfect loop
@show round.(sim.traj.q[end-2] - q0_ref, digits=3)
@show round.(sim.traj.q[end-1] - q1_ref, digits=3)

# Save trajectory
traj = deepcopy(sim.traj)
gait_path = joinpath(@__DIR__, "..", "dynamics", "hopper_3D", "gaits", "gait_in_place.jld2")
@save gait_path traj

# Reload trajectory
res = JLD2.jldopen(gait_path)
loaded_traj = res["traj"]

traj = get_trajectory("hopper_3D", "gait_in_place", load_type=:joint_traj)
plot(hcat(Vector.(traj.q)...)')
