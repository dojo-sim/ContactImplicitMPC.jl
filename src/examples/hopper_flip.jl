include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# qet hopper model
# model = get_model("hopper_2D")
model = get_model("hopper_2D", surf="stairs")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

H = 119
h = 0.024337
κ = 1.0e-8

# Design open-loop control trajectory
q0_ref = SVector{nq,T}([0.75, 1.25, 0.0, 0.5])
q1_ref = SVector{nq,T}([0.75, 1.25, 0.0, 0.5])

ref_traj_ = get_trajectory("hopper_2D", "gait_forward_high", load_type=:joint_traj)
ref_traj = deepcopy(ref_traj_)
u_for = deepcopy(ref_traj.u)

α = 0.016
u_ref = []
push!(u_ref,  [[0.0,  4.00*α*model.g*(model.mb+model.ml)/2] for k=1:9]...);
push!(u_ref,  [[0.0, -0.80*α*model.g*(model.mb+model.ml)/2] for k=1:5]...);
push!(u_ref,  [[0.0, -0.55*α*model.g*(model.mb+model.ml)/2] for k=1:7]...);
push!(u_ref,  [[0.0,  1.20*α*model.g*(model.mb+model.ml)/2] for k=1:16]...);
push!(u_ref,  [[0.0,  1.65*α*model.g*(model.mb+model.ml)/2] for k=1:H-length(u_ref)]...);

for t = 1:H
    u_ref[t][1] = 0.1* u_for[t][1]
end

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
plot(hcat(Vector.(sim.traj.u)...)')
plot(hcat(Vector.(sim.traj.q)...)')
visualize_robot!(vis, model, sim.traj.q)

for t = 1:10
    @show norm(residual(model, sim.traj.z[t], sim.traj.θ[t], [κ]))
end

plot_surface!(vis, model.env, n=400)
# Check ~perfect loop
@show round.(sim.traj.q[end-2] - q0_ref, digits=3)
@show round.(sim.traj.q[end-1] - q1_ref, digits=3)

# Save trajectory
traj = deepcopy(sim.traj)
gait_path = joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "gaits", "gait_in_place_high.jld2")
@save gait_path traj

# Reload trajectory
res = JLD2.jldopen(gait_path)
loaded_traj = res["traj"]

traj = get_trajectory("hopper_2D", "gait_in_place_high", load_type=:joint_traj)
plot(hcat(Vector.(traj.q)...)')
