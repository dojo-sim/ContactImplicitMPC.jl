include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# qet hopper model
model = get_model("hopper_2D")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

H = 119
h = 0.02377
κ = 1.0e-8

# Design open-loop control trajectory
q0_ref = SVector{nq,T}([0.0, 0.50, 0.0, 0.5])
q1_ref = SVector{nq,T}([0.0, 0.50, 0.0, 0.5])

α = 0.02
u_ref = []
push!(u_ref,  [[0.0,  5.00*α*model.g*(model.mb+model.ml)/2] for k=1:8]...);
push!(u_ref,  [[0.0, -1.00*α*model.g*(model.mb+model.ml)/2] for k=1:10]...);
push!(u_ref,  [[0.0, -0.55*α*model.g*(model.mb+model.ml)/2] for k=1:15]...);
push!(u_ref,  [[0.0,  0.14*α*model.g*(model.mb+model.ml)/2] for k=1:35]...);
push!(u_ref,  [[0.0,  2.105*α*model.g*(model.mb+model.ml)/2] for k=1:H-length(u_ref)]...);

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
visualize_robot!(vis, model, sim.traj)

plot_surface!(vis, model.env)
for t = 1:10
    @show norm(residual(model, sim.traj.z[t], sim.traj.θ[t], [κ]))
end

# Check ~perfect loop
@show round.(sim.traj.q[end-2] - q0_ref, digits=3)
@show round.(sim.traj.q[end-1] - q1_ref, digits=3)

# Save trajectory
traj = deepcopy(sim.traj)
gait_path = joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "gaits", "gait_in_place_flip.jld2")
@save gait_path traj

# Reload trajectory
res = JLD2.jldopen(gait_path)
loaded_traj = res["traj"]

traj = get_trajectory("hopper_2D", "gait_in_place_flip", load_type=:joint_traj)
plot(hcat(Vector.(traj.q)...)')
