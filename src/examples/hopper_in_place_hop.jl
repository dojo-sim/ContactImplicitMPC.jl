include(joinpath(module_dir(), "src", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# qet hopper simulation
s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
model = s.model
env = s.env
nq = s.model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = nc * friction_dim(env)
nz = num_var(model, env)
nθ = num_data(model)

H = 92
# h = 0.04
# h = 0.02
h = 0.01
κ = 1.0e-8

# Design open-loop control trajectory
q0_ref = SVector{nq,T}([0.0, 0.50, 0.0, 0.5])
q1_ref = SVector{nq,T}([0.0, 0.50, 0.0, 0.5])

α = 0.0077 # N=46 h=0.01

u_ref = []
push!(u_ref,  [[0.0,  5.0*α*model.g*(model.mb+model.ml)/2] for k=1:2*6]...);
push!(u_ref,  [[0.0, -0.60*α*model.g*(model.mb+model.ml)/2] for k=1:2*10]...);
push!(u_ref,  [[0.0,  0.14*α*model.g*(model.mb+model.ml)/2] for k=1:2*15]...);
push!(u_ref,  [[0.0,  2.19*α*model.g*(model.mb+model.ml)/2] for k=1:H-length(u_ref)]...);

# Simulate
sim = simulator(s, q0_ref, q1_ref, h, H;
    p = open_loop_policy(u_ref),
    d = no_disturbances(s.model),
    ip_opts = InteriorPointOptions(r_tol=1e-9, κ_tol=2κ),
    sim_opts = SimulatorOptions())

simulate!(sim)

plot(hcat(Vector.(sim.traj.u)...)')
plot(hcat(Vector.(sim.traj.q)...)')
visualize_robot!(vis, s.model, sim.traj.q)

for t = 1:10
    @show norm(residual(s.model, s.env, sim.traj.z[t], sim.traj.θ[t], [κ]))
end

# Check ~perfect loop
@show round.(sim.traj.q[end-2] - q0_ref, digits=3)
@show round.(sim.traj.q[end-1] - q1_ref, digits=3)

# Save trajectory
traj = deepcopy(sim.traj)
gait_path = joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "gaits", "gait_in_place.jld2")
@save gait_path traj

# Reload trajectory
res = JLD2.jldopen(gait_path)
loaded_traj = res["traj"]

traj = get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/hopper_2D/gaits/gait_in_place.jld2"),
    load_type = :joint_traj)
plot(hcat(Vector.(traj.q)...)')
