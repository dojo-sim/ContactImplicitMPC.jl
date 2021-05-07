include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# get hopper model
model_sim = get_model("hopper_2D", surf="stairs")
model = get_model("hopper_2D", surf="stairs")
# model = get_model("hopper_2D")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj0 = get_trajectory("hopper_2D", "gait_forward_high", load_type=:joint_traj)
ref_traj = deepcopy(ref_traj0)
for t = 1:ref_traj.H+2
    ref_traj.q[t][1] += 0.25
end

LinearizedStep(model, ref_traj.z[end], ref_traj.θ[end], κ)
LinearizedStep(model, ref_traj.z[end], ref_traj.θ[end], κ)

nz = num_var(model)
nθ = num_data(model)
z0 = SizedVector{nz,T}(ref_traj.z[end])
θ0 = SizedVector{nθ,T}(ref_traj.θ[end])
κ0 = 1e-8
r0 = zeros(SizedVector{nz,T})
model.res.r!(r0, z0, θ0, κ0)
r0

ϕ_fast(model, ref_traj.q[end])
ϕ_func(model, ref_traj.q[end])

H = ref_traj.H
h = ref_traj.h
κ = 1.0e-8

# Cost function
obj = TrackingObjective(H, model.dim,
    q = [[Diagonal(1.0e-5 * [0.01, 1,1,1]) for t = 1:H-1]; [Diagonal(1.0e-4 * ones(nq))]],
    u = [Diagonal(1.0e-5 * [0.1, 1.0]) for t = 1:H],
    γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:H],
    b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:H])
n_opts = NewtonOptions(r_tol=3e-8, max_iter=100)
im_traj = ImplicitTraj(ref_traj, model_sim, κ=κ)
core = Newton(H, h, model_sim, ref_traj, im_traj, obj = obj, opts = n_opts)

ϕ_fast(model_sim, ref_traj.q[end])

q0_dist = deepcopy(ref_traj.q[1] + [-0.00,0.0,0,0])
q1_dist = deepcopy(ref_traj.q[2] + [-0.00,0.0,0,0])
newton_solve!(core, model_sim, im_traj, ref_traj, q0=q0_dist, q1=q1_dist, verbose=true)

plot_surface!(vis, model_sim.env, n=200)
visualize_robot!(vis, model, core.traj.q)
visualize_robot!(vis, model, ref_traj.q)

plot(hcat(Vector.(ref_traj.u)...)')
plot(hcat(Vector.(core.traj.u)...)')
plot(hcat(Vector.(ref_traj.q)...)')
plot(hcat(Vector.(core.traj.q)...)')
# Check ~perfect loop
@show round.(core.traj.q[end-2] - ref_traj.q[1], digits=3)
@show round.(core.traj.q[end-1] - ref_traj.q[2], digits=3)

# Save trajectory
traj = deepcopy(core.traj)
gait_path = joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "gaits", "gait_stair.jld2")
@save gait_path traj

# Reload trajectory
res = JLD2.jldopen(gait_path)
loaded_traj = res["traj"]

traj = get_trajectory("hopper_2D", "gait_stair", load_type=:joint_traj)
plot(hcat(Vector.(traj.q)...)')
