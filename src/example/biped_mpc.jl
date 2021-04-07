include(joinpath(@__DIR__, "..", "dynamics", "biped", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
render(vis)

# get hopper model
model = get_model("biped")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("biped", "gait1", load_type=:split_traj)
ref_traj_copy = deepcopy(ref_traj)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 1
H_mpc = 10
h_sim = h / N_sample
H_sim = 10

# barrier parameter
κ_mpc = 1.0e-4

cost = CostFunction(H_mpc, model.dim,
    q = [Diagonal(1.0 * [1.0, 1.0, 1.0, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]) for t = 1:H_mpc],
    u = [Diagonal(0.1 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, cost,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5))

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-5,
        κ_tol = 1.0e-6),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim)

qq = []
for q in ref_traj_copy.q
    for i = 1:N_sample
        push!(qq, q)
    end
end

L = min(H_sim, length(qq))

plot(hcat(qq...)[1:model.dim.q, 1:L]',
    label = "", color = :black, width = 3.0)
plot!(hcat(sim.traj.q...)[1:model.dim.q, 1:L]',
    label = "", color = :cyan, width = 1.0, legend = :topleft)

visualize!(vis, model, sim.traj.q, Δt = h_sim)#, name = :mpc)
# visualize!(vis, model, ref_traj.q, Δt=10*h/m_opts0.N_sample, name=:mpc)

# filename = "quadruped_mpc_wind"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
ref_traj.q[end-1] - ref_traj.q[1]
