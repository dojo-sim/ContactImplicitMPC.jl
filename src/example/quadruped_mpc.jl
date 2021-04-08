include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
T = Float64
vis = Visualizer()
render(vis)

# get hopper model
model = get_model("quadruped")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("quadruped", "gait1", load_type=:split_traj)
H = ref_traj.H
h = ref_traj.h
κ = 1.0e-4

# initial conditions
q0 = SVector{model.dim.q}(ref_traj.q[1])
q1 = SVector{model.dim.q}(ref_traj.q[1])

# simulator
H_sim = 25
sim = ContactControl.simulator(model, q0, q1, 1.0 * h, H_sim,
    p = no_policy(model), #ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u], N_sample = 1),
    ip_opts = ContactControl.InteriorPointOptions(r_tol = 1.0e-8, κ_init = 1.0e-5, κ_tol = 1.0e-6),
    sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot(hcat(ref_traj.q...)[1:model.dim.q, 1:H_sim]',
    label = ["x" "y" "z"], color = :black, width = 3.0)
plot!(hcat(sim.traj.q...)[1:model.dim.q, 1:H_sim]',
    label = ["x" "y" "z"], color = :red, width = 1.0, legend = :topleft)

vis = Visualizer()
# open(vis)
render(vis)
visualize!(vis, model, sim.traj.q, Δt = h) #, name = :mpc)


ref_traj0 = deepcopy(ref_traj)
n_opts0 = NewtonOptions(r_tol=3e-4, κ_init=κ, κ_tol=2κ, solver_inner_iter=5)
m_opts0 = MPCOptions{T}(
            N_sample=2,
            M=100,
            H_mpc=10,
            κ=κ,
            κ_sim=1e-8,
            r_tol_sim=1e-8,
            open_loop_mpc=false,
            w_amp=0.0*[-0.10, -0.90],
            ip_max_time=0.1,
            live_plotting=false)
cost0 = CostFunction(H, model.dim,
    q = [Diagonal(1e-2 * [0.02,0.02,1,.15,.15,.15,.15,.15,.15,.15,.15,]) for t = 1:m_opts0.H_mpc],
    u = [Diagonal(3e-2 * ones(nu)) for t = 1:m_opts0.H_mpc],
    γ = [Diagonal(1.0e-100 * ones(nc)) for t = 1:m_opts0.H_mpc],
    b = [Diagonal(1.0e-100 * ones(nb)) for t = 1:m_opts0.H_mpc])
core0 = Newton(m_opts0.H_mpc, h, model, cost=cost0, opts=n_opts0)
mpc0 = MPC(model, ref_traj0, m_opts=m_opts0)
@time dummy_mpc(model, core0, mpc0, verbose=true)
# @profiler dummy_mpc(model, core0, mpc0, verbose=true)

# mpc0.impl.ip[1].solver

# mpc0.impl
# lin0 = mpc0.impl.lin[1]
# r  = RLin(model, lin0.z, lin0.θ, lin0.r, lin0.rz, lin0.rθ)
# rz = RZLin(model, lin0.rz)
# rθ = RθLin(model, lin0.rθ)
#
#
# 2.9/(100*h)

using Plots
plt = plot(layout=(2,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], m_opts0.N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.(mpc0.q_sim)...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], m_opts0.N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in mpc0.u_sim]*m_opts0.N_sample)...)', color=:blue, linewidth=1.0)

visualize!(vis, model, mpc0.q_sim[1:10:end], Δt=10*h/m_opts0.N_sample, name=:mpc)

filename = "quadruped_mpc_downwind"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)
11/(100*h)
