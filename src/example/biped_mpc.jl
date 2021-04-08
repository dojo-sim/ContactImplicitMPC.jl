include(joinpath(@__DIR__, "..", "dynamics", "biped", "visuals.jl"))
T = Float64

# get hopper model
model = get_model("biped")

# get trajectory
ref_traj = get_trajectory("biped", "gait13", load_type=:split_traj)
ref_traj_copy = deepcopy(ref_traj)

plot(hcat(ref_traj.u...)', linetype = :steppost)
plot(hcat(ref_traj.γ...)', linetype = :steppost)
plot(hcat([γ[1] for γ in ref_traj.γ]...)', linetype = :steppost)
plot!(hcat([sum(2.0 * b[1:2]) for b in ref_traj.b]...)', linetype = :steppost)
plot!(hcat([γ[2] for γ in ref_traj.γ]...)', linetype = :steppost)
plot!(hcat([sum(2.0 * b[3:4]) for b in ref_traj.b]...)', linetype = :steppost)
# plot(hcat([sum(b[1:2]) for b in ref_traj.b]...)', linetype = :steppost)
plot(hcat([γ[3] for γ in ref_traj.γ]...)', linetype = :steppost)
plot!(hcat([sum(2.0 * b[5:6]) for b in ref_traj.b]...)', linetype = :steppost)
plot!(hcat([γ[4] for γ in ref_traj.γ]...)', linetype = :steppost)
plot!(hcat([sum(2.0 * b[7:8]) for b in ref_traj.b]...)', linetype = :steppost)


ref_traj.h
for t = 1:ref_traj.H
	r = residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)
	@show norm(r)
end

@show maximum([norm(ContactControl.dynamics(model,
	ref_traj.h, ref_traj.q[t], ref_traj.q[t+1], ref_traj.u[t],
	zeros(model.dim.w), ref_traj.γ[t], ref_traj.b[t], ref_traj.q[t+2]), Inf) for t = 1:ref_traj.H])

# initial conditions
q0 = SVector{model.dim.q}(ref_traj.q[1])
q1 = SVector{model.dim.q}(ref_traj.q[2])

# simulator
H_sim = 5
sim = ContactControl.simulator(model, q0, q1, 1.0 * ref_traj.h, H_sim,
    p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u], N_sample = 1),
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
visualize!(vis, model, [ref_traj.q[1]], Δt = h) #, name = :mpc)

# mpc
N_sample = 1
H_mpc = 10
h_sim = h / N_sample
H_sim = 10

# barrier parameter
κ_mpc = 1.0e-4

cost = CostFunction(H_mpc, model.dim,
    q = [Diagonal(1.0e-2 * [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) for t = 1:H_mpc],
    u = [Diagonal(1.0 * ones(model.dim.u)) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-5 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-5 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, cost,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5))

B_func(model, rand(model.dim.q))
B_fast(model, rand(model.dim.q))

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

plot(hcat(qq...)[3:3, 1:L]',
    label = "", color = :black, width = 3.0)
plot!(hcat(sim.traj.q...)[3:3, 1:L]',
    label = "", color = :cyan, width = 1.0, legend = :topleft)


vis = Visualizer()
# open(vis)
render(vis)
visualize!(vis, model, sim.traj.q, Δt = h_sim) #, name = :mpc)
# visualize!(vis, model, ref_traj.q, Δt=10*h/m_opts0.N_sample, name=:mpc)

# filename = "quadruped_mpc_wind"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
