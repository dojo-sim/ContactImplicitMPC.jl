model = get_model("particle", surf = "quadratic")

# time
h = 0.01
T = 1000

## DROP
# initial conditions
q1 = @SVector [1.0, 0.5, 2.0]
q0 = @SVector [1.1, 0.5, 2.0]

# simulator
sim = ContactControl.simulator(model, q0, q1, h, T,
	ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-6, κ_tol = 1.0e-6),
	sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
@time status = ContactControl.simulate!(sim)
@test status

plot(hcat(sim.traj.q[1:3:T]...)', label = ["x" "y" "z"], legend = :bottomleft)
@show ϕ_func(model, sim.traj.q[end])
@show model.env.surf(sim.traj.q[end])
@show sim.traj.q[end]

include(joinpath(pwd(), "src/dynamics/particle/visuals.jl"))

vis = Visualizer()
render(vis)
visualize!(vis, model, sim.traj.q,
	Δt = h, r = 0.1)
plot_surface!(vis, model.env)

using Plots
N = 100
t = range(-2.0, stop = 2.0, length = N)
s = t.^2.0
plot(t, s, ratio = :equal, color = :black, width = 2.0, label = "")

nw = [0.0; 0.0; 1.0]
tw = [1.0; 0.0; 0.0]
env = Environment{R2}(x-> x^2.0, x -> )
t_cand = -1.0
s_cand = t_cand.^2.0
p_cand = [t_cand; 0.0; s_cand]
tan_norm = rotation(model.env, p_cand)' * tw
norm_norm = rotation(model.env, p_cand)' * nw

plot!([p_cand[1], p_cand[1] + tan_norm[1]],
	[p_cand[3], p_cand[3] + tan_norm[3]], color = :green, width = 2.0, label = "")
plot!([p_cand[1], p_cand[1] + norm_norm[1]],
	[p_cand[3], p_cand[3] + norm_norm[3]], color = :cyan, width = 2.0, label = "")
