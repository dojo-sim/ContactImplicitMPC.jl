const ContactControl = Main
include(joinpath(module_dir(), "src/dynamics/particle_2D/visuals.jl"))
vis = Visualizer()
open(vis)

################################################################################
# Parameters
################################################################################
# time
h = 0.01
T = 100

# initial conditions
q0 = SVector{2}([0.00, 0.00])
q1 = SVector{2}([0.00, 0.00])

# Simulation tolerance
tol = 1e-4

################################################################################
# Linearized Cone
################################################################################
s = get_simulation("particle_2D", "flat_2D_lc", "flat_lc")
model = s.model
env = s.env
model.μ_world = 0.5

# Policy
N_sample = 1
p = open_loop_policy([[0e-2, 0.05] for i = 1:T], N_sample = N_sample)

# simulator
sim = ContactControl.simulator(s, deepcopy(q0), deepcopy(q1), h, T,
	p = p,
	ip_opts = ContactControl.InteriorPointOptions(
		undercut = 5.0,
		γ_reg = 0.0,
		κ_tol = tol,
		r_tol = 1e-8,
		diff_sol = true,
	),
	sim_opts = ContactControl.SimulatorOptions(warmstart = false))

# simulate
@time status = ContactControl.simulate!(sim)

plot_surface!(vis, s.env, xlims=[-10,10], ylims=[-0.3,0.3])
plot_lines!(vis, s.model, sim.traj.q, name = :LClines)
anim = visualize_robot!(vis, s.model, sim.traj, name = :LC)
visualize_force!(vis, s.model, s.env, sim.traj, anim=anim, shift=-0.25, name = :LC)

plt = plot(legend=false)
plot!(plt, hcat(Vector.([q[1:end] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot!(plt, hcat(Vector.([u[1:end] for u in sim.traj.u])...)', color=:red, linewidth=4.0)
plot!(plt, hcat(Vector.([γ[1:end] for γ in sim.traj.γ])...)', color=:green, linewidth=4.0)

plot(hcat(Vector.(sim.traj.z)...)')
plot(hcat(Vector.(sim.traj.θ)...)')

plot(hcat(Vector.([q[1:1] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)
plot(hcat(Vector.([q[2:2] for q in sim.traj.q])...)', color=:cyan, linewidth=4.0)

idx = OptimizationIndices(model, env)
iq2 = index_q2(model, env)
iu1 = index_u1(model)

# Interesting smoothed gradients
dxdux = scn(sim.ip.δz[iq2[1], iu1[1]])
dydux = scn(sim.ip.δz[iq2[2], iu1[1]])
dxduy = scn(sim.ip.δz[iq2[1], iu1[2]])
dyduy = scn(sim.ip.δz[iq2[2], iu1[2]])
