# @testset "Simulator: Flamingo" begin
    # Reference trajectory
model = deepcopy(ContactControl.get_model("flamingo", surf = "flat"))
model.μ_world = 0.1
ref_traj = deepcopy(ContactControl.get_trajectory("flamingo", "gait0", load_type = :split_traj_alt))
update_friction_coefficient!(ref_traj, model)

for t = 1:T
	r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)
	@test norm(r) < 1.0e-4
end

# initial conditions
q0 = SVector{model.dim.q}(ref_traj.q[1])
q1 = SVector{model.dim.q}(ref_traj.q[2])

h = ref_traj.h
H = ref_traj.H

# simulator
sim = ContactControl.simulator(model, q0, q1, h, 40,
    p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u]),
    ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8, κ_tol = 1.0e-8, κ_init = 1.0e-5, solver = :lu_solver),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
status = ContactControl.simulate!(sim, verbose = false)
@test status
# @test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3], Inf) < 0.15
# end

include(joinpath(pwd(), "src", "dynamics", "flamingo", "visuals.jl"))
vis = Visualizer()
render(vis)
anim = visualize_robot!(vis, model, sim.traj)
