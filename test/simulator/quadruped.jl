@testset "Simulator: Quadruped" begin
    # Reference trajectory
    model = ContactControl.get_model("quadruped", surf = "flat")
    ref_traj = ContactControl.get_trajectory("quadruped", "gait1")
    T = ref_traj.H
    h = ref_traj.h

    for t = 1:T
    	r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)
    	@test norm(r) < 1.0e-4
    end

    # initial conditions
    q0 = SVector{model.dim.q}(ref_traj.q[1])
    q1 = SVector{model.dim.q}(ref_traj.q[2])

    # simulator
    sim = ContactControl.simulator(model, q0, q1, h, T,
        p = ContactControl.open_loop_policy([SVector{model.dim.u}(ut) for ut in ref_traj.u]),
        ip_opts = ContactControl.InteriorPointOptions(
    		r_tol = 1.0e-8, κ_tol = 1.0e-8, κ_init = 1.0e-5, solver = :lu_solver),
        sim_opts = ContactControl.SimulatorOptions(warmstart = true))

    # simulate
    status = ContactControl.simulate!(sim, verbose = false)
    @test status
    @test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3], Inf) < 1.0e-3
end
