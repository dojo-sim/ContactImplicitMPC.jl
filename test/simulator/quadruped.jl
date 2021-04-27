@testset "Simulator: Quadruped" begin
    # Reference trajectory
    model = deepcopy(ContactControl.get_model("quadruped", surf = "flat"))
    model.μ_world = 0.5

    ref_traj = deepcopy(ContactControl.get_trajectory("quadruped", "gait1", load_type = :split_traj_alt))
    ContactControl.update_friction_coefficient!(ref_traj, model)

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
    @test status = ContactControl.simulate!(sim, verbose = false)
    @show sim.traj.q[end][1:3]
    @test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3], Inf) < 0.025
end

# Reference trajectory
model = deepcopy(ContactControl.get_model("quadruped", surf = "flat"))
model.μ_world = 0.5

ref_traj = deepcopy(ContactControl.get_trajectory("quadruped", "gait2", load_type = :split_traj_alt))
ContactControl.update_friction_coefficient!(ref_traj, model)

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
@test status = ContactControl.simulate!(sim, verbose = false)
@show sim.traj.q[end][1:3]
@test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3]) < 0.01

# time
h = ref_traj.h
H = ref_traj.H

# initial conditions
q0 = SVector{model.dim.q}(ref_traj.q[1])
q1 = SVector{model.dim.q}(ref_traj.q[2])

function ContactControl.z_initialize!(z, model::ContactControl.Quadruped, q1)
    nq = model.dim.q
    z .= 1.0e-1
    z[1:nq] = q1
end

ref_traj1 = deepcopy(ref_traj)
for t = 1:H+2
    ref_traj1.q[t] .+= 1.0e-0 * rand(model.dim.q)
end

im_traj0 = ContactControl.ImplicitTraj(ref_traj, model)
im_traj1 = ContactControl.ImplicitTraj(ref_traj1, model)

ContactControl.implicit_dynamics!(im_traj0, model, ref_traj)

im_traj0.ip[1].r
ref_traj.b[1]
im_traj0.db1[1]
norm(im_traj0.d[1])
@test maximum([norm(d) for d in im_traj0.d]) < 5.0e-3
ContactControl.implicit_dynamics!(im_traj1, model, ref_traj1)
@test maximum(norm.(im_traj1.d, 2)) > 5.0e-1
