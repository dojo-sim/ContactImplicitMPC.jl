# @testset "Simulator: Flamingo" begin
    # Reference trajectory
model = deepcopy(ContactControl.get_model("flamingo", surf = "flat"))
model.μ_world = 0.1
ref_traj = deepcopy(ContactControl.get_trajectory("flamingo", "gait1", load_type = :split_traj_alt))
update_friction_coefficient!(ref_traj, model)

for t = 1:ref_traj.H
	r = ContactControl.residual(model, ref_traj.z[t], ref_traj.θ[t], 0.0)
	@test norm(r) < 1.0e-4
end

# initial conditions
q0 = SVector{model.dim.q}(ref_traj.q[1])
q1 = SVector{model.dim.q}(ref_traj.q[2])

h = ref_traj.h
H = ref_traj.H
model.μ_world = 0.9

"""
    PD tracking policy
"""
mutable struct PD2 <: Policy
	model::ContactModel
    traj::ContactTraj
	q0::AbstractVector
    idx::Int
    cnt::Int
    N_sample::Int
end

function pd_policy(model, traj; N_sample = 1)
    PD2(model, traj, copy(traj.q[1]), 0, N_sample, N_sample)
end

function policy(p::PD2, x, traj, t)
    # reset
    if t == 1
        p.idx = 0
        p.cnt = p.N_sample
    end

    if p.cnt == p.N_sample
        p.idx += 1
        p.cnt = 0
		p.q0 .= copy(x)
    end

    p.cnt += 1

	u = p.traj.u[p.idx]
	# @show u
	# PD
	kp = 1.0 * ones(p.model.dim.u)
	kp[1] *= 100.0
	# kd = 10.0

	# u = Diagonal(kp) * B_func(p.model, x) * (x - p.traj.q[p.idx + 1])
	# @show u
	# u -= kd * B_func(p.model, x) * (x - traj.q[t]) / traj.h

    return u ./ p.N_sample
end

p = pd_policy(model, ref_traj)

# simulator
sim = ContactControl.simulator(model, q0, q1, h, 24,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
		r_tol = 1.0e-8,
		κ_init = 1.0e-8,
		κ_tol = 2.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

# simulate
status = ContactControl.simulate!(sim, verbose = false)
@test status
# @test norm(ref_traj.q[end][1:3] - sim.traj.q[end][1:3], Inf) < 0.15
# end

# include(joinpath(pwd(), "src", "dynamics", "flamingo", "visuals.jl"))
# vis = Visualizer()
# render(vis)
anim = visualize_robot!(vis, model, sim.traj)
