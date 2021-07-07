sim = get_simulation("flamingo", "flat_2D_lc", "flat")
model = sim.model
env = sim.env

H_mpc = 10
s = newton_structure_solver(model.dim.q, model.dim.u, H_mpc)
obj_mpc = quadratic_objective(model, H_mpc)

update_objective!(s, obj_mpc)

# random fill
for t = 1:s.H
	s.qa[t] = randn(s.nq)
	s.qb[t] = randn(s.nq)
	t == s.H && continue
	s.u[t] = rand(s.m)
	s.ν1[t] = randn(s.nq)
	s.ν2[t] = randn(s.nq)
end

for t = 1:s.H+1
	s.q_ref[t] = randn(s.nq)
	t >= s.H && continue
	s.u_ref[t] = randn(s.m)
end

z0 = zeros(s.nz + s.nd)
for t = 1:s.H-1
	z0[s.u_idx[t]] = s.u[t]
	z0[s.qa_idx[t]] = s.qa[t+1]
	z0[s.qb_idx[t]] = s.qb[t+1]
	z0[s.nz .+ (t-1) * 2 * s.nq .+ (1:s.nq)] = s.ν1[t]
	z0[s.nz .+ (t-1) * 2 * s.nq + s.nq .+ (1:s.nq)] = s.ν2[t]
end

ref_traj = deepcopy(ContactControl.get_trajectory(model, flat_2D_lc,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

traj = deepcopy(ref_traj)

ip_opts = eval(interior_point_options(:interior_point))(
			κ_init = 1.0e-4,
			κ_tol = 2.0 * 1.0e-4,
			r_tol = 1.0e-8,
			diff_sol = true,
			solver = :empty_solver)

lci_traj = LCIDynamicsTrajectory(traj, sim,
	ip_type = :interior_point,
	κ = 1.0e-4,
	opts=ip_opts)

linear_contact_implicit_dynamics!(lci_traj, sim, traj, κ = [1.0e-4])

function update_implicit_dynamics!()

end

function update_dynamics_jacobian!(s::NewtonStructureSolver, lci_traj::LCIDynamicsTrajectory)
	for t = 1:s.H-1
		s.Aa[t] = -lci_traj.δq0[t]
		s.Ab[t] = -lci_traj.δq1[t]
		# s.Ac[t] = I # pre-allocated
		s.Ba[t] = -lci_traj.δu1[t]
	end
	return nothing
end

update_dynamics_jacobian!(s, lci_traj)
info = @benchmark update_dynamics_jacobian!($s, $lci_traj)
@code_warntype update_dynamics_jacobian!(s, lci_traj)

@test info.memory == 0
@test info.allocs == 0

function dynamics_constraints_linear!(s::NewtonStructureSolver, lci_traj::LCIDynamicsTrajectory)
	# update implicit dynamics

	for t = 1:s.H-1
		# q^{t+1}_{t-1} - q^t_t = 0
		s.rdyn1[t] = s.qa[t+1] - s.qb[t]

		# q^{t+1}_t - s(q^t_{t-1}, q^t_{t}, u_t)
		s.rdyn2[t] = s.qb[t+1] - s.Aa[t] * s.qa[t] - s.Ab[t] * s.qb[t] - s.Ba[t] * s.u[t]
	end
	return nothing
end

dynamics_constraints_linear!(s, lci_traj)
info = @benchmark dynamics_constraints_linear!($s, $lci_traj)
@code_warntype dynamics_constraints_linear!(s, lci_traj)
@test info.memory == 0
@test info.allocs == 0

function dynamics_constraints!(s::NewtonStructureSolver, lci_traj::ImplicitTraj)
	# update implicit dynamics

	for t = 1:s.H-1
		# q^{t+1}_{t-1} - q^t_t = 0
		s.rdyn1[t] = s.qa[t+1] - s.qb[t]

		# q^{t+1}_t - s(q^t_{t-1}, q^t_{t}, u_t)
		s.rdyn2[t] = s.qb[t+1] - lci_traj.dq2[t] #s.Aa[t] * s.qa[t] - s.Ab[t] * s.qb[t] - s.Ba[t] * s.u[t]
	end
	return nothing
end

function lagrangian_gradient!(s::NewtonStructureSolver)
	# objective terms
	for t = 1:s.H-1
		s.rlagu[t] = s.Ra[t] * (s.u[t] - s.u_ref[t])

		s.rlagqa[t] = s.Qa[t+1] * (s.qa[t+1] - s.q_ref[t+1])
		s.rlagqb[t] = s.Qb[t+1] * (s.qb[t+1] - s.q_ref[t+2])

		s.rlagqa[t] += s.Qv[t+1] * s.qb[t+1]
		s.rlagqb[t] += transpose(s.Qv[t+1]) * s.qa[t+1]
	end

	# configuration equality terms
	for t = 1:s.H-1
		s.rlagqa[t] += s.ν1[t]
		t == 1 && continue
		s.rlagqb[t-1] -= s.ν1[t]
	end

	# dynamics terms
	for t = 1:s.H-1
		s.rlagu[t] -= transpose(s.Ba[t]) * s.ν2[t]
		s.rlagqb[t] += s.ν2[t]
		t == 1 && continue
		s.rlagqa[t-1] -= transpose(s.Aa[t]) * s.ν2[t]
		s.rlagqb[t-1] -= transpose(s.Ab[t]) * s.ν2[t]
	end
	return nothing
end

lagrangian_gradient!(s)
info = @benchmark lagrangian_gradient!($s)
@code_warntype lagrangian_gradient!(s)

@test info.memory == 0
@test info.allocs == 0

function lagrangian_dynamics(z)
	L = 0.0

	for t = 1:s.H-1
		u1 = z[s.u_idx[t]]
		qa1 = z[s.qa_idx[t]]
		qb1 = z[s.qb_idx[t]]

		ν1 = z[s.nz .+ (t-1) * 2 * s.nq .+ (1:s.nq)]
		ν2 = z[s.nz .+ (t-1) * 2 * s.nq + s.nq .+ (1:s.nq)]

		L += 0.5 * transpose(u1 - s.u_ref[t]) * s.Ra[t] * (u1 - s.u_ref[t])
		L += 0.5 * transpose(qa1 - s.q_ref[t+1]) * s.Qa[t+1] * (qa1 - s.q_ref[t+1])
		L += 0.5 * transpose(qb1 - s.q_ref[t+2]) * s.Qb[t+1] * (qb1 - s.q_ref[t+2])

		L += transpose(qa1) * s.Qv[t+1] * qb1

		if t == 1
			L += transpose(ν1) * (qa1 - s.qb[1])
			L += transpose(ν2) * (qb1 - s.Aa[1] * s.qa[1] - s.Ab[1] * s.qb[1] - s.Ba[1] * u1)
		else
			qa0 = z[s.qa_idx[t-1]]
			qb0 = z[s.qb_idx[t-1]]
			L += transpose(ν1) * (qa1 - qb0)
			L += transpose(ν2) * (qb1 - s.Aa[t] * qa0 - s.Ab[t] * qb0 - s.Ba[t] * u1)
		end
	end

	return L
end

lag_grad = vcat(vcat([[s.rlagu[t]; s.rlagqa[t]; s.rlagqb[t]] for t = 1:s.H-1])..., vcat([[s.rdyn1[t]; s.rdyn2[t]] for t = 1:s.H-1])...)
lagrangian_dynamics(z0)
lag_grad_fd = ForwardDiff.gradient(lagrangian_dynamics, z0)

@test norm((lag_grad - lag_grad_fd)[1:s.nz]) < 1.0e-12
@test norm((lag_grad - lag_grad_fd)[s.nz .+ (1:s.nd)]) < 1.0e-12
