# Flamingo
function implicit_dynamics!(im_traj::ImplicitTraj, model, env, u, qa, qb, H,
	traj; κ = traj.κ) where {T, D}

	nq = model.dim.q
	m = model.dim.u
	for t = 1:H-1
		traj.θ[t][1:nq] = qa[t]
		traj.θ[t][nq .+ (1:nq)] = qb[t]
		traj.θ[t][2nq .+ (1:m)] = u[t]

		# initialized solver
		z_initialize!(im_traj.ip[t].z, model, env, copy(qb[t])) #TODO: try alt. schemes
		im_traj.ip[t].θ .= traj.θ[t]

		# solve
		status = interior_point_solve!(im_traj.ip[t])

		# !status && error("implicit dynamics failure (t = $t)")
		!status && (@warn "implicit dynamics failure (t = $t)")
	end
	return nothing
end

include(joinpath(pwd(), "src/controller/lci_dynamics.jl"))
include(joinpath(pwd(), "src/controller/newton_structure_solver/methods.jl"))

sim = get_simulation("flamingo", "flat_2D_lc", "flat")
model = sim.model
env = sim.env

ref_traj = deepcopy(ContactControl.get_trajectory(model, flat_2D_lc,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

copy_traj = deepcopy(ref_traj)

ip_opts = eval(interior_point_options(:interior_point))(
			κ_init = 1.0e-4,
			κ_tol = 2.0 * 1.0e-4,
			r_tol = 1.0e-8,
			diff_sol = true,
			solver = :empty_solver)

im_traj = ImplicitTraj(ref_traj, sim,
	ip_type = :interior_point,
	κ = 1.0e-4,
	mode = :configuration,
	opts=ip_opts)

lci_traj = LCIDynamicsTrajectory(ref_traj, sim,
	ip_type = :interior_point,
	κ = 1.0e-4,
	μ = 0.1,
	opts=ip_opts)

nq = model.dim.q
n = 2 * nq
m = model.dim.u
T = 5

nz = n * (T - 1) + m * (T - 1)
nd = n * (T - 1)

u_idx = [collect((t - 1) * (m + n) .+ (1:m)) for t = 1:T-1]
x_idx = [collect((t - 1) * (m + n) + m .+ (1:n)) for t = 1:T-1]
n_idx = [(t - 1) * n .+ (1:n) for t = 1:T-1]

implicit_dynamics!(im_traj, sim, copy_traj, κ = 1.0e-4)
implicit_dynamics!(im_traj, model, env, ref_traj.u, ref_traj.q[1:end-1], ref_traj.q[2:end], T, copy_traj, κ = 1.0e-4)
# linear_contact_implicit_dynamics!(lci_traj, ref_traj.u, ref_traj.q[1:end-1], ref_traj.q[2:end], T-1)

A = [[zeros(nq, nq) I; im_traj.δq0[t] im_traj.δq1[t]] for t = 1:T-1]
B = [[zeros(nq, m); im_traj.δu1[t]] for t = 1:T-1]
Q = [Diagonal(0.01 * ones(n)) for t = 1:T]
R = [Diagonal(0.001 * ones(m)) for t = 1:T-1]

x_ref = [[ref_traj.q[t+1]; ref_traj.q[t+2]] for t = 1:T-1]
u_ref = [ref_traj.u[t] for t = 1:T-1]

x_init = [ref_traj.q[1]; ref_traj.q[2]]
z0 = rand(nz + nd)

function constraints(z)
	c = zeros(eltype(z), nd)
	x_traj = [z[x_idx[t]] for t = 1:T-1]
	u_traj = [z[u_idx[t]] for t = 1:T-1]
	qa = [x_traj[t][1:nq] for t = 1:T-1]
	qb = [x_traj[t][nq .+ (1:nq)] for t = 1:T-1]

	# linear_contact_implicit_dynamics!(lci_traj, u_traj, qa, qb, T-1)
	implicit_dynamics!(im_traj, model, env, u_traj, qa, qb, T, copy_traj, κ = 1.0e-4)

	for t = 1:T-1
		x2 = z[x_idx[t]]
		x1 = (t == 1 ? x_init : z[x_idx[t-1]])
		u1 = z[u_idx[t]]
		c[n_idx[t]] = x2 - [x1[nq .+ (1:nq)]; im_traj.dq2[t]]#A[t] * x1 - B[t] * u1
		# c[n_idx[t]] = x2 - A[t] * x1 - B[t] * u1
		# c[n_idx[t]] = x2 - [zeros(nq, nq) I; im_traj.δq0[t] im_traj.δq1[t]] * x1 - [zeros(nq, m); im_traj.δu1[t]] * u1

	end

	return c
end

constraints(z0)

function constraints_jacobian(z)
	C = zeros(eltype(z), nd, nz)

	x_traj = [z[x_idx[t]] for t = 1:T-1]
	u_traj = [z[u_idx[t]] for t = 1:T-1]
	qa = [x_traj[t][1:nq] for t = 1:T-1]
	qb = [x_traj[t][nq .+ (1:nq)] for t = 1:T-1]

	# linear_contact_implicit_dynamics!(lci_traj, u_traj, qa, qb, T-1)
	implicit_dynamics!(im_traj, model, env, u_traj, qa, qb, T, copy_traj, κ = 1.0e-4)

	for t = 1:T-1
		C[n_idx[t], u_idx[t]] = -[zeros(nq, m); im_traj.δu1[t]]
		C[n_idx[t], x_idx[t]] = Diagonal(ones(n))
		t == 1 && continue
		C[n_idx[t], x_idx[t-1]] = -[zeros(nq, nq) I; im_traj.δq0[t] im_traj.δq1[t]]
	end

	return C
end

constraints_jacobian(z0)

function objective_gradient(z)
	j = zeros(eltype(z), nz)
	for t = 1:T-1
		u1 = z[u_idx[t]]
		x1 = z[x_idx[t]]
		j[u_idx[t]] = R[t] * (u1 - u_ref[t])
		j[x_idx[t]] = Q[t] * (x1 - x_ref[t])
	end
	return j
end

objective_gradient(z0)

function objective_hessian(z)
	J = zeros(eltype(z), nz, nz)
	for t = 1:T-1
		J[u_idx[t], u_idx[t]] = R[t]
		J[x_idx[t], x_idx[t]] = Q[t]
	end
	return J
end

objective_hessian(z0)

ρ = 1.0e-32
function hessian(z)
	J = objective_hessian(z)
	C = constraints_jacobian(z)

	return [J + ρ * I C'; C -ρ * I]
end

function gradient(z)
	y = z[nz .+ (1:nd)]
	return [objective_gradient(z) + constraints_jacobian(z)' * y; constraints(z)]
end

u = [ref_traj.u[t] for t = 1:T-1]
x = deepcopy(x_ref)
# for t = 1:T-1
# 	push!(x, A[t] * x[end], + B[t] * u[t])
# end
z = [vcat([[u[t]; x[t]] for t = 1:T-1]...); zeros(nd)]

res_norm = norm(gradient(z), 1)

Δ = hessian(z) \ gradient(z)

α = 1.0
z_cand = z - α * Δ

res_cand_norm = norm(gradient(z_cand), 1)

z_cand[nz .+ (1:nd)]

ref_traj.u[1]
z_cand[1:m]

z = z_cand
ref_traj.u[1]
z[1:m]
