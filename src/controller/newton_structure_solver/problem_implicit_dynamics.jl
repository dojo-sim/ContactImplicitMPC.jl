sim = get_simulation("flamingo", "flat_2D_lc", "flat")
model = sim.model
env = sim.env

ref_traj = deepcopy(ContactControl.get_trajectory(model, flat_2D_lc,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

# dimensions
nq = 9
n = 2 * nq
m = 6
T = 15
nz = m * (T - 1) + n * (T - 1)
nd = n * (T - 1)

H_mpc = 10
s = newton_structure_solver(model.dim.q, model.dim.u, H_mpc)
obj_mpc = quadratic_objective(model, H_mpc)

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


# problem data
Aa = [SMatrix{nq, nq}(Array(-1.0 * lci_traj.δq0[t])) for t = 1:T-1]
Ab = [SMatrix{nq, nq}(Array(-1.0 * lci_traj.δq1[t])) for t = 1:T-1]
Ac = [SMatrix{nq, nq}(Array(Diagonal(ones(nq)))) for t = 1:T-1]
Ba = [SMatrix{nq,m}(-1.0 * lci_traj.δu1[t]) for t = 1:T-1]
Qa = [Diagonal(SVector{nq}(0.1 * rand(nq) + ones(nq))) for t = 1:T]
Qb = [Diagonal(SVector{nq}(0.1 * rand(nq) + ones(nq))) for t = 1:T]
Qv = [Diagonal(SVector{nq}(-1.0 * ones(nq) + 0.1 * rand(nq))) for t = 1:T]

# indices
u_idx = [collect((t - 1) * (m + n) .+ (1:m)) for t = 1:T-1]
x_idx = [collect((t - 1) * (m + n) + m .+ (1:n)) for t = 1:T-1]
qa_idx = [collect(x_idx[t][1:nq]) for t = 1:T-1]
qb_idx = [x_idx[t][nq .+ (1:nq)] for t = 1:T-1]
n_idx = [(t - 1) * n .+ (1:n) for t = 1:T-1]

Q = [[Qa[t] Qv[t]; Qv[t]' Qb[t]] for t = 1:T]
R = [Diagonal(SVector{m}(0.1 * rand(m) + ones(m))) for t = 1:T-1]
A = [SMatrix{n,n}([Aa[t] Ab[t]; zeros(nq, nq) -I]) for t = 1:T-1]
B = [SMatrix{n,m}([Ba[t]; zeros(nq, m)]) for t = 1:T-1]
P = [SMatrix{n,n}([zeros(nq, nq) Ac[t]; I zeros(nq, nq)]) for t = 1:T-1]

# direct problem
S = zeros(nz, nz)
C = zeros(nd, nz)

for t = 1:T-1
	S[u_idx[t], u_idx[t]] = R[t]
	S[x_idx[t], x_idx[t]] = Q[t+1]

	C[n_idx[t], u_idx[t]] = B[t]
	C[n_idx[t], x_idx[t]] = P[t]
	t == 1 && continue
	C[n_idx[t], x_idx[t-1]] = A[t]
end

Q̃ = [inv(Array(Q[t])) for t = 1:T]

Q̃a = [Q̃[t][1:nq, 1:nq] for t = 1:T]
Q̃b = [Q̃[t][nq .+ (1:nq), nq .+ (1:nq)] for t = 1:T]
Q̃v = [Q̃[t][1:nq, nq .+ (1:nq)] for t = 1:T]

R̃ = [inv(R[t]) for t = 1:T-1]
S̃ = zeros(nz, nz)

for t = 1:T-1
	S̃[u_idx[t], u_idx[t]] = R̃[t]
	S̃[x_idx[t], x_idx[t]] = Q̃[t+1]
end

norm(S̃ - inv(S))

# build kkt system
J = [S C'; C zeros(nd, nd)]
Js = sparse(J)
r = rand(nz + nd)
Δ = J \ r

Y = C * S̃ * C'
rlag = view(r, 1:nz)
rdyn = view(r, nz .+ (1:nd))
β = -rdyn + C * S̃ * rlag
Δν = Y \ β
Δz = S̃ * (rlag - C' * Δν)

# LU
solver_lu = lu_solver(Js)
Δlu = zeros(nz + nd)
@benchmark linear_solve!($solver_lu, $Δlu, $Js, $r)

solver_lu = lu_solver(J)
Δlu = zeros(nz + nd)
@benchmark linear_solve!($solver_lu, $Δlu, $J, $r)

# # QDLDL
J = sparse([S C'; C -1.0e-16 * I])
solver_ldl = ldl_solver(J)
Δldl = zeros(nz + nd)
@benchmark linear_solve!($solver_ldl, $Δldl, $J, $r)

norm(Δ - Δlu)
norm(Δ - [Δz; Δν])
