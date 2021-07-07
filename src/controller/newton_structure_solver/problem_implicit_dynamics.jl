sim = get_simulation("quadruped", "flat_2D_lc", "flat")
model = sim.model
env = sim.env

ref_traj = deepcopy(ContactControl.get_trajectory(sim.model, sim.env,
    joinpath(module_dir(), "src/dynamics/quadruped/gaits/gait2.jld2"),
    load_type = :split_traj_alt))

# dimensions
nq = model.dim.q
n = 2 * nq
m = model.dim.u
T = 15
nz = m * (T - 1) + n * (T - 1)
nd = n * (T - 1)

H_mpc = T
s = newton_structure_solver(model.dim.q, model.dim.u, H_mpc, ρ = 0.0)
obj_mpc = quadratic_objective(model, H_mpc)

ip_opts = eval(interior_point_options(:interior_point))(
			κ_init = 1.0e-4,
			κ_tol = 2.0 * 1.0e-4,
			r_tol = 1.0e-8,
			diff_sol = true,
			solver = :empty_solver)

lci_traj = LCIDynamicsTrajectory(ref_traj, sim,
	ip_type = :interior_point,
	κ = 1.0e-4,
	opts=ip_opts)

linear_contact_implicit_dynamics!(lci_traj, sim, ref_traj, κ = [1.0e-4])
update_dynamics_jacobian!(s, lci_traj)

# problem data
Aa = [SMatrix{nq, nq}(Array(lci_traj.δq0[t])) for t = 1:T-1]
Ab = [SMatrix{nq, nq}(Array(lci_traj.δq1[t])) for t = 1:T-1]
Ac = [SMatrix{nq, nq}(Array(Diagonal(ones(nq)))) for t = 1:T-1]
Ba = [SMatrix{nq,m}(lci_traj.δu1[t]) for t = 1:T-1]
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
A = [SMatrix{n,n}([zeros(nq, nq) I; Aa[t] Ab[t]]) for t = 1:T-1]
B = [SMatrix{n,m}([zeros(nq, m); Ba[t]]) for t = 1:T-1]
P = [SMatrix{n,n}([I zeros(nq, nq); zeros(nq, nq) Ac[t]]) for t = 1:T-1]

# direct problem
S = zeros(nz, nz)
C = zeros(nd, nz)

for t = 1:T-1
	S[u_idx[t], u_idx[t]] = R[t]
	S[x_idx[t], x_idx[t]] = Q[t+1]

	C[n_idx[t], u_idx[t]] = -B[t]
	C[n_idx[t], x_idx[t]] = P[t]
	t == 1 && continue
	C[n_idx[t], x_idx[t-1]] = -A[t]
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

for t = 1:T
	s.Q̃a[t] = Q̃a[t]
	s.Q̃b[t] = Q̃b[t]
	s.Q̃v[t] = Q̃v[t]

	t == T && continue
	s.Aa[t] = Aa[t]
	s.Ab[t] = Ab[t]
	s.Ac[t] = Ac[t]
	s.Ba[t] = Ba[t]

	s.R̃a[t] = R̃[t]
end

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
J = sparse([S + 1.0e-16 * I C'; C -1.0e-16 * I])
solver_ldl = ldl_solver(J)
Δldl = zeros(nz + nd)
linear_solve!(solver_ldl, Δldl, J, r)
@benchmark linear_solve!($solver_ldl, $Δldl, $J, $r)

norm(Δ - Δlu)
norm(Δ - [Δz; Δν])

norm(Δldl - Δlu)


Y = zeros(n * (T - 1), n * (T - 1))
Yii = [SMatrix{n,n}(zeros(n, n)) for t = 1:T-1]
Yij = [SMatrix{n,n}(zeros(n, n)) for t = 1:T-1]
tmp_nm = SMatrix{n,m}(zeros(n, m))
tmp_nn = SMatrix{n,n}(zeros(n, n))
tmp_nn2 = SMatrix{n,n}(zeros(n, n))

function computeY!(Yii, Yij, A, B, P, Q̃, R̃, tmp_nn, tmp_nn2, tmp_nm, T)
	for t = 1:T-1
		if t == 1
			# Yii[t] = B[t] * R̃[t] * B[t]' + P[t] * Q̃[t+1] * P[t]'

			# mul!(tmp_nn, P[t], Q̃[t+1])
			# tmp_nn .= P[t]
			# rmul!(tmp_nn, Q̃[t+1])
			Yii[t] = Q̃[t+1]
			# mul!(Yii[t], tmp_nn, transpose(P[t]))

			tmp_nm = B[t] * R̃[t]
			# mul!(tmp_nm, B[t], R̃[t])
			# tmp_nm .= B[t]
			# rmul!(tmp_nm, R̃[t])

			tmp_nn = tmp_nm * transpose(B[t])
			Yii[t] += tmp_nn
			# mul!(tmp_nn, tmp_nm, transpose(B[t]))
			# Yii[t] .+= tmp_nn
		else
			# Yii[t] = A[t] * Q̃[t] * A[t]' + B[t] * R̃[t] * B[t]' + P[t] * Q̃[t+1] * P[t]'

			tmp_nn = A[t] * Q̃[t]
			# mul!(tmp_nn, A[t], Q̃[t])
			# tmp_nn .= A[t]
			# rmul!(tmp_nn, Q̃[t])

			Yii[t] = tmp_nn * transpose(A[t])
			# mul!(Yii[t], tmp_nn, transpose(A[t]))

			tmp_nm = B[t] * R̃[t]
			# mul!(tmp_nm, B[t], R̃[t])
			# tmp_nm .= B[t]
			# rmul!(tmp_nm, R̃[t])

			tmp_nn = tmp_nm * transpose(B[t])
			Yii[t] += tmp_nn
			# mul!(tmp_nn, tmp_nm, transpose(B[t]))
			# Yii[t] .+= tmp_nn

			tmp_nn = P[t] * Q̃[t+1]
			# mul!(tmp_nn, P[t], Q̃[t+1])
			# tmp_nn .= P[t]
			# rmul!(tmp_nn, Q̃[t+1])

			tmp_nn2 = tmp_nn * transpose(P[t])
			Yii[t] += tmp_nn2
			# mul!(tmp_nn2, tmp_nn, transpose(P[t]))
			# Yii[t] .+= tmp_nn2
		end

		t == T-1 && continue
		# Yij[t] = P[t] * Q̃[t+1] * A[t+1]'

		# tmp_nn = P[t] *
		# mul!(tmp_nn, P[t], Q̃[t+1])
		# tmp_nn .= P[t]
		# rmul!(tmp_nn, Q̃[t+1])

		Yij[t] = -Q̃[t+1] * transpose(A[t+1])
		# mul!(Yij[t], tmp_nn, transpose(A[t+1]))
	end
	nothing
end

computeY!(Yii, Yij, A, B, P, Q̃, R̃, tmp_nn, tmp_nn, tmp_nm, T)
compute_Y!(s.Yiia, s.Yiib, s.Yiic, s.Yiid, s.Yija, s.Yijb, s.Yijc, s.Yijd,
	s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a,
	s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.Inq, s.H)

update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Yiiav, s.Yiibv, s.Yiicv, s.Yiidv,
	s.Yijav, s.Yijbv, s.Yijcv, s.Yijdv, s.H)

sum([norm(s.Yiis[t] - Yii[t]) for t = 1:T-1])
sum([norm(s.Yijs[t] - Yij[t]) for t = 1:T-1])

for t = 1:T-1
	Y[s.n_idx[t], s.n_idx[t]] = Yii[t]

	t == T-1 && continue

	Y[s.n_idx[t], s.n_idx[t+1]] = Yij[t]
	Y[s.n_idx[t+1], s.n_idx[t]] = Yij[t]'
end

norm((Y - C * S̃ * C'))

Y = zeros(n * (T - 1), n * (T - 1))
for t = 1:T-1
	Y[s.n_idx[t], s.n_idx[t]] = s.Yiis[t]

	t == T-1 && continue

	Y[s.n_idx[t], s.n_idx[t+1]] = s.Yijs[t]
	Y[s.n_idx[t+1], s.n_idx[t]] = s.Yijs[t]'
end

norm((Y - C * S̃ * C'))


L = zeros(n * (T - 1), n * (T - 1))
Lii = [zeros(n, n) for t = 1:T-1]
Lji = [zeros(n, n) for t = 1:T-1]

function computeL!(Lii, Lji, Yii, Yij, T)
	for t = 1:T-1
		# @show t
		if t == 1
			# Lii[t] = cholesky(Hermitian(Yii[t])).L

			Lii[t] .= Yii[t]
			LAPACK.potrf!('L', Lii[t])
		else
			# Lii[t] = cholesky(Hermitian(Yii[t] - Lji[t - 1] * Lji[t - 1]')).L

			Lii[t] .= Yii[t]
			mul!(Lii[t], transpose(Lji[t-1]), Lji[t-1], -1.0, 1.0)
			LAPACK.potrf!('L', Lii[t])

		end
		# Lji[t] = (LowerTriangular(Lii[t]) \ Yij[t])'

		Lji[t] .= Yij[t]
		LAPACK.trtrs!('L', 'N', 'N', Lii[t], Lji[t])
	end
	nothing
end

# @benchmark computeL!($Lii, $Lji, $Yii, $Yij, $T)
# @code_warntype computeL!(Lii, Lji, Yii, Yij, T)

computeL!(Lii, Lji, Yii, Yij, T)

compute_L!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)

sum([norm(s.Liis[t] - LowerTriangular(Lii[t])) for t = 1:T-1])
sum([norm(s.Ljis[t] - Lji[t]) for t = 1:T-1])


























# direct problem
S = zeros(nz, nz)
C = zeros(nd, nz)

ρ = 1.0e-16
for t = 1:T-1
	S[u_idx[t], u_idx[t]] = R[t] + ρ * I
	S[x_idx[t], x_idx[t]] = Q[t+1] + ρ * I

	C[n_idx[t], u_idx[t]] = B[t]
	C[n_idx[t], x_idx[t]] = P[t]
	t == 1 && continue
	C[n_idx[t], x_idx[t-1]] = A[t]
end

Q̃ = [inv(Array(Q[t] + ρ * I)) for t = 1:T]

Q̃a = [Q̃[t][1:nq, 1:nq] for t = 1:T]
Q̃b = [Q̃[t][nq .+ (1:nq), nq .+ (1:nq)] for t = 1:T]
Q̃v = [Q̃[t][1:nq, nq .+ (1:nq)] for t = 1:T]

R̃ = [inv(R[t] + ρ * I) for t = 1:T-1]
S̃ = zeros(nz, nz)

for t = 1:T-1
	S̃[u_idx[t], u_idx[t]] = R̃[t]
	S̃[x_idx[t], x_idx[t]] = Q̃[t+1]
end

Y = C * S̃ * C' + ρ * I
rlag = view(r, 1:nz)
rdyn = view(r, nz .+ (1:nd))
β = -rdyn + C * S̃ * rlag
Δν = Y \ β
Δz = S̃ * (rlag - C' * Δν)

norm(Δ - [Δz; Δν])
