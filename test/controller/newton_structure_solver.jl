@testset "Controller: Newton Structure Solver" begin
	# dimensions
	nq = 10
	n = 2 * nq
	m = 10
	T = 20
	nz = m * (T - 1) + n * (T - 1)
	nd = n * (T - 1)

	# problem data
	Aa = [SMatrix{nq, nq}(Array(1.0 * Diagonal(0.1 * rand(nq) + ones(nq)))) for t = 1:T-1]
	Ab = [SMatrix{nq, nq}(Array(1.0 * Diagonal(0.1 * rand(nq) + ones(nq)))) for t = 1:T-1]
	Ac = [SMatrix{nq, nq}(Array(Diagonal(ones(nq)))) for t = 1:T-1]
	Ba = [SMatrix{nq,m}(1.0 * rand(nq, m)) for t = 1:T-1]
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
	B = [SMatrix{n,m}([zeros(nq, m); Ba[t];]) for t = 1:T-1]
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

	s = ContactControl.newton_structure_solver(nq, m, T,
		opts = NewtonOptions(β_init = 0.0))

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

		rlag = view(r, 1:nz)

		s.rlagu[t] = view(rlag, s.u_idx[t])
		s.rlagqa[t] = view(rlag, s.qa_idx[t])
		s.rlagqb[t] = view(rlag, s.qb_idx[t])

		rdyn = view(r, nz .+ (1:nd))
		s.rdyn1[t] = view(rdyn, s.n_idx[t][1:nq])
		s.rdyn2[t] = view(rdyn, s.n_idx[t][nq .+ (1:nq)])
	end

	ContactControl.compute_Y!(s.Yiia, s.Yiib, s.Yiic, s.Yiid, s.Yija, s.Yijb, s.Yijc, s.Yijd,
		s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a,
		s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.Inq, s.H)

	ContactControl.update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
		s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Yiiav, s.Yiibv, s.Yiicv, s.Yiidv,
		s.Yijav, s.Yijbv, s.Yijcv, s.Yijdv, s.H)

	Y = zeros(n * (T - 1), n * (T - 1))
	for t = 1:T-1
		Y[s.n_idx[t], s.n_idx[t]] = s.Yiis[t]

		t == T-1 && continue

		Y[s.n_idx[t], s.n_idx[t+1]] = s.Yijs[t]
		Y[s.n_idx[t+1], s.n_idx[t]] = s.Yijs[t]'
	end

	@test norm(Y - C * S̃ * C') < 1.0e-12

	ContactControl.compute_L!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)

	L = zeros(n * (T - 1), n * (T - 1))
	for t = 1:T-1
		L[s.n_idx[t], s.n_idx[t]] = LowerTriangular(s.Liis[t])

		t == T-1 && continue

		L[s.n_idx[t+1], s.n_idx[t]] = transpose(s.Ljis[t])
	end

	@test norm(cholesky(Hermitian(Y)).L - L) < 1.0e-12

	ContactControl.compute_β!(s.βn, s.β1, s.β2, s.rlagu, s.rlagqa, s.rlagqb,
		s.rdyn1, s.rdyn2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a, s.H)

	ContactControl.compute_y!(s.y, s.Liis, s.Ljis, s.βn, s.H)

	ContactControl.compute_Δν!(s.Δνn, s.Δν1, s.Δν2, s.Liis, s.Ljis, s.y, s.idx_nq, s.idx_nq2, s.H)

	@test norm(vcat(s.Δνn...) - Δν, Inf) < 1.0e-12

	ContactControl.compute_Δz!(s.Δu, s.Δqa, s.Δqb, s.Δν1, s.Δν2, s.Aa, s.Ab, s.Ac, s.Ba,
		s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a, s.rlagu, s.rlagqa, s.rlagqb, s.H)

	@test norm((Δz - vcat([[s.Δu[t]; s.Δqa[t]; s.Δqb[t]] for t = 1:T-1]...))) < 1.0e-12

	s = ContactControl.newton_structure_solver(nq, m, T,
		opts = NewtonOptions(β_init = 0.0))

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

		rlag = view(r, 1:nz)

		s.rlagu[t] = view(rlag, s.u_idx[t])
		s.rlagqa[t] = view(rlag, s.qa_idx[t])
		s.rlagqb[t] = view(rlag, s.qb_idx[t])

		rdyn = view(r, nz .+ (1:nd))
		s.rdyn1[t] = view(rdyn, s.n_idx[t][1:nq])
		s.rdyn2[t] = view(rdyn, s.n_idx[t][nq .+ (1:nq)])
	end

	ContactControl.factorize!(s)
	ContactControl.solve!(s)

	@test norm(Δ - vcat(vcat([[s.Δu[t]; s.Δqa[t]; s.Δqb[t]] for t = 1:T-1]...), vcat([[s.Δν1[t]; s.Δν2[t]] for t = 1:T-1]...)))  < 1.0e-12
end
