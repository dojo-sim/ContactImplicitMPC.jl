"""
Fast MPC (modified for implicit dynamics (f(xt, ut, xt+1) = 0))
https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

mutable struct NewtonStructureSolver{N,M,NQ,NN,NM,NQNQ,NQM,T,AA,AB,AC,BA,QA,QB,QV,QAI,QBI,QVI,RA}
	"""
	[S C'; C 0] [Δz; Δν] = r
	"""
	# dimensions
	n::Int  # state
	m::Int  # control
	nq::Int # configuration
	H::Int  # horizon
	nz::Int # decision variables
	nd::Int # dual variables

	# indices
	# [u1, x2, u2, ..., xT]
	x_idx::Vector{Vector{Int}}
	u_idx::Vector{Vector{Int}}
	qa_idx::Vector{Vector{Int}}
	qb_idx::Vector{Vector{Int}}
	n_idx::Vector{Vector{Int}}

	# trajectories x = [qa; qb]
	u::Vector{SVector{M,T}}
	qa::Vector{SVector{NQ,T}}
	qb::Vector{SVector{NQ,T}}

	ν1::Vector{SVector{NQ,T}}
	ν2::Vector{SVector{NQ,T}}

	# candidate trajectories
	u_cand::Vector{SVector{M,T}}
	qa_cand::Vector{SVector{NQ,T}}
	qb_cand::Vector{SVector{NQ,T}}

	ν1_cand::Vector{SVector{NQ,T}}
	ν2_cand::Vector{SVector{NQ,T}}

	# dynamics Jacobians
	Aa::Vector{AA} # ∂ft/∂qt-1
	Ab::Vector{AB} # ∂ft/∂qt
	Ac::Vector{AC} # ∂ft/∂qt+1
	Ba::Vector{BA} # ∂ft/∂ut

	# objective
	q_ref::Vector{SVector{NQ,T}}
	u_ref::Vector{SVector{M,T}}

	Qa::Vector{QA}
	Qb::Vector{QB}
	Qv::Vector{QV}
	Ra::Vector{RA}

	# objective inverse
	Q̃a::Vector{QAI}
	Q̃b::Vector{QBI}
	Q̃v::Vector{QVI}
	R̃a::Vector{RA}

	# Y = C * inv(S) * C' = [Y11 Y12...;Y21 Y22...;...;...YT-1T-1]
	Yii::Vector{Matrix{T}}
	Yiis::Vector{SMatrix{N,N,T,NN}}
	Yij::Vector{Matrix{T}}
	Yijs::Vector{SMatrix{N,N,T,NN}}

	# Yii = [Yiia Yiib; Yiic Yiid]
	Yiia::Vector{SMatrix{NQ,NQ,T,NQNQ}}
	Yiib::Vector{SMatrix{NQ,NQ,T,NQNQ}}
	Yiic::Vector{SMatrix{NQ,NQ,T,NQNQ}}
	Yiid::Vector{SMatrix{NQ,NQ,T,NQNQ}}

	Yiiav::Vector{SubArray{T,2,Matrix{T},Tuple{Vector{Int},Vector{Int}},false}}
	Yiibv::Vector{SubArray{T,2,Matrix{T},Tuple{Vector{Int},Vector{Int}},false}}
	Yiicv::Vector{SubArray{T,2,Matrix{T},Tuple{Vector{Int},Vector{Int}},false}}
	Yiidv::Vector{SubArray{T,2,Matrix{T},Tuple{Vector{Int},Vector{Int}},false}}

	# Yij = [Yija Yijb; Yijc Yijd]
	Yija::Vector{SMatrix{NQ,NQ,T,NQNQ}}
	Yijb::Vector{SMatrix{NQ,NQ,T,NQNQ}}
	Yijc::Vector{SMatrix{NQ,NQ,T,NQNQ}}
	Yijd::Vector{SMatrix{NQ,NQ,T,NQNQ}}

	Yijav::Vector{SubArray{T,2,Matrix{T},Tuple{Vector{Int},Vector{Int}},false}}
	Yijbv::Vector{SubArray{T,2,Matrix{T},Tuple{Vector{Int},Vector{Int}},false}}
	Yijcv::Vector{SubArray{T,2,Matrix{T},Tuple{Vector{Int},Vector{Int}},false}}
	Yijdv::Vector{SubArray{T,2,Matrix{T},Tuple{Vector{Int},Vector{Int}},false}}

	# L = cholesky(Y)
	Liis::Vector{LowerTriangular{T,SMatrix{N,N,T,NN}}}
	Ljis::Vector{SMatrix{N,N,T,NN}}

	# residual: r = [rlag; rdyn]
	# rlag = [rlagu; rlagqa; rlagqb]
	rlagu::Vector{SVector{M,T}}
	rlagqa::Vector{SVector{NQ,T}}
	rlagqb::Vector{SVector{NQ,T}}

	# rdyn = [rdyn1; rdyn2]
	rdyn1::Vector{SVector{NQ,T}}
	rdyn2::Vector{SVector{NQ,T}}

	# candidate residual
	rlagu_cand::Vector{SVector{M,T}}
	rlagqa_cand::Vector{SVector{NQ,T}}
	rlagqb_cand::Vector{SVector{NQ,T}}

	# rdyn = [rdyn1; rdyn2]
	rdyn1_cand::Vector{SVector{NQ,T}}
	rdyn2_cand::Vector{SVector{NQ,T}}

	# β = -rdyn + C * inv(S) * rlag = [β1..., βT-1]; β1 = [βd1; βe1], ...
	βn::Vector{SVector{N,T}}
	βd::Vector{SVector{NQ,T}}
	βe::Vector{SVector{NQ,T}}

	# y = L \ β
	y::Vector{SVector{N,T}}

	# Δν = L' \ y
	Δνn::Vector{SVector{N,T}}
	Δν1::Vector{SVector{NQ,T}}
	Δν2::Vector{SVector{NQ,T}}

	# Δz = inv(S) * (rlag - C' * Δν)
	Δu::Vector{SVector{M,T}}
	Δqa::Vector{SVector{NQ,T}}
	Δqb::Vector{SVector{NQ,T}}

	# temporary arrays
	tmp_nn::SMatrix{N,N,T,NN}
	tmp_nn2::SMatrix{N,N,T,NN}
	tmp_nqnq::SMatrix{NQ,NQ,T,NQNQ}
	tmp_nqnq2::SMatrix{NQ,NQ,T,NQNQ}
	tmp_nm::SMatrix{N,M,T,NM}
	tmp_nqm::SMatrix{NQ,M,T,NQM}
	tmp_nq::SVector{NQ,T}

	# idx
	idx_nq::SVector{NQ, Int}
	idx_nq2::SVector{NQ, Int}
end

function newton_structure_solver(nq::Int, m::Int, H::Int)
	# dimensions
	n = 2 * nq
	nz = m * (H - 1) + n * (H - 1)
	nd = n * (H - 1)

	# indices
	u_idx = [collect((t - 1) * (m + n) .+ (1:m)) for t = 1:H-1]
	x_idx = [collect((t - 1) * (m + n) + m .+ (1:n)) for t = 1:H-1]
	qa_idx = [collect(x_idx[t][1:nq]) for t = 1:H-1]
	qb_idx = [collect(x_idx[t][nq .+ (1:nq)]) for t = 1:H-1]
	n_idx = [collect((t - 1) * n .+ (1:n)) for t = 1:H-1]

	# trajectories
	u = [SVector{m}(zeros(m)) for t = 1:H-1]
	qa = [SVector{nq}(zeros(nq)) for t = 1:H]
	qb = [SVector{nq}(zeros(nq)) for t = 1:H]

	ν1 = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	ν2 = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# candidate trajectories
	u_cand = [SVector{m}(zeros(m)) for t = 1:H-1]
	qa_cand = [SVector{nq}(zeros(nq)) for t = 1:H]
	qb_cand = [SVector{nq}(zeros(nq)) for t = 1:H]

	ν1_cand = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	ν2_cand = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# dynamics Jacobians
	Aa = [SMatrix{nq, nq}(zeros(nq, nq)) for t = 1:H-1]
	Ab = [SMatrix{nq, nq}(zeros(nq, nq)) for t = 1:H-1]
	Ac = [SMatrix{nq, nq}(Diagonal(ones(nq))) for t = 1:H-1]
	Ba = [SMatrix{nq,m}(zeros(nq, m)) for t = 1:H-1]

	# objective
	q_ref = [SVector{nq}(zeros(nq)) for t = 1:H+1]
	u_ref = [SVector{m}(zeros(m)) for t = 1:H-1]

	Qa = [SMatrix{nq,nq}(Diagonal(ones(nq))) for t = 1:H]
	Qb = [SMatrix{nq,nq}(Diagonal(ones(nq))) for t = 1:H]
	Qv = [SMatrix{nq,nq}(Diagonal(0.1 * ones(nq))) for t = 1:H]
	Ra = [SMatrix{m,m}(Diagonal(ones(m))) for t = 1:H-1]

	# objective inverse
	Q̃ = [inv(Array([Qa[t] Qv[t]; Qv[t]' Qb[t]])) for t = 1:H]

	Q̃a = [SMatrix{nq,nq}(Q̃[t][1:nq, 1:nq]) for t = 1:H]
	Q̃b = [SMatrix{nq,nq}(Q̃[t][nq .+ (1:nq), nq .+ (1:nq)]) for t = 1:H]
	Q̃v = [SMatrix{nq,nq}(Q̃[t][1:nq, nq .+ (1:nq)]) for t = 1:H]
	R̃a = [inv(Ra[t]) for t = 1:H-1]

	# Y
	Yii = [zeros(n, n) for t = 1:H-1]
	Yij = [zeros(n, n) for t = 1:H-1]
	Yiis = [SMatrix{n,n}(zeros(n, n)) for t = 1:H-1]
	Yijs = [SMatrix{n,n}(zeros(n, n)) for t = 1:H-1]

	Yiiav = [view(Yii[t], collect(1:nq), collect(1:nq)) for t = 1:H-1]
	Yiibv = [view(Yii[t], collect(1:nq), collect(nq .+ (1:nq))) for t = 1:H-1]
	Yiicv = [view(Yii[t], collect(nq .+ (1:nq)), collect(1:nq)) for t = 1:H-1]
	Yiidv = [view(Yii[t], collect(nq .+ (1:nq)), collect(nq .+ (1:nq))) for t = 1:H-1]

	Yijav = [view(Yij[t], collect(1:nq), collect(1:nq)) for t = 1:H-1]
	Yijbv = [view(Yij[t], collect(1:nq), collect(nq .+ (1:nq))) for t = 1:H-1]
	Yijcv = [view(Yij[t], collect(nq .+ (1:nq)), collect(1:nq)) for t = 1:H-1]
	Yijdv = [view(Yij[t], collect(nq .+ (1:nq)), collect(nq .+ (1:nq))) for t = 1:H-1]

	Yiia = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:H-1]
	Yiib = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:H-1]
	Yiic = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:H-1]
	Yiid = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:H-1]

	Yija = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:H-1]
	Yijb = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:H-1]
	Yijc = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:H-1]
	Yijd = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:H-1]

	# L
	Liis = [LowerTriangular(SMatrix{n,n}(zeros(n, n))) for t = 1:H-1]
	Ljis = [SMatrix{n,n}(zeros(n, n)) for t = 1:H-1]

	# r
	rlagu = [SVector{m}(zeros(m)) for t = 1:H-1]
	rlagqa = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	rlagqb = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	rdyn1 = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	rdyn2 = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# r
	rlagu_cand = [SVector{m}(zeros(m)) for t = 1:H-1]
	rlagqa_cand = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	rlagqb_cand = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	rdyn1_cand = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	rdyn2_cand = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# β
	βn = [SVector{n}(zeros(n)) for t = 1:H-1]
	βd = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	βe = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# y
	y = [SVector{n}(zeros(n)) for t = 1:H-1]

	# Δν
	Δνn = [SVector{n}(zeros(n)) for t = 1:H-1]
	Δν1 = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	Δν2 = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# Δz
	Δu = [SVector{m}(zeros(m)) for t = 1:H-1]
	Δqa = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	Δqb = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# temporary arrays
	tmp_nn = SMatrix{n,n}(zeros(n, n))
	tmp_nn2 = SMatrix{n,n}(zeros(n, n))
	tmp_nqnq = SMatrix{nq,nq}(zeros(nq, nq))
	tmp_nqnq2 = SMatrix{nq,nq}(zeros(nq, nq))
	tmp_nm = SMatrix{n,m}(zeros(n, m))
	tmp_nqm = SMatrix{nq,m}(zeros(nq, m))
	tmp_nq = SVector{nq}(zeros(nq))

	# idx
	idx_nq = SVector{nq}(collect(1:nq))
	idx_nq2 = SVector{nq}(collect(nq .+ (1:nq)))

	NewtonStructureSolver(
		n,
		m,
		nq,
		H,
		nz,
		nd,
		x_idx,
		u_idx,
		qa_idx,
		qb_idx,
		n_idx,
		u,
		qa,
		qb,
		ν1,
		ν2,
		u_cand,
		qa_cand,
		qb_cand,
		ν1_cand,
		ν2_cand,
		Aa,
		Ab,
		Ac,
		Ba,
		q_ref,
		u_ref,
		Qa,
		Qb,
		Qv,
		Ra,
		Q̃a,
		Q̃b,
		Q̃v,
		R̃a,
		Yii,
		Yiis,
		Yij,
		Yijs,
		Yiia,
		Yiib,
		Yiic,
		Yiid,
		Yiiav,
		Yiibv,
		Yiicv,
		Yiidv,
		Yija,
		Yijb,
		Yijc,
		Yijd,
		Yijav,
		Yijbv,
		Yijcv,
		Yijdv,
		Liis,
		Ljis,
		rlagu,
		rlagqa,
		rlagqb,
		rdyn1,
		rdyn2,
		rlagu_cand,
		rlagqa_cand,
		rlagqb_cand,
		rdyn1_cand,
		rdyn2_cand,
		βn,
		βd,
		βe,
		y,
		Δνn,
		Δν1,
		Δν2,
		Δu,
		Δqa,
		Δqb,
		tmp_nn,
		tmp_nn2,
		tmp_nqnq,
		tmp_nqnq2,
		tmp_nm,
		tmp_nqm,
		tmp_nq,
		idx_nq,
		idx_nq2)
end

function compute_Y!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Aa, Ab, Ac, Ba, Q̃a, Q̃b, Q̃v, R̃a, tmp_q_nn, tmp_q_nn2, tmp_q_nm, T)

	for t = 1:T-1
		if t == 1
			# Yiia = Ac[t] * Q̃b[t+1] Ac[t]' + Ba[t] * R̃[t] * Ba[t]'

			tmp_q_nn = Ac[t] * Q̃b[t+1]
			Yiia[t] = tmp_q_nn * transpose(Ac[t])

			tmp_q_nm = Ba[t] * R̃a[t]
			tmp_q_nn = tmp_q_nm * transpose(Ba[t])
			Yiia[t] += tmp_q_nn

			# Yiib = Ac[t] * Q̃v[t+1]'
			Yiib[t] = Ac[t] * transpose(Q̃v[t+1])

			# Yiic = Q̃v[t+1] * Ac[t]'
			Yiic[t] = transpose(Yiib[t])

			# Yiid
			Yiid[t] = Q̃a[t+1]

		else
			#Yiia[t] = Aa[t] * Q̃a[t] * Aa[t]' + Ab[t] * Q̃b[t] * Ab[t]' + Ba[t] * R̃[t] * Ba[t]' + Ac[t] * Q̃b[t+1] * Ac[t]' + Aa[t] * Q̃v[t+1] * Ab[t]' + Ab[t] * Q̃v[t+1] * Aa[t+1]'

			tmp_q_nn = Aa[t] * Q̃a[t]
			Yiia[t] = tmp_q_nn * transpose(Aa[t])
			tmp_q_nn = Ab[t] * Q̃b[t]
			tmp_q_nn2 = tmp_q_nn * transpose(Ab[t])
			Yiia[t] += tmp_q_nn2

			tmp_q_nm = Ba[t] * R̃a[t]
			tmp_q_nn = tmp_q_nm * transpose(Ba[t])
			Yiia[t] += tmp_q_nn

			tmp_q_nn = Ac[t] * Q̃b[t+1]
			tmp_q_nn2 = tmp_q_nn * transpose(Ac[t])
			Yiia[t] += tmp_q_nn2

			tmp_q_nn = Aa[t] * Q̃v[t]
			tmp_q_nn2 = tmp_q_nn * transpose(Ab[t])
			Yiia[t] += tmp_q_nn2
			tmp_q_nn = Ab[t] * transpose(Q̃v[t])
			tmp_q_nn2 = tmp_q_nn * transpose(Aa[t])
			Yiia[t] += tmp_q_nn2

			# Yiib[t] .= -Ab[t] * Q̃b[t] - Aa[t] * Q̃v[t] + Ac[t] * Q̃v[t+1]
			Yiib[t] = -1.0 * Ab[t] * Q̃b[t]
			Yiib[t] += -1.0 * Aa[t] * Q̃v[t]
			Yiib[t] += Ac[t] * Q̃v[t+1]

			# Yiic[t] .= -Q̃b[t] * Ab[t]' - Q̃v[t]' * Ac[t]'
			Yiic[t] = transpose(Yiib[t])

			# Yiid[t] .= Q̃b[t] + Q̃a[t+1]
			Yiid[t] = Q̃b[t]
			Yiid[t] += Q̃a[t+1]
		end

		t == T-1 && continue
		# Yija[t] = Ac[t] * Q̃b[t+1] * Ab[t+1] + Ac[t] * Q̃v[t+1] * Aa[t+1]'
		tmp_q_nn = Ac[t] * Q̃b[t+1]
		Yija[t] = tmp_q_nn * transpose(Ab[t+1])
		tmp_q_nn = Ac[t] * Q̃v[t+1]
		tmp_q_nn2 = tmp_q_nn * transpose(Aa[t+1])
		Yija[t] += tmp_q_nn2

		# Yijb[t] = -Ac[t] * Q̃b[t+1]
		Yijb[t] = -1.0 * Ac[t] * Q̃b[t+1]

		# Yijc[t] = Q̃a[t+1] * Aa[t+1] + Q̃v[t+1] * Ab[t+1]'
		Yijc[t] = Q̃a[t+1] * Aa[t+1]
		Yijc[t] += Q̃v[t+1] * transpose(Ab[t+1])

		# Yijd[t] = -Q̃v[t+1]
		Yijd[t] = -1.0 * Q̃v[t+1]

	end
	nothing
end

function update_Y!(Yiis, Yijs, Yii, Yij, Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Yiiav, Yiibv, Yiicv, Yiidv, Yijav, Yijbv, Yijcv, Yijdv, T)
	for t = 1:T-1
		Yiiav[t] .= Yiia[t]
		Yiibv[t] .= Yiib[t]
		Yiicv[t] .= Yiic[t]
		Yiidv[t] .= Yiid[t]

		Yijav[t] .= Yija[t]
		Yijbv[t] .= Yijb[t]
		Yijcv[t] .= Yijc[t]
		Yijdv[t] .= Yijd[t]

		Yiis[t] = Yii[t]
		Yijs[t] = Yij[t]
	end
	nothing
end

function compute_L!(Lii, Lji, Yii, Yij, tmp_nn, tmp_nn2, T)
	for t = 1:T-1
		if t == 1
			Lii[t] = cholesky(Hermitian(Yii[t])).L

			# Lii[t] = cholesky(Yii[t]).L
		else
			# Lii[t] = cholesky(Hermitian(Yii[t] - Lji[t - 1] * Lji[t - 1]')).L

			tmp_nn = Yii[t]
			tmp_nn2 = transpose(Lji[t-1]) * Lji[t-1]
			tmp_nn -= tmp_nn2
			Lii[t] = cholesky(Hermitian(tmp_nn)).L
		end
		# Lji[t] = (LowerTriangular(Lii[t]) \ Yij[t])'

		Lji[t] = Lii[t] \ Yij[t]
	end
	nothing
end

function compute_β!(βn, βd, βe, rlagu, rlagqa, rlagqb, rdyn1, rdyn2, Aa, Ab, Ac, Ba, Q̃a, Q̃b, Q̃v, R̃, T)
	for t = 1:T-1
		if t == 1
			βd[1] = -1.0 * rdyn1[1] + Ba[1] * R̃[1] * rlagu[1] + Ac[1] * Q̃b[2] * rlagqb[1]
			βd[1] += Ac[1] * Q̃v[2] * rlagqa[1]

			βe[1] = -1.0 * rdyn2[1] + Q̃a[2] * rlagqa[1]
			βe[1] += Q̃v[2] * rlagqb[1]
		else
			βd[t] = -1.0 * rdyn1[t] + Ba[t] * R̃[t] * rlagu[t] + Ac[t] * Q̃b[t+1] * rlagqb[t] + Aa[t] * Q̃a[t] * rlagqa[t-1] + Ab[t] * Q̃b[t] * rlagqb[t-1]
			βd[t] += Aa[t] * Q̃v[t] * rlagqb[t-1] + Ab[t] * transpose(Q̃v[t]) * rlagqa[t-1] + Ac[t] * Q̃v[t+1] * rlagqa[t]

			βe[t] = -1.0 * rdyn2[t] + Q̃a[t+1] * rlagqa[t] - Q̃b[t] * rlagqb[t-1]
			βe[t] += -1.0 * Q̃v[t] * rlagqa[t-1] + Q̃v[t+1] * rlagqb[t]
		end
		βn[t] = [βd[t]; βe[t]]
	end
end

function compute_y!(y, Lii, Lji, b, T)
	for t = 1:T-1
		if t == 1
			# y[1] .= Lii[1] \ β[n_idx[1]]

			y[1] = Lii[1] \ b[1]
		else
			# y[t] .= Lii[t] \ (β[n_idx[t]] - transpose(Lji[t - 1]) * y[t - 1])

			y[t] = Lii[t] \ (b[t] - transpose(Lji[t - 1]) * y[t - 1])
		end
	end
	nothing
end

function compute_Δν!(x, x1, x2, Lii, Lji, y, idx1, idx2, T)
	for t = T-1:-1:1
		if t == T-1
			# x[t] = LowerTriangular(Lii[t])' \ y[t]

			x[t] = transpose(Lii[t]) \ y[t]
		else
			# x[t] = LowerTriangular(Lii[t])' \ (y[t] - Lji[t] * x[t+1])

			x[t] = transpose(Lii[t]) \ (y[t] - Lji[t] * x[t + 1])
		end

		x1[t] = x[t][idx1]
		x2[t] = x[t][idx2]
	end
	nothing
end

function compute_Δz!(Δu, Δqa, Δqb, Δν1, Δν2, Aa, Ab, Ac, Ba, Q̃a, Q̃b, Q̃v, R̃, rlagu, rlagqa, rlagqb, T)
	for t = 1:T-1
		Δu[t] = R̃[t] * (rlagu[t] - transpose(Ba[t]) * Δν1[t])

		if t < T-1
			Δqa[t] = Q̃a[t+1] * (rlagqa[t] - Δν2[t] - transpose(Aa[t+1]) * Δν1[t+1])
			Δqa[t] += Q̃v[t+1] * (rlagqb[t] - transpose(Ac[t]) * Δν1[t] - transpose(Ab[t+1]) * Δν1[t+1] + Δν2[t+1])

			Δqb[t] = Q̃b[t+1] * (rlagqb[t] - transpose(Ac[t]) * Δν1[t] - transpose(Ab[t+1]) * Δν1[t+1] + Δν2[t+1])
			Δqb[t] += transpose(Q̃v[t+1]) * (rlagqa[t] - Δν2[t] - transpose(Aa[t+1]) * Δν1[t+1])
		else
			Δqa[t] = Q̃a[t+1] * (rlagqa[t] - Δν2[t])
			Δqa[t] += Q̃v[t+1] * (rlagqb[t] - transpose(Ac[t]) * Δν1[t])

			Δqb[t] = Q̃b[t+1] * (rlagqb[t] - transpose(Ac[t]) * Δν1[t])
			Δqb[t] += transpose(Q̃v[t+1]) * (rlagqa[t] - Δν2[t])
		end
	end
end

function ContactControl.solve!(s::NewtonStructureSolver)
	compute_Y!(s.Yiia, s.Yiib, s.Yiic, s.Yiid, s.Yija, s.Yijb, s.Yijc, s.Yijd,
		s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a,
		s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.H)

	update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
		s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Yiiav, s.Yiibv, s.Yiicv, s.Yiidv,
		s.Yijav, s.Yijbv, s.Yijcv, s.Yijdv, s.H)

	compute_L!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)

	compute_β!(s.βn, s.βd, s.βe, s.rlagu, s.rlagqa, s.rlagqb,
		s.rdyn1, s.rdyn2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a, s.H)

	compute_y!(s.y, s.Liis, s.Ljis, s.βn, s.H)

	compute_Δν!(s.Δνn, s.Δν1, s.Δν2, s.Liis, s.Ljis, s.y, s.idx_nq, s.idx_nq2, s.H)

	compute_Δz!(s.Δu, s.Δqa, s.Δqb, s.Δν1,  s.Δν2, s.Aa, s.Ab, s.Ac, s.Ba,
		s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a, s.rlagu, s.rlagqa, s.rlagqb, s.H)
end

# objective
mutable struct QuadraticObjective{Q,V,U} <: Objective
    q::Vector{Q}
	v::Vector{V}
    u::Vector{U}
end

function quadratic_objective(model::ContactModel, H::Int;
    q = [Diagonal(SVector{model.dim.q}(ones(model.dim.q))) for t = 1:H+1],
	v = [Diagonal(SVector{model.dim.q}(0.1 * ones(model.dim.q))) for t = 1:H],
    u = [Diagonal(SVector{model.dim.q}(ones(model.dim.q))) for t = 1:H-1])
    return QuadraticObjective(q, v, u)
end

function update_objective!(s::NewtonStructureSolver, obj::QuadraticObjective)
	n = s.n
	m = s.m
	nq = s.nq
	H = s.H

	Q = [[(t > 1 ? 0.5 : 1.0) * obj.q[t] + obj.v[t]  -1.0 * obj.v[t];
	      -1.0 * obj.v[t]' (t < H ? 0.5 : 1.0) * obj.q[t] + obj.v[t]] for t = 1:H]
	Q̃ = [inv(Array(Q[t])) for t = 1:H]

	R̃ = [inv(obj.u[t]) for t = 1:H-1]

	for t = 1:H
		s.Qa[t] = Q[t][1:nq, 1:nq]
		s.Qb[t] = Q[t][nq .+ (1:nq), nq .+ (1:nq)]
		s.Qv[t] = Q[t][1:nq, nq .+ (1:nq)]

		s.Q̃a[t] = Q̃[t][1:nq, 1:nq]
		s.Q̃b[t] = Q̃[t][nq .+ (1:nq), nq .+ (1:nq)]
		s.Q̃v[t] = Q̃[t][1:nq, nq .+ (1:nq)]
	end
end

function update_dynamics_jacobian!(s::NewtonStructureSolver, lci_traj::LCIDynamicsTrajectory)
	for t = 1:s.H-1
		s.Aa[t] = lci_traj.δq0[t]
		s.Ab[t] = lci_traj.δq1[t]
		# s.Ac[t] = I # pre-allocated
		s.Ba[t] = lci_traj.δu1[t]
	end
	return nothing
end

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

function compute_residual!(s::NewtonStructureSolver, lci_traj::LCIDynamicsTrajectory)

	# update lci_traj

	# update static Jacobians
	update_dynamics_jacobian!(s, lci_traj)

	# update constraints
	dynamics_constraints!(s, lci_traj)

	# update gradient of Lagrangian
	lagrangian_gradient!(s)
end
