"""
Fast MPC (modified for implicit dynamics (f(xt, ut, xt+1) = 0))
https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

mutable struct NewtonStructureSolver{N,M,NQ,NN,NM,NQNQ,NQM,T,AA,AB,AC,BA,QA,QB,RA}
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

	# dynamics Jacobians
	Aa::Vector{AA} # ∂ft/∂qt-1
	Ab::Vector{AB} # ∂ft/∂qt
	Ac::Vector{AC} # ∂ft/∂qt+1
	Ba::Vector{BA} # ∂ft/∂ut

	# objective
	Qa::Vector{QA}
	Qb::Vector{QB}
	Ra::Vector{RA}

	# objective inverse
	Q̃a::Vector{QA}
	Q̃b::Vector{QB}
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

	# r = [rlag; rdyn]
	# rlag = [rlagu; rlagqa; rlagqb]
	rlagu::Vector{SVector{M,T}}
	rlagqa::Vector{SVector{NQ,T}}
	rlagqb::Vector{SVector{NQ,T}}

	# rdyn = [rdyn1; rdyn2]
	rdyn1::Vector{SVector{NQ,T}}
	rdyn2::Vector{SVector{NQ,T}}

	# β = -rdyn + C * inv(S) * rlag = [β1..., βT-1]; β1 = [βd1; βe1], ...
	βn::Vector{SVector{N,T}}
	βd::Vector{SVector{NQ,T}}
	βe::Vector{SVector{NQ,T}}

	# y = L \ β
	y::Vector{SVector{N,T}}

	# Δν = L' \ y
	Δνn::Vector{SVector{N,T}}
	Δνd::Vector{SVector{NQ,T}}
	Δνe::Vector{SVector{NQ,T}}

	# Δz = inv(S) * (rlag - C' * Δν)
	Δzu::Vector{SVector{M,T}}
	Δzqa::Vector{SVector{NQ,T}}
	Δzqb::Vector{SVector{NQ,T}}

	# temporary arrays
	tmp_nn::SMatrix{N,N,T,NN}
	tmp_nn2::SMatrix{N,N,T,NN}
	tmp_nqnq::SMatrix{NQ,NQ,T,NQNQ}
	tmp_nqnq2::SMatrix{NQ,NQ,T,NQNQ}
	tmp_nm::SMatrix{N,M,T,NM}
	tmp_nqm::SMatrix{NQ,M,T,NQM}
	tmp_nq::SVector{NQ,T}
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

	# dynamics Jacobians
	Aa = [SMatrix{nq, nq}(zeros(nq, nq)) for t = 1:H-1]
	Ab = [SMatrix{nq, nq}(zeros(nq, nq)) for t = 1:H-1]
	Ac = [SMatrix{nq, nq}(Diagonal(ones(nq))) for t = 1:H-1]
	Ba = [SMatrix{nq,m}(zeros(nq, m)) for t = 1:H-1]

	# objective # TODO fix objective
	Qa = [SMatrix{nq,nq}(Diagonal(ones(nq))) for t = 1:H]
	Qb = [SMatrix{nq,nq}(Diagonal(ones(nq))) for t = 1:H]
	Ra = [SMatrix{m,m}(Diagonal(ones(m))) for t = 1:H-1]

	# objective inverse
	Q̃a = [inv(Qa[t]) for t = 1:H]
	Q̃b = [inv(Qb[t]) for t = 1:H]
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

	# β
	βn = [SVector{n}(zeros(n)) for t = 1:H-1]
	βd = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	βe = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# y
	y = [SVector{n}(zeros(n)) for t = 1:H-1]

	# Δν
	Δνn = [SVector{n}(zeros(n)) for t = 1:H-1]
	Δνd = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	Δνe = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# Δz
	Δzu = [SVector{m}(zeros(m)) for t = 1:H-1]
	Δzqa = [SVector{nq}(zeros(nq)) for t = 1:H-1]
	Δzqb = [SVector{nq}(zeros(nq)) for t = 1:H-1]

	# temporary arrays
	tmp_nn = SMatrix{n,n}(zeros(n, n))
	tmp_nn2 = SMatrix{n,n}(zeros(n, n))
	tmp_nqnq = SMatrix{nq,nq}(zeros(nq, nq))
	tmp_nqnq2 = SMatrix{nq,nq}(zeros(nq, nq))
	tmp_nm = SMatrix{n,m}(zeros(n, m))
	tmp_nqm = SMatrix{nq,m}(zeros(nq, m))
	tmp_nq = SVector{nq}(zeros(nq))

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
		Aa,
		Ab,
		Ac,
		Ba,
		Qa,
		Qb,
		Ra,
		Q̃a,
		Q̃b,
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
		βn,
		βd,
		βe,
		y,
		Δνn,
		Δνd,
		Δνe,
		Δzu,
		Δzqa,
		Δzqb,
		tmp_nn,
		tmp_nn2,
		tmp_nqnq,
		tmp_nqnq2,
		tmp_nm,
		tmp_nqm,
		tmp_nq)
end

function computeYs!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃a, tmp_q_nn, tmp_q_nn2, tmp_q_nm, T)

	for t = 1:T-1
		if t == 1
			# Yiia = Ac[t] * Q̃b[t+1] Ac[t]' + Ba[t] * R̃[t] * Ba[t]'

			tmp_q_nn = Ac[t] * Q̃b[t+1]
			Yiia[t] = tmp_q_nn * transpose(Ac[t])

			tmp_q_nm = Ba[t] * R̃a[t]
			tmp_q_nn = tmp_q_nm * transpose(Ba[t])
			Yiia[t] += tmp_q_nn

			Yiid[t] = Q̃a[t+1]
		else
			#Yiia[t] = Aa[t] * Q̃a[t] * Aa[t]' + Ab[t] * Q̃b[t] * Ab[t]' + Ba[t] * R̃[t] * Ba[t]' + Ac[t] * Q̃[t+1] * Ac[t]'

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

			# Yiib[t] .= -Ab[t] * Q̃b[t]
			Yiib[t] = -1.0 * Ab[t] * Q̃b[t]

			# Yiic[t] .= -Q̃b[t] * Ab[t]'
			Yiic[t] = transpose(Yiib[t])

			# Yiid[t] .= Q̃b[t] + Q̃a[t+1]
			Yiid[t] = Q̃b[t]
			Yiid[t] += Q̃a[t+1]
		end

		t == T-1 && continue
		# Yija[t] = Ac[t] * Q̃b[t+1] * Ab[t+1]
		tmp_q_nn = Ac[t] * Q̃b[t+1]
		Yija[t] = tmp_q_nn * transpose(Ab[t+1])

		# Yijb[t] = -Ac[t] * Q̃b[t+1]
		Yijb[t] = -1.0 * Ac[t] * Q̃b[t+1]

		# Yijc[t] = Q̃a[t+1] * Aa[t+1]
		Yijc[t] = Q̃a[t+1] * Aa[t+1]
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

function computeLs!(Lii, Lji, Yii, Yij, tmp_nn, tmp_nn2, T)
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

function update_β!(βn, βd, βe, rlagu, rlagqa, rlagqb, rdyn1, rdyn2, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃, T)
	for t = 1:T-1
		if t == 1
			βd[1] = -1.0 * rdyn1[1] + Ba[1] * R̃[1] * rlagu[1] + Ac[1] * Q̃b[2] * rlagqb[1]
			βe[1] = -1.0 * rdyn2[1] + Q̃a[2] * rlagqa[1]
		else
			βd[t] = -1.0 * rdyn1[t] + Ba[t] * R̃[t] * rlagu[t] + Ac[t] * Q̃b[t+1] * rlagqb[t] + Aa[t] * Q̃a[t] * rlagqa[t-1] + Ab[t] * Q̃b[t] * rlagqb[t-1]
			βe[t] = -1.0 * rdyn2[t] + Q̃a[t+1] * rlagqa[t] - Q̃b[t] * rlagqb[t-1]
		end
		βn[t] = [βd[t]; βe[t]]
	end
end

function forward_substitution_s!(y, Lii, Lji, b, T)
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

function backward_substitution_s!(x, Lii, Lji, y, T)
	for t = T-1:-1:1
		if t == T-1
			# x[t] = LowerTriangular(Lii[t])' \ y[t]

			x[t] = transpose(Lii[t]) \ y[t]
		else
			# x[t] = LowerTriangular(Lii[t])' \ (y[t] - Lji[t] * x[t+1])

			x[t] = transpose(Lii[t]) \ (y[t] - Lji[t] * x[t + 1])
		end
	end
	nothing
end

function update_Δz!(Δzu, Δzqa, Δzqb, Δνd, Δνe, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃, rlagu, rlagqa, rlagqb, T)
	for t = 1:T-1
		Δzu[t] = R̃[t] * (rlagu[t] - transpose(Ba[t]) * Δνd[t])
		if t < T-1
			Δzqa[t] = Q̃a[t+1] * (rlagqa[t] - Δνe[t] - transpose(Aa[t+1]) * Δνd[t+1])
			Δzqb[t] = Q̃b[t+1] * (rlagqb[t] - transpose(Ac[t]) * Δνd[t] - transpose(Ab[t+1]) * Δνd[t+1] + Δνe[t+1])
		else
			Δzqa[t] = Q̃a[t+1] * (rlagqa[t] - Δνe[t])
			Δzqb[t] = Q̃b[t+1] * (rlagqb[t] - transpose(Ac[t]) * Δνd[t])
		end
	end
end
