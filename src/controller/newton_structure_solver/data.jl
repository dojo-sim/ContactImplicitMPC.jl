"""
Fast MPC (modified for implicit dynamics (f(xt, ut, xt+1) = 0))
https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

using BenchmarkTools
using InteractiveUtils

struct NewtonStructureSolver{N,M,NQ,NN,NM,NQNQ,NQM,T,AA,AB,AC,BA,QA,QB,RA}
	"""
	[P C'; C 0] [Δz; Δν] = r
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

	# Y = C * inv(H) * C' = [Y11 Y12...;Y21 Y22...;...;...YT-1T-1]
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

	# β = -rdyn + C * inv(H) * rlag = [β1..., βT-1]; β1 = [βd1; βe1], ...
	βn::Vector{SVector{N,T}}
	βd::Vector{SVector{NQ,T}}
	βe::Vector{SVector{NQ,T}}

	# y = L \ β
	y::Vector{SVector{N,T}}

	# Δν = L' \ y
	Δνn::Vector{SVector{N,T}}
	Δνd::Vector{SVector{NQ,T}}
	Δνe::Vector{SVector{NQ,T}}

	# Δz = inv(H) * (rlag - C' * Δν)
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

function newton_structure_solver(nq::Int, m::Int, H::Int, Q, R)
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
	Qa = [SMatrix{nq,nq}(Q[t]) for t = 1:H]
	Qb = [SMatrix{nq,nq}(Q[t]) for t = 1:H]
	Ra = [SMatrix{m,m}(R[t]) for t = 1:H-1]

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

# dimensions
nq = 10
n = 2 * nq
m = 10
T = 20
nz = m * (T - 1) + n * (T - 1)
nd = n * (T - 1)

# indices
u_idx = [collect((t - 1) * (m + n) .+ (1:m)) for t = 1:T-1]
x_idx = [collect((t - 1) * (m + n) + m .+ (1:n)) for t = 1:T-1]
qa_idx = [collect(x_idx[t][1:nq]) for t = 1:T-1]
qb_idx = [x_idx[t][nq .+ (1:nq)] for t = 1:T-1]
n_idx = [(t - 1) * n .+ (1:n) for t = 1:T-1]

# problem data
Aa = [SMatrix{nq, nq}(Array(-1.0 * Diagonal(0.1 * rand(nq) + ones(nq)))) for t = 1:T-1]
Ab = [SMatrix{nq, nq}(Array(-1.0 * Diagonal(0.1 * rand(nq) + ones(nq)))) for t = 1:T-1]
Ac = [SMatrix{nq, nq}(Array(Diagonal(ones(nq)))) for t = 1:T-1]
Ba = [SMatrix{nq,m}(-1.0 * rand(nq, m)) for t = 1:T-1]
Qa = [Diagonal(SVector{nq}(0.1 * rand(nq) + ones(nq))) for t = 1:T]
Qb = [Diagonal(SVector{nq}(0.1 * rand(nq) + ones(nq))) for t = 1:T]

Q = [Diagonal(SVector{n}(Array(diag(cat(Qa[t], Qb[t], dims=(1,2)))))) for t = 1:T]
R = [Diagonal(SVector{m}(0.1 * rand(m) + ones(m))) for t = 1:T-1]
A = [SMatrix{n,n}([Aa[t] Ab[t]; zeros(nq, nq) -I]) for t = 1:T-1]
B = [SMatrix{n,m}([Ba[t]; zeros(nq, m)]) for t = 1:T-1]
P = [SMatrix{n,n}([zeros(nq, nq) Ac[t]; I zeros(nq, nq)]) for t = 1:T-1]

# direct problem
H = zeros(nz, nz)
C = zeros(nd, nz)

for t = 1:T-1
	H[u_idx[t], u_idx[t]] = R[t]
	H[x_idx[t], x_idx[t]] = Q[t+1]

	C[n_idx[t], u_idx[t]] = B[t]
	C[n_idx[t], x_idx[t]] = P[t]
	t == 1 && continue
	C[n_idx[t], x_idx[t-1]] = A[t]
end

# build kkt system
J = sparse([H C'; C -1.0e-16 * I])
r = rand(nz + nd)

# QDLDL
solver_ldl = ldl_solver(J)
Δldl = zeros(nz + nd)
@benchmark linear_solve!($solver_ldl, $Δldl, $J, $r)

# custom
# objective inverse
Q̃a = [inv(Qa[t]) for t = 1:T]
Q̃b = [inv(Qb[t]) for t = 1:T]

Q̃ = [inv(Q[t]) for t = 1:T]
R̃ = [inv(R[t]) for t = 1:T-1]
H̃ = zeros(nz, nz)

for t = 1:T-1
	H̃[u_idx[t], u_idx[t]] = R̃[t]
	H̃[x_idx[t], x_idx[t]] = Q̃[t+1]
end

norm(H̃ - inv(H))


###
s = newton_structure_solver(nq, m, T, Qa, R)
for t = 1:T
	s.Q̃a[t] = Q̃a[t]
	s.Q̃b[t] = Q̃b[t]

	t == T && continue
	s.Aa[t] = Aa[t]
	s.Ab[t] = Ab[t]
	s.Ac[t] = Ac[t]
	s.Ba[t] = Ba[t]

	s.R̃a[t] = R̃[t]
end

computeYs!(s.Yiia, s.Yiib, s.Yiic, s.Yiid, s.Yija, s.Yijb, s.Yijc, s.Yijd,
	s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a,
	s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.H)

for t = 1:T-1
	@show t
	@show norm(Yiia[t] - s.Yiia[t])
	@show norm(Yiib[t] - s.Yiib[t])
	@show norm(Yiic[t] - s.Yiic[t])
	@show norm(Yiid[t] - s.Yiid[t])

	@show norm(Yija[t] - s.Yija[t])
	@show norm(Yijb[t] - s.Yijb[t])
	@show norm(Yijc[t] - s.Yijc[t])
	@show norm(Yijd[t] - s.Yijd[t])

	@show norm(Yiiav[t] - s.Yiia[t])
	@show norm(Yiibv[t] - s.Yiib[t])
	@show norm(Yiicv[t] - s.Yiic[t])
	@show norm(Yiidv[t] - s.Yiid[t])

	@show norm(Yijav[t] - s.Yija[t])
	@show norm(Yijbv[t] - s.Yijb[t])
	@show norm(Yijcv[t] - s.Yijc[t])
	@show norm(Yijdv[t] - s.Yijd[t])
end

###

Ys = zeros(n * (T - 1), n * (T - 1))
Ysii = [zeros(n, n) for t = 1:T-1]
Ysij = [zeros(n, n) for t = 1:T-1]
tmp_nm = zeros(n, m)
tmp_nn = zeros(n, n)
tmp_nn2 = zeros(n, n)

Yiiav = [view(Ysii[t], collect(1:nq), collect(1:nq)) for t = 1:T-1]
Yiibv = [view(Ysii[t], collect(1:nq), collect(nq .+ (1:nq))) for t = 1:T-1]
Yiicv = [view(Ysii[t], collect(nq .+ (1:nq)), collect(1:nq)) for t = 1:T-1]
Yiidv = [view(Ysii[t], collect(nq .+ (1:nq)), collect(nq .+ (1:nq))) for t = 1:T-1]

Yijav = [view(Ysij[t], collect(1:nq), collect(1:nq)) for t = 1:T-1]
Yijbv = [view(Ysij[t], collect(1:nq), collect(nq .+ (1:nq))) for t = 1:T-1]
Yijcv = [view(Ysij[t], collect(nq .+ (1:nq)), collect(1:nq)) for t = 1:T-1]
Yijdv = [view(Ysij[t], collect(nq .+ (1:nq)), collect(nq .+ (1:nq))) for t = 1:T-1]

# tmp_q_nm = zeros(nq, m)
# tmp_q_nn = zeros(nq, nq)
# tmp_q_nn2 = zeros(nq, nq)

Yiia = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:T-1]
Yiib = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:T-1]
Yiic = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:T-1]
Yiid = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:T-1]

Yija = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:T-1]
Yijb = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:T-1]
Yijc = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:T-1]
Yijd = [SMatrix{nq,nq}(zeros(nq,nq)) for t = 1:T-1]

tmp_q_nm = SMatrix{nq,m}(zeros(nq, m))
tmp_q_nn = SMatrix{nq,nq}(zeros(nq, nq))
tmp_q_nn2 = SMatrix{nq,nq}(zeros(nq, nq))

function computeYs!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃, tmp_q_nn, tmp_q_nn2, tmp_q_nm, T)
	for t = 1:T-1
		if t == 1
			# Yiia = Ac[t] * Q̃b[t+1] Ac[t]' + Ba[t] * R̃[t] * Ba[t]'

			tmp_q_nn = Ac[t] * Q̃b[t+1]
			Yiia[t] = tmp_q_nn * transpose(Ac[t])

			tmp_q_nm = Ba[t] * R̃[t]
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

			tmp_q_nm = Ba[t] * R̃[t]
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
Yiia[1]
computeYs!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃, tmp_q_nn, tmp_q_nn2, tmp_q_nm, T)
@benchmark computeYs!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃, tmp_q_nn, tmp_q_nn2, tmp_q_nm, T)
@code_warntype computeYs!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃, tmp_q_nn, tmp_q_nn2, tmp_q_nm, T)

function update_Y!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Yiiav, Yiibv, Yiicv, Yiidv, Yijav, Yijbv, Yijcv, Yijdv, T)
	for t = 1:T-1
		Yiiav[t] .= Yiia[t]
		Yiibv[t] .= Yiib[t]
		Yiicv[t] .= Yiic[t]
		Yiidv[t] .= Yiid[t]

		Yijav[t] .= Yija[t]
		Yijbv[t] .= Yijb[t]
		Yijcv[t] .= Yijc[t]
		Yijdv[t] .= Yijd[t]
	end
	nothing
end

update_Y!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Yiiav, Yiibv, Yiicv, Yiidv, Yijav, Yijbv, Yijcv, Yijdv, T)
@benchmark update_Y!($Yiia, $Yiib, $Yiic, $Yiid, $Yija, $Yijb, $Yijc, $Yijd, $Yiiav, $Yiibv, $Yiicv, $Yiidv, $Yijav, $Yijbv, $Yijcv, $Yijdv, $T)
@code_warntype update_Y!(Yiia, Yiib, Yiic, Yiid, Yija, Yijb, Yijc, Yijd, Yiiav, Yiibv, Yiicv, Yiidv, Yijav, Yijbv, Yijcv, Yijdv, T)

idx_n = n_idx
for t = 1:T-1
	Ys[idx_n[t], idx_n[t]] = Ysii[t]

	t == T-1 && continue

	Ys[idx_n[t], idx_n[t+1]] = Ysij[t]
	Ys[idx_n[t+1], idx_n[t]] = Ysij[t]'
end

norm(Ys - C * H̃ * C')#[1:2n, 1:2n])

###
L = zeros(n * (T - 1), n * (T - 1))
Liis = [LowerTriangular(SMatrix{n,n}(zeros(n, n))) for t = 1:T-1]
Ljis = [SMatrix{n,n}(zeros(n, n)) for t = 1:T-1]

Yiis = [SMatrix{n,n}(Hermitian(Ysii[t])) for t = 1:T-1]
Yijs = [SMatrix{n,n}(Ysij[t]) for t = 1:T-1]

tmp_nn = SMatrix{n,n}(zeros(n,n))
tmp_nn2 = SMatrix{n,n}(zeros(n,n))

function computeLs!(Lii, Lji, Yii, Yij, tmp_nn, tmp_nn2, tmp_nn_a, T)
	for t = 1:T-1
		# @show t
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

@benchmark computeLs!($Liis, $Ljis, $Yiis, $Yijs, $tmp_nn, $tmp_nn2, $tmp_nn_a, $T)
@code_warntype computeLs!(Liis, Ljis, Yiis, Yijs, tmp_nn, tmp_nn2, tmp_nn_a, T)

for t = 1:T-1
	L[n_idx[t], n_idx[t]] = LowerTriangular(Liis[t])

	t == T-1 && continue

	L[n_idx[t+1], n_idx[t]] = transpose(Ljis[t])
end
norm(cholesky(Hermitian(Ys)).L - L)

###

#solve system
Δz = zeros(nz)
Δν = zeros(nd)

rlag = view(r, 1:nz)
_rlagu = [view(rlag, u_idx[t]) for t = 1:T-1]
_rlagx = [view(rlag, x_idx[t]) for t = 1:T-1]
_rlagqa = [view(rlag, qa_idx[t]) for t = 1:T-1]
_rlagqb = [view(rlag, qb_idx[t]) for t = 1:T-1]

rlagu = [SVector{m}(zeros(m)) for t = 1:T-1]
rlagqa = [SVector{nq}(zeros(nq)) for t = 1:T-1]
rlagqb = [SVector{nq}(zeros(nq)) for t = 1:T-1]

for t = 1:T-1
	rlagu[t] = _rlagu[t]
	rlagqa[t] = _rlagqa[t]
	rlagqb[t] = _rlagqb[t]
end


rdyn = view(r, nz .+ (1:nd))
_rdyn1 = [view(rdyn, n_idx[t][1:nq]) for t = 1:T-1]
_rdyn2 = [view(rdyn, n_idx[t][nq .+ (1:nq)]) for t = 1:T-1]

rdyn1 = [SVector{nq}(zeros(nq)) for t = 1:T-1]
rdyn2 = [SVector{nq}(zeros(nq)) for t = 1:T-1]

for t = 1:T-1
	rdyn1[t] = _rdyn1[t]
	rdyn2[t] = _rdyn2[t]
end

β = -rdyn + C * H̃ * rlag

βn = [SVector{n}(zeros(n)) for t = 1:T-1]
βd = [SVector{nq}(zeros(nq)) for t = 1:T-1]
βe = [SVector{nq}(zeros(nq)) for t = 1:T-1]
tmp_q_n = SVector{nq}(zeros(nq))

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

@benchmark update_β!($βn, $βd, $βe, $rlagu, $rlagqa, $rlagqb, $rdyn1, $rdyn2, $Aa, $Ab, $Ac, $Ba, $Q̃a, $Q̃b, $R̃, $T)
@code_warntype update_β!(βn, βd, βe, rlagu, rlagqa, rlagqb, rdyn1, rdyn2, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃, T)
norm(β - vcat([[βd[t]; βe[t]] for t = 1:T-1]...))
norm((β - vcat(βn...)))


# forward subsitution
y = [zeros(n) for t = 1:T-1]
b = [SVector{n}(view(β, idx_n[t])) for t = 1:T-1]

function forward_substitution!(y, Lii, Lji, b, T)
	for t = 1:T-1
		if t == 1
			# y[1] .= Lii[1] \ β[idx_n[1]]

			y[1] .= b[t]
			LAPACK.trtrs!('L', 'N', 'N', Array(Lii[t]), y[1])
		else
			# y[t] .= Lii[t] \ (β[idx_n[t]] - Lji[t - 1] * y[t - 1])

			y[t] .= b[t]
			mul!(y[t], transpose(Array(Lji[t - 1])), y[t - 1], -1.0, 1.0)
			LAPACK.trtrs!('L', 'N', 'N', Array(Lii[t]), y[t])
		end
	end
	nothing
end

@benchmark forward_substitution!($y, $Liis, $Ljis, $b, $T)
@code_warntype forward_substitution!(y, Liis, Ljis, b, T)
norm(vcat(y...) - L \ β, Inf)

ys = [@SVector zeros(n) for t = 1:T-1]
b = [SVector{n}(view(β, idx_n[t])) for t = 1:T-1]

function forward_substitution_s!(y, Lii, Lji, b, T)
	for t = 1:T-1
		if t == 1
			# y[1] .= Lii[1] \ β[n_idx[1]]

			y[1] = Lii[1] \ b[1]
		else
			# y[t] .= Lii[t] \ (β[n_idx[t]] - transpose(Lji[t - 1]) * y[t - 1])
			# y[t] = b[t]
			# b[t] -= transpose(Lji[t - 1]) * y[t - 1]
			# y[t] = Lii[t] \ b[t]
			y[t] = Lii[t] \ (b[t] - transpose(Lji[t - 1]) * y[t - 1])
		end
	end
	nothing
end

@benchmark forward_substitution_s!($ys, $Liis, $Ljis, $b, $T)
@code_warntype forward_substitution_s!(ys, Liis, Ljis, b, T)
norm(vcat(ys...) - L \ β, Inf)

# backward substitution
x_vec = zeros(n * (T-1))
x = [view(x_vec, idx_n[t]) for t = 1:T-1]
function backward_substitution!(x, Lii, Lji, y, T)
	for t = T-1:-1:1
		if t == T-1
			# x[t] = LowerTriangular(Lii[t])' \ y[t]

			x[t] .= y[t]
			LAPACK.trtrs!('L', 'T', 'N', Array(Lii[t]), x[t])
		else
			# x[t] = LowerTriangular(Lii[t])' \ (y[t] - Lji[t]' * x[t+1])

			x[t] .= y[t]
			mul!(x[t], Array(Lji[t]), x[t+1], -1.0, 1.0)
			LAPACK.trtrs!('L', 'T', 'N', Array(Lii[t]), x[t])
		end
	end
	nothing
end

@benchmark backward_substitution!($x, $Liis, $Ljis, $y, $T)
@code_warntype backward_substitution!(x, Liis, Ljis, y, T)

norm(vcat(x...) - Ys \ β, Inf)
norm(vcat(x...) - L' \ (L \ β), Inf)
norm(x_vec - L' \ (L \ β), Inf)
norm(x_vec - L' \ (L \ β), Inf)

Δνn = [SVector{n}(zeros(n)) for t = 1:T-1]
function backward_substitution_s!(x, Lii, Lji, y, T)
	for t = T-1:-1:1
		if t == T-1
			# x[t] = LowerTriangular(Lii[t])' \ y[t]

			x[t] = transpose(Lii[t]) \ y[t]
		else
			# x[t] = LowerTriangular(Lii[t])' \ (y[t] - Lji[t] * x[t+1])

			# y[t] -= Lji[t] * x[t + 1]
			# x[t] = transpose(Lii[t]) \ y[t]
			x[t] = transpose(Lii[t]) \ (y[t] - Lji[t] * x[t + 1])
		end
	end
	nothing
end

@benchmark backward_substitution_s!($Δνn, $Liis, $Ljis, $ys, $T)
@code_warntype backward_substitution_s!(Δνn, Liis, Ljis, ys, T)

norm(vcat(Δνn...) - Ys \ β, Inf)
norm(vcat(x...) - L' \ (L \ β), Inf)
norm(x_vec - L' \ (L \ β), Inf)
norm(x_vec - L' \ (L \ β), Inf)

Δν = zeros(nd)
Δν .= Ys \ β

Δν = vcat(Δνn...)

Δνd = [SVector{nq}(view(Δν, n_idx[t][1:nq])) for t = 1:T-1]
Δνe = [SVector{nq}(view(Δν, n_idx[t][nq .+ (1:nq)])) for t = 1:T-1]

Δzu = [@SVector zeros(m) for t = 1:T-1]
Δzqa = [@SVector zeros(nq) for t = 1:T-1]
Δzqb = [@SVector zeros(nq) for t = 1:T-1]

# Δν .= β
# LAPACK.potrs!('L', L, Δν)
Δz .= H̃ * (rlag - C' * Δν)

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

@benchmark update_Δz!($Δzu, $Δzqa, $Δzqb, $Δνd, $Δνe, $Aa, $Ab, $Ac, $Ba, $Q̃a, $Q̃b, $R̃, $rlagu, $rlagqa, $rlagqb, $T)
@code_warntype update_Δz!(Δzu, Δzqa, Δzqb, Δνd, Δνe, Aa, Ab, Ac, Ba, Q̃a, Q̃b, R̃, rlagu, rlagqa, rlagqb, T)

norm((Δz - vcat([[Δzu[t]; Δzqa[t]; Δzqb[t]] for t = 1:T-1]...))[1:end])
