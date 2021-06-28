"""
Fast MPC (modified for implicit dynamics (f(xt, ut, xt+1) = 0))
https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

using BenchmarkTools
using InteractiveUtils

# dimensions
nq = 10
n = 2 * nq
m = 10
T = 20
nz = m * (T - 1) + n * (T - 1)
nd = n * (T - 1)

# indices
u_idx = [(t - 1) * (m + n) .+ (1:m) for t = 1:T-1]
x_idx = [(t - 1) * (m + n) + m .+ (1:n) for t = 1:T-1]
qa_idx = [x_idx[t][1:nq] for t = 1:T-1]
qb_idx = [x_idx[t][nq .+ (1:nq)] for t = 1:T-1]
n_idx = [(t - 1) * n .+ (1:n) for t = 1:T-1]

# problem data
_Aa = [Array(-1.0 * Diagonal(0.1 * rand(nq) + ones(nq))) for t = 1:T-1]
_Ab = [Array(-1.0 * Diagonal(0.1 * rand(nq) + ones(nq))) for t = 1:T-1]
_Ac = [Array(Diagonal(ones(nq))) for t = 1:T-1]
_Ba = [-1.0 * rand(nq, m) for t = 1:T-1]

Aa = [SMatrix{nq, nq}(zeros(nq, nq)) for t = 1:T-1]
Ab = [SMatrix{nq, nq}(zeros(nq, nq)) for t = 1:T-1]
Ac = [SMatrix{nq, nq}(zeros(nq, nq)) for t = 1:T-1]
Ba = [SMatrix{nq,m}(zeros(nq, m)) for t = 1:T-1]

Qa = [Diagonal(0.1 * rand(nq) + ones(nq)) for t = 1:T]
Qb = [Diagonal(0.1 * rand(nq) + ones(nq)) for t = 1:T]

function static_dynamics_jacobians!(Aas, Abs, Acs, Aa, Ab, Ac, T)
	for t = 1:T-1
		Aas[t] = Aa[t]
		Abs[t] = Ab[t]
		Acs[t] = Ac[t]
	end
end

@benchmark static_dynamics_jacobians!($Aa, $Ab, $Ac, $_Aa, $_Ab, $_Ac, $T)
@code_warntype static_dynamics_jacobians!(Aa, Ab, Ac, _Aa, _Ab, _Ac, T)


Q = [Diagonal(SVector{n}(Array(diag(cat(Qa[t], Qb[t], dims=(1,2)))))) for t = 1:T]
R = [Diagonal(SVector{m}(0.1 * rand(m) + ones(m))) for t = 1:T-1]
A = [SMatrix{n,n}([Aa[t] Ab[t]; zeros(nq, nq) -I]) for t = 1:T-1]
B = [SMatrix{n,m}([Ba[t]; zeros(nq, m)]) for t = 1:T-1]
P = [SMatrix{n,n}([zeros(nq, nq) Ac[t]; I zeros(nq, nq)]) for t = 1:T-1]

# direct problem
H = spzeros(nz, nz)
for t = 1:T-1
	H[u_idx[t], u_idx[t]] = R[t]
	H[qa_idx[t], qa_idx[t]] = Qa[t+1]
	H[qb_idx[t], qb_idx[t]] = Qb[t+1]
end

# objective inverse
Q̃a = [inv(Qa[t]) for t = 1:T]
Q̃b = [inv(Qb[t]) for t = 1:T]

Q̃ = [inv(Q[t]) for t = 1:T]
R̃ = [inv(R[t]) for t = 1:T-1]
H̃ = zeros(nz, nz)

for t = 1:T-1
	H̃[u_idx[t], u_idx[t]] = R̃[t]
	H̃[qa_idx[t], qa_idx[t]] = Q̃a[t+1]
	H̃[qb_idx[t], qb_idx[t]] = Q̃b[t+1]
end

norm(H̃ - inv(Array(H)))

C = spzeros(nd, nz)
Cdu = [view(C, n_idx[t][1:nq], u_idx[t]) for t = 1:T-1]
Cdqa = [t == 1 ? view(C, n_idx[t][1:nq], 1:0) : view(C, n_idx[t][1:nq], x_idx[t-1][1:nq]) for t = 1:T-1]
Cdqb = [t == 1 ? view(C, n_idx[t][1:nq], 1:0) : view(C, n_idx[t][1:nq], x_idx[t-1][nq .+ (1:nq)]) for t = 1:T-1]
Cdqc = [view(C, n_idx[t][1:nq], x_idx[t][nq .+ (1:nq)]) for t = 1:T-1]

Ceqb = [t == 1 ? view(C, n_idx[t][nq .+ (1:nq)], 1:0) : view(C, n_idx[t][nq .+ (1:nq)], x_idx[t-1][nq .+ (1:nq)]) for t = 1:T-1]
Ceqa = [view(C, n_idx[t][nq .+ (1:nq)], x_idx[t][1:nq]) for t = 1:T-1]

function build_C!(Cdu, Cdqa, Cdqb, Cdqc, Ceqa, Ceqb, Aa, Ab, Ac, Ba, T)
	for t = 1:T-1
		Cdu[t] .= Ba[t]
		Cdqc[t] .= Ac[t]

		Ceqa[t] .= Diagonal(ones(nq))
		t == 1 && continue
		Cdqa[t] .= Aa[t]
		Cdqb[t] .= Ab[t]

		Ceqb[t] .= Diagonal(-1.0 * ones(nq))
	end
end

@benchmark build_C!($Cdu, $Cdqa, $Cdqb, $Cdqc, $Ceqa, $Ceqb, $_Aa, $_Ab, $_Ac, $_Ba, $T)

C1 = spzeros(nd, nz)
for t = 1:T-1
	C1[n_idx[t], u_idx[t]] = B[t]
	C1[n_idx[t], x_idx[t]] = P[t]
	t == 1 && continue
	C1[n_idx[t], x_idx[t-1]] = A[t]
end

t = T-1
C[n_idx[t][nq .+ (1:nq)], x_idx[t-1][nq .+ (1:nq)]] #- Ab[t]
C1[n_idx[t][nq .+ (1:nq)], x_idx[t-1][nq .+ (1:nq)]]

C = C1
norm(C - C1)

# build kkt system
J = sparse([H C'; C -1.0e-16 * I])
r = rand(nz + nd)

# custom


### structured Y
Ys = zeros(n * (T - 1), n * (T - 1))
Ysii = [zeros(n, n) for t = 1:T-1]
Ysij = [zeros(n, n) for t = 1:T-1]

Yiiav = [view(Ysii[t], 1:nq, 1:nq) for t = 1:T-1]
Yiibv = [view(Ysii[t], 1:nq, nq .+ (1:nq)) for t = 1:T-1]
Yiicv = [view(Ysii[t], nq .+ (1:nq), 1:nq) for t = 1:T-1]
Yiidv = [view(Ysii[t], nq .+ (1:nq), nq .+ (1:nq)) for t = 1:T-1]

Yijav = [view(Ysij[t], 1:nq, 1:nq) for t = 1:T-1]
Yijbv = [view(Ysij[t], 1:nq, nq .+ (1:nq)) for t = 1:T-1]
Yijcv = [view(Ysij[t], nq .+ (1:nq), 1:nq) for t = 1:T-1]
Yijdv = [view(Ysij[t], nq .+ (1:nq), nq .+ (1:nq)) for t = 1:T-1]

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

for t = 1:T-1
	Ys[idx_n[t], idx_n[t]] = Ysii[t]

	t == T-1 && continue

	Ys[idx_n[t], idx_n[t+1]] = Ysij[t]
	Ys[idx_n[t+1], idx_n[t]] = Ysij[t]'
end

norm((Ys - C * H̃ * C'))#[1:2n, 1:2n])
#####

L = zeros(n * (T - 1), n * (T - 1))
Liis = [LowerTriangular(SMatrix{n,n}(zeros(n, n))) for t = 1:T-1]
Ljis = [SMatrix{n,n}(zeros(n, n)) for t = 1:T-1]

Yiis = [SMatrix{n,n}(Hermitian(Yii[t])) for t = 1:T-1]
Yijs = [SMatrix{n,n}(Yij[t]) for t = 1:T-1]

tmp_nn = SMatrix{n,n}(zeros(n,n))
tmp_nn2 = SMatrix{n,n}(zeros(n,n))

function computeLs!(Lii, Lji, Yii, Yij, tmp_nn, tmp_nn2, T)
	for t = 1:T-1
		# @show t
		if t == 1
			# Lii[t] = cholesky(Hermitian(Yii[t])).L

			Lii[t] = cholesky(Yii[t]).L
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

@benchmark computeLs!($Liis, $Ljis, $Yiis, $Yijs, $tmp_nn, $tmp_nn2, $T)
@code_warntype computeLs!(Liis, Ljis, Yiis, Yijs, tmp_nn, tmp_nn2, T)

for t = 1:T-1
	L[idx_n[t], idx_n[t]] = LowerTriangular(Liis[t])

	t == T-1 && continue

	L[idx_n[t+1], idx_n[t]] = transpose(Ljis[t])
end
norm(cholesky(Hermitian(Y)).L - L)

###

# solve system
Δz = zeros(nz)
Δν = zeros(nd)

rd = view(r, 1:nz)
rdu = [view(rd, u_idx[t]) for t = 1:T-1]
rdx = [view(rd, x_idx[t]) for t = 1:T-1]
rdqa = [view(rd, qa_idx[t]) for t = 1:T-1]
rdqb = [view(rd, qb_idx[t]) for t = 1:T-1]

rp = view(r, nz .+ (1:nd))
rpd = [view(rp, n_idx[t]) for t = 1:T-1]
rpdd = [view(rpd[t], 1:nq) for t = 1:T-1]
rpde = [view(rpd[t], nq .+ (1:nq)) for t = 1:T-1]

n_idx

β = -rp + C * H̃ * rd

βn = [zeros(n) for t = 1:T-1]
βnd = [view(βn[t], 1:nq) for t = 1:T-1]
βne = [view(βn[t], nq .+ (1:nq)) for t = 1:T-1]

for t = 1:T-1
	if t == 1
		βnd[1] .= -1.0 * rpdd[1] + Ba[1] * R̃[1] * rdu[1] + Ac[1] * Q̃b[2] * rdqb[1]
		βne[1] .= -1.0 * rpde[1] + Q̃a[2] * rdqa[1]
	else
		βnd[t] .= -1.0 * rpdd[t] + Ba[t] * R̃[t] * rdu[t] + Ac[t] * Q̃b[t+1] * rdqb[t] + Aa[t] * Q̃a[t] * rdqa[t-1] + Ab[t] * Q̃b[t] * rdqb[t-1]
		βne[t] .= -1.0 * rpde[t] + Q̃a[t+1] * rdqa[t] - Q̃b[t] * rdqb[t-1]
	end
end

norm((β - vcat(βn...)))


b = [view(β, idx_n[t]) for t = 1:T-1]
Δν .= Y \ β

Δνd = [view(Δν, n_idx[t][1:nq]) for t = 1:T-1]
Δνe = [view(Δν, n_idx[t][nq .+ (1:nq)]) for t = 1:T-1]

Δzu = [zeros(m) for t = 1:T-1]
Δzqa = [zeros(nq) for t = 1:T-1]
Δzqb = [zeros(nq) for t = 1:T-1]

# Δν .= β
# LAPACK.potrs!('L', L, Δν)
Δz .= H̃ * (rd - C' * Δν)

for t = 1:T-1
	Δzu[t] .= R̃[t] * (rdu[t] - transpose(Ba[t]) * Δνd[t])
	if t < T-1
		Δzqa[t] .= Q̃a[t+1] * (rdqa[t] - Δνe[t] - transpose(Aa[t+1]) * Δνd[t+1])
		Δzqb[t] .= Q̃b[t+1] * (rdqb[t] - transpose(Ac[t]) * Δνd[t] - transpose(Ab[t+1]) * Δνd[t+1] + Δνe[t+1])
	else
		Δzqa[t] .= Q̃a[t+1] * (rdqa[t] - Δνe[t])
		Δzqb[t] .= Q̃b[t+1] * (rdqb[t] - transpose(Ac[t]) * Δνd[t])
	end
end

norm((Δz - vcat([[Δzu[t]; Δzqa[t]; Δzqb[t]] for t = 1:T-1]...))[1:end])

ys = [@SVector zeros(n) for t = 1:T-1]
bs = [SVector{n}(β[idx_n[t]]) for t = 1:T-1]

function forward_substitution_s!(y, Lii, Lji, b, T)
	for t = 1:T-1
		if t == 1
			# y[1] .= Lii[1] \ β[idx_n[1]]

			y[1] = Lii[1] \ b[1]
		else
			# y[t] .= Lii[t] \ (β[idx_n[t]] - Lji[t - 1] * y[t - 1])

			b[t] -= Lji[t - 1] * y[t - 1]
			y[t] = Lii[t] \ b[t]
		end
	end
	nothing
end

@benchmark forward_substitution_s!($ys, $Liis, $Ljis, $bs, $T)
@code_warntype forward_substitution_s!(ys, Liis, Ljis, bs, T)
norm(vcat(y...) - L \ β, Inf)

xs = [SVector{n}(zeros(n)) for t = 1:T-1]
function backward_substitution_s!(x, Lii, Lji, y, T)
	for t = T-1:-1:1
		if t == T-1
			# x[t] = LowerTriangular(Lii[t])' \ y[t]

			x[t] = transpose(Lii[t]) \ y[t]
		else
			# x[t] = LowerTriangular(Lii[t])' \ (y[t] - Lji[t]' * x[t+1])

			y[t] -= transpose(Lji[t]) * x[t + 1]
			x[t] = transpose(Lii[t]) \ y[t]
		end
	end
	nothing
end

@benchmark backward_substitution_s!($xs, $Liis, $Ljis, $ys, $T)
@code_warntype backward_substitution!(x, Lii, Lji, y, T)
