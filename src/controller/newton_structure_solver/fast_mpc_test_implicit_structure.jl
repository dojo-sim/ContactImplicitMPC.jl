"""
Fast MPC (modified for implicit dynamics (f(xt, ut, xt+1) = 0))
https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

using BenchmarkTools
using InteractiveUtils

# dimensions
nq = 3
n = 2 * nq
m = 2
T = 3
nz = m * (T - 1) + n * (T - 1)
nd = n * (T - 1)

# indices
u_idx = [(t - 1) * (m + n) .+ (1:m) for t = 1:T-1]
x_idx = [(t - 1) * (m + n) + m .+ (1:n) for t = 1:T-1]
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

idx_n = [(t - 1) * n .+ (1:n) for t = 1:T-1]

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

			tmp_nn = P[t] * Q̃[t+1]
			# mul!(tmp_nn, P[t], Q̃[t+1])
			# tmp_nn .= P[t]
			# rmul!(tmp_nn, Q̃[t+1])
			Yii[t] = tmp_nn * transpose(P[t])
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

		tmp_nn = P[t] * Q̃[t+1]
		# mul!(tmp_nn, P[t], Q̃[t+1])
		# tmp_nn .= P[t]
		# rmul!(tmp_nn, Q̃[t+1])

		Yij[t] = tmp_nn * transpose(A[t+1])
		# mul!(Yij[t], tmp_nn, transpose(A[t+1]))
	end
	nothing
end

@benchmark computeY!($Yii, $Yij, $A, $B, $P, $Q̃, $R̃, $tmp_nn, $tmp_nn2, $tmp_nm, $T)
@code_warntype computeY!(Yii, Yij, A, B, P, Q̃, R̃, tmp_nn, tmp_nn, tmp_nm, T)
for t = 1:T-1
	Y[idx_n[t], idx_n[t]] = Yii[t]

	t == T-1 && continue

	Y[idx_n[t], idx_n[t+1]] = Yij[t]
	Y[idx_n[t+1], idx_n[t]] = Yij[t]'
end

norm((Y - C * H̃ * C'))#[1:2n, 1:2n])

### structured Y
Ys = zeros(n * (T - 1), n * (T - 1))
Ysii = [zeros(n, n) for t = 1:T-1]
Ysij = [zeros(n, n) for t = 1:T-1]
tmp_nm = zeros(n, m)
tmp_nn = zeros(n, n)
tmp_nn2 = zeros(n, n)

Yiiav = [view(Ysii[t], 1:nq, 1:nq) for t = 1:T-1]
Yiibv = [view(Ysii[t], 1:nq, nq .+ (1:nq)) for t = 1:T-1]
Yiicv = [view(Ysii[t], nq .+ (1:nq), 1:nq) for t = 1:T-1]
Yiidv = [view(Ysii[t], nq .+ (1:nq), nq .+ (1:nq)) for t = 1:T-1]

Yijav = [view(Ysij[t], 1:nq, 1:nq) for t = 1:T-1]
Yijbv = [view(Ysij[t], 1:nq, nq .+ (1:nq)) for t = 1:T-1]
Yijcv = [view(Ysij[t], nq .+ (1:nq), 1:nq) for t = 1:T-1]
Yijdv = [view(Ysij[t], nq .+ (1:nq), nq .+ (1:nq)) for t = 1:T-1]
#
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
idx_n
for t = 1:T-1
	Ys[idx_n[t], idx_n[t]] = Ysii[t]

	t == T-1 && continue

	Ys[idx_n[t], idx_n[t+1]] = Ysij[t]
	Ys[idx_n[t+1], idx_n[t]] = Ysij[t]'
end

norm(Ys - Y)
norm((Ys - C * H̃ * C'))#[1:2n, 1:2n])
#####

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

function test_cholesky1(A)
	B = copy(A)
	cholesky(Hermitian(B)).L
end

function test_cholesky2(A)
	B = copy(A)
	LAPACK.potrf!('L', B)
end

cholesky(Hermitian(Yii[1])).L
@benchmark test_cholesky1($Yii[1])

_Y = Array(Yii[1])
@benchmark test_cholesky2($_Y)

@benchmark computeL!($Lii, $Lji, $Yii, $Yij, $T)
@code_warntype computeL!(Lii, Lji, Yii, Yij, T)

computeL!(Lii, Lji, Yii, Yij, T)

for t = 1:T-1
	L[idx_n[t], idx_n[t]] = LowerTriangular(Lii[t])

	t == T-1 && continue

	L[idx_n[t+1], idx_n[t]] = transpose(Lji[t])
end

norm(cholesky(Hermitian(Ys)).L - L)

###
L = zeros(n * (T - 1), n * (T - 1))
Liis = [LowerTriangular(SMatrix{n,n}(zeros(n, n))) for t = 1:T-1]
Ljis = [SMatrix{n,n}(zeros(n, n)) for t = 1:T-1]

Yiis = [SMatrix{n,n}(Hermitian(Yiis[t])) for t = 1:T-1]
Yijs = [SMatrix{n,n}(Yijs[t]) for t = 1:T-1]

tmp_nn = SMatrix{n,n}(zeros(n,n))
tmp_nn2 = SMatrix{n,n}(zeros(n,n))
tmp_nn_a = zeros(n,n)
function computeLs!(Lii, Lji, Yii, Yij, tmp_nn, tmp_nn2, tmp_nn_a, T)
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

@benchmark computeLs!($Liis, $Ljis, $Yiis, $Yijs, $tmp_nn, $tmp_nn2, $tmp_nn_a, $T)
@code_warntype computeLs!(Liis, Ljis, Yiis, Yijs, tmp_nn, tmp_nn2, tmp_nn_a, T)

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

rd = r[1:nz]
rp = r[nz .+ (1:nd)]
β = -rp + C * H̃ * rd
b = [view(β, idx_n[t]) for t = 1:T-1]
Δν .= Y \ β
# Δν .= β
# LAPACK.potrs!('L', L, Δν)
Δz .= H̃ * (rd - C' * Δν)
norm(Δldl - [Δz; Δν])

# forward subsitution
y = [zeros(n) for t = 1:T-1]

function forward_substitution!(y, Lii, Lji, b, T)
	for t = 1:T-1
		if t == 1
			# y[1] .= Lii[1] \ β[idx_n[1]]

			y[1] .= b[t]
			LAPACK.trtrs!('L', 'N', 'N', Lii[t], y[1])
		else
			# y[t] .= Lii[t] \ (β[idx_n[t]] - Lji[t - 1] * y[t - 1])

			y[t] .= b[t]
			mul!(y[t], transpose(Lji[t - 1]), y[t - 1], -1.0, 1.0)
			LAPACK.trtrs!('L', 'N', 'N', Lii[t], y[t])
		end
	end
	nothing
end

@benchmark forward_substitution!($y, $Lii, $Lji, $b, $T)
@code_warntype forward_substitution!(y, Lii, Lji, b, T)
norm(vcat(y...) - L \ β, Inf)

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

# backward substitution
x_vec = zeros(n * (T-1))
x = [view(x_vec, idx_n[t]) for t = 1:T-1]
function backward_substitution!(x, Lii, Lji, y, T)
	for t = T-1:-1:1
		if t == T-1
			# x[t] = LowerTriangular(Lii[t])' \ y[t]

			x[t] .= y[t]
			LAPACK.trtrs!('L', 'T', 'N', Lii[t], x[t])
		else
			# x[t] = LowerTriangular(Lii[t])' \ (y[t] - Lji[t]' * x[t+1])

			x[t] .= y[t]
			mul!(x[t], Lji[t], x[t+1], -1.0, 1.0)
			LAPACK.trtrs!('L', 'T', 'N', Lii[t], x[t])
		end
	end
	nothing
end

@benchmark backward_substitution!($x, $Lii, $Lji, $y, $T)
@code_warntype backward_substitution!(x, Lii, Lji, y, T)

norm(vcat(x...) - Y \ β, Inf)
norm(vcat(x...) - L' \ (L \ β), Inf)
norm(x_vec - L' \ (L \ β), Inf)
norm(x_vec - L' \ (L \ β), Inf)

xs = [SVector{n}(zeros(n)) for t = 1:T-1]
function backward_substitution_s!(x, Lii, Lji, y, T)
	for t = T-1:-1:1
		if t == T-1
			# x[t] = LowerTriangular(Lii[t])' \ y[t]

			x[t] = transpose(Lii[t]) \ y[t]
		else
			# x[t] = LowerTriangular(Lii[t])' \ (y[t] - Lji[t]' * x[t+1])

			y[t] -= transpose(Lji[t]) * x[t+1]
			x[t] = transpose(Lii[t]) \ y[t]
		end
	end
	nothing
end

@benchmark backward_substitution_s!($xs, $Liis, $Ljis, $ys, $T)
@code_warntype backward_substitution!(x, Lii, Lji, y, T)
