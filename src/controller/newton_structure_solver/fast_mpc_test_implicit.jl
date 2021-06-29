"""
Fast MPC (modified for implicit dynamics (f(xt, ut, xt+1) = 0))
https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

using BenchmarkTools
using InteractiveUtils

# dimensions
n = 20
m = 10
T = 25
nz = m * (T - 1) + n * (T - 1)
nd = n * (T - 1)

# indices
u_idx = [(t - 1) * (m + n) .+ (1:m) for t = 1:T-1]
x_idx = [(t - 1) * (m + n) + m .+ (1:n) for t = 1:T-1]
d_idx = [(t - 1) * n .+ (1:n) for t = 1:T-1]

# problem data
Q = [SMatrix{n, n}(Diagonal(0.1 * rand(n) + ones(n))) for t = 1:T]
R = [SMatrix{m, m}(Diagonal(0.1 * rand(m) + ones(m))) for t = 1:T-1]
A = [SMatrix{n, n}(Diagonal(0.1 * rand(n) + ones(n))) for t = 1:T-1]
B = [SMatrix{n, m}(rand(n, m)) for t = 1:T-1]
P = [SMatrix{n, n}(Diagonal(0.01 * rand(n) + ones(n))) for t = 1:T-1]

# direct problem
H = zeros(nz, nz)
C = zeros(nd, nz)

for t = 1:T-1
	H[u_idx[t], u_idx[t]] = R[t]
	H[x_idx[t], x_idx[t]] = Q[t+1]

	C[d_idx[t], u_idx[t]] = B[t]
	C[d_idx[t], x_idx[t]] = P[t]
	t == 1 && continue
	C[d_idx[t], x_idx[t-1]] = A[t]
end

# build kkt system
J = sparse([H C'; C zeros(nd, nd)])
r = rand(nz + nd)

# QDLDL
solver_ldl = ldl_solver(J)
Δldl = zeros(nz + nd)
@benchmark linear_solve!($solver_ldl, $Δldl, $J, $r)

# Custom

# objective inverse
Q̃ = [SMatrix{n, n}(inv(Q[t])) for t = 1:T]
R̃ = [SMatrix{m, m}(inv(R[t])) for t = 1:T-1]
H̃ = zeros(nz, nz)

for t = 1:T-1
	H̃[u_idx[t], u_idx[t]] = R̃[t]
	H̃[x_idx[t], x_idx[t]] = Q̃[t+1]
end

norm(H̃ - inv(H))

idx_n = [(t - 1) * n .+ (1:n) for t = 1:T-1]

Y = zeros(n * (T - 1), n * (T - 1))
Yii = [@SMatrix zeros(n, n) for t = 1:T-1]
Yij = [@SMatrix zeros(n, n) for t = 1:T-1]
tmp_nm = @SMatrix zeros(n, m)
tmp_nn = @SMatrix zeros(n, n)
tmp_nn2 = @SMatrix zeros(n, n)

function computeY!(Yii, Yij, A, B, P, Q̃, R̃, tmp_nn, tmp_nn2, tmp_nm, T)
	for t = 1:T-1
		if t == 1
			# Yii[t] = B[t] * R̃[t] * B[t]' + P[t] * Q̃[t+1] * P[t]'

			tmp_nn = P[t] * Q̃[t+1]
			Yii[t] = tmp_nn * transpose(P[t])
			tmp_nm = B[t] * R̃[t]
			tmp_nn = tmp_nm * transpose(B[t])
			Yii[t] += tmp_nn
		else
			# Yii[t] = A[t] * Q̃[t] * A[t]' + B[t] * R̃[t] * B[t]' + P[t] * Q̃[t+1] * P[t]'

			tmp_nn = A[t] * Q̃[t]
			Yii[t] = tmp_nn * transpose(A[t])
			tmp_nm = B[t] * R̃[t]
			tmp_nn = tmp_nm * transpose(B[t])
			Yii[t] += tmp_nn
			tmp_nn = P[t] * Q̃[t+1]
			tmp_nn2 = tmp_nn * transpose(P[t])
			Yii[t] += tmp_nn2
		end

		t == T-1 && continue
		# Yij[t] = P[t] * Q̃[t+1] * A[t+1]'

		tmp_nn = P[t] * Q̃[t+1]
		Yij[t] = tmp_nn * transpose(A[t+1])
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

@benchmark computeL!($Lii, $Lji, $Yii, $Yij, $T)
@code_warntype computeL!(Lii, Lji, Yii, Yij, T)

computeL!(Lii, Lji, Yii, Yij, T)

for t = 1:T-1
	L[idx_n[t], idx_n[t]] = LowerTriangular(Lii[t])

	t == T-1 && continue

	L[idx_n[t+1], idx_n[t]] = transpose(Lji[t])
end

norm(cholesky(Hermitian(Y)).L - L)

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
			# # y[t] .= Lii[t] \ (β[idx_n[t]] - Lji[t - 1] * y[t - 1])

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
