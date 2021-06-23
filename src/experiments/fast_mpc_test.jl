"""
Fast MPC
https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

using BenchmarkTools

# dimensions
n = 5
m = 3
T = 25
nz = m * (T - 1) + n * (T - 1)
nd = n * (T - 1)

# indices
u_idx = [(t - 1) * (m + n) .+ (1:m) for t = 1:T-1]
x_idx = [(t - 1) * (m + n) + m .+ (1:n) for t = 1:T-1]
d_idx = [(t - 1) * n .+ (1:n) for t = 1:T-1]

# problem data
Q = [Diagonal(0.1 * rand(n) + ones(n)) for t = 1:T]
R = [Diagonal(0.1 * rand(m) + ones(m)) for t = 1:T-1]
A = [Diagonal(0.1 * rand(n) + ones(n)) for t = 1:T-1]
B = [rand(n, m) for t = 1:T-1]

# direct problem
H = zeros(nz, nz)
C = zeros(nd, nz)

for t = 1:T-1
	H[u_idx[t], u_idx[t]] = R[t]
	H[x_idx[t], x_idx[t]] = Q[t+1]

	C[d_idx[t], u_idx[t]] = -B[t]
	C[d_idx[t], x_idx[t]] = Diagonal(ones(n))
	t == 1 && continue
	C[d_idx[t], x_idx[t-1]] = -A[t]
end

# build kkt system
J = sparse([H C'; C zeros(nd, nd)])
r = rand(nz + nd)

# QDLDL
solver_ldl = ldl_solver(J)
Δldl = zeros(nz + nd)
@benchmark linear_solve!($solver_ldl, $Δldl, $J, $r)

# # LU
# solver_lu = lu_solver(J)
# Δlu = zeros(nz + nd)
# @benchmark linear_solve!($solver_lu, $Δlu, $J, $r)

# Custom

# objective inverse
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
Yii = [zeros(n, n) for t = 1:T-1]
Yij = [zeros(n, n) for t = 1:T-1]
tmp_nm = zeros(n, m)
tmp_nn = zeros(n, n)

for t = 1:T-1
	if t == 1
		# Yii[t] = B[t] * R̃[t] * B[t]' + Q̃[t+1]

		Yii[t] .= Q̃[t+1]
		mul!(tmp_nm, B[t], R̃[t])
		mul!(Yii[t], tmp_nm, transpose(B[t]), 1.0, 1.0)
	else
		# Yii[t] = A[t] * Q̃[t] * A[t]' + B[t] * R̃[t] * B[t]' + Q̃[t+1]

		mul!(tmp_nn, A[t], Q̃[t])
		mul!(Yii[t], tmp_nn, transpose(A[t]))
		mul!(tmp_nm, B[t], R̃[t])
		mul!(Yii[t], tmp_nm, transpose(B[t]), 1.0, 1.0)
		Yii[t] .+= Q̃[t+1]
	end

	t == T-1 && continue
	# Yij[t] = -Q̃[t+1] * A[t+1]'

	Yij[t] .= 0.0
	mul!(Yij[t], Q̃[t+1], transpose(A[t+1]), -1.0, 1.0)
end

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

for t = 1:T-1
	@show t
	if t == 1
		# Lii[t] = cholesky(Hermitian(Yii[t])).L

		Lii[t] .= copy(Yii[t])
		LAPACK.potrf!('L', Lii[t])
		# Lii[t] .= LowerTriangular(Lii[t])
	else
		# Lii[t] = cholesky(Hermitian(Yii[t] - Lji[t - 1] * Lji[t - 1]')).L
		Lii[t] .= copy(Yii[t])
		mul!(Lii[t], Lji[t-1], transpose(Lji[t-1]), -1.0, 1.0)
		LAPACK.potrf!('L', Lii[t])
		# Lii[t] .= LowerTriangular(Lii[t])
	end
	# Lji[t] = (LowerTriangular(Lii[t]) \ Yij[t])'

	Lji[t] .= copy(Yij[t])
	LAPACK.trtrs!('L', 'N', 'N', Lii[t], Lji[t])
	Lji[t] .= transpose(Lji[t])
end

for t = 1:T-1
	L[idx_n[t], idx_n[t]] = LowerTriangular(Lii[t])

	t == T-1 && continue

	L[idx_n[t+1], idx_n[t]] = Lji[t]
end

norm(cholesky(Hermitian(Y)).L - L)

# solve system
Δz = zeros(nz)
Δν = zeros(nd)

rd = r[1:nz]
rp = r[nz .+ (1:nd)]
β = -rp + C * H̃ * rd
Δν .= Y \ β
# Δν .= β
# LAPACK.potrs!('L', L, Δν)
Δz .= H̃ * (rd - C' * Δν)
norm(Δldl - [Δz; Δν])

# forward subsitution
y = [zeros(n) for t = 1:T-1]
for t = 1:T-1
	if t == 1
		# y[1] .= Lii[1] \ β[idx_n[1]]
		y[1] .= copy(β[idx_n[1]])
		LAPACK.trtrs!('L', 'N', 'N', Lii[t], y[1])
	else
		# y[t] .= Lii[t] \ (β[idx_n[t]] - Lji[t - 1] * y[t - 1])
		y[t] .= β[idx_n[t]]
		mul!(y[t], Lji[t - 1], y[t - 1], -1.0, 1.0)
		LAPACK.trtrs!('L', 'N', 'N', Lii[t], y[t])
	end
end

norm(vcat(y...) - L \ β, Inf)

# backward substitution
x = [zeros(n) for t = 1:T-1]
for t = T-1:-1:1
	@show t
	if t == T-1
		# x[t] = LowerTriangular(Lii[t])' \ y[t]
		x[t] .= y[t]
		LAPACK.trtrs!('L', 'T', 'N', Lii[t], x[t])
	else
		# x[t] = LowerTriangular(Lii[t])' \ (y[t] - Lji[t]' * x[t+1])
		x[t] .= y[t]
		mul!(x[t], transpose(Lji[t]), x[t+1], -1.0, 1.0)
		LAPACK.trtrs!('L', 'T', 'N', Lii[t], x[t])
	end
end

s = T-1
# s = T-1
norm(vcat(x...) - Y \ β, Inf)#[end-s*n+1:end-(s-1)], Inf)
norm(vcat(x...) - L' \ (L \ β), Inf)#[end-s*n+1:end-(s-1)], Inf)
