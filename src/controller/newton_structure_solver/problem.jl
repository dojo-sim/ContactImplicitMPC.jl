# dimensions
nq = 10
n = 2 * nq
m = 10
T = 20
nz = m * (T - 1) + n * (T - 1)
nd = n * (T - 1)

# problem data
Aa = [SMatrix{nq, nq}(Array(-1.0 * Diagonal(0.1 * rand(nq) + ones(nq)))) for t = 1:T-1]
Ab = [SMatrix{nq, nq}(Array(-1.0 * Diagonal(0.1 * rand(nq) + ones(nq)))) for t = 1:T-1]
Ac = [SMatrix{nq, nq}(Array(Diagonal(ones(nq)))) for t = 1:T-1]
Ba = [SMatrix{nq,m}(-1.0 * rand(nq, m)) for t = 1:T-1]
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
@benchmark linear_solve!($solver_lu, $Δlu, $J, $r)

# # QDLDL
# solver_ldl = ldl_solver(J)
# Δldl = zeros(nz + nd)
# @benchmark linear_solve!($solver_ldl, $Δldl, $J, $r)
# Δldl

norm(Δ - Δlu)
norm(Δ - [Δz; Δν])
