# A = Diagonal(5.0 * ones(2))
# B = Diagonal(ones(2))
# z = zeros(2,2)
#
# H = [A B z z z;
#      B' A B z z
# 	 z B' A B z;
# 	 z z B' A B;
# 	 z z z B' A]
#
# H
# cholesky(H)
# inv(Array(H))
#
# H = [A B z z z;
#      B' A z z z
# 	 z z A B z;
# 	 z z B' A z;
# 	 z z z z A]
#
# inv(Array(H))
#
#
#
# R = Diagonal(2.0 * ones(2))
# W = Diagonal(5.0 * ones(3))
# N = zeros(3, 2)
# M = zeros(2, 2)
# P = zeros(3, 3)
#
# H = [R N' M N' M N';
#      N W N W N P;
# 	 M N' R N' M N';
# 	 N W N W N W;
# 	 M N' M N' R N';
# 	 N P N W N W]
#
#
# H + inv(Array(H))

n = 2
m = 1
T = 10
nz = n * T + m * T

Q = [Diagonal(ones(n)) for t = 1:T]
R = [Diagonal(ones(m)) for t = 1:T]
V = [Diagonal(ones(n)) for t = 1:T]

u_idx = [(t - 1) * (n + m) .+ (1:m) for t = 1:T]
x_idx = [(t - 1) * (n + m) + m .+ (1:n) for t = 1:T]

H = zeros(nz, nz)

for t = 1:T
	# control cost
	H[u_idx[t], u_idx[t]] = R[t]

	# configuration cost + velocity components
	H[x_idx[t], x_idx[t]] = Q[t] + V[t] + (t == T ? zero(V[t]) : V[t+1])

	# velocity components
	t == T && continue
	H[x_idx[t], x_idx[t+1]] = -V[t]
	H[x_idx[t+1], x_idx[t]] = -V[t]
end

using UnicodePlots
UnicodePlots.spy(H)
UnicodePlots.spy(inv(H))

p = n * T

C = zeros(p, nz)
dfdq0 = ones(n, n)
dfdq1 = ones(n, n)
dfdu1 = ones(n, m)

for t = 1:T
	t_idx = (t - 1) * n .+ (1:n)

	C[t_idx, u_idx[t]] = -dfdu1
	C[t_idx, x_idx[t]] = Diagonal(ones(n))

	t == 1 && continue
	C[t_idx, x_idx[t-1]] = -dfdq1

	t == 2 && continue
	C[t_idx, x_idx[t-2]] = -dfdq0
end

UnicodePlots.spy(H)












UnicodePlots.spy(C)










UnicodePlots.spy(C * inv(H) * C')
