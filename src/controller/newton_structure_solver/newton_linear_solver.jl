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

inv(Array(cholesky(H).U)) * inv(Array(cholesky(H).U))'







# run Flamingo MPC example

using UnicodePlots
UnicodePlots.show(p.newton.jac.R)

# @show "newton_solve!" #@@@
reset!(p.newton, p.ref_traj, warm_start = false,
	initial_offset = false, q0 = q0_sim, q1 = q1_sim)

# Compute implicit dynamics about traj
implicit_dynamics!(p.im_traj, s, p.newton.traj, κ = p.im_traj.ip[1].κ)

# Compute residual
residual!(p.newton.res, p.newton, p.newton.ν, p.im_traj, p.newton.traj, p.ref_traj)
update_jacobian!(p.newton.jac, p.im_traj, p.newton.obj, p.newton.traj.H, p.newton.β)

using BenchmarkTools
@benchmark linear_solve!($p.newton.solver, $p.newton.Δ.r, $p.newton.jac.R, $p.newton.res.r)
nq = model.dim.q # configuration
nu = model.dim.u # control
nc = model.dim.c # contact
nb = nc * friction_dim(env) # linear friction
nd = nq + nc + nb # implicit dynamics constraint
nr = nq + nu + nc + nb# + nd # size of a one-time-step block

idx_pr = collect(1:nr * H_mpc)
idx_du = collect(nr * H_mpc .+ (1:nd * H_mpc))
H = Array(@view Array(p.newton.jac.R)[idx_pr, idx_pr])
C = Array(@view Array(p.newton.jac.R)[idx_du, idx_pr])
invH = inv(Array(H))
rd = Array(@view p.newton.res.r[idx_du])
rp = Array(@view p.newton.res.r[idx_pr])
Δz = zeros(nr * H_mpc)
Δν = zeros(nd * H_mpc)
D = zeros(nd * H_mpc, nr * H_mpc)
Y = zeros(nd * H_mpc, nd * H_mpc)

function solve_schur!(Δz, Δν, H, C, invH, rd, rp, D, Y)
	mul!(D, C, invH)
	mul!(Y, D, transpose(C))
	LAPACK.potrf!('U', Y)
	mul!(Δν, D, rp)
	Δν .-= rd
	LAPACK.potrs!('U', Y, Δν)

	mul!(Δz, invH, rp)
	Δz .-= transpose(D) * Δν
end

@code_warntype solve_schur!(Δz, Δν, H, C, invH, rd, rp, D, Y)
@benchmark solve_schur!($Δz, $Δν, $H, $C, $invH, $rd, $rp, $D, $Y)

function solve_schur!(Δz, Δν, R, invH, r, D, Y)
	mul!(D, p.newton.jac.R[idx_du, idx_pr], invH)
	mul!(Y, D, transpose(p.newton.jac.R[idx_du, idx_pr]))
	LAPACK.potrf!('U', Y)
	mul!(Δν, D, p.newton.res.r[idx_pr])
	Δν .-= p.newton.res.r[idx_du]
	LAPACK.potrs!('U', Y, Δν)

	mul!(Δz, invH, p.newton.res.r[idx_pr])
	Δz .-= transpose(D) * Δν
end

@benchmark solve_schur!($Δz, $Δν, $p.newton.jac.R, $invH, $p.newton.res.r, $D, $Y)

norm(p.newton.Δ.r - [Δz; Δν], Inf)

# save trajectory
# @save joinpath(module_dir(), "src/dynamics/flamingo/simulations/flat.jld2") sim
# @load joinpath(module_dir(), "src/dynamics/flamingo/simulations/flat.jld2") sim

idx_pr = 1:nr * H_mpc
idx_du = nr * H_mpc .+ (1:nd * H_mpc)
H = view(p.newton.jac.R, idx_pr, idx_pr)
C = view(p.newton.jac.R, idx_du, idx_pr)
invH = inv(Array(H))
D = zeros(nd * H_mpc, nr * H_mpc)
rd = view(p.newton.res.r, idx_du)
rp = view(p.newton.res.r, idx_pr)
Δz = zeros(nr * H_mpc)
Δν = zeros(nd * H_mpc)
mul!(D, C, invH)
mul!(Y, D, transpose(C))
Δν .= -rd + D * rp
Δν .= Y \ (-rd + D * rp)

LAPACK.potrf!('U', Y)
LAPACK.potrs!('U', Y, Δν)
Δz .= invH * rp - transpose(C) * Δν

norm(p.newton.Δ.r - [Δz; Δν], Inf)
