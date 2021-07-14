"""
based on https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

using BenchmarkTools
using InteractiveUtils

include(joinpath(module_dir(),
	"src/controller/newton_structure_solver/problem.jl"))
include(joinpath(module_dir(),
	"src/controller/newton_structure_solver/methods.jl"))

s = ContactControl.newton_structure_solver(nq, m, T,
	opts = NewtonOptions(β_init = 0.0))

for t = 1:T
	s.Q̃a[t] = Q̃a[t]
	s.Q̃b[t] = Q̃b[t]
	s.Q̃v[t] = Q̃v[t]

	t == T && continue
	s.Aa[t] = Aa[t]
	s.Ab[t] = Ab[t]
	s.Ac[t] = Ac[t]
	s.Ba[t] = Ba[t]

	s.R̃a[t] = R̃[t]

	rlag = view(r, 1:nz)

	s.rlagu[t] = view(rlag, s.u_idx[t])
	s.rlagqa[t] = view(rlag, s.qa_idx[t])
	s.rlagqb[t] = view(rlag, s.qb_idx[t])

	rdyn = view(r, nz .+ (1:nd))
	s.rdyn1[t] = view(rdyn, s.n_idx[t][1:nq])
	s.rdyn2[t] = view(rdyn, s.n_idx[t][nq .+ (1:nq)])
end

ContactControl.compute_Y!(s.Yiia, s.Yiib, s.Yiic, s.Yiid, s.Yija, s.Yijb, s.Yijc, s.Yijd,
	s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a,
	s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.Inq, s.H)

info = @benchmark ContactControl.compute_Y!($s.Yiia, $s.Yiib, $s.Yiic, $s.Yiid,
	$s.Yija, $s.Yijb, $s.Yijc, $s.Yijd, $s.Aa, $s.Ab, $s.Ac, $s.Ba,
	$s.Q̃a, $s.Q̃b, $s.Q̃v, $s.R̃a, $s.tmp_nqnq, $s.tmp_nqnq2, $s.tmp_nqm, $s.Inq, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype ContactControl.compute_Y!(s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a,
	s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.Inq, s.H)

ContactControl.update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Yiiav, s.Yiibv, s.Yiicv, s.Yiidv,
	s.Yijav, s.Yijbv, s.Yijcv, s.Yijdv, s.H)
info = @benchmark ContactControl.update_Y!($s.Yiis, $s.Yijs, $s.Yii, $s.Yij, $s.Yiia, $s.Yiib, $s.Yiic, $s.Yiid,
	$s.Yija, $s.Yijb, $s.Yijc, $s.Yijd,
	$s.Yiiav, $s.Yiibv, $s.Yiicv, $s.Yiidv, $s.Yijav, $s.Yijbv, $s.Yijcv, $s.Yijdv, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype ContactControl.update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Yiiav, s.Yiibv, s.Yiicv, s.Yiidv,
	s.Yijav, s.Yijbv, s.Yijcv, s.Yijdv, s.H)

Y = zeros(n * (T - 1), n * (T - 1))
for t = 1:T-1
	Y[s.n_idx[t], s.n_idx[t]] = s.Yiis[t]

	t == T-1 && continue

	Y[s.n_idx[t], s.n_idx[t+1]] = s.Yijs[t]
	Y[s.n_idx[t+1], s.n_idx[t]] = s.Yijs[t]'
end

@test norm(Y - C * S̃ * C') < 1.0e-12

ContactControl.compute_L!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)

info = @benchmark ContactControl.compute_L!($s.Liis, $s.Ljis, $s.Yiis, $s.Yijs, $s.tmp_nn, $s.tmp_nn2, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype ContactControl.compute_L!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)

L = zeros(n * (T - 1), n * (T - 1))
for t = 1:T-1
	L[s.n_idx[t], s.n_idx[t]] = LowerTriangular(s.Liis[t])

	t == T-1 && continue

	L[s.n_idx[t+1], s.n_idx[t]] = transpose(s.Ljis[t])
end

@test norm(cholesky(Hermitian(Y)).L - L) < 1.0e-12

ContactControl.compute_β!(s.βn, s.β1, s.β2, s.rlagu, s.rlagqa, s.rlagqb,
	s.rdyn1, s.rdyn2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a, s.H)

info = @benchmark ContactControl.compute_β!($s.βn, $s.β1, $s.β2, $s.rlagu, $s.rlagqa, $s.rlagqb,
	$s.rdyn1, $s.rdyn2, $s.Aa, $s.Ab, $s.Ac, $s.Ba, $s.Q̃a, $s.Q̃b, $s.Q̃v, $s.R̃a, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype ContactControl.compute_β!(s.βn, s.β1, s.β2, s.rlagu, s.rlagqa, s.rlagqb,
	s.rdyn1, s.rdyn2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a, s.H)

@test norm((β - vcat([[s.β1[t]; s.β2[t]] for t = 1:T-1]...))) < 1.0e-12
@test norm((β - vcat(s.βn...))[1:s.nq]) < 1.0e-12

ContactControl.compute_y!(s.y, s.Liis, s.Ljis, s.βn, s.H)

info = @benchmark ContactControl.compute_y!($s.y, $s.Liis, $s.Ljis, $s.βn, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype ContactControl.compute_y!(s.y, s.Liis, s.Ljis, s.βn, s.H)

@test norm(vcat(s.y...) - L \ β, Inf) < 1.0e-12

ContactControl.compute_Δν!(s.Δνn, s.Δν1, s.Δν2, s.Liis, s.Ljis, s.y, s.idx_nq, s.idx_nq2, s.H)

info = @benchmark ContactControl.compute_Δν!($s.Δνn, $s.Δν1, $s.Δν2, $s.Liis, $s.Ljis,
	$s.y, $s.idx_nq, $s.idx_nq2, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype ContactControl.compute_Δν!(s.Δνn, s.Δν1, s.Δν2, s.Liis, s.Ljis, s.y, s.idx_nq, s.idx_nq2, s.H)

@test norm(vcat(s.Δνn...) - Δν, Inf) < 1.0e-12
@test norm(vcat(s.Δνn...) - Y \ β, Inf) < 1.0e-12
@test norm(vcat(s.Δνn...) - L' \ (L \ β), Inf) < 1.0e-12
@test norm(vcat([[s.Δν1[t]; s.Δν2[t]] for t = 1:T-1]...) - Y \ β, Inf) < 1.0e-12

ContactControl.compute_Δz!(s.Δu, s.Δqa, s.Δqb, s.Δν1, s.Δν2, s.Aa, s.Ab, s.Ac, s.Ba,
	s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a, s.rlagu, s.rlagqa, s.rlagqb, s.H)

info = @benchmark ContactControl.compute_Δz!($s.Δu, $s.Δqa, $s.Δqb, $s.Δν1, $s.Δν2,
	$s.Aa, $s.Ab, $s.Ac, $s.Ba, $s.Q̃a, $s.Q̃b, $s.Q̃v, $s.R̃a,
	$s.rlagu, $s.rlagqa, $s.rlagqb, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype ContactControl.compute_Δz!(s.Δu, s.Δqa, s.Δqb, s.Δν1, s.Δν2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.Q̃v, s.R̃a, s.rlagu, s.rlagqa, s.rlagqb, s.H)

@test norm((Δz - vcat([[s.Δu[t]; s.Δqa[t]; s.Δqb[t]] for t = 1:T-1]...))) < 1.0e-12

ContactControl.factorize!(s)
info = @benchmark ContactControl.factorize!($s)
@code_warntype ContactControl.factorize!(s)

ContactControl.solve!(s)
info = @benchmark ContactControl.solve!($s)
@code_warntype ContactControl.solve!(s)

@test norm(Δ - vcat(vcat([[s.Δu[t]; s.Δqa[t]; s.Δqb[t]] for t = 1:T-1]...), vcat([[s.Δν1[t]; s.Δν2[t]] for t = 1:T-1]...)))  < 1.0e-12

# # test residual
# sim = ContactControl.get_simulation("flamingo", "flat_2D_lc", "flat")
# model = sim.model
# env = sim.env
#
# ref_traj = deepcopy(ContactControl.get_trajectory(model, flat_2D_lc,
#     joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
#     load_type = :split_traj_alt))
#
# H_mpc = 10
# s = ContactControl.newton_structure_solver(model.dim.q, model.dim.u, H_mpc, ρ = 1.0e-5)
# obj_mpc = ContactControl.quadratic_objective(model, H_mpc)
# ContactControl.update_objective!(s, obj_mpc)
#
# # random fill
# for t = 1:s.H
# 	s.qa[t] = ref_traj.q[t]
# 	s.qb[t] = ref_traj.q[t+1]
# 	t == s.H && continue
# 	s.u[t] = ref_traj.u[t]
# 	s.ν1[t] = zeros(s.nq)
# 	s.ν2[t] = zeros(s.nq)
# end
#
# for t = 1:s.H+1
# 	s.q_ref[t] = ref_traj.q[t]
# 	t >= s.H && continue
# 	s.u_ref[t] = ref_traj.u[t]
# end
#
# z0 = zeros(s.nz + s.nd)
# for t = 1:s.H-1
# 	z0[s.u_idx[t]] = s.u[t]
# 	z0[s.qa_idx[t]] = s.qa[t+1]
# 	z0[s.qb_idx[t]] = s.qb[t+1]
# 	z0[s.nz .+ (t-1) * 2 * s.nq .+ (1:s.nq)] = s.ν1[t]
# 	z0[s.nz .+ (t-1) * 2 * s.nq + s.nq .+ (1:s.nq)] = s.ν2[t]
# end
#
# ip_opts = eval(ContactControl.interior_point_options(:interior_point))(
# 			κ_init = 1.0e-4,
# 			κ_tol = 2.0 * 1.0e-4,
# 			r_tol = 1.0e-8,
# 			diff_sol = true,
# 			solver = :empty_solver)
#
# update_dynamics_jacobian!(s, )
# info = @benchmark update_dynamics_jacobian!($s, $)
# @code_warntype update_dynamics_jacobian!(s, )
#
# @test info.memory == 0
# @test info.allocs == 0
#
# dynamics_constraints!(s, s.qa, s.qb, )
#
# info = @benchmark dynamics_constraints!($s, $s.qa, $s.qb, $lci_traj.dq2)
# @code_warntype dynamics_constraints!(s, s.qa, s.qb, lci_traj.dq2)
# @test info.memory == 0
# @test info.allocs == 0
#
#
# lagrangian_gradient!(s, s.u, s.qa, s.qb, s.ν1, s.ν2)
# info = @benchmark lagrangian_gradient!($s, $s.u, $s.qa, $s.qb, $s.ν1, $s.ν2)
# @code_warntype lagrangian_gradient!(s, s.u, s.qa, s.qb, s.ν1, s.ν2)
#
# @test info.memory == 0
# @test info.allocs == 0
#
# function lagrangian_dynamics(z)
# 	L = 0.0
#
# 	for t = 1:s.H-1
# 		u1 = z[s.u_idx[t]]
# 		qa1 = z[s.qa_idx[t]]
# 		qb1 = z[s.qb_idx[t]]
#
# 		ν1 = z[s.nz .+ (t-1) * 2 * s.nq .+ (1:s.nq)]
# 		ν2 = z[s.nz .+ (t-1) * 2 * s.nq + s.nq .+ (1:s.nq)]
#
# 		L += 0.5 * transpose(u1 - s.u_ref[t]) * s.Ra[t] * (u1 - s.u_ref[t])
# 		L += 0.5 * transpose(qa1 - s.q_ref[t+1]) * s.Qa[t+1] * (qa1 - s.q_ref[t+1])
# 		L += 0.5 * transpose(qb1 - s.q_ref[t+2]) * s.Qb[t+1] * (qb1 - s.q_ref[t+2])
#
# 		L += transpose(qa1) * s.Qv[t+1] * qb1
#
# 		if t == 1
# 			L += transpose(ν1) * (qa1 - s.qb[1])
# 			L += transpose(ν2) * (qb1 - s.Aa[1] * s.qa[1] - s.Ab[1] * s.qb[1] - s.Ba[1] * u1)
# 		else
# 			qa0 = z[s.qa_idx[t-1]]
# 			qb0 = z[s.qb_idx[t-1]]
# 			L += transpose(ν1) * (qa1 - qb0)
# 			L += transpose(ν2) * (qb1 - s.Aa[t] * qa0 - s.Ab[t] * qb0 - s.Ba[t] * u1)
# 		end
# 	end
#
# 	return L
# end
#
# lag_grad = vcat(vcat([[s.rlagu[t]; s.rlagqa[t]; s.rlagqb[t]] for t = 1:s.H-1])..., vcat([[s.rdyn1[t]; s.rdyn2[t]] for t = 1:s.H-1])...)
# lag_grad_fd = ForwardDiff.gradient(lagrangian_dynamics, z0)
#
# @test norm((lag_grad - lag_grad_fd)[1:s.nz]) < 1.0e-12
# @test norm((lag_grad - lag_grad_fd)[s.nz .+ (1:s.nd)]) < 1.0e-12
#
# compute_residual!(s, s.u, s.qa, s.qb, s.ν1, s.ν2, lci_traj)
# info = @benchmark compute_residual!($s, $s.u, $s.qa, $s.qb, $s.ν1, $s.ν2,  $lci_traj)
# @code_warntype compute_residual!(s, s.u, s.qa, s.qb, s.ν1, s.ν2, lci_traj)
#
# factorize!(s)
# ContactControl.solve!(s)
#
# α = 1.0
# step!(s, α)
# @code_warntype step!(s, α)
# info = @benchmark step!($s, $α)
