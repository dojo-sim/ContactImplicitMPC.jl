"""
based on https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
"""

using BenchmarkTools
using InteractiveUtils

include(joinpath(module_dir(),
	"src/controller/newton_structure_solver/problem.jl"))
include(joinpath(module_dir(),
	"src/controller/newton_structure_solver/methods.jl"))

s = newton_structure_solver(nq, m, T)

for t = 1:T
	s.Q̃a[t] = Q̃a[t]
	s.Q̃b[t] = Q̃b[t]

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

compute_Y!(s.Yiia, s.Yiib, s.Yiic, s.Yiid, s.Yija, s.Yijb, s.Yijc, s.Yijd,
	s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a,
	s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.H)

info = @benchmark compute_Y!($s.Yiia, $s.Yiib, $s.Yiic, $s.Yiid,
	$s.Yija, $s.Yijb, $s.Yijc, $s.Yijd, $s.Aa, $s.Ab, $s.Ac, $s.Ba,
	$s.Q̃a, $s.Q̃b, $s.R̃a, $s.tmp_nqnq, $s.tmp_nqnq2, $s.tmp_nqm, $s.H)
@test info.memory == 0
@test info.allocs == 0

@code_warntype compute_Y!(s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a,
	s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.H)


update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Yiiav, s.Yiibv, s.Yiicv, s.Yiidv,
	s.Yijav, s.Yijbv, s.Yijcv, s.Yijdv, s.H)
info = @benchmark update_Y!($s.Yiis, $s.Yijs, $s.Yii, $s.Yij, $s.Yiia, $s.Yiib, $s.Yiic, $s.Yiid,
	$s.Yija, $s.Yijb, $s.Yijc, $s.Yijd,
	$s.Yiiav, $s.Yiibv, $s.Yiicv, $s.Yiidv, $s.Yijav, $s.Yijbv, $s.Yijcv, $s.Yijdv, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
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

compute_L!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)

info = @benchmark compute_L!($s.Liis, $s.Ljis, $s.Yiis, $s.Yijs, $s.tmp_nn, $s.tmp_nn2, $s.H)
@test info.memory == 0
@test info.allocs == 0

@code_warntype compute_L!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)

L = zeros(n * (T - 1), n * (T - 1))
for t = 1:T-1
	L[s.n_idx[t], s.n_idx[t]] = LowerTriangular(s.Liis[t])

	t == T-1 && continue

	L[s.n_idx[t+1], s.n_idx[t]] = transpose(s.Ljis[t])
end

@test norm(cholesky(Hermitian(Y)).L - L) < 1.0e-12

compute_β!(s.βn, s.βd, s.βe, s.rlagu, s.rlagqa, s.rlagqb,
	s.rdyn1, s.rdyn2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a, s.H)

info = @benchmark compute_β!($s.βn, $s.βd, $s.βe, $s.rlagu, $s.rlagqa, $s.rlagqb,
	$s.rdyn1, $s.rdyn2, $s.Aa, $s.Ab, $s.Ac, $s.Ba, $s.Q̃a, $s.Q̃b, $s.R̃a, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype compute_β!(s.βn, s.βd, s.βe, s.rlagu, s.rlagqa, s.rlagqb,
	s.rdyn1, s.rdyn2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a, s.H)

@test norm(β - vcat([[s.βd[t]; s.βe[t]] for t = 1:T-1]...)) < 1.0e-12
@test norm((β - vcat(s.βn...))) < 1.0e-12

compute_y!(s.y, s.Liis, s.Ljis, s.βn, s.H)

info = @benchmark compute_y!($s.y, $s.Liis, $s.Ljis, $s.βn, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype compute_y!(s.y, s.Liis, s.Ljis, s.βn, s.H)

@test norm(vcat(s.y...) - L \ β, Inf) < 1.0e-12

compute_Δν!(s.Δνn, s.Δνd, s.Δνe, s.Liis, s.Ljis, s.y, s.idx_nq, s.idx_nq2, s.H)

info = @benchmark compute_Δν!($s.Δνn, $s.Δνd, $s.Δνe, $s.Liis, $s.Ljis, $s.y, $s.idx_nq, $s.idx_nq2, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype compute_Δν!(s.Δνn, s.Δνd, s.Δνe, s.Liis, s.Ljis, s.y, s.idx_nq, s.idx_nq2, s.H)

@test norm(vcat(s.Δνn...) - Y \ β, Inf) < 1.0e-12
@test norm(vcat(s.Δνn...) - L' \ (L \ β), Inf) < 1.0e-12

compute_Δz!(s.Δzu, s.Δzqa, s.Δzqb, s.Δνd, s.Δνe, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a, s.rlagu, s.rlagqa, s.rlagqb, s.H)
info = @benchmark compute_Δz!($s.Δzu, $s.Δzqa, $s.Δzqb, $s.Δνd, $s.Δνe, $s.Aa, $s.Ab, $s.Ac, $s.Ba, $s.Q̃a, $s.Q̃b, $s.R̃a, $s.rlagu, $s.rlagqa, $s.rlagqb, $s.H)

@test info.memory == 0
@test info.allocs == 0

@code_warntype compute_Δz!(s.Δzu, s.Δzqa, s.Δzqb, s.Δνd, s.Δνe, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a, s.rlagu, s.rlagqa, s.rlagqb, s.H)

@test norm((Δz - vcat([[s.Δzu[t]; s.Δzqa[t]; s.Δzqb[t]] for t = 1:T-1]...))) < 1.0e-12

ContactControl.solve!(s)
info = @benchmark ContactControl.solve!(s)
@code_warntype ContactControl.solve!(s)

@test norm(Δ - vcat(vcat([[s.Δzu[t]; s.Δzqa[t]; s.Δzqb[t]] for t = 1:T-1]...), vcat([[s.Δνd[t]; s.Δνe[t]] for t = 1:T-1]...)))  < 1.0e-12
