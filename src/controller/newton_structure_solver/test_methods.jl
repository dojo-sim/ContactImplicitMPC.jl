"""
Fast MPC (modified for implicit dynamics (f(xt, ut, xt+1) = 0))
https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf
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

Y = zeros(n * (T - 1), n * (T - 1))
Y
computeYs!(s.Yiia, s.Yiib, s.Yiic, s.Yiid, s.Yija, s.Yijb, s.Yijc, s.Yijd,
	s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a,
	s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.H)
@benchmark computeYs!($s.Yiia, $s.Yiib, $s.Yiic, $s.Yiid,
	$s.Yija, $s.Yijb, $s.Yijc, $s.Yijd, $s.Aa, $s.Ab, $s.Ac, $s.Ba,
	$s.Q̃a, $s.Q̃b, $s.R̃a, $s.tmp_nqnq, $s.tmp_nqnq2, $s.tmp_nqm, $s.H)
@code_warntype computeYs!(s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a,
	s.tmp_nqnq, s.tmp_nqnq2, s.tmp_nqm, s.H)

update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Yiiav, s.Yiibv, s.Yiicv, s.Yiidv,
	s.Yijav, s.Yijbv, s.Yijcv, s.Yijdv, s.H)
@benchmark update_Y!($s.Yiis, $s.Yijs, $s.Yii, $s.Yij, $s.Yiia, $s.Yiib, $s.Yiic, $s.Yiid,
	$s.Yija, $s.Yijb, $s.Yijc, $s.Yijd,
	$s.Yiiav, $s.Yiibv, $s.Yiicv, $s.Yiidv, $s.Yijav, $s.Yijbv, $s.Yijcv, $s.Yijdv, $s.H)
@code_warntype update_Y!(s.Yiis, s.Yijs, s.Yii, s.Yij, s.Yiia, s.Yiib, s.Yiic, s.Yiid,
	s.Yija, s.Yijb, s.Yijc, s.Yijd, s.Yiiav, s.Yiibv, s.Yiicv, s.Yiidv,
	s.Yijav, s.Yijbv, s.Yijcv, s.Yijdv, s.H)

for t = 1:T-1
	Y[s.n_idx[t], s.n_idx[t]] = s.Yiis[t]

	t == T-1 && continue

	Y[s.n_idx[t], s.n_idx[t+1]] = s.Yijs[t]
	Y[s.n_idx[t+1], s.n_idx[t]] = s.Yijs[t]'
end

norm(Y - C * S̃ * C')

###
L = zeros(n * (T - 1), n * (T - 1))

computeLs!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)
@benchmark computeLs!($s.Liis, $s.Ljis, $s.Yiis, $s.Yijs, $s.tmp_nn, $s.tmp_nn2, $s.H)
@code_warntype computeLs!(s.Liis, s.Ljis, s.Yiis, s.Yijs, s.tmp_nn, s.tmp_nn2, s.H)


for t = 1:T-1
	L[s.n_idx[t], s.n_idx[t]] = LowerTriangular(s.Liis[t])

	t == T-1 && continue

	L[s.n_idx[t+1], s.n_idx[t]] = transpose(s.Ljis[t])
end

norm(cholesky(Hermitian(Y)).L - L)

update_β!(s.βn, s.βd, s.βe, s.rlagu, s.rlagqa, s.rlagqb,
	s.rdyn1, s.rdyn2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a, s.H)
@benchmark update_β!($s.βn, $s.βd, $s.βe, $s.rlagu, $s.rlagqa, $s.rlagqb,
	$s.rdyn1, $s.rdyn2, $s.Aa, $s.Ab, $s.Ac, $s.Ba, $s.Q̃a, $s.Q̃b, $s.R̃a, $s.H)
@code_warntype update_β!(s.βn, s.βd, s.βe, s.rlagu, s.rlagqa, s.rlagqb,
	s.rdyn1, s.rdyn2, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a, s.H)
norm(β - vcat([[s.βd[t]; s.βe[t]] for t = 1:T-1]...))
norm((β - vcat(s.βn...)))


forward_substitution_s!(s.y, s.Liis, s.Ljis, s.βn, s.H)
@benchmark forward_substitution_s!($s.y, $s.Liis, $s.Ljis, $s.βn, $s.H)
@code_warntype forward_substitution_s!(s.y, s.Liis, s.Ljis, s.βn, s.H)
norm(vcat(s.y...) - L \ β, Inf)

backward_substitution_s!(s.Δνn, s.Liis, s.Ljis, s.y, s.H)
@benchmark backward_substitution_s!($s.Δνn, $s.Liis, $s.Ljis, $s.y, $s.H)
@code_warntype backward_substitution_s!(s.Δνn, s.Liis, s.Ljis, s.y, s.H)

norm(vcat(s.Δνn...) - Y \ β, Inf)
norm(vcat(s.Δνn...) - L' \ (L \ β), Inf)

for t = 1:T-1
	s.Δνd[t] = s.Δνn[t][1:s.nq]
	s.Δνe[t] = s.Δνn[t][s.nq .+ (1:nq)]
end

update_Δz!(s.Δzu, s.Δzqa, s.Δzqb, s.Δνd, s.Δνe, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a, s.rlagu, s.rlagqa, s.rlagqb, s.H)
@benchmark update_Δz!($s.Δzu, $s.Δzqa, $s.Δzqb, $s.Δνd, $s.Δνe, $s.Aa, $s.Ab, $s.Ac, $s.Ba, $s.Q̃a, $s.Q̃b, $s.R̃a, $s.rlagu, $s.rlagqa, $s.rlagqb, $s.H)
@code_warntype update_Δz!(s.Δzu, s.Δzqa, s.Δzqb, s.Δνd, s.Δνe, s.Aa, s.Ab, s.Ac, s.Ba, s.Q̃a, s.Q̃b, s.R̃a, s.rlagu, s.rlagqa, s.rlagqb, s.H)

norm((Δz - vcat([[s.Δzu[t]; s.Δzqa[t]; s.Δzqb[t]] for t = 1:T-1]...)))
