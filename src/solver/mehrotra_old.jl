#
#
# # interior point solver
# function interior_point_solve_old!(ip::Mehrotra{T}) where T
#
#     # space
#     s = ip.s
#
#     # methods
#     r! = ip.methods.r!
#     rm! = ip.methods.rm!
#     rz! = ip.methods.rz!
#     rθ! = ip.methods.rθ!
#
#     # options
#     opts = ip.opts
#     r_tol = opts.r_tol
#     κ_tol = opts.κ_tol
#     max_iter_inner = opts.max_iter_inner
#     max_time = opts.max_time
#     diff_sol = opts.diff_sol
#     res_norm = opts.res_norm
#     reg = opts.reg
#     ϵ_min = opts.ϵ_min
#     verbose = opts.verbose
#
#     # unpack pre-allocated data
#     z = ip.z
#     Δaff = ip.Δaff
#     Δ = ip.Δ
#     r = ip.r
#     rz = ip.rz
#     idx_ineq = ip.idx_ineq
#     idx_soc = ip.idx_soc
#     θ = ip.θ
#     κ = ip.κ
#     v_pr = ip.v_pr
#     v_du = ip.v_du
#     z_y1 = ip.z_y1
#     z_y2 = ip.z_y2
#     Δaff_y1 = ip.Δaff_y1
#     Δaff_y2 = ip.Δaff_y2
#     Δ_y1 = ip.Δ_y1
#     Δ_y2 = ip.Δ_y2
#     iy1 = ip.iy1
#     iy2 = ip.iy2
#     reg_pr = ip.reg_pr
#     reg_du = ip.reg_du
#     solver = ip.solver
#     ip.iterations = 0
#
#     # initialize regularization
#     reg_pr[1] = opts.reg_pr_init
#     reg_du[1] = opts.reg_du_init
#
#     # z .= 1e-3
#     println("z: ", scn.(z[[10,12,15]], digits=4))
#
#     # compute residual, residual Jacobian
#     ip.methods.rm!(r, z, 0.0 .* Δaff, θ, 0.0) # here we set κ = 0, Δ = 0
#     r_merit = norm(r, res_norm)
#
#     elapsed_time = 0.0
#
#     for j = 1:max_iter_inner
#         elapsed_time >= max_time && break
#         elapsed_time += @elapsed begin
#
#             # @show j
#             # @show norm(z)
#             # @show norm(θ)
#             # @show minimum(abs.(z))
#             # plt = plot()
#             # plot!(z)
#             # display(plt)
#
#             # check for converged residual
#             if r_merit < r_tol
#                 break
#             end
#             ip.iterations += 1
#             # compute residual Jacobian
#             rz!(rz, z, θ)
#             # @show norm(rz)
#
#             # regularize (fixed, TODO: adaptive)
#             reg && regularize!(v_pr, v_du, reg_pr[1], reg_du[1])
#
#             # compute affine search direction
#             linear_solve!(solver, Δaff, rz, r)
#             # @show norm(Δaff)
#             # plt = plot()
#             # plot!(Δaff)
#             # display(plt)
#             # αaff = step_length(z_y1, z_y2, Δaff_y1, Δaff_y2, τ=1.0)
#             αaff = step_length(z, Δaff, iy1, iy2; τ = 1.0)
#             println("αaff: ", scn(αaff, digits=6))
#
#             μaff = (z_y1 - αaff * Δaff[iy1])' * (z_y2 - αaff * Δaff[iy2]) / length(iy1)
#             println("μaff: ", scn(μaff, digits=6))
#             # @show nbil
#             # @show norm(αaff)
#             # @show norm(μaff)
#
#             μ = z_y1'*z_y2 / length(z_y1)
#             σ = (μaff / μ)^3
#             println("μ: ", scn(μ, digits=6))
#             println("σ: ", scn(σ, digits=6))
#
#
#             # Compute corrector residual
#             # rm!(rm, z, Δaff, θ, σ*μ) # here we set κ = σ*μ, Δ = Δaff
#             ip.methods.rm!(r, z, Δaff, θ, σ*μ)
#
#
#             # Compute corrector search direction
#             # linear_solve!(solver, Δ, rz, rm)
#             linear_solve!(solver, Δ, rz, r)
#             τ = progress(r_merit, ϵ_min=ϵ_min)
#             # α = step_length(z_y1, z_y2, Δ_y1, Δ_y2, τ=τ)
#             α = step_length(z, Δ, iy1, iy2; τ = τ)
#
#             # candidate point
#             candidate_point!(z, s, z, Δ, α)
#
#             # update
#             r!(r, z, θ, 0.0) # we set κ= 0.0 to measure the bilinear constraint violation
#             r_merit = norm(r, res_norm)
#             verbose && println("iter: ", j, "   res∞: ", scn(r_merit))
#         end
#     end
#
#     if r_merit > r_tol
#         @error "Mehrotra solver failed to reduce residual below r_tol."
#         return false
#     end
#
#     # differentiate solution
#     diff_sol && differentiate_solution!(ip)
#     return true
# end
#
# # TODO maybe we will need to implement this merit function to use κ_tol > b and r_tol > a
# # function merit(rlin::AbstractVector, rbil::AbstractVector, t::Real)
# # 	a = norm(rlin, t)
# # 	b = norm(rbil, t)
# # 	return a, b
# # end
#
# function progress(merit; ϵ_min=0.05)
#     ϵ = min(ϵ_min, merit^2)
#     τ = 1 - ϵ
#     return τ
# end
