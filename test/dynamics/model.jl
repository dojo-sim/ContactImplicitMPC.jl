# @testset "Indices" begin
#     s = ContactImplicitMPC.get_simulation("hopper_2D", "flat_2D_lc", "flat")
#
#     model = s.model
#     env = s.env
#     for quat ∈ (false, true)
#         terms, vars = ContactImplicitMPC.get_bilinear_indices(model, env, quat = quat)
#         idyn, irst, ibil, ialt = ContactImplicitMPC.linearization_term_index(model, env, quat = quat)
#         ix, iy1, iy2 = ContactImplicitMPC.linearization_var_index(model, env, quat = quat)
#         nquat = 0
#         @test iy1 == vcat(Vector.([vars[1][1], vars[2][1], vars[3][1]])...)
#         @test iy2 == vcat(Vector.([vars[1][2], vars[2][2], vars[3][2]])...)
#         @test ix == setdiff(1:ContactImplicitMPC.num_var(model, env) - nquat, vcat(iy1, iy2))
#         @test ibil == vcat(Vector.([terms[1], terms[2], terms[3]])...)
#     end
# end
