@testset "Indices" begin
    s = ContactControl.get_simulation("hopper_2D", "flat_2D_lc", "flat")

    model = s.model
    env = s.env
    for nquat = 0:2
        terms, vars = ContactControl.get_bilinear_indices(model, env, nquat = nquat)
        idyn, irst, ibil, ialt = ContactControl.linearization_term_index(model, env, nquat = nquat)
        ix, iy1, iy2 = ContactControl.linearization_var_index(model, env, nquat = nquat)

        @test iy1 == vcat(Vector.([vars[1][1], vars[2][1], vars[3][1]])...)
        @test iy2 == vcat(Vector.([vars[1][2], vars[2][2], vars[3][2]])...)
        @test ix == setdiff(1:ContactControl.num_var(model, env) - nquat, vcat(iy1, iy2))
        @test ibil == vcat(Vector.([terms[1], terms[2], terms[3]])...)
    end
end
