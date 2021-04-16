@testset "Indices" begin
    model = get_model("hopper_2D")
    terms, vars = ContactControl.get_bilinear_indices(model)
    idyn, irst, ibil, ialt = ContactControl.linearization_term_index(model)
    ix, iy1, iy2 = ContactControl.linearization_var_index(model)

    @test iy1 == vcat(Vector.([vars[1][1], vars[2][1], vars[3][1]])...)
    @test iy2 == vcat(Vector.([vars[1][2], vars[2][2], vars[3][2]])...)
    @test ix == setdiff(1:nz, vcat(iy1, iy2))

    @test ibil == vcat(Vector.([terms[1], terms[2], terms[3]])...)
end
