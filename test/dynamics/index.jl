@testset "Indices Linear Friction Cone" begin
    s = ContactControl.get_simulation("hopper_2D", "flat_2D_lc", "flat")
    model = s.model
    env = s.env

    # Test basic indices

    nq = model.dim.q
    nc = model.dim.c
    nquat = 1
    off = 0

    iq2 = index_q2(model, env, nquat = nquat)
    @test iq2 == Vector(off .+ (1:nq - nquat))
    off += nq - nquat
    index_γ1(model, env, nquat = nquat) == Vector(off .+ (1:nc))
    off += nc
    index_b1(model, env, nquat = nquat) == Vector(off .+ (1:nb))
    off += nb
    index_ψ1(model, env, nquat = nquat) == Vector(off .+ (1:nc))
    off += nc
    index_s1(model, env, nquat = nquat) == Vector(off .+ (1:nc))
    off += nc
    index_η1(model, env, nquat = nquat) == Vector(off .+ (1:nb))
    off += nb
    index_s2(model, env, nquat = nquat) == Vector(off .+ (1:nc))
    off += nc

    # Test aggregated indices
    for nquat = 0:2
        terms, vars = get_bilinear_indices(model, env, nquat = nquat)
        idyn, irst, ibil, ialt = linearization_term_index(model, env, nquat = nquat)
        ix, iy1, iy2 = linearization_var_index(model, env, nquat = nquat)

        @test iy1 == vcat(Vector.([vars[1][1], vars[2][1], vars[3][1]])...)
        @test iy2 == vcat(Vector.([vars[1][2], vars[2][2], vars[3][2]])...)
        @test ix == setdiff(1:ContactControl.num_var(model, env) - nquat, vcat(iy1, iy2))
        @test ibil == vcat(Vector.([terms[1], terms[2], terms[3]])...)
    end

    nquat = 1
    inequality_indices(model, env, nquat = nquat) == Vector(nq:num_var(model, env) .- nquat)
    soc_indices(model, env, nquat = nquat) == []
end


@testset "Indices Nonlinear friction cone" begin
    s = ContactControl.get_simulation("hopper_2D", "flat_2D_nc", "flat")
    model = s.model
    env = s.env

    # Test basic indices

    nq = model.dim.q
    nc = model.dim.c
    nb = nc * friction_dim(env)
    nquat = 1
    off = 0

    iq2 = index_q2(model, env, nquat = nquat)
    @test iq2 == Vector(off .+ (1:nq - nquat))
    off += nq - nquat
    index_γ1(model, env, nquat = nquat) == Vector(off .+ (1:nc))
    off += nc
    index_b1(model, env, nquat = nquat) == Vector(off .+ (1:nb))
    off += nb
    index_ψ1(model, env, nquat = nquat) == Vector(off .+ (1:nc))
    off += nc
    index_s1(model, env, nquat = nquat) == Vector(off .+ (1:nc))
    off += nc
    index_η1(model, env, nquat = nquat) == Vector(off .+ (1:nb))
    off += nb
    index_s2(model, env, nquat = nquat) == Vector(off .+ (1:nc))
    off += nc

    # Test aggregated indices
    for nquat = 0:2
        terms, vars = get_bilinear_indices(model, env, nquat = nquat)
        idyn, irst, ibil, ialt = linearization_term_index(model, env, nquat = nquat)
        ix, iy1, iy2 = linearization_var_index(model, env, nquat = nquat)

        @test iy1 == vcat(Vector.([vars[1][1], vars[2][1], vars[3][1]])...)
        @test iy2 == vcat(Vector.([vars[1][2], vars[2][2], vars[3][2]])...)
        @test ix == setdiff(1:ContactControl.num_var(model, env) - nquat, vcat(iy1, iy2))
        @test ibil == vcat(Vector.([terms[1], terms[2], terms[3]])...)
    end

    nquat = 1
    inequality_indices(model, env, nquat = nquat) == [(nq - nquat) .+ (1:nc); (nq - nquat + nc + nb + nc) .+ (1:nc)]
    soc_indices(model, env, nquat = nquat) == [(nq - nquat + nc) .+ [nb .+ (1:nc); 1:nb],
        (nq - nquat + nc + nb + nc + nc) .+ [nb .+ (1:nc); 1:nb]]
end
