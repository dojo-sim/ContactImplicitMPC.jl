@testset "Cost Function" begin
    model = ContactControl.get_model("quadruped")
    nq = model.dim.q
    nu = model.dim.u
    nw = model.dim.w
    nc = model.dim.c
    nb = model.dim.b

    H = 20
    cost0 = ContactControl.CostFunction(H, model.dim)
    @test cost0.H == H
    @test length(cost0.Qq) == H
    @test length(cost0.Qu) == H
    @test length(cost0.Qγ) == H
    @test length(cost0.Qb) == H
end
