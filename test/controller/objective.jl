@testset "Controller: Objective" begin
    s = get_simulation("quadruped", "flat_2D_lc", "flat")
    model = s.model
    env = s.env

    H = 10
    obj = ContactImplicitMPC.TrackingObjective(model, env, H)
    @test length(obj.q) == H
    @test length(obj.u) == H
    @test length(obj.γ) == H
    @test length(obj.b) == H
end
