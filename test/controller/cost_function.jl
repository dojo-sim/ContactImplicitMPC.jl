@testset "Cost Function" begin
    model = ContactControl.get_model("quadruped")

    H = 10
    cost = ContactControl.cost_function(H, model.dim)
    @test length(cost.q) == H
    @test length(cost.u) == H
    @test length(cost.γ) == H
    @test length(cost.b) == H
end
