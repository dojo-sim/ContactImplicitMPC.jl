@testset "Objective" begin
    model = ContactControl.get_model("quadruped")

    H = 10
    obj = ContactControl.TrackingObjective(H, model.dim)
    @test length(obj.q) == H
    @test length(obj.u) == H
    @test length(obj.γ) == H
    @test length(obj.b) == H
end
