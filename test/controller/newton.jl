@testset "Newton" begin

    # Test set_traj!
    T = Float64
    H = 10
    h = 0.1
    model = ContactControl.get_model("quadruped")
    target = ContactControl.ContactTraj(H, h, model.dim)
    source = ContactControl.ContactTraj(H, h, model.dim)
    nd = model.dim.q + model.dim.c + model.dim.b
    νtarget = [-30*ones(SizedVector{nd,T}) for t=1:H]
    νsource = [+30*ones(SizedVector{nd,T}) for t=1:H]
    for t = 1:H
        source.q[t] .= +1.0
        source.u[t] .= +2.0
        source.w[t] .= +3.0
        source.γ[t] .= +4.0
        source.b[t] .= +5.0
        source.z[t] .= +6.0
        source.θ[t] .= +7.0

        target.q[t] .= -1.0
        target.u[t] .= -2.0
        target.w[t] .= -3.0
        target.γ[t] .= -4.0
        target.b[t] .= -5.0
        target.z[t] .= -6.0
        target.θ[t] .= -7.0
    end
    Δ0 = ContactControl.Residual11(H, model.dim)
    Δ0.r .+= 100*ones(ind0.nr*H)
    ContactControl.set_traj!(target, source, νtarget, νsource, Δ0, ind0)
    source.q[1] .+= 1000.0
    source.u[1] .+= 1000.0
    @test target.q[1][1] == 101.0
    @test target.u[1][1] == 102.0
    @test target.γ[1][1] == 104.0
    @test target.b[1][1] == 105.0
    @test νtarget[1][1]  == 130.0

end
