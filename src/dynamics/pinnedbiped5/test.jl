
@testset "h parameterization" begin
    model = get_model("biped5", surf="flat")
    pinnedmodel = get_model("pinnedbiped5", surf="flat")

    q0 = [0.0, 0.0, 0.0, 0.5, 0.5, -0.5, -0.5]
    # set_robot!(vis, model, q0, r=0.004)
    p_com = com_func(model, q0)
    @test p_com[1] == 0.0

    q1 = [0.0, 0.0, 0.0, π/2, π/2, -π/2, -π/2]
    # set_robot!(vis, model, q1, r=0.004)
    p_com = com_func(model, q1)
    @test p_com[1] == 0.0
    @test p_com[2] >= 0.0
    h1 = h_func(model, q1)
    r_com1, r_foot1, θ_foot1, θ_trunk1 = h1
    θ_foot1 < 0.0
    θ_foot1 > - pi
    θ_trunk1

    q2 = [0.0, 0.0, 0.0, 0.5, 0.5, -0.5, -0.5]
    # set_robot!(vis, model, q2, r=0.004)
    h2 = h_func(model, q2)
    r_com2, r_foot2, θ_foot2, θ_trunk2 = h2
    θcom2 = θcom_func(model, q2)
    @test r_com2 > 0.0
    @test r_foot2 > 0.0
    @test norm(r_com2 - r_foot2) < 1e-10
    @test norm(θ_foot2 - 2θ_trunk2) < 1e-10
    @test θ_foot2 < -1.0
    @test θcom2 > π/2

    nqp = pinnedmodel.dim.q
    q  = [0.0; 0.0; rand(nqp)]
    qp = rand(nqp)
    @test q  == unpin_state(pin_state(q))
    @test qp == pin_state(unpin_state(qp))

    q_inip = [-0.228, 0.228, -0.1, -0.1, -0.3]
    q_midp = [-0.35,  0.2,   -0.1,  0.3, -0.4]
    q_endp = [-0.3, -0.1, -0.1, 0.228, -0.228]
    @test s_func(pinnedmodel, q_inip, q_ini=q_inip, q_mid=q_midp, q_end=q_endp) == 0.0
    @test s_func(pinnedmodel, q_endp, q_ini=q_inip, q_mid=q_midp, q_end=q_endp) == 1.0
    @test 0.0 < s_func(pinnedmodel, q_midp, q_ini=q_inip, q_mid=q_midp, q_end=q_endp) < 1.0

    qp = q_inip
    qdp = (q_endp - q_inip)/10
    norm(sd_func(pinnedmodel, qp, qdp, q_ini=q_inip, q_mid=q_midp, q_end=q_endp) - 0.10) < 3e-2

end
