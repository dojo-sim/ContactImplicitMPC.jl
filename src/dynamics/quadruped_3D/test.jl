include(joinpath(@__DIR__, "model.jl"))

model = quadruped_3d

x = [1., 0, 0]
y = [0., 1, 0]
z = [0., 0, 1]

@testset "Quadruped 3D Kinematics" begin
    # test kinematics
    q = [rand(3); 0.0; 0.4142138; 0.0; π/2;0;0; π/2;0;0; π/2;0;0; π/2;0;0;]
    p = q[1:3]

    # TORSO
    ptorso_ee_ = p + model.l_torso*x
    ptorso_ee = kinematics_1(model, q, body = :torso, mode = :ee)
    @test norm(ptorso_ee - ptorso_ee_) < 1e-5
    ptorso_com_ = p + model.d_torso*x
    ptorso_com = kinematics_1(model, q, body = :torso, mode = :com)
    @test norm(ptorso_com - ptorso_com_) < 1e-5

    # SHOULDER1
    pshoulder_1ee_ = p + model.shoulder_lateral_offset*y + model.l_shoulder1*z
    pshoulder_1ee = kinematics_1(model, q, body=:shoulder_1, mode=:ee)
    @test norm(pshoulder_1ee - pshoulder_1ee_) < 1e-5
    pshoulder_1com_ = p + model.shoulder_lateral_offset*y + model.d_shoulder1*z
    pshoulder_1com = kinematics_1(model, q, body=:shoulder_1, mode=:com)
    @test norm(pshoulder_1com - pshoulder_1com_) < 1e-5

    # SHOULDER2
    pshoulder_2ee_ = p - model.shoulder_lateral_offset*y + model.l_shoulder1*z
    pshoulder_2ee = kinematics_1(model, q, body=:shoulder_2, mode=:ee)
    @test norm(pshoulder_2ee - pshoulder_2ee_) < 1e-5
    pshoulder_2com_ = p - model.shoulder_lateral_offset*y + model.d_shoulder1*z
    pshoulder_2com = kinematics_1(model, q, body=:shoulder_2, mode=:com)
    @test norm(pshoulder_2com - pshoulder_2com_) < 1e-5

    # SHOULDER3
    pshoulder_3ee_ = ptorso_ee + model.shoulder_lateral_offset*y + model.l_shoulder1*z
    pshoulder_3ee = kinematics_2(model, q, body=:shoulder_3, mode=:ee)
    @test norm(pshoulder_3ee - pshoulder_3ee_) < 1e-5
    pshoulder_3com_ = ptorso_ee + model.shoulder_lateral_offset*y + model.d_shoulder1*z
    pshoulder_3com = kinematics_2(model, q, body=:shoulder_3, mode=:com)
    @test norm(pshoulder_3com - pshoulder_3com_) < 1e-5

    # SHOULDER4
    pshoulder_4ee_ = ptorso_ee - model.shoulder_lateral_offset*y + model.l_shoulder1*z
    pshoulder_4ee = kinematics_2(model, q, body=:shoulder_4, mode=:ee)
    @test norm(pshoulder_4ee - pshoulder_4ee_) < 1e-5
    pshoulder_4com_ = ptorso_ee - model.shoulder_lateral_offset*y + model.d_shoulder1*z
    pshoulder_4com = kinematics_2(model, q, body=:shoulder_4, mode=:com)
    @test norm(pshoulder_4com - pshoulder_4com_) < 1e-5

    # THIGH1
    pthigh_1ee_ = pshoulder_1ee + model.l_thigh1*y
    pthigh_1ee = kinematics_2(model, q, body=:thigh_1, mode=:ee)
    @test norm(pthigh_1ee - pthigh_1ee_) < 1e-5
    pthigh_1com_ = pshoulder_1ee + model.d_thigh1*y
    pthigh_1com = kinematics_2(model, q, body=:thigh_1, mode=:com)
    @test norm(pthigh_1com - pthigh_1com_) < 1e-5

    # THIGH2
    pthigh_2ee_ = pshoulder_2ee - model.l_thigh1*y
    pthigh_2ee = kinematics_2(model, q, body=:thigh_2, mode=:ee)
    @test norm(pthigh_2ee - pthigh_2ee_) < 1e-5
    pthigh_2com_ = pshoulder_2ee - model.d_thigh1*y
    pthigh_2com = kinematics_2(model, q, body=:thigh_2, mode=:com)
    @test norm(pthigh_2com - pthigh_2com_) < 1e-5

    # THIGH3
    pthigh_3ee_ = pshoulder_3ee + model.l_thigh3*y
    pthigh_3ee = kinematics_3(model, q, body=:thigh_3, mode=:ee)
    @test norm(pthigh_3ee - pthigh_3ee_) < 1e-5
    pthigh_3com_ = pshoulder_3ee + model.d_thigh3*y
    pthigh_3com = kinematics_3(model, q, body=:thigh_3, mode=:com)
    @test norm(pthigh_3com - pthigh_3com_) < 1e-5

    # THIGH4
    pthigh_4ee_ = pshoulder_4ee - model.l_thigh4*y
    pthigh_4ee = kinematics_3(model, q, body=:thigh_4, mode=:ee)
    @test norm(pthigh_4ee - pthigh_4ee_) < 1e-5
    pthigh_4com_ = pshoulder_4ee - model.d_thigh4*y
    pthigh_4com = kinematics_3(model, q, body=:thigh_4, mode=:com)
    @test norm(pthigh_4com - pthigh_4com_) < 1e-5

    # CALF1
    pcalf_1ee_ = pthigh_1ee + model.l_calf1*y
    pcalf_1ee = kinematics_3(model, q, body=:calf_1, mode=:ee)
    @test norm(pcalf_1ee - pcalf_1ee_) < 1e-5
    pcalf_1com_ = pthigh_1ee + model.d_calf1*y
    pcalf_1com = kinematics_3(model, q, body=:calf_1, mode=:com)
    @test norm(pcalf_1com - pcalf_1com_) < 1e-5

    # CALF2
    pcalf_2ee_ = pthigh_2ee - model.l_calf1*y
    pcalf_2ee = kinematics_3(model, q, body=:calf_2, mode=:ee)
    @test norm(pcalf_2ee - pcalf_2ee_) < 1e-5
    pcalf_2com_ = pthigh_2ee - model.d_calf1*y
    pcalf_2com = kinematics_3(model, q, body=:calf_2, mode=:com)
    @test norm(pcalf_2com - pcalf_2com_) < 1e-5

    # CALF3
    pcalf_3ee_ = pthigh_3ee + model.l_calf1*y
    pcalf_3ee = kinematics_4(model, q, body=:calf_3, mode=:ee)
    @test norm(pcalf_3ee - pcalf_3ee_) < 1e-5
    pcalf_3com_ = pthigh_3ee + model.d_calf1*y
    pcalf_3com = kinematics_4(model, q, body=:calf_3, mode=:com)
    @test norm(pcalf_3com - pcalf_3com_) < 1e-5

    # CALF4
    pcalf_4ee_ = pthigh_4ee - model.l_calf1*y
    pcalf_4ee = kinematics_4(model, q, body=:calf_4, mode=:ee)
    @test norm(pcalf_4ee - pcalf_4ee_) < 1e-5
    pcalf_4com_ = pthigh_4ee - model.d_calf1*y
    pcalf_4com = kinematics_4(model, q, body=:calf_4, mode=:com)
    @test norm(pcalf_4com - pcalf_4com_) < 1e-5
end

@testset "Quadruped 3D Jacobians" begin
    # test Jacobians
    q = rand(model.dim.q)

    # TORSO
    Jtorso_ee_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body=:torso, mode=:ee), q)
    Jtorso_ee = jacobian_1(model, q, body = :torso, mode = :ee)
    @test norm(Jtorso_ee - Jtorso_ee_) < 1e-5
    Jtorso_com_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body=:torso, mode=:com), q)
    Jtorso_com = jacobian_1(model, q, body = :torso, mode = :com)
    @test norm(Jtorso_com - Jtorso_com_) < 1e-5

    # SHOULDER1
    Jshoulder_1ee_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body=:shoulder_1, mode=:ee), q)
    Jshoulder_1ee = jacobian_1(model, q, body = :shoulder_1, mode = :ee)
    @test norm(Jshoulder_1ee - Jshoulder_1ee_) < 1e-5
    Jshoulder_1com_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body=:shoulder_1, mode=:com), q)
    Jshoulder_1com = jacobian_1(model, q, body = :shoulder_1, mode = :com)
    @test norm(Jshoulder_1com - Jshoulder_1com_) < 1e-5

    # SHOULDER2
    Jshoulder_2ee_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body=:shoulder_2, mode=:ee), q)
    Jshoulder_2ee = jacobian_1(model, q, body = :shoulder_2, mode = :ee)
    @test norm(Jshoulder_2ee - Jshoulder_2ee_) < 1e-5
    Jshoulder_2com_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body=:shoulder_2, mode=:com), q)
    Jshoulder_2com = jacobian_1(model, q, body = :shoulder_2, mode = :com)
    @test norm(Jshoulder_2com - Jshoulder_2com_) < 1e-5

    # SHOULDER3
    Jshoulder_3ee_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body=:shoulder_3, mode=:ee), q)
    Jshoulder_3ee = jacobian_2(model, q, body = :shoulder_3, mode = :ee)
    @test norm(Jshoulder_3ee - Jshoulder_3ee_) < 1e-5
    Jshoulder_3com_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body=:shoulder_3, mode=:com), q)
    Jshoulder_3com = jacobian_2(model, q, body = :shoulder_3, mode = :com)
    @test norm(Jshoulder_3com - Jshoulder_3com_) < 1e-5

    # SHOULDER4
    Jshoulder_4ee_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body=:shoulder_4, mode=:ee), q)
    Jshoulder_4ee = jacobian_2(model, q, body = :shoulder_4, mode = :ee)
    @test norm(Jshoulder_4ee - Jshoulder_4ee_) < 1e-5
    Jshoulder_4com_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body=:shoulder_4, mode=:com), q)
    Jshoulder_4com = jacobian_2(model, q, body = :shoulder_4, mode = :com)
    @test norm(Jshoulder_4com - Jshoulder_4com_) < 1e-5

    # THIGH1
    Jthigh_1ee_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body=:thigh_1, mode=:ee), q)
    Jthigh_1ee = jacobian_2(model, q, body = :thigh_1, mode = :ee)
    @test norm(Jthigh_1ee - Jthigh_1ee_) < 1e-5
    Jthigh_1com_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body=:thigh_1, mode=:com), q)
    Jthigh_1com = jacobian_2(model, q, body = :thigh_1, mode = :com)
    @test norm(Jthigh_1com - Jthigh_1com_) < 1e-5

    # THIGH2
    Jthigh_2ee_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body=:thigh_2, mode=:ee), q)
    Jthigh_2ee = jacobian_2(model, q, body = :thigh_2, mode = :ee)
    @test norm(Jthigh_2ee - Jthigh_2ee_) < 1e-5
    Jthigh_2com_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body=:thigh_2, mode=:com), q)
    Jthigh_2com = jacobian_2(model, q, body = :thigh_2, mode = :com)
    @test norm(Jthigh_2com - Jthigh_2com_) < 1e-5

    # THIGH3
    Jthigh_3ee_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body=:thigh_3, mode=:ee), q)
    Jthigh_3ee = jacobian_3(model, q, body = :thigh_3, mode = :ee)
    @test norm(Jthigh_3ee - Jthigh_3ee_) < 1e-5
    Jthigh_3com_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body=:thigh_3, mode=:com), q)
    Jthigh_3com = jacobian_3(model, q, body = :thigh_3, mode = :com)
    @test norm(Jthigh_3com - Jthigh_3com_) < 1e-5

    # THIGH4
    Jthigh_4ee_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body=:thigh_4, mode=:ee), q)
    Jthigh_4ee = jacobian_3(model, q, body = :thigh_4, mode = :ee)
    @test norm(Jthigh_4ee - Jthigh_4ee_) < 1e-5
    Jthigh_4com_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body=:thigh_4, mode=:com), q)
    Jthigh_4com = jacobian_3(model, q, body = :thigh_4, mode = :com)
    @test norm(Jthigh_4com - Jthigh_4com_) < 1e-5


    # CALF1
    Jcalf_1ee_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body=:calf_1, mode=:ee), q)
    Jcalf_1ee = jacobian_3(model, q, body = :calf_1, mode = :ee)
    @test norm(Jcalf_1ee - Jcalf_1ee_) < 1e-5
    Jcalf_1com_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body=:calf_1, mode=:com), q)
    Jcalf_1com = jacobian_3(model, q, body = :calf_1, mode = :com)
    @test norm(Jcalf_1com - Jcalf_1com_) < 1e-5

    # CALF2
    Jcalf_2ee_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body=:calf_2, mode=:ee), q)
    Jcalf_2ee = jacobian_3(model, q, body = :calf_2, mode = :ee)
    @test norm(Jcalf_2ee - Jcalf_2ee_) < 1e-5
    Jcalf_2com_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body=:calf_2, mode=:com), q)
    Jcalf_2com = jacobian_3(model, q, body = :calf_2, mode = :com)
    @test norm(Jcalf_2com - Jcalf_2com_) < 1e-5

    # CALF3
    Jcalf_3ee_ = ForwardDiff.jacobian(q -> kinematics_4(model, q, body=:calf_3, mode=:ee), q)
    Jcalf_3ee = jacobian_4(model, q, body = :calf_3, mode = :ee)
    @test norm(Jcalf_3ee - Jcalf_3ee_) < 1e-5
    Jcalf_3com_ = ForwardDiff.jacobian(q -> kinematics_4(model, q, body=:calf_3, mode=:com), q)
    Jcalf_3com = jacobian_4(model, q, body = :calf_3, mode = :com)
    @test norm(Jcalf_3com - Jcalf_3com_) < 1e-5

    # CALF4
    Jcalf_4ee_ = ForwardDiff.jacobian(q -> kinematics_4(model, q, body=:calf_4, mode=:ee), q)
    Jcalf_4ee = jacobian_4(model, q, body = :calf_4, mode = :ee)
    @test norm(Jcalf_4ee - Jcalf_4ee_) < 1e-5
    Jcalf_4com_ = ForwardDiff.jacobian(q -> kinematics_4(model, q, body=:calf_4, mode=:com), q)
    Jcalf_4com = jacobian_4(model, q, body = :calf_4, mode = :com)
    @test norm(Jcalf_4com - Jcalf_4com_) < 1e-5
end

q̇ = rand(nq)


a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
a = 10



k1(z) = kinematics_1(model, z, body = :torso, mode = :ee)
@assert norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :torso, mode = :ee)) ≈ 0.0
k1(z) = kinematics_1(model, z, body = :torso, mode = :com)
@assert norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :torso, mode = :com)) ≈ 0.0
k1(z) = kinematics_1(model, z, body = :thigh_1, mode = :ee)
@assert norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :thigh_1, mode = :ee)) ≈ 0.0
k1(z) = kinematics_1(model, z, body = :thigh_1, mode = :com)
@assert norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :thigh_1, mode = :com)) ≈ 0.0
k1(z) = kinematics_1(model, z, body = :thigh_2, mode = :ee)
@assert norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :thigh_2, mode = :ee)) ≈ 0.0
k1(z) = kinematics_1(model, z, body = :thigh_2, mode = :com)
@assert norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :thigh_2, mode = :com)) ≈ 0.0

k2(z) = kinematics_2(model, z, body = :calf_1, mode = :ee)
@assert norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :calf_1, mode = :ee)) ≈ 0.0
k2(z) = kinematics_2(model, z, body = :calf_1, mode = :com)
@assert norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :calf_1, mode = :com)) ≈ 0.0
k2(z) = kinematics_2(model, z, body = :calf_2, mode = :ee)
@assert norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :calf_2, mode = :ee)) ≈ 0.0
k2(z) = kinematics_2(model, z, body = :calf_2, mode = :com)
@assert norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :calf_2, mode = :com)) ≈ 0.0
k2(z) = kinematics_2(model, z, body = :thigh_3, mode = :ee)
@assert norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :thigh_3, mode = :ee)) ≈ 0.0
k2(z) = kinematics_2(model, z, body = :thigh_3, mode = :com)
@assert norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :thigh_3, mode = :com)) ≈ 0.0
k2(z) = kinematics_2(model, z, body = :thigh_4, mode = :ee)
@assert norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :thigh_4, mode = :ee)) ≈ 0.0
k2(z) = kinematics_2(model, z, body = :thigh_4, mode = :com)
@assert norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :thigh_4, mode = :com)) ≈ 0.0

k3(z) = kinematics_3(model, z, body = :calf_3, mode = :ee)
@assert norm(ForwardDiff.jacobian(k3, q) - jacobian_3(model, q, body = :calf_3, mode = :ee)) ≈ 0.0
k3(z) = kinematics_3(model, z, body = :calf_3, mode = :com)
@assert norm(ForwardDiff.jacobian(k3, q) - jacobian_3(model, q, body = :calf_3, mode = :com)) ≈ 0.0
k3(z) = kinematics_3(model, z, body = :calf_4, mode = :ee)
@assert norm(ForwardDiff.jacobian(k3, q) - jacobian_3(model, q, body = :calf_4, mode = :ee)) ≈ 0.0
k3(z) = kinematics_3(model, z, body = :calf_4, mode = :com)
@assert norm(ForwardDiff.jacobian(k3, q) - jacobian_3(model, q, body = :calf_4, mode = :com)) ≈ 0.0

q̇ = rand(nq)
eigen(M_func(model, q))

tmp_q(z) = _dLdq̇(model, z, q̇)
tmp_q̇(z) = _dLdq̇(model, q, z)
norm(ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))

norm((ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))[:,1:3])
norm((ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))[:,4:6])
norm((ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))[:,7:9])
norm((ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))[:,10:12])
norm((ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))[:,13:15])
norm((ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))[:,16:18])
