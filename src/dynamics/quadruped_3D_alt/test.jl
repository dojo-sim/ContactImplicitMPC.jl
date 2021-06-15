include(joinpath("..", "quadruped.jl"))

# test kinematics
q = rand(nq)

@assert norm(kinematics_1(model, q, body = :torso, mode = :ee) - [q[1] + model.l_torso * sin(q[3]); q[2] - model.l_torso * cos(q[3])]) ≈ 0.0
@assert norm(kinematics_1(model, q, body = :torso, mode = :com) - [q[1] + model.d_torso * sin(q[3]); q[2] - model.d_torso * cos(q[3])]) ≈ 0.0
@assert norm(kinematics_1(model, q, body = :thigh_1, mode = :ee) - [q[1] + model.l_thigh1 * sin(q[4]); q[2] - model.l_thigh1 * cos(q[4])]) ≈ 0.0
@assert norm(kinematics_1(model, q, body = :thigh_1, mode = :com) - [q[1] + model.d_thigh1 * sin(q[4]); q[2] - model.d_thigh1 * cos(q[4])]) ≈ 0.0
@assert norm(kinematics_1(model, q, body = :thigh_2, mode = :ee) - [q[1] + model.l_thigh2 * sin(q[6]); q[2] - model.l_thigh2 * cos(q[6])]) ≈ 0.0
@assert norm(kinematics_1(model, q, body = :thigh_2, mode = :com) - [q[1] + model.d_thigh2 * sin(q[6]); q[2] - model.d_thigh2 * cos(q[6])]) ≈ 0.0

@assert norm(kinematics_2(model, q, body = :calf_1, mode = :ee) - [q[1] + model.l_thigh1 * sin(q[4]) + model.l_calf1 * sin(q[5]); q[2] - model.l_thigh1 * cos(q[4]) - model.l_calf1 * cos(q[5])]) ≈ 0.0
@assert norm(kinematics_2(model, q, body = :calf_1, mode = :com) - [q[1] + model.l_thigh1 * sin(q[4]) + model.d_calf1 * sin(q[5]); q[2] - model.l_thigh1 * cos(q[4]) - model.d_calf1 * cos(q[5])]) ≈ 0.0
@assert norm(kinematics_2(model, q, body = :calf_2, mode = :ee) - [q[1] + model.l_thigh2 * sin(q[6]) + model.l_calf2 * sin(q[7]); q[2] - model.l_thigh2 * cos(q[6]) - model.l_calf2 * cos(q[7])]) ≈ 0.0
@assert norm(kinematics_2(model, q, body = :calf_2, mode = :com) - [q[1] + model.l_thigh2 * sin(q[6]) + model.d_calf2 * sin(q[7]); q[2] - model.l_thigh2 * cos(q[6]) - model.d_calf2 * cos(q[7])]) ≈ 0.0
@assert norm(kinematics_2(model, q, body = :thigh_3, mode = :ee) - [q[1] + model.l_torso * sin(q[3]) + model.l_thigh3 * sin(q[8]); q[2] - model.l_torso * cos(q[3]) - model.l_thigh3 * cos(q[8])]) ≈ 0.0
@assert norm(kinematics_2(model, q, body = :thigh_3, mode = :com) - [q[1] + model.l_torso * sin(q[3]) + model.d_thigh3 * sin(q[8]); q[2] - model.l_torso * cos(q[3]) - model.d_thigh3 * cos(q[8])]) ≈ 0.0
@assert norm(kinematics_2(model, q, body = :thigh_4, mode = :ee) - [q[1] + model.l_torso * sin(q[3]) + model.l_thigh4 * sin(q[10]); q[2] - model.l_torso * cos(q[3]) - model.l_thigh4 * cos(q[10])]) ≈ 0.0
@assert norm(kinematics_2(model, q, body = :thigh_4, mode = :com) - [q[1] + model.l_torso * sin(q[3]) + model.d_thigh4 * sin(q[10]); q[2] - model.l_torso * cos(q[3]) - model.d_thigh4 * cos(q[10])]) ≈ 0.0

@assert norm(kinematics_3(model, q, body = :calf_3, mode = :ee) - [q[1] + model.l_torso * sin(q[3]) + model.l_thigh3 * sin(q[8]) + model.l_calf3 * sin(q[9]); q[2] - model.l_torso * cos(q[3]) - model.l_thigh3 * cos(q[8]) - model.l_calf3 * cos(q[9])]) ≈ 0.0
@assert norm(kinematics_3(model, q, body = :calf_3, mode = :com) - [q[1] + model.l_torso * sin(q[3]) + model.l_thigh3 * sin(q[8]) + model.d_calf3 * sin(q[9]); q[2] - model.l_torso * cos(q[3]) - model.l_thigh3 * cos(q[8]) - model.d_calf3 * cos(q[9])]) ≈ 0.0
@assert norm(kinematics_3(model, q, body = :calf_4, mode = :ee) - [q[1] + model.l_torso * sin(q[3]) + model.l_thigh4 * sin(q[10]) + model.l_calf4 * sin(q[11]); q[2] - model.l_torso * cos(q[3]) - model.l_thigh4 * cos(q[10]) - model.l_calf4 * cos(q[11])]) ≈ 0.0
@assert norm(kinematics_3(model, q, body = :calf_4, mode = :com) - [q[1] + model.l_torso * sin(q[3]) + model.l_thigh4 * sin(q[10]) + model.d_calf4 * sin(q[11]); q[2] - model.l_torso * cos(q[3]) - model.l_thigh4 * cos(q[10]) - model.d_calf4 * cos(q[11])]) ≈ 0.0

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
