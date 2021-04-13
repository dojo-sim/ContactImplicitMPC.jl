include(joinpath("..", "biped.jl"))

# test kinematics
q = rand(nq)

norm(kinematics_1(model, q, body = :torso, mode = :ee)
    - [q[1] - model.l_torso * sin(q[3]); q[2] + model.l_torso * cos(q[3])])
norm(kinematics_1(model, q, body = :torso, mode = :com)
    - [q[1] - model.d_torso * sin(q[3]); q[2] + model.d_torso * cos(q[3])])
norm(kinematics_1(model, q, body = :thigh_1, mode = :ee)
    - [q[1] + model.l_thigh1 * sin(q[4]); q[2] - model.l_thigh1 * cos(q[4])])
norm(kinematics_1(model, q, body = :thigh_1, mode = :com)
    - [q[1] + model.d_thigh1 * sin(q[4]); q[2] - model.d_thigh1 * cos(q[4])])
norm(kinematics_1(model, q, body = :thigh_2, mode = :ee)
    - [q[1] + model.l_thigh2 * sin(q[6]); q[2] - model.l_thigh2 * cos(q[6])])
norm(kinematics_1(model, q, body = :thigh_2, mode = :com)
    - [q[1] + model.d_thigh2 * sin(q[6]); q[2] - model.d_thigh2 * cos(q[6])])

norm(kinematics_2(model, q, body = :calf_1, mode = :ee)
    - [q[1] + model.l_thigh1 * sin(q[4]) + model.l_calf1 * sin(q[5]);
    q[2] - model.l_thigh1 * cos(q[4]) - model.l_calf1 * cos(q[5])])
norm(kinematics_2(model, q, body = :calf_1, mode = :com)
    - [q[1] + model.l_thigh1 * sin(q[4]) + model.d_calf1 * sin(q[5]);
    q[2] - model.l_thigh1 * cos(q[4]) - model.d_calf1 * cos(q[5])])
norm(kinematics_2(model, q, body = :calf_2, mode = :ee)
    - [q[1] + model.l_thigh2 * sin(q[6]) + model.l_calf2 * sin(q[7]);
    q[2] - model.l_thigh2 * cos(q[6]) - model.l_calf2 * cos(q[7])])
norm(kinematics_2(model, q, body = :calf_2, mode = :com)
    - [q[1] + model.l_thigh2 * sin(q[6]) + model.d_calf2 * sin(q[7]);
    q[2] - model.l_thigh2 * cos(q[6]) - model.d_calf2 * cos(q[7])])

k1(z) = kinematics_1(model, z, body = :torso, mode = :ee)
norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :torso, mode = :ee))

k1(z) = kinematics_1(model, z, body = :torso, mode = :com)
norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :torso, mode = :com))

k1(z) = kinematics_1(model, z, body = :thigh_1, mode = :ee)
norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :thigh_1, mode = :ee))

k1(z) = kinematics_1(model, z, body = :thigh_1, mode = :com)
norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :thigh_1, mode = :com))

k1(z) = kinematics_1(model, z, body = :thigh_2, mode = :ee)
norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :thigh_2, mode = :ee))

k1(z) = kinematics_1(model, z, body = :thigh_2, mode = :com)
norm(ForwardDiff.jacobian(k1,q) - jacobian_1(model, q, body = :thigh_2, mode = :com))

k2(z) = kinematics_2(model, z, body = :calf_1, mode = :ee)
norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :calf_1, mode = :ee))

k2(z) = kinematics_2(model, z, body = :calf_1, mode = :com)
norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :calf_1, mode = :com))

k2(z) = kinematics_2(model, z, body = :calf_2, mode = :ee)
norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :calf_2, mode = :ee))

k2(z) = kinematics_2(model, z, body = :calf_2, mode = :com)
norm(ForwardDiff.jacobian(k2,q) - jacobian_2(model, q, body = :calf_2, mode = :com))

q̇ = rand(nq)
lagrangian(model, q, q̇)
dLdq(model, q, q̇)
dLdq̇(model, q, q̇)

eigen(M_func(model, q))

tmp_q(z) = dLdq̇(model, z, q̇)
tmp_q̇(z) = dLdq̇(model, q, z)
norm(ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))
