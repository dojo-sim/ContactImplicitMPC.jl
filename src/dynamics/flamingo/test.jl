include(joinpath("..", "walker.jl"))

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

norm(kinematics_3(model, q, body = :foot_1, mode = :toe)
	- [q[1] + model.l_thigh1 * sin(q[4]) + model.l_calf1 * sin(q[5]) + model.l_foot1 * sin(q[8]);
	q[2] - model.l_thigh1 * cos(q[4]) - model.l_calf1 * cos(q[5]) - model.l_foot1 * cos(q[8])])
norm(kinematics_3(model, q, body = :foot_1, mode = :heel)
	- [q[1] + model.l_thigh1 * sin(q[4]) + model.l_calf1 * sin(q[5]) - model.d_foot1 * sin(q[8]);
	q[2] - model.l_thigh1 * cos(q[4]) - model.l_calf1 * cos(q[5]) + model.d_foot1 * cos(q[8])])
norm(kinematics_3(model, q, body = :foot_1, mode = :com)
	- [q[1] + model.l_thigh1 * sin(q[4]) + model.l_calf1 * sin(q[5]) + 0.5 * (model.l_foot1 - model.d_foot1) * sin(q[8]);
	q[2] - model.l_thigh1 * cos(q[4]) - model.l_calf1 * cos(q[5]) - 0.5 * (model.l_foot1 - model.d_foot1) * cos(q[8])])

norm(kinematics_3(model, q, body = :foot_2, mode = :toe)
	- [q[1] + model.l_thigh2 * sin(q[6]) + model.l_calf2 * sin(q[7]) + model.l_foot2 * sin(q[9]);
	 q[2] - model.l_thigh2 * cos(q[6]) - model.l_calf2 * cos(q[7]) - model.l_foot2 * cos(q[9])])
norm(kinematics_3(model, q, body = :foot_2, mode = :heel)
	- [q[1] + model.l_thigh2 * sin(q[6]) + model.l_calf2 * sin(q[7]) - model.d_foot2 * sin(q[9]);
	q[2] - model.l_thigh2 * cos(q[6]) - model.l_calf2 * cos(q[7]) + model.d_foot2 * cos(q[9])])
norm(kinematics_3(model, q, body = :foot_2, mode = :com)
	- [q[1] + model.l_thigh2 * sin(q[6]) + model.l_calf2 * sin(q[7]) + 0.5 * (model.l_foot2 -  model.d_foot2) * sin(q[9]);
	q[2] - model.l_thigh2 * cos(q[6]) - model.l_calf2 * cos(q[7]) - 0.5 * (model.l_foot2 - model.d_foot2) * cos(q[9])])

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

k3(z) = kinematics_3(model, z, body = :foot_1, mode = :toe)
norm(ForwardDiff.jacobian(k3,q) - jacobian_3(model, q, body = :foot_1, mode = :toe))

k3(z) = kinematics_3(model, z, body = :foot_1, mode = :heel)
norm(ForwardDiff.jacobian(k3,q) - jacobian_3(model, q, body = :foot_1, mode = :heel))

k3(z) = kinematics_3(model, z, body = :foot_1, mode = :com)
norm(ForwardDiff.jacobian(k3,q) - jacobian_3(model, q, body = :foot_1, mode = :com))

k3(z) = kinematics_3(model, z, body = :foot_2, mode = :toe)
norm(ForwardDiff.jacobian(k3,q) - jacobian_3(model, q, body = :foot_2, mode = :toe))

k3(z) = kinematics_3(model, z, body = :foot_2, mode = :heel)
norm(ForwardDiff.jacobian(k3,q) - jacobian_3(model, q, body = :foot_2, mode = :heel))

k3(z) = kinematics_3(model, z, body = :foot_2, mode = :com)
norm(ForwardDiff.jacobian(k3,q) - jacobian_3(model, q, body = :foot_2, mode = :com))

q̇ = rand(nq)
lagrangian(model, q, q̇)
dLdq(model, q, q̇)
dLdq̇(model, q, q̇)

eigen(M_func(model, q))

tmp_q(z) = dLdq̇(model, z, q̇)
tmp_q̇(z) = dLdq̇(model, q, z)
norm(ForwardDiff.jacobian(tmp_q̇,q̇) - M_func(model, q))
