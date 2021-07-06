model = deepcopy(racecar)
q = rand(model.dim.q)

k1_1 = kinematics_1(model, q, body = :wheel1)
k1_2 = kinematics_1(model, q, body = :wheel2)
k1_3 = kinematics_1(model, q, body = :wheel3)
k1_4 = kinematics_1(model, q, body = :wheel4)

k2_1 = kinematics_2(model, q, body = :wheel_hub1)
k2_2 = kinematics_2(model, q, body = :wheel_hub2)
k2_3 = kinematics_2(model, q, body = :wheel_hub3)
k2_4 = kinematics_2(model, q, body = :wheel_hub4)

k3_1 = kinematics_3(model, q, body = :wheel_contact1)
k3_2 = kinematics_3(model, q, body = :wheel_contact2)
k3_3 = kinematics_3(model, q, body = :wheel_contact3)
k3_4 = kinematics_3(model, q, body = :wheel_contact4)

J1_1 = jacobian_1(model, q, body = :wheel1)
J1_2 = jacobian_1(model, q, body = :wheel2)
J1_3 = jacobian_1(model, q, body = :wheel3)
J1_4 = jacobian_1(model, q, body = :wheel4)

J2_1 = jacobian_2(model, q, body = :wheel_hub1)
J2_2 = jacobian_2(model, q, body = :wheel_hub2)
J2_3 = jacobian_2(model, q, body = :wheel_hub3)
J2_4 = jacobian_2(model, q, body = :wheel_hub4)

J3_1 = jacobian_3(model, q, body = :wheel_contact1)
J3_2 = jacobian_3(model, q, body = :wheel_contact2)
J3_3 = jacobian_3(model, q, body = :wheel_contact3)
J3_4 = jacobian_3(model, q, body = :wheel_contact4)

J1_1_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body = :wheel1), q)
J1_2_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body = :wheel2), q)
J1_3_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body = :wheel3), q)
J1_4_ = ForwardDiff.jacobian(q -> kinematics_1(model, q, body = :wheel4), q)
@test norm(J1_1 - J1_1_, Inf) < 1e-10
@test norm(J1_2 - J1_2_, Inf) < 1e-10
@test norm(J1_3 - J1_3_, Inf) < 1e-10
@test norm(J1_4 - J1_4_, Inf) < 1e-10

J2_1_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body = :wheel_hub1), q)
J2_2_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body = :wheel_hub2), q)
J2_3_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body = :wheel_hub3), q)
J2_4_ = ForwardDiff.jacobian(q -> kinematics_2(model, q, body = :wheel_hub4), q)
@test norm(J2_1 - J2_1_, Inf) < 1e-10
@test norm(J2_2 - J2_2_, Inf) < 1e-10
@test norm(J2_3 - J2_3_, Inf) < 1e-10
@test norm(J2_4 - J2_4_, Inf) < 1e-10

J3_1_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body = :wheel_contact1), q)
J3_2_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body = :wheel_contact2), q)
J3_3_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body = :wheel_contact3), q)
J3_4_ = ForwardDiff.jacobian(q -> kinematics_3(model, q, body = :wheel_contact4), q)
@test norm(J3_1 - J3_1_, Inf) < 1e-10
@test norm(J3_2 - J3_2_, Inf) < 1e-10
@test norm(J3_3 - J3_3_, Inf) < 1e-10
@test norm(J3_4 - J3_4_, Inf) < 1e-10

J3_1_[:,1:7]  - J3_1[:,1:7]
J3_1_[:,8:15] - J3_1[:,8:15]

J3_2_[:,1:7]  - J3_2[:,1:7]
J3_2_[:,8:15] - J3_2[:,8:15]

J3_3_[:,1:7]  - J3_3[:,1:7]
J3_3_[:,8:15] - J3_3[:,8:15]

J3_4_[:,1:7]  - J3_4[:,1:7]
J3_4_[:,8:15] - J3_4[:,8:15]
