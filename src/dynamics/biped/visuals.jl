# visualization
function visualize!(vis, model::Biped, q;
      r = 0.040, Δt = 0.1, name::Symbol=:biped, camera::Bool=false)

	default_background!(vis)

	torso = Cylinder(Point3f0(0.0, 0.0, 0.0), Point3f0(0.0, 0.0, model.l_torso),
		convert(Float32, 0.035))
	setobject!(vis[name]["torso"], torso,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_1 = Cylinder(Point3f0(0.0,0.0,0.0), Point3f0(0.0, 0.0, model.l_thigh1),
		convert(Float32, 0.035))
	setobject!(vis[name]["thigh1"], thigh_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_1 = Cylinder(Point3f0(0.0,0.0,0.0), Point3f0(0.0, 0.0, model.l_calf1),
		convert(Float32, 0.035))
	setobject!(vis[name]["calf1"], calf_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	foot_1 = Cylinder(Point3f0(0.0,0.0,0.0),
		Point3f0(0.0, 0.0, model.l_foot1 + model.d_foot1),
		convert(Float32, 0.035))
	setobject!(vis[name]["foot1"], foot_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_2 = Cylinder(Point3f0(0.0,0.0,0.0), Point3f0(0.0, 0.0, model.l_thigh2),
		convert(Float32, 0.035))
	setobject!(vis[name]["thigh2"], thigh_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_2 = Cylinder(Point3f0(0.0,0.0,0.0), Point3f0(0.0, 0.0, model.l_calf2),
		convert(Float32, 0.035))
	setobject!(vis[name]["calf2"], calf_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	foot_2 = Cylinder(Point3f0(0.0,0.0,0.0),
		Point3f0(0.0, 0.0, model.l_foot2 + model.d_foot2),
		convert(Float32, 0.035))
	setobject!(vis[name]["foot2"], foot_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	setobject!(vis[name]["heel1"], Sphere(Point3f0(0.0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))
	setobject!(vis[name]["heel2"], Sphere(Point3f0(0.0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))
	setobject!(vis[name]["toe1"], Sphere(Point3f0(0.0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))
	setobject!(vis[name]["toe2"], Sphere(Point3f0(0.0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))
	setobject!(vis[name]["knee1"], Sphere(Point3f0(0.0),
		convert(Float32, 0.035)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))
	setobject!(vis[name]["knee2"], Sphere(Point3f0(0.0),
		convert(Float32, 0.035)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))
	setobject!(vis[name]["hip"], Sphere(Point3f0(0.0),
		convert(Float32, 0.035)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))
	setobject!(vis[name]["torso_top"], Sphere(Point3f0(0.0),
		convert(Float32, 0.035)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

	T = length(q)
	p_shift = [0.0; 0.0; r]
	for t = 1:T
		MeshCat.atframe(anim, t) do
			p = [q[t][1]; 0.0; q[t][2]] + p_shift

			k_torso = kinematics_1(model, q[t], body = :torso, mode = :ee)
			p_torso = [k_torso[1], 0.0, k_torso[2]] + p_shift

			k_thigh_1 = kinematics_1(model, q[t], body = :thigh_1, mode = :ee)
			p_thigh_1 = [k_thigh_1[1], 0.0, k_thigh_1[2]] + p_shift

			k_calf_1 = kinematics_2(model, q[t], body = :calf_1, mode = :ee)
			p_calf_1 = [k_calf_1[1], 0.0, k_calf_1[2]] + p_shift

			k_thigh_2 = kinematics_1(model, q[t], body = :thigh_2, mode = :ee)
			p_thigh_2 = [k_thigh_2[1], 0.0, k_thigh_2[2]] + p_shift

			k_calf_2 = kinematics_2(model, q[t], body = :calf_2, mode = :ee)
			p_calf_2 = [k_calf_2[1], 0.0, k_calf_2[2]] + p_shift

			k_toe_1 = kinematics_3(model, q[t], body = :foot_1, mode = :toe)
			p_toe_1 = [k_toe_1[1], 0.0, k_toe_1[2]] + p_shift

			k_heel_1 = kinematics_3(model, q[t], body = :foot_1, mode = :heel)
			p_heel_1 = [k_heel_1[1], 0.0, k_heel_1[2]] + p_shift

			k_toe_2 = kinematics_3(model, q[t], body = :foot_2, mode = :toe)
			p_toe_2 = [k_toe_2[1], 0.0, k_toe_2[2]] + p_shift

			k_heel_2 = kinematics_3(model, q[t], body = :foot_2, mode = :heel)
			p_heel_2 = [k_heel_2[1], 0.0, k_heel_2[2]] + p_shift

			settransform!(vis[name]["thigh1"], cable_transform(p, p_thigh_1))
			settransform!(vis[name]["calf1"], cable_transform(p_thigh_1, p_calf_1))
			settransform!(vis[name]["foot1"], cable_transform(p_toe_1, p_heel_1))

			settransform!(vis[name]["thigh2"], cable_transform(p, p_thigh_2))
			settransform!(vis[name]["calf2"], cable_transform(p_thigh_2, p_calf_2))
			settransform!(vis[name]["foot2"], cable_transform(p_toe_2, p_heel_2))

			settransform!(vis[name]["torso"], cable_transform(p_torso,p))
			settransform!(vis[name]["heel1"], Translation(p_heel_1))
			settransform!(vis[name]["heel2"], Translation(p_heel_2))
			settransform!(vis[name]["toe1"], Translation(p_toe_1))
			settransform!(vis[name]["toe2"], Translation(p_toe_2))
			settransform!(vis[name]["knee1"], Translation(p_thigh_1))
			settransform!(vis[name]["knee2"], Translation(p_thigh_2))
			settransform!(vis[name]["hip"], Translation(p))
			settransform!(vis[name]["torso_top"], Translation(p_torso))
		end
	end
	if camera
		settransform!(vis["/Cameras/default"],
		    compose(Translation(0.0, 0.0, -1.0), LinearMap(RotZ(-pi / 2.0))))
	end

	MeshCat.setanimation!(vis, anim)
end
