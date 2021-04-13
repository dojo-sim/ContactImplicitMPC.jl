# visualization
function visualize!(vis, model::Biped5, q;
      r = 0.035, Δt = 0.1)

	default_background!(vis)

	torso = Cylinder(Point3f0(0.0, 0.0, 0.0), Point3f0(0.0, 0.0, model.l_torso),
		convert(Float32, 0.025))
	setobject!(vis["torso"], torso,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_1 = Cylinder(Point3f0(0.0,0.0,0.0), Point3f0(0.0, 0.0, model.l_thigh1),
		convert(Float32, 0.025))
	setobject!(vis["thigh1"], thigh_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_1 = Cylinder(Point3f0(0.0,0.0,0.0), Point3f0(0.0, 0.0, model.l_calf1),
		convert(Float32, 0.025))
	setobject!(vis["calf1"], calf_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_2 = Cylinder(Point3f0(0.0,0.0,0.0), Point3f0(0.0, 0.0, model.l_thigh2),
		convert(Float32, 0.025))
	setobject!(vis["thigh2"], thigh_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_2 = Cylinder(Point3f0(0.0,0.0,0.0), Point3f0(0.0, 0.0, model.l_calf2),
		convert(Float32, 0.025))
	setobject!(vis["calf2"], calf_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	setobject!(vis["foot1"], Sphere(Point3f0(0.0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))
	setobject!(vis["foot2"], Sphere(Point3f0(0.0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))
	setobject!(vis["knee1"], Sphere(Point3f0(0.0),
		convert(Float32, 0.025)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))
	setobject!(vis["knee2"], Sphere(Point3f0(0.0),
		convert(Float32, 0.025)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))
	setobject!(vis["hip"], Sphere(Point3f0(0.0),
		convert(Float32, 0.025)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))
	setobject!(vis["torso_top"], Sphere(Point3f0(0.0),
		convert(Float32, 0.025)),
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

			settransform!(vis["thigh1"], cable_transform(p, p_thigh_1))
			settransform!(vis["calf1"], cable_transform(p_thigh_1, p_calf_1))
			settransform!(vis["thigh2"], cable_transform(p, p_thigh_2))
			settransform!(vis["calf2"], cable_transform(p_thigh_2, p_calf_2))
			settransform!(vis["torso"], cable_transform(p_torso,p))
			settransform!(vis["foot1"], Translation(p_calf_1))
			settransform!(vis["foot2"], Translation(p_calf_2))
			settransform!(vis["knee1"], Translation(p_thigh_1))
			settransform!(vis["knee2"], Translation(p_thigh_2))
			settransform!(vis["hip"], Translation(p))
			settransform!(vis["torso_top"], Translation(p_torso))


		end
	end

	settransform!(vis["/Cameras/default"],
	   compose(Translation(0.0 , 1.0 , -1.0), LinearMap(RotZ(-pi / 2.0))))

	MeshCat.setanimation!(vis, anim)
end
