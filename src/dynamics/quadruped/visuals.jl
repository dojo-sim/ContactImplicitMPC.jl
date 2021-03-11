# visualization
function visualize!(vis, model::Quadruped, q;
	r = 0.025, Δt = 0.1, name::String="quadruped")

	default_background!(vis)

	torso = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0, 0.0, 0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_torso),
		convert(Float32, 0.035))
	setobject!(vis[name*"/torso"], torso,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_1 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_thigh1),
		convert(Float32, 0.0175))
	setobject!(vis[name*"/thigh1"], thigh_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_1 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_calf1),
		convert(Float32, 0.0125))
	setobject!(vis[name*"/leg1"], calf_1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_2 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_thigh2),
		convert(Float32, 0.0175))
	setobject!(vis[name*"/thigh2"], thigh_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_2 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_calf2),
		convert(Float32, 0.0125))
	setobject!(vis[name*"/leg2"], calf_2,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_3 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_thigh3),
		convert(Float32, 0.0175))
	setobject!(vis[name*"/thigh3"], thigh_3,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_3 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_calf3),
		convert(Float32, 0.0125))
	setobject!(vis[name*"/leg3"], calf_3,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	thigh_4 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_thigh4),
		convert(Float32, 0.0175))
	setobject!(vis[name*"/thigh4"], thigh_4,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	calf_4 = GeometryBasics.Cylinder(GeometryBasics.Point3f0(0.0,0.0,0.0),
		GeometryBasics.Point3f0(0.0, 0.0, model.l_calf4),
		convert(Float32, 0.0125))
	setobject!(vis[name*"/leg4"], calf_4,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

	hip1 = setobject!(vis[name*"/hip1"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, 0.035)),
        MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	hip2 = setobject!(vis[name*"/hip2"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, 0.035)),
        MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	knee1 = setobject!(vis[name*"/knee1"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, 0.025)),
        MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	knee2 = setobject!(vis[name*"/knee2"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, 0.025)),
		MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	knee3 = setobject!(vis[name*"/knee3"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, 0.025)),
		MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	knee4 = setobject!(vis[name*"/knee4"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, 0.025)),
        MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))

	feet1 = setobject!(vis[name*"/feet1"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet2 = setobject!(vis[name*"/feet2"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet3 = setobject!(vis[name*"/feet3"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet4 = setobject!(vis[name*"/feet4"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	T = length(q)
	p_shift = [0.0, 0.0, r]
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


			k_thigh_3 = kinematics_2(model, q[t], body = :thigh_3, mode = :ee)
			p_thigh_3 = [k_thigh_3[1], 0.0, k_thigh_3[2]] + p_shift

			k_calf_3 = kinematics_3(model, q[t], body = :calf_3, mode = :ee)
			p_calf_3 = [k_calf_3[1], 0.0, k_calf_3[2]] + p_shift

			k_thigh_4 = kinematics_2(model, q[t], body = :thigh_4, mode = :ee)
			p_thigh_4 = [k_thigh_4[1], 0.0, k_thigh_4[2]] + p_shift

			k_calf_4 = kinematics_3(model, q[t], body = :calf_4, mode = :ee)
			p_calf_4 = [k_calf_4[1], 0.0, k_calf_4[2]] + p_shift

			settransform!(vis[name*"/thigh1"], cable_transform(p, p_thigh_1))
			settransform!(vis[name*"/leg1"], cable_transform(p_thigh_1, p_calf_1))
			settransform!(vis[name*"/thigh2"], cable_transform(p, p_thigh_2))
			settransform!(vis[name*"/leg2"], cable_transform(p_thigh_2, p_calf_2))
			settransform!(vis[name*"/thigh3"], cable_transform(p_torso, p_thigh_3))
			settransform!(vis[name*"/leg3"], cable_transform(p_thigh_3, p_calf_3))
			settransform!(vis[name*"/thigh4"], cable_transform(p_torso, p_thigh_4))
			settransform!(vis[name*"/leg4"], cable_transform(p_thigh_4, p_calf_4))
			settransform!(vis[name*"/torso"], cable_transform(p, p_torso))
			settransform!(vis[name*"/hip1"], MeshCat.Translation(p))
			settransform!(vis[name*"/hip2"], MeshCat.Translation(p_torso))
			settransform!(vis[name*"/knee1"], MeshCat.Translation(p_thigh_1))
			settransform!(vis[name*"/knee2"], MeshCat.Translation(p_thigh_2))
			settransform!(vis[name*"/knee3"], MeshCat.Translation(p_thigh_3))
			settransform!(vis[name*"/knee4"], MeshCat.Translation(p_thigh_4))
			settransform!(vis[name*"/feet1"], MeshCat.Translation(p_calf_1))
			settransform!(vis[name*"/feet2"], MeshCat.Translation(p_calf_2))
			settransform!(vis[name*"/feet3"], MeshCat.Translation(p_calf_3))
			settransform!(vis[name*"/feet4"], MeshCat.Translation(p_calf_4))
		end
	end

	settransform!(vis["/Cameras/default"],
	    compose(Translation(0.0, 0.0, -1.0), LinearMap(RotZ(-pi / 2.0))))

	MeshCat.setanimation!(vis, anim)
end
