function visualize!(vis, model::Hopper3D, q;
	Δt = 0.1)

	default_background!(vis)

	r_foot = 0.05
	setobject!(vis["body"], Sphere(Point3f0(0),
	   convert(Float32, 0.1)),
	   MeshPhongMaterial(color = RGBA(0, 1, 0, 1.0)))
	setobject!(vis["foot"], Sphere(Point3f0(0),
	   convert(Float32, r_foot)),
	   MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))

	r_leg = 0.5 * r_foot
	n_leg = 100
	for i = 1:n_leg
	   setobject!(vis["leg$i"], Sphere(Point3f0(0),
	       convert(Float32, r_leg)),
	       MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))
	end
	p_leg = [zeros(3) for i = 1:n_leg]
	anim = MeshCat.Animation(convert(Int,floor(1.0 / Δt)))

	z_shift = [0.0 ; 0.0; r_foot]
	for t = 1:length(q)
	   p_body = q[t][1:3]

	   q_tmp = Array(copy(q[t]))
	   r_range = range(0.0, stop = q[t][7], length = n_leg)
	   for i = 1:n_leg
	       q_tmp[7] = r_range[i]
	       p_leg[i] = kinematics(model, q_tmp)
	   end
	   q_tmp[7] = q[t][7]
	   p_foot = kinematics(model, q_tmp)

	   MeshCat.atframe(anim, t) do
	       settransform!(vis["body"], Translation(p_body + z_shift))
	       settransform!(vis["foot"], Translation(p_foot + z_shift))

	       for i = 1:n_leg
	           settransform!(vis["leg$i"], Translation(p_leg[i] + z_shift))
	       end
	   end
	end
	MeshCat.setanimation!(vis, anim)
end
