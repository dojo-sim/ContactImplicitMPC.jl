function visualize!(vis, model, q;
       Δt = 0.1,
	   color = RGBA(1,0,0,1.0))

	default_background!(vis)

	l2 = Cylinder(Point3f0(-model.l * 4.0, 0.0, 0.0),
		Point3f0(model.l * 4.0, 0.0, 0.0),
		convert(Float32, 0.025))

	setobject!(vis["slider"], l2, MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    l1 = Cylinder(Point3f0(0.0, 0.0, 0.0),
		Point3f0(0.0, 0.0, model.l),
		convert(Float32, 0.025))

    setobject!(vis["arm"], l1,
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    setobject!(vis["base"], HyperSphere(Point3f0(0.0),
        convert(Float32, 0.1)),
        MeshPhongMaterial(color = color))

    setobject!(vis["ee"], HyperSphere(Point3f0(0.0),
        convert(Float32, 0.05)),
        MeshPhongMaterial(color = color))

	anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

	for t = 1:length(q)
	    MeshCat.atframe(anim,t) do
			x = q[t]
			px = x[1] + model.l * sin(x[2])
			pz = -model.l * cos(x[2])
            settransform!(vis["arm"], cable_transform([x[1]; 0;0], [px; 0.0; pz]))
            settransform!(vis["base"], Translation([x[1]; 0.0; 0.0]))
            settransform!(vis["ee"], Translation([px; 0.0; pz]))
	    end
	end

	settransform!(vis["/Cameras/default"],
		compose(Translation(0.0, 0.0, -1.0), LinearMap(RotZ(- pi / 2))))

	MeshCat.setanimation!(vis,anim)
end
