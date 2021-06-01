function visualize!(vis, model::Box, q;
        Δt = 0.1)

	default_background!(vis)

	r = abs(model.corner_offset[1][1])
    setobject!(vis["box"], GeometryBasics.Rect(Vec(-1.0 * r,
		-1.0 * r,
		-1.0 * r),
		Vec(2.0 * r, 2.0 * r, 2.0 * r)),
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    for i = 1:model.n_corners
        setobject!(vis["corner$i"], GeometryBasics.Sphere(Point3f0(0),
            convert(Float32, 0.05)),
            MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))
    end

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["box"],
				compose(Translation(q[t][1:3]...), LinearMap(UnitQuaternion(q[t][4:7]...))))

            for i = 1:model.n_corners
                settransform!(vis["corner$i"],
                    Translation((q[t][1:3] + UnitQuaternion(q[t][4:7]...) * (corner_offset[i]))...))
            end
        end
    end
    # settransform!(vis["/Cameras/default"], compose(Translation(-1, -1, 0),LinearMap(RotZ(pi/2))))
    MeshCat.setanimation!(vis, anim)
end
