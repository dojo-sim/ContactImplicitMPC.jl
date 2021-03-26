# Visualization
function visualize!(vis, model::Hopper2D, q;
		Δt = 0.1, scenario = :vertical, name::Symbol=:hopper_2D)

	default_background!(vis)

    r_foot = 0.05
    r_leg = 0.5 * r_foot

    setobject!(vis[name]["body"], Sphere(Point3f0(0),
        convert(Float32, 0.1)),
        MeshPhongMaterial(color = RGBA(0, 1, 0, 1.0)))

    setobject!(vis[name]["foot"], Sphere(Point3f0(0),
        convert(Float32, r_foot)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

    n_leg = 100
    for i = 1:n_leg
        setobject!(vis[name]["leg$i"], Sphere(Point3f0(0),
            convert(Float32, r_leg)),
            MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))
    end

    p_leg = [zeros(3) for i = 1:n_leg]
    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        p_body = [q[t][1], 0.0, q[t][2]]
        p_foot = [kinematics(model, q[t])[1], 0.0, kinematics(model, q[t])[2]]

        q_tmp = Array(copy(q[t]))
        r_range = range(0, stop = q[t][4], length = n_leg)
        for i = 1:n_leg
            q_tmp[4] = r_range[i]
            p_leg[i] = [kinematics(model, q_tmp)[1], 0.0, kinematics(model, q_tmp)[2]]
        end
        q_tmp[4] = q[t][4]
        p_foot = [kinematics(model, q_tmp)[1], 0.0, kinematics(model, q_tmp)[2]]

        z_shift = [0.0; 0.0; r_foot]

        MeshCat.atframe(anim, t) do
            settransform!(vis[name]["body"], Translation(p_body + z_shift))
            settransform!(vis[name]["foot"], Translation(p_foot + z_shift))

            for i = 1:n_leg
                settransform!(vis[name]["leg$i"], Translation(p_leg[i] + z_shift))
            end
        end
    end

	settransform!(vis["/Cameras/default"],
		compose(Translation(0.0, 0.5, -1.0),LinearMap(RotZ(-pi / 2.0))))

    MeshCat.setanimation!(vis, anim)
end
