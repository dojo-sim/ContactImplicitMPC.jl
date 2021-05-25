# visuals
function visualize!(vis, p, q; Δt = 0.1)
	setvisible!(vis["/Background"], true)
	setprop!(vis["/Background"], "top_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setprop!(vis["/Background"], "bottom_color", RGBA(1.0, 1.0, 1.0, 1.0))
	setvisible!(vis["/Axes"], false)
	setvisible!(vis["/Grid"], false)

    setobject!(vis[:satellite],
    	Rect(Vec(-0.25, -0.25, -0.25),Vec(0.5, 0.5, 0.5)),
    	MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    arrow_x = ArrowVisualizer(vis[:satellite][:arrow_x])
    mat = MeshPhongMaterial(color=RGBA(1.0, 0.0, 0.0, 1.0))
    setobject!(arrow_x, mat)
    settransform!(arrow_x,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.75, 0.0, 0.0),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    arrow_y = ArrowVisualizer(vis[:satellite][:arrow_y])
    mat = MeshPhongMaterial(color=RGBA(0.0, 1.0, 0.0, 1.0))
    setobject!(arrow_y, mat)
    settransform!(arrow_y,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.0, 0.75, 0.0),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    arrow_z = ArrowVisualizer(vis[:satellite][:arrow_z])
    mat = MeshPhongMaterial(color=RGBA(0.0, 0.0, 1.0, 1.0))
    setobject!(arrow_z, mat)
    settransform!(arrow_z,
    	Point(0.0, 0.0, 0.0),
    	Vec(0.0, 0.0, 0.75),
    	shaft_radius=0.05,
    	max_head_radius=0.1)

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

     # for t = 1:length(q)
    	#  MeshCat.atframe(anim, t) do
    	# 	 settransform!(vis["satellite"],
    	# 		   compose(Translation((q[t][1:3] + [-0.25; -0.25; -0.25])...),
    	# 				 LinearMap(UnitQuaternion(q[t][4:7]...))))
    	#  end
     # end

	 for t = 1:length(q)
		MeshCat.atframe(anim, t) do
			settransform!(vis["satellite"],
				  compose(Translation(([-0.0; -0.0; -0.0])...),
						LinearMap(UnitQuaternion(q[t][1:4]...))))
		end
	 end
	#
	#
    MeshCat.setanimation!(vis, anim)
end
