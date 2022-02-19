function visualize!(vis, model::BoxMRP, q;
        Δt = 0.1)

	default_background!(vis)

	r = abs(model.corner_offset[1][1])
    setobject!(vis["box"], GeometryBasics.Rect(Vec(-1.0 * r,
		-1.0 * r,
		-1.0 * r),
		Vec(2.0 * r, 2.0 * r, 2.0 * r)),
		MeshPhongMaterial(color = RGBA(13/255, 152/255, 186/255, 1.0)))

    for i = 1:model.n_corners
        setobject!(vis["corner$i"], GeometryBasics.Sphere(Point3f0(0),
            convert(Float32, 0.05)),
            MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0)))
    end

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["box"],
				compose(Translation(q[t][1:3]...), LinearMap(MRP(q[t][4:6]...))))

            for i = 1:model.n_corners
                settransform!(vis["corner$i"],
                    Translation((q[t][1:3] + MRP(q[t][4:6]...) * (corner_offset[i]))...))
            end
        end
    end
    # settransform!(vis["/Cameras/default"], compose(Translation(-1, -1, 0),LinearMap(RotZ(pi/2))))
    MeshCat.setanimation!(vis, anim)
end




function plot_lines!(vis::Visualizer, model::BoxMRP, q::AbstractVector;
		r_corner=0.04, size=10, name::Symbol=:BoxMRP, col::Bool=true)

	nc = model.dim.c
	# Point Traj
	com_point = Vector{Point{3,Float64}}()
	corner_points = [Vector{Point{3,Float64}}() for i=1:nc]

	for qi in q
		push!(com_point, Point(qi[1], qi[2], qi[3]))
		k = kinematics(model, qi)
		for i = 1:nc
			push!(corner_points[i], Point(k[(i-1)*3 .+ (1:3)]...))
		end
	end

	# Set lines
	orange_mat, blue_mat, black_mat = get_line_material(size)
	if col
		setobject!(vis[name]["/lines/com"], MeshCat.Line(com_point, orange_mat))
		for i = 1:nc
			setobject!(vis[name]["/lines/corners/corner$i"], MeshCat.Line(corner_points[i], blue_mat))
		end
	else
		setobject!(vis[name]["/lines/com"], MeshCat.Line(com_point, black_mat))
		for i = 1:nc
			setobject!(vis[name]["/lines/corners/corner$i"], MeshCat.Line(corner_points[i], black_mat))
		end
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::BoxMRP; name::Symbol=:BoxMRP, r_contact=0.04, α=1.0)
	nc = model.dim.c
	r_contact = convert(Float32, r_contact)
	body_mat = MeshPhongMaterial(color = RGBA(13/255, 152/255, 186/255, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))

	default_background!(vis)

	r = abs(model.corner_offset[1][1])
    setobject!(vis[name][:robot]["box"], GeometryBasics.Rect(
		Vec(-1.0 * r, -1.0 * r, -1.0 * r),
		Vec(2.0 * r, 2.0 * r, 2.0 * r)),
		body_mat)
	for i = 1:nc
		setobject!(vis[name][:robot]["corners"]["corner$i"], Sphere(Point3f0(0), r_contact), contact_mat)
	end
	return nothing
end

function set_robot!(vis::Visualizer, model::BoxMRP, q::AbstractVector; name::Symbol=:BoxMRP)
	nc = model.dim.c

	settransform!(vis[name][:robot]["box"],
		compose(Translation(q[1:3]...), LinearMap(MRP(q[4:6]...))))
	k = kinematics(model, q)
	for i = 1:nc
		pi = k[(i-1)*3 .+ (1:3)]
		settransform!(vis[name][:robot]["corners"]["corner$i"], Translation(pi...))
	end
	return nothing
end

function contact_point(model::BoxMRP, q::AbstractVector)
	nc = model.dim.c
	k = kinematics(model, q)
	pc = [k[(i-1)*3 .+ (1:3)] for i=1:nc]
	return pc
end
