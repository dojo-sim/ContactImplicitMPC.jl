function plot_lines!(vis::Visualizer, model::Box2D, q::AbstractVector;
		r_corner=0.04, size=10, name::Symbol=:Box2D, col::Bool=true)

	nc = model.dim.c
	# Point Traj
	com_point = Vector{Point{3,Float64}}()
	corner_points = [Vector{Point{3,Float64}}() for i=1:nc]

	for qi in q
		push!(com_point, Point(qi[1], 0.0, qi[2]))
		k = kinematics(model, qi)
		for i = 1:nc
			push!(corner_points[i], Point(k[(i-1)*2 + 1], 0.0, k[(i-1)*2 + 2]))
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

function build_robot!(vis::Visualizer, model::Box2D; name::Symbol=:Box2D, r_contact=0.04, α=1.0)
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

function set_robot!(vis::Visualizer, model::Box2D, q::AbstractVector; name::Symbol=:Box2D)
	nc = model.dim.c

	settransform!(vis[name][:robot]["box"],
		compose(Translation(q[1], 0, q[2]), LinearMap(RotY(q[3]))))
	k = kinematics(model, q)
	for i = 1:nc
		pi = [k[(i-1)*2 + 1], 0.0, k[(i-1)*2 + 2]]
		settransform!(vis[name][:robot]["corners"]["corner$i"], Translation(pi...))
	end
	return nothing
end

function contact_point(model::Box2D, q::AbstractVector)
	nc = model.dim.c
	k = kinematics(model, q)
	pc = [k[(i-1)*3 .+ (1:3)] for i=1:nc]
	return pc
end
