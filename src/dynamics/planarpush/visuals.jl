
function plot_lines!(vis::Visualizer, model::PlanarPush, q::AbstractVector;
		offset=0.15, size=10, name::Symbol=:PlanarPush, col::Bool=true)
	p_shift = [0.0, 0.0, offset]

	# Point Traj
	object_point = Vector{Point{3,Float64}}()
	pusher_point = Vector{Point{3,Float64}}()

	for qi in q
		push!(object_point, Point(qi[1:3]...) + p_shift)
		push!(pusher_point, Point(qi[5:7]...) + p_shift)
	end

	# Set lines
	orange_mat, blue_mat, black_mat = get_line_material(size)
	if col
		setobject!(vis[name]["/lines/object"], MeshCat.Line(object_point, orange_mat))
		setobject!(vis[name]["/lines/pusher"], MeshCat.Line(pusher_point, blue_mat))
	else
		setobject!(vis[name]["/lines/object"], MeshCat.Line(object_point, black_mat))
		setobject!(vis[name]["/lines/pusher"], MeshCat.Line(pusher_point, black_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::PlanarPush; name::Symbol=:PlanarPush, height=0.20, α=1.0)
	height = convert(Float32, height)
	r = convert(Float32, model.r)
	rp = convert(Float32, model.rp)
	rc = convert(Float32, model.r/sqrt(2))
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(51/255, 1.0, 1.0, α))

	default_background!(vis)

	Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, height), r)
	obj_p = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, height), rp)
	setobject!(vis[name][:robot]["object"]["cyl"], Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, height/5), r), body_mat)
	setobject!(vis[name][:robot]["object"]["point1"], Cylinder(Point3f0(rc/2, 0, 0), Point3f0(rc/2, 0.0, height/4), r/10), contact_mat)
	setobject!(vis[name][:robot]["object"]["point2"], Cylinder(Point3f0(-rc/2, 0, 0), Point3f0(-rc/2, 0.0, height/4), r/10), contact_mat)
    setobject!(vis[name][:robot]["pusher"], Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, height), rp), contact_mat)
	return nothing
end

function set_robot!(vis::Visualizer, model::PlanarPush, q::AbstractVector;
		name::Symbol=:PlanarPush)
	r = convert(Float32, model.r)
	rp = convert(Float32, model.rp)

	p_object = q[1:3]
    p_pusher = q[5:7]

    settransform!(vis[name][:robot]["object"], compose(Translation(p_object), LinearMap(RotZ(q[4]))))
    settransform!(vis[name][:robot]["pusher"], Translation(p_pusher))
	return nothing
end

function contact_point(model::PlanarPush, q::AbstractVector)
	p_floor = kinematics_1(model, q, body = :floor)
	p_object = kinematics_1(model, q, body = :object)
	p_pusher = kinematics_1(model, q, body = :pusher)
	pc = [p_floor, p_object, p_pusher]
	return pc
end
