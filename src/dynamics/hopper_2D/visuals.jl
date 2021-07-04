
function plot_lines!(vis::Visualizer, model::Hopper2D, q::AbstractVector;
		r_foot=0.04, offset=0.04, size=10, name::Symbol=:hopper_2D, col::Bool=true)
	p_shift = [0.0, 0.0, r_foot]

	# Point Traj
	top_point = Vector{Point{3,Float64}}()
	bot_point = Vector{Point{3,Float64}}()

	for qi in q
		push!(top_point, Point(qi[1], -offset, qi[2])+p_shift)
		push!(bot_point, Point(kinematics(model, qi)[1], -offset, kinematics(model, qi)[2])+p_shift)
	end

	# Set lines
	orange_mat, blue_mat, black_mat = get_line_material(size)
	if col
		setobject!(vis[name]["/lines/top"], MeshCat.Line(top_point, orange_mat))
		setobject!(vis[name]["/lines/bot"], MeshCat.Line(bot_point, blue_mat))
	else
		setobject!(vis[name]["/lines/top"], MeshCat.Line(top_point, black_mat))
		setobject!(vis[name]["/lines/bot"], MeshCat.Line(bot_point, black_mat))
	end
	return nothing
end

function build_robot!(vis::Visualizer, model::Hopper2D; name::Symbol=:Hopper2D, r=0.04, α=1.0)
	r = convert(Float32, r)
	r_leg = convert(Float32, 0.5 * r)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))

	default_background!(vis)

    setobject!(vis[name][:robot]["body"], Sphere(Point3f0(0),
        convert(Float32, 0.1)), body_mat)
    setobject!(vis[name][:robot]["foot"], Sphere(Point3f0(0), r), contact_mat)
    setobject!(vis[name][:robot]["leg"], Cylinder(Point3f0(0),
		Point3f0(0, 0, 1.0), r_leg), body_mat)
	return nothing
end

function set_robot!(vis::Visualizer, model::Hopper2D, q::AbstractVector;
		name::Symbol=:Hopper2D, r=0.04)
	r = convert(Float32, r)
	r_leg = convert(Float32, 0.5 * r)
	p_shift = [0.0; 0.0; r]

    p_body = [q[1], 0.0, q[2]]
    p_foot = [kinematics(model, q)[1], 0.0, kinematics(model, q)[2]]
	leg_length = norm(p_body - p_foot)

    settransform!(vis[name][:robot]["body"], Translation(p_body + p_shift))
    settransform!(vis[name][:robot]["foot"], Translation(p_foot + p_shift))
	settransform!(vis[name][:robot]["leg"],
		compose(
			cable_transform(p_foot + p_shift, p_body + p_shift),
			LinearMap(Diagonal([1,1,leg_length])))
		)
	return nothing
end

function contact_point(model::Hopper2D, q::AbstractVector)
	p_foot = [kinematics(model, q)[1], 0.0, kinematics(model, q)[2]]
	pc = [p_foot]
	return pc
end

function stairs!(vis)
	setobject!(vis["box1"], GeometryBasics.HyperRectangle(Vec(0.0, 0.0, 0.0),
		Vec(0.25, 0.5, 0.25)), MeshPhongMaterial(color = RGBA(0.5, 0.5, 0.5, 1.0)))
	settransform!(vis["box1"], Translation(0.125, -0.25, 0))

	setobject!(vis["box2"], GeometryBasics.HyperRectangle(Vec(0.0, 0.0, 0.0),
		Vec(0.25, 0.5, 2 * 0.25)), MeshPhongMaterial(color = RGBA(0.5, 0.5, 0.5, 1.0)))
	settransform!(vis["box2"], Translation(0.125 + 0.25, -0.25, 0))

	setobject!(vis["box3"], GeometryBasics.HyperRectangle(Vec(0.0, 0.0, 0.0),
		Vec(0.25, 0.5, 3 * 0.25)), MeshPhongMaterial(color = RGBA(0.5, 0.5, 0.5, 1.0)))
	settransform!(vis["box3"], Translation(0.125 + 2 * 0.25, -0.25, 0))
end

function hopper_parkour_ghost!(vis, sim, traj, ref_traj;
    idx = [1, 35, 50, 110, 130, 190, 210, 265, 270, 284, 295, 300, 305, 320],
    α = [convert(Float64, i / length(idx)) for i = 1:length(idx)],
    _name = "_stairs")

    plot_lines!(vis, sim.s.model, ref_traj.q[1:1:end], size = 5, offset = -0.5)
    stairs!(vis)
    # settransform!(vis["/Cameras/default"],
    #         compose(Translation(0.0, -95.0, -1.0), LinearMap(RotY(0.0 * π) * RotZ(-π / 2.0))))
    # setprop!(vis["/Cameras/default/rotated/<object>"], "zoom", 20)

    for (i, t) in enumerate(idx)
        name = Symbol("Hopper" * "$i" * _name)
        build_robot!(vis, sim.s.model, name=name, α = α[i])
        set_robot!(vis, sim.s.model, traj.q[t], name = name)
    end
end
