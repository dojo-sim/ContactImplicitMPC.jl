
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
	n_leg = 100
	r = convert(Float32, r)
	r_leg = convert(Float32, 0.5 * r)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))

	default_background!(vis)

    setobject!(vis[name][:robot]["body"], Sphere(Point3f0(0),
        convert(Float32, 0.1)), body_mat)
    setobject!(vis[name][:robot]["foot"], Sphere(Point3f0(0), r), contact_mat)
    for i = 1:n_leg
        setobject!(vis[name][:robot]["leg$i"], Sphere(Point3f0(0), r_leg), body_mat)
    end
	return nothing
end

function set_robot!(vis::Visualizer, model::Hopper2D, q::AbstractVector;
		name::Symbol=:Hopper2D, r=0.04)
	n_leg = 100
	r = convert(Float32, r)
	r_leg = convert(Float32, 0.5 * r)
	p_shift = [0.0; 0.0; r]

    p_body = [q[1], 0.0, q[2]]
    p_foot = [kinematics(model, q)[1], 0.0, kinematics(model, q)[2]]

    settransform!(vis[name][:robot]["body"], Translation(p_body + p_shift))
    settransform!(vis[name][:robot]["foot"], Translation(p_foot + p_shift))
	r_range = range(0, stop = q[4], length = n_leg)
    for i = 1:n_leg
		q_tmp = Array(copy(q))
		q_tmp[4] = r_range[i]
		p_leg = [kinematics(model, q_tmp)[1], 0.0, kinematics(model, q_tmp)[2]]
        settransform!(vis[name][:robot]["leg$i"], Translation(p_leg + p_shift))
    end
	return nothing
end

function contact_point(model::Hopper2D, q::AbstractVector)
	p_foot = [kinematics(model, q)[1], 0.0, kinematics(model, q)[2]]
	pc = [p_foot]
	return pc
end
