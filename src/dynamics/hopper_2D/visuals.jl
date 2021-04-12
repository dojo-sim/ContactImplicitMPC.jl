# # Visualization
# function visualize!(vis, model::Hopper2D, q;
# 		Δt = 0.1, scenario = :vertical, name::Symbol=:hopper_2D, r_foot = 0.04, camera::Bool=false)
#
# 	default_background!(vis)
#
#     r_leg = 0.5 * r_foot
#
#     setobject!(vis[name]["body"], Sphere(Point3f0(0),
#         convert(Float32, 0.1)),
#         MeshPhongMaterial(color = RGBA(0, 1, 0, 1.0)))
#
#     setobject!(vis[name]["foot"], Sphere(Point3f0(0),
#         convert(Float32, r_foot)),
#         MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))
#
#     n_leg = 100
#     for i = 1:n_leg
#         setobject!(vis[name]["leg$i"], Sphere(Point3f0(0),
#             convert(Float32, r_leg)),
#             MeshPhongMaterial(color = RGBA(0, 0, 0, 1.0)))
#     end
#
#     p_leg = [zeros(3) for i = 1:n_leg]
#     anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))
#
#     for t = 1:length(q)
#         p_body = [q[t][1], 0.0, q[t][2]]
#         p_foot = [kinematics(model, q[t])[1], 0.0, kinematics(model, q[t])[2]]
#
#         q_tmp = Array(copy(q[t]))
#         r_range = range(0, stop = q[t][4], length = n_leg)
#         for i = 1:n_leg
#             q_tmp[4] = r_range[i]
#             p_leg[i] = [kinematics(model, q_tmp)[1], 0.0, kinematics(model, q_tmp)[2]]
#         end
#         q_tmp[4] = q[t][4]
#         p_foot = [kinematics(model, q_tmp)[1], 0.0, kinematics(model, q_tmp)[2]]
#
#         z_shift = [0.0; 0.0; r_foot]
#
#         MeshCat.atframe(anim, t) do
#             settransform!(vis[name]["body"], Translation(p_body + z_shift))
#             settransform!(vis[name]["foot"], Translation(p_foot + z_shift))
#
#             for i = 1:n_leg
#                 settransform!(vis[name]["leg$i"], Translation(p_leg[i] + z_shift))
#             end
#         end
#     end
# 	if camera
# 		settransform!(vis["/Cameras/default"],
# 			compose(Translation(0.0, 0.5, -1.0),LinearMap(RotZ(-pi / 2.0))))
# 	end
#
#     MeshCat.setanimation!(vis, anim)
# end


function plot_lines!(vis::Visualizer, model::Hopper2D, q::AbstractVector;
		r_foot=0.04, offset=0.04, size=10, name::Symbol=:hopper_2D, col::Bool=true)
	p_shift = [0.0, 0.0, r_foot]

	kinematics(model::Hopper2D, q) = [q[1] + q[4] * sin(q[3]), q[2] - q[4] * cos(q[3])]

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

function build_robot!(vis::Visualizer, model::Hopper2D; name::Symbol=:Hopper2D, r=0.04)
	n_leg = 100
	r = convert(Float32, r)
	r_leg = convert(Float32, 0.5 * r)
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, 1.0))

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
