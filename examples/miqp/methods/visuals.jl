function build_robot!(vis::Visualizer, model::WallPendulum; name::Symbol=:WallPendulum,
		r = 0.01, l = 1.0, n_leg::Int = 30, α = 1.0)

	l = model.l
	d = model.d

	r = convert(Float32, r)
	r_contact = convert(Float32, r * 2)
	# r_arm = convert(Float32, r * 0.3)
	wall_thickness = 0.15
	# l_joint = 0.16
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))
	wall_mat = MeshPhongMaterial(color = RGBA(0.7, 0.7, 0.7, 1.0))

	default_background!(vis)

	link = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, l), r)
	setobject!(vis[name][:robot]["link"], link, body_mat)
	setobject!(vis[name][:robot]["endeffector"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	setobject!(vis[name][:robot]["linkbase"], Sphere(Point3f0(0.0), r), body_mat)

	# setobject!(vis[name][:robot]["contact1"], Sphere(Point3f0(0.0), r_contact), contact_mat)
	#
	# for i = 1:n_leg
	# 	setobject!(vis[name][:robot]["arm$i"], Sphere(Point3f0(0), r_arm), body_mat)
	# end

	setobject!(vis[name][:env]["fixed_wall1"],
		Rect(Vec(0, 0, 0),Vec(wall_thickness, wall_thickness, 1.5)), wall_mat)
	setobject!(vis[name][:env]["fixed_wall2"],
		Rect(Vec(0, 0, 0),Vec(wall_thickness, wall_thickness, 1.5)), wall_mat)

	setobject!(vis[name][:env][:spring_bundle1]["spring_wall1"],
		Rect(Vec(0, 0, 0),Vec(wall_thickness/10, wall_thickness, 0.25)), wall_mat)
	setobject!(vis[name][:env][:spring_bundle2]["spring_wall2"],
		Rect(Vec(0, 0, 0),Vec(wall_thickness/10, wall_thickness, 0.25)), wall_mat)

	spring = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, wall_thickness), r)
	setobject!(vis[name][:env][:spring_bundle1]["spring1"], spring, body_mat)
	setobject!(vis[name][:env][:spring_bundle2]["spring2"], spring, body_mat)

	p_spring1 = [                     0*l*sin(d) + r_contact; -wall_thickness/2; l-0.125]
	p_spring2 = [-wall_thickness/10 - 0*l*sin(d) - r_contact; -wall_thickness/2; l-0.125]
	settransform!(vis[name][:env]["fixed_wall1"], Translation([                  0.2 + 1.5 * r; -wall_thickness/2; 0.0]))
	settransform!(vis[name][:env]["fixed_wall2"], Translation([-wall_thickness - 0.2 - 1.5 * r; -wall_thickness/2; 0.0]))
	settransform!(vis[name][:env][:spring_bundle1]["spring_wall1"], Translation(p_spring1))
	settransform!(vis[name][:env][:spring_bundle2]["spring_wall2"], Translation(p_spring2))
	settransform!(vis[name][:env][:spring_bundle1]["spring1"], compose(Translation(p_spring1 + [wall_thickness/10, wall_thickness/2, 0.125]), LinearMap(RotY(π/2))))
	settransform!(vis[name][:env][:spring_bundle2]["spring2"], compose(Translation(p_spring2 + [0,                 wall_thickness/2, 0.125]), LinearMap(RotY(-π/2))))

	return nothing
end


function set_robot!(vis::Visualizer, model::WallPendulum, x::AbstractVector;
		name::Symbol=:WallPendulum, r = 0.01, n_leg::Int=30)

	l = model.l
	d = model.d

	r = convert(Float32, r)
	r_contact = convert(Float32, r * 2)

	p = [kinematics(x)[1], 0.0, kinematics(x)[2]]

	# settransform!(vis[name][:robot]["joint1"], compose(cable_transform(zeros(3), pee), Translation(0,0,0)))
	# settransform!(vis[name][:robot]["joint2"], cable_transform(zeros(3), pee))
	settransform!(vis[name][:robot]["link"], cable_transform(p, zeros(3)))
	settransform!(vis[name][:robot]["endeffector"], Translation(p))
	if p[1] < -l * sin(d) # left
		settransform!(vis[name][:env][:spring_bundle2], Translation(p[1], 0, 0))
	else
		settransform!(vis[name][:env][:spring_bundle2], Translation(-l*sin(d), 0, 0))
	end
	if p[1] > l * sin(d) # right
		# settransform!(vis[name][:env][:spring_bundle1], Translation(p[1], 0, 0))
		settransform!(vis[name][:env][:spring_bundle1], Translation(l*sin(d), 0, 0))
	else
		settransform!(vis[name][:env][:spring_bundle1], Translation(l*sin(d), 0, 0))
	end

	# settransform!(vis[name][:robot]["contact1"], Translation(p1))

	# r_range = range(0, stop = q[2], length = n_leg)
	# for i = 1:n_leg
	# 	q_tmp = Array(copy(q))
	# 	q_tmp[2] = r_range[i]
	# 	p_leg = [_kinematics(model, q_tmp, mode = :d)[1], 0.0, _kinematics(model, q_tmp, mode = :d)[2]]
    #     settransform!(vis[name][:robot]["arm$i"], Translation(p_leg))
    # end

	return nothing
end

function pusher_arc(center::AbstractVector, radius::Real, α::Real; sense::Symbol=:trig)
	p = center + radius*[cos(α), 0, sin(α)]
	if sense == :trig
		θ = α + π/2
	elseif sense == :antitrig
		θ = α - π/2
	else
		error("Invalid sense.")
	end
	pθ = [p; θ]
	return pθ
end

function pusher_arc_traj(center::AbstractVector, radius::Real, α_start::Real, α_end::Real,
	 	N::Int; sense::Symbol=:trig)
	pθ = []
	for α in range(α_start, stop=α_end, length=N)
		pθ_ = pusher_arc(center, radius, α, sense=sense)
		push!(pθ, pθ_)
	end
	return pθ
end

function generate_pusher_traj(ind::Vector{Int}, w::Vector, x::Vector; side::Symbol=:right,
		center=zeros(3), radius::Real=0.8, N_arc::Int=17)
	H = length(x)
	if side == :right
		sense = :trig
		θ_rest = 0.0
	elseif side == :left
		sense = :antitrig
		θ_rest = π
	else
		error("Invalid side.")
	end

	pθ = [pusher_arc(center, radius, θ_rest; sense=sense) for t=1:H]

	for (i, idx) in enumerate(ind)
		wi = w[i][1]
		θ_impact = x[idx+1][1] + π/2
		θ_next = x[idx+2][1] + π/2

		if ((wi > 0.0) && (side == :right)) || ((wi < 0.0) && (side == :left))
			forward_arc = pusher_arc_traj(center, radius, θ_rest, θ_impact, N_arc; sense=sense)
			backward_arc = pusher_arc_traj(center, radius, θ_next, θ_rest, N_arc; sense=sense)
			arc = [forward_arc; backward_arc] # 2N
			pθ[idx-N_arc-0:idx+N_arc-1] .= arc
		end
	end
	return pθ
end

function build_disturbance!(vis::Visualizer, model::WallPendulum; name::Symbol=:Pusher, r=0.01, α=1.0)
	r1 = convert(Float32, r)
	r2 = convert(Float32, 6r)
	pusher_mat = MeshPhongMaterial(color = RGBA(1.0, 0.0, 0.0, α))

	# build pusher
	p1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, 0.150), r1)
	p2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, 0.03), r2)
	setobject!(vis[name, "pusher1"], p1, pusher_mat)
	setobject!(vis[name, "pusher2"], p2, pusher_mat)
	return anim
end

function set_disturbance!(vis::Visualizer, model::WallPendulum, pθ::AbstractVector; name::Symbol=:Pusher, offset::Real=0.0)
	p = pθ[1:3]
	θ = pθ[4]
	settransform!(vis[name],
		compose(Translation(p...),
		compose(LinearMap(RotY(-θ-pi/2)),
				Translation(0,0,offset))))
	return nothing
end

function animate_disturbance!(vis::Visualizer, anim::MeshCat.Animation, model::WallPendulum,
		pθ::AbstractVector; name::Symbol=:Pusher, offset::Real=0.0)
	H = length(pθ)
	for t in 1:H
		MeshCat.atframe(anim, t) do
			set_disturbance!(vis, model, pθ[t], name=name, offset=offset)
		end
	end
	setanimation!(vis, anim)
	return anim
end

function visualize_disturbance!(vis::Visualizer, model::ContactModel, pθ::AbstractVector;
		h=0.01, α=1.0,
		sample=max(1, Int(floor(length(pθ) / 100))),
		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
		name::Symbol=:Pusher, offset::Real=0.0)

	build_disturbance!(vis, model, name=name, α=α)
	animate_disturbance!(vis, anim, model, pθ[1:sample:end], name=name, offset=offset)
	return anim
end
