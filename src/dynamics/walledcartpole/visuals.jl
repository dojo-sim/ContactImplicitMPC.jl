function build_robot!(vis::Visualizer, model::WalledCartpole13; name::Symbol=:WalledCartpole13,
		r=0.05, α=1.0)
	r = convert(Float32, r)
	r_contact = convert(Float32, r * 1.5)

	wall_thickness = 0.15
	body_mat = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, α))
	contact_mat = MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0.0, α))
	wall_mat = MeshPhongMaterial(color = RGBA(0.7, 0.7, 0.7, 1.0))

	default_background!(vis)

	link = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, model.l), r/2)
	setobject!(vis[name][:robot]["link"], link, body_mat)
	vec_base = Vec(2r_contact, 4r, 1.5r)
	vec_shaft = Vec(2model.w + 2r_contact + 2wall_thickness, 0.5r, 0.5r)
	setobject!(vis[name][:robot]["linkbase"], Rect(-0.5 * vec_base, vec_base), body_mat)
	setobject!(vis[name][:robot]["shaft"], Rect(-0.5 * vec_shaft, vec_shaft), body_mat)
	setobject!(vis[name][:robot]["contact"], Sphere(Point3f0(0.0), r_contact), contact_mat)

	setobject!(vis[name][:env]["wall1"],
		Rect(Vec(0, 0, 0),Vec(wall_thickness, wall_thickness, 0.75)), wall_mat)
	setobject!(vis[name][:env]["wall2"],
		Rect(Vec(0, 0, 0),Vec(wall_thickness, wall_thickness, 0.75)), wall_mat)
	settransform!(vis[name][:env]["wall1"], Translation([model.w + 1.5 * r;                   -wall_thickness/2; -0.25r]))
	settransform!(vis[name][:env]["wall2"], Translation([-wall_thickness - model.w - 1.5 * r; -wall_thickness/2; -0.25r]))
	return nothing
end


function set_robot!(vis::Visualizer, model::WalledCartpole13, q::AbstractVector;
		name::Symbol=:WalledCartpole13)

	pt = [_kinematics(model, q, mode =  :tip)[1], 0.0, _kinematics(model, q, mode =  :tip)[2]]
    pb = [_kinematics(model, q, mode = :base)[1], 0.0, _kinematics(model, q, mode = :base)[2]]

	settransform!(vis[name][:robot]["link"], cable_transform(pt, pb))
	settransform!(vis[name][:robot]["linkbase"], Translation(pb))
	settransform!(vis[name][:robot]["contact"], Translation(pt))

	xw1 = q[3]
	xw2 = q[4]
	settransform!(vis[name][:env]["wall1"], Translation([xw1, 0.0, 0.0]))
	settransform!(vis[name][:env]["wall2"], Translation([xw2, 0.0, 0.0]))

	return nothing
end


# function pusher_arc(center::AbstractVector, radius::Real, α::Real; sense::Symbol=:trig)
# 	p = center + radius*[cos(α), 0, sin(α)]
# 	if sense == :trig
# 		θ = α + π/2
# 	elseif sense == :antitrig
# 		θ = α - π/2
# 	else
# 		error("Invalid sense.")
# 	end
# 	pθ = [p; θ]
# 	return pθ
# end
#
# function pusher_arc_traj(center::AbstractVector, radius::Real, α_start::Real, α_end::Real,
# 	 	N::Int; sense::Symbol=:trig)
# 	pθ = []
# 	for α in range(α_start, stop=α_end, length=N)
# 		pθ_ = pusher_arc(center, radius, α, sense=sense)
# 		push!(pθ, pθ_)
# 	end
# 	return pθ
# end
#
# function generate_pusher_traj(d::ImpulseDisturbance, traj::ContactTraj; side::Symbol=:right,
# 		center=zeros(3), radius::Real=0.8, N_arc::Int=17)
# 	H = traj.H
# 	if side == :right
# 		sense = :trig
# 		θ_rest = 0.0
# 	elseif side == :left
# 		sense = :antitrig
# 		θ_rest = π
# 	else
# 		error("Invalid side.")
# 	end
#
# 	pθ = [pusher_arc(center, radius, θ_rest; sense=sense) for t=1:H_sim]
#
# 	for (i,idx) in enumerate(d.idx)
# 		w = d.w[i][1]
# 		θ_impact = traj.q[idx+1][1] + π/2
# 		θ_next = traj.q[idx+2][1] + π/2
#
# 		if ((w > 0.0) && (side == :right)) || ((w < 0.0) && (side == :left))
# 			forward_arc = pusher_arc_traj(center, radius, θ_rest, θ_impact, N_arc; sense=sense)
# 			backward_arc = pusher_arc_traj(center, radius, θ_next, θ_rest, N_arc; sense=sense)
# 			arc = [forward_arc; backward_arc] # 2N
# 			pθ[idx-N_arc-0:idx+N_arc-1] .= arc
# 		end
# 	end
# 	return pθ
# end
#
# function build_disturbance!(vis::Visualizer, model::WalledCartpole13; name::Symbol=:Pusher, r=0.025, α=1.0)
# 	r1 = convert(Float32, r)
# 	r2 = convert(Float32, 4r)
# 	pusher_mat = MeshPhongMaterial(color = RGBA(1.0, 0.0, 0.0, α))
#
# 	# build pusher
# 	p1 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, 0.250), r1)
# 	p2 = Cylinder(Point3f0(0.0), Point3f0(0.0, 0.0, 0.075), r2)
# 	setobject!(vis[name, "pusher1"], p1, pusher_mat)
# 	setobject!(vis[name, "pusher2"], p2, pusher_mat)
# 	return anim
# end
#
# function set_disturbance!(vis::Visualizer, model::WalledCartpole13, pθ::AbstractVector; name::Symbol=:Pusher, offset::Real=0.0)
# 	p = pθ[1:3]
# 	θ = pθ[4]
# 	settransform!(vis[name],
# 		compose(Translation(p...),
# 		compose(LinearMap(RotY(-θ-pi/2)),
# 				Translation(0,0,offset))))
# 	return nothing
# end
#
# function animate_disturbance!(vis::Visualizer, anim::MeshCat.Animation, model::WalledCartpole13,
# 		pθ::AbstractVector; name::Symbol=:Pusher, offset::Real=0.0)
# 	H = length(pθ)
# 	for t in 1:H
# 		MeshCat.atframe(anim, t) do
# 			set_disturbance!(vis, model, pθ[t], name=name, offset=offset)
# 		end
# 	end
# 	setanimation!(vis, anim)
# 	return anim
# end
#
# function visualize_disturbance!(vis::Visualizer, model::ContactModel, pθ::AbstractVector;
# 		h=0.01, α=1.0,
# 		sample=max(1, Int(floor(length(pθ) / 100))),
# 		anim::MeshCat.Animation=MeshCat.Animation(Int(floor(1/h))),
# 		name::Symbol=:Pusher, offset::Real=0.0)
#
# 	build_disturbance!(vis, model, name=name, α=α)
# 	animate_disturbance!(vis, anim, model, pθ[1:sample:end], name=name, offset=offset)
# 	return anim
# end
