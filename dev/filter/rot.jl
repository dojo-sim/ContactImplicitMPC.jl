function visualize!(vis, model::Particle, q;
	Δt = 0.1, r = 0.25, name=:particle)

	default_background!(vis)

    setobject!(vis[name][:body],
		# GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		# convert(Float32, r)),
		MeshCat.HyperRectangle(MeshCat.Vec(0,0,0),MeshCat.Vec(r,r,r)),
		MeshPhongMaterial(color = RGBA(165.0 / 255.0, 0, 1, 1.0)))

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis[name][:body], MeshCat.Translation(q[t][1:3]...))
        end
    end

	# settransform!(vis["/Cameras/default"],
	#     compose(Translation(-2.5, 7.5, 1.0),LinearMap(RotZ(0.0))))

    MeshCat.setanimation!(vis, anim)
end

function visualize_uncertainty!(vis, model::Particle, q, P;
	Δt = 0.1, r = 0.25, name=:particle)

	H = length(q)
	ve = [covariance_processing(P[t])[1] for t=1:H]
	va = [covariance_processing(P[t])[2] for t=1:H]
	default_background!(vis)

    setobject!(vis[name][:body],
		MeshCat.HyperRectangle(MeshCat.Vec(0,0,0),MeshCat.Vec(r,r,r)),
		MeshPhongMaterial(color = RGBA(0.0, 0.9, 0.7, 1.0)))
	for t = 1:H
		setobject!(vis[name][:ellip][Symbol(t)],
			HyperEllipsoid(Point(0,0,0.), Vec(5.0*va[t]...)),
			MeshPhongMaterial(color = RGBA(1.0, 0.0, 0.5, 0.3)))
	end
    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:H
        MeshCat.atframe(anim, t) do
			for l in setdiff(1:H,t)
				setvisible!(vis[name][:ellip][Symbol(l)], false)
			end
			setvisible!(vis[name][:ellip][Symbol(t)], true)

			settransform!(vis[name][:body], MeshCat.Translation(q[t][1:3]...))
			for l = 1:H
					settransform!(
					vis[name][:ellip][Symbol(l)],
					compose(
						Translation(r/2*ones(3)+q[t][1:3]...),
						LinearMap(ve[l]),
					)
				)
			end
        end
    end

	# settransform!(vis["/Cameras/default"],
	#     compose(Translation(-2.5, 7.5, 1.0),LinearMap(RotZ(0.0))))

    MeshCat.setanimation!(vis, anim)
end


function vis_3d_gaussian2!(vis::Visualizer, μ::AbstractVector, ve, va; r=0.25,
	name=:ellip, ind::Int=0)
	default_background!(vis)

    setobject!(vis[name][Symbol(ind)],
		HyperEllipsoid(Point(0,0,0.), Vec(1.0*va...)),
		MeshPhongMaterial(color = RGBA(165.0 / 255.0, 0, 1, 0.4)))
	settransform!(
		vis[name][Symbol(ind)],
		compose(
			Translation(r/2*ones(3)+μ...),
			LinearMap(ve),
		)
	)
	return nothing
end

function covariance_processing(P)
	P_ = sqrt(P)
	ve = real.(eigvecs(P_))
	va = real.(eigvals(P_))
	if real.(det(ve)) < 0.0
		ve[:,1] = -ve[:,1]
	end
	return ve, va
end

# vis = Visualizer()
# open(vis)
#
# function vvv(vis, )
# 	anim = MeshCat.Animation(convert(Int, floor(1.0 / 0.1)))
#
# 	for t = 1:H
# 		MeshCat.atframe(anim, t) do
# 			# setobject!(vis[:test],
# 			# 	HyperEllipsoid(Point(0,0,0.), Vec(0.1, 0.2, 0.4)),
# 			# 	MeshPhongMaterial(color = RGBA(165.0 / 255.0, 0, 1, 0.4)))
#
#
# 			# settransform!(vis[name][:body], MeshCat.Translation(q[t][1:3]...))
# 			# ve, va = covariance_processing(P[t])
# 			# # settransform!(vis[name][:ellip][Symbol(t)], MeshCat.Translation(q[t][1:3]...))
# 			# settransform!(vis[name][:ellip], MeshCat.Translation(q[t][1:3]...))
# 		end
# 	end
#
# 	return nothing
# end
#
# vvv(vis)
# vis.core
