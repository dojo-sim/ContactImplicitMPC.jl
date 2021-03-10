using Colors
using CoordinateTransformations
using FileIO
using GeometryBasics
using MeshCat, MeshIO, Meshing
using Rotations

function visualize!(vis, model::Particle, q;
	Δt = 0.1, r = 0.25)

	default_background!(vis)

    setobject!(vis["particle"],
		GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["particle"], Translation(q[t][1:3]...))
        end
    end

	# settransform!(vis["/Cameras/default"],
	#     compose(Translation(-2.5, 7.5, 1.0),LinearMap(RotZ(0.0))))

    MeshCat.setanimation!(vis, anim)
end
