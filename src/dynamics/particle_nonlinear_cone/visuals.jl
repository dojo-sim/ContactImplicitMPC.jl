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
		MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["particle"], MeshCat.Translation(q[t][1:3]...))
        end
    end

    MeshCat.setanimation!(vis, anim)
end
