function visual_element(b::Body{BoxGeometry})
    GeometryBasics.Rect(
        Vec(-1.0 * 0.5 * b.geometry[:length],
    		-1.0 * 0.5 * b.geometry[:width],
    		-1.0 * 0.5 * b.geometry[:height]),
	    Vec(2.0 * 0.5 * b.geometry[:length],
            2.0 * 0.5 * b.geometry[:width],
            2.0 * 0.5 * b.geometry[:height]))
end

function visual_element(b::Body{RodGeometry})
    GeometryBasics.Cylinder(
        Point3f0(0.0, 0.0, -0.5 * b.geometry[:length]),
		Point3f0(0.0, 0.0,  0.5 * b.geometry[:length]),
		convert(Float32, b.geometry[:radius]))
end

function visual_element(b::Body{SphereGeometry})
    GeometryBasics.Sphere(Point3f0(0), convert(Float32, b.geometry[:radius]))
end

function visualize!(vis, mechanism::Mechanism, q;
    colors = [RGBA(1.0, 0.0, 0.0, 1.0) for i = 1:length(mechanism.bodies)])

    # set default background
    default_background!(vis)


    # visual elements for each body in mechanism
    visual_elements = [visual_element(body) for body in mechanism.bodies]

    # add elements to the visualizer
    for (i, body) in enumerate(mechanism.bodies)
        setobject!(vis[body.id],
            visual_elements[i],
            MeshPhongMaterial(color = colors[i]))
    end

end
