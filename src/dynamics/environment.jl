abstract type R2 end
abstract type R3 end

struct Environment{T}
    surf::Any
	surf_grad::Any
end

function environment_2D_flat()
	Environment{R2}(x -> 0.0, x-> zero(x))
end

function environment_3D_flat()
	Environment{R3}(x -> 0.0, x-> zero(x))
end

function environment_2D(surf)
	@variables q[1:1]

	s = surf(q)
	s = ModelingToolkit.simplify.(s)
	ds = ModelingToolkit.gradient(s, q, simplify = true)

	surf_fast = eval(ModelingToolkit.build_function(s, q))
	surf_grad_fast = eval(ModelingToolkit.build_function(ds, q)[1])

	Environment{R2}(surf_fast, surf_grad_fast)
end

function environment_3D(surf)
	@variables q[1:2]

	s = surf(q)
	s = ModelingToolkit.simplify.(s)
	ds = ModelingToolkit.gradient(s, q, simplify = true)

	surf_fast = eval(ModelingToolkit.build_function(s, q))
	surf_grad_fast = eval(ModelingToolkit.build_function(ds, q)[1])

	Environment{R3}(surf_fast, surf_grad_fast)
end

function skew(x)
	SMatrix{3,3}([0.0 -x[3] x[2];
	               x[3] 0.0 -x[1];
				   -x[2] x[1] 0.0])
end

# rotation matrix rotating unit vector a onto unit vector b
function rot(a, b)
	v = cross(a, b)
	s = sqrt(transpose(v) * v)
	c = transpose(a) * b

	R = Diagonal(@SVector ones(3)) + skew(v) + 1.0 / (1.0 + c) * skew(v) * skew(v)
end

function rotation(env::Environment{R3}, q)
	# unit surface normal (3D)
	n = [-1.0 * env.surf_grad(q[1:2]); 1.0]
	ns = n ./ sqrt(transpose(n) * n)

	# world-frame normal
	nw = @SVector [0.0, 0.0, 1.0]

	rot(ns, nw)
end

function rotation(env::Environment{R2}, q)
	# unit surface normal (3D)
	n = [-1.0 * env.surf_grad(q[1]); 1.0]
	ns = n ./ sqrt(transpose(n) * n)

	# world-frame normal
	nw = @SVector [0.0, 1.0]

	ang = atan(nw[2], nw[1]) - atan(ns[2], ns[1])

	SMatrix{2,2}([cos(ang) -sin(ang); sin(ang) cos(ang)])
end

function friction_mapping(env::Environment{R2})
    SMatrix{1, 2}([1.0 -1.0])
end

function friction_mapping(env::Environment{R3})
    SMatrix{2, 4}([1.0 0.0 -1.0 0.0;
                   0.0 1.0 0.0 -1.0])
end

dim(env::Environment{R2}) = 2
dim(env::Environment{R3}) = 3
