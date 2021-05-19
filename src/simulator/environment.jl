# world
abstract type World end
abstract type R2 <: World end
abstract type R3 <: World end

# friction cone
abstract type FrictionCone end
abstract type LinearizedCone <: FrictionCone end
abstract type NonlinearCone <:FrictionCone end

# environment
struct Environment{W <: World, F <: FrictionCone}
    surf::Any
	surf_grad::Any
end

function environment_2D_flat(; cone = LinearizedCone)
	Environment{R2,cone}(x -> 0.0, x-> zero(x))
end

function environment_3D_flat(; cone = LinearizedCone)
	Environment{R3,cone}(x -> 0.0, x-> zero(x))
end

function environment_2D(surf; cone = LinearizedCone)
	# Generate two functions: they both take as input a vector and return a vector.
	@variables q[1:1]
	s = surf(q)
	s = Symbolics.simplify.(s)
	ds = Symbolics.jacobian([s], q, simplify = true)
	ds = reshape(ds, 1)

	surf_fast = eval(Symbolics.build_function(s, q))
	surf_grad_fast = eval(Symbolics.build_function(ds, q)[1])

	Environment{R2,cone}(surf_fast, surf_grad_fast)
end

function environment_3D(surf; cone = LinearizedCone)
	# Generate two functions: they both take as input a vector and return a vector.
	@variables q[1:2]
	s = surf(q)
	s = Symbolics.simplify.(s)
	ds = Symbolics.jacobian([s], q, simplify = true)
	ds = reshape(ds, 2)

	surf_fast = eval(Symbolics.build_function(s, q))
	surf_grad_fast = eval(Symbolics.build_function(ds, q)[1])

	Environment{R3,cone}(surf_fast, surf_grad_fast)
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

function rotation(env::Environment{R3,FrictionCone}, q)
	# unit surface normal (3D)
	n = [-1.0 * env.surf_grad(q[1:2]); 1.0]
	ns = n ./ sqrt(transpose(n) * n)

	# world-frame normal
	nw = @SVector [0.0, 0.0, 1.0]

	rot(ns, nw)
end

function rotation(env::Environment{R2,FrictionCone}, q)
	# unit surface normal (3D)
	sg = env.surf_grad(q[1:1])[1]
	n = [-1.0 * sg; 1.0]
	ns = n ./ sqrt(transpose(n) * n)

	# world-frame normal
	nw = @SVector [0.0, 1.0]

	ang = atan(nw[2], nw[1]) - atan(ns[2], ns[1])

	# Rw->s
	SMatrix{2,2}([cos(ang) -sin(ang); sin(ang) cos(ang)]) # this the rotation from world frame to surface frame
end

function rotation_3d(env::Environment{R3,FrictionCone}, q)
	r = rotation(env, q)
end

function rotation_3d(env::Environment{R2,FrictionCone}, q)
	r = rotation(env, q)
	r3d = [r[1,1] 0 r[1,2];
		     0    1   0   ;
		   r[2,1] 0 r[2,2]]
end

function friction_mapping(env::Environment{R2,LinearizedCone})
    SMatrix{1, 2}([1.0 -1.0])
end

function friction_mapping(env::Environment{R3,LinearizedCone})
    SMatrix{2, 4}([1.0 0.0 -1.0 0.0;
                   0.0 1.0 0.0 -1.0])
end

function friction_mapping(env::Environment{R2,NonlinearCone})
    SMatrix{1, 1}([1.0])
end

function friction_mapping(env::Environment{R3,NonlinearCone})
    SMatrix{2, 2}([1.0 0.0;
                   0.0 1.0])
end

dim(env::Environment{R2,FrictionCone}) = 2
dim(env::Environment{R3,FrictionCone}) = 3

friction_dim(env::Environment{R2,LinearizedCone}) = 2
friction_dim(env::Environment{R3,LinearizedCone}) = 4

friction_dim(env::Environment{R2,NonlinearCone}) = 1
friction_dim(env::Environment{R3,NonlinearCone}) = 2

function verify_2D_surface(env::Environment{R2,FrictionCone}, x;
	x_range = range(0.0, stop = 5.0, length = 1000))

	nw = [0.0; 1.0]
	tw = [1.0; 0.0]

	s = [env.surf([xx])[1] for xx in x_range]
	plt = plot(x_range, s, label = "", aspect_ratio = :equal)

	for t_cand in x
		s_cand = env.surf([t_cand])[1]
		p_cand = [t_cand; s_cand]
		tan_norm = rotation(env, [t_cand])' * tw
		norm_norm = rotation(env, [t_cand])' * nw

		plt = plot!([p_cand[1], p_cand[1] + tan_norm[1]],
			[p_cand[2], p_cand[2] + tan_norm[2]], color = :green, width = 2.0, label = "")
		plt = plot!([p_cand[1], p_cand[1] + norm_norm[1]],
			[p_cand[2], p_cand[2] + norm_norm[2]], color = :cyan, width = 2.0, label = "")
	end

	return plt
end

# TODO: verify this method
function verify_3D_surface(env::Environment{R3,FrictionCone}, xl;
	x_range = range(0.0, stop = 5.0, length = 1000))

	nw = [0.0; 0.0; 1.0]
	tw = [1.0; 1.0; 0.0]

	sx = [env.surf([xx, 0])[1] for xx in x_range]
	sy = [env.surf([0, xx])[1] for xx in x_range]
	plt = plot(layout=(2,1), legend=false)

	plot!(plt[1,1], x_range, sx, label = "", aspect_ratio = :equal)
	plot!(plt[2,1], x_range, sy, label = "", aspect_ratio = :equal)

	for x in xl
		sx_cand = env.surf([x[1], 0])[1]
		sy_cand = env.surf([0, x[2]])[1]
		px_cand = [x[1]; sx_cand]
		py_cand = [x[2]; sy_cand]
		rotx = rotation(env, [x[1], 0])
		roty = rotation(env, [0, x[2]])
		t = ([1,1,0])/norm([1,1,0])
		xw = [1,0,0]
		yw = [0,1,0]
		zw = [0,0,1]
		tanx = (rotx'*xw)[[1,3]]
		tany = (roty'*yw)[[2,3]]
		normx = (rotx'*zw)[[1,3]]
		normy = (roty'*zw)[[2,3]]
		@show rotx
		@show roty
		@show tanx
		@show normx
		plot!(plt[1,1], [x[1], x[1]+tanx[1]], [sx_cand, sx_cand+tanx[2]])
		plot!(plt[1,1], [x[1], x[1]+normx[1]], [sx_cand, sx_cand+normx[2]])
		plot!(plt[2,1], [x[2], x[2]+tany[1]], [sy_cand, sy_cand+tany[2]])
		plot!(plt[2,1], [x[2], x[2]+normy[1]], [sy_cand, sy_cand+normy[2]])
		# tan_norm = rotation(env, [t_cand])' * tw
		# norm_norm = rotation(env, [t_cand])' * nw
		#
		# plt = plot!([p_cand[1], p_cand[1] + tan_norm[1]],
		# 	[p_cand[2], p_cand[2] + tan_norm[2]], color = :green, width = 2.0, label = "")
		# plt = plot!([p_cand[1], p_cand[1] + norm_norm[1]],
		# 	[p_cand[2], p_cand[2] + norm_norm[2]], color = :cyan, width = 2.0, label = "")
	end

	return plt
end
