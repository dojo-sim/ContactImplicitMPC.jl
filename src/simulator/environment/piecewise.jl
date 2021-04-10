# https://math.stackexchange.com/questions/2692572/smooth-approximation-of-three-phased-linear-models
poly_func(x, b, s) = 0.5 * sqrt(b * (4 * s + b * x^2.0))

function comp_func(x, pts, slopes, dirs, sf)
	a = dirs[1] * slopes[1] * x[1]
	for i = 2:length(slopes)
		a += dirs[i] * poly_func(x[1] - pts[i-1], slopes[i] - slopes[i-1], sf)
	end
	return a
end

function comp_func_rot(x, pts, slopes, dirs, sf)
	comp_func(x, pts, slopes, dirs, sf) - (comp_func(0.0, pts, slopes, dirs, sf) + x[1] * (comp_func(pts[1], pts, slopes, dirs, sf) - comp_func(0.0, pts, slopes, dirs, sf)) / (pts[1] - 0.0))
end

function generate_fast_terrain_2D(pts, slopes, dirs, sf)
	@variables x[1:1]
	s = comp_func_rot(x, pts, slopes, dirs, sf)
	s = Symbolics.simplify.(s)

	fast_terrain = eval(Symbolics.build_function(s, x))

	return fast_terrain
end

# example

#pts
x1 = 1.0
x2 = 2.0
x3 = 3.0
x4 = 4.0
pts = [x1, x2, x3, x4]

#slopes
s1 = 0.5
s2 = 0.65
s3 = 0.75
s4 = 1.0
s5 = 1.125
slopes = [s1, s2, s3, s4, s5]

# directions
dirs = [1.0, -1.0, 1.0, 1.0, -1.0]

# smoothness
sf = 0.0001

fast_terrain = generate_fast_terrain_2D(pts, slopes, dirs, sf)
# env = environment_2D(fast_terrain)
# x_terrain = range(-1.0, stop = 5.0, length = 1000)
# z_terrain = [env.surf([x])[1] for x in x_terrain]
#
# plot(x_terrain, z_terrain, aspect_ratio = :equal)
# verify_2D_surface(env, [0.0, 1.0, 2.0, 3.0, 4.0], x_range = x_terrain)
