using Plots
# https://math.stackexchange.com/questions/2692572/smooth-approximation-of-three-phased-linear-models
B1 = 1.0
B2 = 2.0
B3 = 3.0
B4 = 4.0
b1 = 2.0
b2 = 2.5
b3 = 2.5
b4 = 3.5
b5 = 4.0
a = 0.0
q(x, b, s) = 0.5 * sqrt(b * (4 * s + b * x^2.0))

function f1(x, a, b1, b2, b3, b4, b5, B1, B2, B3, B4, s)
	a + b1 * x - q(x - B1, b2 - b1, s) + q(x - B2, b3 - b2, s) + q(x - B3, b4 - b3, s) - q(x - B4, b5 - b4, s)
end

s = 0.01

x_first = -1.0
x_last = 5.0
t = range(x_first, stop = x_last, length = 1000)
y_first = f1.(x_first, a, b1, b2, b3, b4, b5, B1, B2, B3, B4, s)
y_last = f1.(x_last, a, b1, b2, b3, b4, b5, B1, B2, B3, B4, s)
(y_last - y_first) / (x_last - x_first)

function f2(x)
	# y_first + f1(x, a, b1, b2, b3, B1, B2, s) * (y_last - y_first) / (x_last - x_first)
	# f1(x, a, b1, b2, b3, b4, b5, B1, B2, B3, B4, s) - (y_first + x * (y_last - y_first) / (x_last - x_first))
	f1(x[1], a, b1, b2, b3, b4, b5, B1, B2, B3, B4, s) - (y_first + x[1] * (f1(B1, a, b1, b2, b3, b4, b5, B1, B2, B3, B4, s) - f1(0, a, b1, b2, b3, b4, b5, B1, B2, B3, B4, s)) / (B1 - 0.0))

end

plot(t, f1.(t, a, b1, b2, b3, b4, b5, B1, B2, B3, B4, s), label = "", aspect_ratio = :equal)
plot(t, f2.(t), label = "", aspect_ratio = :equal)

env = environment_2D(f2)
verify_2D_surface(env, [0.0, 1.0, 2.0, 3.0, 4.0], x_range = range(-1.0, stop = 5.0, length = 1000))

##
# https://math.stackexchange.com/questions/2692572/smooth-approximation-of-three-phased-linear-models
x1 = 1.0
x2 = 2.0
x3 = 3.0
x4 = 4.0
s1 = 2.0
s2 = 2.5
s3 = 2.5
s4 = 3.5
s5 = 4.0
poly_func(x, b, s) = 0.5 * sqrt(b * (4 * s + b * x^2.0))

pts = [x1, x2, x3, x4]
slopes = [s1, s2, s3, s4, s5]
dirs = [1.0, -1.0, 1.0, 1.0, -1.0]
sf = 0.01

function comp_func(x, pts, slopes, dirs, sf)
	a = dirs[1] * slopes[1] * x[1]
	for i = 2:length(slopes)
		a += dirs[i] * poly_func(x[1] - pts[i-1], slopes[i] - slopes[i-1], sf)
	end
	return a
end
comp_func(0.0, pts, slopes, dirs, sf)

function comp_func_rot(x, pts, slopes, dirs, sf)
	comp_func(x, pts, slopes, dirs, sf) - (comp_func(0.0, pts, slopes, dirs, sf) + x[1] * (comp_func(pts[1], pts, slopes, dirs, sf) - comp_func(0.0, pts, slopes, dirs, sf)) / (pts[1] - 0.0))
end

cf = [comp_func(tt, pts, slopes, dirs, sf) for tt in t]
cfr = [comp_func_rot(tt, pts, slopes, dirs, sf) for tt in t]
plot(t, cf, label = "", aspect_ratio = :equal)
plot(t, cfr, label = "", aspect_ratio = :equal)

@variables x[1:1]
s = comp_func_rot(x, pts, slopes, dirs, sf)
s = Symbolics.simplify.(s)

fast_terrain = eval(Symbolics.build_function(s, x))
fast_terrain(10.0)
