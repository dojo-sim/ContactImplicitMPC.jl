# using Interpolations, Plots
#
# N = 100
# x = range(-1.0, stop = 1.0, length = N)
# function gen_z(x, N)
# 	z = Float64[]
# 	for i = 1:N
# 		if x[i] < -0.5
# 			push!(z, 0.0)
# 		elseif x[i] >= -0.5 && x[i] < 0.5
# 			push!(z, 1.0 * x[i] + 0.5)
# 		else
# 			push!(z, 1.0)
# 		end
# 	end
# 	return z
# end
# z = gen_z(x, N)
# # z = [xi < 0.0 ? 0.0 : 1.0 for xi in x]
# interp = CubicSplineInterpolation(x, z)
#
# plot(x, z, color = :black, width = 2.0)
#
# x_fine = range(-1.0, stop = 1.0, length = 10 * N)
# plot!(x_fine, interp.(x_fine), color = :cyan, width = 1.0)
#
# @variables x_sym[1:1]
# dx = interp(x_sym)
#
# t = range(-π, stop = π, length = N)
# s = sin.(t)
# plot(t, s, ratio = :equal, color = :black, width = 2.0, label = "")
#
# function skew(x)
# 	SMatrix{3,3}([0.0 -x[3] x[2];
# 	               x[3] 0.0 -x[1];
# 				   -x[2] x[1] 0.0])
# end
#
# function rot(a, b)
# 	v = cross(a, b)
# 	s = norm(v)
# 	c = a' * b
#
# 	R = Diagonal(@SVector ones(3)) + skew(v) + 1.0 / (1.0 + c) * skew(v) * skew(v)
# end
using Plots
N = 100
t = range(-π, stop = π, length = N)
s = sin.(t)
plot(t, s, ratio = :equal, color = :black, width = 2.0, label = "")

nw = [0.0; 1.0]
tw = [1.0; 0.0]
env = Environment{R2}(sin, rotation)
t_cand = 2.0
s_cand = sin(t_cand)
p_cand = [t_cand; s_cand]
tan_norm = env.R(env, t_cand)' * tw
norm_norm = env.R(env, t_cand)' * nw

plot!([p_cand[1], p_cand[1] + tan_norm[1]],
	[p_cand[2], p_cand[2] + tan_norm[2]], color = :green, width = 2.0, label = "")
plot!([p_cand[1], p_cand[1] + norm_norm[1]],
	[p_cand[2], p_cand[2] + norm_norm[2]], color = :cyan, width = 2.0, label = "")


t = range(-1.0, stop = 1.0, length = N)
s = t.^2.0
plot(t, s, ratio = :equal, color = :black, width = 2.0, label = "")

env = Environment{R3}(x -> x[1:2]' * x[1:2], rotation)
_p = [0.25; 0.25]
p_cand = [_p; env.surf(_p)]

nw = [0.0; 0.0; 1.0]
tw = [1.0; 0.0; 0.0]
tan_norm = env.R(env, p_cand)' * tw
norm_norm = env.R(env, p_cand)' * nw


plot!([p_cand[1], p_cand[1] + tan_norm[1]],
	[p_cand[3], p_cand[3] + tan_norm[3]], color = :green, width = 2.0, label = "")
plot!([p_cand[1], p_cand[1] + norm_norm[1]],
	[p_cand[3], p_cand[3] + norm_norm[3]], color = :cyan, width = 2.0, label = "")
