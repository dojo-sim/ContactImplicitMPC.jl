# Smooth ReLU
t_ss = 25.0
m_ss = tan(deg2rad(10.0)) # 20 degree slope
x_off_ss = 0.5
slope_smooth(x) = m_ss * 1.0 / t_ss .* log.(1.0 .+ exp.(t_ss .* (x[1] .- x_off_ss)))

# x = range(-1.0, stop = 4.0, length = 1000)
# plot(x, ss.(x), aspect_ratio = :equal)

slope_smooth_2D_lc = environment_2D(slope_smooth)

slope1_2D_lc = environment_2D(x -> 0.5 * x[1])
