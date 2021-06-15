sine1_3D_lc = environment_3D(x -> sin(x[1]) + sin(x[2]))
sine2_3D_lc = environment_3D(x -> 0.075 * sin(2π * x[1]))
sine3_3D_lc = environment_3D(x -> 0.075 * sin(2π * x[1]) * sin(2π * x[2]))

sine1_2D_lc = environment_2D(x -> 0.05 * (cos(pi * x[1]) - 1.0))
sine2_2D_lc = environment_2D(x -> 0.10 * sin(2π * x[1]))
sine3_2D_lc = environment_2D(x -> 0.03 * (cos(pi * x[1]) - 1.0))
