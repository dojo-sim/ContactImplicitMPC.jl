quadratic_bowl_3D_lc = environment_3D(x -> transpose(x[1:2]) * x[1:2], cone = LinearizedCone)
quadratic_bowl_3D_nc = environment_3D(x -> transpose(x[1:2]) * x[1:2], cone = NonlinearCone)
