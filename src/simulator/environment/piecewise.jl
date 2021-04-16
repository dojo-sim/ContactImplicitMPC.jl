function terrain(x)
	if x[1] < 1.0
		return [0.0]
	elseif x[1] < 2.0
		return [-0.125 * x[1] + 0.125]
	elseif x[1] < 3.0
		return [-0.075 * x[1] + 0.025]
	elseif x[1] < 4.0
		return [0.3 * x[1] - 1.1]
	else
		return [0.1]
	end
end

function terrain_sym(x)
	IfElse.ifelse(x[1] < 1.0, 0.0,
		IfElse.ifelse(x[1] < 2.0, -0.125 * x[1] + 0.125,
			IfElse.ifelse(x[1] < 3.0, -0.075 * x[1] + 0.025,
				IfElse.ifelse(x[1] < 4.0, 0.3 * x[1] - 1.1,
					0.1))))

	# IfElse.ifelse(x[1] < 1.0, 0.0, 1.0 * x[1])
	# [0.0]
end

function d_terrain(x)
	if x[1] < 1.0
		return [0.0]
	elseif x[1] < 2.0
		return [-0.125]
	elseif x[1] < 3.0
		return [-0.075]
	elseif x[1] < 4.0
		return [0.3]
	else
		return [0.0]
	end
end

function d_terrain_sym(x)
	IfElse.ifelse(x[1] < 1.0, 0.0,
		IfElse.ifelse(x[1] < 2.0, -0.125,
			IfElse.ifelse(x[1] < 3.0, -0.075,
				IfElse.ifelse(x[1] < 4.0, 0.3,
					0.0))))
	# IfElse.ifelse(x[1] < 1.0, 0.0, 1.0)
	# [0.0]
end

# hopper test
# include(joinpath(pwd(), "src/dynamics/hopper_2D/model.jl"))
# model = Hopper2D(Dimensions(nq, nu, nw, nc, nb),
# 			   mb, ml, Jb, Jl,
# 			   μ_world, μ_joint, g,
# 			   BaseMethods(), DynamicsMethods(), ResidualMethods(), ResidualMethods(),
# 			   SparseStructure(spzeros(0, 0), spzeros(0, 0)),
# 			   SVector{4}(zeros(4)),
# 			   Environment{R2}(terrain_sym, d_terrain_sym))

# quadruped test
# include(joinpath(pwd(), "src/dynamics/quadruped/model.jl"))

# model = Quadruped(Dimensions(nq, nu, nw, nc, nb),
# 				g, μ_world, μ_joint,
# 				l_torso, d_torso, m_torso, J_torso,
# 				l_thigh, d_thigh, m_thigh, J_thigh,
# 				l_leg, d_leg, m_leg, J_leg,
# 				l_thigh, d_thigh, m_thigh, J_thigh,
# 				l_leg, d_leg, m_leg, J_leg,
# 				l_thigh, d_thigh, m_thigh, J_thigh,
# 				l_leg, d_leg, m_leg, J_leg,
# 				l_thigh, d_thigh, m_thigh, J_thigh,
# 				l_leg, d_leg, m_leg, J_leg,
# 				zeros(nc),
# 				BaseMethods(), DynamicsMethods(), ResidualMethods(), ResidualMethods(),
# 				SparseStructure(spzeros(0, 0), spzeros(0, 0)),
# 				SVector{nq}([zeros(3); μ_joint * ones(nq - 3)]),
# 				Environment{R2}(terrain_sym, d_terrain_sym))

# dir = joinpath(pwd(), "src/dynamics/quadruped")
# dir = joinpath(pwd(), "src/dynamics/hopper_2D")
# path_base = joinpath(dir, "piecewise/base.jld2")
# path_dyn = joinpath(dir, "piecewise/dynamics.jld2")
# path_res = joinpath(dir, "piecewise/residual.jld2")
# path_jac = joinpath(dir, "piecewise/sparse_jacobians.jld2")
# path_linearized = joinpath(dir, "piecewise/linearized.jld2")
#
# function contact_forces(model::Hopper2D, γ1, b1, q2)
# 	k = kinematics(model, q2)
# 	m = friction_mapping(model.env)
# 	SVector{2}(transpose(rotation(model.env, k)) * [m * b1; γ1])
# end
#
# function velocity_stack(model::Hopper2D, q1, q2, h)
# 	k = kinematics(model, q2)
# 	v = J_func(model, q2) * (q2 - q1) / h[1]
# 	v_surf = rotation(model.env, k) * v
# 	SVector{2}([v_surf[1]; -v_surf[1]])
# end
#
# function contact_forces(model::Quadruped, γ1, b1, q2)
# 	k = kinematics(model, q2)
# 	m = friction_mapping(model.env)
#
# 	SVector{8}([transpose(rotation(model.env, k[1:2])) * [m * b1[1:2]; γ1[1]];
# 				transpose(rotation(model.env, k[3:4])) * [m * b1[3:4]; γ1[2]];
# 				transpose(rotation(model.env, k[5:6])) * [m * b1[5:6]; γ1[3]];
# 				transpose(rotation(model.env, k[7:8])) * [m * b1[7:8]; γ1[4]]])
# end
#
# function velocity_stack(model::Quadruped, q1, q2, h)
# 	k = kinematics(model, q2)
# 	v = J_func(model, q2) * (q2 - q1) / h[1]
#
# 	v1_surf = rotation(model.env, k[1:2]) * v[1:2]
# 	v2_surf = rotation(model.env, k[3:4]) * v[3:4]
# 	v3_surf = rotation(model.env, k[5:6]) * v[5:6]
# 	v4_surf = rotation(model.env, k[7:8]) * v[7:8]
#
# 	SVector{8}([v1_surf[1]; -v1_surf[1];
# 				v2_surf[1]; -v2_surf[1];
# 				v3_surf[1]; -v3_surf[1];
# 				v4_surf[1]; -v4_surf[1]])
# end
#
#
# # expr_base = generate_base_expressions(model)
# expr_base = generate_base_expressions_analytical(model)
# save_expressions(expr_base, path_base, overwrite=true)
# instantiate_base!(model, path_base)
#
# expr_dyn = generate_dynamics_expressions(model)
# save_expressions(expr_dyn, path_dyn, overwrite=true)
# instantiate_dynamics!(model, path_dyn)
#
# nq = model.dim.q
# nu = model.dim.u
# nc = model.dim.c
# nb = model.dim.b
# nz = num_var(model)
# nθ = num_data(model)
#
# # Declare variables
# @variables z[1:nz]
# @variables θ[1:nθ]
# @variables κ
#
# # Residual
# r = residual(model, z, θ, κ)
# r = Symbolics.simplify.(r)
# rz = Symbolics.jacobian(r, z)#, simplify = true)
# # rzs = Symbolics.sparsejacobian(r, z)#, simplify = true)
# # rzs.nzval
# # rθ = Symbolics.jacobian(r, θ, simplify = true) # TODO: sparse version
#
# rz_sp = similar(rz, Float64)
# # rzs_sp = similar(rzs, Float64)
# # rθ_sp = similar(rθ, T)
# # rzs_sp.nzval
#
# # Build function
# expr = Dict{Symbol, Expr}()
# expr[:r]  = build_function(r, z, θ, κ)[2]
# expr[:rz] = build_function(rz, z, θ)[2]
# # expr[:rzs] = build_function(rzs.nzval, z, θ)[2]
# # expr[:rθ] = build_function(rθ, z, θ)[2]
#
# r0 = zeros(num_var(model))
# z0 = rand(num_var(model))
# θ0 = rand(num_data(model))
# eval(expr[:r])(r0, z0, θ0, 0.1)
# @time eval(expr[:rz])(rz_sp, z0, θ0)
# # @time eval(expr[:rzs])(rzs_sp.nzval, z0, θ0)
