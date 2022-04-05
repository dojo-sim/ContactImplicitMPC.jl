
# # ## model
# include("ci_model.jl") 
# model = RoboDojo.centroidal_quadruped

# # ## horizon
# T = 101
# Tm = 51
# T_fix = 10

# # ## Time step
# tf = 1.0
# h = 0.01

# # ## feet templates
# include("template.jl") 

# function configuration(model; 
#     body_position=[0.0; 0.0; 0.15],
#     body_orientation=[0.0; 0.0; 0.0],
#     foot1_position=[0.15 - 0.025; 0.1; 0.0],
#     foot2_position=[0.1 - 0.025;-0.1; 0.0],
#     foot3_position=[-0.1 - 0.025; 0.1; 0.0],
#     foot4_position=[-0.05 - 0.025;-0.1; 0.0]) 
#     [
#         body_position; 
#         body_orientation;
#         foot1_position;
#         foot2_position;
#         foot3_position;
#         foot4_position;
#     ]
# end

# q1 = configuration(model)

# # ## visualize 
# vis = Visualizer() 
# render(vis)
# # q_sol = state_to_configuration([x[1:nx] for x in x_sol])
# RoboDojo.visualize!(vis, RoboDojo.centroidal_quadruped, [q1], Δt=h);


# # feet positions
# pf1 = q1[6 .+ (1:3)]
# pf2 = q1[9 .+ (1:3)]

# pf3 = q1[12 .+ (1:3)]
# pf4 = q1[15 .+ (1:3)]

# strd = 2 * (pf1 - pf2)[1]

# q_shift1 = zeros(model.nq)
# q_shift1[1] = 0.5 * strd
# q_shift1[10] = strd
# q_shift1[13] = strd

# qM = copy(q1) + q_shift1

# q_shift2 = zeros(model.nq)
# q_shift2[1] = 0.5 * strd
# q_shift2[7] = strd
# q_shift2[16] = strd

# qT = copy(qM) + q_shift2

# q_ref = [q1, linear_interpolation(q1, qM, 30)..., linear_interpolation(qM, qT, 30)...]
# visualize!(vis, model, q_ref)

# zh = 0.05
# xf1_el, zf1_el = ellipse_traj(pf1[1], pf1[1] + strd, zh, Tm - T_fix)
# xf1 = [[pf1[1] for t = 1:Tm + T_fix + 1]..., xf1_el[2:end]...]
# zf1 = [[pf1[3] for t = 1:Tm + T_fix + 1]..., zf1_el[2:end]...]
# pf1_ref = [[xf1[t]; pf1[2];  zf1[t]] for t = 1:T + 1]

# xf4_el, zf4_el = ellipse_traj(pf4[1], pf4[1] + strd, zh, Tm - T_fix)
# xf4 = [[pf4[1] for t = 1:Tm + T_fix + 1]..., xf4_el[2:end]...]
# zf4 = [[pf4[3] for t = 1:Tm + T_fix + 1]..., zf4_el[2:end]...]
# pf4_ref = [[xf4[t]; pf4[2]; zf4[t]] for t = 1:T + 1]

# xf2_el, zf2_el = ellipse_traj(pf2[1], pf2[1] + strd, zh, Tm - T_fix)
# xf2 = [[xf2_el[1] for t = 1:T_fix + 1]..., xf2_el..., [xf2_el[end] for t = 1:Tm-1 + T_fix]...]
# zf2 = [[zf2_el[1] for t = 1:T_fix + 1]..., zf2_el..., [zf2_el[end] for t = 1:Tm-1 + T_fix]...]
# pf2_ref = [[xf2[t]; pf2[2]; zf2[t]] for t = 1:T + 1]

# xf3_el, zf3_el = ellipse_traj(pf3[1], pf3[1] + strd, zh, Tm - T_fix)
# xf3 = [[xf3_el[1] for t = 1:T_fix + 1]..., xf3_el..., [xf3_el[end] for t = 1:Tm-1]...]
# zf3 = [[zf3_el[1] for t = 1:T_fix + 1]..., zf3_el..., [zf3_el[end] for t = 1:Tm-1]...]
# pf3_ref = [[xf3[t]; pf3[2]; zf3[t]] for t = 1:T + 1]

# using Plots
# tr = range(0, stop = tf, length = T+1)
# plot(tr, hcat(pf1_ref...)')
# plot!(tr, hcat(pf4_ref...)')

# plot!(tr, hcat(pf2_ref...)')
# plot!(tr, hcat(pf3_ref...)')

# q_ref = [q1, linear_interpolation(q1, qM, Tm)[1:end-1]..., linear_interpolation(qM, qT, Tm)...,]

# for (t, q) in enumerate(q_ref) 
#     q[6 .+ (1:3)] = pf1_ref[t]
#     q[9 .+ (1:3)] = pf2_ref[t] 
#     q[12 .+ (1:3)] = pf3_ref[t] 
#     q[15 .+ (1:3)] = pf4_ref[t] 
# end

# x_ref = [[q_ref[t]; q_ref[t+1]] for t = 1:length(q_ref)-1]
# visualize!(vis, model, q_ref, Δt=h)

# # ## centroidal_quadruped 
# model = RoboDojo.centroidal_quadruped
# nx = 2 * model.nq
# nc = 4 # model.nc
# nu = model.nu + nc + 4 * nc + nc + 4 * nc + 1
# nθ = 5 

# RoboDojo.mass_matrix(model, ones(model.nq))
# RoboDojo.dynamics_bias(model, ones(model.nq), ones(model.nq))
# RoboDojo.contact_jacobian(model, ones(model.nq))[1:12, :]

# # ## model
# d1 = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dyn1(model, RoboDojo.mass_matrix, RoboDojo.dynamics_bias, [h], y, x, u, w), nx + nθ + nx + model.nu, nx, nu)
# dt = DTO.Dynamics((y, x, u, w) -> centroidal_quadruped_dynt(model, RoboDojo.mass_matrix, RoboDojo.dynamics_bias, [h], y, x, u, w), nx + nθ + nx + model.nu, nx + nθ + nx + model.nu, nu)

# dyn = [d1, [dt for t = 2:T-1]...]

# # ## objective
# obj = DTO.Cost{Float64}[]

# for t = 1:T 

#     if t == 1 
#         function obj1(x, u, w)
#             J = 0.0 
#             v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
#             J += 0.5 * 1.0 * dot(v, v)
#             J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(100.0 * ones(nx)) * (x[1:nx] - x_ref[t]) 
#             J += 0.5 * transpose(u) * Diagonal([ones(model.nu); zeros(nu - model.nu)]) * u
#             J += 1000.0 * u[end] # slack
#             return J
#         end
#         push!(obj, DTO.Cost(obj1, nx, nu))
#     elseif t == T 
#         function objT(x, u, w)
#             J = 0.0 
#             v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
#             J += 0.5 * 1.0 * dot(v, v)
#             J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(100.0 * ones(nx)) * (x[1:nx] - x_ref[t]) 
#             return J
#         end
#         push!(obj, DTO.Cost(objT, nx + nθ + nx + model.nu, 0))
#     else 
#         function objt(x, u, w)
#             J = 0.0 
#             v = (x[model.nq .+ (1:model.nq)] - x[1:model.nq]) ./ h
#             J += 0.5 * 1.0 * dot(v, v)
#             J += 0.5 * transpose(x[1:nx] - x_ref[t]) * Diagonal(100.0 * ones(nx)) * (x[1:nx] - x_ref[t]) 
#             J += 0.5 * transpose(u) * Diagonal([ones(model.nu); zeros(nu - model.nu)]) * u
#             J += 1000.0 * u[end] # slack
#             return J
#         end
#         push!(obj, DTO.Cost(objt, nx + nθ + nx + model.nu, nu))
#     end
# end

# # ## constraints
# # initial condition
# xl1 = x_ref[1] 
# xu1 = [q_ref[1]; q_ref[1]]

# # stage
# xlt = [-Inf * ones(nx); -Inf * ones(nθ); -Inf * ones(nx); -Inf * ones(model.nu)] 
# xut = [Inf * ones(nx); Inf * ones(nθ); Inf * ones(nx); Inf * ones(model.nu)]

# # final condition
# xlT = [q_ref[end]; q_ref[end]; -Inf * ones(nθ); -Inf * ones(nx); -Inf * ones(model.nu)] 
# xuT = [q_ref[end]; q_ref[end]; Inf * ones(nθ); Inf * ones(nx); Inf * ones(model.nu)]

# ul = [-Inf * ones(model.nu); zeros(nu - model.nu)]
# uu = [Inf * ones(model.nu); Inf * ones(nu - model.nu)]

# bnd1 = DTO.Bound(nx, nu, xl=xl1, xu=xu1, ul=ul, uu=uu)
# bndt = DTO.Bound(nx + nθ + nx + model.nu, nu, xl=xlt, xu=xut, ul=ul, uu=uu)
# bndT = DTO.Bound(nx + nθ + nx + model.nu, 0, xl=xlT, xu=xuT)
# bnds = [bnd1, [bndt for t = 2:T-1]..., bndT];

# function constraints_1(x, u, w) 
#     [
#      # equality (16)
#      contact_constraints_equality(model, h, x, u, w); 
#      # inequality (28)
#      contact_constraints_inequality_1(model, h, x, u, w);
#     ]
# end

# function constraints_t(x, u, w) 
#     [
#      # equality (16)
#      contact_constraints_equality(model, h, x, u, w); 
#      # inequality (32)
#      contact_constraints_inequality_t(model, h, x, u, w);
#     ]
# end

# function constraints_T(x, u, w) 
#     [
#      # inequality (8)
#      contact_constraints_inequality_T(model, h, x, u, w);
#     ]
# end

# con1 = DTO.Constraint(constraints_1, nx, nu, idx_ineq=collect(16 .+ (1:28))) 
# cont = DTO.Constraint(constraints_t, nx + nθ + nx + model.nu, nu, idx_ineq=collect(16 .+ (1:32))) 
# conT = DTO.Constraint(constraints_T, nx + nθ + nx + model.nu, nu, idx_ineq=collect(0 .+ (1:8))) 
# cons = [con1, [cont for t = 2:T-1]..., conT];

# # ## problem 
# p = DTO.solver(dyn, obj, cons, bnds, 
#     options=DTO.Options(
#         tol=1.0e-2,
#         constr_viol_tol=1.0e-2,
#         ))

# # ## initialize
# x_interpolation = [x1, [[x1; zeros(nθ); zeros(nx); zeros(model.nu)] for t = 2:T]...]
# u_guess = [1.0e-4 * rand(nu) for t = 1:T-1] # may need to run more than once to get good trajectory
# DTO.initialize_states!(p, x_interpolation)
# DTO.initialize_controls!(p, u_guess)

# # ## solve
# @time DTO.solve!(p)

# # ## solution
# x_sol, u_sol = DTO.get_trajectory(p)
# @show x_sol[1]
# @show x_sol[T]
# sum([u[end] for u in u_sol[1:end-1]])

# # ## visualize 
# vis = Visualizer() 
# render(vis)
# # q_sol = state_to_configuration([x[1:nx] for x in x_sol])
# RoboDojo.visualize!(vis, RoboDojo.centroidal_quadruped, x_sol, Δt=h);
