T = 100
n = 2
m = 1

function f(x1, u1, x2)
	h = 0.01

	function dyn(x, u)
		[x[2];
        (u[1] / 0.25
        - 9.81 * sin(x[1]) / 0.5
        - 0.1 * x[2] / 0.25)]
	end

	x2 - (x1 + h * dyn(0.5 * (x1 + x2), u1))
end

Q = [1.0 0.0; 0.0 1.0]
R = [1.0]

x_init = [0.0; 0.0]
x_final = [π; 0.0]
x_ref = linear_interpolation(x_init, x_final, T)
u_ref = [randn(m) for t = 1:T-1]
y_ref = [zeros(n) for t = 1:T+1]

nz = T * n + (T - 1) * m + (T + 1) * n

x_idx = [collect((t - 1) * n .+ (1:n)) for t = 1:T]
u_idx = [collect(n * T + (t - 1) * m .+ (1:m)) for t = 1:T-1]
y_idx = [collect(n * T + m * (T - 1) + (t - 1) * n .+ (1:n)) for t = 1:T + 1]

function obj(x, u)
	J = 0.0

	J += transpose(x - x_final) * Q * (x - x_final)
	J += (transpose(u) * R * u)[1]

	return J
end

function initial_condition(x, y)
	transpose(y) * (x_init - x)
end

function final_condition(x, y)
	transpose(y) * (x_final - x)
end

function dynamics(x1, u1, x2, y2)
	transpose(y2) * f(x1, u1, x2)
end

z0 = rand(nz)

@variables z_sym[1:nz]
@variables x_sym[1:n]
@variables x2_sym[1:n]
@variables u_sym[1:m]
@variables y_sym[1:n]
@variables θ_sym
@variables κ_sym

# L = lagrangian(z_sym);
# L = simplify.(L);
#
# dL = Symbolics.gradient(L, z_sym)
# ddL = Symbolics.sparsehessian(L, z_sym)
#
# L_grad = eval(Symbolics.build_function(dL, z_sym, θ_sym, κ_sym)[1])
# L_hess = eval(Symbolics.build_function(ddL, z_sym, θ_sym)[1])
#
# L_grad! = eval(Symbolics.build_function(dL, z_sym, θ_sym, κ_sym)[2])
# L_hess! = eval(Symbolics.build_function(ddL, z_sym, θ_sym)[2])
#
# L_grad(z0, nothing, nothing)
# L_hess(z0, nothing)
J = obj(x_sym, u_sym)
J = simplify(J)
dJ = Symbolics.gradient(J, [x_sym; u_sym], simplify = true)
ddJ = Symbolics.hessian(J, [x_sym; u_sym], simplify = true)
J_grad = eval(Symbolics.build_function(dJ, x_sym, u_sym)[1])
J_hess = eval(Symbolics.build_function(ddJ, x_sym, u_sym)[1])

ic = initial_condition(x_sym, y_sym)
ic = simplify(ic)
dic = Symbolics.gradient(ic, [x_sym; y_sym], simplify = true)
ddic = Symbolics.hessian(ic, [x_sym; y_sym], simplify = true)
ic_grad = eval(Symbolics.build_function(dic, x_sym, y_sym)[1])
ic_hess = eval(Symbolics.build_function(ddic, x_sym, y_sym)[1])

dyn = dynamics(x_sym, u_sym, x2_sym, y_sym)
ddyn = Symbolics.gradient(dyn, [x_sym; u_sym; x2_sym; y_sym], simplify = true)
dddyn = Symbolics.hessian(dyn, [x_sym; u_sym; x2_sym; y_sym], simplify = true)
dyn_grad = eval(Symbolics.build_function(ddyn, x_sym, u_sym, x2_sym, y_sym)[1])
dyn_hess = eval(Symbolics.build_function(dddyn, x_sym, u_sym, x2_sym, y_sym)[1])

fc = final_condition(x_sym, y_sym)
fc = simplify(fc)
dfc = Symbolics.gradient(fc, [x_sym; y_sym], simplify = true)
ddfc = Symbolics.hessian(fc, [x_sym; y_sym], simplify = true)
fc_grad = eval(Symbolics.build_function(dfc, x_sym, y_sym)[1])
fc_hess = eval(Symbolics.build_function(ddfc, x_sym, y_sym)[1])


function r!(r, z, θ, κ)
	r .= 0.0

	x1 = view(z, x_idx[1])
	y1 = view(z, y_idx[1])

	ic_idx = collect([x_idx[1]; y_idx[1]])
	r[ic_idx] += ic_grad(x1, y1)

	for t = 1:T-1
		xt = view(z, x_idx[t])
		ut = view(z, u_idx[t])
		yt⁺ = view(z, y_idx[t+1])
		xt⁺ = view(z, x_idx[t+1])

		J_idx = collect([x_idx[t]; u_idx[t]])
		r[J_idx] += J_grad(xt, ut)

		dyn_idx = collect([x_idx[t]; u_idx[t]; x_idx[t+1]; y_idx[t+1]])
		r[dyn_idx] += dyn_grad(xt, ut, xt⁺, yt⁺)
	end

	xT = view(z, x_idx[T])
	yT⁺ = view(z, y_idx[T+1])

	fc_idx = collect([x_idx[T]; y_idx[T+1]])
	r[fc_idx] += fc_grad(xT, yT⁺)

	nothing
end

# r0 = zeros(nz)
# z0 = rand(nz)
# r!(r0, z0, nothing, nothing)

function rz!(rz, z, θ)
	rz .= 0.0

	x1 = view(z, x_idx[1])
	y1 = view(z, y_idx[1])

	ic_idx = collect([x_idx[1]; y_idx[1]])
	rz[ic_idx, ic_idx] += ic_hess(x1, y1)

	for t = 1:T-1
		xt = view(z, x_idx[t])
		ut = view(z, u_idx[t])
		yt⁺ = view(z, y_idx[t+1])
		xt⁺ = view(z, x_idx[t+1])

		J_idx = collect([x_idx[t]; u_idx[t]])
		rz[J_idx, J_idx] += J_hess(xt, ut)

		dyn_idx = collect([x_idx[t]; u_idx[t]; x_idx[t+1]; y_idx[t+1]])
		rz[dyn_idx, dyn_idx] += dyn_hess(xt, ut, xt⁺, yt⁺)
	end

	xT = view(z, x_idx[T])
	yT⁺ = view(z, y_idx[T+1])

	fc_idx = collect([x_idx[T]; y_idx[T+1]])
	rz[fc_idx, fc_idx] += fc_hess(xT, yT⁺)

	nothing
end
# rz0 = zeros(nz, nz)
# rz!(rz0, z0, nothing)

# options
opts = ContactControl.InteriorPointOptions(
	κ_init = 1.0,
	κ_tol = 1.0,
	r_tol = 1.0e-5,
	diff_sol = false)

# solver
zz = vcat(vcat(x_ref...), vcat(u_ref...), vcat(y_ref...))
ip = ContactControl.interior_point(z0, zeros(0),
	r! = r!, rz! = rz!,
	rz = zeros(nz, nz),
	opts = opts)

# solve
status = ContactControl.interior_point_solve!(ip)

# test
@test status
@test norm(ip.z[x_idx[1]] - x_init) < 1.0e-8

ip.z[x_idx[T]]

x_traj = [ip.z[x_idx[t]] for t = 1:T]
u_traj = [ip.z[u_idx[t]] for t = 1:T-1]

plot(hcat(x_traj...)')
plot(hcat(u_traj...)')
