T = 40
n = 2
m = 1

function f(x1, u1, x2)
	h = 0.05

	function dyn(x, u)
		[x[2];
        (u[1] / 0.25
        - 9.81 * sin(x[1]) / 0.5
        - 0.1 * x[2] / 0.25)]
	end

	x2 - (x1 + h * dyn(0.5 * (x1 + x2), u1))
end

Q = [1.0e-3 0.0; 0.0 1.0e-3]
R = [1.0e-1]

x_init = [0.0; 0.0]
x_final = [π; 0.0]

nz = T * n + (T - 1) * m + (T + 1) * n

x_idx = [collect((t - 1) * n .+ (1:n)) for t = 1:T]
u_idx = [collect(n * T + (t - 1) * m .+ (1:m)) for t = 1:T-1]
y_idx = [collect(n * T + m * (T - 1) + (t - 1) * n .+ (1:n)) for t = 1:T + 1]

function lagrangian(z)
	L = 0.0

	x1 = view(z, x_idx[1])
	y1 = view(z, y_idx[1])

	L += transpose(y1) * (x_init - x1)

	for t = 1:T-1
		xt = view(z, x_idx[t])
		ut = view(z, u_idx[t])
		yt⁺ = view(z, y_idx[t+1])
		xt⁺ = view(z, x_idx[t+1])

		L += transpose(xt - x_final) * Q * (xt - x_final)
		L += (transpose(ut) * R * ut)[1]

		L += transpose(yt⁺) * f(xt, ut, xt⁺)
	end

	xT = view(z, x_idx[T])
	yT⁺ = view(z, y_idx[T+1])

	L += transpose(xT - x_final) * Q * (xT - x_final)
	L += transpose(yT⁺) * (x_final - xT)

	return L
end

z0 = rand(nz)

lagrangian(z0)

@variables z_sym[1:nz]
@variables θ_sym
@variables κ_sym

L = lagrangian(z_sym);
L = simplify.(L);

dL = Symbolics.gradient(L, z_sym)
ddL = Symbolics.sparsehessian(L, z_sym)

L_grad = eval(Symbolics.build_function(dL, z_sym, θ_sym, κ_sym)[1])
L_hess = eval(Symbolics.build_function(ddL, z_sym, θ_sym)[1])

L_grad! = eval(Symbolics.build_function(dL, z_sym, θ_sym, κ_sym)[2])
L_hess! = eval(Symbolics.build_function(ddL, z_sym, θ_sym)[2])

L_grad(z0, nothing, nothing)
L_hess(z0, nothing)

# options
opts = ContactControl.InteriorPointOptions(
	κ_init = 1.0,
	κ_tol = 1.0,
	r_tol = 1.0e-8,
	diff_sol = false)

# solver
ip = ContactControl.interior_point(z0, zeros(0),
	r! = L_grad!, rz! = L_hess!,
	rz = similar(ddL, Float64),
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
