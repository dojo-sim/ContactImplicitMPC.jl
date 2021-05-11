m_ss = tan(deg2rad(10.0)) # 10 degree slope

function piecewise(x)
	IfElse.ifelse(x[1] < 0.5, 0.0,
		IfElse.ifelse(x[1] < 2.0, m_ss * x[1] - 0.5 * m_ss,
			-0.25 * m_ss * (x[1] - 2.0) + 1.5 * m_ss))
end

function d_piecewise(x)
	IfElse.ifelse(x[1] < 0.5, 0.0,
		IfElse.ifelse(x[1] < 2.0, m_ss,
			-0.25 * m_ss))
end

# smooth kinks w/ cubic polynomial
poly(a, z) = a[4] + a[3] * z + a[2] * z^2.0 + a[1] * z^3.0
d_poly(a, z) = a[3] + 2.0 * a[2] * z + 3.0 * a[1] * z^2.0

# piece 1
m1 = 0.0
x1 = 0.4
y1 = m1 * x1

m2 = m_ss
x2 = 0.6
y2 = m2 * 0.1

A1 = [x1^3.0 x1^2.0 x1 1.0;
      x2^3.0 x2^2.0 x2 1.0;
	  3.0 * x1^2.0 2.0 * x1 1.0 0.0;
	  3.0 * x2^2.0 2.0 * x2 1.0 0.0]
b1 = [y1; y2; m1; m2]

a1 = A1 \ b1

@assert isapprox(poly(a1, x1) - y1, 0.0, atol = 1.0e-8)
@assert isapprox(poly(a1, x2) - y2, 0.0, atol = 1.0e-8)

# piece 2
m1 = m_ss
x1 = 1.4
y1 = m_ss * x1

m2 = -0.25 * m_ss
x2 = 1.6
y2 = m_ss * 1.5 + m2 * 0.1

A2 = [x1^3.0 x1^2.0 x1 1.0;
      x2^3.0 x2^2.0 x2 1.0;
	  3.0 * x1^2.0 2.0 * x1 1.0 0.0;
	  3.0 * x2^2.0 2.0 * x2 1.0 0.0]
b2 = [y1; y2; m1; m2]

a2 = A2 \ b2

@assert isapprox(poly(a2, x1) - y1, 0.0, atol = 1.0e-8)
@assert isapprox(poly(a2, x2) - y2, 0.0, atol = 1.0e-8)


function piecewise_smoothed(x)
	IfElse.ifelse(x[1] < 0.4, 0.0,
		IfElse.ifelse(x[1] < 0.6, poly(a1, x[1]),
			IfElse.ifelse(x[1] < 1.9, m_ss * x[1] - 0.5 * m_ss,
				IfElse.ifelse(x[1] < 2.1, poly(a2, x[1] - 0.5),
				-0.25 * m_ss * (x[1] - 2.0) + 1.5 * m_ss))))
end

function d_piecewise_smoothed(x)
	IfElse.ifelse(x[1] < 0.4, 0.0,
		IfElse.ifelse(x[1] < 0.6, d_poly(a1, x[1]),
			IfElse.ifelse(x[1] < 1.9, m_ss,
				IfElse.ifelse(x[1] < 2.1, d_poly(a2, x[1] - 0.5),
				-0.25 * m_ss))))
end

x = range(-1.0, stop = 4.0, length = 1000)
plot(x, piecewise.(x))#, aspect_ratio = :equal)
plot!(x, piecewise_smoothed.(x))#, aspect_ratio = :equal)

plot(x, d_piecewise.(x), aspect_ratio = :equal)
plot!(x, d_piecewise_smoothed.(x), aspect_ratio = :equal)
