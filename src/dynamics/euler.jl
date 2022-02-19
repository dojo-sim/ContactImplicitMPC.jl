#https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions

function euler_rotation_matrix(θ)
	a = θ[1]
	b = θ[2]
	c = θ[3]

	[cos(a) * cos(b) (cos(a) * sin(b) * sin(c) - sin(a) * cos(c)) (cos(a) * sin(b) * cos(c) + sin(a) * sin(c));
	 sin(a) * cos(b) (sin(a) * sin(b) * sin(c) + cos(a) * cos(c)) (sin(a) * sin(b) * cos(c) - cos(a) * sin(c));
	 -sin(b) cos(b) * sin(c) cos(b) * cos(c)]
end
