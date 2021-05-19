function terrain(x)
	IfElse.ifelse(x[1] < 1.0, 0.0,
		IfElse.ifelse(x[1] < 2.0, -0.125 * x[1] + 0.125,
			IfElse.ifelse(x[1] < 3.0, -0.075 * x[1] + 0.025,
				IfElse.ifelse(x[1] < 4.0, 0.3 * x[1] - 1.1,
					0.1))))
end

function d_terrain(x)
	IfElse.ifelse(x[1] < 1.0, 0.0,
		IfElse.ifelse(x[1] < 2.0, -0.125,
			IfElse.ifelse(x[1] < 3.0, -0.075,
				IfElse.ifelse(x[1] < 4.0, 0.3,
					0.0))))
end
