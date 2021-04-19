function s_func(model::PinnedBiped513, q;
		# q_ini=[-0.40,  0.70, 0.0, 0.53, -0.584],
		# q_mid=[-0.44,  0.54, 0.0, 0.66, -0.56 ],
		# q_end=[-0.584, 0.53, 0.0, 0.70, -0.40 ],
		ref=REF,
		)
	q_ini, q_mid, q_end = ref
	θcom_ini = θcom_func(model, q_ini)
	θcom_end = θcom_func(model, q_end)

	θcom = θcom_func(model, q)
	s = (θcom - θcom_ini)/(θcom_end - θcom_ini)
	# s = clamp(s, 0.0, 1.0)
	return s
end

function sd_func(model::PinnedBiped513, q, qd;
		# q_ini=[-0.40,  0.70, 0.0, 0.53, -0.584],
		# q_mid=[-0.44,  0.54, 0.0, 0.66, -0.56 ],
		# q_end=[-0.584, 0.53, 0.0, 0.70, -0.40 ],
		ref=REF,
		)
	q_ini, q_mid, q_end = ref
	θcom_(q) = θcom_func(model, q)
	∇qθcom = ForwardDiff.gradient(θcom_, q)
	θcomd = ∇qθcom'*qd

	θcom_ini = θcom_func(model, q_ini)
	θcom_end = θcom_func(model, q_end)
	sd = θcomd/(θcom_end - θcom_ini)

	return sd
end

function p_func(model::PinnedBiped513, s;
		# q_ini=[-0.40,  0.70, 0.0, 0.53, -0.584],
		# q_mid=[-0.44,  0.54, 0.0, 0.66, -0.56 ],
		# q_end=[-0.584, 0.53, 0.0, 0.70, -0.40 ],
		ref=REF,
		)
	nh = 4
	q_ini, q_mid, q_end = ref
	θcom_ini = θcom_func(model, q_ini)
	θcom_mid = θcom_func(model, q_mid)
	θcom_end = θcom_func(model, q_end)

	h_ini = h_func(model, q_ini)
	h_mid = h_func(model, q_mid)
	h_end = h_func(model, q_end)

	s_ini = 0.0
	s_mid = (θcom_mid - θcom_ini)/(θcom_end - θcom_ini)
	s_end = 1.0

	p = zeros(nh)
	if s < s_mid
		α = (s - s_ini)/(s_mid - s_ini)
		p = h_mid*α + (1-α)*h_ini
	elseif s >= 0*s+s_mid
		α = (s - s_mid)/(s_end - s_mid)
		p = h_end*α + (1-α)*h_mid
	end
	return p
end

function pd_func(model::PinnedBiped513, s;
		# q_ini=[-0.40,  0.70, 0.0, 0.53, -0.584],
		# q_mid=[-0.44,  0.54, 0.0, 0.66, -0.56 ],
		# q_end=[-0.584, 0.53, 0.0, 0.70, -0.40 ],
		ref=REF,
		)
	nh = 4
	q_ini, q_mid, q_end = ref
	θcom_ini = θcom_func(model, q_ini)
	θcom_mid = θcom_func(model, q_mid)
	θcom_end = θcom_func(model, q_end)

	h_ini = h_func(model, q_ini)
	h_mid = h_func(model, q_mid)
	h_end = h_func(model, q_end)

	s_ini = 0.0
	s_mid = (θcom_mid - θcom_ini)/(θcom_end - θcom_ini)
	s_end = 1.0

	pd = zeros(nh)
	if s < s_mid
		pd = (h_mid - h_ini)/(s_mid - s_ini)
	elseif s >= s_mid
		pd = (h_end - h_mid)/(s_end - s_mid)
	end
	return pd
end

function pdd_func(model::PinnedBiped513, s;
		# q_ini=[-0.40,  0.70, 0.0, 0.53, -0.584],
		# q_mid=[-0.44,  0.54, 0.0, 0.66, -0.56 ],
		# q_end=[-0.584, 0.53, 0.0, 0.70, -0.40 ],
		ref=REF,
		)
	nh = 4
	q_ini, q_mid, q_end = ref
	pdd = zeros(nh)
	return pdd
end
