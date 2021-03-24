function meas(x)
	q1, q2, γ1, b1 = unpack_x(x)
	y = q2
	return y
end

function unpack_x(x)
	off = 0
	q1 = x[off .+ (1:nq)]; off += nq
	q2 = x[off .+ (1:nq)]; off += nq
	γ1 = x[off .+ (1:nc)]; off += nc
	b1 = x[off .+ (1:nb)]; off += nb
	return q1, q2, γ1, b1
end
