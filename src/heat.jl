"""
ms_in -- входящий материальный поток
W -- входящий тепловой поток
Δp -- перепад давления
"""
function heater_W(ms_in::MaterialStream, W::EnergyStream, Δp::Float64 = 0.0) #W поток в ваттах
	h_in = mstream_H_T(ms_in, ms_in.T)    #начальная удельная энтальпия
    h_k = h_in + W.Q / ms_in.G	#посчитать конечную энтальпия
	T_k = find_zero(T -> mstream_H_T(ms_in, T) - h_k, ms_in.T)
	#возвращаем streamout
	ms_out = copy_stream(ms_in)
	ms_out.T = T_k
    A = @. ms_in.Q * ms_in.y + (1.0-ms_in.Q) * ms_in.x
    p_k = ms_in.p
    return mstream_TpA(ms_in.G, T_k, p_k, A, ms_in.model)
end

"""
ms_in -- входящий материальный поток
W -- входящий тепловой поток
Δp -- перепад давления
"""
function heater_Q(ms_in::MaterialStream,Q::Float64, Δp::Float64 = 0.0)
	h_n = mstream_H_T(ms_in, ms_in.T)	 #начальная удельная энтальпия
	A = @. ms_in.Q * ms_in.y + (1.0-ms_in.Q) * ms_in.x	#исходный состав
	ms_out = mstream_pQA(ms_in.G, p, A, Q, ms_in.model)	#конечный состав
	h_k = mstream_H_T(ms_out, ms_out.T)		#конечная энтальпия
	dQ = (h_k - h_n) * sout.G		#подведенная теплота
	return ms_out, EnergyStream(dQ)
end

