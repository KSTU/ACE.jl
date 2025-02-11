"""
ms_in -- входящий материальный поток
W -- входящий тепловой поток
Δp -- перепад давления

"""
function heater_W(ms_in::MaterialStream, W::EnergyStream, Δp::Float64 = 0.0) #W поток в ваттах
	h_in = mstream_H_T(ms_in, ms_in.T)    #начальная
    #посчитать конечную энтальпия
    h_k = h_n + W.Q / ms_in.G
	T_k = find_zero(T -> mstream_H_T(ms_in, T) - h_k, ms_in.T)
	#возвращаем streamout
	ms_out = copy_stream(ms_in)
	ms_out.T = T_k
    A = @. ms_in.Q * ms_in.y + (1.0-ms_in.Q) * ms_in.x
    p_k = ms_in.p - Δp
    return mstream_TpA(ms_in.G, T_k, p_k, A, ms_in.model)
end

"""
ms_in -- входящий материальный поток
W -- входящий тепловой поток
Δp -- перепад давления

"""
function heater_Q(ms_in::MaterialStream,Q::Float64, Δp::Float64 = 0.0)
	h_n = mstream_H_T(ms_in, ms_in.T)
	A = @. ms_in.Q * ms_in.y + (1.0-ms_in.Q) * ms_in.x

	ms_out = mstream_pQA(ms_in.G, p, A, Q, ms_in.model)
	H_k = myH(sout, sout.T)
	dQ = (H_k - H_n) * sout.G
	return [sout, dQ] 
end