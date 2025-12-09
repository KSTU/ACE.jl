"""
Делит поток на 2 части где задается мольная доля первого потока
"""
function splitter_dol(ms_in::MaterialStream, dol::Float64)
    if dol >= 0 && dol <= 1
        ms_out1 = material_stream_copy(ms_in)
        ms_out1.N = ms_in.N * dol
        ms_out2 = material_stream_copy(ms_in)
        ms_out2.N = ms_in.N - ms_out1.N
    return ms_out1, ms_out2
    else 
        println("Доля должна быть в диапазоне от 0 до 1")
    end
end

"""
Делит поток на 2 части где задается мольный расход первого потока
"""
function splitter_N(ms_in::MaterialStream, N1::Float64)
    if 0.0 <= N1 <= ms_in.N 
        ms_out1 = material_stream_copy(ms_in)
        ms_out1.N = N1
        ms_out2 = material_stream_copy(ms_in)
        ms_out2.N = ms_in.N - ms_out1.N
        return ms_out1, ms_out2
    else
        println("Расход выходящего потока должен быть меньше входящего расхода")
    end
end