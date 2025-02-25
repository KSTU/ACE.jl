"""
Делит поток на 2 части где задается массовая доля первого потока
"""
function splitter_dol(ms_in::MaterialStream, dol::Float64)
    if dol >= 0 && dol <= 1
        ms_out1 = material_stream_copy(ms_in)
        ms_out1.G = ms_in.G * dol
        ms_out2 = material_stream_copy(ms_in)
        ms_out2.G = ms_in.G - ms_out1.G
    return ms_out1, ms_out2
    else 
        println("Доля должна быть в диапазоне от 0 до 1")
    end
end

"""
Делит поток на 2 части где расход первого потока
"""
function splitter_G(ms_in::MaterialStream, G1::Float64)
    if G1 <= ms_in.G 
        ms_out1 = material_stream_copy(ms_in)
        ms_out1.G = G1
        ms_out2 = material_stream_copy(ms_in)
        ms_out2.G = ms_in.G - ms_out1.G
        return ms_out1, ms_out2
    else
        println("Расход выходящего потока должен быть меньше входящего расхода")
    end
end