"""
Mix 2 stream
Pressure of 2 streams must be equal
    ms_in₁ - 1 stream
    ms_in₂ - 2 stream
"""
function app_mixer(ms_in₁::MaterialStream, ms_in₂::MaterialStream)
    if ms_in₁.p == ms_in₂.p
        a₁ = @. ms_in₁.Q * ms_in₁.y + (1.0-ms_in₁.Q) * ms_in₁.x #суммарны доли компонентов в фазах
        a₂ = @. ms_in₂.Q * ms_in₂.y + (1.0-ms_in₂.Q) * ms_in₂.x
        aout = @. (a₁ * ms_in₁.N + a₂ * ms_in₂.N) / (ms_in₁.N + ms_in₂.N)

        h_in₁ = mstream_H_T(ms_in₁, ms_in₁.T)   #удельная энтельпия первого потока
        h_in2 = mstream_H_T(ms_in₂, ms_in₂.T)   #удельная энтальпия второго потока
        h_out = (ms_in₁.N * h_in₁ + ms_in₂.N * h_in2)/(ms_in₁.N + ms_in₂.N)    #удельная энтальпия выходящего потока

        Tinit = (ms_in₁.T * ms_in₁.N + ms_in₂.T * ms_in₂.N) / (ms_in₁.N + ms_in₂.N)  #начальное приближение для температуры
        ms_out = mstream_TpA(ms_in₁.N + ms_in₂.N, ms_in₁.T, ms_in₁.p, aout, ms_in₁.model)    #определяется поток без заданной температуры
        Troot = find_zero(T -> mstream_H_T(ms_out, T) - h_out, Tinit)   #определяем конечную температуру
        return mstream_TpA(ms_in₁.N + ms_in₂.N, Troot, ms_in₁.p, aout, ms_in₁.model)    #подставляем найденное значение и определяем состояние выходящего потока
    else
        ace_error("Давление входящих потоков должно быть одинаковым")
    end
end