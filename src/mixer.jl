"""
Функция для расчета смесителя двух потоков
Давление потоков должно быть одинаковым
Рассчитываются энтальпии для двух входящий потоков (h1, h2), находится общая энтальпия (h)
Рассчитывается общее количество молей в каждом потоке (A1, A2)
N - количество вещества в молях
"""
function mixer(ms_in1::MaterialStream, ms_in2::MaterialStream)
    if ms_in1.p == ms_in2.p
        h_in1 = mstream_H_T(ms_in1, ms_in1.T)   #энтельпия первого потока
        h_in2 = mstream_H_T(ms_in2, ms_in2.T)   #энтальпия второго потока
        h_out = (ms_in1.G * h_in1 + ms_in2.G * h_in2)/(ms_in1.G + ms_in2.G)    #удельная энтальпия выходящего потока
        A1 = @. ms_in1.Q * ms_in1.y + (1.0-ms_in1.Q) * ms_in1.x #суммарны доли компонентов в фазах
        A2 = @. ms_in2.Q * ms_in2.y + (1.0-ms_in2.Q) * ms_in2.x
        M_sm1 = sum(ms_in1.model.params.Mw .* A1) / 1000.0  #молярная масса смеси
        M_sm2 = sum(ms_in2.model.params.Mw .* A2) / 1000.0
        N1 = @. A1 * ms_in1.G / M_sm1     #количество молей
        N2 = @. A2 * ms_in2.G / M_sm2
        N = @. N1 + N2  #суммарное количество молей
        A = N / sum(N)  #мольная доля в двух потоках
        ms_out = mstream_TpA(ms_in1.G + ms_in2.G, ms_in1.T, ms_in1.p, A, ms_in1.model)    #определяется поток без заданной температуры
        Tinit = (ms_in1.T * ms_in1.G + ms_in2.T * ms_in2.G) / (ms_in1.G + ms_in2.G)  #начальное приближение для температуры
        ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - h_out, Tinit)   #определяем конечную температуру
        return ms_out
    else
        println("Давление входящих потоков должно быть одинаковым")
    end
end