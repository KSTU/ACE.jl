"""
Функция для расчета смесителя двух потоков
Давление потоков должно быть одинаковым
Рассчитываются энтальпии для двух входящий потоков (h1, h2), находится общая энтальпия (h)
Рассчитывается общее количество молей в каждом потоке (A1, A2)
N - количество вещества в молях
"""
function mixer(ms_in1::MaterialStream, ms_in2::MaterialStream)
    if ms_in1.p == ms_in2.p
        a1 = @. ms_in1.Q * ms_in1.y + (1.0-ms_in1.Q) * ms_in1.x #суммарны доли компонентов в фазах
        a2 = @. ms_in2.Q * ms_in2.y + (1.0-ms_in2.Q) * ms_in2.x
        aout = @. (a1 * ms_in1.N + a2 * ms_in2.N) / (ms_in1.N + ms_in2.N)

        h_in1 = mstream_H_T(ms_in1, ms_in1.T)   #удельная энтельпия первого потока
        h_in2 = mstream_H_T(ms_in2, ms_in2.T)   #удельная энтальпия второго потока
        h_out = (ms_in1.N * h_in1 + ms_in2.N * h_in2)/(ms_in1.N + ms_in2.N)    #удельная энтальпия выходящего потока

        Tinit = (ms_in1.T * ms_in1.N + ms_in2.T * ms_in2.N) / (ms_in1.N + ms_in2.N)  #начальное приближение для температуры
        ms_out = mstream_TpA(ms_in1.N + ms_in2.N, ms_in1.T, ms_in1.p, aout, ms_in1.model)    #определяется поток без заданной температуры
        Troot = find_zero(T -> mstream_H_T(ms_out, T) - h_out, Tinit)   #определяем конечную температуру
        return mstream_TpA(ms_in1.N + ms_in2.N, Troot, ms_in1.p, aout, ms_in1.model)    #подставляем найденное значение и определяем состояние выходящего потока
    else
        println("Давление входящих потоков должно быть одинаковым")
    end
end