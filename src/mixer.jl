"""
Функция для расчета смесителя двух потоков
Рассчитываются энтальпии для двух входящий потоков (h1, h2), находится общая энтальпия (h)
Рассчитывается общее количество молей в каждом потоке (A1, A2)
N - количество вещества в молях
"""
function mixer(sin1::Stream, sin2::Stream)
    sout = copy_stream(sin1)
    m_sm = sum(sin1.sys.model.params.Mw .* sin1.x) / 1000.0
    hl1 = enthalpy(sin1.sys.model, sin1.p, sin1.T, sin1.x) / m_sm
    m_sm = sum(sin1.sys.model.params.Mw .* sin1.y) / 1000.0
    hv1 = enthalpy(sin1.sys.model, sin1.p, sin1.T, sin1.y) / m_sm
    h1 = hl1+hv1
    m_sm = sum(sin2.sys.model.params.Mw .* sin2.x) / 1000.0
    hl2 = enthalpy(sin2.sys.model, sin2.p, sin2.T, sin2.x) / m_sm
    m_sm = sum(sin2.sys.model.params.Mw .* sin2.y) / 1000.0
    hv2 = enthalpy(sin2.sys.model, sin2.p, sin2.T, sin2.y) / m_sm
    h2 = hl2+hv2
    h = (sin1.G * h1 + sin2.G * h2)/(sin1.G + sin2.G)
    sout.G = sin1.G + sin2.G
    A1 = @. sin1.Q * sin1.y + (1-sin1.Q) * sin1.x
    A2 = @. sin2.Q * sin2.y + (1-sin2.Q) * sin2.x
    m_sm1 = sum(sin1.sys.model.params.Mw .* A1) / 1000.0
    m_sm2 = sum(sin2.sys.model.params.Mw .* A2) / 1000.0
    N1 = @. A1 * sin1.G / m_sm1
    N2 = @. A2 * sin2.G / m_sm2
    N = @. N1 + N2
    A = N / sum(N)
    sout.T = find_zero(T -> myH(sout, T) - h, sout.T)
    return sout
end