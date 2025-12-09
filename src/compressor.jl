"""
фуникция что то делает

"""
function compress(ms, Δp)
    ms_out = material_stream_copy(ms)
    hᵢ = mstream_H(ms_out)
    N = 1000
    pₒᵤₜ = ms.p + Δp
    δp = Δp / N
    A =  @. ms_out.Q * ms_out.y + (1.0-ms_out.Q) * ms_out.x         #расчет суммарной доли по фазам

    #1) сверхкритика
    #2) если пар и давление <
    #3) если жидкость и давление >

    CP = crit_mix(ms_out.model, A)

    if ms_out.T > CP[1] || (ms_out.Q == 1 && Δp < 0) || (ms_out.Q == 0 && Δp > 0)
        for i = 1:N
            ρ = mstream_D(ms_out)
            hᵢ₊₁ = hᵢ + δp / ρ        
            ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
            ms_out.p = ms_out.p + δp
            ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
        end
    elseif ms_out.Q == 1 && Δp > 0.0
        #1 интеграл
        println(ms_out.p)
        for i = 1:N
            pₖ = dew_pressure(ms_out.model, ms_out.T, ms_out.y)[1]
            if ms_out.p + δp < pₖ
                ρ = mstream_Dm(ms_out)
                hᵢ₊₁ = hᵢ + δp / ρ * 1000
                ms_out.p = ms_out.p + δp
                ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
                ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
            else
                ρ = mstream_Dm(ms_out)
                hᵢ₊₁ = hᵢ + (pₖ - ms_out.p) / ρ
                ms_out.p = ms_out.p + (pₖ - ms_out.p)
                ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
                ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
                break
            end
        end
        println("p int 1 ", ms_out.p, " T ", ms_out.T)
            #2 интеграл
        if ms_out.p < pₒᵤₜ
            for q = 1:-0.05:0
                ρ = mstream_Dm(ms_out)
                pᵢ₊₁ = mstream_TQA(ms_out.N, ms_out.T, q, A, ms_out.model).p
                hᵢ₊₁ = hᵢ + (pᵢ₊₁ - ms_out.p) / ρ
                if pᵢ₊₁ > pₒᵤₜ
                    q = find_zero(q -> mstream_H_pQ(ms.out, pₒᵤₜ, q) - hᵢ₊₁, ms_out.q)
                    ms_out = mstream_pQA(ms_out.N, pₒᵤₜ, q, A, ms_out.model)
                    break
                else
                    T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
                    ms_out = mstream_TpA(ms_out.N, T, pᵢ₊₁, A, ms_out.model)
                end
            end
        end
        #третий интеграл
        if ms_out.p < pₒᵤₜ
            δp₃ = (pₒᵤₜ - ms_out.p) / N
            for i = 1:N
                ρ = mstream_Dm(ms_out)
                hᵢ₊₁ = hᵢ + δp₃ / ρ
                ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
                ms_out.p = ms_out.p + δp₃
                ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
            end
        end
    elseif ms_out.Q == 0 && Δp < 0.0
        #1 интеграл
        for i = 1:N
            pₖ = bubble_pressure(ms_out.model, ms_out.T, ms_out.x)[1]
            if ms_out.p + δp > pₖ
                ρ = mstream_Dm(ms_out)
                hᵢ₊₁ = hᵢ + δp / ρ
                ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
                ms_out.p = ms_out.p + δp
                ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
            else
                ρ = mstream_Dm(ms_out)
                hᵢ₊₁ = hᵢ + (pₖ - ms_out.p) / ρ
                ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
                ms_out.p = ms_out.p + (pₖ - ms_out.p)
                ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
                break
            end
        end
            #2 интеграл
        if ms_out.p > pₒᵤₜ
            for q = 0:0.05:1
                ρ = mstream_Dm(ms_out)
                pᵢ₊₁ = mstream_TQA(ms_out.N, ms_out.T, q, A, ms_out.model).p
                hᵢ₊₁ = hᵢ + (pᵢ₊₁ - ms_out.p) / ρ
                if pᵢ₊₁ < pₒᵤₜ
                    q = find_zero(q -> mstream_H_pQ(ms.out, pₒᵤₜ, q) - hᵢ₊₁, ms_out.q)
                    ms_out = mstream_pQA(ms_out.N, pₒᵤₜ, q, A, ms_out.model)
                    break
                else
                    T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
                    ms_out = mstream_TpA(ms_out.N, T, pᵢ₊₁, A, ms_out.model)
                end
            end
        end
        #третий интеграл
        if ms_out.p > pₒᵤₜ
            δp₃ = (pₒᵤₜ - ms_out.p) / N
            for i = 1:N
                ρ = mstream_Dm(ms_out)
                hᵢ₊₁ = hᵢ + δp₃ / ρ
                ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ₊₁, ms_out.T)
                ms_out.p = ms_out.p + δp₃
                ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
            end
        end
    end
    
    return ms_out
end