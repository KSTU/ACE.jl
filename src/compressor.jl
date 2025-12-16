"""
фуникция что то делает

"""
function compress(ms, Δp)
    ms_out = material_stream_copy(ms)
    hᵢ = mstream_H(ms_out)
    N = 100
    pₒᵤₜ = ms.p + Δp
    δp = Δp / N
    A =  @. ms_out.Q * ms_out.y + (1.0-ms_out.Q) * ms_out.x         #расчет суммарной доли по фазам

    # println(hᵢ)
    # for i = 1:N
    #     ρ = mstream_Dm(ms_out)
    #     hᵢ = hᵢ + δp / ρ
    #     println("ρm ", ρ, " h ", hᵢ, " ρ ", mstream_D(ms_out), " δp / ρ ", δp / ρ)
    #     ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ, ms_out.T)
    #     ms_out.p = ms_out.p + δp
    #     ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
    #     println("")
    # end

    # #1) сверхкритика
    # #2) если пар и давление <
    # #3) если жидкость и давление >

    CP = crit_mix(ms_out.model, A)

    if ms_out.T > CP[1] || (ms_out.Q == 1 && Δp < 0) || (ms_out.Q == 0 && Δp > 0)
        for i = 1:N
            ρ = mstream_Dm(ms_out)
            hᵢ = hᵢ + δp / ρ        
            ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ, ms_out.T)
            ms_out.p = ms_out.p + δp
            ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
        end
    elseif ms_out.Q == 1 && Δp > 0.0
        #1 интеграл
        #println(ms_out.p)
        for i = 1:N
            pₖ = dew_pressure(ms_out.model, ms_out.T, ms_out.y)[1]
            if ms_out.p + δp < pₖ
                ρ = mstream_Dm(ms_out)
                hᵢ = hᵢ + δp / ρ
                ms_out.p = ms_out.p + δp
                ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ, ms_out.T)
                ms_out = mstream_TpA(ms_out.N, ms_out.T, ms_out.p, A, ms_out.model)
            else
                ρ = mstream_Dm(ms_out)
                #hᵢ = hᵢ + (pₖ - ms_out.p) / ρ
                ms_out.p = pₖ
                #ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ, ms_out.T)
                ms_out = mstream_pQA(ms_out.N, pₖ, 1.0, A, ms_out.model)
                hᵢ = mstream_H(ms_out)
                break
            end
        end
        #println("p int 1 ", ms_out.p, " T ", ms_out.T)
        #2 интеграл
        if ms_out.p < pₒᵤₜ
            for q = 1:-0.01:0
                ρ = mstream_Dm(ms_out)
                pᵢ₊₁ = mstream_TQA(ms_out.N, ms_out.T, q, A, ms_out.model).p
                hᵢ = hᵢ + (pᵢ₊₁ - ms_out.p) / ρ
                println("H ", mstream_H_pQ(ms_out, pᵢ₊₁, q), " p ", pᵢ₊₁, " q ", q, "ms out Q ", ms_out.Q, " p ", pₒᵤₜ, " T ", ms_out.T)
                if pᵢ₊₁ > pₒᵤₜ
                    println(" err ", mstream_H_pQ(ms_out, pₒᵤₜ, q), " ms out Q", ms_out.Q)
                    q = find_zero(q -> mstream_H_pQ(ms_out, pₒᵤₜ, q) - hᵢ, q)
                    ms_out = mstream_pQA(ms_out.N, pₒᵤₜ, q, A, ms_out.model)
                    break
                else
                    T = find_zero(T -> mstream_H_TQ(ms_out, T, q) - hᵢ, ms_out.T)
                    #ms_out = mstream_pQA(ms_out.N, pᵢ₊₁, q, A, ms_out.model)
                    #ms_out = mstream_TQA(ms_out.N, T, q, A, ms_out.model)
                    ms_out = mstream_TpA(ms_out.N, T, pᵢ₊₁, A, ms_out.model)
                end
            end
        end
        #третий интеграл
        if ms_out.p < pₒᵤₜ
            δp₃ = (pₒᵤₜ - ms_out.p) / N
            for i = 1:N
                ρ = mstream_Dm(ms_out)
                hᵢ = hᵢ + δp₃ / ρ
                ms_out.T = find_zero(T -> mstream_H_T(ms_out, T) - hᵢ, ms_out.T)
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