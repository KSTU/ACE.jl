"""
Set material stream at p - pressure
h - enthalpy
A - initial mole fraction
"""
function mstream_phA(N::Float64, p::Float64, h::Float64, A::Vector{Float64}, model::EoSModel, name::String)
    ms = mstream_pQA(N, p, 0.5, A, model, name)
    if length(A) > 1
        T = find_zero(T -> mstream_H_Tp(ms, T, p) - h, ms.T)
        return mstream_TpA(N, T, p, A, model, name)
    else
        Q = find_zero(Q -> mstream_H_pQ(ms, p, Q) - h, 0.5)
        return mstream_pQA(N, p, Q, A, model, name)
    end
end

function mstream_phA(N::Float64, p::Float64, h::Float64, A::Vector{Float64},
     model::EoSModel, name::String)
    return mstream_phA(N, p, h, A, model, "none")
end

"""
Set material stream at p - pressure
h - enthalpy
A - initial mole fraction

Tinit - initial temperature
"""
function mstream_phA(N::Float64, p::Float64, h::Float64, A::Vector{Float64}, model::EoSModel, Tinit::Float64)
    ms = mstream_pQA(N, p, 1.0, A, model)
    T = find_zero(T -> mstream_H_Tp(ms, T, p) - h, Tinit)
    return mstream_TpA(N, T, p, A, model)
end


function ms_support_D_p(model, p, T, A)
    ρ = 0
end

"""
Set material stream at 
    D - density kg/m^3
    T - temperature
    A - initial mole fraction
"""
function mstream_DTA(N::Float64, D::Float64, T::Float64, A::Vector{Float64}, model::EoSModel, name::String)
    pᵥ, _, _, _ = bubble_pressure(model, T, A)
    pₗ, _, _, _ = dew_pressure(model,T, A)
    println("p v ", pᵥ)
    ρᵥ = mass_density(model, pᵥ, T, A; phase = :vapor)
    ρₗ = mass_density(model, pₗ, T, A; phase = :liquid)
    println("ρ v ", ρᵥ)
    println("ρ l ", ρₗ)
    if ρᵥ < D < ρₗ
        function nl_sum(B)
            x = similar(B)

        end
        println("test")
    elseif D <= ρᵥ
        p = find_zero(p -> mass_density(model, p, T, A; phase = :vapor) - D, pᵥ)
        Q = 1.0
        y = A
        x = zeros(length(A))
    elseif D >= ρₗ
        p = find_zero(p -> mass_density(model, p, T, A; phase = :liquid) - D, pᵥ)
        Q = 0.0
        y = zeros(length(A))
        x = A
    end
    println( "p ",  p)
    return MaterialStream(N, T, p, x, y, Q, model, name)
end