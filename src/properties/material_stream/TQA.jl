"""
Set stream at T -- temperature
    Q -- vapor mole fraction
    A -- initial mole fraction
"""
function mstream_TQA(N::Float64, T::Float64, Q::Float64,
     A::Vector{Float64}, model::Props, name::String)
    cp = crit_mix(model.td, A)     #критическая точка
    nsub = length(A)
    if T < cp[1]    #если температура меньше критичсекой
        if Q == 1.0 #частный случай - все в паровой фазе
            y = A
            prop =dew_pressure(model.td, T, y)
            p = prop[1]
            x = zeros(Float64, nsub)
        elseif Q == 0.0     #все в жидкой фазе
            x = A
            prop = bubble_pressure(model.td, T, x)
            p = prop[1]
            y = zeros(Float64, nsub)
        else    #двухфазная система
            function nl_dew(y)
                prop = dew_pressure(model.td, T, y)
                eq = collect(Q * y[i] + (1.0-Q) * prop[4][i] - A[i] for i = 1:nsub)
                return eq
            end
            y = nlsolve(nl_dew, A).zero
            prop = dew_pressure(model.td, T, y)
            p = prop[1]
            x = collect(prop[4])
        end
        return MaterialStream(N, T, p, x, y, Q, model, name)
    else
        ace_error("Сверхкритическое состояние")
    end
end

function mstream_TQA(N::Float64, T::Float64, Q::Float64, A::Vector{Float64}, model::Props)
    return MaterialStream(N, T, Q, A, model, "none")
end