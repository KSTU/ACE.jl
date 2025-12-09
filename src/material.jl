"""
Set stream at T -- temperature
Q -- vapor mole fraction
A -- initial mole fraction
"""
function mstream_TQA(N::Float64, T::Float64, Q::Float64, A::Vector{Float64}, model::EoSModel)
    CP = crit_mix(model, A)     #критическая точка
    nsub = length(A)
    if T < CP[1]    #если температура меньше критичсекой
        if Q == 1.0 #частный случай - все в паровой фазе
            y = A
            prop =dew_pressure(model, T, y)
            p = prop[1]
            x = zeros(Float64, nsub)
        elseif Q == 0.0     #все в жидкой фазе
            x = A
            prop = bubble_pressure(model, T, x)
            p = prop[1]
            y = zeros(Float64, nsub)
        else    #двухфазная система
            function nl_dew(y)
                prop = dew_pressure(model, T, y)
                eq = collect(Q * y[i] + (1.0-Q) * prop[4][i] - A[i] for i = 1: nsub)
                return eq
            end
            y = nlsolve(nl_dew, A).zero
            prop =dew_pressure(model, T, y)
            p = prop[1]
            x = collect(prop[4])
        end
        return MaterialStream(N, T, p, x, y, Q, model)
    else
        println("Сверхкритическое состояние")
    end
end

"""
Set stream at p -- pressure
Q -- vapor mole fraction
A -- initial mole fraction
"""
function mstream_pQA(N::Float64, p::Float64,  Q::Float64, A::Vector{Float64}, model::EoSModel)
    CP = crit_mix(model, A)
    nsub = length(A)
    if Q == 1.0
        T = dew_temperature(model, p, A)[1]
        x = zeros(Float64, nsub)
        y = A
    elseif Q == 0.0        
        T = bubble_temperature(model, p, A)[1]
        x = A
        y = zeros(Float64, nsub)
    else
        function nl_dew(B)
            eq = []
            y = zeros(Float64, length(B))
            for i = 1:nsub-1
                y[i] = B[i]
            end
            y[nsub] = 1.0 - sum(y[1:nsub-1])
            T = B[nsub]
            dp = dew_pressure(model, T, y)
            for i = 1:nsub-1
                push!(eq, Q * y[i] + (1.0-Q) * dp[4][i] - A[i])
            end
            push!(eq, dp[1] - p)
            return eq
        end
        #начальное приближение
        guess = zeros(Float64, nsub)
        guess[1:nsub-1] = A[1:nsub-1]
        guess[nsub] = dew_temperature(model, p, A)[1]    #начальное приближенеи для смеси
        sol = nlsolve(nl_dew, guess).zero
        y = zeros(Float64, nsub)
        y[1:nsub-1] = sol[1:nsub-1]
        y[nsub] = 1.0 - sum(sol[1:nsub-1])
        T = sol[nsub]
        x = dew_pressure(model, T, y)[4]
    end
    return MaterialStream(N, T, p, x, y, Q, model)
end

"""
Set stream at T -- temperature
p -- pressure
A -- initial mole fraction
"""
function mstream_TpA(N::Float64, T::Float64, p::Float64, A::Vector{Float64}, model::EoSModel)
    CP = crit_mix(model, A)
    if T < CP[1]
        pkip = bubble_pressure(model, T, A)[1]
        pkon = dew_pressure(model, T, A)[1]
        if p <= pkon
            #одна фаза паровая
            Q = 1.0
            y = A
            x = zeros(Float64, length(A))
            #println("Паровая фаза")
        elseif p >= pkip
            #одна фаза жидкая
            Q = 0.0
            x = A
            y = zeros(Float64, length(A))
            # println("Жидкая фаза")
        else
            #двухфазная система
            nsub = length(A)
            function nl_sym(B)
                # B[1] - y[1]
                # B[2] - y[2]
                # B[end] - Q
                y = zeros(Float64, length(B))
                for i = 1:nsub-1
                    y[i] = B[i]
                end
                y[nsub] = 1.0 - sum(B[1:end-1])
                eq = []
                dp = dew_pressure(model, T, y)
                for i = 1:nsub-1
                    push!(eq, B[nsub] * y[i] + (1.0-B[nsub]) * dp[4][i] - A[i])
                end
                push!(eq, dp[1] - p)
                return eq
            end
            #initial guess
            guess = zeros(Float64, nsub)
            guess[1:end-1] = A[1:end-1]
            guess[end] = 0.4
            #solve system of equation
            sol = nlsolve(nl_sym, guess).zero
            y = zeros(Float64, nsub)
            y[1:end-1] = sol[1:end-1]
            y[end] = 1.0 - sum(sol[1:end-1])
            x = dew_pressure(model, T, y)[4]
            Q = sol[end]
        end
    else
        Q = 1.0
        y = A
        x = zeros(Float64, length(A))
    end
    return MaterialStream(N, T, p, x, y, Q, model)
end

function mstream_phA(N::Float64, p::Float64, h::Float64, A::Vector{Float64}, model::EoSModel)
    ms = mstream_pQA(N, p, 0.5, A, model)
    if length(A) > 1
        T = find_zero(T -> mstream_H_Tp(ms, T, p) - h, ms.T)
        return mstream_TpA(N, T, p, A, model)
    else
        Q = find_zero(Q -> mstream_H_pQ(ms, p, Q) - h, 0.5)
        return mstream_pQA(N, p, Q, A, model)
    end
end

function mstream_phA(N::Float64, p::Float64, h::Float64, A::Vector{Float64}, model::EoSModel, Tinit::Float64)
    ms = mstream_pQA(N, p, 1.0, A, model)
    T = find_zero(T -> mstream_H_Tp(ms, T, p) - h, Tinit)
    return mstream_TpA(N, T, p, A, model)
end
