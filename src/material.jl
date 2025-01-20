"""
Set stream at T -- temperature
Q -- vapor mole fraction
A -- initial mole fraction
"""
function mstream_TQA(G::Float64, T::Float64, A::Vector{Float64}, Q::Float64, model::EoSModel)
    CP = crit_mix(model, A)
    nsub = length(A)
    if T < CP[1]
        function nl_dew(y)
            eq = []
            for i = 1:nsub
                push!(eq, Q * y[i] + (1-Q) * dew_pressure(model, T, y)[4][i] - A[i])
            end
            return eq
        end
        y = nlsolve(nl_dew, A).zero
        p = dew_pressure(model, T, y)[1]
        x = dew_pressure(model, T, y)[4]
    else
        println("Сверхкритическое состояние")
    end
    return MaterialStream(G, T, p, x, y, Q, model)
end

"""
Set stream at T -- temperature
p -- pressure
A -- initial mole fraction
"""
function mstream_TpA(G::Float64, T::Float64, p::Float64, A::Vector{Float64}, model::EoSModel)
    CP = crit_mix(model, A)
    if T < CP[1]
        pkip = bubble_pressure(model, T, A)[1]
        pkon = dew_pressure(model, T, A)[1]
        if p <= pkon
            #одна фаза паровая
            Q = 1.0
            y = A
            x = zeros(Float64, length(A))
            println("Паровая фаза")
        elseif p >= pkip
            #одна фаза жидкая
            Q = 0.0
            x = A
            y = zeros(Float64, length(A))
            println("Жидкая фаза")
        else
            #двухфазная система
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
            x = dew_pressure(sys.model, T, y)[4]
            Q = sol[end]
        end
    else
        Q = 1.0
        y = A
        x = zeros(Float64, length(A))
    end
    return Stream(G, T, p, x, y, Q, model)
end

