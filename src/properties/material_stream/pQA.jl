"""
Set stream at p -- pressure
    Q -- vapor mole fraction
    A -- initial mole fraction
"""
function mstream_pQA(N::Float64, p::Float64,  Q::Float64, A::Vector{Float64},
    model::Props, name::String)
    cp = crit_mix(model.td, A) #critial point
    nsub = length(A)
    if p < cp[2]
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
                y = zeros(Float64, length(B))
                for i = 1:nsub-1
                    y[i] = B[i]
                end
                y[nsub] = 1.0 - sum(y[1:nsub-1])
                y = frac_norm(y)
                T = B[nsub]
                dp = dew_pressure(model.td, T, y)
                eq = zeros(Float64, nsub)
                for i = 1:nsub-1
                    eq[i] = Q * y[i] + (1.0-Q) * dp[4][i] - A[i]
                end
                eq[nsub] = dp[1] - p
                println("y ", y)
                println("p ", p)
                println("eq ", eq)
                return eq
            end
            #начальное приближение
            guess = zeros(Float64, nsub)
            guess[1:nsub-1] = A[1:nsub-1]
            guess[nsub] = dew_temperature(model.td, p, A)[1]    #начальное приближенеи для смеси
            sol = nlsolve(nl_dew, guess).zero
            y = zeros(Float64, nsub)
            y[1:nsub-1] = sol[1:nsub-1]
            y[nsub] = 1.0 - sum(sol[1:nsub-1])
            y = frac_norm(y)
            T = sol[nsub]
            x = dew_pressure(model.td, T, y)[4]
        end
    else
        ace_error("Сверхкритическое состояние")
    end
    return MaterialStream(N, T, p, x, y, Q, model, name)
end

function mstream_pQA(N::Float64, p::Float64,  Q::Float64, A::Vector{Float64}, model::Props)
    return mstream_pQA(N, p, Q, A, model, "none")
end
