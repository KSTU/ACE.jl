
"""
Рассчитывается колонна
R - флегмовое число
Nw - расход куба
p - давление
Δp - перепад давления
N - количество тарелок
Nf - номер тарелки питания
Ey - кпд по Мерфри
...
"""
function rect_R_Nw(ms::MaterialStream, R::Float64, Nw::Float64, p::Float64, Δp::Float64, N::Int64, Nf::Int64, Ey::Array{Float64})
    K = length(ms.x)    #количество компонентов

    N = N+2
    Nd = ms.N - Nw
    Nr = Nd * R

    pₜ = collect(p + Δp * (i-1) / (N-1) for i = 1:N) #давление на тарелке

    println("ND ", Nd)
    println("Nr ", Nr)
    println("p ", pₜ)

    function f(B)
        Nₗ = B[1:N] #мольный расход жидкой фазы на тарелке
        Nᵥ = B[N+1:2*N] #мольный расход паровой фазы на тарелке
        x = zeros(Float64, N, K)
        y = zeros(Float64, N, K)
        for i = 1:N
            for j = 1:K-1
                x[i,j] = B[2*N+(i-1)*(K-1)+j]
                y[i,j] = B[2*N+(N*(K-1))+(i-1)*(K-1)+j]
            end
            x[i,K] = 1.0 - sum(x[i,1:K-1])
            y[i,K] = 1.0 - sum(y[i,1:K-1])
        end
        println("Nl ", Nₗ)
        println("Nv ", Nᵥ)
        println("x ", x)
        println("y ", y)

        eq = zeros(Float64, length(B))
        for i = 1:N
            if i == 1
                eq[i] = Nᵥ[2] - Nd - Nr
            elseif i == N
                eq[N] = Nₗ[N-1] - Nᵥ[N] - Nw
            elseif i == Nf
                eq[i] = ms.N + Nₗ[i-1] + Nᵥ[i+1] - Nₗ[i] - Nᵥ[i]
            else
                eq[i] = Nₗ[i-1] + Nᵥ[i+1] - Nₗ[i] - Nᵥ[i]
            end
        end
        A =  @. ms.Q * ms.y + (1.0-ms.Q) * ms.x
        ndx = N
        for i = 1:N
            for j = 1:K-1
                ndx += 1
                if i == 1
                    eq[ndx] = Nᵥ[2] * y[2,j] - Nd * x[1,j] - Nr * x[1,j]
                elseif i == N
                    eq[ndx] = Nₗ[N-1] * x[N-1,j] - Nᵥ[N] * y[N,j] - Nw * x[N,j]
                elseif i == Nf
                    eq[ndx] = ms.N * A[j] + Nₗ[i-1] * x[i-1,j] + Nᵥ[i+1] * y[i+1,j] - Nₗ[i] * x[i,j] - Nᵥ[i] * y[i,j]
                else
                    eq[ndx] = Nₗ[i-1] * x[i-1,j] + Nᵥ[i+1] * y[i+1,j] - Nₗ[i] * x[i,j] - Nᵥ[i] * y[i,j]
                end
            end
        end
        Tₜ = zeros(Float64, N)
        for i = 1:N
            println("x ", x[i,:])
            Tₜ[i], _, _, yr = bubble_temperature(ms.model, pₜ[i], x[i,:])
            if i == N    #условие в ребойлере
                for j = 1:K-1
                    ndx += 1
                    eq[ndx] = y[i,j] - yr[j]
                end
            else    #не последняя тарелка с учетом эффективности
                for j = 1:K-1
                    ndx += 1
                    eq[ndx] = Ey[i] * (yr[j] - y[i+1,j]) - (y[i,j] - y[i+1,j])
                end
            end
        end
        Hₗ = collect(enthalpy(ms.model, pₜ[i], Tₜ[i], Nₗ[i] .* x[i,:], phase =:liquid) for i = 1:N) #энтальпия жидкой фазы
        Hᵥ = collect(enthalpy(ms.model, pₜ[i], Tₜ[i], Nᵥ[i] .* y[i,:], phase =:vapor) for i = 1:N)
        Hf = mstream_H(ms) * ms.N
        Hd = enthalpy(ms.model, pₜ[1], Tₜ[1], Nd .* x[1,:], phase =:liquid)
        Hr = enthalpy(ms.model, pₜ[1], Tₜ[1], Nr .* x[1,:], phase =:liquid)
        Hw = enthalpy(ms.model, pₜ[N], Tₜ[N], Nw .* x[N,:], phase =:liquid)
        
        for i = 1:N    #тепловой баланс
            ndx += 1
            if i == 1
                eq[ndx] = Nᵥ[2] - Nd - Nr
            elseif i == N
                eq[ndx] = Nₗ[N-1] - Nᵥ[N] - Nw
            elseif i == Nf
                eq[ndx] = ms.N + Nₗ[i-1] + Nᵥ[i+1] - Nₗ[i] - Nᵥ[i]
            else
                eq[ndx] = Nₗ[i-1] + Nᵥ[i+1] - Nₗ[i] - Nᵥ[i]
            end
        end

        return eq
    end

    Nₗ = fill(0.3, N)
    Nᵥ = fill(0.4, N)
    x = fill(0.5, N, K-1)
    y = fill(0.5, N, K-1)

    B = zeros(Float64, 2*N*K)


    B[1:N] = Nₗ
    B[N+1:2*N] = Nᵥ
    for i = 1:N
        for j = 1:K-1
            # println("i ", i, " j ", j)
            # println("id x ", 2*N+(i-1)*(K-1)+j)
            # println("id y ", 2*N+(N*(K-1))+(i-1)*(K-1)+j)
            B[2*N+(i-1)*(K-1)+j] = x[i,j]
            B[2*N+(N*(K-1))+(i-1)*(K-1)+j] = y[i,j]
        end
    end

    println("B test", B)

    println("f(b) ", f(B))

end