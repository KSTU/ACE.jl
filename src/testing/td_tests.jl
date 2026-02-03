"""
Generate list of random substance
"""
function test_rand_subs(N::Int64)
    list = ["cycloheptane", "hexane", "octane", "methane", "cyclohexane"]
    return sample(list,N; replace=false)
end
"""
Generate random fractions
"""
function test_rand_frac(N::Int64)
    x = frac_norm(rand(N))
    while NaN in x
        x = frac_norm(rand(N))
    end
    return x
end

"""
Test pQA
"""
function test_pQA(Ntest::Int64, Nsub::Int64)
    time = zeros(Float64, Ntest)
    time2 = zeros(Float64, Ntest)
    for i = 1:Ntest
        sub = test_rand_subs(Nsub)
        frac = test_rand_frac(Nsub)
        println(sub)
        println(frac)
        time2[i] = @elapsed td_model = PR(sub, idealmodel = ReidIdeal)
        v_model = viscosity_empty()
        t_model = term_empty()
        model = Props(td_model, v_model, t_model)
        println(model)
        crit = crit_mix(model.td, frac)
        println("critical ", crit)
        cur_p = rand() * crit[2] * 0.8
        println("p ", cur_p)
        cur_Q = rand()
        println("Q ", cur_Q)
        time[i] = @elapsed mstream_pQA(1.0, cur_p, cur_Q, frac, model, "test1")
    end
    scatter(time .* 1000, time2 .* 1000)
    savefig("pQA.pdf")
    #t = @elapsed 
end
