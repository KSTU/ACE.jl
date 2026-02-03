"""
normalize fractions to 1
"""
function frac_norm(A)
    Σ = sum(abs.(A))
    return abs.(A) / Σ
end


