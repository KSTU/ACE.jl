"""
Split material stream for 2
    frac - molar fraction of 1 strem
"""
function app_splitter_frac(ms_in::MaterialStream, frac::Float64)
    if frac >= 0 && frac <= 1
        ms_out₁ = material_stream_copy(ms_in)
        ms_out₁.N = ms_in.N * frac
        ms_out₂ = material_stream_copy(ms_in)
        ms_out₂.N = ms_in.N - ms_out₁.N
    return ms_out₁, ms_out₂
    else 
        ace_error("The fraction must be in the range from 0 to 1")
    end
end

"""
Split material stream for 2
    N₁ - molar flow of 1 stream
"""
function app_splitter_N(ms_in::MaterialStream, N₁::Float64)
    if 0.0 <= N₁ <= ms_in.N 
        ms_out₁ = material_stream_copy(ms_in)
        ms_out₁.N = N₁
        ms_out₂ = material_stream_copy(ms_in)
        ms_out₂.N = ms_in.N - ms_out₁.N
        return ms_out₁, ms_out₂
    else
        ace_error("The outlet flow rate must be lower than the inlet flow rate.")
    end
end