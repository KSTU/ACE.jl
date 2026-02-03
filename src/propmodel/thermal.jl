mutable struct TermModelDef <: TermModel
    p::Float64
end

function term_empty()
    return TermModelDef(0.1)
end