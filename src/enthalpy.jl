"""
calculates enthalpy of stream
at temperature T
"""
function mstream_HT(s::Stream, T)
    teststream = copy_stream(s)
    teststream.T = T
    if issingle(teststream)
        return myH_single(teststream)
    else
        return myH_double(teststream)
    end
end
