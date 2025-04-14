

"""
For calculation of propagation constant
"""
function get_postfix(ω_specs::Tuple{Number, Number, Number}, ρf_specs::Tuple{Number, Number, Number}, n_specs::Tuple{Number, Number, Number})
    postfix_components = [
        "omega_$(join(ω_specs, ","))",
        "rhof_$(join(ρf_specs, ","))",
        "n_$(join(n_specs, ","))"
    ]
    return join(postfix_components, "_")
end


"""
For Fiber struct postfix
"""
function get_postfix(ρf::Number, n::Number, ω::Number)
    postfix_components = [
        "rhof_$ρf",
        "n_$n",
        "omega_$ω"
    ]
    return join(postfix_components, "_")
end


"""
For calculation of steady state σ and Bα
"""
function get_postfix(Δ, d, να, ηα, arrayDescription, fiberPostfix)
    if typeof(d) == String
        dipole_moment_string = "d_$d"
    else
        dipole_moment_string = "d_$(join(format_Complex_to_String.(d), ","))"
    end
    
    postfix_components = [
        "Delta_$(Δ)",
        dipole_moment_string,
        "trapFreqs_$(join(να, ","))",
        "LamDic_$(join(ηα, ","))",
        arrayDescription,
        fiberPostfix
    ]
    return join(postfix_components, "_")
end


"""
For time evolution
"""
function get_postfix(Δ, d, να, ηα, arrayDescription, fiberPostfix, initialStateDescription, tspan, dtmax)
    if typeof(d) == String
        dipole_moment_string = "d_$d"
    else
        dipole_moment_string = "d_$(join(format_Complex_to_String.(d), ","))"
    end
    
    postfix_components = [
        "Delta_$(Δ)",
        dipole_moment_string,
        "trapFreqs_$(join(να, ","))",
        "LamDic_$(join(ηα, ","))",
        arrayDescription,
        fiberPostfix,
        initialStateDescription,
        "tspan_$(join(tspan, ","))",
        "dtm_$dtmax"
    ]
    return join(postfix_components, "_")
end


"""
For calculation of the imaginary part of the transverse part of radiation mode Green's function or its derivatives 
"""
function get_postfix(ω, coords, derivOrder, α, fiberPostfix)
    postfix_components = [
        "omega_$ω",
        "coords_$(join(coords, ","))",
        "derivOrder_$(join(derivOrder, ","))",
        "alpha_$α",
        fiberPostfix
    ]
    return join(postfix_components, "_")
end


function save_as_jld2(data, saveDir, filename)
    JLD2.save(saveDir * filename * ".jld2", "data", data)
end


function load_as_jld2(saveDir, filename)
    return load(saveDir * filename * ".jld2", "data")
end


function save_as_txt(data, saveDir, filename)
    writedlm(saveDir * filename * ".txt", data, ' ')
end


function load_as_txt(saveDir, filename)
    return readdlm(saveDir * filename * ".txt", ' ', Float64, '\n')
end