

"""
For calculation of propagation constant
"""
function get_postfix_scan_propConst(ω_specs, ρf_specs, n_specs)
    postfix_components = [
        "omega_$(join((ro(ω_specs[1]), ro(ω_specs[2]), ω_specs[3]), ","))",
        "rhof_$(join((ro(ρf_specs[1]), ro(ρf_specs[2]), ρf_specs[3]), ","))",
        "n_$(join((ro(n_specs[1]), ro(n_specs[2]), n_specs[3]), ","))"
    ]
    return join(postfix_components, "_")
end


"""
For Fiber struct postfix
"""
function get_postfix_Fiber(ρf, n, ω)
    postfix_components = [
        "rhof_$(ro(ρf))",
        "n_$(ro(n))",
        "omega_$(ro(ω))"
    ]
    return join(postfix_components, "_")
end


"""
For calculation of steady state σ and Bα
"""
function get_postfix_steadyState(Δ, ΔvariDescription, dDescription, να, ηα, incField_wlf, arrayDescription, fiberPostfix)
    postfix_components = [
        "Delta_$(ro(Δ))",
        ΔvariDescription,
        dDescription,
        "trapFreqs_$(join(ro.(να), ","))",
        "LamDic_$(join(ro.(ηα), ","))",
        "wlf_$("[" * join(["(" * join([format_Complex_to_String(wlf[1]), wlf[2], wlf[3]], ",") * ")" for wlf in incField_wlf], ",") * "]")",
        arrayDescription,
        fiberPostfix
    ]
    return join(postfix_components, "_")
end


"""
For calculation of transmission over classically disordered arrays
"""
function get_postfix_classDisorder_transmission(Δ_specs, ΔvariDescription, dDescription, να, ηα, incField_wlf, n_inst, arrayDescription, fiberPostfix)
    postfix_components = [
        "Delta_$(join((ro(Δ_specs[1]), ro(Δ_specs[2]), Δ_specs[3]), ","))",
        ΔvariDescription,
        dDescription,
        "trapFreqs_$(join(ro.(να), ","))",
        "LamDic_$(join(ro.(ηα), ","))",
        "wlf_$("[" * join(["(" * join([format_Complex_to_String(wlf[1]), wlf[2], wlf[3]], ",") * ")" for wlf in incField_wlf], ",") * "]")",
        "nInst_$(n_inst)",
        arrayDescription,
        fiberPostfix
    ]
    return join(postfix_components, "_")
end


"""
For time evolution
"""
function get_postfix_timeEvolution(Δ, ΔvariDescription, dDescription, να, ηα, incField_wlf, arrayDescription, fiberPostfix, initialStateDescription, tspan, dtmax)
    postfix_components = [
        "Delta_$(ro(Δ))",
        ΔvariDescription,
        dDescription,
        "trapFreqs_$(join(ro.(να), ","))",
        "LamDic_$(join(ro.(ηα), ","))",
        "wlf_$("[" * join(["(" * join([format_Complex_to_String(wlf[1]), wlf[2], wlf[3]], ",") * ")" for wlf in incField_wlf], ",") * "]")",
        arrayDescription,
        fiberPostfix,
        initialStateDescription,
        "tspan_$(join(ro.(tspan), ","))",
        "dtm_$(ro(dtmax))"
    ]
    return join(postfix_components, "_")
end


"""
For calculation of the imaginary part of the transverse part of radiation mode Green's function or its derivatives 
"""
function get_postfix_Im_Grm_trans(ω, coords, derivOrder, α, abstol, fiberPostfix)
    postfix_components = [
        "omega_$(ro(ω))",
        "coords_$(join(ro.(coords), ","))",
        "derivOrder_$(join(derivOrder, ","))",
        "alpha_$α",
        "abstol_$abstol",
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