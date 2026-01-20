

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
        "n_$(ro(n))"
    ]
    if ω == ωa push!(postfix_components, "o_2pi") else push!(postfix_components, "o_$(ro(ω))") end
    return join(postfix_components, "_")
end


"""
For calculation of steady state σ and Bα
"""
function get_postfix_steadyState(Δ, ΔvariDescription, dDescription, να, ηα, noPhonons, incField_wlf, tildeG_flags, arrayDescription, fiberPostfix, include3rdLevel, cDriveDescription, Δc, Ωc, cDriveArgs)
    postfix_components = [
        "D_$(ro(Δ, Int(4 + ceil(log10(ceil(abs(Δ + eps(Δ))))))))",
        ΔvariDescription,
        dDescription
    ]
    if noPhonons push!(postfix_components, "noPh") else push!(postfix_components, "tFr_$(join(ro.(να), ","))") end
    if all(ηα .== 0) push!(postfix_components, "LD_0.0") else push!(postfix_components, "LD_$(join(ro.(ηα), ","))") end
    push!(postfix_components, "wlf_$("[" * join(["(" * join([format_Complex_to_String(wlf[1]), wlf[2], wlf[3]], ",") * ")" for wlf in incField_wlf], ",") * "]")")
    if any(tildeG_flags .!= 1) push!(postfix_components, "tGfl_$(join(Int.(tildeG_flags), ","))") end
    append!(postfix_components, [
            arrayDescription,
            fiberPostfix,
    ])
    if include3rdLevel
        append!(postfix_components, [
          "cDr_$cDriveDescription", 
          "cDe_$Δc", 
          "cOm_$Ωc"
        ])
        if cDriveDescription == "plW" push!(postfix_components, "kc_$(join(ro.(cDriveArgs.kc), ","))") end
    end
    return "_" * join(postfix_components, "_")
end


"""
For calculation of transmission over classically disordered arrays
"""
function get_postfix_imperfectArray_transmission(Δ_specs, ΔvariDescription, dDescription, να, ηα, noPhonons, incField_wlf, n_inst, tildeG_flags, arrayDescription, fiberPostfix)
    postfix_components = [
        "D_$(join((ro(Δ_specs[1]), ro(Δ_specs[2]), Δ_specs[3]), ","))",
        ΔvariDescription,
        dDescription
    ]
    if noPhonons push!(postfix_components, "noPh") else push!(postfix_components, "tFr_$(join(ro.(να), ","))") end
    if all(ηα .== 0) push!(postfix_components, "LD_0.0") else push!(postfix_components, "LD_$(join(ro.(ηα), ","))") end
    append!(postfix_components, [
        "wlf_$("[" * join(["(" * join([format_Complex_to_String(wlf[1]), wlf[2], wlf[3]], ",") * ")" for wlf in incField_wlf], ",") * "]")",
        "nInst_$(n_inst)"
    ])
    if any(tildeG_flags .!= 1) push!(postfix_components, "tGfl_$(join(Int.(tildeG_flags), ","))") end
    append!(postfix_components, [
            arrayDescription,
            fiberPostfix,
    ])
    return join(postfix_components, "_")
end


"""
For time evolution
"""
function get_postfix_timeEvolution(Δ, ΔvariDescription, dDescription, να, ηα, noPhonons, incField_wlf, tildeG_flags, arrayDescription, fiberPostfix, initialStateDescription, tspan, dtmax, include3rdLevel, cDriveDescription, Δc, Ωc, cDriveArgs)
    postfix_components = [
        "D_$(ro(Δ, Int(4 + ceil(log10(ceil(abs(Δ + eps(Δ))))))))",
        ΔvariDescription,
        dDescription
    ]
    if noPhonons push!(postfix_components, "noPh") else push!(postfix_components, "tFr_$(join(ro.(να), ","))") end
    if all(ηα .== 0) push!(postfix_components, "LD_0.0") else push!(postfix_components, "LD_$(join(ro.(ηα), ","))") end
    push!(postfix_components, "wlf_$("[" * join(["(" * join([format_Complex_to_String(wlf[1]), wlf[2], wlf[3]], ",") * ")" for wlf in incField_wlf], ",") * "]")")
    if any(tildeG_flags .!= 1) push!(postfix_components, "tGfl_$(join(Int.(tildeG_flags), ","))") end
    append!(postfix_components, [
            arrayDescription,
            fiberPostfix,
            initialStateDescription,
            "t_$(join(ro.(tspan), ",")),$(ro(dtmax))"
    ])
    if include3rdLevel
        append!(postfix_components, [
          "cDr_$cDriveDescription", 
          "cDe_$Δc", 
          "cOm_$Ωc"
        ])
        if cDriveDescription == "plW" push!(postfix_components, "kc_$(join(ro.(cDriveArgs.kc), ","))") end
    end
    return join(postfix_components, "_")
end


"""
For memory efficiency 
"""
function get_postfix_memoryEfficiency(Δ, ΔvariDescription, dDescription, να, ηα, noPhonons, incField_wlf, tildeG_flags, arrayDescription, fiberPostfix, initialStateDescription, tspan, dtmax, radDecayRateAndStateNorm_LowerTol, cDriveDescription, Δc, Ωc, cDriveArgs)
    postfix_components = [
        "D_$(ro(Δ, Int(4 + ceil(log10(ceil(abs(Δ + eps(Δ))))))))",
        ΔvariDescription,
        dDescription
    ]
    if noPhonons push!(postfix_components, "noPh") else push!(postfix_components, "tFr_$(join(ro.(να), ","))") end
    if all(ηα .== 0) push!(postfix_components, "LD_0.0") else push!(postfix_components, "LD_$(join(ro.(ηα), ","))") end
    push!(postfix_components, "wlf_$("[" * join(["(" * join([format_Complex_to_String(wlf[1]), wlf[2], wlf[3]], ",") * ")" for wlf in incField_wlf], ",") * "]")")
    if any(tildeG_flags .!= 1) push!(postfix_components, "tGfl_$(join(Int.(tildeG_flags), ","))") end
    append!(postfix_components, [
        arrayDescription,
        fiberPostfix,
        initialStateDescription,
        "t_$(join(ro.(tspan), ",")),$(ro(dtmax))",
        "lTol_$(join(ro.(radDecayRateAndStateNorm_LowerTol), ","))",
        "cDr_$cDriveDescription", 
        "cDe_$Δc", 
        "cOm_$Ωc"
    ])
    if cDriveDescription == "plW" push!(postfix_components, "kc_$(join(ro.(cDriveArgs.kc), ","))") end
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


"""
Function to rename existing data files
"""
function rename()
    dataFolder = saveDir * "memoryEfficiency/"
    replacementPairs = ["__" => "_", ]
    
    for oldFilename in readdir(dataFolder)
        # newFilename = replace(oldFilename, replacementPairs...)
        # if newFilename != oldFilename mv(dataFolder * oldFilename, dataFolder * newFilename) end
        
        # if occursin("_LD_0.0_", oldFilename)
        #     indices = findfirst.(["tFr_", "_LD_0.0_"], oldFilename)
        #     startIndex = indices[1][1]
        #     endIndex   = indices[2][1] - 1
        #     replacementPairs = [oldFilename[startIndex:endIndex] => "noPh", ]
        #     newFilename = replace(oldFilename, replacementPairs...)
        #     mv(dataFolder * oldFilename, dataFolder * newFilename)
        # end
    end
end