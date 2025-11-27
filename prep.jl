

"""
Prepare the atomic coherences (in the case of no phonons)
"""
function prep_times_σTrajectories(xTrajectories, N)
    times = xTrajectories[:, 1]
    
    σTrajectories_t = unpack_σFromx.(eachrow(xTrajectories[:, 2:end]))
    
    σTrajectories = [[σTrajectories_t[t][i] for t in eachindex(times)] for i in 1:N]
    return times, σTrajectories
end


"""
Prepare the atomic coherences and the atom-phonon correlations 
"""
function prep_times_σBαTrajectories(xTrajectories, N)
    times = xTrajectories[:, 1]
    
    σBαTrajectories = unpack_σBαFromx.(eachrow(xTrajectories[:, 2:end]))
    
    σTrajectories  =  [[σBαTrajectories[t][1][i]    for t in eachindex(times)] for i in 1:N]
    BαTrajectories = [[[σBαTrajectories[t][2][α][i] for t in eachindex(times)] for i in 1:N^2] for α in 1:3]
    return times, σTrajectories, BαTrajectories
end



"""
Prepare the squared magnitude and phase of some vector
"""
function prep_squaredNorm_phase(t)
    return abs2.(t), angle.(t)
end



"""
Prepare the squared magnitude, phase, unwrapped phase, phase per atom, and slope of phase
with respect to Δ_range of some vector
"""
function prep_squaredNorm_phase_unwrappedPhase_phasePerAtom_phaseSlope(SP, t)
    absSquared = abs2.(t)
    phase = angle.(t)
    unwrappedPhase = vcat([angle(t[1])], angle(t[1]) .+ cumsum(angle.(t[2:end]./t[1:end-1])))
    phasePerAtom = unwrappedPhase/SP.N
    phaseSlope = diff(unwrappedPhase)./diff(SP.Δ_range)   
    return absSquared, phase, unwrappedPhase, phasePerAtom, phaseSlope
end


"""
Prepare the means and standard deviations of the transmission as calculated over
a set of imperfect atomic arrays
"""
function prep_imperfectArray_transmission(ts)
    ts_mat = vectorOfRows2Matrix(ts)
    return squeeze(mean(real.(ts_mat), dims=1)),
           squeeze( std(real.(ts_mat), dims=1)), 
           squeeze(mean(imag.(ts_mat), dims=1)),
           squeeze( std(imag.(ts_mat), dims=1))
    
    # Ts     = [abs2.(t) for t in ts]
    # phases = [angle.(t) for t in ts]
    # # phases = [vcat([angle(t[1])], angle(t[1]) .+ cumsum(angle.(t[2:end]./t[1:end-1]))) for t in ts]
    
    # T_mat      = vectorOfRows2Matrix(Ts)
    # phases_mat = vectorOfRows2Matrix(phases)
    # return           squeeze(mean(T_mat, dims=1)), 
    #                  squeeze( std(T_mat, dims=1)), 
    #                  squeeze(mean(phases_mat, dims=1)), 
    #     #    wrapPhase(squeeze(mean(phases_mat, dims=1))), 
    #                  squeeze( std(phases_mat, dims=1))
end


"""
Prepare the means and standard deviations of T and arg(t) from the same quantities of t
"""
function prep_T_argt_statistics(t_real_means, t_real_stds, t_imag_means, t_imag_stds)
    T_means     = @. t_real_means^2 + t_imag_means^2
    T_stds      = @. sqrt(4*t_real_means^2*t_real_stds^2 + 4*t_imag_means^2*t_imag_stds^2)
    phase_means = @. atan(t_imag_means, t_real_means)
    phase_stds  = @. sqrt(  (-t_imag_means/(t_real_means^2 + t_imag_means^2))^2*t_real_stds^2
                          + ( t_real_means/(t_real_means^2 + t_imag_means^2))^2*t_imag_stds^2 )
    return T_means, T_stds, phase_means, phase_stds
end


"""
Prepare the means and standard deviations of the transmission as calculated over
a set of imperfect atomic arrays
"""
function prep_compareImperfectArray_transmission_vs_X(Δ_index, T_meanss, T_stdss, phase_meanss, phase_stdss, T_indepDecayss, phase_indepDecayss, refrIndex_realss, refrIndex_imagss, refrIndex_real_indepDecayss, refrIndex_imag_indepDecayss)
    T_means                    = [T_means[Δ_index]                    for T_means                    in T_meanss]
    T_stds                     = [T_stds[Δ_index]                     for T_stds                     in T_stdss]
    phase_means                = [phase_means[Δ_index]                for phase_means                in phase_meanss]
    phase_stds                 = [phase_stds[Δ_index]                 for phase_stds                 in phase_stdss]
    T_indepDecays              = [T_indepDecays[Δ_index]              for T_indepDecays              in T_indepDecayss]
    phase_indepDecays          = [phase_indepDecays[Δ_index]          for phase_indepDecays          in phase_indepDecayss]
    refrIndex_reals            = [refrIndex_reals[Δ_index]            for refrIndex_reals            in refrIndex_realss]
    refrIndex_imags            = [refrIndex_imags[Δ_index]            for refrIndex_imags            in refrIndex_imagss]
    refrIndex_real_indepDecays = [refrIndex_real_indepDecays[Δ_index] for refrIndex_real_indepDecays in refrIndex_real_indepDecayss]
    refrIndex_imag_indepDecays = [refrIndex_imag_indepDecays[Δ_index] for refrIndex_imag_indepDecays in refrIndex_imag_indepDecayss]
    return T_means, T_stds, phase_means, phase_stds, T_indepDecays, phase_indepDecays, refrIndex_reals, refrIndex_imags, refrIndex_real_indepDecays, refrIndex_imag_indepDecays
end


"""
Prepare the loss, absolute value of weights, and absolute value of resonances
"""
function prep_loss_weights_resonances(t, r, weights, resonances)
    loss = 1 .- abs2.(t) .- abs2.(r)
    weights_abs = abs.(weights)
    resonances_abs = broadcast(x -> abs.(x), resonances)
    return loss, weights_abs, resonances_abs
end


"""
Prepare the loss, absolute value of resonances, maximum absolute value of resonances, 
and excitation-sector populations
"""
function prep_loss_resonances_pops(t, r, resonances, eigenModesMatrix, noPhonons)
    loss = 1 .- abs2.(t) .- abs2.(r)
    resonances_abs = broadcast(x -> abs.(x), resonances)
    resonances_abs_max = maximum.(resonances_abs)
    exci_pops = [x[1] for x in statePopulations.(eachcol(eigenModesMatrix), Ref(noPhonons))]
    return loss, resonances_abs, resonances_abs_max, exci_pops
end


"""
Prepare title for the transmission plot
"""
function prep_transmission_title(SP)
    title_components = [
        "Δvari: " * SP.ΔvariDescription,
        "Array: " * SP.arrayDescription,
        "tildeG flags: $(join(SP.tildeG_flags, ", "))",
        "να = $(ro.(SP.να)), ηα = $(ro.(SP.ηα)),",
        "Dipole moments: " * SP.dDescription,
        "No phonons: " * string(SP.noPhonons)
    ]
    if SP.include3rdLevel
        append!(title_components, [
          "Control drive, detuning, Rabi freq.: $(SP.cDriveDescription), $(SP.Δc), $(SP.Ωc)",
        ])
        if SP.cDriveDescription == "plW"
            append!(title_components, [
                "Control drive wavenumber: $(join(ro.(SP.cDriveArgs.kc), ","))"
                ])
        end
    end
    return join(title_components, "\n")
end



"""
Prepare title for the transmission plot with imperfect arrays
"""
function prep_imperfectArray_transmission_title(SP)
    title_components = [
        "Δvari: " * SP.ΔvariDescription,
        "Array: " * SP.arrayDescription,
        "tildeG flags: $(join(SP.tildeG_flags, ", "))",
        "n_inst: $(SP.n_inst)",
        "να = $(ro.(SP.να)), ηα = $(ro.(SP.ηα)),",
        "Dipole moments: " * SP.dDescription,
        "No phonons: " * string(SP.noPhonons)
    ]
    if SP.include3rdLevel
        append!(title_components, [
          "Control drive, detuning, Rabi freq.: $(SP.cDriveDescription), $(SP.Δc), $(SP.Ωc)",
        ])
        if SP.cDriveDescription == "plW"
            append!(title_components, [
                "Control drive wavenumber: $(join(ro.(SP.cDriveArgs.kc), ","))"
                ])
        end
    end
    return join(title_components, "\n")
end


"""
Prepare title for the state plot
"""
function prep_state_title(SP, Δ)
    title_components = [
        "Δvari: " * SP.ΔvariDescription,
        "Array: " * SP.arrayDescription,
        "tildeG flags: $(join(SP.tildeG_flags, ", "))",
        "να = $(ro.(SP.να)), ηα = $(ro.(SP.ηα)),",
        "Dipole moments: " * SP.dDescription,
        "No phonons: " * string(SP.noPhonons),
        "Δ = $(ro(Δ))"
    ]
    if SP.include3rdLevel
        append!(title_components, [
          "Control drive, detuning, Rabi freq.: $(SP.cDriveDescription), $(SP.Δc), $(SP.Ωc)",
        ])
        if SP.cDriveDescription == "plW"
            append!(title_components, [
                "Control drive wavenumber: $(join(ro.(SP.cDriveArgs.kc), ","))"
                ])
        end
    end
    return join(title_components, "\n")
end


"""
Prepare title for the coupling matrix eigenmodes plot
"""
function prep_GnmEigenModes_title(SP)
    title_components = [
        "Δvari: " * SP.ΔvariDescription,
        "Array: " * SP.arrayDescription,
        "tildeG flags: $(join(SP.tildeG_flags, ", "))",
        "να = $(ro.(SP.να)), ηα = $(ro.(SP.ηα)),",
        "Dipole moments: " * SP.dDescription,
        "No phonons: " * string(SP.noPhonons)
    ]
    if SP.include3rdLevel
        append!(title_components, [
          "Control drive, detuning, Rabi freq.: $(SP.cDriveDescription), $(SP.Δc), $(SP.Ωc)",
        ])
        if SP.cDriveDescription == "plW"
            append!(title_components, [
                "Control drive wavenumber: $(join(ro.(SP.cDriveArgs.kc), ","))"
                ])
        end
    end
    return join(title_components, "\n")
end


"""
Prepare title for the coupling matrix eigenenergies plot
"""
function prep_GnmEigenEnergies_title(SP)
    title_components = [
        "Δvari: " * SP.ΔvariDescription,
        "Array: " * SP.arrayDescription,
        "tildeG flags: $(join(SP.tildeG_flags, ", "))",
        "να = $(ro.(SP.να)), ηα = $(ro.(SP.ηα)),",
        "Dipole moments: " * SP.dDescription,
        "No phonons: " * string(SP.noPhonons)
    ]
    if SP.include3rdLevel
        append!(title_components, [
          "Control drive, detuning, Rabi freq.: $(SP.cDriveDescription), $(SP.Δc), $(SP.Ωc)",
        ])
        if SP.cDriveDescription == "plW"
            append!(title_components, [
                "Control drive wavenumber: $(join(ro.(SP.cDriveArgs.kc), ","))"
                ])
        end
    end
    return join(title_components, "\n")
end


"""
Prepare title for the transmission with coupling matrix eigenenergies etc. plot
"""
function prep_transmissionWithGnmEigenEnergies_title(SP)
    title_components = [
        "Δvari: " * SP.ΔvariDescription,
        "Array: " * SP.arrayDescription,
        "tildeG flags: $(join(SP.tildeG_flags, ", "))",
        "να = $(ro.(SP.να)), ηα = $(ro.(SP.ηα)),",
        "Dipole moments: " * SP.dDescription,
        "No phonons: " * string(SP.noPhonons)
    ]
    if SP.include3rdLevel
        append!(title_components, [
          "Control drive, detuning, Rabi freq.: $(SP.cDriveDescription), $(SP.Δc), $(SP.Ωc)",
        ])
        if SP.cDriveDescription == "plW"
            append!(title_components, [
                "Control drive wavenumber: $(join(ro.(SP.cDriveArgs.kc), ","))"
                ])
        end
    end
    return join(title_components, "\n")
end





