

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
Prepare the squared magnitude and phase of the transmission
"""
function prep_squaredNorm_phase(t)
    return abs2.(t), angle.(t)
end


"""
Prepare the means and standard deviations of the transmission as calculated over
a set of imperfect atomic arrays
"""
function prep_imperfectArray_transmission(ts)
    Ts     = [abs2.(t) for t in ts]
    phases = [vcat([angle(t[1])], angle(t[1]) .+ cumsum(angle.(t[2:end]./t[1:end-1]))) for t in ts]
    
    T_mat      = vectorOfRows2Matrix(Ts)
    phases_mat = vectorOfRows2Matrix(phases)
    return           squeeze(mean(T_mat, dims=1)), 
                     squeeze( std(T_mat, dims=1)), 
           wrapPhase(squeeze(mean(phases_mat, dims=1))), 
                     squeeze( std(phases_mat, dims=1))
end


"""
Prepare the means and standard deviations of the transmission as calculated over
a set of imperfect atomic arrays
"""
function prep_compareImperfectArray_transmission_vs_ffOrηα(Δ_index, T_meanss, T_stdss, phase_meanss, phase_stdss)
    T_means     = [T_means[Δ_index]     for T_means     in T_meanss]
    T_stds      = [T_stds[Δ_index]      for T_stds      in T_stdss]
    phase_means = [phase_means[Δ_index] for phase_means in phase_meanss]
    phase_stds  = [phase_stds[Δ_index]  for phase_stds  in phase_stdss]
    return T_means, T_stds, phase_means, phase_stds
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
Prepare title for the transmission plot
"""
function prep_transmission_title(SP)
    title_components = [
        "Δvari: " * SP.ΔvariDescription,
        "Array: " * SP.arrayDescription,
        "tildeG flags: $(join(SP.tildeG_flags, ", "))",
        "να = $(ro.(SP.να)), ηα = $(ro.(SP.ηα)),",
        "Dipole moments: " * SP.dDescription
    ]
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
        "Dipole moments: " * SP.dDescription
    ]
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
        "Δ = $(ro(Δ))"
    ]
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
        "Dipole moments: " * SP.dDescription
    ]
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
        "Dipole moments: " * SP.dDescription
    ]
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
        "Dipole moments: " * SP.dDescription
    ]
    return join(title_components, "\n")
end





