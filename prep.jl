

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
function prep_transmission(t)
    return abs2.(t), angle.(t)
end


"""
Prepare the means and standard deviations of the transmission as calculated over
a set of imperfect atomic arrays
"""
function prep_imperfectArray_transmission(ts)
    ts_mat = vectorOfRows2Matrix(ts)
    
    T_mat     = abs2.(ts_mat)
    phase_mat = unwrapPhase.(angle.(ts_mat))
    
    return            squeeze(mean(T_mat, dims=1)), 
                      squeeze( std(T_mat, dims=1)), 
           wrapPhase.(squeeze(mean(phase_mat, dims=1))), 
                      squeeze( std(phase_mat, dims=1))
end


"""
Prepare the loss, absolute value of weights, and absolute value of resonances
"""
function prep_loss_weights_resonances(t, weights, resonances)
    loss = 1 .- abs2.(t)
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
        "να = $(ro.(SP.να)), ηα = $(ro.(SP.ηα)),",
        "Dipole moments: " * SP.dDescription
    ]
    return join(title_components, "\n")
end





