
include("preamble.jl")


# ================================================
#   Main functions
# ================================================
function define_ω_ρf_n_ranges()
    λ0  = 852       #nm, guided mode wavelength, transition frequency of cs133
    ω0  = 2π/λ0     #nm^-1, guided mode angular frequency
    ρf0 = 200       #nm, fiber radius
    n0  = 1.45      #unitless, index of refraction
    
    ω_specs  = (0.1*ω0, 2*ω0, 100) 
    ρf_specs = (ρf0, ρf0, 1)
    n_specs  = (n0, n0, 1)
    return ω_ρf_n_ranges(ω_specs, ρf_specs, n_specs)
end 


function define_SP_BerlinCS()
    # Fiber specs from "Magic-wavelength nanofiber-based two-color dipole trap with sub-λ/2 spacing"
    λ0  = 852       #nm, guided mode wavelength, transition frequency of cs133
    ω0  = 2π/λ0     #nm^-1, guided mode angular frequency
    γ0  = 2π*5.22e3 #kHz, free decay rate of cs133
    ρf0 = 200       #nm, fiber radius
    n0  = 1.45      #unitless, index of refraction
    
    # Atomic array specs
    ρa0 = 550   #nm, atomic array radial coordinate
    a0  = 300   #nm, atomic array lattice constant
    
    # Trap specs
    ν0_radial    = 2π*109 #kHz, radial atomic trap angular frequency
    ν0_axial     = 2π*139 #kHz, axial atomic trap angular frequency
    ν0_azimuthal = 2π*62  #kHz, azimuthal atomic trap angular frequency (estimate from graph: 18 kHz, but usually half of the others)
    να0 = [ν0_radial, ν0_azimuthal, ν0_axial] #trap frequencies in a Cartesian basis (x, y, z) which matches with (radial, azimuthal, axial) if the position is taken to be on the x-axis
    
    # Recoil energy
    νR0 = 2π*2.0663 #kHz, recoil energy of cesium atoms (as an angular frequency)
    
    # Lamb-Dicke parameters in a Cartesian basis (x, y, z)
    ηα0 = @. sqrt(νR0/να0)
    
    # Unitless versions
    ρf0_ul = ρf0/λ0 #unitless version of ρf0
    ρa0_ul = ρa0/λ0 #unitless version of ρa0
    a0_ul  = a0/λ0  #unitless version of a0
    να0_ul = να0/γ0 #unitless version of να0
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-0.5, 0.5, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "Gaussian"
    Δvari_args = -3, 5*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = ηα0 .* [0.1, 0.2, 0.1]
    # ηα = ηα0 * 0.4
    # ηα = [0.01, 0.01, 0.01]
    ηα = [0., 0., 0.]
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set array specs and generate array, as well as description for postfix
    N = 50
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 1.0
    pos_unc = any(ηα .!= 0) ? 0.0 : 0.0
    pos_unc = any(ηα .!= 0) ? 0.0 : ηα0/ωa * 0.4
    n_inst  = any(ηα .!= 0) ?   1 : 10
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Prepare initial state for time evolution, as well as description for postfix
    N_eff = Int(floor(N*ff))
    initialState = groundstate(N_eff, all(ηα .== 0))
    initialStateDescription = "gs"
    
    # Atomic dipole moment
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul)
    d = "chiral"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    # incField_wlf = [(1, 1, 1), (1, -1, 1)]
    incField_wlf = []
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively, usually (true, false))
    approx_Grm_trans = (true, false)
    
    # Whether to save individual results (Im_Grm_trans, steady states, time evolutions)
    save_individual_res = n_inst == 1 && ff == 1 && pos_unc == 0
    # save_individual_res = pos_unc == 0
    
    # Ranges of z and x values to define r_field for calculating the radiated E-field
    arrayL = (N - 1)*a0_ul
    # z_range = range(-0.5*arrayL, 1.5*arrayL, 60)
    # x_range = range(ρa0_ul - 0.3*arrayL, ρa0_ul + 0.3*arrayL, 60)
    z_range = range(-10, arrayL + 10, 60)
    x_range = range(-ρf0_ul - 10, ρf0_ul + ρa0_ul + 10, 60)
    y_fix   = ρa0_ul
    
    
    if n_inst != 1
        return [SysPar(ρf0_ul, n0, ωa,
                       Δ_specs,
                       ΔvariDependence, Δvari_args, ΔvariDescription,
                       tspan, dtmax, initialState, initialStateDescription,
                       arrayType, N, ρa0_ul, a0_ul, ff, pos_unc,
                       να0_ul, ηα,
                       d, incField_wlf, save_individual_res, approx_Grm_trans,
                       z_range, x_range, y_fix) for _ in 1:n_inst]
    else
        return SysPar(ρf0_ul, n0, ωa,
                      Δ_specs,
                      ΔvariDependence, Δvari_args, ΔvariDescription,
                      tspan, dtmax, initialState, initialStateDescription,
                      arrayType, N, ρa0_ul, a0_ul, ff, pos_unc,
                      να0_ul, ηα,
                      d, incField_wlf, save_individual_res, approx_Grm_trans,
                      z_range, x_range, y_fix)
    end
end


function define_SP_Olmos()
    # Fiber specs from "Modified dipole-dipole interactions in the presence of a nanophotonic waveguide"
    λ0 = 852  #nm, guided mode wavelength, transition frequency of cs133
    n  = 1.45 #unitless, index of refraction
    # n = 1.452467
    ρf = 250  #nm, Fiber radius
    
    # Unitless fiber radius
    ρf_ul = ρf/λ0 
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-10, 10, 1000)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = nothing
    ΔvariDescription = ΔvariDescription(ΔvariDependence, Δvari_args)
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
    arrayType = "1Dchain"
    N  = 10
    ρa = ρf_ul + 50/λ0
    a  = 0.1
    
    # Phonon bare energies, i.e. trap frequencies
    να = [0, 0, 0]
    
    # Lamb-Dicke parameters
    ηα = [0, 0, 0]
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, all(ηα .== 0))
    initialStateDescription = "gs"
     
    # Atomic dipole moment
    d = [1, 0, 0]
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively)
    approx_Grm_trans = (true, false)
    
    return SysPar(ρf_ul, n, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N, ρa, a,
                  να, ηα,
                  d, incField_wlf, approx_Grm_trans)
end


function define_SP_Rauschenbeutel()
    # Fiber specs from "Modified dipole-dipole interactions in the presence of a nanophotonic waveguide"
    λ0 = 852  #nm, guided mode wavelength, transition frequency of cs133
    n  = 1.45 #unitless, index of refraction
    ρf = 250  #nm, Fiber radius
    
    # Unitless fiber radius
    ρf_ul = ρf/λ0
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-30, 30, 100)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = nothing
    ΔvariDescription = ΔvariDescription(ΔvariDependence, Δvari_args)
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 5)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
    arrayType = "1Dchain"
    N  = 2
    ρa = ρf_ul + 100/λ0
    a  = 0.1
    
    # Phonon bare energies, i.e. trap frequencies
    να = [0, 0, 0]
    
    # Lamb-Dicke parameters
    ηα = [0, 0, 0]
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, all(ηα .== 0))
    initialStateDescription = "gs"
     
    # Atomic dipole moment
    d = conj([1im, 0, -1]/sqrt(2))
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively)
    approx_Grm_trans = (true, false)
    
    return SysPar(ρf_ul, n, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N, ρa, a,
                  να, ηα,
                  d, incField_wlf, approx_Grm_trans)
end


function define_SP_Chang()
    # Fiber specs from "Modified dipole-dipole interactions in the presence of a nanophotonic waveguide"
    n  = 2 #unitless, index of refraction
    ρf = 1.2/(2π)  #unitless, fiber radius
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-10, 10, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = nothing
    ΔvariDescription = ΔvariDescription(ΔvariDependence, Δvari_args)
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 5)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
    arrayType = "1Dchain"
    N  = 20
    ρa = 1.5*ρf
    a  = 0.25
    
    # Phonon bare energies, i.e. trap frequencies
    να = [0, 0, 0]
    
    # Lamb-Dicke parameters
    ηα = [0, 0, 0]
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, all(ηα .== 0))
    initialStateDescription = "gs"
     
    # Atomic dipole moment
    d = [1, 0, 0]
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively)
    approx_Grm_trans = (true, false)
    
    return SysPar(ρf, n, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N, ρa, a,
                  να, ηα,
                  d, incField_wlf, approx_Grm_trans)
end


function main()
    # Define system parameters
    # ωρfn_ranges = define_ω_ρf_n_ranges()
    SP = define_SP_BerlinCS()
    # SP = define_SP_Olmos()
    # SP = define_SP_Rauschenbeutel()
    # SP = define_SP_Chang()
    # show(SP)
     
    
    
    # plot_propConst_inOutMom(ωρfn_ranges)
    # plot_coupling_strengths(SP)
    # plot_σBαTrajectories_σBαSS(SP)
    # plot_transmission_vs_Δ(SP)
    plot_classDisorder_transmission_vs_Δ(SP)
    # plot_steadyState_radiation_Efield(SP)
    # plot_radiation_Efield(SP)
    # plot_GnmEigenModes(SP)
    # plot_emissionPatternOfGnmeigenModes(SP)
    # plot_GnmEigenEnergies(SP)
    # plot_lossWithGnmEigenEnergies(SP)
    
    
    return nothing
end


# ================================================
#   Generate figures
# ================================================
function plot_propConst_inOutMom(ω_ρf_n_ranges)
    κ = scan_propConst(ω_ρf_n_ranges)
    
    # Plot propConst with a specific choice of ρf and n from the ranges
    ρf_ind = 1
    n_ind  = 1
    fig_propConst_vs_ω(ω_ρf_n_ranges.ω_range, κ[:, ρf_ind, n_ind], ω_ρf_n_ranges.ρf_range[ρf_ind], ω_ρf_n_ranges.n_range[n_ind])
    
    # Plot the momenta inside and outside the fiber
    h = in_momentum.(κ[:, ρf_ind, n_ind], ω_ρf_n_ranges.ω_range, ω_ρf_n_ranges.n_range[n_ind])
    q = out_momentum.(κ[:, ρf_ind, n_ind], ω_ρf_n_ranges.ω_range)
    fig_inout_momenta_vs_ω(ω_ρf_n_ranges.ω_range, h, q, ω_ρf_n_ranges.n_range[n_ind])
end


function plot_coupling_strengths(SP)
    
    # PRESENTLY HARDCODED for define_SP_Rauschenbeutel
    
    # Guided mode local decay as a function of distance to fiber
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γgm = get_Γgm.(Ref(SP.fiber), Ref(SP.d), r_range, r_range)
    x_label = L"$ \rho - \rho_f $"
    y_label = L"$ \Gamma_{gm, nn} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γgm, x_label, y_label)
    
    # Guided mode dissipative interaction as a function of radial distance
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γgm = get_Γgm.(Ref(SP.fiber), Ref(SP.d), Ref(r_range[1]), r_range)
    x_label = L"$ \rho_2 - \rho_f $"
    y_label = L"$ \Gamma_{gm, 12} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γgm, x_label, y_label)
    
    # Radiation mode local decay as a function of distance to fiber
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γrm = get_Γrm.(Ref(SP.fiber), Ref(SP.d), r_range, r_range, SP.approx_Grm_trans)
    x_label = L"$ \rho - \rho_f $"
    y_label = L"$ \Gamma_{rm, nn} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γrm, x_label, y_label)
    
    # Radation mode dissipative interaction as a function of radial distance
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γrm = get_Γrm.(Ref(SP.fiber), Ref(SP.d), Ref(r_range[1]), r_range, SP.approx_Grm_trans)
    x_label = L"$ \rho_2 - \rho_f $"
    y_label = L"$ \Gamma_{rm, 12} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γrm, x_label, y_label)
    
end


function plot_σBαTrajectories_σBαSS(SP)
    Δ = 0.0
    σBα_SS = calc_σBα_steadyState(SP, Δ)
    xTrajectories = timeEvolution(SP, Δ)
    # xTrajectories = timeEvolution_eigenmodes(SP, Δ)
    
    if all(SP.ηα .== 0)
        times, σTrajectories = prep_times_σTrajectories(xTrajectories, SP.N)
        fig_σTrajectories_σSS(times, σTrajectories, σBα_SS)
    else
        times, σTrajectories, BαTrajectories = prep_times_σBαTrajectories(xTrajectories, SP.N)
        fig_σBαTrajectories_σBαSS(times, σTrajectories, BαTrajectories, σBα_SS...)
    end
end


function plot_transmission_vs_Δ(SP)
    σBα_scan = scan_σBα_steadyState(SP)
    t = calc_transmission.(Ref(SP), σBα_scan)
    # t = scan_transmission_eigenmodes(SP)
    
    T, phase = prep_transmission(t)
    titl = prep_transmission_title(SP)
    fig_transmission_vs_Δ(SP.Δ_range, T, phase, titl)
end


function plot_classDisorder_transmission_vs_Δ(SPs)
    if typeof(SPs) == SysPar throw(ArgumentError("plot_classDisorder_transmission_vs_Δ requires n_inst > 1")) end
    
    n_inst = length(SPs)
    postfix = get_postfix(SPs[1].Δ_specs, SPs[1].ΔvariDescription, SPs[1].d, SPs[1].να, SPs[1].ηα, SPs[1].incField_wlf, n_inst, SPs[1].arrayDescription, SPs[1].fiber.postfix)
    filename = "T_phase" * postfix
    folder = "classDisorder_T_phase/"
    
    if isfile(saveDir * folder * filename * ".txt") 
        T_means, T_stds, phase_means, phase_stds = eachrow(load_as_txt(saveDir * folder, filename))
    else
        ts = []
        for SP in SPs
            σBα_scan = scan_σBα_steadyState(SP)
            push!(ts, calc_transmission.(Ref(SP), σBα_scan))
        end
            
        # Prepare means and standard deviations of (squared) magnitudes and phases
        T_means, T_stds, phase_means, phase_stds = prep_classDisorder_transmission(ts)
        formattedResult = vectorOfRows2Matrix([T_means, T_stds, phase_means, phase_stds])
        save_as_txt(formattedResult, saveDir * folder, filename)
    end
    
    fig_classDisorder_transmission_vs_Δ(SPs[1].Δ_range, T_means, T_stds, phase_means, phase_stds)
end


function plot_steadyState_radiation_Efield(SP)
    Δ = -0.25
    rs = [site[3] for site in SP.array]
    if all(SP.ηα .== 0) σ_SS = calc_σBα_steadyState(SP, Δ)
    else                σ_SS, Bα_SS = calc_σBα_steadyState(SP, Δ)
    end
    ks, σ_SS_FT = discFourierTransform(σ_SS, SP.a, true, 1000)
    
    E = calc_radiation_Efield.(Ref(SP), Ref(σ_SS), SP.r_field)
    intensity = norm.(E).^2
    # intensity = zeros(size(SP.r_field))
    
    titl = prep_state_title(SP, Δ)
    fig_state(rs, σ_SS, ks, σ_SS_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, SP.fiber.propagation_constant, titl)
end


function plot_radiation_Efield(SP)
    if any(SP.ηα .!= 0) throw(ArgumentError("plot_radiation_Efield is not implemented for the case of including phonons")) end
    
    Δ = 0.0
    σ_SS = calc_σBα_steadyState(SP, Δ)
    E = calc_radiation_Efield.(Ref(SP), Ref(σ_SS), SP.r_field)
    intensity = norm.(E).^2
    
    fig_radiation_Efield(SP.z_range, SP.x_range, intensity, SP.ρf, SP.array)
end


function plot_GnmEigenModes(SP)
    if any(SP.ηα .!= 0) throw(ArgumentError("plot_GnmEigenModes is not implemented for the case of including phonons")) end
    
    Gnm = get_tildeGs(SP.fiber, SP.d, SP.array, SP.save_individual_res, SP.approx_Grm_trans)
    
    eigenEnergies, eigenModes, dominant_ks = spectrum(Gnm, SP.a)
    eigenModes_FT = discFourierTransform.(eigenModes, SP.a, true, 1000)
    
    rs = [site[3] for site in SP.array]
    iter_list = collect(zip(eigenModes, eigenModes_FT, eigenEnergies, dominant_ks))[sortperm(dominant_ks)]
    for (mode, (ks, mode_FT), eigenEnergy, dom_k) in iter_list
        # if -7 < dom_k < -5
            E = calc_radiation_Efield.(Ref(SP), Ref(mode), SP.r_field)
            intensity = norm.(E).^2
            titl = prep_GnmEigenModes_title(SP)
            fig_GnmEigenModes(rs, mode, ks, mode_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, eigenEnergy, SP.fiber.propagation_constant, titl)
        # end
    end
end


function plot_emissionPatternOfGnmeigenModes(SP)
    if any(SP.ηα .!= 0) throw(ArgumentError("plot_emissionPatternOfGnmeigenModes is not implemented for the case of including phonons")) end
    
    Gnm = get_tildeGs(SP.fiber, SP.d, SP.array, SP.save_individual_res, (true, false))
    eigenEnergies, eigenModes = eigbasis(Gnm)
    
    for evec in eigenModes
        E = calc_radiation_Efield.(Ref(SP), Ref(evec), SP.r_field)
        intensity = norm.(E).^2
        fig_radiation_Efield(SP.z_range, SP.x_range, intensity, SP.ρf, SP.array)
    end
end


function plot_GnmEigenEnergies(SP)
    if any(SP.ηα .!= 0) throw(ArgumentError("plot_GnmEigenEnergies is not implemented for the case of including phonons")) end
    
    Δvari, tildeΩ, tildeG = get_parameterMatrices(SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.save_individual_res, SP.approx_Grm_trans)
    # tildeG = get_tildeG0(SP.fiber, SP.d, SP.array)
    
    eigenEnergies, eigenModes, dominant_ks = spectrum(Δvari + tildeG, SP.a)
    collΔ, collΓ = collEnergies_from_eigenEnergies(eigenEnergies)
    weights, resonances = transmission_eigenmodes_weights_resonances(SP.Δ_range, tildeΩ, eigenEnergies, eigenModes, SP.fiber.propagation_constant_derivative)
    weights_abs = abs.(weights)
    
    titl = prep_GnmEigenEnergies_title(SP)
    fig_eigenEnergies_vs_k(dominant_ks, collΔ, collΓ, weights_abs, SP.fiber.propagation_constant, titl) 
end


function plot_lossWithGnmEigenEnergies(SP)
    if any(SP.ηα .!= 0) throw(ArgumentError("plot_lossWithGnmEigenEnergies is not implemented for the case of including phonons")) end
    
    σBα_scan = scan_σBα_steadyState(SP)
    t = calc_transmission.(Ref(SP), σBα_scan)
    # t = scan_transmission_eigenmodes(SP)
    
    Δvari, tildeΩ, tildeG = get_parameterMatrices(SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.save_individual_res, SP.approx_Grm_trans)
    eigenEnergies, eigenModes = eigbasis(Δvari + tildeG)
    collΔ, collΓ = collEnergies_from_eigenEnergies(eigenEnergies)
    weights, resonances = transmission_eigenmodes_weights_resonances(SP.Δ_range, tildeΩ, eigenEnergies, eigenModes, SP.fiber.propagation_constant_derivative)
    
    loss, weights_abs, resonances_abs = prep_loss_weights_resonances(t, weights, resonances)
    titl = prep_transmissionWithGnmEigenEnergies_title(SP)
    fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, collΓ, weights_abs, titl)
end





println("\n -- Running main() -- \n")
@time main()




# TODO list:

# Make chiral dipole moment work for atoms on the opposite side of the array

# Maybe calculate tildeG, etc., in SP? They are not expected to change anyway
    # Only calculate them if needed..?

# Implement titles for figures that show parameters etc.

# Implement that Im_Grm_trans_ simply returns zero for a certain size of relative_z/fiber_radius or so? (Will go to zero for large inter-atom distance)

# Implement Gnm functions for the case of including phonons

# Implement radiation_Efield for the case of including phonons

# Calculate Poynting vector and plot instead of/together with emission pattern?

# Get it to work on the cluster
    # Use MPI?
    # clean up before moving to the cluster: comments, naming, minor issues below, checks of derivatives

# Systematically scan the effect of η and ff
    # Pick a value of Δ with T ~ 1 and phase ~ pi, and plot these things as a function of η and ff?
    
# Calculate reflection and loss

# Implement n_inst in a better way? That is, dont have a list of SPs, but include the many instantiations in the same SP?

# Implement saving and loading of the parameter matrices?

# Clean up and split up files? physics.jl and utility.jl in particular

# Implement non-lazy version of get_tildeGs(fiber, d::String... ? Presumably significantly faster when exploiting knowledge of which components etc. are actually needed, but also very messy...

# Implement the correct real part of the radiation GF
    # But maybe not relevant since atoms are far-ish from the fiber, such that using vacuum-approximation should be good
    
# When calculating Im_Grm_trans, only calculate needed components? And exploit that some entries appear to be zero depending on derivOrder and α,
    # presumably because of the simple array structure

# Reconsider factoring of code, particularly whether calculation of parameters etc. 
    # should take place outside of scanning loops
    # One step deeper, the calculation of mode components should perhaps also be done 
    # outside of loops over atomic array
# Perhaps make a struct that contains all information needed to make a scan of some type?
    # Prepare for a scan by optimally preparing all parameters
    # then perform scan with the prepared quantities as input
    # Each type of scan can have its own struct-thing?
    # Functions that calculate different things can either take SP-like parameters
    # and calculate parameters on its own, or it takes the already prepared parameters?
    # or some wrapper which can direct to do either one thing or the other. Such that
    # different calculation functions can still be called without making a whole scan...
# When calculating GFs (and maybe other quantities), go through which calculations actually depend on the input that is varied (positions, frequencies, ...) 
    # and refactor to optimize this (no need to repeat the same position-independent calculation for each point in the array, etc.)

# Use only zeroth order tildeFα in σBα_steadyState? Reduce to second order in eta in other ways (i.e. calculate for eta=0 and eta≠0 to extract exact dependencies)? Time evolution anyways also finds results that are effectively higher order
    
# Update interface comments and make a description of notation/conventions somewhere

# Make small ω version of fiber_equation work? Implement HomotopyContinuation solution?

# Consider making structs for x-vector and σBα, to make their structure easier to get an overview of

# Split utility.jl according to theme or which file has need of the functions
    # simply code utility/no physics concatenated
    # physics functions according to which file has need of them
    
# Consider never using JLD2
    # Invent some packing/unpacking scheme to put any kind of data into real matrices

# Consider using StaticArrays in some places?

# Consider writing tests...
    # Test/check validity of SP (d is normalized, some parameters are non-negative), which could also serve as a reminder of the implicit assumptions regarding parameters
    # Test that transmission is only a phase for 1 atom and no coupling to radiation spectrum
    # Test that long-time time evolution matches with steady state
    # Test calculations in a limit where it can be done analytically? (One atom...)