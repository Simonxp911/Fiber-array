

include("preamble.jl")
const saveDir = "C:/Users/Simon/Forskning/Data/fiber_array_data/"

using GLMakie #plotting, specialized interactive plots


# ================================================
#   Main functions
# ================================================
function define_SP_BerlinCs()
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
    ηα0 = @. sqrt(νR0/να0) #[0.1377, 0.1826, 0.1219]
    
    # Unitless versions
    ρf0_ul = ρf0/λ0 #unitless version of ρf0, 0.2347
    ρa0_ul = ρa0/λ0 #unitless version of ρa0, 0.6455
    a0_ul  = a0/λ0  #unitless version of a0 , 0.3521
    να0_ul = να0/γ0 #unitless version of να0, [0.0209, 0.0119, 0.0266]
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-1, 0.0, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "Gaussian"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = ηα0 * 0.4
    ηα = [0., 0., 0.]
    
    # Whether phonons are excluded or not from the calculations (a finite ηα but noPhonons = true will result in including ground state motion into tildeG)
    noPhonons = all(ηα .== 0)
    # noPhonons = true
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set number of atomic sites 
    N_sites = 750
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 1.0
    pos_unc = 0.0
    # pos_unc = ηα0/ωa
    n_inst = 1
    
    # Generate the array, its description, and the number of atoms
    array, arrayDescription, N = get_array(arrayType, N_sites, ρa0_ul, a0_ul, ff, pos_unc, n_inst)
    
    # Time span and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, noPhonons)
    initialStateDescription = "gs"
    
    # Atomic dipole moment
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul, array)
    d = "chiral"
    dDescription = "chiral"
    # d = rightCircularDipoleMoment(array)
    # dDescription = "rgtCrc"
    # d = [[1, 0, 0] for site in array]
    # dDescription = "xPol"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    if typeof(d) == String incField_wlf = [] end
    
    # Whether to include the guided contribution, the radiated contribution, and the radiated interactions when calculating tildeG
    tildeG_flags = (true, true, true)
    
    # Absolute tolerance in the calculations of Im_Grm_trans
    abstol_Im_Grm_trans = 1e-5
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively, usually (true, false))
    approx_Grm_trans = (true, false)
    
    # Whether to interpolate Im_Grm_trans
    interpolate_Im_Grm_trans = arrayType == "randomZ"
    
    # Whether to save individual results (Im_Grm_trans, steady states, time evolutions)
    save_Im_Grm_trans = pos_unc == 0 && arrayType != "randomZ" && !interpolate_Im_Grm_trans
    save_steadyState  = n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    save_timeEvol     = n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a0_ul
    margin = 0.1*arrayL
    z_range = range(-margin, arrayL + margin, 30)
    x_range = range(-ρf0_ul - margin, ρf0_ul + ρa0_ul + margin, 30)
    y_fix   = ρa0_ul
    
    # Get the interpolation function for the imaginary, transverse part of the radiation Green's function, if needed
    if interpolate_Im_Grm_trans interpolation_Im_Grm_trans = interpolation1D_Im_Grm_trans(Fiber(ρf0_ul, n0, ωa), Int(ceil(arrayL/0.1)) + 1, ρa0_ul, 0.1, noPhonons) else interpolation_Im_Grm_trans = nothing end
    
    
    return SysPar(ρf0_ul, n0, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N_sites, ρa0_ul, a0_ul, ff, pos_unc, n_inst, array, arrayDescription, N,
                  να0_ul, ηα, noPhonons,
                  d, dDescription, incField_wlf, tildeG_flags, 
                  interpolate_Im_Grm_trans, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                  save_steadyState, save_timeEvol,
                  interpolation_Im_Grm_trans,
                  z_range, x_range, y_fix) 
end


function define_SP_BerlinSr()
    # Fiber specs 
    λ0  = 689       #nm, guided mode wavelength, transition frequency of Sr
    ω0  = 2π/λ0     #nm^-1, guided mode angular frequency
    γ0  = 2π*7.4    #kHz, free decay rate of Sr88
    ρf0 = 115       #nm, fiber radius
    n0  = 1.45      #unitless, index of refraction
    
    # Atomic array specs
    ρa0 = 275   #nm, atomic array radial coordinate
    a0  = 200   #nm, atomic array lattice constant
    
    # Trap specs (might be lower than the Cs experiment)
    ν0_radial    = 2π*109 #kHz, radial atomic trap angular frequency
    ν0_axial     = 2π*139 #kHz, axial atomic trap angular frequency
    ν0_azimuthal = 2π*62  #kHz, azimuthal atomic trap angular frequency (estimate from graph: 18 kHz, but usually half of the others)
    να0 = [ν0_radial, ν0_azimuthal, ν0_axial] #trap frequencies in a Cartesian basis (x, y, z) which matches with (radial, azimuthal, axial) if the position is taken to be on the x-axis
    
    # Recoil energy
    νR0 = 2π*4.775 #kHz, recoil energy of strontium atoms (as an angular frequency)
    
    # Lamb-Dicke parameters in a Cartesian basis (x, y, z)
    ηα0 = @. sqrt(νR0/να0) # [0.2093, 0.2775, 0.1853]
    
    # Unitless versions
    ρf0_ul = ρf0/λ0 #unitless version of ρf0, 0.1669
    ρa0_ul = ρa0/λ0 #unitless version of ρa0, 0.3991
    a0_ul  = a0/λ0  #unitless version of a0 , 0.2903
    να0_ul = να0/γ0 #unitless version of να0, [14.7297, 8.3784, 18.7838]
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (18.0, 20.0, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = ηα0 * 0.8
    
    # Whether phonons are excluded or not from the calculations (a finite ηα but noPhonons = true will result in including ground state motion into tildeG)
    noPhonons = all(ηα .== 0)
    # noPhonons = true
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set number of atomic sites 
    N_sites = 50
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 1.0
    pos_unc = 0.0
    # pos_unc = ηα0/ωa
    n_inst = 1
    
    # Generate the array, its description, and the number of atoms
    array, arrayDescription, N = get_array(arrayType, N_sites, ρa0_ul, a0_ul, ff, pos_unc, n_inst)
    
    # Time span and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, noPhonons)
    initialStateDescription = "gs"
    
    # Atomic dipole moment
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul, array)
    d = "chiral"
    dDescription = "chiral"
    # d = rightCircularDipoleMoment(array)
    # dDescription = "rgtCrc"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    if typeof(d) == String incField_wlf = [] end
    
    # Whether to include the guided contribution, the radiated contribution, and the radiated interactions when calculating tildeG
    tildeG_flags = (true, true, true)
    
    # Absolute tolerance in the calculations of Im_Grm_trans
    abstol_Im_Grm_trans = 1e-5
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively, usually (true, false))
    approx_Grm_trans = (true, false)
    
    
    # Whether to interpolate Im_Grm_trans
    interpolate_Im_Grm_trans = arrayType == "randomZ"
    
    # Whether to save individual results (Im_Grm_trans, steady states, time evolutions)
    save_Im_Grm_trans = pos_unc == 0 && arrayType != "randomZ" && !interpolate_Im_Grm_trans
    save_steadyState  = n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    save_timeEvol     = n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a0_ul
    margin = 0.1*arrayL
    z_range = range(-margin, arrayL + margin, 100)
    x_range = range(-ρf0_ul - margin, ρf0_ul + ρa0_ul + margin, 100)
    y_fix   = ρa0_ul
    
    # Get the interpolation function for the imaginary, transverse part of the radiation Green's function, if needed
    if interpolate_Im_Grm_trans interpolation_Im_Grm_trans = interpolation1D_Im_Grm_trans(Fiber(ρf0_ul, n0, ωa), Int(ceil(arrayL/0.1)) + 1, ρa0_ul, 0.1, noPhonons) else interpolation_Im_Grm_trans = nothing end
    
    
    return SysPar(ρf0_ul, n0, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N_sites, ρa0_ul, a0_ul, ff, pos_unc, n_inst, array, arrayDescription, N,
                  να0_ul, ηα, noPhonons,
                  d, dDescription, incField_wlf, tildeG_flags, 
                  interpolate_Im_Grm_trans, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                  save_steadyState, save_timeEvol,
                  interpolation_Im_Grm_trans,
                  z_range, x_range, y_fix)
end


function define_SP_artificial()
    # Fiber and array parameters
    n  = 1.45
    ρf = 0.2347
    ρa = 0.6455 - 0.2347 + ρf
    a  = 0.352112676056
    # να = [0.0209, 0.0119, 0.0266]
    να = [5, 12, 7]
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (6.5, 7.5, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = [0.1377, 0.1826, 0.1219] #Cs
    # ηα = [0.2093, 0.2775, 0.1853] #Sr
    ηα = [0.1784, 0.2366, 0.1580]
    # ηα = [0.0, 0.0, 0.0]
    
    # Whether phonons are excluded or not from the calculations (a finite ηα but noPhonons = true will result in including ground state motion into tildeG)
    noPhonons = all(ηα .== 0)
    # noPhonons = false
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set number of atomic sites 
    N_sites = 500
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 1.0
    pos_unc = 0.0
    # pos_unc = ηα0/ωa
    n_inst = 1
    
    # Generate the array, its description, and the number of atoms
    array, arrayDescription, N = get_array(arrayType, N_sites, ρa, a, ff, pos_unc, n_inst)
    
    # Time span and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, noPhonons)
    initialStateDescription = "gs"
    
    # Atomic dipole moment
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul, array)
    d = "chiral"
    dDescription = "chiral"
    # d = rightCircularDipoleMoment(array)
    # dDescription = "rgtCrc"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    if typeof(d) == String incField_wlf = [] end
    
    # Whether to include the guided contribution, the radiated contribution, and the radiated interactions when calculating tildeG
    tildeG_flags = (true, true, true)
    
    # Absolute tolerance in the calculations of Im_Grm_trans
    abstol_Im_Grm_trans = 1e-5
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively, usually (true, false))
    approx_Grm_trans = (true, false)
    
    # Whether to interpolate Im_Grm_trans
    interpolate_Im_Grm_trans = arrayType == "randomZ"
    
    # Whether to save individual results (Im_Grm_trans, steady states, time evolutions)
    save_Im_Grm_trans = pos_unc == 0 && arrayType != "randomZ" && !interpolate_Im_Grm_trans
    save_steadyState  = n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    save_timeEvol     = n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a
    margin = 0.1*arrayL
    z_range = range(-margin, arrayL + margin, 100)
    x_range = range(-ρf - margin, ρf + ρa + margin, 100)
    y_fix   = ρa
    
    # Get the interpolation function for the imaginary, transverse part of the radiation Green's function, if needed
    if interpolate_Im_Grm_trans interpolation_Im_Grm_trans = interpolation1D_Im_Grm_trans(Fiber(ρf, n, ωa), Int(ceil(arrayL/0.1)) + 1, ρa, 0.1, noPhonons) else interpolation_Im_Grm_trans = nothing end
    
    
    return SysPar(ρf, n, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N_sites, ρa, a, ff, pos_unc, n_inst, array, arrayDescription, N,
                  να, ηα, noPhonons,
                  d, dDescription, incField_wlf, tildeG_flags, 
                  interpolate_Im_Grm_trans, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                  save_steadyState, save_timeEvol,
                  interpolation_Im_Grm_trans,
                  z_range, x_range, y_fix)
end


function main()
    # Define system parameters
    # SP = define_SP_BerlinCs()
    SP = define_SP_BerlinSr()
    # SP = define_SP_artificial()
    # show(SP)
    
    
    # plot_propConst_inOutMom(ωρfn_ranges)
    # plot_coupling_strengths(SP)
    # plot_arrayIn3D(SP)
    # plot_σBαTrajectories_σBαSS(SP)
    # plot_transmission_vs_Δ(SP)
    # plot_imperfectArray_transmission_vs_Δ(SP)
    # plot_compareImperfectArray_transmission_vs_Δ(SP)
    plot_effectiveBetaFactor(SP)
    # plot_effectiveBetaFactor_perfectArray(SP)
    # plot_steadyState_radiation_Efield(SP)
    # plot_radiation_Efield(SP)
    # plot_GnmEigenModes(SP)
    # plot_emissionPatternOfGnmeigenModes(SP)
    # plot_GnmEigenEnergies(SP)
    # plot_GnmFourierTransformed(SP)
    # plot_compareGnmEigenEnergies(SP)
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
    
    
    # ρa_rel = 0.411
    # # d_type = "rgtCrc"
    # d_type = "chiral"
    
    # ρfs = range(0.05, 0.5, 100)
    # ρas = ρfs .+ ρa_rel
    # fibers = Fiber.(ρfs, SP.n, ωa)
    # if d_type == "rgtCrc"
    #     ds = fill([1im, 0, 1]/sqrt(2), length(ρfs))
    # elseif d_type == "chiral"
    #     ds = [chiralDipoleMoment(fiber, ρa) for (fiber, ρa) in zip(fibers, ρas)]
    # end
    
    # κs = [fiber.propagation_constant for fiber in fibers]
    # γs = [get_γs(fiber, d, [ρa, 0, 0], SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans) for (ρa, fiber, d) in zip(ρas, fibers, ds)]
    # βs = [γ_gm/(γ_gm + γ_rm) for (γ_gm, γ_rm) in γs]
    
    # # display(GLMakie.Screen(), lines(ρfs, κs, label="ρa_rel=$ρa_rel, d_type=$d_type"))
    # display(GLMakie.Screen(), lines(ρfs, βs))
    
end


function plot_arrayIn3D(SP)
    fig_arrayIn3D(SP.array, SP.x_range, SP.z_range, SP.ρf)
end


function plot_σBαTrajectories_σBαSS(SP)
    Δ = 0.0
    σBα_SS = calc_steadyState(SP, Δ)
    xTrajectories = calc_timeEvolution(SP, Δ)
    # xTrajectories = calc_timeEvolution_eigenmodes(SP, Δ)
    
    if SP.noPhonons
        times, σTrajectories = prep_times_σTrajectories(xTrajectories, SP.N)
        fig_σTrajectories_σSS(times, σTrajectories, σBα_SS)
    else
        times, σTrajectories, BαTrajectories = prep_times_σBαTrajectories(xTrajectories, SP.N)
        fig_σBαTrajectories_σBαSS(times, σTrajectories, BαTrajectories, σBα_SS...)
    end
end


function plot_transmission_vs_Δ(SP)
    titl = prep_transmission_title(SP)
    σBα_scan = scan_steadyState(SP)
    
    t = calc_transmission.(Ref(SP), σBα_scan)
    # t = scan_transmission_eigenmodes(SP)
    # t = scan_transmission_indepDecay(SP)
    
    # T, tPhase = prep_squaredNorm_phase(t)
    # fig_transmission_vs_Δ(SP.Δ_range, T, tPhase, titl)
    
    T, tPhase, unwrappedPhase, phasePerAtom, phaseSlope = prep_squaredNorm_phase_unwrappedPhase_phasePerAtom_phaseSlope(SP, t)
    # fig_transmission_vs_Δ_phaseDetails(SP.Δ_range, T, tPhase, unwrappedPhase, phasePerAtom, phaseSlope, titl)
    # fig_transmission_polar(SP.Δ_range, t, titl)
    fig_transmission_vs_Δ_phaseDetails_polar(SP.Δ_range, T, tPhase, t, unwrappedPhase, phaseSlope, titl)
    
    # r = calc_reflection.(Ref(SP), σBα_scan)
    # R, rPhase = prep_squaredNorm_phase(r)
    # fig_transmissionAndReflection_vs_Δ(SP.Δ_range, T, tPhase, R, rPhase, titl)
    
    # titl = L"$ N = %$(SP.N) $"
    # fig = fig_presentation_transmission_vs_Δ(SP.Δ_range, T, phase, titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\transmission_N500_linear_negative_zoom1.png", fig, px_per_unit=2)

end


function plot_imperfectArray_transmission_vs_Δ(SP)
    if SP.n_inst == 1 throw(ArgumentError("plot_imperfectArray_transmission_vs_Δ requires n_inst > 1")) end
    
    postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, SP.arrayDescription, SP.fiber.postfix)
    filename = "T_phase_" * postfix
    folder = "imperfectArray_T_phase/"
    
    if isfile(saveDir * folder * filename * ".txt") 
        t_real_means, t_real_stds, t_imag_means, t_imag_stds = eachrow(load_as_txt(saveDir * folder, filename))
    else
        ts = []
        if typeof(SP.d) == String
            for array in SP.array
                σBα_scan = scan_steadyState(SP, SP.d, array)
                push!(ts, calc_transmission.(Ref(SP), σBα_scan, Ref(SP.d), Ref(array)))
            end
        else
            for (d, array) in zip(SP.d, SP.array)
                σBα_scan = scan_steadyState(SP, d, array)
                push!(ts, calc_transmission.(Ref(SP), σBα_scan, Ref(d), Ref(array)))
            end
        end
        
        # Prepare means and standard deviations of (squared) magnitudes and phases
        t_real_means, t_real_stds, t_imag_means, t_imag_stds = prep_imperfectArray_transmission(ts)
        formattedResult = vectorOfRows2Matrix([t_real_means, t_real_stds, t_imag_means, t_imag_stds])
        save_as_txt(formattedResult, saveDir * folder, filename)
    end
    
    T_means, T_stds, phase_means, phase_stds = prep_T_argt_statistics(t_real_means, t_real_stds, t_imag_means, t_imag_stds)
    titl = prep_imperfectArray_transmission_title(SP)
    fig_imperfectArray_transmission_vs_Δ(SP.Δ_range, T_means, T_stds, phase_means, phase_stds, titl)
end


function plot_compareImperfectArray_transmission_vs_Δ(SP)
    if SP.n_inst == 1 throw(ArgumentError("plot_imperfectArray_transmission_vs_Δ requires n_inst > 1")) end
    
    # ff_list = [1.0]
    # ff_list = [0.8, 0.85, 0.9, 0.95, 1.0]
    # ff_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    # ff_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7]
    # ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    # ηαFactor_list = [0.0, 0.1, 0.4, 0.7, 1.0]
    ηαFactor_list = [1.0]
    # Nsites_list = [2000, 1000, 667, 500, 400, 334, 250, 200, 167, 143] # N = 100
    # Nsites_list = [6000, 3000, 2000, 1500, 1200, 1000, 750, 600, 500, 429] # N = 300
    Nsites_list = [ 200   100    67    50    40    34    25    20   17   15
                    800   400   267   200   160   134   100    80   67   58
                   1400   700   467   350   280   234   175   140  117  100
                   2000  1000   667   500   400   334   250   200  167  143
                   4000  2000  1334  1000   800   667   500   400  334  286
                   6000  3000  2000  1500  1200  1000   750   600  500  429
                   8000  4000  2667  2000  1600  1334  1000   800  667  572
                  10000  5000  3334  2500  2000  1667  1250  1000  834  715]
    # Nsites_list = [ 200   100   67   50   40   34   25   20   17   15   13   12   10;
    #                 400   200  134  100   80   67   50   40   34   29   25   23   20;
    #                 600   300  200  150  120  100   75   60   50   43   38   34   30;
    #                 800   400  267  200  160  134  100   80   67   58   50   45   40;
    #                1000   500  334  250  200  167  125  100   84   72   63   56   50;
    #                1200   600  400  300  240  200  150  120  100   86   75   67   60;
    #                1400   700  467  350  280  234  175  140  117  100   88   78   70;
    #                1600   800  534  400  320  267  200  160  134  115  100   89   80;
    #                1800   900  600  450  360  300  225  180  150  129  113  100   90;
    #                2000  1000  667  500  400  334  250  200  167  143  125  112  100]
    Ns = [10, 40, 70, 100, 200, 300, 400, 500]
    # Ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    arrayType_list = ["1Dchain", "randomZ"]
    # for ηαFactor in ηαFactor_list
    for ηαFactor in ηαFactor_list, (i, ff) in enumerate(ff_list)
    # for ff in ff_list
        T_meanss, T_stdss, phase_meanss, phase_stdss = [], [], [], []
        T_indepDecayss, phase_indepDecayss = [], []
        labels = []
        # for ηαFactor in ηαFactor_list
        # for ff in ff_list
        # for (ff, N_sites) in zip(ff_list, Nsites_list)
        for N_sites in Nsites_list[:, i]
        # for arrayType in arrayType_list, ff in ff_list
        # for arrayType in arrayType_list, (ff, N_sites) in zip(ff_list, Nsites_list)
        # for ff in ff_list, ηαFactor in ηαFactor_list
            ηα = SP.ηα * ηαFactor
            pos_unc = SP.pos_unc * ηαFactor
            # arrayDescription = arrayDescript(SP.arrayType, SP.N_sites, SP.ρa, SP.a, ff, pos_unc)
            # arrayDescription = arrayDescript(arrayType, SP.N_sites, SP.ρa, SP.a, ff, pos_unc)
            arrayDescription = arrayDescript(SP.arrayType, N_sites, SP.ρa, SP.a, ff, pos_unc)
            # arrayDescription = arrayDescript(arrayType, N_sites, SP.ρa, SP.a, ff, pos_unc)
            
            postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, arrayDescription, SP.fiber.postfix)
            filename = "T_phase_" * postfix
            folder = "imperfectArray_T_phase/"
        
            if isfile(saveDir * folder * filename * ".txt") 
                push!.([T_meanss, T_stdss, phase_meanss, phase_stdss], prep_T_argt_statistics(eachrow(load_as_txt(saveDir * folder, filename))...))
                # push!(labels, L"$ ff = %$(ff) $, $ ηα = %$(ηαFactor) \cdot ηα0 $")
                # push!(labels, L"$ ff = %$(ff) $")
                push!(labels, L"$ N_{sites} = %$(N_sites) $, $ ff = %$(ff) $")
                # push!(labels, L"$ η_{α} = %$(ηαFactor) \cdot η_{α}^{(0)} $")                
                
                # γ_gm, γ_rm = get_γs(SP)
                # N = Int(floor(N_sites*ff))
                # t_indepDecay = transmission_indepDecay.(SP.Δ_range, γ_gm, γ_rm, N)
                # push!.([T_indepDecayss, phase_indepDecayss], prep_squaredNorm_phase(t_indepDecay))
            else
                throw(ArgumentError("The following file can not be found: " * filename))
            end
        end
        
        
        titl = prep_imperfectArray_transmission_title(SP)
        fig_compareImperfectArray_transmission_vs_Δ(SP.Δ_range, T_meanss, T_stdss, phase_meanss, phase_stdss, labels, titl)
        
        # Δ_index = 100
        # T_means, T_stds, phase_means, phase_stds, T_indepDecays, phase_indepDecays = prep_compareImperfectArray_transmission_vs_X(Δ_index, T_meanss, T_stdss, phase_meanss, phase_stdss, T_indepDecayss, phase_indepDecayss)
        # # titl = titl * "\nηα = $(ηαFactor) * ηα0" * "\nΔ = $(SP.Δ_range[Δ_index])"
        # titl = titl * "\nηα = $(ηαFactor) * ηα0, ff = $(ff)" * "\nΔ = $(SP.Δ_range[Δ_index])"
        # # fig_compareImperfectArray_transmission_vs_X(ff_list, L"Filling fraction $$", T_means, T_stds, phase_means, phase_stds, titl)
        # fig_compareImperfectArray_transmission_vs_X(Ns, L"$ N $", T_means, T_stds, phase_means, phase_stds, T_indepDecays, phase_indepDecays, titl)
        
        # titl = L"$ N = %$(SP.N) $, random $ z_{n} $"
        # titl = L"$ N_{sites} = 1000 $, ordered $ z_{n} $"
        # titl = L"$ N_{sites} = 1000 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # titl = L"$ N = 100 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # labels = ff_list
        # labels = ηαFactor_list
        # fig = fig_presentation_compareImperfectArray_transmission_vs_Δ(SP.Δ_range, T_meanss, T_stdss, phase_meanss, phase_stdss, labels, titl)
        # save("C:\\Users\\Simon\\Downloads\\compareff_N100_random.png", fig, px_per_unit=2)
        
        # titl = L"$ N = 100 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # fig = fig_presentation_compareImperfectArray_transmission_vs_ffOrηα(ff_list, L"Filling fraction $$", T_means, T_stds, phase_means, phase_stds, T_noRadInt, phase_noRadInt, titl)
        # save("C:\\Users\\Simon\\Downloads\\compareff_N100.png", fig, px_per_unit=2)
    end
end


function plot_effectiveBetaFactor(SP)
    if SP.n_inst == 1 throw(ArgumentError("plot_effectiveBetaFactor requires n_inst > 1")) end
    
    # Parameters
    # ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7]
    # Nsites_list = [ 200   100    67    50    40    34    25    20   17   15
    #                 800   400   267   200   160   134   100    80   67   58
    #                1400   700   467   350   280   234   175   140  117  100
    #                2000  1000   667   500   400   334   250   200  167  143
    #                4000  2000  1334  1000   800   667   500   400  334  286
    #                6000  3000  2000  1500  1200  1000   750   600  500  429
    #                8000  4000  2667  2000  1600  1334  1000   800  667  572
    #               10000  5000  3334  2500  2000  1667  1250  1000  834  715]
    # Ns = [10, 40, 70, 100, 200, 300, 400, 500]
    ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    Nsites_list = [ 200   100   67   50   40   34   25   20   17   15   13   12   10;
                    400   200  134  100   80   67   50   40   34   29   25   23   20;
                    600   300  200  150  120  100   75   60   50   43   38   34   30;
                    800   400  267  200  160  134  100   80   67   58   50   45   40;
                   1000   500  334  250  200  167  125  100   84   72   63   56   50;
                   1200   600  400  300  240  200  150  120  100   86   75   67   60;
                   1400   700  467  350  280  234  175  140  117  100   88   78   70;
                   1600   800  534  400  320  267  200  160  134  115  100   89   80;
                   1800   900  600  450  360  300  225  180  150  129  113  100   90;
                   2000  1000  667  500  400  334  250  200  167  143  125  112  100]
    Ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    
    arrayType_list = ["1Dchain", "randomZ"][1:1]
    
    # Get independent decay transmissions
    T_indepDecays = zeros(length(Ns), length(SP.Δ_range))
    phase_indepDecays = deepcopy(T_indepDecays)
    γ_gm, γ_rm = get_γs(SP)
    β_indepDecay = γ_gm/(γ_gm + γ_rm)
    for (i, N) in enumerate(Ns)
        t_indepDecay = transmission_indepDecay.(SP.Δ_range, γ_gm, γ_rm, N)
        # t_indepDecay = scan_transmission_indepDecay(SP)
        T_indepDecays[i, :], phase_indepDecays[i, :] = prep_squaredNorm_phase(t_indepDecay)
    end
    
    # Load pre-calculated transmissions
    t_real_means = zeros(length(arrayType_list), length(ff_list), length(Ns), length(SP.Δ_range))    
    t_real_stds, t_imag_means, t_imag_stds, T_means, T_stds, phase_means, phase_stds = deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means)
    for (i, arrayType) in enumerate(arrayType_list), (j, ff) in enumerate(ff_list), (k, N_sites) in enumerate(Nsites_list[:, j])
        arrayDescription = arrayDescript(arrayType, N_sites, SP.ρa, SP.a, ff, SP.pos_unc)
        postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, arrayDescription, SP.fiber.postfix)
        filename = "T_phase_" * postfix
        folder = "imperfectArray_T_phase/"
    
        if isfile(saveDir * folder * filename * ".txt")
            t_real_means[i, j, k, :], t_real_stds[i, j, k, :], t_imag_means[i, j, k, :], t_imag_stds[i, j, k, :] = eachrow(load_as_txt(saveDir * folder, filename))
            T_means[i, j, k, :], T_stds[i, j, k, :], phase_means[i, j, k, :], phase_stds[i, j, k, :] = prep_T_argt_statistics(t_real_means[i, j, k, :], t_real_stds[i, j, k, :], t_imag_means[i, j, k, :], t_imag_stds[i, j, k, :])
        else
            throw(ArgumentError("The following file can not be found: " * filename))
        end
    end
    
    # Fit effective decay rates
    T_fits     = fill(NaN, length(arrayType_list), length(ff_list), length(Ns), length(SP.Δ_range))
    phase_fits = deepcopy(T_fits)
    β_effs  = fill(NaN, length(arrayType_list), length(ff_list), length(SP.Δ_range))
    Δ_effs  = deepcopy(β_effs)
    for (i, arrayType) in enumerate(arrayType_list), (j, ff) in enumerate(ff_list), (l, Δ) in enumerate(SP.Δ_range)
        if l < 400
            # println("i = ", i, ", j = ", j, ", l = ", l)
            model(x, (β_eff, Δ_eff)) = (1 - 2*β_eff/(1 - 2im*Δ_eff))^x
            t           = t_real_means[i, j, :, l] + 1im*t_imag_means[i, j, :, l]
            tDeviations = t_real_stds[i, j, :, l]  + 1im*t_imag_stds[i, j, :, l]
            # pmin = fitComplexData(Ns, t, model, [β_indepDecay*(1 + 5*ff^2), Δ*(1 + 5*ff^2)]; ydataDeviations=tDeviations)
            pmin = fitComplexData(Ns, t, model, [β_indepDecay*(1 + 5*ff^2), Δ*(1 + 5*ff^2)])
            
            T_fits[i, j, :, l], phase_fits[i, j, :, l] = prep_squaredNorm_phase(model.(Ns, Ref(pmin)))
            β_effs[i, j, l], Δ_effs[i, j, l] = pmin
        end
    end
    
    # Unwrap phase for clarity in plots
    phase_means       = mapslices(x -> unwrapPhase(x, "forcePositiveSlope"), phase_means      , dims=3)
    phase_indepDecays = mapslices(x -> unwrapPhase(x, "forcePositiveSlope"), phase_indepDecays, dims=1)
    phase_fits        = mapslices(x -> unwrapPhase(x, "forcePositiveSlope"), phase_fits       , dims=3)
    
    # Plot example of fits
    Δ_index = 75
    for (i, arrayType) in enumerate(arrayType_list)
        labels = [L" ff$ = %$(ff) $, $ β = %$(round(β_effs[i, j, Δ_index], digits=3)) $" for (j, ff) in enumerate(ff_list)]
        titl = prep_imperfectArray_transmission_title(SP) * "\narrayType: $arrayType" * "\nΔ = $(SP.Δ_range[Δ_index])"
        fig_compareImperfectArray_transmission_vs_N(Ns, T_means[i, :, :, Δ_index], T_stds[i, :, :, Δ_index], phase_means[i, :, :, Δ_index], phase_stds[i, :, :, Δ_index], T_indepDecays[:, Δ_index], phase_indepDecays[:, Δ_index], T_fits[i, :, :, Δ_index], phase_fits[i, :, :, Δ_index], β_indepDecay, labels, titl)
    end
    
    # Plot effective β-factors and detunings
    for (i, arrayType) in enumerate(arrayType_list)
        labels = [L" ff$ = %$(ff) $" for (j, ff) in enumerate(ff_list)]
        titl = prep_imperfectArray_transmission_title(SP) * "\narrayType: $arrayType"
        fig_effectiveβΔ_vs_Δ(SP.Δ_range, β_effs[i, :, :], Δ_effs[i, :, :], β_indepDecay, labels, titl)
    end
end


function plot_effectiveBetaFactor_perfectArray(SP)
    Ns = 200:50:750
    
    # Load pre-calculated steady states and calculate transmissions
    ts = zeros(ComplexF64, length(Ns), length(SP.Δ_range))
    Ts, phases = zeros(size(ts)), zeros(size(ts))
    for (i, N) in enumerate(Ns)
        array = get_array(SP.arrayType, N, SP.ρa, SP.a, 1.0, 0.0)
        arrayDescription = arrayDescript(SP.arrayType, N, SP.ρa, SP.a, 1.0, 0.0)
        postfixes = get_postfix_steadyState.(SP.Δ_range, SP.ΔvariDescription, SP.dDescription, Ref(SP.να), Ref(SP.ηα), Ref(SP.incField_wlf), Ref(SP.tildeG_flags), arrayDescription, SP.fiber.postfix)
        for (j, postfix) in enumerate(postfixes)
            if SP.noPhonons filename = "sigma_" * postfix else filename = "sigmaBalpha_" * postfix end
            folder = "steadyStates/"
            if isfile(saveDir * folder * filename * ".jld2")
                σBα = load_as_jld2(saveDir * folder, filename)
                
                if SP.noPhonons
                    tildeΩ = get_tildeΩs(SP.fiber, SP.d, SP.incField_wlf, array)
                    ts[i, j] = transmission(σBα, tildeΩ, SP.fiber)
                else
                    tildeΩ, tildeΩα = get_tildeΩs(SP.fiber, SP.d, SP.ηα, SP.incField_wlf, array)
                    ts[i, j] = transmission(σBα..., tildeΩ, tildeΩα, SP.fiber)
                end
                Ts[i, j], phases[i, j] = prep_squaredNorm_phase(ts[i, j])
            else
                throw(ArgumentError("The following file can not be found: " * filename))
            end
        end
    end
    
    # Fit effective decay rates
    T_fits          = fill(NaN, length(Ns), length(SP.Δ_range))
    phase_fits      = deepcopy(T_fits)
    β_effs          = fill(NaN, length(SP.Δ_range))
    Δ_effs, As, ϕ0s = deepcopy(β_effs), deepcopy(β_effs), deepcopy(β_effs)
    for (j, Δ) in enumerate(SP.Δ_range)
        t = ts[:, j]
        
        # Fit with β, Δ, overall attenuation, and overall phase as fitting parameters
        model(x, (β_eff, Δ_eff, A, ϕ0)) = (1 - 2*β_eff/(1 - 2im*Δ_eff))^x * A*exp(1im*ϕ0)
        if j == 1
            pmin = fitComplexData(Ns, t, model, [0.9999, 7.5*Δ, 0.999, 5.8])
        else
            pmin = fitComplexData(Ns, t, model, [β_effs[j - 1], Δ_effs[j - 1], As[j - 1], ϕ0s[j - 1]])
        end
        β_effs[j], Δ_effs[j], As[j], ϕ0s[j] = pmin[1], pmin[2], pmin[3], mod2pi(pmin[4])
        
        # # Fit with Δ and overall phase as fitting parameters
        # model(x, (Δ_eff, ϕ0)) = (1 - 2/(1 - 2im*Δ_eff))^x * exp(1im*ϕ0)
        # if j == 1
        #     pmin = fitComplexData(Ns, t, model, [7.5*Δ, 5.8])
        # else
        #     pmin = fitComplexData(Ns, t, model, [Δ_effs[j - 1], ϕ0s[j - 1]])
        # end
        # β_effs[j], Δ_effs[j], As[j], ϕ0s[j] = 1.0, pmin[1], 1.0, mod2pi(pmin[2])
        
        T_fits[:, j], phase_fits[:, j] = prep_squaredNorm_phase(model.(Ns, Ref(pmin)))
    end
    
    # Unwrap phase for clarity in plots
    phases     = mapslices(x -> unwrapPhase(x, "forcePositiveSlope"), phases    , dims=1)
    phase_fits = mapslices(x -> unwrapPhase(x, "forcePositiveSlope"), phase_fits, dims=1)
    
    # Plot example of fits
    Δ_index = 100
    label = [L"$ β_{eff} = %$(round(β_effs[Δ_index], digits=5)) $, $ Δ_{eff} = %$(round(Δ_effs[Δ_index], digits=3)) $, $ A = %$(round(As[Δ_index], digits=5)) $, $ ϕ_{0} = %$(round(ϕ0s[Δ_index], digits=3)) $"]
    # label = [L"$ β = 1.0 $, $ Δ = %$(round(Δ_effs[Δ_index], digits=3)) $, $ ϕ_{0} = %$(round(ϕ0s[Δ_index], digits=3)) $"]
    titl = prep_transmission_title(SP) * "\nΔ = $(SP.Δ_range[Δ_index])"
    fig_transmission_vs_N(Ns, Ts[:, Δ_index], phases[:, Δ_index], T_fits[:, Δ_index], phase_fits[:, Δ_index], label, titl)
    
    # Plot effective β-factors and detunings
    titl = prep_transmission_title(SP)
    fig_effectiveβΔ_vs_Δ_perfectArray(SP.Δ_range, β_effs, Δ_effs, titl)
end


function plot_steadyState_radiation_Efield(SP)
    # Get steady state and its Fourier transform
    # Δ = 0.1775
    Δ = 0.201
    σBα_SS = calc_steadyState(SP, Δ)
    if SP.noPhonons σ_SS = σBα_SS else σ_SS = σBα_SS[1] end
    ks, σ_SS_FT = discFourierTransform(σ_SS, SP.a)
    
    # Get emission intensity pattern for plotting
    E = scan_radiation_Efield(SP, σBα_SS)
    intensity = norm.(E).^2
    
    # Get parameter matrices (only considers the excitation part of the Hilbert space)
    Δvari = get_Δvari(SP.ΔvariDependence, SP.Δvari_args, SP.array)
    tildeG_gm, tildeG_rm = get_tildeGs_split(SP)
    tildeG = tildeG_gm + tildeG_rm
    
    # Find eigenmodes etc.
    eigenEnergies, dominant_ks, eigenModesMatrix, eigenModesMatrix_inv = spectrum_dominant_ks_basisMatrices(Δvari + tildeG, SP.a)
    σ_SS_eigenModes = eigenModesMatrix_inv*σ_SS
    
    # Get separate expectation values and total eigenvalues
    collΔ_Δv, collΓ_Δv, collΔ_gm, collΓ_gm, collΔ_rm, collΓ_rm = splitCollEnergies(eigenModesMatrix, Δvari, tildeG_gm, tildeG_rm)
    
    zs = [site[3] for site in SP.array]
    titl = prep_state_title(SP, Δ)
    fig_state(zs, σ_SS, ks, σ_SS_FT, dominant_ks, collΓ_gm, collΓ_rm, σ_SS_eigenModes, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, SP.fiber.propagation_constant, titl)
    
    # fig = fig_presentation_state(zs, σ_SS, ks, σ_SS_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, SP.fiber.propagation_constant, titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\steadyState_N500_linear_positive.png", fig, px_per_unit=2)
end


function plot_radiation_Efield(SP)
    Δ = 0.0
    σ_SS = calc_steadyState(SP, Δ)
    E = scan_radiation_Efield(SP, σ_SS)
    intensity = norm.(E).^2
    
    fig_radiation_Efield(SP.z_range, SP.x_range, intensity, SP.ρf, SP.array)
    
    # fig = fig_presentation_radiation_Efield(SP.z_range, SP.x_range, intensity, SP.ρf, SP.array)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\radiation_Efield_N100_ff0.9.png", fig, px_per_unit=2)
end


function plot_GnmEigenModes(SP)
    zs = [site[3] for site in SP.array]
    
    if SP.noPhonons
        # Get coupling matrix and its spectrum
        Δvari = get_Δvari(SP.ΔvariDependence, SP.Δvari_args, SP.array)
        tildeG_gm, tildeG_rm = get_tildeGs_split(SP)
        tildeG = tildeG_gm + tildeG_rm
        eigenEnergies, dominant_ks, eigenModesMatrix, eigenModesMatrix_inv = spectrum_dominant_ks_basisMatrices(Δvari + tildeG, SP.a)
        eigen_σs = eachcol(eigenModesMatrix)
        collΔ_Δvs, collΓ_Δvs, collΔ_gms, collΓ_gms, collΔ_rms, collΓ_rms = splitCollEnergies(eigenModesMatrix, Δvari, tildeG_gm, tildeG_rm)
        collΔs = collΔ_Δvs + collΔ_gms + collΔ_rms
    
        # Pack the eigenmodes
        eigen_σs_FT = discFourierTransform.(collect.(eigen_σs), SP.a)
        
        # Prepare iteration list to facilitate sorting according to dominant_ks
        iter_list = collect(zip(eigen_σs, eigen_σs_FT, collΔs, collΓ_gms, collΓ_rms, dominant_ks))[sortperm(dominant_ks)]
        
        # Plot
        for (eigen_σ, (ks, eigen_σ_FT), collΔ, collΓ_gm, collΓ_rm, dom_k) in iter_list
            # if collΓ < 10^-2.7
            # if -6.3 < dom_k < -5.5
                E = scan_radiation_Efield(SP, eigen_σ)
                intensity = norm.(E).^2
                # intensity = zeros(size(SP.r_fields))
                
                titl = prep_GnmEigenModes_title(SP)
                fig_GnmEigenModes(zs, eigen_σ, ks, eigen_σ_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, collΔ, collΓ_gm, collΓ_rm, SP.fiber.propagation_constant, titl)
                
                # titl = L"$ \tilde{Δ}_{i}/γ_{a} = %$(round(collΔ, digits=2)) $, $ \tilde{Γ}_{i}/γ_{a} = %$(round(collΓ, digits=4)) $, dominant $ λ_{a}k_z = %$(round(dom_k, digits=2)) $"
                # fig = fig_presentation_GnmEigenModes(zs, eigen_σ, ks, eigen_σ_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, eigenEnergy, SP.fiber.propagation_constant, titl)
                # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\GnmEigenmode_N1000_$(round(dom_k, digits=2)).png", fig, px_per_unit=2)
            # end
        end
        
    else
        # Get coupling matrix and its spectrum
        Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2 = get_parameterMatrices(SP)
        fullCoupling = get_fullCouplingMatrix(Δvari, tildeG, tildeFα, tildeGα1, tildeGα2)
        eigenEnergies, eigenModes = spectrum(fullCoupling)
        
        # Pack the eigenmodes
        eigen_σBαs = unpack_σBαFromσBαVec.(eigenModes)
        eigen_σs   = [eigen_σBα[1] for eigen_σBα in eigen_σBαs]
        eigen_Bαs  = [eigen_σBα[2] for eigen_σBα in eigen_σBαs]
        eigen_diagBαs = [diag.(eigen_Bα) for eigen_Bα in eigen_Bαs]
        
        # Fourier transform
        eigen_σs_FT      =   discFourierTransform.(eigen_σs, SP.a)
        eigen_diagBαs_FT = [[discFourierTransform(eigen_diagBα[α], SP.a)[2] for α in 1:3] for eigen_diagBα in eigen_diagBαs]
        
        # Plot
        for (eigen_σBα, eigen_σ, (ks, eigen_σ_FT), eigen_diagBα, eigen_diagBα_FT, eigenEnergy) in 
            zip(eigen_σBαs, eigen_σs, eigen_σs_FT, eigen_diagBαs, eigen_diagBαs_FT, eigenEnergies)
            # collect(zip(eigen_σBαs, eigen_σs, eigen_σs_FT, eigen_diagBαs, eigen_diagBαs_FT, eigenEnergies))[1:30]
            collΔ, collΓ = collEnergies(eigenEnergy)
            # if collΓ < 10^-2.7
                E = scan_radiation_Efield(SP, eigen_σBα)
                intensity = norm.(E).^2
                
                titl = prep_GnmEigenModes_title(SP)
                # fig_GnmEigenModes(zs, eigen_σ, eigen_diagBα, ks, eigen_σ_FT, eigen_diagBα_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, collΔ, collΓ_gm, collΓ_rm, SP.fiber.propagation_constant, titl)
            # end
        end
        # for (eigen_Bα, eigenEnergy) in zip(eigen_Bαs, eigenEnergies)
        #     collΔ, collΓ = collEnergies(eigenEnergy)
        #     if collΓ < 10^-2.7
        #         for α in 1:3
        #             eigen_Bα_FT = discFourierTransform(eigen_Bα[α], SP.a)
        #             fig_fOnSquare(zs, eigen_Bα[α], eigen_Bα_FT...)
        #         end
        #     end
        # end
    end
end


function plot_emissionPatternOfGnmeigenModes(SP)
    if SP.noPhonons
        Δvari, tildeΩ, tildeG = get_parameterMatrices(SP)
        eigenEnergies, eigen_σs = spectrum(Δvari + tildeG)
        eigen_σBαs = eigenModes
    else
        Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2 = get_parameterMatrices(SP)
        fullCoupling = get_fullCouplingMatrix(Δvari, tildeG, tildeFα, tildeGα1, tildeGα2)
        eigenEnergies, eigenModes = spectrum(fullCoupling)
        eigen_σBαs = unpack_σBαFromσBαVec.(eigenModes)
        eigen_σs   = [eigen_σBα[1] for eigen_σBα in eigen_σBαs]
        eigen_Bαs  = [eigen_σBα[2] for eigen_σBα in eigen_σBαs]
    end
    
    for eigen_σBα in eigen_σBαs
        E = scan_radiation_Efield.(SP, eigen_σBα)
        intensity = norm.(E).^2
        fig_radiation_Efield(SP.z_range, SP.x_range, intensity, SP.ρf, SP.array)
    end
end


function plot_GnmEigenEnergies(SP)
    # The eigenmodes when including phonons do not have a clear band structure
    # if any(SP.ηα .!= 0) throw(ArgumentError("plot_GnmEigenEnergies is not implemented for the case of including phonons")) end
    # but we can still consider the band structure of just the sigma block of the full coupling matrix, which does experience a second order shift due to ground state motion
    
    # Get parameters (explicitly getting guided and radiated part of Gnm)
    Δvari = get_Δvari(SP.ΔvariDependence, SP.Δvari_args, SP.array)
    if SP.noPhonons
        tildeΩ    = get_tildeΩs(SP.fiber, SP.d, SP.incField_wlf, SP.array)
    else 
        tildeΩ, _ = get_tildeΩs(SP.fiber, SP.d, SP.ηα, SP.incField_wlf, SP.array)
    end
    tildeG_gm, tildeG_rm = get_tildeGs_split(SP)
    tildeG = tildeG_gm + tildeG_rm
    # tildeG = get_tildeG0(SP.fiber, SP.d, SP.array)
    # Δvari, tildeΩ, tildeG = get_parameterMatrices(SP)
    
    # Find eigenmodes etc.
    eigenEnergies, dominant_ks, eigenModesMatrix, eigenModesMatrix_inv = spectrum_dominant_ks_basisMatrices(Δvari + tildeG, SP.a)
    
    # Get separate expectation values and total eigenvalues
    collΔ_Δv, collΓ_Δv, collΔ_gm, collΓ_gm, collΔ_rm, collΓ_rm = splitCollEnergies(eigenModesMatrix, Δvari, tildeG_gm, tildeG_rm)
    collΔ = collΔ_Δv + collΔ_gm + collΔ_rm
    collΓ = collΓ_Δv + collΓ_gm + collΓ_rm
    # collΔ, collΓ = collEnergies(eigenEnergies)
    
    # Get weight and resonances
    weights, resonances = transmission_eigenmodes_weights_resonances(SP.Δ_range, tildeΩ, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, SP.fiber.propagation_constant_derivative)
    weights_abs = abs.(weights)
    
    
    titl = prep_GnmEigenEnergies_title(SP)
    # fig_eigenEnergies_vs_k(dominant_ks, collΔ, collΓ, weights_abs, SP.fiber.propagation_constant, titl) 
    fig_eigenEnergies_vs_k(dominant_ks, collΔ, abs.(collΓ), weights_abs, SP.fiber.propagation_constant, titl) 
    # fig_eigenEnergies_vs_k(dominant_ks, collΔ_gm, collΓ_gm, weights_abs, SP.fiber.propagation_constant, titl * "\nGuided contribution") 
    # fig_eigenEnergies_vs_k(dominant_ks, collΔ_rm, collΓ_rm, weights_abs, SP.fiber.propagation_constant, titl * "\nRadiated contribution") 
    # if SP.ΔvariDependence != "flat"
    #     fig_eigenEnergies_vs_k(dominant_ks, collΔ_Δv, collΓ_Δv, weights_abs, SP.fiber.propagation_constant, titl * "\nΔ-variation contribution") 
    # end
    # beta = collΓ_gm./abs.(collΓ)
    # fig_eigenEnergies_vs_k(dominant_ks, collΔ, beta, weights_abs, SP.fiber.propagation_constant, titl * "\nβ-factor") 
    
    # titl = L"$ N = %$(SP.N) $, including $ η_{α}^2 $ shift"
    # titl = L"$ N = %$(SP.N) $, without guided modes"
    # fig = fig_presentation_eigenEnergies_vs_k(dominant_ks, collΔ, abs.(collΓ), weights_abs, SP.fiber.propagation_constant, titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\band_N200_log_shifted.png", fig, px_per_unit=2)
end


function plot_GnmFourierTransformed(SP)
    # N = SP.N÷2
    N = 1000
    kz_range = range(-π/SP.a, π/SP.a, 1000)
    ρ_field = SP.array[1][1:2]
    ρ_source = ρ_field
    zs = -N*SP.a:SP.a:N*SP.a
    r_source = [ρ_source..., 0]
    r_fields = [[ρ_field..., zn] for zn in zs]
    
    if typeof(SP.d) == String
        if SP.d == "chiral"
            d = chiralDipoleMoment(SP.fiber, SP.ρa)
        else
            throw(ArgumentError("plot_GnmFourierTransformed is not implemented for any String dipole moments other than 'chiral'"))
        end
    else
        d = SP.d[1]
    end
        
    # Guided part
    Ggm_ = Ggm.(Ref(SP.fiber), r_fields, Ref(r_source))
    Ggm_ = 3*π/ωa*(Ref(d').*Ggm_.*Ref(d))
    Ggm_kz = [sum(Ggm_.*exp.(-1im*kz*zs)) for kz in kz_range]
    
    # Radiation part
    Grm_ = Grm.(Ref(SP.fiber), ωa, r_fields, Ref(r_source), Ref((0, 0)), 1, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, Ref(SP.approx_Grm_trans), SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)
    Grm_ = 3*π/ωa*(Ref(d').*Grm_.*Ref(d))
    Grm_kz = [sum(Grm_.*exp.(-1im*kz*zs)) for kz in kz_range]
    
    # Get energies
    collΔ_gm, collΓ_gm = collEnergies(Ggm_kz)
    collΔ_rm, collΓ_rm = collEnergies(Grm_kz)
    collΔ = collΔ_gm + collΔ_rm
    collΓ = collΓ_gm + collΓ_rm
    
    
    titl = prep_GnmEigenEnergies_title(SP)
    titl = "Fourier transformed coupling\n" * titl
    weights_abs = zeros(size(collΔ))
    # fig_eigenEnergies_vs_k(kz_range, collΔ, abs.(collΓ), weights_abs, SP.fiber.propagation_constant, titl) 
    fig_eigenEnergies_vs_k(kz_range, collΔ_gm, collΓ_gm, weights_abs, SP.fiber.propagation_constant, titl * "\nGuided contribution") 
    # fig_eigenEnergies_vs_k(kz_range, collΔ_rm, collΓ_rm, weights_abs, SP.fiber.propagation_constant, titl * "\nRadiated contribution") 
end


function plot_compareGnmEigenEnergies(SP)
    # The eigenmodes when including phonons do not have a clear band structure
    if any(SP.ηα .!= 0) throw(ArgumentError("plot_compareGnmEigenEnergies is not implemented for the case of including phonons")) end
    
    # Prepare finite array spectra
    Ns = [1000, 200, 50]
    collΔs = []
    collΓs = []
    dominant_kss = []
    for N in Ns
        array, _, _ = get_array(SP.arrayType, N, SP.ρa, SP.a, SP.ff, SP.pos_unc, 1)
        d = chiralDipoleMoment(SP.fiber, SP.ρa, array)
        tildeG = get_tildeGs(SP.fiber, d, array, SP.tildeG_flags, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)
    
        eigenEnergies, eigenModes, dominant_ks = spectrum_dominant_ks(tildeG, SP.a)
        collΔ, collΓ = collEnergies(eigenEnergies)
        push!(collΔs, collΔ)
        push!(collΓs, collΓ)
        push!(dominant_kss, dominant_ks)
    end
    
    # Prepare infinite case spectrum
    kz_range = range(-π/SP.a, π/SP.a, 300)
    d = chiralDipoleMoment(SP.fiber, SP.ρa)
    spectrum_infArray = Ref(d').*G0_1DFT.(ωa, SP.a, kz_range).*Ref(d)
    collΔ_inf, collΓ_inf = collEnergies(spectrum_infArray)
    
    titl = prep_GnmEigenEnergies_title(SP)
    fig_compareEigenEnergies_vs_k(dominant_kss, collΔs, collΓs, kz_range, collΔ_inf, collΓ_inf, Ns, SP.fiber.propagation_constant, titl) 
end


function plot_lossWithGnmEigenEnergies(SP)
    σBα_scan = scan_steadyState(SP)
    t = calc_transmission.(Ref(SP), σBα_scan)
    r = calc_reflection.(Ref(SP), σBα_scan)
    
    if SP.noPhonons
        Δvari, tildeΩ, tildeG = get_parameterMatrices(SP)
        drive = tildeΩ
        fullCoupling = Δvari + tildeG
    else
        Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2 = get_parameterMatrices(SP)
        drive = get_fullDriveVector(tildeΩ, tildeΩα)
        fullCoupling = get_fullCouplingMatrix(Δvari, tildeG, tildeFα, tildeGα1, tildeGα2)
    end
    
    eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv = spectrum_basisMatrices(fullCoupling)
    collΔ, collΓ = collEnergies(eigenEnergies)
    
    weights, resonances = transmission_eigenmodes_weights_resonances(SP.Δ_range, drive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, SP.fiber.propagation_constant_derivative)
    loss, resonances_abs, resonances_abs_max, exci_pops = prep_loss_resonances_pops(t, r, resonances, eigenModesMatrix, SP.noPhonons)
    titl = prep_transmissionWithGnmEigenEnergies_title(SP)
    
    # flt = resonances_abs_max .> 0.01
    # resonances_abs = resonances_abs[flt]
    # resonances_abs_max = resonances_abs_max[flt]
    # collΔ = collΔ[flt]
    # collΓ = collΓ[flt]
    # exci_pops = exci_pops[flt]
    
    # fig_transAndResonances_polar(SP.Δ_range, t, resonances, titl)
    fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, abs.(collΓ), resonances_abs_max, exci_pops, titl)
    
    # flt = weights_abs .> 1e-20
    # titl = L"$ γ_{a} = 0.01 \cdot γ_{a}^{(0)} $"
    # fig = fig_presentation_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs[flt], collΔ[flt], collΓ[flt], weights_abs[flt], titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\loss_N10_withPhonons_e1.png", fig, px_per_unit=2)
end





@time begin
    println("\n -- Running main() -- \n")
    main()
    println(" -- -- ")
    println("@time of main():")
end




### Major things TODO:
# Consider a pulse/localized excitation in a very long chain to see their dynamics before they reach the ends of the chain, do they decay before the it hits the end?
    # Use transmission spectrum to predict the dynamics of a pulse that is narrow enough in momentum to live within the range of detuning where the transmission is close to 1 when doing Δvari
    # The pulse, implemented via time-dependent driving, will be spatially very very long in order to fit inside the window of high transmission with non-zero phase
    # It should be spatially compressed in the chain, due to the phase gradient (?)
    # Maybe a pulse will not see the narrow loss resonances, because they take a long time to populate?
    # Somehow calculate the loss of the pulse? I.e. does the pulse stay in the chain more than steady state light..?

# Find parameters that result in a strong effective β-factor
    # Plot independent decay transmission for different "strong" βs
    # Scan/explore system parameters to see if any can bring closer to strong β
    # Effect of propagation constant coming close to light cone?

# Driving phonon states for high N
    # Find out whether these states are driven directly by Ωα or indirectly via single-excitation states (that is, a second order process which never populates the single-excitation modes)
        # High density of states counters the fact the jump/interaction is first (or second?) order in the Lamb-Dicke parameter
    # Find out whether T > 1 is due to bug/error or because of ηα expansion
        # Compare time evolution, steady state via basic variables, and steady state via vectorized variables

# Calculate disordered case by having atoms in a tube, rather than randomly positioned traps
    # Keep traps in x and y, but sample fixed positions along z
    # This corresponds to setting ηz = 0 and sampling random positions, as we already did before

# Calculate β_eff for larger and smaller ηα (do this by changing the trap frequencies)
    # Consider +- 50% for example

# Compare with classical sampling for more cases, including Sr

# Fixed N varying Nsites ordered case, the spectra cross, for small detuning higher ff means lower trans, otherwise higher ff increases trans. Why?
    # Is this an artefact of the resonances or the band edge?

# Understand small features in the transmission spectra for random positioned traps (in figures of Summary document)
    # Transmission goes slightly up for random case compared to ordered case for ff = 0.1
    # A slight kink or shoulder for Δ = 0.6
        # Edge of energy band? But why is it not there for ordered case then?
        # General wobbliness due to number of instantiations?

# How does β look for Sr?


### Minor things TODO:
# Check if the decay rates when having eta nonzero and excluding or including phonons is very different
    # Does the addition of phonons increase or decrease the decay rates
    
# Implement cutting out the y-block from the full coupling matrix when it is anyways analytically decoupled
    # Somewhat complicated though as many places have 'α in 1:3' or similar and would (probably) need to be changed to ensure consistency
    # 
    # Whether to include the y-phonon block (the azimuthal atomic motion) in the full coupling matrix
    # include_y_block = arrayType ∉ ("1Dchain", "randomZ")

# Metrics of interaction strength
    # Read up on slow light (it comes from a lot of phase per atom?), maybe we have it here without three-level/EIT setup?
    # t = exp(χN), look at ratio of imaginary and real part of susceptibility χ (these give the phase and loss, depending on whether an i is included in the exponent)
    # Calculate phase per loss for transmission vs. N?
    # Alternative metrics: Slope of phase close to T=1, phase per loss, phase per atom close to T=1, integral of phase over interval with T=1

# Consider a = 500 nm - make t vs ff plots and talk to Philip because he thought that it explained some change in beta?

# Consider small lattice spacing to see if holes matter less? (For small lattice spacing a single hole is not resolved by the wavelength)

# Band structure quickly converges as N increases, but fraction that is from guided or radiation GF changes a lot?

# Zoom in on resonances in T (for no motion or imperfections) and consider phase shift across these
    # To see what could be achieved in the ideal case (in terms of a narrow resonance with near unity transmission and strong light-matter coupling)

# Make randomZ array but with a minimum distance allowed (use the lattice spacing as minimum distance? Then it would only be different from 1Dchain for low ff)
    # Maybe t vs ff plots will level out/find a plateau 
    # This is relevant for comparing with the ordered case, where there is a minimum distance allowed

# Effect of individual eta - scale each individually to see their effect
    # We previously tentatively concluded that axial eta is the most significant
    # Argue which of the ηα is the most significant by looking at the parameter matrices
        # If certain derivatives of tildeΩ or tildeG are large, the corresponding ηα would have a greater effect
        # With this we can say which aspect of the atomic trap is the most important (radial, azimuthal, or axial trapping)
        # Six curves for Ωnα and Ωnαα, twelve heatmaps for Gnmα1, Gnmα2, Gnmαα11, and Gnmαα22 (potentially Gnmα1 and Gnmα2 will show the same plot and likewise for Gnmαα11 and Gnmαα22 - so only six heatmaps)

# Is it relevant for our case that excitations with a certain kz (if this parameter even makes sense to use)
    # can "jump" to another kz' by creating a phonon with momentum kz - kz'?
    # If the kz excitation mode has an energy that matches that of kz' plus the phonon mode this jump is resonant
    # But we only have momentum at kappa, or?

# One can think of the band structure of sigma block of the interaction matrix, 
    # and then the bands of the bsigma blocks (diagonalizing the blocks without taking interaction between them into account)
    # Maybe this is a helpful way of imagining things?
    # Coupling between blocks is small because it scales with eta
    # The bsigma block is given by the sigma block plus a shift and a Kronecker product with the identity
        # so, naively, the band structure is identical but shifted by the trap frequencies
        # and the states are identical but Kronecker multiplied with any basis of phonon states, fx simply having a phonon at site n
        # Hence, the eigenstates and -energies of the full coupling matrix, will indeed generally just show the original spectrum
        # plus a shifted state for each x,y,z
        # and a small perturbation due to the interaction between the excitation sector and the phonon sector,
        # which is small as it scales with eta

# Figure out the calculation of Ggm in the limit of z1=z2
    # an overall sign depending on the sign of Δz = ±eps?
    # if the limit depends on from which direction you approach it, that is problematic?

# Implement the correct real part of the radiation GF
    # But maybe not relevant since atoms are far-ish from the fiber, such that using vacuum-approximation should be good



### Fixes and new features TODO:
# Implement proper rounding of parameters in postfixes
    # Presently there is a unique (and somewhat complicated) rounding for Δ

# Look at calculation of reflection/emission in other guided modes - is it working properly? Plot all four emissions

# Figure out T>1 error for doubleChain 
    # Something with the coupling..?
    # Or dipole moment cannot be such that both sides of the fiber couple to the same mode?
    # Should be sigma+ above and sigma- below (or the other way around), or the corresponding matching pair for the chiral case
    # The atoms are driven by the guided light, so they enter a state that matches the polarization of that light
    # so they should indeed be sigma+ and sigma- polarized above and below respectively

# Plot the Poynting vector instead of the radiation intensity to see where the light is leaked out of the atomic array

# We include terms that go as eta^3 in terms like tildeG*B_α or tildeGα1*σ
    # Removing them can be messy (the eta^2 contribution in σ is not easily removed...)
    # By construction this error should be negligible
    # Use only zeroth order tildeFα in steadyState? Reduce to second order in eta in other ways (i.e. calculate for eta=0 and eta≠0 to extract exact dependencies)? Time evolution anyways also finds results that are effectively higher order

# Optimize calculation of Im_Grm_trans for the case of classical disorder

# Implement plot_GnmEigenEnergies for the case of including phonons
    # Figure out how to order the eigenvalues of the fullCouplingMatrix
    # real part vs imag part?
    # separate into dominated by σ or Bα? 
    # 2D momentum?

# Implement non-lazy version of get_tildeGs(fiber, d::String... for the case of including phonons? 
    # Presumably significantly faster when exploiting knowledge of which components etc. are actually needed, but also very messy...
    # Not needed for classical disorder calculations, so the need is not so great, since the calculations without classical disorder can exploit the z-translational invariance to reduce number of calculations
    # Possibly only implement optimized calculation of integral
        # Less work, and this is obviously the slow part of the overall calculation

# Implement saving and loading of the parameter matrices?

# Clean up and split up files? physics.jl and utility.jl in particular
    
# When calculating Im_Grm_trans, exploit that some entries appear to be zero depending on derivOrder and α,
    # presumably because of the simple array structure

# When scanning over different variables (detuning, (relative) atomic positions, ...) consider using structs as input for the scans
    # In some optimal way calculate all things needed for the scan
    # Put in a struct and use as input for scan
    # To optimize calculation time when scanning

# Update interface comments and make a description of notation/conventions somewhere

# Make small ω version of fiber_equation work? Implement HomotopyContinuation solution?

# Consider making structs for x-vector and σBα, to make their structure easier to get an overview of

# Consider never using JLD2
    # Invent some packing/unpacking scheme to put any kind of data into real matrices

# Consider using StaticArrays in some places?
    # make structs for 3-vectors and 3,3-matrices

# Consider writing tests...
    # Test/check validity of SP (d is normalized, some parameters are non-negative), which could also serve as a reminder of the implicit assumptions regarding parameters
    # Test that transmission is only a phase for 1 atom and no coupling to radiation spectrum
    # Test that long-time time evolution matches with steady state
    # Test calculations in a limit where it can be done analytically? (One atom...)