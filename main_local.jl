

include("preamble.jl")
const saveDir = "C:/Users/Simon/Forskning/Data/fiber_array_data/"

using GLMakie #plotting, specialized interactive plots

# Random.seed!(1234)


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
    Δ_specs = (-1.0, 1.0, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = ηα0 .* [0.1, 0.2, 0.1]
    # ηα = ηα0 * 0.4
    # ηα = [0.01, 0.01, 0.01]
    # ηα = [0., 0., 0.]
    
    # Whether phonons are excluded or not from the calculations
    noPhonons = all(ηα .== 0)
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set number of atomic sites 
    N_sites = 100
    
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
    margin = 4
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


function define_SP_BerlinCS_mod()
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
    # να0 = να0/0.7^2
    ηα0 = @. sqrt(νR0/να0)
    
    # Unitless versions
    ρf0_ul = ρf0/λ0 #unitless version of ρf0
    ρa0_ul = ρa0/λ0 #unitless version of ρa0
    a0_ul  = a0/λ0  #unitless version of a0
    να0_ul = να0/γ0 #unitless version of να0
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-2.0, 2.0, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = [0., 0., 0.]
    
    # Whether phonons are excluded or not from the calculations
    noPhonons = all(ηα .== 0)
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set number of atomic sites 
    N_sites = 1000
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 0.4
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
    save_steadyState  = n_inst == 1 && pos_unc == 0 && ((ff == 1 && arrayType != "randomZ") || !tildeG_flags[3])
    save_timeEvol     = save_steadyState
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a0_ul
    margin = 4
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
    initialState = groundstate(N, true)
    initialStateDescription = "gs"
     
    # Atomic dipole moment
    d = fill([1, 0, 0], N)
    dDescription = "xPol"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively)
    approx_Grm_trans = (true, false)
    
    return SysPar(ρf_ul, n, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N, ρa, a,
                  να, ηα, true,
                  d, dDescription, incField_wlf, approx_Grm_trans)
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
    initialState = groundstate(N, true)
    initialStateDescription = "gs"
     
    # Atomic dipole moment
    d = fill(conj([1im, 0, -1]/sqrt(2)), N)
    dDescription = "rightCirc"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively)
    approx_Grm_trans = (true, false)
    
    return SysPar(ρf_ul, n, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N, ρa, a,
                  να, ηα, true,
                  d, dDescription, incField_wlf, approx_Grm_trans)
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
    initialState = groundstate(N, true)
    initialStateDescription = "gs"
     
    # Atomic dipole moment
    d = fill([1, 0, 0], N)
    dDescription = "xPol"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively)
    approx_Grm_trans = (true, false)
    
    return SysPar(ρf, n, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N, ρa, a,
                  να, ηα, true,
                  d, dDescription, incField_wlf, approx_Grm_trans)
end


function main()
    # Define system parameters
    # ωρfn_ranges = define_ω_ρf_n_ranges()
    # SP = define_SP_BerlinCS()
    SP = define_SP_BerlinCS_mod()
    # SP = define_SP_Olmos()
    # SP = define_SP_Rauschenbeutel()
    # SP = define_SP_Chang()
    # show(SP)
    
    
    # plot_propConst_inOutMom(ωρfn_ranges)
    # plot_coupling_strengths(SP)
    # plot_arrayIn3D(SP)
    # plot_σBαTrajectories_σBαSS(SP)
    plot_transmission_vs_Δ(SP)
    # plot_imperfectArray_transmission_vs_Δ(SP)
    # plot_compareImperfectArray_transmission_vs_Δ(SP)
    # plot_steadyState_radiation_Efield(SP)
    # plot_radiation_Efield(SP)
    # plot_GnmEigenModes(SP)
    # plot_emissionPatternOfGnmeigenModes(SP)
    # plot_GnmEigenEnergies(SP)
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
end


function plot_coupling_strengths(SP)
    
    # PRESENTLY HARDCODED for define_SP_Rauschenbeutel
    
    # Guided mode local decay as a function of distance to fiber
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γgm = get_Γgm.(Ref(SP.fiber), Ref(SP.d[1]), r_range, r_range)
    x_label = L"$ \rho - \rho_f $"
    y_label = L"$ \Gamma_{gm, nn} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γgm, x_label, y_label)
    
    # Guided mode dissipative interaction as a function of radial distance
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γgm = get_Γgm.(Ref(SP.fiber), Ref(SP.d[1]), Ref(r_range[1]), r_range)
    x_label = L"$ \rho_2 - \rho_f $"
    y_label = L"$ \Gamma_{gm, 12} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γgm, x_label, y_label)
    
    # Radiation mode local decay as a function of distance to fiber
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γrm = get_Γrm.(Ref(SP.fiber), Ref(SP.d[1]), r_range, r_range, SP.approx_Grm_trans)
    x_label = L"$ \rho - \rho_f $"
    y_label = L"$ \Gamma_{rm, nn} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γrm, x_label, y_label)
    
    # Radation mode dissipative interaction as a function of radial distance
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γrm = get_Γrm.(Ref(SP.fiber), Ref(SP.d[1]), Ref(r_range[1]), r_range, SP.approx_Grm_trans)
    x_label = L"$ \rho_2 - \rho_f $"
    y_label = L"$ \Gamma_{rm, 12} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γrm, x_label, y_label)
    
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
    σBα_scan = scan_steadyState(SP)
    t = calc_transmission.(Ref(SP), σBα_scan)
    # t = scan_transmission_eigenmodes(SP)
    
    T, phase = prep_transmission(t)
    titl = prep_transmission_title(SP)
    fig_transmission_vs_Δ(SP.Δ_range, T, phase, titl)
    
    # titl = L"$ N = %$(SP.N) $"
    # fig = fig_presentation_transmission_vs_Δ(SP.Δ_range, T, phase, titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\transmission_N500_linear_negative_zoom1.png", fig, px_per_unit=2)
end


function plot_imperfectArray_transmission_vs_Δ(SP)
    if SP.n_inst == 1 throw(ArgumentError("plot_imperfectArray_transmission_vs_Δ requires n_inst > 1")) end
    
    postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, SP.arrayDescription, SP.fiber.postfix)
    filename = "T_phase" * postfix
    folder = "imperfectArray_T_phase/"
    
    if isfile(saveDir * folder * filename * ".txt") 
        T_means, T_stds, phase_means, phase_stds = eachrow(load_as_txt(saveDir * folder, filename))
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
        T_means, T_stds, phase_means, phase_stds = prep_imperfectArray_transmission(ts)
        formattedResult = vectorOfRows2Matrix([T_means, T_stds, phase_means, phase_stds])
        save_as_txt(formattedResult, saveDir * folder, filename)
    end
    
    titl = prep_imperfectArray_transmission_title(SP)
    fig_imperfectArray_transmission_vs_Δ(SP.Δ_range, T_means, T_stds, phase_means, phase_stds, titl)
end


function plot_compareImperfectArray_transmission_vs_Δ(SP)
    if SP.n_inst == 1 throw(ArgumentError("plot_imperfectArray_transmission_vs_Δ requires n_inst > 1")) end
    
    # ff_list = [1.0]
    # ff_list = [0.8, 0.85, 0.9, 0.95, 1.0]
    # ff_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    # ff_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7][1:1]
    # ηαFactor_list = [0.0, 0.1, 0.4, 0.7, 1.0]
    ηαFactor_list = [1.0]
    Nsites_list = [2000, 1000, 667, 500, 400, 334, 250, 200, 167, 143]
    # Nsites_list = [6000, 3000, 2000, 1500, 1200, 1000, 750, 600, 500, 429]
    arrayType_list = ["1Dchain", "randomZ"]
    for ηαFactor in ηαFactor_list
    # for ff in ff_list
        T_meanss, T_stdss, phase_meanss, phase_stdss = [], [], [], []
        labels = []
        # for ηαFactor in ηαFactor_list
        # for ff in ff_list
        # for (ff, N_sites) in zip(ff_list, Nsites_list)
        for arrayType in arrayType_list, ff in ff_list
        # for arrayType in arrayType_list, (ff, N_sites) in zip(ff_list, Nsites_list)
        # for ff in ff_list, ηαFactor in ηαFactor_list
            ηα = SP.ηα * ηαFactor
            pos_unc = SP.pos_unc * ηαFactor
            # arrayDescription = arrayDescript(SP.arrayType, SP.N_sites, SP.ρa, SP.a, ff, pos_unc)
            arrayDescription = arrayDescript(arrayType, SP.N_sites, SP.ρa, SP.a, ff, pos_unc)
            # arrayDescription = arrayDescript(SP.arrayType, N_sites, SP.ρa, SP.a, ff, pos_unc)
            # arrayDescription = arrayDescript(arrayType, N_sites, SP.ρa, SP.a, ff, pos_unc)
            postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, arrayDescription, SP.fiber.postfix)
            filename = "T_phase" * postfix
            folder = "imperfectArray_T_phase/"
        
            if isfile(saveDir * folder * filename * ".txt") 
                push!.([T_meanss, T_stdss, phase_meanss, phase_stdss], eachrow(load_as_txt(saveDir * folder, filename)))
                # push!(labels, L"$ ff = %$(ff) $, $ ηα = %$(ηαFactor) \cdot ηα0 $")
                push!(labels, L"$ ff = %$(ff) $")
                # push!(labels, L"$ η_{α} = %$(ηαFactor) \cdot η_{α}^{(0)} $")
            else
                throw(ArgumentError("The following file can not be found: " * filename))
            end
        end
        
        titl = prep_imperfectArray_transmission_title(SP)
        fig_compareImperfectArray_transmission_vs_Δ(SP.Δ_range, T_meanss, T_stdss, phase_meanss, phase_stdss, labels, titl)
        
        Δ_index = 1
        T_means, T_stds, phase_means, phase_stds = prep_compareImperfectArray_transmission_vs_ffOrηα(Δ_index, T_meanss, T_stdss, phase_meanss, phase_stdss)
        # titl = titl * "\nηα = $(ηαFactor) * ηα0" * "\nΔ = $(SP.Δ_range[Δ_index])"
        # fig_compareImperfectArray_transmission_vs_ffOrηα(ff_list, L"Filling fraction $$", T_means, T_stds, phase_means, phase_stds, titl)
        
        # titl = L"$ N = %$(SP.N) $, random $ z_{n} $"
        # titl = L"$ N_{sites} = 1000 $, ordered $ z_{n} $"
        titl = L"$ N_{sites} = 1000 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # titl = L"$ N = 100 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # labels = ff_list
        # labels = ηαFactor_list
        # fig = fig_presentation_compareImperfectArray_transmission_vs_Δ(SP.Δ_range, T_meanss, T_stdss, phase_meanss, phase_stdss, labels, titl)
        # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\compareff_eta1.0_NsitesFixed.png", fig, px_per_unit=2)
        
        # fig = fig_presentation_compareImperfectArray_transmission_vs_ffOrηα(ff_list, L"Filling fraction $$", T_means, T_stds, phase_means, phase_stds, titl)
        # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\compareff_eta1.0_NsitesFixed_vs_ff.png", fig, px_per_unit=2)
    end
end


function plot_steadyState_radiation_Efield(SP)
    Δ = -0.5
    zs = [site[3] for site in SP.array]
    σBα_SS = calc_steadyState(SP, Δ)
    if SP.noPhonons σ_SS = σBα_SS else σ_SS = σBα_SS[1] end
    ks, σ_SS_FT = discFourierTransform(σ_SS, SP.a, true, 1000)
    
    E = scan_radiation_Efield(SP, σBα_SS)
    intensity = norm.(E).^2
    
    titl = prep_state_title(SP, Δ)
    fig_state(zs, σ_SS, ks, σ_SS_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, SP.fiber.propagation_constant, titl)
    
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
        Δvari, tildeΩ, tildeG = get_parameterMatrices(SP)
        eigenEnergies, eigen_σs, dominant_ks = spectrum_dominant_ks(Δvari + tildeG, SP.a)
        
        # Pack the eigenmodes
        eigen_σs_FT = discFourierTransform.(eigen_σs, SP.a, true, 1000)
        
        # Prepare iteration list to facilitate sorting according to dominant_ks
        iter_list = collect(zip(eigen_σs, eigen_σs_FT, eigenEnergies, dominant_ks))[sortperm(dominant_ks)]
        
        # Plot
        for (eigen_σ, (ks, eigen_σ_FT), eigenEnergy, dom_k) in iter_list
            # collΔ, collΓ = collEnergies_from_eigenEnergies(eigenEnergy)
            # if collΓ < 10^-2.7
            if 6.4 < dom_k < 7
                E = scan_radiation_Efield(SP, eigen_σ)
                intensity = norm.(E).^2
                # intensity = zeros(size(SP.r_fields))
                
                titl = prep_GnmEigenModes_title(SP)
                fig_GnmEigenModes(zs, eigen_σ, ks, eigen_σ_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, eigenEnergy, SP.fiber.propagation_constant, titl)
            end
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
        eigen_σs_FT      =   discFourierTransform.(eigen_σs, SP.a, true, 1000)
        eigen_diagBαs_FT = [[discFourierTransform(eigen_diagBα[α], SP.a, true, 1000)[2] for α in 1:3] for eigen_diagBα in eigen_diagBαs]
        
        # Plot
        for (eigen_σBα, eigen_σ, (ks, eigen_σ_FT), eigen_diagBα, eigen_diagBα_FT, eigenEnergy) in 
            zip(eigen_σBαs, eigen_σs, eigen_σs_FT, eigen_diagBαs, eigen_diagBαs_FT, eigenEnergies)
            # collect(zip(eigen_σBαs, eigen_σs, eigen_σs_FT, eigen_diagBαs, eigen_diagBαs_FT, eigenEnergies))[1:30]
            collΔ, collΓ = collEnergies_from_eigenEnergies(eigenEnergy)
            if collΓ < 10^-2.7
                E = scan_radiation_Efield(SP, eigen_σBα)
                intensity = norm.(E).^2
                
                titl = prep_GnmEigenModes_title(SP)
                fig_GnmEigenModes(zs, eigen_σ, eigen_diagBα, ks, eigen_σ_FT, eigen_diagBα_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, eigenEnergy, SP.fiber.propagation_constant, titl)
            end
        end
        # for (eigen_Bα, eigenEnergy) in zip(eigen_Bαs, eigenEnergies)
        #     collΔ, collΓ = collEnergies_from_eigenEnergies(eigenEnergy)
        #     if collΓ < 10^-2.7
        #         for α in 1:3
        #             eigen_Bα_FT = discFourierTransform(eigen_Bα[α], SP.a, true, 1000)
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
    if any(SP.ηα .!= 0) throw(ArgumentError("plot_GnmEigenEnergies is not implemented for the case of including phonons")) end
    
    Δvari, tildeΩ, tildeG = get_parameterMatrices(SP)
    # tildeG = get_tildeG0(SP.fiber, SP.d, SP.array)
    
    eigenEnergies, dominant_ks, eigenModesMatrix, eigenModesMatrix_inv = spectrum_dominant_ks_basisMatrices(Δvari + tildeG, SP.a)
    collΔ, collΓ = collEnergies_from_eigenEnergies(eigenEnergies)
    weights, resonances = transmission_eigenmodes_weights_resonances(SP.Δ_range, tildeΩ, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, SP.fiber.propagation_constant_derivative)
    weights_abs = abs.(weights)
    
    titl = prep_GnmEigenEnergies_title(SP)
    fig_eigenEnergies_vs_k(dominant_ks, collΔ, collΓ, weights_abs, SP.fiber.propagation_constant, titl) 
    
    # titl = L"$ N = %$(SP.N) $"
    # titl = L"$ N = %$(SP.N) $, without guided modes"
    # fig = fig_presentation_eigenEnergies_vs_k(dominant_ks, collΔ, collΓ, weights_abs, SP.fiber.propagation_constant, titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\band_N50.png", fig, px_per_unit=2)
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
        collΔ, collΓ = collEnergies_from_eigenEnergies(eigenEnergies)
        push!(collΔs, collΔ)
        push!(collΓs, collΓ)
        push!(dominant_kss, dominant_ks)
    end
    
    # Prepare infinite case spectrum
    kz_range = range(-π/SP.a, π/SP.a, 300)
    d = chiralDipoleMoment(SP.fiber, SP.ρa)
    spectrum_infArray = Ref(d').*G0_1DFT.(ωa, SP.a, kz_range).*Ref(d)
    collΔ_inf, collΓ_inf = collEnergies_from_eigenEnergies(spectrum_infArray)
    
    titl = prep_GnmEigenEnergies_title(SP)
    fig_compareEigenEnergies_vs_k(dominant_kss, collΔs, collΓs, kz_range, collΔ_inf, collΓ_inf, Ns, SP.fiber.propagation_constant, titl) 
end


function plot_lossWithGnmEigenEnergies(SP)
    σBα_scan = scan_steadyState(SP)
    t = calc_transmission.(Ref(SP), σBα_scan)
    
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
    collΔ, collΓ = collEnergies_from_eigenEnergies(eigenEnergies)
    weights, resonances = transmission_eigenmodes_weights_resonances(SP.Δ_range, drive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, SP.fiber.propagation_constant_derivative)
    
    loss, weights_abs, resonances_abs = prep_loss_weights_resonances(t, weights, resonances)
    titl = prep_transmissionWithGnmEigenEnergies_title(SP)
    fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, collΓ, weights_abs, titl)
    # fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, abs.(collΓ), weights_abs, titl)
    # fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, abs.(collΓ), 2*weights_abs./abs.(collΓ), titl)
    
    # flt = weights_abs .> 1e-20
    # fig = fig_presentation_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs[flt], collΔ[flt], collΓ[flt], weights_abs[flt], titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\loss_N10_withPhonons.png", fig, px_per_unit=2)
end





@time begin
    println("\n -- Running main() -- \n")
    main() 
    println(" -- -- ")
    println("@time of main():")
end



# TODO list:

# Set up code to systematically calculate with
    # Only Grm
    # Only Ggm
    # Independent decay (only diagonal from Grm)
# Make sure it is indicated in the postfixes...


# Look at modes of Gnm with only eta^2 corrections to Hamiltonian but without bsigma part of Hilbert space (find eigenbasis of only first block of coupling matrix with eta^2 included)
# When looking at loss spectrum, and we have extra significant peaks when including phonons, these are exactly polaron modes which actually have an effect on the transmission
    # whereas for smaller gamma the phonon modes (polaron modes) are pushed far away
# Plot total "population" in sigma part of wave function and in bsigma part of wave function - find some metric for how much weight each wave function has
    # plot this weight when doing loss and eigenmodes - maybe color markers of decay rate or have another line of plotting
# Effect of individual eta - scale each individually to show that axial is most important

# In slides:
    # Emphasize that even with phonon energies far away, there is still a second order contribution to the parameters (write it out more explicitly)
        # In particular, there is a contribution to the decay rate, such that modes outside the light cone start to decay
    # Show loss and modes for smaller gamma (same eta)
    # Talk about width of transmission features coming from gamma, and whether these widths then include the spread of phonons due to trap energies or not
    # Phonons pushed away is because they are only created by photon-mediated interactions and with gamma small, these are small, such that we dont excite phonons
    # Show which eta is most important?
    # Compare t vs ff for fixed N with result for independent decay (no interaction through radiation)
    


# Look at gamma much smaller than trap frequencies (presently we have gamma 50 times greater than the traps)
    # With gamma much greater than trap frequencies, the atoms decay before they move, so it should match with classical disorder (which is indeed what we see)
    # If gamma is smaller we will see more quantum effects due to motion, and ther should be a difference between the phonon calculation and simply including classical disorder

# Consider a pulse/localized excitation in a very long chain to see their dynamics before they reach the ends of the chain, do they decay before the it hits the end?
    # Use transmission spectrum to predict the dynamics of a pulse that is narrow enough in momentum to live within the range of detuning where the transmission is close to 1 when doing Δvari
    # The pulse, implemented via time-dependent driving, will be spatially very very long in order to fit inside the window of high transmission with non-zero phase
    # It should be spatially compressed in the chain, due to the phase gradient (?)

# Make plots of effect of ff without motion, i.e. eta=0, focus on T = 1 peak right after flat plateau at negative detuning when including Deltavari 
    # This peak is not always there (depending on the number of atoms that is), so it's a bit arbitrary to focus on it
    # These peaks are result of a complicated interference of many modes, so the appearance of such a peak is somewhat random/arbitrary
    
# Decompose state in terms of eigenmodes
    # Decompose steady state
    # See how the steady state distribution is over momentum, energies, decay rates, and weights
    # Does it make sense?
    # Is it consistent with expectations from Fourier transforming steady state?

# The eigenmodes which are due to the guided mode GF are not strongly peaked in kz
    # Does it really make sense to plot them on the kz spectrum as we have?
    # Is it misleading to think of them as concentrated around kz = kappa?

# Is it relevant for our case that excitations with a certain kz (if this parameter even makes sense to use)
# can "jump" to another kz' by creating a phonon with momentum kz - kz'?
    # If the kz excitation mode has an energy that matches that of kz' plus the phonon mode this jump is resonant
    # But we only have momentum at kappa, or?

# One can think of the band structure of sigma block of the interaction matrix, 
    # and then the bands of the bsigma blocks (diagonalizing the blocks without taking interaction between them into account)
    # Maybe this is a helpful way of imagining things?


# Figure out T>1 error for doubleChain 
    # Something with the coupling..?
    # Or dipole moment cannot be such that both sides of the fiber couple to the same mode?


# Be careful that large phase is not just N times a small phase per atom

# FT Im_Grm_trans in z?

# Determine effective β-factor from the phase?
    # Compare phase with that of a single atom in the case of no radiation GF
    # Find analytical expression for phase to see how
    # ... but in that case there is no radiation and thus β is just 1..?
    # Well, something can be done perhaps...
    # Compare with model of 1Dchain with some decay rate into free space, but no interaction through free space (i.e. independent decay into free space, but still interacting through fiber)

# Find out why transmission phase sometimes has an extra dip/swing..? (i.e. whether the winding number is different for some random instantiations)

# Make scans where only one etaalpha is varied away from zero, to see their individual effects and significances

# Argue which of the ηα is the most significant by looking at the parameter matrices
    # If certain derivatives of tildeΩ or tildeG are large, the corresponding ηα would have a greater effect
    # With this we can say which aspect of the atomic trap is the most important (radial, azimuthal, or axial trapping)
    # Six curves for Ωnα and Ωnαα, twelve heatmaps for Gnmα1, Gnmα2, Gnmαα11, and Gnmαα22 (potentially Gnmα1 and Gnmα2 will show the same plot and likewise for Gnmαα11 and Gnmαα22 - so only six heatmaps)

# While high-weight modes give the big overarching contributions to the loss, the actual visible peaks of the loss are given by relatively low-weight low-decay modes
    # What's up with that?

# Implement analytic calculation of steady state and transmission for the case of no radiation interactions

# Plot the Poynting vector instead of the radiation intensity to see where the light is leaked out of the atomic array

# We include terms that go as eta^3 in terms like tildeG*B_α or tildeGα1*σ
    # Removing them can be messy (the eta^2 contribution in σ is not easily removed...)
    # By construction this error should be negligible

# Figure out the calculation of Ggm in the limit of z1=z2
    # an overall sign depending on the sign of Δz = ±eps?
    # if the limit depends on from which direction you approach it, that is problematic?

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

# Implement the correct real part of the radiation GF
    # But maybe not relevant since atoms are far-ish from the fiber, such that using vacuum-approximation should be good
    
# When calculating Im_Grm_trans, exploit that some entries appear to be zero depending on derivOrder and α,
    # presumably because of the simple array structure

# When scanning over different variables (detuning, (relative) atomic positions, ...) consider using structs as input for the scans
    # In some optimal way calculate all things needed for the scan
    # Put in a struct and use as input for scan
    # To optimize calculation time when scanning

# Use only zeroth order tildeFα in steadyState? Reduce to second order in eta in other ways (i.e. calculate for eta=0 and eta≠0 to extract exact dependencies)? Time evolution anyways also finds results that are effectively higher order

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