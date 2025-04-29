
include("preamble.jl")

    
#================================================
    Main functions
================================================#
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
    ν0_radial    = 109 #kHz, radial atomic trap frequency
    ν0_axial     = 139 #kHz, axial atomic trap frequency
    ν0_azimuthal = 18  #kHz, azimuthal atomic trap frequency (estimated from graph)
    να0 = [ν0_radial, ν0_azimuthal, ν0_axial] #trap frequencies in a Cartesian basis (x, y, z) which matches with (radial, azimuthal, axial) if the position is taken to be on the x-axis
    
    # Recoil energy
    ħ = 1.054e-34                                 #m^2*kg/s, Planck's reduced constant
    cesium_mass = 132.9*1.66e-27                  #kg, mass of cesium-133
    νR0 = ħ*(2π/(λ0*1e-9))^2/(2*cesium_mass)*1e-3 #kHz, recoil energy of cesium atoms
    
    # Lamb-Dicke parameters in a Cartesian basis (x, y, z)
    ηα0 = @. sqrt(νR0/να0)
    
    # Unitless versions
    ρf0_ul = ρf0/λ0 #unitless version of ρf0
    ρa0_ul = ρa0/λ0 #unitless version of ρa0
    a0_ul  = a0/λ0  #unitless version of a0
    να0_ul = να0/γ0 #unitless version of να0
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-30, 30, 300)
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 5)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
    N  = 5
    
    # Lamb-Dicke parameters
    ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    ηα = [0., 0., 0.]
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, all(ηα .== 0))
    initialStateDescription = "gs"
    
    # Atomic dipole moment
    # d = [1im, 0, -1]/sqrt(2)
    # d = chiralDipoleMoment(fiber, ρa)
    d = "chiral"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
    # Whether to approximate real, transverse part of radiation GF
    approx_Re_Grm_trans = true
    
    return SysPar(ρf0_ul, n0, ωa,
                  Δ_specs,
                  tspan, dtmax, initialState, initialStateDescription,
                  N, ρa0_ul, a0_ul,
                  να0_ul, ηα,
                  d, incField_wlf, approx_Re_Grm_trans)
end


function define_SP_Olmos()
    # Fiber specs from "Modified dipole-dipole interactions in the presence of a nanophotonic waveguide"
    λ0 = 852  #nm, guided mode wavelength, transition frequency of cs133
    n  = 1.45 #unitless, index of refraction
    ρf = 250  #nm, Fiber radius
    
    # Unitless fiber radius
    ρf_ul = ρf/λ0 
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-30, 30, 300)
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 5)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
    N  = 5
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
    d = [1, 0, 0]
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1)]
    
    # Whether to approximate real, transverse part of radiation GF
    approx_Re_Grm_trans = true
    
    return SysPar(ρf_ul, n, ωa,
                  Δ_specs,
                  tspan, dtmax, initialState, initialStateDescription,
                  N, ρa, a,
                  να, ηα,
                  d, incField_wlf, approx_Re_Grm_trans)
end


function define_SP_Rauschenbeutel()
    # Fiber specs from "Modified dipole-dipole interactions in the presence of a nanophotonic waveguide"
    λ0 = 852  #nm, guided mode wavelength, transition frequency of cs133
    n  = 1.45 #unitless, index of refraction
    ρf = 250  #nm, Fiber radius
    
    # Unitless fiber radius
    ρf_ul = ρf/λ0
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-10, 10, 100)
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 5)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
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
    incField_wlf = [(1, 1, 1)]
    
    # Whether to approximate real, transverse part of radiation GF
    approx_Re_Grm_trans = true
    
    return SysPar(ρf_ul, n, ωa,
                  Δ_specs,
                  tspan, dtmax, initialState, initialStateDescription,
                  N, ρa, a,
                  να, ηα,
                  d, incField_wlf, approx_Re_Grm_trans)
end


function main()
    # Define system parameters
    # ωρfn_ranges = define_ω_ρf_n_ranges()
    SP = define_SP_BerlinCS()
    # SP = define_SP_Olmos()
    # SP = define_SP_Rauschenbeutel()
    # show(SP)
    
    
    
    # plot_propConst_inOutMom(ωρfn_ranges)
    # plot_coupling_strengths(SP)
    # plot_σBαTrajectories_σBαSS(SP)
    plot_transmission_vs_Δ(SP)
    
    return nothing
end


#================================================
    Generate figures
================================================#
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
    Γrm = get_Γrm.(Ref(SP.fiber), Ref(SP.d), r_range, r_range, SP.approx_Re_Grm_trans)
    x_label = L"$ \rho - \rho_f $"
    y_label = L"$ \Gamma_{rm, nn} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γrm, x_label, y_label)
    
    # Radation mode dissipative interaction as a function of radial distance
    ρ_range = range(SP.ρf, SP.ρf + 1000/852, 100)
    r_range = [[ρ, 0, 0] for ρ in ρ_range]
    Γrm = get_Γrm.(Ref(SP.fiber), Ref(SP.d), Ref(r_range[1]), r_range, SP.approx_Re_Grm_trans)
    x_label = L"$ \rho_2 - \rho_f $"
    y_label = L"$ \Gamma_{rm, 12} $"
    x_range = ρ_range .- SP.ρf
    fig_coupling_vs_x(x_range, Γrm, x_label, y_label)
    
end


function plot_σBαTrajectories_σBαSS(SP)
    Δ = 0.0
    σBα_SS = calc_σBα_steadyState(SP, Δ)
    xTrajectories = timeEvolution(SP, Δ)
    
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
    
    fig_transmission_vs_Δ(SP.Δ_range, t)
end





println("\n -- Running main() -- \n")
@time main()




# TODO list:

# Figure out how to add (super-)titles to figures
    # Possibly define function that takes a plot and adds the title to it, by moving existing plots around
# Related to this, figure out how to define proper, nice-looking, robust layouts of plots

# Calculate the free-space emitted E-field and plot it, to see where the atoms "leak"

# Consider making structs for x-vector and σBα, to make their structure easier to get an overview of

# Split utility.jl according to theme or which file has need of the functions
    # simply code utility/no physics concatenated
    # physics functions according to which file has need of them
    
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
    
# Consider never using JLD2
    # Invent some packing/unpacking scheme to put any kind of data into real matrices

# Consider using StaticArrays in some places?

# Consider writing tests...
    # Test/check validity of SP (d is normalized, some parameters are non-negative), which could also serve as a reminder of the implicit assumptions regarding parameters
    # Test that transmission is only a phase for 1 atom and no coupling to radiation spectrum
    # Test that long-time time evolution matches with steady state
    # Test calculations in a limit where it can be done analytically? (One atom...)