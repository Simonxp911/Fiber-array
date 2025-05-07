
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
    Δ_specs = (-5, 5, 300)
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
    N = 100
    
    # Lamb-Dicke parameters
    ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    ηα = [0., 0., 0.]
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, all(ηα .== 0))
    initialStateDescription = "gs"
    
    # Atomic dipole moment
    # d = conj([1im, 0, -1]/sqrt(2))
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul)
    d = "chiral"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    # incField_wlf = [(1, 1, 1), (1, -1, 1)]
    incField_wlf = []
    
    # Whether to approximate real, transverse part of radiation GF
    approx_Re_Grm_trans = true
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 1.0
    pos_unc = any(ηα .!= 0) ? 0.0 : 0.0#ηα0/ωa
    n_inst  = any(ηα .!= 0) ?   1 : 1
    
    # Whether to save individual results (Im_Grm_trans, steady states, time evolutions)
    save_individual_res = n_inst == 1
    
    # Ranges of z and x values to define r_field for calculating the radiated E-field
    arrayL = N*a0_ul
    z_range = range(-arrayL, 2*arrayL, 60)
    x_range = range(ρa0_ul - arrayL, ρa0_ul + arrayL, 60)
    
    
    if n_inst == 1
        return SysPar(ρf0_ul, n0, ωa,
                      Δ_specs,
                      tspan, dtmax, initialState, initialStateDescription,
                      N, ρa0_ul, a0_ul, ff, pos_unc,
                      να0_ul, ηα,
                      d, incField_wlf, save_individual_res, approx_Re_Grm_trans,
                      z_range, x_range)
    else
        return [SysPar(ρf0_ul, n0, ωa,
                       Δ_specs,
                       tspan, dtmax, initialState, initialStateDescription,
                       N, ρa0_ul, a0_ul, ff, pos_unc,
                       να0_ul, ηα,
                       d, incField_wlf, save_individual_res, approx_Re_Grm_trans,
                       z_range, x_range) for _ in 1:n_inst]
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
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
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
    Δ_specs = (-30, 30, 100)
    
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


function define_SP_Chang()
    # Fiber specs from "Modified dipole-dipole interactions in the presence of a nanophotonic waveguide"
    n  = 2 #unitless, index of refraction
    ρf = 1.2/(2π)  #unitless, fiber radius
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-10, 10, 300)
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 5)
    dtmax = 0.01
    
    # Set array specs and generate array, as well as description for postfix
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
    
    # Whether to approximate real, transverse part of radiation GF
    approx_Re_Grm_trans = true
    
    return SysPar(ρf, n, ωa,
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
    # SP = define_SP_Chang()
    # show(SP)
    
    
    
    # plot_propConst_inOutMom(ωρfn_ranges)
    # plot_coupling_strengths(SP)
    # plot_σBαTrajectories_σBαSS(SP)
    # plot_transmission_vs_Δ(SP)
    # plot_classDisorder_transmission_vs_Δ(SP)
    plot_E_radiation(SP)
    
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
    
    T, phase = prep_transmission(t)
    fig_transmission_vs_Δ(SP.Δ_range, T, phase)
end


function plot_classDisorder_transmission_vs_Δ(SPs)
    if typeof(SPs) == SysPar throw(ArgumentError("plot_classDisorder_transmission_vs_Δ requires n_inst > 1")) end
    
    n_inst = length(SPs)
    postfix = get_postfix(SPs[1].Δ_specs, SPs[1].d, SPs[1].να, SPs[1].ηα, SPs[1].incField_wlf, n_inst, SPs[1].arrayDescription, SPs[1].fiber.postfix)
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


function plot_E_radiation(SP)
    if any(SP.ηα .!= 0) throw(ArgumentError("plot_E_radiation is not implemented for the case of including phonons")) end
        
    Δ = -2.76
    σ = calc_σBα_steadyState(SP, Δ)
    E = calc_E_radiation.(Ref(SP), Ref(σ), SP.r_field)
    intensity = norm.(E).^2
    
    fig_E_radiation(SP.z_range, SP.x_range, intensity, SP.ρf, SP.array)
end





println("\n -- Running main() -- \n")
@time main()




# TODO list:

# Get it to work on the cluster
    # Use MPI?
    
# Calculate reflection and loss

# Implement n_inst in a better way? That is, dont have a list of SPs, but include the many instantiations in the same SP?

# Check that the derivatives of the modes are correctly implemented by comparing with numerically calculated derivatives

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

# Figure out how to add (super-)titles to figures
    # Possibly define function that takes a plot and adds the title to it, by moving existing plots around
# Related to this, figure out how to define proper, nice-looking, robust layouts of plots

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