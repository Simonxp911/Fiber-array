

include("preamble.jl")
const saveDir = "C:/Users/Simon/Forskning/Data/fiber_array_data/"

using GLMakie #plotting, specialized interactive plots

# Random.seed!(1234)


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
    # ηα = ηα0 * 0.4
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
    # d = "chiral"
    # dDescription = "chiral"
    d = rightCircularDipoleMoment(array)
    dDescription = "rgtCrc"
    
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
    # Fiber specs from "Magic-wavelength nanofiber-based two-color dipole trap with sub-λ/2 spacing"
    λ0  = 689       #nm, guided mode wavelength, transition frequency of sr
    ω0  = 2π/λ0     #nm^-1, guided mode angular frequency
    γ0  = 2π*7.4    #kHz, free decay rate of cs133
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
    νR0 = nothing #kHz, recoil energy of strontium atoms (as an angular frequency)
    # νR0 = 2π*3.47 #kHz, recoil energy of strontium atoms (as an angular frequency)
    # THIS VALUE MAY BE WRONG - calculate from transition wavelength and mass (ask which isotope of Strontium is used, maybe 87?)
    
    # Lamb-Dicke parameters in a Cartesian basis (x, y, z)
    ηα0 = @. sqrt(νR0/να0)
    
    # Unitless versions
    ρf0_ul = ρf0/λ0 #unitless version of ρf0
    ρa0_ul = ρa0/λ0 #unitless version of ρa0
    a0_ul  = a0/λ0  #unitless version of a0
    να0_ul = να0/γ0 #unitless version of να0
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-0.2, 0.4, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = ηα0 * 0.0
    
    # Whether phonons are excluded or not from the calculations
    noPhonons = all(ηα .== 0)
    
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


function main()
    # Define system parameters
    SP = define_SP_BerlinCs()
    # SP = define_SP_BerlinSr()
    # show(SP)
    
    
    # plot_propConst_inOutMom(ωρfn_ranges)
    # plot_coupling_strengths(SP)
    # plot_arrayIn3D(SP)
    # plot_σBαTrajectories_σBαSS(SP)
    plot_transmission_vs_Δ(SP)
    # plot_imperfectArray_transmission_vs_Δ(SP)
    # plot_compareImperfectArray_transmission_vs_Δ(SP)
    # plot_effectiveBetaFactor(SP)
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
    titl = prep_transmission_title(SP)
    σBα_scan = scan_steadyState(SP)
    
    t = calc_transmission.(Ref(SP), σBα_scan)
    # t = scan_transmission_eigenmodes(SP)
    # t = scan_transmission_indepDecay(SP)
    
    # T, tPhase = prep_squaredNorm_phase(t)
    # fig_transmission_vs_Δ(SP.Δ_range, T, tPhase, titl)
    
    T, tPhase, unwrappedPhase, phasePerAtom, phaseSlope = prep_squaredNorm_phase_unwrappedPhase_phasePerAtom_phaseSlope(SP, t)
    fig_transmission_vs_Δ_phaseDetails(SP.Δ_range, T, tPhase, unwrappedPhase, phasePerAtom, phaseSlope, titl)
    
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
    ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7]
    # ηαFactor_list = [0.0, 0.1, 0.4, 0.7, 1.0]
    ηαFactor_list = [1.0]
    Nsites_list = [2000, 1000, 667, 500, 400, 334, 250, 200, 167, 143] # N = 100
    # Nsites_list = [6000, 3000, 2000, 1500, 1200, 1000, 750, 600, 500, 429] # N = 300
    Nsites_list = reverse([ 15   17    20    25    34    40    50    67   100    200;
                            58   67    80   100   134   160   200   267   400    800;
                           100  117   140   175   234   280   350   467   700   1400;
                           143  167   200   250   334   400   500   667  1000   2000;
                           286  334   400   500   667   800  1000  1334  2000   4000;
                           429  500   600   750  1000  1200  1500  2000  3000   6000;
                           572  667   800  1000  1334  1600  2000  2667  4000   8000;
                           715  834  1000  1250  1667  2000  2500  3334  5000  10000], dims=2) # N = 10, 40, 70, 100, 200, 300, 400, 500
    Ns = [10, 40, 70, 100, 200, 300, 400, 500]
    arrayType_list = ["1Dchain", "randomZ"]
    # for ηαFactor in ηαFactor_list
    for ηαFactor in ηαFactor_list, (i, ff) in enumerate(ff_list)
        # if i != 10 continue end
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
                push!.([T_meanss, T_stdss, phase_meanss, phase_stdss], eachrow(load_as_txt(saveDir * folder, filename)))
                # push!(labels, L"$ ff = %$(ff) $, $ ηα = %$(ηαFactor) \cdot ηα0 $")
                # push!(labels, L"$ ff = %$(ff) $")
                push!(labels, L"$ N_{sites} = %$(N_sites) $, $ ff = %$(ff) $")
                # push!(labels, L"$ η_{α} = %$(ηαFactor) \cdot η_{α}^{(0)} $")                
                
                γ_gm = 2*diag(imag(get_tildeGs(SP.fiber, SP.d, SP.array[1][1:1], (true, false, false), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)))
                γ_rm = 2*diag(imag(get_tildeGs(SP.fiber, SP.d, SP.array[1][1:1], (false, true, false), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)))
                N = Int(floor(N_sites*ff))
                t_indepDecay = transmission_indepDecay.(SP.Δ_range, γ_gm, γ_rm, N)
                push!.([T_indepDecayss, phase_indepDecayss], prep_squaredNorm_phase(t_indepDecay))
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
    if SP.n_inst == 1 throw(ArgumentError("plot_imperfectArray_transmission_vs_Δ requires n_inst > 1")) end
    
    # Parameters
    ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7]
    Nsites_list = reverse([ 15   17    20    25    34    40    50    67   100    200;
                            58   67    80   100   134   160   200   267   400    800;
                           100  117   140   175   234   280   350   467   700   1400;
                           143  167   200   250   334   400   500   667  1000   2000;
                           286  334   400   500   667   800  1000  1334  2000   4000;
                           429  500   600   750  1000  1200  1500  2000  3000   6000;
                           572  667   800  1000  1334  1600  2000  2667  4000   8000;
                           715  834  1000  1250  1667  2000  2500  3334  5000  10000], dims=2) # N = 10, 40, 70, 100, 200, 300, 400, 500
    Ns = [10, 40, 70, 100, 200, 300, 400, 500]
    arrayType_list = ["1Dchain", "randomZ"]
    
    # Get independent decay transmissions
    T_indepDecays = zeros(length(Ns), length(SP.Δ_range))
    phase_indepDecays = deepcopy(T_indepDecays)
    γ_gm = 2*imag(get_tildeGs(SP.fiber, SP.d, SP.array[1][1:1], (true, false, false), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)[1])
    γ_rm = 2*imag(get_tildeGs(SP.fiber, SP.d, SP.array[1][1:1], (false, true, false), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)[1])
    for (i, N) in enumerate(Ns)
        t_indepDecay = transmission_indepDecay.(SP.Δ_range, γ_gm, γ_rm, N)
        # t_indepDecay = scan_transmission_indepDecay(SP)
        T_indepDecays[i, :], phase_indepDecays[i, :] = prep_squaredNorm_phase(t_indepDecay)
    end
    
    # Load pre-calculated transmissions
    T_means = zeros(length(arrayType_list), length(ff_list), length(Ns), length(SP.Δ_range))
    T_stds, phase_means, phase_stds = deepcopy(T_means), deepcopy(T_means), deepcopy(T_means)
    for (i, arrayType) in enumerate(arrayType_list), (j, ff) in enumerate(ff_list), (k, N_sites) in enumerate(Nsites_list[:, j])
        arrayDescription = arrayDescript(arrayType, N_sites, SP.ρa, SP.a, ff, SP.pos_unc)
        postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, arrayDescription, SP.fiber.postfix)
        filename = "T_phase_" * postfix
        folder = "imperfectArray_T_phase/"
    
        if isfile(saveDir * folder * filename * ".txt")
            T_means[i, j, k, :], T_stds[i, j, k, :], phase_means[i, j, k, :], phase_stds[i, j, k, :] = eachrow(load_as_txt(saveDir * folder, filename))
        else
            throw(ArgumentError("The following file can not be found: " * filename))
        end
    end
    
    # Fit effective decay rates
    T_fits     = fill(NaN, length(arrayType_list), length(ff_list), length(Ns), length(SP.Δ_range))
    phase_fits = deepcopy(T_fits)
    γ_gm_effs  = fill(NaN, length(arrayType_list), length(ff_list), length(SP.Δ_range))
    γ_rm_effs  = deepcopy(γ_gm_effs)
    for (i, arrayType) in enumerate(arrayType_list), (j, ff) in enumerate(ff_list), (l, Δ) in enumerate(SP.Δ_range)
        if l < 100
            model(x, p) = transmission_indepDecay(Δ, p..., x)
            ydata = sqrt.(T_means[i, j, :, l]).*exp.(1im*phase_means[i, j, :, l])
            pmin = fitComplexData(Ns, ydata, model, [γ_gm, γ_rm])
            
            T_fits[i, j, :, l], phase_fits[i, j, :, l] = prep_squaredNorm_phase(model.(Ns, Ref(pmin)))
            γ_gm_effs[i, j, l], γ_rm_effs[i, j, l] = pmin
        end
    end
    βs = γ_gm_effs./(γ_gm_effs + γ_rm_effs)
    
    # Unwrap phase for clarity in plots
    phase_means = mapslices(unwrapPhase, phase_means, dims=4)
    phase_indepDecays = mapslices(unwrapPhase, phase_indepDecays, dims=2)
    phase_fits = mapslices(unwrapPhase, phase_fits, dims=3)
    
    # Plot example of fits
    Δ_index = 75
    for (i, arrayType) in enumerate(arrayType_list)
        labels = [L" ff$ = %$(ff) $, $ β = %$(round(βs[i, j, Δ_index], digits=3)) $" for (j, ff) in enumerate(ff_list)]
        titl = prep_imperfectArray_transmission_title(SP) * "\narrayType: $arrayType" * "\nΔ = $(SP.Δ_range[Δ_index])"
        fig_compareImperfectArray_transmission_vs_N(Ns, labels, T_means[i, :, :, Δ_index], T_stds[i, :, :, Δ_index], phase_means[i, :, :, Δ_index], phase_stds[i, :, :, Δ_index], T_indepDecays[:, Δ_index], phase_indepDecays[:, Δ_index], T_fits[i, :, :, Δ_index], phase_fits[i, :, :, Δ_index], titl)
    end
    
    # Plot β-factors
    β_indepDecay = γ_gm/(γ_gm + γ_rm)
    for (i, arrayType) in enumerate(arrayType_list)
        labels = [L" ff$ = %$(ff) $" for (j, ff) in enumerate(ff_list)]
        titl = prep_imperfectArray_transmission_title(SP) * "\narrayType: $arrayType"
        fig_βfactor_vs_Δ(SP.Δ_range, labels, βs[i, :, :], β_indepDecay, titl)
    end
    
    # Plot effective decay rates
    for (i, arrayType) in enumerate(arrayType_list)
        labels = [L" ff$ = %$(ff) $" for (j, ff) in enumerate(ff_list)]
        titl = prep_imperfectArray_transmission_title(SP) * "\narrayType: $arrayType"
        fig_effectiveDecayRates_vs_Δ(SP.Δ_range, labels, γ_gm_effs[i, :, :], γ_rm_effs[i, :, :], γ_gm, γ_rm, titl)
    end
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
        # eigen_σs_FT = discFourierTransform.(eigen_σs, SP.a)
        eigen_σs_FT = discFourierTransform.(eigen_σs, SP.a)
        
        # Prepare iteration list to facilitate sorting according to dominant_ks
        iter_list = collect(zip(eigen_σs, eigen_σs_FT, collΔs, collΓ_gms, collΓ_rms, dominant_ks))[sortperm(dominant_ks)]
        
        # Plot
        for (eigen_σ, (ks, eigen_σ_FT), collΔ, collΓ_gm, collΓ_rm, dom_k) in iter_list
            # if collΓ < 10^-2.7
            if -6.3 < dom_k < -5.5
                E = scan_radiation_Efield(SP, eigen_σ)
                intensity = norm.(E).^2
                # intensity = zeros(size(SP.r_fields))
                
                titl = prep_GnmEigenModes_title(SP)
                fig_GnmEigenModes(zs, eigen_σ, ks, eigen_σ_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, collΔ, collΓ_gm, collΓ_rm, SP.fiber.propagation_constant, titl)
                
                # titl = L"$ \tilde{Δ}_{i}/γ_{a} = %$(round(collΔ, digits=2)) $, $ \tilde{Γ}_{i}/γ_{a} = %$(round(collΓ, digits=4)) $, dominant $ λ_{a}k_z = %$(round(dom_k, digits=2)) $"
                # fig = fig_presentation_GnmEigenModes(zs, eigen_σ, ks, eigen_σ_FT, SP.z_range, SP.x_range, intensity, SP.ρf, SP.array, eigenEnergy, SP.fiber.propagation_constant, titl)
                # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\GnmEigenmode_N1000_$(round(dom_k, digits=2)).png", fig, px_per_unit=2)
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
        eigen_σs_FT      =   discFourierTransform.(eigen_σs, SP.a)
        eigen_diagBαs_FT = [[discFourierTransform(eigen_diagBα[α], SP.a)[2] for α in 1:3] for eigen_diagBα in eigen_diagBαs]
        
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
    fig_eigenEnergies_vs_k(dominant_ks, collΔ_gm, collΓ_gm, weights_abs, SP.fiber.propagation_constant, titl * "\nGuided contribution") 
    fig_eigenEnergies_vs_k(dominant_ks, collΔ_rm, collΓ_rm, weights_abs, SP.fiber.propagation_constant, titl * "\nRadiated contribution") 
    if SP.ΔvariDependence != "flat"
        fig_eigenEnergies_vs_k(dominant_ks, collΔ_Δv, collΓ_Δv, weights_abs, SP.fiber.propagation_constant, titl * "\nΔ-variation contribution") 
    end
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
    fig_eigenEnergies_vs_k(kz_range, collΔ, abs.(collΓ), weights_abs, SP.fiber.propagation_constant, titl) 
    fig_eigenEnergies_vs_k(kz_range, collΔ_gm, collΓ_gm, weights_abs, SP.fiber.propagation_constant, titl * "\nGuided contribution") 
    fig_eigenEnergies_vs_k(kz_range, collΔ_rm, collΓ_rm, weights_abs, SP.fiber.propagation_constant, titl * "\nRadiated contribution") 
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
    
    loss, weights_abs, resonances_abs = prep_loss_weights_resonances(t, r, weights, resonances)
    titl = prep_transmissionWithGnmEigenEnergies_title(SP)
    # fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, collΓ, weights_abs, titl)
    fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, abs.(collΓ), weights_abs, titl)
    # fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, abs.(collΓ), 2*weights_abs./abs.(collΓ), titl)
    
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



# TODO list:

# Make trans vs. ff plots for right circular dipole moments

# Metric of interaction strength
# Quantum motion regime
# Pulse driving
# Strontium parameters
# Write notes with plots


# Talk with Thomas about:
    # Arno presented an argument that the increase in decay rate due to ground state motion is due to suddenly having "which path"-information
        # The perfect subradiance of subwavelength arrays can be said to be due to not knowing which atom emits what light
        # If you somehow get information of where light is emitted, the perfect destructive interference disappears
        # And now, if some atom occasionally also generates a phonon when it emits a photon you would be able to know which atom emitted, by looking at where the phonon is
        # But I don't understand this "which path" stuff generally, and the calculation of the decay rates including the eta^2 shift does not involve the phonon part of the state space?

# Understand peak in decay rate vs kz
    # Are these modes always "spread" into the light cone and that is why they can radiate?
    # Could a mode that lives entirely outside of the light cone be radiating because of the finite size of the system and thus the breaking of translational symmetry?
    # The plotted decay rates are the TOTAL decay rates - could it be that their decay into free space is small but their decay into the fiber is large?
    # Divide decay rate (and energies?) by taking expectation value of gm or rm part of Gnm with respect to eigenmodes of full Gnm
    # FT Im_Grm_trans in z?
    # Understand kink in decay rate band at negative light cone boundary
    # Plot bands for non-chiral case or at least sigma+ dipole moment 
    # Calculate emitted guided and radiation light through cylindrical surface around a collection of atoms? 
        # Intensity on surface or flux through surface or?
        # See if these quantities are related simply to the expectation value of the Im_Ggm and Im_Grm, i.e. do the split decay rates really tell us about the ratio of light emitted in the guided or radiated modes?
    # Band structure quickly converges as N increases, but fraction that is from guided or radiation GF changes a lot?

# Metrics of interaction strength
    # Slope of phase in area where T is close to 1
    # Total amount of phase accumulated in area where T is close to 1 (look at unfolded phase plots)
    # T = e^(-OD) = e^(4*β_eff*N) - effective definition of β-factor?
        # Here, the β-factor becomes larger when there is a lot of loss?
        # Also, this expression must be only for a certain regime - namely small β, as e^(4*β_eff*N) easily becomes greater than 1 otherwise
    # Consider phase per atom
        # Is the large accumulated phase we have simply due to lots of atoms? 
    # Read up on slow light (it comes from a lot of phase per atom?), maybe we have it here without three-level/EIT setup?
    # Compare with β-factor for a single atom or analytic expressions for independent atoms?
    # Winding number of phase?
    # Compare transmission of N and N+1 atoms - deduce an effective single-atom β-factor from this, and compare with actual β-factor?
        # Derive how transmission for independent decay is power of single atom transmission, and how this transmission depends on β - from this we can define a β from the ratio of transmissions for N and N+1 atoms
    
        # See what factor you get each time you add an atom
        # For independent atoms you get a constant factor on the transmission for each atom
        # Calculate that number when including the motion and ff and show that the number is better than independent case
        # I.e. you have less loss per phase
        # Fix ff and change number of sites, make transmission vs. ff as before
        # t = exp(χN)
        # Fix ff, change N, see if χ converges to something
        # Look at ratio of imaginary and real part (these give the phase and loss)
        # Look at independent case whether this ratio is "best" on resonance Δ = 0 (it's not. at Δ = 0, the phase is zero. best ratio is found slightly off resonance)
        # Compare with independent case, analytically
        # Fit transmission (inlc. motion and ff) with the expression from independent case to get a effective beta
        # Try to plot χ vs detuning instead of t?
    
    # ... Compare with independent atoms

# Consider a = 500 nm - make t vs ff plots and talk to Philip because he thought that it explained some change in beta?

# Try to calculate with circular dipole moment (circular in xz plane) - it's closer to the experiment
    # Also calculate reflection (emission in all four guided modes)

# Consider small lattice spacing to see if holes matter less? (For small lattice spacing a single hole is not resolved by the wavelength)

# Zoom in on resonances in T (for no motion or imperfections) and consider phase shift across these
    # To see what could be achieved in the ideal case (in terms of a narrow resonance with near unity transmission and strong light-matter coupling)

# Make randomZ array but with a minimum distance allowed (use the lattice spacing as minimum distance? Then it would only be different from 1Dchain for low ff)
    # Maybe t vs ff plots will level out/find a plateau 
    # This is relevant for comparing with the ordered case, where there is a minimum distance allowed

# Look at modes of Gnm with only eta^2 corrections to Hamiltonian but without bsigma part of Hilbert space (find eigenbasis of only first block of coupling matrix with eta^2 included)
    # Make notes about the pleateu'ing of the decay rate
    # See how plateau scales with eta^2? I guess it's obvious that is scales as eta^2

# Look at gamma much smaller than trap frequencies (presently we have gamma 50 times greater than the traps)
    # With gamma much greater than trap frequencies, the atoms decay before they move, so it should match with classical disorder (which is indeed what we see)
    # If gamma is smaller we will see more quantum effects due to motion, and ther should be a difference between the phonon calculation and simply including classical disorder
    # Make notes about the phonon modes moving away and becoming irrelevant due to weak driving

# Plot total "population" in sigma part of wave function and in bsigma part of wave function - find some metric for how much weight each wave function has
    # plot this weight when doing loss and eigenmodes - maybe color markers of decay rate or have another line of plotting

# When looking at loss spectrum, and we have extra significant peaks appearing when including phonons, 
    # these are exactly polaron modes which actually have an effect on the transmission
    # whereas for smaller gamma the phonon modes (polaron modes) are pushed far away

# Effect of individual eta - scale each individually to see their effect
    # We previously tentatively concluded that axial eta is the most significant
    # Argue which of the ηα is the most significant by looking at the parameter matrices
        # If certain derivatives of tildeΩ or tildeG are large, the corresponding ηα would have a greater effect
        # With this we can say which aspect of the atomic trap is the most important (radial, azimuthal, or axial trapping)
        # Six curves for Ωnα and Ωnαα, twelve heatmaps for Gnmα1, Gnmα2, Gnmαα11, and Gnmαα22 (potentially Gnmα1 and Gnmα2 will show the same plot and likewise for Gnmαα11 and Gnmαα22 - so only six heatmaps)

# Consider a pulse/localized excitation in a very long chain to see their dynamics before they reach the ends of the chain, do they decay before the it hits the end?
    # Use transmission spectrum to predict the dynamics of a pulse that is narrow enough in momentum to live within the range of detuning where the transmission is close to 1 when doing Δvari
    # The pulse, implemented via time-dependent driving, will be spatially very very long in order to fit inside the window of high transmission with non-zero phase
    # It should be spatially compressed in the chain, due to the phase gradient (?)
    # Maybe a pulse will not see the narrow loss resonances, because they take a long time to populate?
    # Somehow calculate the loss of the pulse? I.e. does the pulse stay in the chain more than steady state light..?

# Decompose state in terms of eigenmodes
    # Decompose steady state
    # See how the steady state distribution is over momentum, energies, decay rates, and weights
    # Does it make sense?
    # Is it consistent with expectations from Fourier transforming steady state?
    # Look at dynamics of kz components
    # Does the final state ever have significant components within the light cone?
    # Writing state in terms of eigenmodes shows transmission is sum of eigenmode coefficients with weights given by driving and mode overlap
        # Transmission is generally a sum of many significant contributions, giving some unremarkable transmission
        # Changing detuning even slightly can bring in and out a narrow resonance that may thus locally change the transmission dramatically
        # All other contributions are pretty much constant there, just giving some complex background for this single contribution to interfere with
        # Thus even a very narrow, very subradiant, mode can suddenly change the transmission significantly
        # Thus even though the steady as written in terms of the eigenmodes only changes one or two of its entries as the detuning is varied slightly, the transmission may change significantly

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

# Figure out T>1 error for doubleChain 
    # Something with the coupling..?
    # Or dipole moment cannot be such that both sides of the fiber couple to the same mode?
    # Should be sigma+ above and sigma- below (or the other way around), or the corresponding matching pair for the chiral case
    # The atoms are driven by the guided light, so they enter a state that matches the polarization of that light
    # so they should indeed be sigma+ and sigma- polarized above and below respectively

# Find out why transmission phase sometimes has an extra dip/swing..? (i.e. whether the winding number is different for some random instantiations)

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