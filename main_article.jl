
include("preamble.jl")
const saveDir = "C:/Users/Simon/Forskning/Data/fiber_array_data/"

using GLMakie #plotting, specialized interactive plots
using StatProfilerHTML #profiling the code to see which parts take the most time to run


# const column_width  = 246
# const page_width    = 510
const column_width  = 320
const page_width    = 650
const figure_height = 250
const fontsize      = 10


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
    # a0  = 410   #nm, atomic array lattice constant
    
    # Trap specs
    ν0_radial    = 2π*109 #kHz, radial atomic trap angular frequency
    ν0_axial     = 2π*139 #kHz, axial atomic trap angular frequency
    ν0_azimuthal = 2π*62  #kHz, azimuthal atomic trap angular frequency (estimate from graph: 18 kHz, but usually half of the others)
    να0 = [ν0_radial, ν0_azimuthal, ν0_axial] #trap frequencies in a Cartesian basis (x, y, z) which matches with (radial, azimuthal, axial) if the position is taken to be on the x-axis
    # να0 *= 0.444 # 0.444 0.64 1.777 4.0
    # να0 = να0 .* [4, 4, 4]
    
    # Recoil energy
    νR0 = 2π*2.0663 #kHz, recoil energy of cesium atoms (as an angular frequency)
    
    # Lamb-Dicke parameters in a Cartesian basis (x, y, z)
    ηα0 = @. sqrt(νR0/να0) #[0.1377, 0.1826, 0.1219]
    
    # Unitless versions
    ρf0_ul = ρf0/λ0 #unitless version of ρf0, 0.2347
    ρa0_ul = ρa0/λ0 #unitless version of ρa0, 0.6455
    a0_ul  = a0/λ0  #unitless version of a0 , 0.3521
    να0_ul = να0/γ0 #unitless version of να0, [0.0209, 0.0119, 0.0266]
    
    # Define the fiber
    fiber = Fiber(ρf0_ul, n0, ωa)
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-1.0, 1.0, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = [0., 0., 0.]
    
    # Whether phonons are excluded or not from the calculations (a finite ηα but noPhonons = true will result in including ground state motion into tildeG)
    noPhonons = all(ηα .== 0)
    # noPhonons = true
    
    # Whether to include a third (metastable) level to facilitate EIT
    include3rdLevel = true
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    # arrayType = "randomZ"
    
    # Set number of atomic sites 
    N_sites = 30
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 1.0
    pos_unc = 0.0
    # pos_unc = ηα0/ωa
    n_inst = 1
    
    # Generate the array, its description, and the number of atoms
    array, arrayDescription, N = get_array(arrayType, N_sites, ρa0_ul, a0_ul, ff, pos_unc, n_inst)
    
    # Which time evolver to use ("OrdinaryDiffEq", "simple")
    whichTimeEvolver = "simple"
    
    # Time span and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Prepare initial state for time evolution, as well as description for postfix
    # initialState = groundstate(N, noPhonons, include3rdLevel)
    # initialStateDescription = "gs"
    initialState = Gaussian_sState(N, array, fiber, sqrt(N)*a0_ul, noPhonons, include3rdLevel)
    initialStateDescription = "Ga"
    
    # Whether to have driving on the g-e transition or not
    ΩDriveOn = false
    
    # Atomic dipole moment
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul, array)
    # dDescription = "chiral"
    d = rightCircularDipoleMoment(array)
    dDescription = "rgtCrc"
    # d = [[1, 0, 0] for site in array]
    # dDescription = "xPol"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
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
    save_steadyState  = false # n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    save_timeEvol     = false # n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a0_ul
    margin = 0.1*arrayL
    z_range = range(-margin, arrayL + margin, 30)
    x_range = range(-ρf0_ul - margin, ρf0_ul + ρa0_ul + margin, 30)
    y_fix   = ρa0_ul
    
    # Get the interpolation function for the imaginary, transverse part of the radiation Green's function, if needed
    if interpolate_Im_Grm_trans interpolation_Im_Grm_trans = interpolation1D_Im_Grm_trans(fiber, Int(ceil(arrayL/0.1)) + 1, ρa0_ul, 0.1, ηα) else interpolation_Im_Grm_trans = nothing end
    
    # Type of control drive the third level transition
    cDriveType = "constant" # "constant", "planeWave"
    cDriveDescription = "cst" # "cst", "plW"
    
    # Detuning of the control drive with respect to the e-s transition
    Δc = 0.0
    
    # Rabi frequency of the control drive with respect to the e-s transition
    Ωc = 0.1
    
    # Additional arguments for the control drive ("planeWave" requires a momentum vector)
    cDriveArgs = (kc = ωa*[-0.9, 0, 0], )
    
    # Lower tolerances of the radiative decay rate and state norm used to stop time evolution for calculation of memory retrieval error 
    radDecayRateAndStateNorm_LowerTol = (1e-6, 0.01)
    
    
    return SysPar(ρf0_ul, n0, ωa,
                  fiber,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  whichTimeEvolver, tspan, dtmax,
                  initialState, initialStateDescription,
                  ΩDriveOn,
                  arrayType, N_sites, ρa0_ul, a0_ul, ff, pos_unc, n_inst, array, arrayDescription, N,
                  να0_ul, ηα, noPhonons,
                  d, dDescription, incField_wlf, tildeG_flags, 
                  interpolate_Im_Grm_trans, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                  save_steadyState, save_timeEvol,
                  interpolation_Im_Grm_trans,
                  z_range, x_range, y_fix,
                  include3rdLevel, cDriveType, cDriveDescription, Δc, Ωc, cDriveArgs,
                  radDecayRateAndStateNorm_LowerTol)
end


function define_SP_BerlinSr()
    # Fiber specs 
    λ0  = 689       #nm, guided mode wavelength, transition frequency of Sr
    ω0  = 2π/λ0     #nm^-1, guided mode angular frequency
    γ0  = 2π*7.4    #kHz, free decay rate of Sr88
    ρf0 = 200       #nm, fiber radius, 115, 150, 100, 150, 200, 250, 300, 350, 400
    n0  = 1.45      #unitless, index of refraction
    
    # Atomic array specs
    ρa0 = ρf0 + 100   #nm, atomic array radial coordinate, 275 (115, 150), ρf0 + 100 (100, 150, 200, 250, 300, 350, 400, 450)
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
    
    
    # # Change trap/phonon frequencies to those of Caesium
    # γ0_Cs  = 2π*5.22e3
    # νR0_Cs = 2π*2.0663
    # να0 *= νR0_Cs/νR0 * γ0/γ0_Cs
    
    
    # Unitless versions
    ρf0_ul = ρf0/λ0 #unitless version of ρf0, 0.1669
    ρa0_ul = ρa0/λ0 #unitless version of ρa0, 0.3991
    a0_ul  = a0/λ0  #unitless version of a0 , 0.2903
    να0_ul = να0/γ0 #unitless version of να0, [14.7297, 8.3784, 18.7838]
    
    # Define the fiber
    fiber = Fiber(ρf0_ul, n0, ωa)
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-2, 2, 1000)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    ηα = [0, 0, 0]
    
    # Whether phonons are excluded or not from the calculations (a finite ηα but noPhonons = true will result in including ground state motion into tildeG)
    # noPhonons = all(ηα .== 0)
    noPhonons = true
    
    # Whether to include a third (metastable) level to facilitate EIT
    include3rdLevel = true
    
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
    
    # Which time evolver to use ("OrdinaryDiffEq", "simple")
    whichTimeEvolver = "OrdinaryDiffEq"
    
    # Time span and maximum time step allowed in time evolution
    tspan = (0, 100)
    dtmax = 0.01
    
    # Whether to have driving on the g-e transition or not
    ΩDriveOn = false
    
    # Atomic dipole moment
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul, array)
    # dDescription = "chiral"
    d = rightCircularDipoleMoment(array)
    dDescription = "rgtCrc"
    # d = [[1, 0, 0] for site in array]
    # dDescription = "xPol"
    
    # Incoming field, described by a set of (w, l, f) corresponding to relative weigth, polarization index, and propagation direction index
    # While the relevant part of the drive naively is given simply its dot product with the dipole moment's orientation, the actual polarization/composition of the drive becomes relevant when you include motion
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    
    # Whether to include the guided contribution, the radiated contribution, and the radiated interactions when calculating tildeG
    tildeG_flags = (true, true, true)
    
    # Absolute tolerance in the calculations of Im_Grm_trans
    abstol_Im_Grm_trans = 1e-7
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively, usually (true, false))
    approx_Grm_trans = (true, false)
    
    # Whether to interpolate Im_Grm_trans
    interpolate_Im_Grm_trans = arrayType == "randomZ"
    
    # Whether to save individual results (Im_Grm_trans, steady states, time evolutions)
    save_Im_Grm_trans = pos_unc == 0 && arrayType != "randomZ" && !interpolate_Im_Grm_trans
    save_steadyState  = false # n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    save_timeEvol     = false # n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a0_ul
    margin = 0.1*arrayL
    z_range = range(-margin, arrayL + margin, 100)
    x_range = range(-ρf0_ul - margin, ρf0_ul + ρa0_ul + margin, 100)
    y_fix   = ρa0_ul
    
    # Get the interpolation function for the imaginary, transverse part of the radiation Green's function, if needed
    if interpolate_Im_Grm_trans interpolation_Im_Grm_trans = interpolation1D_Im_Grm_trans(fiber, Int(ceil(arrayL/0.1)) + 1, ρa0_ul, 0.1, ηα) else interpolation_Im_Grm_trans = nothing end
    
    # Type of control drive the third level transition
    cDriveType = "hyperbolic" # "constant", "planeWave", "hyperbolic"
    cDriveDescription = "hyp" # "cst", "plW", "hyp"
    
    # Detuning of the control drive with respect to the e-s transition
    Δc = 0
    
    # Rabi frequency of the control drive with respect to the e-s transition
    Ωc = 0.005
    
    # Additional arguments for the control drive ("planeWave" requires a momentum vector)
    cDriveArgs = (kc = ωa*[-1, 0, 0], N_sites=N_sites, a=a0_ul)
    
    # Lower tolerances of the radiative decay rate and state norm used to stop time evolution for calculation of memory retrieval error 
    radDecayRateAndStateNorm_LowerTol = (1e-6, 1e-2)
    
    # Prepare initial state for time evolution, as well as description for postfix
    # initialState = groundstate(N, noPhonons, include3rdLevel)
    # initialStateDescription = "gs"
    initialState = Gaussian_sState(N, array, fiber, sqrt(N)*a0_ul, noPhonons, include3rdLevel)
    initialStateDescription = "Ga"
    # initialState = triangle_sState(N, array, fiber, noPhonons, include3rdLevel)
    # initialStateDescription = "tr"
    # initialState = interpolatedOptimalMemoryEigenstate(ΔvariDescription, dDescription, arrayType, N_sites, ρa0_ul, a0_ul, να0_ul, ηα, noPhonons, include3rdLevel, tildeG_flags, fiber, cDriveDescription, Δc, Ωc, cDriveArgs)
    # initialStateDescription = "op"
    
    
    return SysPar(ρf0_ul, n0, ωa,
                  fiber,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  whichTimeEvolver, tspan, dtmax,
                  initialState, initialStateDescription,
                  ΩDriveOn,
                  arrayType, N_sites, ρa0_ul, a0_ul, ff, pos_unc, n_inst, array, arrayDescription, N,
                  να0_ul, ηα, noPhonons,
                  d, dDescription, incField_wlf, tildeG_flags, 
                  interpolate_Im_Grm_trans, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                  save_steadyState, save_timeEvol,
                  interpolation_Im_Grm_trans,
                  z_range, x_range, y_fix,
                  include3rdLevel, cDriveType, cDriveDescription, Δc, Ωc, cDriveArgs,
                  radDecayRateAndStateNorm_LowerTol) 
end


function main()
    # Define system parameters
    # SP = define_SP_BerlinCs()
    SP = define_SP_BerlinSr()
    # show(SP)
    
    
    # fig_evanescent_field_intensity()
    # fig_states_comparison(SP)
    # fig_scaling_noMotion_comparison(SP)
    # fig_error_vs_fiberradius(SP)
    
    return nothing
end


# ================================================
#   Generate figures
# ================================================
function fig_evanescent_field_intensity()
    # System parameters
    incField_wlf = [(1, 1, 1), (1, -1, 1)]
    ρf = 0.2
    n = 1.45
    fiber = Fiber(ρf, n, ωa)
        
    # Define radial positions
    rs = range(-1.5, 1.5, 1000)
    r_fields = [[0, r, 0] for r in rs]
    
    # Calculate field
    En = fill(zeros(ComplexF64, 3), size(r_fields))
    for (w, l, f) in incField_wlf
        En += w*Egm.(Ref(fiber), l, f, r_fields)
    end
    En /= sqrt(sum([abs2(w) for (w, l, f) in incField_wlf], init=0.0))
    
    # Calculate intensity
    intensity = real.(adjoint.(En).*En)
    
    
    # Start figure 
    fig = Figure(size=(200, 200))
    
    # Color
    pink = RGB([255, 102, 165]/255...)
    
    # Make axis
    Axis(fig[1, 1], limits=(0, 1.05*maximum(intensity), extrema(rs)...), 
        # xlabel=L"Evanescent field intensity$$", 
        # ylabel=L"$ \rho/\lambda_{a} $",
        topspinevisible=false,
        rightspinevisible=false)
    
    # Plot
    lines!(intensity, rs, color=pink, linewidth=2)
    hlines!([-ρf, ρf], color=:black, linewidth=1, linestyle=:dash)
    
    # Finish figure
    hidedecorations!(label=false)
    # hidedecorations!(label=false, ticklabels=false, ticks=false)
    display(GLMakie.Screen(), fig)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Article 3, fiber memory\\figures\\evanescent_field_intensity.png", fig, px_per_unit=4)
end


function fig_states_comparison(SP)
    # Requires SP from BerlinSr, with N = 100, noPhonons = true, include3rdLevel = true
    
    # Prepare optimal state
    ϵ_eigvals, ϵ_eigmods = calc_memoryRetrievalErrorMatrixEigenmodes(SP)
    optimal_state = unpack_σvarFromσvarVec(ϵ_eigmods[1], SP.N, SP.noPhonons, SP.include3rdLevel)[2]
    
    # Prepare Gaussian state
    zs = [site[3] for site in SP.array]
    zc = (maximum(zs) - minimum(zs))/2
    w = sqrt(SP.N)*SP.a
    Gaussian_state = GaussianState(SP.N, SP.array, SP.fiber.propagation_constant, zc, w, "e", true, false)[1]
    
    # Prepare triangular state
    triangular_state = triangleState(SP.N, SP.array, SP.fiber.propagation_constant, "e", true, false)[1]
    
    # Gather the states
    states = [optimal_state, Gaussian_state, triangular_state]
    labels = [L"Optimal$$", L"Gaussian$$", L"Triangular$$"]
    
    
    # Start figure 
    fig = Figure(size=(column_width, figure_height), figure_padding=5, fontsize=fontsize)
    
    # Colors
    colors = distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Make titles and axes
    ax1 = Axis(fig[1, 1], xticklabelsvisible=false, ylabel=L"Magnitude$$", xlabelsize=fontsize, xticklabelsize=fontsize, ylabelsize=fontsize, yticklabelsize=fontsize)
    ax2 = Axis(fig[2, 1], xlabel=L"$ z_{n}/L $", ylabel=L"Phase$/2\pi$", xlabelsize=fontsize, xticklabelsize=fontsize, ylabelsize=fontsize, yticklabelsize=fontsize)
    
    # Plot the real and imaginary part separately
    for (state, color, label) in zip(states, colors, labels)
        # lines!(ax1, zs/maximum(zs), abs.(state), color=color)
        # lines!(ax2, zs/maximum(zs), unwrapPhase(angle.(state/state[1]), "forcePositiveSlope")/2π, color=color, label=label)
        scatter!(ax1, zs/maximum(zs), abs.(state), color=color, markersize=3)
        scatter!(ax2, zs/maximum(zs), unwrapPhase(angle.(state/state[1]), "forcePositiveSlope")/2π, color=color, markersize=3)
        scatter!(ax2, [-1], [NaN], color=color, label=label)
    end
    
    # Labels
    text!(fig.scene, (0.03, 0.95, 0), text=L"\textbf{a)}$$", space=:relative, fontsize=fontsize+2)
    text!(fig.scene, (0.03, 0.5, 0), text=L"\textbf{b)}$$", space=:relative, fontsize=fontsize+2)
    
    # Finish figure
    axislegend(ax2, position=:lt, fontsize=fontsize, rowgap=1, patchsize=(15.0f0, 10.0f0))
    display(GLMakie.Screen(), fig)
    save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Article 3, fiber memory\\figures\\states_comparison.png", fig, px_per_unit=4)
end


function fig_scaling_noMotion_comparison(SP)
    # Requires SP from BerlinSr
    
    
    params_list = [("eigbasis", SP.να, zeros(3), ""), ("timeEvol", SP.να, zeros(3), "Ga"), ("timeEvol", SP.να, zeros(3), "tr")]
    labels = [L"Optimal$$", L"Gaussian, transient$$", L"Triangular$$", L"Gaussian$$"]
    # Ns = 10:10:200
    Ns = vcat(10:10:200, 220:20:500)
    # Ns = 20:20:440
    ϵs = zeros(length(params_list), length(Ns))
    
    # Load data
    for (i, params) in enumerate(params_list)
        timeEvol_or_eigbasis, να, ηα, initialStateDescription = params
        for (j, N_sites) in enumerate(Ns)
            arrayDescription = arrayDescript(SP.arrayType, N_sites, SP.ρa, SP.a, SP.ff, SP.pos_unc)
            
            folder = "memoryEfficiency/"
            if timeEvol_or_eigbasis == "timeEvol"
                postfix = get_postfix_memoryEfficiency(SP.ΔvariDescription, SP.dDescription, να, ηα, SP.noPhonons, SP.tildeG_flags, arrayDescription, SP.fiber.postfix, initialStateDescription, SP.tspan, SP.dtmax, SP.radDecayRateAndStateNorm_LowerTol, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
                filename = "memEff_" * postfix
                
                if isfile(saveDir * folder * filename * ".txt")
                    ϵs[i, j] = load_as_txt(saveDir * folder, filename)[1]
                else
                    ϵs[i, j] = NaN
                    # throw(ArgumentError("The following file can not be found: " * filename))
                end
            elseif timeEvol_or_eigbasis == "eigbasis"
                postfix = get_postfix_memoryRetrievalErrorMatrixEigenmodes(SP.ΔvariDescription, SP.dDescription, να, ηα, SP.noPhonons, SP.tildeG_flags, arrayDescription, SP.fiber.postfix, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
                filename_eigvals = "memEff_eigvals_" * postfix
            
                if isfile(saveDir * folder * filename_eigvals * ".txt")
                    ϵs[i, j] = minimum(load_as_txt(saveDir * folder, filename_eigvals))
                else
                    ϵs[i, j] = NaN
                    # throw(ArgumentError("The following file can not be found: " * filename_eigvals))
                end
            end
            if ϵs[i, j] < 0 ϵs[i, j] = NaN end
        end
    end
    
    # Remove outliers from optimal state error
    ϵs[1, 15:16] .= NaN
    
    # Make fits
    model_pol(N, (a, α)) = a*N^α
    model_pol1(N, (a,))  = a/N
    model_pol4(N, (a,))  = a/N^4
    model_exp(N, (a, α)) = a*exp(-α*N)
    model_exp_asymp(N, (a, α, c)) = a*exp(-α*N) + c
    p0_pol  = [1000.0, -4.0]
    p0_pol1 = [0.1]
    p0_pol4 = [1000.0]
    p0_exp  = [0.1, 0.01]
    p0_exp_asymp  = [0.1, 0.01, 1e-1]
    # label_pol(p)  = L", $ %$(round(p[1], digits=2)) * N^{%$(round(p[2], digits=2))} $"
    # label_pol1(p) = L", $ %$(round(p[1], digits=3)) * N^{-1} $"
    # label_pol4(p) = L", $ %$(round(p[1], digits=3)) * N^{-4} $"
    # label_exp(p)  = L", $ %$(round(p[1], digits=2)) * e^{-%$(round(p[2], digits=2)) * N} $"
    label_pol(p)  = L", $ \propto N^{%$(round(p[2], digits=2))} $"
    label_pol1(p) = L", $ \propto N^{-1} $"
    label_pol4(p) = L", $ \propto N^{-4} $"
    label_exp(p)  = L", $ \propto e^{-%$(round(p[2], digits=2)) * N} $"
    label_exp_asymp(p) = L", $ %$(round(p[1], digits=3)) * e^{-%$(round(p[2], digits=3)) * N} %$(formatNumberWithSign(round(p[3], digits=3))) $"
    setup_pol  = [model_pol , p0_pol , label_pol]
    setup_pol1 = [model_pol1, p0_pol1, label_pol1]
    setup_pol4 = [model_pol4, p0_pol4, label_pol4]
    setup_exp  = [model_exp , p0_exp , label_exp]
    setup_exp_asymp  = [model_exp_asymp , p0_exp_asymp , label_exp_asymp]
    
    
    fitting_intervals = [5:14, 5:15, 5:35]
    # setups = [setup_pol, setup_exp, setup_pol]
    setups = [setup_pol4, setup_exp, setup_pol1]
    ϵ_fits = fill(NaN, size(ϵs))
    for i in eachindex(params_list)
        model, p0, label = setups[i]
        pmin = fitComplexData(Ns[fitting_intervals[i]], log.(ϵs[i, fitting_intervals[i]]), (N, p) -> log(abs(model(N, p))), p0)
        ϵ_fits[i, :] = abs.(model.(Ns, Ref(pmin)))
        labels[i] *= label(pmin)
    end
    
    fitting_interval = 25:35
    # model, p0, label = setup_pol
    model, p0, label = setup_pol4
    pmin = fitComplexData(Ns[fitting_interval], log.(ϵs[2, fitting_interval]), (N, p) -> log(abs(model(N, p))), p0)
    tail_fit = abs.(model.(Ns, Ref(pmin)))
    labels[4] *= label(pmin)
    
    # Start figure 
    fig = Figure(size=(column_width, figure_height), figure_padding=5, fontsize=fontsize)
    
    # Colors
    colors = distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Make titles and axes
    ax1 = Axis(fig[1, 1], xlabel=L"$ N $", ylabel=L"Memory error, $ ϵ $", yscale=log10, 
        limits=(0, 520, 10^(-8.5), 10^(1.7)),
        xlabelsize=fontsize, xticklabelsize=fontsize, ylabelsize=fontsize, yticklabelsize=fontsize)
    
    # Plot
    scatter!(ax1, Ns, ϵs[1, :], color=colors[1])
    # scatter!(ax1, Ns[vcat(2:2:20, 21:35)], ϵs[2, vcat(2:2:20, 21:35)], color=colors[2])
    scatter!(ax1, Ns, ϵs[2, :], color=colors[2])
    scatter!(ax1, Ns[vcat(2:2:20, 21:35)], ϵs[3, vcat(2:2:20, 21:35)], color=colors[3])
    lines!(ax1, Ns       , ϵ_fits[1, :]   , color=colors[1], label=labels[1])
    lines!(ax1, Ns[16:35], tail_fit[16:35], color=colors[2], label=labels[4])
    lines!(ax1, Ns       , ϵ_fits[3, :]   , color=colors[3], label=labels[3])
    lines!(ax1, Ns[1:25] , ϵ_fits[2, 1:25], color=:black, label=labels[2])
    
    # Finish figure
    axislegend(ax1, position=:rt, fontsize=fontsize, rowgap=1, patchsize=(15.0f0, 10.0f0))
    display(GLMakie.Screen(), fig)
    save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Article 3, fiber memory\\figures\\scaling_noMotion_comparison.png", fig, px_per_unit=4)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Article 3, fiber memory\\figures\\scaling_noMotion_comparison_variableExponents.png", fig, px_per_unit=4)
end


function fig_error_vs_fiberradius(SP)
    # Requires SP from BerlinSr, with 
    
    # Set parameters
    ηα = zeros(3)
    params_list = [("eigbasis", ""), ("timeEvol", "Ga"), ("timeEvol", "tr")]
    labels = [L"Optimal$$", L"Gaussian$$", L"Triangular$$"]
    N_sites = 140
    ρf_list = [0.1451, 0.1814, 0.2177, 0.2540, 0.2903, 0.3266, 0.3628, 0.3991, 0.4354, 0.4717, 0.508, 0.5443, 0.5806]
    ρf_list_nm = 100:25:400
    ρa_list = [0.2903, 0.3266, 0.3628, 0.3991, 0.4354, 0.4717, 0.508, 0.5443, 0.5806, 0.6168, 0.6531, 0.6894, 0.7257]
    
    # Arrays for storing data
    ϵs = zeros(length(params_list), length(ρf_list))
    ls = zeros(length(ρf_list))
    Γs = zeros(length(ρf_list))
    
    # Load data
    for (i, params) in enumerate(params_list)
        timeEvol_or_eigbasis, initialStateDescription = params
        for (j, (ρf, ρa)) in enumerate(zip(ρf_list, ρa_list))
            fiber = Fiber(ρf, SP.n, ωa)
            array, arrayDescription, N = get_array(SP.arrayType, N_sites, ρa, SP.a, SP.ff, SP.pos_unc, 1)
            
            # Get memory efficiency
            folder = "memoryEfficiency/"
            if timeEvol_or_eigbasis == "timeEvol"
                postfix = get_postfix_memoryEfficiency(SP.ΔvariDescription, SP.dDescription, SP.να, ηα, SP.noPhonons, SP.tildeG_flags, arrayDescription, fiber.postfix, initialStateDescription, SP.tspan, SP.dtmax, SP.radDecayRateAndStateNorm_LowerTol, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
                filename = "memEff_" * postfix

                if isfile(saveDir * folder * filename * ".txt")
                    ϵs[i, j] = load_as_txt(saveDir * folder, filename)[1]
                else
                    throw(ArgumentError("The following file can not be found: " * filename))
                end
            elseif timeEvol_or_eigbasis == "eigbasis"
                postfix = get_postfix_memoryRetrievalErrorMatrixEigenmodes(SP.ΔvariDescription, SP.dDescription, SP.να, ηα, SP.noPhonons, SP.tildeG_flags, arrayDescription, fiber.postfix, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
                filename_eigvals = "memEff_eigvals_" * postfix
            
                if isfile(saveDir * folder * filename_eigvals * ".txt")
                    ϵs[i, j] = minimum(load_as_txt(saveDir * folder, filename_eigvals))
                else
                    throw(ArgumentError("The following file can not be found: " * filename_eigvals))
                end
            end
            if ϵs[i, j] < 0 ϵs[i, j] = NaN end
            
            # if i == 1
            #     # Get the evanescent field length scale
            #     ls[j] = 1/sqrt(fiber.propagation_constant^2 - ωa^2)
                
            #     # Get the collective decay near the propagation constant
            #     Δvari = get_Δvari(SP.ΔvariDependence, SP.Δvari_args, array)
            #     d = rightCircularDipoleMoment(array)
            #     tildeG_gm, tildeG_rm = get_tildeGs_split(fiber, d, array, SP.tildeG_flags[3], SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)
            #     tildeG = tildeG_gm + tildeG_rm
            #     eigenEnergies, dominant_ks, eigenModesMatrix, eigenModesMatrix_inv = spectrum_dominant_ks_basisMatrices(Δvari + tildeG, SP.a)
            #     collΔ_Δv, collΓ_Δv, collΔ_gm, collΓ_gm, collΔ_rm, collΓ_rm = splitCollEnergies(eigenModesMatrix, Δvari, tildeG_gm, tildeG_rm)
            #     Γs[j] = collΓ_rm[argmin(abs.(dominant_ks .- fiber.propagation_constant))]
            #     # Γs[j] = collΓ_gm[argmin(abs.(dominant_ks .- fiber.propagation_constant))]
            #     # Γs[j] = collΓ_rm[argmin(abs.(dominant_ks .- fiber.propagation_constant))]/collΓ_gm[argmin(abs.(dominant_ks .- fiber.propagation_constant))]
            # end
    
        end
    end
    
    # ls = exp.(-ρa_list./ls)./sqrt.(ρa_list)
    
    # Start figure 
    fig = Figure(size=(column_width, figure_height), figure_padding=10, fontsize=fontsize)
    # fig = Figure(size=(page_width, figure_height), figure_padding=10, fontsize=fontsize)
    
    # Colors
    colors = distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Make titles and axes
    # ax1 = Axis(fig[1, 1], xlabel=L"$ ρ_{f} $ [nm]", ylabel=L"$ Γ_{rm} $", 
        # xlabelsize=fontsize, xticklabelsize=fontsize, ylabelsize=fontsize, yticklabelsize=fontsize)
    # ax2 = Axis(fig[1, 2], xlabel=L"$ ρ_{f} $ [nm]", ylabel=L"$ l $", 
        # xlabelsize=fontsize, xticklabelsize=fontsize, ylabelsize=fontsize, yticklabelsize=fontsize)
    ax3 = Axis(fig[1, 1], xlabel=L"$ ρ_{f} $ [nm]", ylabel=L"Memory error, $ ϵ $", yscale=log10, 
    # ax3 = Axis(fig[1, 3], xlabel=L"$ ρ_{f} $ [nm]", ylabel=L"$ ϵ $", yscale=log10, 
        limits=(nothing, nothing, 10^(-6.5), 10^(0.5)),
        xlabelsize=fontsize, xticklabelsize=fontsize, ylabelsize=fontsize, yticklabelsize=fontsize)
    
    # Plot
    # lines!(ax1, ρf_list_nm, Γs, color=:black)
    # lines!(ax2, ρf_list_nm, ls, color=:black)
    for (ϵ, color, label) in zip(eachrow(ϵs), colors, labels)
        lines!(ax3, ρf_list_nm, ϵ, color=color, label=label)
    end
    
    # Finish figure
    axislegend(ax3, position=:ct, fontsize=fontsize, rowgap=1, patchsize=(15.0f0, 10.0f0))
    display(GLMakie.Screen(), fig)
    save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Article 3, fiber memory\\figures\\error_vs_fiberradius.png", fig, px_per_unit=4)
end





@time begin
    println("\n -- Running main() -- \n")
    main()
    println(" -- -- ")
    println("@time of main():")
end