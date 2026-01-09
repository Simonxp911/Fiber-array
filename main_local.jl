

include("preamble.jl")
const saveDir = "C:/Users/Simon/Forskning/Data/fiber_array_data/"

using GLMakie #plotting, specialized interactive plots
using StatProfilerHTML #profiling the code to see which parts take the most time to run

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
    N_sites = 5
    
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
    initialState = groundstate(N, noPhonons, include3rdLevel)
    initialStateDescription = "gs"
    
    # Whether to have driving on the g-e transition or not
    ΩDriveOn = true
    
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
    cDriveType = "planeWave" # "constant", "planeWave"
    cDriveDescription = "plW" # "cst", "plW"
    
    # Detuning of the control drive with respect to the e-s transition
    Δc = 0.0
    
    # Rabi frequency of the control drive with respect to the e-s transition
    Ωc = 0.5
    
    # Additional arguments for the control drive ("planeWave" requires a momentum vector)
    cDriveArgs = (kc = ωa*[-0.9, 0, 0], )
    
    
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
                  include3rdLevel, cDriveType, cDriveDescription, Δc, Ωc, cDriveArgs) 
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
    
    # Define the fiber
    fiber = Fiber(ρf0_ul, n0, ωa)
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-1, 1, 300)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = ηα0 .* [1, 1, 0]
    ηα = [0, 0, 0]
    
    # Whether phonons are excluded or not from the calculations (a finite ηα but noPhonons = true will result in including ground state motion into tildeG)
    noPhonons = all(ηα .== 0)
    # noPhonons = true
    
    # Whether to include a third (metastable) level to facilitate EIT
    include3rdLevel = true
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set number of atomic sites 
    N_sites = 200
    
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
    tspan = (0, 1000)
    dtmax = 0.01
    
    # Prepare initial state for time evolution, as well as description for postfix
    # initialState = groundstate(N, noPhonons, include3rdLevel)
    # initialStateDescription = "gs"
    GaussWidth = sqrt(N)*a0_ul
    initialState = Gaussian_sState(N, array, fiber, GaussWidth, noPhonons, include3rdLevel)
    # initialState = GaussianState(N, array, 0, N/2*a0_ul, GaussWidth, "e", noPhonons, include3rdLevel)
    initialStateDescription = "Ga"
    
    # Whether to have driving on the g-e transition or not
    ΩDriveOn = false
    
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
    cDriveType = "planeWave" # "constant", "planeWave"
    cDriveDescription = "plW" # "cst", "plW"
    # cDriveType = "constant"
    # cDriveDescription = "cst"
    
    # Detuning of the control drive with respect to the e-s transition
    Δc = 0
    
    # Rabi frequency of the control drive with respect to the e-s transition
    Ωc = 0.1
    
    # Additional arguments for the control drive ("planeWave" requires a momentum vector)
    cDriveArgs = (kc = ωa*[0.0, 0.9, 0.0], )
    
    
    
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
                  include3rdLevel, cDriveType, cDriveDescription, Δc, Ωc, cDriveArgs)
end


function define_SP_ChangExponential()
    # Fiber and array parameters
    n  = 2
    ρf = 1.2/ωa
    ρa = 1.5*ρf
    a  = 0.25
    να = [1, 2, 3]
    
    # Define the fiber
    fiber = Fiber(ρf, n, ωa)
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-10, 10, 3000)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = [0.1377, 0.1826, 0.1219] #Cs
    # ηα = [0.2093, 0.2775, 0.1853] #Sr
    ηα = [0.0, 0.0, 0.0]
    
    # Whether phonons are excluded or not from the calculations (a finite ηα but noPhonons = true will result in including ground state motion into tildeG)
    noPhonons = all(ηα .== 0)
    # noPhonons = false
    
    # Whether to include a third (metastable) level to facilitate EIT
    include3rdLevel = true
    
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
    array, arrayDescription, N = get_array(arrayType, N_sites, ρa, a, ff, pos_unc, n_inst)
    
    # Which time evolver to use ("OrdinaryDiffEq", "simple")
    whichTimeEvolver = "simple"
    
    # Time span and maximum time step allowed in time evolution
    tspan = (0, 1000)
    dtmax = 0.01
    
    # Prepare initial state for time evolution, as well as description for postfix
    # initialState = groundstate(N, noPhonons, include3rdLevel)
    # initialStateDescription = "gs"
    GaussWidth = sqrt(N)*a
    initialState = Gaussian_sState(N, array, fiber, GaussWidth, noPhonons, include3rdLevel)
    # initialState = GaussianState(N, array, 0, N/2*a, GaussWidth, "e", noPhonons, include3rdLevel)
    initialStateDescription = "Ga"
    
    # Whether to have driving on the g-e transition or not
    ΩDriveOn = false
    
    # Atomic dipole moment
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul, array)
    # d = "chiral"
    # dDescription = "chiral"
    # d = rightCircularDipoleMoment(array)
    # dDescription = "rgtCrc"
    d = [[1, 0, 0] for site in array]
    dDescription = "xPol"
    
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
    save_steadyState  = false # n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    save_timeEvol     = false # n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a
    margin = 0.1*arrayL
    z_range = range(-margin, arrayL + margin, 100)
    x_range = range(-ρf - margin, ρf + ρa + margin, 100)
    y_fix   = ρa
    
    # Get the interpolation function for the imaginary, transverse part of the radiation Green's function, if needed
    if interpolate_Im_Grm_trans interpolation_Im_Grm_trans = interpolation1D_Im_Grm_trans(fiber, Int(ceil(arrayL/0.1)) + 1, ρa, 0.1, ηα) else interpolation_Im_Grm_trans = nothing end
    
    # Type of control drive the third level transition
    cDriveType = "hyperbolic" # "constant", "planeWave", "hyperbolic"
    cDriveDescription = "hyp" # "cst", "plW", "hyp"
    
    # Detuning of the control drive with respect to the e-s transition
    Δc = 0
    
    # Rabi frequency of the control drive with respect to the e-s transition
    Ωc = 0.1
    
    # Additional arguments for the control drive ("planeWave" requires a momentum vector)
    cDriveArgs = (kc = ωa*[-1, 0, 0], N_sites=N_sites, a=a)
    
    
    return SysPar(ρf, n, ωa,
                  fiber,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  whichTimeEvolver, tspan, dtmax,
                  initialState, initialStateDescription,
                  ΩDriveOn,
                  arrayType, N_sites, ρa, a, ff, pos_unc, n_inst, array, arrayDescription, N,
                  να, ηα, noPhonons,
                  d, dDescription, incField_wlf, tildeG_flags, 
                  interpolate_Im_Grm_trans, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                  save_steadyState, save_timeEvol,
                  interpolation_Im_Grm_trans,
                  z_range, x_range, y_fix,
                  include3rdLevel, cDriveType, cDriveDescription, Δc, Ωc, cDriveArgs)
end


function define_SP_artificial()
    # Fiber and array parameters
    n  = 1.45
    ρf = 0.167
    ρa = 0.3991
    a  = 0.2903
    να = [1, 2, 3]
    
    # Define the fiber
    fiber = Fiber(ρf, n, ωa)
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs = (-10, 10, 3000)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = [0.1377, 0.1826, 0.1219] #Cs
    # ηα = [0.2093, 0.2775, 0.1853] #Sr
    # ηα = [0.1784, 0.2366, 0.1580]
    ηα = [0.0, 0.0, 0.0]
    
    # Whether phonons are excluded or not from the calculations (a finite ηα but noPhonons = true will result in including ground state motion into tildeG)
    noPhonons = all(ηα .== 0)
    # noPhonons = false
    
    # Whether to include a third (metastable) level to facilitate EIT
    include3rdLevel = false
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set number of atomic sites 
    N_sites = 20
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 1.0
    pos_unc = 0.0
    # pos_unc = ηα0/ωa
    n_inst = 1
    
    # Generate the array, its description, and the number of atoms
    array, arrayDescription, N = get_array(arrayType, N_sites, ρa, a, ff, pos_unc, n_inst)
    
    # Which time evolver to use ("OrdinaryDiffEq", "simple")
    whichTimeEvolver = "simple"
    
    # Time span and maximum time step allowed in time evolution
    tspan = (0, 20)
    dtmax = 0.01
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N, noPhonons, include3rdLevel)
    initialStateDescription = "gs"
    # GaussWidth = sqrt(N)*a
    # initialState = Gaussian_sState(N, array, fiber, GaussWidth, noPhonons, include3rdLevel)
    # # initialState = GaussianState(N, array, 0, N/2*a, GaussWidth, "e", noPhonons, include3rdLevel)
    # initialStateDescription = "Ga"
    
    # Whether to have driving on the g-e transition or not
    ΩDriveOn = true
    
    # Atomic dipole moment
    # d = chiralDipoleMoment(Fiber(ρf0_ul, n0, ωa), ρa0_ul, array)
    # d = "chiral"
    # dDescription = "chiral"
    # d = rightCircularDipoleMoment(array)
    # dDescription = "rgtCrc"
    d = [[1, 0, 0] for site in array]
    dDescription = "xPol"
    
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
    save_steadyState  = false # n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    save_timeEvol     = false # n_inst == 1 && ff == 1 && pos_unc == 0 && arrayType != "randomZ"
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a
    margin = 0.1*arrayL
    z_range = range(-margin, arrayL + margin, 100)
    x_range = range(-ρf - margin, ρf + ρa + margin, 100)
    y_fix   = ρa
    
    # Get the interpolation function for the imaginary, transverse part of the radiation Green's function, if needed
    if interpolate_Im_Grm_trans interpolation_Im_Grm_trans = interpolation1D_Im_Grm_trans(fiber, Int(ceil(arrayL/0.1)) + 1, ρa, 0.1, ηα) else interpolation_Im_Grm_trans = nothing end
    
    # Type of control drive the third level transition
    cDriveType = "constant" # "constant", "planeWave"
    cDriveDescription = "cst" # "cst", "plW"
    
    # Detuning of the control drive with respect to the e-s transition
    Δc = 0
    
    # Rabi frequency of the control drive with respect to the e-s transition
    Ωc = 0.1
    
    # Additional arguments for the control drive ("planeWave" requires a momentum vector)
    cDriveArgs = (kc = ωa*[-1, 0, 0], )
    
    
    return SysPar(ρf, n, ωa,
                  fiber,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  whichTimeEvolver, tspan, dtmax,
                  initialState, initialStateDescription,
                  ΩDriveOn,
                  arrayType, N_sites, ρa, a, ff, pos_unc, n_inst, array, arrayDescription, N,
                  να, ηα, noPhonons,
                  d, dDescription, incField_wlf, tildeG_flags, 
                  interpolate_Im_Grm_trans, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                  save_steadyState, save_timeEvol,
                  interpolation_Im_Grm_trans,
                  z_range, x_range, y_fix,
                  include3rdLevel, cDriveType, cDriveDescription, Δc, Ωc, cDriveArgs)
end

    
function main()
    # Define system parameters
    # SP = define_SP_BerlinCs()
    # SP = define_SP_BerlinSr()
    SP = define_SP_ChangExponential()
    # SP = define_SP_artificial()
    # show(SP)
    
    
    
    # plot_propConst_inOutMom(ωρfn_ranges)
    # plot_coupling_strengths(SP)
    # plot_arrayIn3D(SP)
    # plot_interPairEnergiesWeights(SP)
    # plot_σBαTrajectories_σBαSS(SP)
    # plot_transmission_vs_Δ(SP)
    # plot_imperfectArray_transmission_vs_Δ(SP)
    # plot_compareImperfectArray_transmission_vs_Δ(SP)
    # plot_effectiveBetaFactor(SP)
    # plot_effectiveBetaFactor_vs_eta(SP)
    # plot_effectiveBetaFactor_perfectArray(SP)
    # plot_steadyState_radiation_Efield(SP)
    # plot_radiation_Efield(SP)
    # plot_GnmEigenModes(SP)
    # plot_emissionPatternOfGnmeigenModes(SP)
    # plot_GnmEigenEnergies(SP)
    # plot_GnmFourierTransformed(SP)
    # plot_compareGnmEigenEnergies(SP)
    # plot_lossWithGnmEigenEnergies(SP)
    plot_memoryEfficiency(SP)
    
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


function plot_interPairEnergiesWeights(SP)
    # We calculate the eigenenergies for a single pair of atoms with a certain distance between them
    # Likewise, we calculate their effective Rabi frequency (their weigth in the expression for the transmission),
    # and the amplitude of their contribution to the transmission
    # We also find the probability that any pair of atoms has a certain distance between them, 
    # given a number of atoms randomly distributed in a certain interval
    # With these we can analyse what strong contributions there would be to the transmission, 
    # when we have a randomZ array
    
    # Set number of atoms and length of interval (for calculation of probabilities)
    N_atoms = 100
    ff = 0.05
    L = N_atoms*SP.a/ff
    # L = 30
    
    # Set whether to include phonon modes
    noPhonons = true
    
    # Set the inter-pair distances we use
    n = 100
    interPairDistance = range(0.06, 0.5, n)
    # dx = L/n
    # interPairDistance = range(dx, L, n)
    
    # Find arrays for these distances
    arrays = [[[SP.ρa, 0, 0], [SP.ρa, 0, r]]  for r in interPairDistance]
    if noPhonons N_modes = 2 else N_modes = 2 + 3*2^2 end
    
    
    # Numerically calculate probability distributions
    # bins = zeros(Int, n + 1)
    
    # Probability for distances between any pair of atoms
    # for _ in 1:1000
    #     positions = rand(N_atoms)*L
    #     pairwise_distances = [abs(p1 - p2) for p1 in positions, p2 in positions]
    #     pairwise_distances_flat = pairwise_distances[pairwise_distances .!= 0]
        
    #     bins[1] += sum(pairwise_distances_flat .< interPairDistance[1])
    #     for i in 1:n-1
    #         bins[i + 1] += sum(interPairDistance[i] .< pairwise_distances_flat .< interPairDistance[i + 1])
    #     end
    #     bins[end] += sum(interPairDistance[end] .< pairwise_distances_flat)
    # end
    # probabilities = bins/sum(bins*dx)
    # theoretical_line = @. 2*(L - interPairDistance)/L^2
    
    # Probability for distances between consecutive atoms
    # for _ in 1:10000
    #     positions = rand(N_atoms)*L
    #     pairwise_distances = diff(sort(positions))
        
    #     bins[1] += sum(pairwise_distances .< interPairDistance[1])
    #     for i in 1:n-1
    #         bins[i + 1] += sum(interPairDistance[i] .< pairwise_distances .< interPairDistance[i + 1])
    #     end
    #     bins[end] += sum(interPairDistance[end] .< pairwise_distances)
    # end
    # probabilities = bins/sum(bins*dx)
    # theoretical_line = @. exp(-N_atoms/L*interPairDistance)/(L/N_atoms*(1 - exp(-N_atoms)))
    
    # # Start figure for plotting numerically calculated probability distribution
    # fig = Figure(size=(600, 400))
    # ax1 = Axis(fig[1, 1])
    # lines!(ax1, vcat([0], interPairDistance), probabilities)
    # lines!(ax1, interPairDistance, theoretical_line)
    # display(GLMakie.Screen(), fig)
    
    
    # Find eigenenergies for a pair of atoms, as well as the weigth of the corresponding contributions to the transmission, the amplitude of those contributions, and the amplitude weigthed with the probability of that specific inter-atom distance
    collΔs, collΓs, weights, amplitudes, amplitudes_prob = zeros(n, N_modes), zeros(n, N_modes), zeros(n, N_modes), zeros(n, N_modes), zeros(n, N_modes)
    for (i, array) in enumerate(arrays)
        if noPhonons
            Δvari, tildeΩ, tildeG = get_parameterMatrices(noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, array, SP.ΩDriveOn, SP.tildeG_flags, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
            drive = tildeΩ
            fullCoupling = Δvari + tildeG
        else
            Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2 = get_parameterMatrices(noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, array, SP.ΩDriveOn, SP.tildeG_flags, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
            drive = get_fullDriveVector(tildeΩ, tildeΩα)
            fullCoupling = get_fullCouplingMatrix(Δvari, tildeG, tildeFα, tildeGα1, tildeGα2)
        end
        
        eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv = spectrum_basisMatrices(fullCoupling)
        collΔ, collΓ = collEnergies(eigenEnergies)
        
        weight, resonances = transmission_eigenmodes_weights_resonances(SP.Δ_range, drive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, SP.fiber.propagation_constant_derivative)
        
        # Probability for a certain distance between any pair of atoms 
        # probability = @. 2*(L - interPairDistance[i])/L^2
        # Probability for a certain distance between consecutive atoms
        probability = @. exp(-N_atoms/L*interPairDistance[i])/(L/N_atoms*(1 - exp(-N_atoms)))
        
        collΔs[i, :]          = collΔ
        collΓs[i, :]          = collΓ
        weights[i, :]         = abs.(weight)
        amplitudes[i, :]      = abs.(2*weight./collΓ)
        amplitudes_prob[i, :] = amplitudes[i, :]*probability
        
        # println(i)
    end
    
    
    # Start figure for plotting the energies and weigths
    fig = Figure(size=(600, 800))
    ax1 = Axis(fig[1, 1], ylabel=L"Eigenenergies$$")
    ax2 = Axis(fig[2, 1], ylabel=L"Decay rates$$")
    ax3 = Axis(fig[3, 1], ylabel=L"Effective Rabi frequencies$$", xlabel=L"Inter-atom distance, $ r/\lambda_{a} $")
    for i in 1:N_modes
        scatter!(ax1, interPairDistance, collΔs[:, i], color=:blue)
        scatter!(ax2, interPairDistance, collΓs[:, i], color=:red)
        scatter!(ax3, interPairDistance, weights[:, i], color=:purple)
    end
    
    # collΔs1  = [collΔs[i, findfirst(collΓs[i, :] .> 1.05)] for i in 1:n]
    # collΔs2  = [collΔs[i, findfirst(collΓs[i, :] .< 1.05)] for i in 1:n]
    # collΓs1  = [collΓs[i, findfirst(collΓs[i, :] .> 1.05)] for i in 1:n]
    # collΓs2  = [collΓs[i, findfirst(collΓs[i, :] .< 1.05)] for i in 1:n]
    # weights1 = [weights[i, findfirst(collΓs[i, :] .> 1.05)] for i in 1:n]
    # weights2 = [weights[i, findfirst(collΓs[i, :] .< 1.05)] for i in 1:n]
    # for (collΔs_i, collΓs_i, weights_i, clr) in zip((collΔs1, collΔs2), (collΓs1, collΓs2), (weights1, weights2), (:blue, :red))
    #     lines!(ax1, interPairDistance, collΔs_i, color=clr)
    #     lines!(ax2, interPairDistance, collΓs_i, color=clr)
    #     lines!(ax3, interPairDistance, weights_i, color=clr)
    # end
    
    display(GLMakie.Screen(), fig)
    # save("C:\\Users\\Simon\\Downloads\\atomPairEnergies.png", fig)
    
    # Start figure for the transmission contribution amplitudes weigthed by distance probability
    fig = Figure(size=(600, 400))
    ax1 = Axis(fig[1, 1], limits=(-2, 2, nothing, nothing), xlabel=L"Eigenenergies$$", ylabel=L"transmission contribution amplitude$$")
    for i in 1:N_modes
        scatter!(ax1, collΔs[:, i], amplitudes_prob[:, i], color=:blue)
    end
    
    # amplitudes1 = [amplitudes[i, findfirst(collΓs[i, :] .> 1.05)] for i in 1:n]
    # amplitudes2 = [amplitudes[i, findfirst(collΓs[i, :] .< 1.05)] for i in 1:n]
    # for (collΔs_i, amplitudes_i, clr) in zip((collΔs1, collΔs2), (amplitudes1, amplitudes2), (:blue, :red))
    #     lines!(ax1, collΔs_i, amplitudes_i, color=clr)
    # end
    
    display(GLMakie.Screen(), fig)
    # save("C:\\Users\\Simon\\Downloads\\atomPairTransmissionContributions.png", fig)
end


function plot_σBαTrajectories_σBαSS(SP)
    # Δ = 1.0
    Δ = SP.Δc + eps(Float64)
    σBα_SS = calc_steadyState(SP, Δ)
    xTrajectories = calc_timeEvolution(SP, Δ)
    # xTrajectories = calc_timeEvolution_eigenmodes(SP, Δ)
    
    if SP.noPhonons
        if SP.include3rdLevel 
            times, σgeTrajectories, σgsTrajectories = prep_times_σgeσgsTrajectories(xTrajectories, SP.N)
            fig_σTrajectories_σSS(times, σgeTrajectories, σBα_SS[1])
            fig_σTrajectories_σSS(times, σgsTrajectories, σBα_SS[2])
        else
            times, σTrajectories = prep_times_σTrajectories(xTrajectories, SP.N)
            fig_σTrajectories_σSS(times, σTrajectories, σBα_SS)
        end
    else
        if SP.include3rdLevel
            times, σgeTrajectories, σgsTrajectories, BαgeTrajectories, BαgsTrajectories = prep_times_σgeσgsBαgeBαgsTrajectories(xTrajectories, SP.N)
            fig_σBαTrajectories_σBαSS(times, σgeTrajectories, BαgeTrajectories, σBα_SS[[1, 3]]...)
            fig_σBαTrajectories_σBαSS(times, σgsTrajectories, BαgsTrajectories, σBα_SS[[2, 4]]...)
        else
            times, σTrajectories, BαTrajectories = prep_times_σBαTrajectories(xTrajectories, SP.N)
            fig_σBαTrajectories_σBαSS(times, σTrajectories, BαTrajectories, σBα_SS...)
        end
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
    # fig_transmission_vs_Δ_phaseDetails_polar(SP.Δ_range, T, tPhase, t, unwrappedPhase, phaseSlope, titl)
    
    r = calc_reflection.(Ref(SP), σBα_scan)
    R, rPhase = prep_squaredNorm_phase(r)
    fig_transmissionAndReflection_vs_Δ(SP.Δ_range, T, tPhase, R, rPhase, titl)
    
    # titl = L"$ N = %$(SP.N) $"
    # fig = fig_presentation_transmission_vs_Δ(SP.Δ_range, T, phase, titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\transmission_N500_linear_negative_zoom1.png", fig, px_per_unit=2)
    
    # n = refractiveIndex.(t, SP.arrayType, SP.N_sites, SP.a)
    # n_real, n_imag = real.(n), imag.(n)
    # groupVel = groupVelocity(n, SP.Δ_range)
    # fig_RefrIndexGroupVelocity_vs_Δ(SP.Δ_range, n_real, n_imag, groupVel, titl)
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
    # ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7]
    ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    # ηαFactor_list = [0.0, 0.1, 0.4, 0.7, 1.0]
    # ηαFactor_list = [1.0]
    # ηαFactor_list = [1.0, [1, 1, 1], [1, 1, 0]]
    # ηαFactor_list = [1.0, [1, 1, 0]]
    # ηα_list = [SP.ηα .* [1, 1, 0], [0.0, 0.0, 0.0]]
    ηα_list = [SP.ηα, [0.0, 0.0, 0.0]]
    # pos_unc_list = [0.0, SP.ηα/ωa .* [1, 1, 0]]
    pos_unc_list = [0.0, SP.ηα/ωa]
    # Nsites_list = [2000, 1000, 667, 500, 400, 334, 250, 200, 167, 143] # N = 100
    Nsites_list = [2000, 1000, 667, 500, 400, 334, 250, 200, 167, 143, 125, 112, 100] # N = 100
    # Nsites_list = [6000, 3000, 2000, 1500, 1200, 1000, 750, 600, 500, 429] # N = 300
    # Nsites_list = [ 200   100    67    50    40    34    25    20   17   15
    #                 800   400   267   200   160   134   100    80   67   58
    #                1400   700   467   350   280   234   175   140  117  100
    #                2000  1000   667   500   400   334   250   200  167  143
    #                4000  2000  1334  1000   800   667   500   400  334  286
    #                6000  3000  2000  1500  1200  1000   750   600  500  429
    #                8000  4000  2667  2000  1600  1334  1000   800  667  572
    #               10000  5000  3334  2500  2000  1667  1250  1000  834  715]
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
    # Ns = [10, 40, 70, 100, 200, 300, 400, 500]
    # Ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # arrayType_list = ["1Dchain", "randomZ", "randomZ"]
    arrayType_list = ["1Dchain", "randomZ"]
    # for ηαFactor in ηαFactor_list
    # for ηαFactor in ηαFactor_list, (i, ff) in enumerate(ff_list)
    # for (i, ff) in enumerate(ff_list)
    # for (arrayType, ηαFactor) in zip(arrayType_list, ηαFactor_list)
    for (ff, N_sites) in zip(ff_list, Nsites_list)
    # for (arrayType, ηαFactor) in zip(arrayType_list, ηαFactor_list), (i, ff) in enumerate(ff_list)
        T_meanss, T_stdss, phase_meanss, phase_stdss = [], [], [], []
        T_indepDecayss, phase_indepDecayss = [], []
        refrIndex_realss, refrIndex_imagss = [], []
        refrIndex_real_indepDecayss, refrIndex_imag_indepDecayss = [], []
        labels = []
        # for ηαFactor in ηαFactor_list
        # for ff in ff_list
        # for (ff, N_sites) in zip(ff_list, Nsites_list)
        # for N_sites in Nsites_list[:, i]
        # for arrayType in arrayType_list, ff in ff_list
        # for (arrayType, ηαFactor) in zip(arrayType_list, ηαFactor_list), ff in ff_list
        # for (arrayType, ηαFactor) in zip(arrayType_list, ηαFactor_list), N_sites in Nsites_list[:, i]
        # for (arrayType, ηαFactor) in zip(arrayType_list, ηαFactor_list), (ff, N_sites) in zip(ff_list, Nsites_list)
        for (ηα, pos_unc) in zip(ηα_list, pos_unc_list)
        # for (ηα, pos_unc) in zip(ηα_list, pos_unc_list), (ff, N_sites) in zip(ff_list, Nsites_list)
        # for arrayType in arrayType_list, (ff, N_sites) in zip(ff_list, Nsites_list)
        # for ff in ff_list, ηαFactor in ηαFactor_list
            # ηα = SP.ηα .* ηαFactor
            # pos_unc = SP.pos_unc #.* ηαFactor
            # arrayDescription = arrayDescript(SP.arrayType, SP.N_sites, SP.ρa, SP.a, ff, pos_unc)
            # arrayDescription = arrayDescript(arrayType, SP.N_sites, SP.ρa, SP.a, ff, pos_unc)
            arrayDescription = arrayDescript(SP.arrayType, N_sites, SP.ρa, SP.a, ff, pos_unc)
            # arrayDescription = arrayDescript(arrayType, N_sites, SP.ρa, SP.a, ff, pos_unc)
            
            postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, arrayDescription, SP.fiber.postfix)
            filename = "T_phase_" * postfix
            folder = "imperfectArray_T_phase/"
        
            if isfile(saveDir * folder * filename * ".txt") 
                # Gather labels
                # push!(labels, L"$ ff = %$(ff) $, $ ηα = %$(ηαFactor) \cdot ηα0 $")
                # push!(labels, L"$ ff = %$(ff) $")
                # push!(labels, L"$ N_{atoms} = %$(Int(floor(N_sites*ff))) $")
                # push!(labels, L"$ N_{sites} = %$(N_sites) $, $ ff = %$(ff) $")
                # push!(labels, L"$ N_{sites} = %$(N_sites) $, $ ff = %$(ff) $, %$(arrayType), %$(ηαFactor)")
                # push!(labels, L"$ N_{sites} = %$(N_sites) $, $ ff = %$(ff) $, $ ηα = %$(round.(ηα, digits=2)) $, pu $ = %$(round.(pos_unc, digits=2)) $")
                push!(labels, L"$ ff = %$(ff) $, pu $ = %$(round.(pos_unc, digits=2)) $")
                # push!(labels, L"$ η_{α} = %$(ηαFactor) \cdot η_{α}^{(0)} $")                
                
                # Gather transmissions
                push!.([T_meanss, T_stdss, phase_meanss, phase_stdss], prep_T_argt_statistics(eachrow(load_as_txt(saveDir * folder, filename))...))
                
                # # Gather independent decay transmissions
                # γ_gm, γ_rm = get_γs(SP)
                # N = Int(floor(SP.N_sites*ff))
                # # N = Int(floor(N_sites*ff))
                # t_indepDecay = transmission_indepDecay.(SP.Δ_range, γ_gm, γ_rm, N)
                # push!.([T_indepDecayss, phase_indepDecayss], prep_squaredNorm_phase(t_indepDecay))
                
                # # Gather refractive indices
                # refrIndex_complex = refractiveIndex.(sqrt.(T_meanss[end]).*exp.(1im*phase_meanss[end]), SP.arrayType, SP.N_sites, SP.a)
                # refrIndex_indepDecay_complex = refractiveIndex.(sqrt.(T_indepDecayss[end]).*exp.(1im*phase_indepDecayss[end]), SP.arrayType, SP.N_sites, SP.a)
                # # refrIndex_complex = refractiveIndex.(sqrt.(T_meanss[end]).*exp.(1im*phase_meanss[end]), SP.arrayType, N_sites, SP.a)
                # # refrIndex_indepDecay_complex = refractiveIndex.(sqrt.(T_indepDecayss[end]).*exp.(1im*phase_indepDecayss[end]), SP.arrayType, N_sites, SP.a)
                # push!.([refrIndex_realss, refrIndex_imagss], [real.(refrIndex_complex), imag.(refrIndex_complex)])
                # push!.([refrIndex_real_indepDecayss, refrIndex_imag_indepDecayss], [real.(refrIndex_indepDecay_complex), imag.(refrIndex_indepDecay_complex)])
            else
                throw(ArgumentError("The following file can not be found: " * filename))
            end
        end
        
        
        # titl = prep_imperfectArray_transmission_title(SP)
        arrayTypeDescription = SP.arrayType == "1Dchain" ? L"ordered 1D chain$$" : L"random $ z $-coordinate"
        # titl = L"Fixed $ N_{sites} = %$(SP.N_sites) $, array type: %$(arrayTypeDescription)"
        # titl = L"Fixed $ N_{atoms} = %$(SP.N_sites) $, array type: %$(arrayTypeDescription)"
        # titl = L"Fixed ff$ = %$(ff) $, array type: %$(arrayTypeDescription)"
        titl = L"$ N_{atoms} = 100 $, ff$ = %$(ff) $, array type: %$(arrayTypeDescription)"
        labels = ["w. phonons", "clas. disorder"]
        fig = fig_compareImperfectArray_transmission_vs_Δ(SP.Δ_range, T_meanss, T_stdss, phase_meanss, phase_stdss, labels, titl)
        # save("C:\\Users\\Simon\\Downloads\\phononVSclassDisorder_ff$(ff)_Sr_rcirc.png", fig)
        # fig_compareImperfectArray_refrIndex_vs_Δ(SP.Δ_range, refrIndex_realss, refrIndex_imagss, labels, titl)
        # fig_compareImperfectArray_transmission_vs_Δ(SP.Δ_range, T_indepDecayss, fill(zeros(size(T_indepDecayss[1])), size(T_indepDecayss)), phase_indepDecayss, fill(zeros(size(phase_indepDecayss[1])), size(phase_indepDecayss)), labels, titl)
        # fig_compareImperfectArray_refrIndex_vs_Δ(SP.Δ_range, refrIndex_real_indepDecayss, refrIndex_imag_indepDecayss, labels, titl)
        
        # Δ_index = 100
        # T_means, T_stds, phase_means, phase_stds, T_indepDecays, phase_indepDecays, refrIndex_reals, refrIndex_imags, refrIndex_real_indepDecays, refrIndex_imag_indepDecays = prep_compareImperfectArray_transmission_vs_X(Δ_index, T_meanss, T_stdss, phase_meanss, phase_stdss, T_indepDecayss, phase_indepDecayss, refrIndex_realss, refrIndex_imagss, refrIndex_real_indepDecayss, refrIndex_imag_indepDecayss)
        # titl = titl * "\nηα = $(ηαFactor) * ηα0" * "\nΔ = $(SP.Δ_range[Δ_index])"
        # titl = titl * "\nηα = $(ηαFactor) * ηα0, ff = $(ff)" * "\nΔ = $(SP.Δ_range[Δ_index])"
        # titl = titl * "\nΔ = $(SP.Δ_range[Δ_index])"
        # fig_compareImperfectArray_transmission_vs_X(ff_list, L"Filling fraction $$", T_means, T_stds, phase_means, phase_stds, T_indepDecays, phase_indepDecays, titl)
        # fig_compareImperfectArray_transmission_vs_X(Ns, L"$ N $", T_means, T_stds, phase_means, phase_stds, T_indepDecays, phase_indepDecays, titl)
        
        # titl = L"$ N = %$(SP.N) $, random $ z_{n} $"
        # titl = L"$ N_{sites} = 1000 $, ordered $ z_{n} $"
        # titl = L"$ N_{sites} = 1000 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # titl = L"$ N = 100 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # labels = ff_list
        # labels = ηαFactor_list
        # fig = fig_presentation_compareImperfectArray_transmission_vs_Δ(SP.Δ_range, T_meanss, T_stdss, phase_meanss, phase_stdss, labels, titl)
        # save("C:\\Users\\Simon\\Downloads\\compareff_N100_random.png", fig, px_per_unit=2)
        
        # phase_means = vcat(unwrapPhase(phase_means[1:length(ff_list)], "forcePositiveSlope"), unwrapPhase(phase_means[length(ff_list)+1:end], "forcePositiveSlope"))
        # phase_indepDecays = unwrapPhase(phase_indepDecays, "forcePositiveSlope")
        
        # titl = L"$ N_{atoms} = 100 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # # titl = L"$ N_{sites} = 1000 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # # titl = L"ff$ = 0.4 $, $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
        # # fig = fig_presentation_compareImperfectArray_transmission_vs_ffOrηα(Ns, L"$ N_{atoms} $", T_means, T_stds, phase_means, phase_stds, T_indepDecays[1:length(T_indepDecays)÷2], phase_indepDecays[1:length(phase_indepDecays)÷2], titl)
        # fig = fig_presentation_compareImperfectArray_transmission_vs_ffOrηα(ff_list, L"Filling fraction $$", T_means, T_stds, phase_means, phase_stds, T_indepDecays[1:length(T_indepDecays)÷2], phase_indepDecays[1:length(phase_indepDecays)÷2], titl)
        # fig = fig_presentation_compareImperfectArray_refrIndex_vs_ffOrηα(ff_list, L"Filling fraction $$", refrIndex_reals, refrIndex_imags, refrIndex_real_indepDecays[1:length(ff_list)], refrIndex_imag_indepDecays[1:length(ff_list)], titl)
        # # save("C:\\Users\\Simon\\Downloads\\compareff_refrIndex_N100_D$(round(SP.Δ_range[Δ_index], digits=2))_etaz0.png", fig, px_per_unit=2)
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
    ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0][1:7]
    Nsites_list = [ 200   100   67   50   40   34   25   20   17   15   13   12   10;
                    400   200  134  100   80   67   50   40   34   29   25   23   20;
                    600   300  200  150  120  100   75   60   50   43   38   34   30;
                    800   400  267  200  160  134  100   80   67   58   50   45   40;
                   1000   500  334  250  200  167  125  100   84   72   63   56   50;
                   1200   600  400  300  240  200  150  120  100   86   75   67   60;
                   1400   700  467  350  280  234  175  140  117  100   88   78   70;
                   1600   800  534  400  320  267  200  160  134  115  100   89   80;
                   1800   900  600  450  360  300  225  180  150  129  113  100   90;
                   2000  1000  667  500  400  334  250  200  167  143  125  112  100][1:7, :]
    Ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100][1:7]
    
    arrayType_list = ["1Dchain", "randomZ"]
    ηαFactor_list = [[1, 1, 1], [1, 1, 0]]
    
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
        ηα = SP.ηα .* ηαFactor_list[i]
        postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, arrayDescription, SP.fiber.postfix)
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
    n_effs  = fill(0.0im, length(arrayType_list), length(ff_list), length(SP.Δ_range))
    for (i, arrayType) in enumerate(arrayType_list), (j, ff) in enumerate(ff_list), (l, Δ) in enumerate(SP.Δ_range)
        # println("i = ", i, ", j = ", j, ", l = ", l)
        model(x, (β_eff, Δ_eff)) = (1 - 2*β_eff/(1 - 2im*Δ_eff))^x
        t           = t_real_means[i, j, :, l] + 1im*t_imag_means[i, j, :, l]
        tDeviations = t_real_stds[i, j, :, l]  + 1im*t_imag_stds[i, j, :, l]
        # pmin = fitComplexData(Ns, t, model, [β_indepDecay*(1 + 5*ff^2), Δ*(1 + 5*ff^2)]; ydataDeviations=tDeviations)
        pmin = fitComplexData(Ns, t, model, [β_indepDecay*(1 + 1*ff^2), Δ*(1 + 1*ff^2)])
        
        T_fits[i, j, :, l], phase_fits[i, j, :, l] = prep_squaredNorm_phase(model.(Ns, Ref(pmin)))
        β_effs[i, j, l], Δ_effs[i, j, l] = pmin
        n_effs[i, j, l] = 1 + ff/(1im*ωa*SP.a)*log(1 - 2*β_effs[i, j, l]/(1 - 2im*Δ_effs[i, j, l]))
    end
    
    # Unwrap phase for clarity in plots
    phase_means       = mapslices(x -> unwrapPhase(x, "forcePositiveSlope"), phase_means      , dims=3)
    phase_indepDecays = mapslices(x -> unwrapPhase(x, "forcePositiveSlope"), phase_indepDecays, dims=1)
    phase_fits        = mapslices(x -> unwrapPhase(x, "forcePositiveSlope"), phase_fits       , dims=3)
    
    # # Plot example of fits
    # Δ_index = 50
    # for (i, arrayType) in enumerate(arrayType_list)
    #     labels = [L" ff$ = %$(ff) $, $ β = %$(round(β_effs[i, j, Δ_index], digits=3)) $" for (j, ff) in enumerate(ff_list)]
    #     titl = prep_imperfectArray_transmission_title(SP) * "\narrayType: $arrayType" * "\nΔ = $(SP.Δ_range[Δ_index])"
    #     fig_compareImperfectArray_transmission_vs_N(Ns, T_means[i, :, :, Δ_index], T_stds[i, :, :, Δ_index], phase_means[i, :, :, Δ_index], phase_stds[i, :, :, Δ_index], T_indepDecays[:, Δ_index], phase_indepDecays[:, Δ_index], T_fits[i, :, :, Δ_index], phase_fits[i, :, :, Δ_index], β_indepDecay, labels, titl)
    # end
    
    # Plot effective β-factors and detunings as well as corresponding refractive indices
    for (i, arrayType) in enumerate(arrayType_list)
        labels = [L" ff$ = %$(ff) $" for (j, ff) in enumerate(ff_list)]
        # titl = prep_imperfectArray_transmission_title(SP) * "\narrayType: $arrayType"
        arrayTypeDescription = arrayType == "1Dchain" ? L"ordered 1D chain$$" : L"random $ z $-coordinate"
        titl = L"Array type: %$(arrayTypeDescription)"
        fig = fig_effectiveβΔ_vs_Δ(SP.Δ_range, β_effs[i, :, :], Δ_effs[i, :, :], β_indepDecay, labels, titl)
        # save("C:\\Users\\Simon\\Downloads\\effectiveBetaDelta_$(arrayType)_Sr_rcirc.png", fig)
        fig = fig_effectiveRefrIndex_vs_Δ(SP.Δ_range, real(n_effs[i, :, :]), imag(n_effs[i, :, :]), labels, titl)
        # save("C:\\Users\\Simon\\Downloads\\effectiveRefrIndex_$(arrayType)_Sr_rcirc.png", fig)
    end
    
    # # arrayTypeDescription = arrayType == "1Dchain" ? L"ordered 1D chain$$" : L"random $ z $-coordinate"
    # # titl = L"Array type: %$(arrayTypeDescription)"
    # titl = L"$ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
    # fig = fig_presentation_effectiveβΔ_vs_ff(ff_list, β_effs[:, :, Δ_index], Δ_effs[:, :, Δ_index], β_indepDecay, SP.Δ_range[Δ_index], titl)
    # # save("C:\\Users\\Simon\\Downloads\\effectiveBetaVSff_D$(round(SP.Δ_range[Δ_index], digits=2)).png", fig, px_per_unit=2)
    
end


function plot_effectiveBetaFactor_perfectArray(SP)
    Ns = 200:50:750
    
    # Load pre-calculated steady states and calculate transmissions
    ts = zeros(ComplexF64, length(Ns), length(SP.Δ_range))
    Ts, phases = zeros(size(ts)), zeros(size(ts))
    for (i, N) in enumerate(Ns)
        array = get_array(SP.arrayType, N, SP.ρa, SP.a, 1.0, 0.0)
        arrayDescription = arrayDescript(SP.arrayType, N, SP.ρa, SP.a, 1.0, 0.0)
        postfixes = get_postfix_steadyState.(SP.Δ_range, SP.ΔvariDescription, SP.dDescription, Ref(SP.να), Ref(SP.ηα), Ref(SP.incField_wlf), Ref(SP.tildeG_flags), arrayDescription, SP.fiber.postfix, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, Ref(SP.cDriveArgs))
        for (j, postfix) in enumerate(postfixes)
            if SP.noPhonons filename = "sigma_" * postfix else filename = "sigmaBalpha_" * postfix end
            folder = "steadyStates/"
            if isfile(saveDir * folder * filename * ".jld2")
                σBα = load_as_jld2(saveDir * folder, filename)
                
                if SP.noPhonons
                    tildeΩ = get_tildeΩs(SP.fiber, SP.d, SP.incField_wlf, array, SP.ΩDriveOn)
                    ts[i, j] = transmission(σBα, tildeΩ, SP.fiber)
                else
                    tildeΩ, tildeΩα = get_tildeΩs(SP.fiber, SP.d, SP.ηα, SP.incField_wlf, array, SP.ΩDriveOn)
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


function plot_effectiveBetaFactor_vs_eta(SP)
    if SP.n_inst == 1 throw(ArgumentError("plot_effectiveBetaFactor requires n_inst > 1")) end
    
    ff_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]#[1:7]
    Nsites_list = [ 200   100   67   50   40   34   25   20   17   15   13   12   10;
                    400   200  134  100   80   67   50   40   34   29   25   23   20;
                    600   300  200  150  120  100   75   60   50   43   38   34   30;
                    800   400  267  200  160  134  100   80   67   58   50   45   40;
                   1000   500  334  250  200  167  125  100   84   72   63   56   50;
                   1200   600  400  300  240  200  150  120  100   86   75   67   60;
                   1400   700  467  350  280  234  175  140  117  100   88   78   70;
                   1600   800  534  400  320  267  200  160  134  115  100   89   80;
                   1800   900  600  450  360  300  225  180  150  129  113  100   90;
                   2000  1000  667  500  400  334  250  200  167  143  125  112  100]#[1:7, :]
    Ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]#[1:7]
    
    ναFactor_list = [0.444, 0.64, 1.0, 1.777, 4.0]
    ηαFactor_list = [1.5, 1.25, 1.0, 0.75, 0.5]
    
    # Load pre-calculated transmissions
    t_real_means = zeros(length(ναFactor_list), length(ff_list), length(Ns), length(SP.Δ_range))    
    t_real_stds, t_imag_means, t_imag_stds, T_means, T_stds, phase_means, phase_stds = deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means), deepcopy(t_real_means)
    for (i, ναFactor) in enumerate(ναFactor_list), (j, ff) in enumerate(ff_list), (k, N_sites) in enumerate(Nsites_list[:, j])
        arrayDescription = arrayDescript(SP.arrayType, N_sites, SP.ρa, SP.a, ff, SP.pos_unc)
        να = SP.να * ναFactor
        ηα = SP.ηα / sqrt(ναFactor)
        postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, να, ηα, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, arrayDescription, SP.fiber.postfix)
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
    γ_gm, γ_rm = get_γs(SP)
    β_indepDecay = γ_gm/(γ_gm + γ_rm)
    β_effs  = fill(NaN, length(ναFactor_list), length(ff_list), length(SP.Δ_range))
    Δ_effs  = deepcopy(β_effs)
    n_effs  = fill(0.0im, length(ναFactor_list), length(ff_list), length(SP.Δ_range))
    for (i, ναFactor) in enumerate(ναFactor_list), (j, ff) in enumerate(ff_list), (l, Δ) in enumerate(SP.Δ_range)
        model(x, (β_eff, Δ_eff)) = (1 - 2*β_eff/(1 - 2im*Δ_eff))^x
        t           = t_real_means[i, j, :, l] + 1im*t_imag_means[i, j, :, l]
        tDeviations = t_real_stds[i, j, :, l]  + 1im*t_imag_stds[i, j, :, l]
        # pmin = fitComplexData(Ns, t, model, [β_indepDecay*(1 + 5*ff^2), Δ*(1 + 5*ff^2)]; ydataDeviations=tDeviations)
        pmin = fitComplexData(Ns, t, model, [β_indepDecay*(1 + 5*ff^2), Δ*(1 + 5*ff^2)])
        
        β_effs[i, j, l], Δ_effs[i, j, l] = pmin
        n_effs[i, j, l] = 1 + ff/(1im*ωa*SP.a)*log(1 - 2*β_effs[i, j, l]/(1 - 2im*Δ_effs[i, j, l]))
    end
    
    # Plot effective β-factors and detunings
    Δ_index = 175
    labels = [L" ff$ = %$(ff) $" for (j, ff) in enumerate(ff_list)]
    # titl = prep_imperfectArray_transmission_title(SP) * "\narrayType: $SP.arrayType" * "\nΔ = $(SP.Δ_range[Δ_index])"
    arrayTypeDescription = SP.arrayType == "1Dchain" ? L"ordered 1D chain$$" : L"random $ z $-coordinate"
    titl = L"Array type: %$(arrayTypeDescription), detuning $ Δ = %$(round(SP.Δ_range[Δ_index], digits=2)) $"
    fig = fig_effectiveβΔ_vs_η(ηαFactor_list, β_effs[:, :, Δ_index], Δ_effs[:, :, Δ_index], labels, titl)
    # save("C:\\Users\\Simon\\Downloads\\effectiveBetaVSeta_$(SP.arrayType)_D$(round(SP.Δ_range[Δ_index], digits=2)).png", fig)
end


function plot_steadyState_radiation_Efield(SP)
    # Get steady state and its Fourier transform
    # Δ = 0.1775
    Δ = 0.326
    σBα_SS = calc_steadyState(SP, Δ)
    if SP.noPhonons σ_SS = σBα_SS else σ_SS = σBα_SS[1] end
    ks, σ_SS_FT = discFourierTransform(σ_SS, SP.a)
    
    # Get emission intensity pattern for plotting
    E = scan_radiation_Efield(SP, σBα_SS)
    intensity = norm.(E).^2
    
    # Get parameter matrices (only considers the excitation part of the Hilbert space)
    Δvari = get_Δvari(SP.ΔvariDependence, SP.Δvari_args, SP.array)
    if all(SP.ηα .== 0)
        tildeG_gm, tildeG_rm = get_tildeGs_split(SP)
    else 
        tildeG_gm, tildeGα1_gm, tildeGα2_gm, tildeG_rm, tildeGα1_rm, tildeGα2_rm = get_tildeGs_split(SP)
    end
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
        if all(SP.ηα .== 0)
            tildeG_gm, tildeG_rm = get_tildeGs_split(SP)
        else 
            tildeG_gm, tildeGα1_gm, tildeGα2_gm, tildeG_rm, tildeGα1_rm, tildeGα2_rm = get_tildeGs_split(SP)
        end
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
        tildeΩ    = get_tildeΩs(SP.fiber, SP.d, SP.incField_wlf, SP.array, SP.ΩDriveOn)
    else 
        tildeΩ, _ = get_tildeΩs(SP.fiber, SP.d, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn)
    end
    if all(SP.ηα .== 0)
        tildeG_gm, tildeG_rm = get_tildeGs_split(SP)
    else 
        tildeG_gm, tildeGα1_gm, tildeGα2_gm, tildeG_rm, tildeGα1_rm, tildeGα2_rm = get_tildeGs_split(SP)
    end
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
    # σBα_scan = scan_steadyState(SP)
    # t = calc_transmission.(Ref(SP), σBα_scan)
    # r = calc_reflection.(Ref(SP), σBα_scan)
    
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
    t = 1 .- sum(resonances)
    r = zeros(ComplexF64, size(t))
    loss, resonances_abs, resonances_abs_max, exci_pops = prep_loss_resonances_pops(t, r, resonances, eigenModesMatrix, SP.noPhonons)
    titl = prep_transmissionWithGnmEigenEnergies_title(SP)
    
    # flt = resonances_abs_max .> 0.01
    # resonances_abs = resonances_abs[flt]
    # resonances_abs_max = resonances_abs_max[flt]
    # collΔ = collΔ[flt]
    # collΓ = collΓ[flt]
    # exci_pops = exci_pops[flt]
    
    # fig_transAndResonances_polar(SP.Δ_range, t, resonances, titl)
    arrayTypeDescription = SP.arrayType == "1Dchain" ? L"ordered 1D chain$$" : L"random $ z $-coordinate"
    titl = L"$ N_{atoms} = 10 $, array type: %$(arrayTypeDescription)"
    fig_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs, collΔ, abs.(collΓ), resonances_abs_max, exci_pops, titl)
    
    # flt = weights_abs .> 1e-20
    # titl = L"$ γ_{a} = 0.01 \cdot γ_{a}^{(0)} $"
    # fig = fig_presentation_loss_withGnmeigenEnergies(SP.Δ_range, loss, resonances_abs[flt], collΔ[flt], collΓ[flt], weights_abs[flt], titl)
    # save("C:\\Users\\Simon\\Forskning\\Dokumenter\\Conferences and visits\\Berlin 2025\\talk\\figures\\loss_N10_withPhonons_e1.png", fig, px_per_unit=2)
end


function plot_memoryEfficiency(SP)
    if !SP.include3rdLevel                throw(ArgumentError("plot_memoryEfficiency assumes the third level (s) is included")) end
    if SP.initialStateDescription != "Ga" throw(ArgumentError("plot_memoryEfficiency assumes a Gaussian initial state")) end
    if SP.ΩDriveOn                        throw(ArgumentError("plot_memoryEfficiency assumes the driving on the g-e transition is off")) end
    
    # # Prepare parameters
    Δ = SP.Δc + eps(Float64)
    # fullCoupling_rm = calc_fullCoupling_rm(SP)
    radDecayRateAndStateNorm_LowerTol = (1e-6, 0.01)
    
    # # Perform time-evolution and calculate memory retrieval error
    # times, states, radDecayRatesAndStateNorm = calc_timeEvolution_forMemoryRetrievalError(SP, Δ, fullCoupling_rm, radDecayRateAndStateNorm_LowerTol)
    # radiativeDecayRates = [x[1] for x in radDecayRatesAndStateNorm]
    # ϵ = calc_memoryRetrievalError(times, radiativeDecayRates)
    
    
    N_sites_list = 10:10:200
    ϵs = []
    for N_sites in N_sites_list
        arrayDescription = arrayDescript(SP.arrayType, N_sites, SP.ρa, SP.a, SP.ff, SP.pos_unc)
        
        postfix = get_postfix_memoryEfficiency(Δ, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.incField_wlf, SP.tildeG_flags, arrayDescription, SP.fiber.postfix, SP.initialStateDescription, SP.tspan, SP.dtmax, radDecayRateAndStateNorm_LowerTol, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
        filename = "memEff_" * postfix
        folder = "memoryEfficiency/"

        if isfile(saveDir * folder * filename * ".txt")
            push!(ϵs, load_as_txt(saveDir * folder, filename)[1])
        else
            throw(ArgumentError("The following file can not be found: " * filename))
        end
    end
    
    titl = prep_memoryRetrievalError_title(SP, Δ)
    fig_memoryRetrievalError(N_sites_list, ϵs, titl)
    
    
    # TEMP
    # println(ϵ)
    
    # σgeσgsTraj = unpack_σgeσgsFromx.(states)
    # σgeσgsBαgeBαgsTraj = unpack_σgeσgsBαgeBαgsFromx.(states)
    # zs = [site[3] for site in SP.array]
    # for i in Int.(floor.(1 .+ [0.1, 0.3, 0.5, 0.7, 0.9].*length(times)))
    # # for i in Int.(floor.(1 .+ [0.5, 0.9].*length(times)))
    #     σge, σgs = σgeσgsTraj[i]
    #     # σge, σgs, Bαge, Bαgs = σgeσgsBαgeBαgsTraj[i]
    #     titl = "time = $(ro(times[i]))"
    #     fig_σgeσgs_vs_z(zs, σge, σgs, titl)
    # end
    
    # fig_complexFunction(times, radiativeDecayRates)
    # fig_complexFunction(times, norm.(states).^2)
    
    # fig_complexFunction(N_sites_list, ϵs)
    
    # TEMP
    
    
    
    # Figure out why the plots of the excitation distribution don't match at all 
    # with Chang's article (also giving us incorrect values for ϵ?)
        # Real part of GF is exact in their calculations?
    
    # FIGURE OUT where the related functions in calcs.jl should be - in calcs.jl or in physics.jl?
        # maybe this function should also be a calcs.jl?
    
    # Implement a postfix, a folder, and save ϵ
    # do the above for different N in parallel (probably quite slow for large N...)
    # and plot that
    # compare with and without phonons, different Ωc, other parameters
    
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

# Driving phonon states for high N
    # Find out whether these states are driven directly by Ωα or indirectly via single-excitation states (that is, a second order process which never populates the single-excitation modes)
        # High density of states counters the fact the jump/interaction is first (or second?) order in the Lamb-Dicke parameter
    # Find out whether T > 1 is due to bug/error or because of ηα expansion
        # Compare time evolution, steady state via basic variables, and steady state via vectorized variables
    # Find out why states with higher energy contribute more? 
        # What is the relation to driving at a certain momentum?

# Consider thermal state as a initial state (instead of total ground state)
    # Still at most one excitation and at most one phonon on top of this
    # How?

# Implement include3rdLevel calculations for 
    # Mode calculations (fullCouplingMatrix)
    # Others?
    
# Study the effect of motion on EIT (seems to be robust against motion)
    # Motion indeed includes an extra term in the transmission which may be not be zero when the original is, but its contribution is small
# Understand dark state of the system 
    # Start without motion to get a simple picture
    # A single dark state on the two-photon resonance (supposedly you also have dark states off resonance?)
    # When including motion, a perfect dark state does not exist?
    # Check steady state distribution of population right on the two-photon resonance, with and without motion, also on the two-photon+phonon resonance
# Understand group velocity around the EIT window
# Compare the classical calculation motional effects with the quantum study
    # Is there a difference? Are bad effects perhaps less when doing the quantum calculations?
    # Compare Strontium and Caesium
# Look at effect on quantum memory quality metrics
    # Implement spatial variation of Ωc
    # Implement scanning over range of Ωc? Maybe just do a few chosen points?
    # Presently initializing in a Gaussian state without phonons
        # how was that state achieved? any driving also would create phonons?
        # initialize in a Gaussian with phonons? 
        # give Bαgs same weigth as just \sigmags



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
# Consider order of for loops in packing/unpacking x and maybe other places (calculation of tildeGs?)
    # Should be faster if order of indices in loops is opposite order of indices in array (column-major something something)

# If we use "direct" for loops instead of matrix products, etc., everywhere, we can implement Δvari as a Vector instead of a Matrix?

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