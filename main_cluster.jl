

include("preamble.jl")
const saveDir = ARGS[1] * "/fiber_array_data/"

using MPI
MPI.Init()

const comm     = MPI.COMM_WORLD
const root     = 0
const myRank   = MPI.Comm_rank(comm)
const commSize = MPI.Comm_size(comm)


# ================================================
#   Main functions
# ================================================
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
    Δ_specs = (-0.5, 0.5, 30)
    
    # Set up the spatial dependence of the detuning ("flat" (nothing), "Gaussian" (amp, edge_width), "linear" (amp, edge_width), "parabolic" (amp))
    ΔvariDependence = "flat"
    Δvari_args = -3, 50*a0_ul
    ΔvariDescription = ΔvariDescript(ΔvariDependence, Δvari_args)
    
    # Lamb-Dicke parameters
    # ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    # ηα = ηα0 .* [0.1, 0.2, 0.1]
    ηα = ηα0 * 0.4
    # ηα = [0.01, 0.01, 0.01]
    # ηα = [0., 0., 0.]
    
    # Whether phonons are excluded or not from the calculations
    noPhonons = all(ηα .== 0)
    
    # Set which kind of array to use ("1Dchain", "doubleChain", "randomZ")
    arrayType = "1Dchain"
    
    # Set number of atomic sites 
    N_sites = 5
    
    # Set filling fraction, positional uncertainty, and number of instantiations 
    ff = 0.8
    pos_unc = 0.0 #ηα0/ωa * 0.4
    n_inst = 10
    
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
    
    # Absolute tolerance in the calculations of Im_Grm_trans
    abstol_Im_Grm_trans = 1e-5
    
    # Whether to approximate transverse part of radiation GF (real part and imaginary part respectively, usually (true, false))
    approx_Grm_trans = (true, false)
    
    # Whether to save individual results (Im_Grm_trans, steady states, time evolutions)
    # save_individual_res = n_inst == 1 && ff == 1 && pos_unc == 0
    # save_individual_res = n_inst == 1
    save_individual_res = pos_unc == 0 && arrayType != "randomZ"
    
    # Ranges of z and x values to define r_fields for calculating the radiated E-field
    arrayL = (N_sites - 1)*a0_ul
    # z_range = range(-0.5*arrayL, 1.5*arrayL, 60)
    # x_range = range(ρa0_ul - 0.3*arrayL, ρa0_ul + 0.3*arrayL, 60)
    z_range = range(-10, arrayL + 10, 60)
    x_range = range(-ρf0_ul - 10, ρf0_ul + ρa0_ul + 10, 60)
    y_fix   = ρa0_ul
    
    
    return SysPar(ρf0_ul, n0, ωa,
                  Δ_specs,
                  ΔvariDependence, Δvari_args, ΔvariDescription,
                  tspan, dtmax, initialState, initialStateDescription,
                  arrayType, N_sites, ρa0_ul, a0_ul, ff, pos_unc, n_inst, array, arrayDescription, N,
                  να0_ul, ηα, noPhonons,
                  d, dDescription, incField_wlf, save_individual_res, abstol_Im_Grm_trans, approx_Grm_trans,
                  z_range, x_range, y_fix)
end


function main()
    # Define system parameters
    SP = define_SP_BerlinCS()
    # show(SP)
        
    
    
    
    # prepare_Im_Grm_trans(SP)
    # plot_imperfectArray_transmission_vs_Δ(SP)
    
    return nothing
end


# ================================================
#   Generate figures
# ================================================
function prepare_Im_Grm_trans(SP)
    totalNumberOfJobs = SP.N^2
    myStartIndex, myEndIndex = myStartIndex_and_myEndIndex(totalNumberOfJobs)
    
    index = 1
    for j in 1:SP.N, i in 1:SP.N
        if myStartIndex <= index <= myEndIndex
            if i <= j     Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[j], (0, 0), 1, true, SP.abstol) end
            for α in 1:3
                          Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[j], (1, 0), α, true, SP.abstol)
                          Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[j], (2, 0), α, true, SP.abstol)
                if i == j Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[j], (1, 1), α, true, SP.abstol) end
            end
        end
        index += 1
    end
end


function plot_imperfectArray_transmission_vs_Δ(SP)
    if SP.n_inst == 1 throw(ArgumentError("plot_imperfectArray_transmission_vs_Δ requires n_inst > 1")) end
    
    totalNumberOfJobs = SP.n_inst
    myStartIndex, myEndIndex = myStartIndex_and_myEndIndex(totalNumberOfJobs)
    
    ts = []
    for index in myStartIndex:myEndIndex
        if typeof(SP.d) == String
            σBα_scan = scan_steadyState(SP, SP.d, SP.array[index])
            push!(ts, calc_transmission.(Ref(SP), σBα_scan, Ref(SP.d), Ref(SP.array[index])))
        else
            σBα_scan = scan_steadyState(SP, SP.d[index], SP.array[index])
            push!(ts, calc_transmission.(Ref(SP), σBα_scan, Ref(SP.d[index]), Ref(SP.array[index])))
        end
    end
    
    ts = gather(ts, comm, root=root)
    
    if myRank == root
        postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.incField_wlf, n_inst, SP.arrayDescription, SP.fiber.postfix)
        filename = "T_phase" * postfix
        folder = "imperfectArray_T_phase/"
        
        # Prepare means and standard deviations of (squared) magnitudes and phases
        T_means, T_stds, phase_means, phase_stds = prep_imperfectArray_transmission(ts)
        formattedResult = vectorOfRows2Matrix([T_means, T_stds, phase_means, phase_stds])
        save_as_txt(formattedResult, saveDir * folder, filename)
    end
end


function myStartIndex_and_myEndIndex(totalNumberOfJobs)
    # First divide the totalNumberOfJobs evenly among the processes, rounded down
    # Then divide the remainingJobs among the lowest-rank processes 
    # In case totalNumberOfJobs < commSize, some processes will have zero jobs to do
    
    jobsPerProcessFloor = totalNumberOfJobs ÷ commSize
    remainingJobs = totalNumberOfJobs - jobsPerProcessFloor*commSize
    if myRank < remainingJobs
        myNumberOfJobs = jobsPerProcessFloor + 1
        myStartIndex   = 1 + myRank*myNumberOfJobs
    elseif jobsPerProcessFloor > 0
        myNumberOfJobs = jobsPerProcessFloor
        myStartIndex   = 1 + remainingJobs + myRank*myNumberOfJobs
    else
        myNumberOfJobs = 0
        myStartIndex   = 0
    end
    myEndIndex = myStartIndex + myNumberOfJobs - 1
    return myStartIndex, myEndIndex
end 





@time begin
    println("\n -- Running main() -- \n")
    main() 
    println(" -- -- ")
    println("@time of main():")
end


