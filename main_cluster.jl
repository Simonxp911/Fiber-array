

include("preamble.jl")
BLAS.set_num_threads(1) # Necessary on the workstation for LinearAlgebra to work properly
const saveDir = ARGS[1] * "/fiber_array_data/"

using MPI
MPI.Init()
const comm     = MPI.COMM_WORLD
const root     = 0
const myRank   = MPI.Comm_rank(comm)
const commSize = MPI.Comm_size(comm)

# Include the input file (which defines define_SP())
include(ARGS[2])


# ================================================
#   Main functions
# ================================================
function main()
    # Define system parameters
    SP = define_SP()
    MPI.Barrier(comm)
    if myRank == root show(SP) end
    MPI.Barrier(comm)
    
    
    task = ARGS[3]
    if task == "prepare_Im_Grm_trans"
        prepare_Im_Grm_trans(SP)
    elseif task == "imperfectArray_transmission_vs_Δ"
        imperfectArray_transmission_vs_Δ(SP)
    elseif task == "steadyState_vs_Δ"
        steadyState_vs_Δ(SP)
    elseif task == "memoryEfficiency"
        memoryEfficiency(SP)
    end
    return nothing
end


# ================================================
#   Generate figures
# ================================================
function prepare_Im_Grm_trans(SP)
    # For 1Dchain it is enough to only consider the interaction between the first atom and the others
    if SP.arrayType == "1Dchain"
        totalNumberOfJobs = SP.N
        myStartIndex, myEndIndex = myStartIndex_and_myEndIndex(totalNumberOfJobs)
        index = 1
        for i in 1:SP.N
            if myStartIndex <= index <= myEndIndex
                println("My rank is $myRank, I am working on index = $index")
                              Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[1], (0, 0), 1, true, SP.abstol_Im_Grm_trans)
                for α in 1:3
                              Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[1], (1, 0), α, true, SP.abstol_Im_Grm_trans)
                              Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[1], (0, 1), α, true, SP.abstol_Im_Grm_trans)
                              Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[1], (2, 0), α, true, SP.abstol_Im_Grm_trans)
                              Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[1], (0, 2), α, true, SP.abstol_Im_Grm_trans)
                    if i == 1 Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[1], (1, 1), α, true, SP.abstol_Im_Grm_trans) end
                end
            end
            index += 1
        end
    
    else
        totalNumberOfJobs = SP.N^2
        myStartIndex, myEndIndex = myStartIndex_and_myEndIndex(totalNumberOfJobs)
        index = 1
        for j in 1:SP.N, i in 1:SP.N
            if myStartIndex <= index <= myEndIndex
                println("My rank is $myRank, I am working on index = $index")
                if i <= j     Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[j], (0, 0), 1, true, SP.abstol_Im_Grm_trans) end
                for α in 1:3
                              Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[j], (1, 0), α, true, SP.abstol_Im_Grm_trans)
                              Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[j], (2, 0), α, true, SP.abstol_Im_Grm_trans)
                    if i == j Im_Grm_trans(SP.fiber, ωa, SP.array[i], SP.array[j], (1, 1), α, true, SP.abstol_Im_Grm_trans) end
                end
            end
            index += 1
        end
    end
end


function imperfectArray_transmission_vs_Δ(SP)
    if SP.n_inst == 1 throw(ArgumentError("imperfectArray_transmission_vs_Δ requires n_inst > 1")) end
    
    totalNumberOfJobs = SP.n_inst
    myStartIndex, myEndIndex = myStartIndex_and_myEndIndex(totalNumberOfJobs)
    
    ts = []
    for index in myStartIndex:myEndIndex
        println("My rank is $myRank, I am working on index = $index")
        if typeof(SP.d) == String
            σBα_scan = scan_steadyState(SP, SP.d, SP.array[index])
            push!(ts, calc_transmission.(Ref(SP), σBα_scan, Ref(SP.d), Ref(SP.array[index])))
        else
            σBα_scan = scan_steadyState(SP, SP.d[index], SP.array[index])
            push!(ts, calc_transmission.(Ref(SP), σBα_scan, Ref(SP.d[index]), Ref(SP.array[index])))
        end
    end
    
    ts = MPI.gather(ts, comm, root=root)
    
    if myRank == root
        postfix = get_postfix_imperfectArray_transmission(SP.Δ_specs, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.noPhonons, SP.incField_wlf, SP.n_inst, SP.tildeG_flags, SP.arrayDescription, SP.fiber.postfix)
        filename = "T_phase_" * postfix
        folder = "imperfectArray_T_phase/"
        
        # Prepare means and standard deviations of (squared) magnitudes and phases
        ts = vcat(ts...)
        formattedResult = vectorOfRows2Matrix([prep_imperfectArray_transmission(ts)...])
        save_as_txt(formattedResult, saveDir * folder, filename)
    end
end


function steadyState_vs_Δ(SP)
    if myRank == root
        postfixes = get_postfix_steadyState.(SP.Δ_range, SP.ΔvariDescription, SP.dDescription, Ref(SP.να), Ref(SP.ηα), SP.noPhonons, Ref(SP.incField_wlf), Ref(SP.tildeG_flags), SP.arrayDescription, SP.fiber.postfix)
        params    = get_parameterMatrices(SP)
    else
        postfixes = nothing
        params    = nothing
    end
    
    postfixes = MPI.bcast(postfixes, root, comm)
    params    = MPI.bcast(params   , root, comm)
    
    totalNumberOfJobs = SP.Δ_specs[3]
    myStartIndex, myEndIndex = myStartIndex_and_myEndIndex(totalNumberOfJobs)
    for index in myStartIndex:myEndIndex
        println("My rank is $myRank, I am working on index = $index")
        calc_steadyState(SP.Δ_range[index], params, postfixes[index])
    end
end


function memoryEfficiency(SP)
    if !SP.include3rdLevel                throw(ArgumentError("plot_memoryEfficiency assumes the third level (s) is included")) end
    if SP.initialStateDescription != "Ga" throw(ArgumentError("plot_memoryEfficiency assumes a Gaussian initial state")) end
    if SP.ΩDriveOn                        throw(ArgumentError("plot_memoryEfficiency assumes the driving on the g-e transition is off")) end
    
    # Prepare parameters
    Δ = SP.Δc + eps(Float64)
    fullCoupling_rm = calc_fullCoupling_rm(SP)
    radDecayRateAndStateNorm_LowerTol = (1e-6, 0.01)
    
    # Perform time-evolution and calculate memory retrieval error
    times, states, radDecayRatesAndStateNorm = calc_timeEvolution_forMemoryRetrievalError(SP, Δ, fullCoupling_rm, radDecayRateAndStateNorm_LowerTol)
    radiativeDecayRates = [x[1] for x in radDecayRatesAndStateNorm]
    ϵ = calc_memoryRetrievalError(times, radiativeDecayRates)
    
    
    postfix = get_postfix_memoryEfficiency(Δ, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.noPhonons, SP.incField_wlf, SP.tildeG_flags, SP.arrayDescription, SP.fiber.postfix, SP.initialStateDescription, SP.tspan, SP.dtmax, radDecayRateAndStateNorm_LowerTol, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
    filename = "memEff_" * postfix
    folder = "memoryEfficiency/"
    save_as_txt(ϵ, saveDir * folder, filename)
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
    println("\n -- Running main() on rank $myRank -- \n")
    main() 
    println(" -- -- ")
    println("@time of main() on rank $myRank:")
end


