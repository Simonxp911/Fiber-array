
using OrdinaryDiffEq #solve differential equation for time-evolution or steady state
using NonlinearSolve #addition to OrdinaryDiffEq for nonlinear EoMs
using Integrals #for computing integrals
using LinearAlgebra #norm of vectors and other standard linear algebra
using JLD2 #saving and loading
using DelimitedFiles #read/write simple text data files
using Plots; pythonplot() #plot using Python-Matplotlib as backend
using Colors #for generating distinguishable colors
using LaTeXStrings #LaTeX formatting in string in plots
# using Polylogarithms #for calculating the linear array FT GF [POSSIBLY FATALLY BUGGED/DEPRECATED]
# using HomotopyContinuation #for finding fixed points (i.e. solving systems of polynomial equations)
using Bessels #Bessel functions for fiber equation and modes
# using StaticArrays #implements arrays with static size, which are faster for matrix manipulations of small matrices
using Printf #for formatting strings

include("calcs.jl")
include("prep.jl")
include("figs.jl")
include("physics.jl")
include("utility.jl")
include("save_load.jl")

const ωa = ro(2π)
const saveDir = "C:/Users/Simon/Forskning/Data/fiber_array_data/"

# ON UNITS:
# We have set μ0 = ϵ0 = ħ = 1 (and thus c = 1)
# Furthermore, we have set the driving strength Ω = d*Ein = 1 and the atomic decay rate γ = 1
# Most quantities we care about become independent of Ω (e.g. the transmission) 
# or are linear in it (we are focussing on the linear physics of the system),
# such that we can express them in terms of Ω (e.g. X/Ω)
# Likewise, γ enters as an overall scale of certain quantities (e.g. the collective energies)
# such that we can express them in terms of γ (e.g. X/γ)
# Finally, we will generally take the wavelength of the atomic transition to also be λ = 1
# This is again because it either cancels or quantities may be proportional to it



#================================================
    Main functions
================================================#
struct Fiber{T<:Real}
    radius::T                           # radius of the fiber 
    refractive_index::T                 # refractive index of the fiber material
    frequency::T                        # frequency of the guided mode chosen as drive
    propagation_constant::T             # propagation constant of that guided mode
    propagation_constant_derivative::T  # derivative of propagation constant of that guided mode at its frequency
    inside_momentum::T                  # momentum parameter pertaining to the inside of the fiber
    outside_momentum::T                 # momentum parameter pertaining to the outside of the fiber
    s_parameter::T                      # so-called s-parameter, calculated for later ease 
    normalization_constant::T           # guided mode normalization constant
    postfix::String                     # postfix for saving/loading pertaining to the fiber
    
    function Fiber(radius::T, refractive_index::T, frequency::T) where {T<:Real}
        if radius < 0 throw(ArgumentError("The fiber radius must be greater than zero.")) end
        if refractive_index < 0 throw(ArgumentError("The refractive index must be greater than zero.")) end
        if frequency < 0 throw(ArgumentError("The wavelength must be greater than zero.")) end

        κ  = calc_propConst(frequency, radius, refractive_index)
        dκ = calc_propConstDerivative(frequency, radius, refractive_index)
        h  = in_momentum(κ, frequency, refractive_index)
        q  = out_momentum(κ, frequency)
        s  = s_parameter(h, q, radius)
        C  = norm_constant(frequency, κ, h, q, s, radius)
        postfix = get_postfix(radius, refractive_index, frequency)
    
        return new{T}(radius, refractive_index, frequency, κ, dκ, h, q, s, C, postfix)
    end
end


function Base.show(io::IO, fiber::Fiber)
    println(io, "Optical fiber with parameters:")
    println(io, "Radius: $(fiber.radius)")
    println(io, "Refractive index: $(fiber.refractive_index)")
    println(io, "Angular frequency: $(fiber.frequency)")
    println(io, "Propagation constant: $(fiber.propagation_constant)")
    println(io, "Postfix: $(fiber.postfix)")
end


function define_SP()
    # The System Parameters
    
    # TODO: check that all numbers are correct, including estimates and units; write to them to confirm numbers/get updated numbers
    #       -estimates/calculations of azimuthal atomic trap frequency
    #       -estimates/calculations of the Lamb-Dicke parameters
    
    # Fiber specs from "Magic-wavelength nanofiber-based two-color dipole trap with sub-λ/2 spacing"
    λ0  = 852               #nm, guided mode wavelength, transition frequency of cs133
    ω0  = ro(2π/λ0)         #nm^-1, guided mode angular frequency
    γ0  = ro(2π*5.22e3)     #kHz, free decay rate of cs133
    ρf0 = 200               #nm, fiber radius
    n0  = 1.45              #unitless, index of refraction
    
    # Atomic array specs
    ρa0 = 550   #nm, atomic array radial coordinate
    a0  = 300   #nm, atomic array lattice constant
    
    # Trap specs
    ν0_radial    = 109 #kHz, radial atomic trap frequency
    ν0_axial     = 139 #kHz, axial atomic trap frequency
    ν0_azimuthal = 18  #kHz, azimuthal atomic trap frequency (estimated from graph)
    να0 = [ν0_radial, ν0_azimuthal, ν0_axial] #trap frequencies in a Cartesian basis (x, y, z) which matches with (radial, azimuthal, axial) if the position is taken to be on the x-axis
    
    # Recoil energy
    ħ = 1.054e-34 #m^2*kg/s, Planck's reduced constant
    cesium_mass = 132.9*1.66e-27 #kg, mass of cesium-133
    νR0 = ro(ħ*(2π/(λ0*1e-9))^2/(2*cesium_mass)*1e-3) #kHz, recoil energy of cesium atoms
    
    
    # # TEMP
    # R87_λ0 = 780 #nm
    # R87_mass = 86.9*1.66e-27 #kg
    # R87_νR = ħ*(2π/(R87_λ0*1e-9))^2/(2*R87_mass)*1e-3 #kHz (matches with table value)
    
    # R87_νz = 2π*70 #kHz
    # println(sqrt(R87_νR/R87_νz))
    # R87_νT1 = 2π*187 #kHz
    # println(sqrt(R87_νR/R87_νT1))
    # R87_νT2 = 2π*375 #kHz
    # println(sqrt(R87_νR/R87_νT2))
    # # TEMP
    
    
    # Lamb-Dicke parameters in a Cartesian basis (x, y, z)
    ηα0 = @. ro(sqrt(νR0/να0))
    
    # Unitless versions
    ρf0_ul = ro(ρf0/λ0)  #unitless version of ρf0
    ρa0_ul = ro(ρa0/λ0)  #unitless version of ρa0
    a0_ul  = ro(a0/λ0)   #unitless version of a0
    να0_ul = ro.(να0/γ0) #unitless version of να0
    
    # Initialize the fiber (the frequency of the chosen driving mode should be ωa; if ωa is used the radius should be unitless)
    fiber = Fiber(ρf0_ul, n0, ωa)
    
    # Set specs and ranges for scanning the fiber propagation constant (expects dimensionfull quantities)
    ω_specs  = (0.1*ω0, 2*ω0, 100) 
    ρf_specs = (ρf0, ρf0, 1)
    n_specs  = (n0, n0, 1)
    ω_range  = range(ω_specs...)
    ρf_range = range(ρf_specs...)
    n_range  = range(n_specs...)
    
    # Set specs and ranges for time evolution and related calculations (expects dimensionless quantities)
    Δ_specs  = (-30, 30, 100)
    Δ_range  = ro.(range(Δ_specs...))
    
    # Time spand and maximum time step allowed in time evolution
    tspan = (0, 2)
    dtmax = 0.1
    
    # Set array specs and generate array, as well as description for postfix
    N  = 1
    ρa = ρa0_ul
    a  = a0_ul
    array = get_array(N, ρa, a)
    arrayDescription = standardArrayDescription(N, ρa, a)
    
    # Prepare initial state for time evolution, as well as description for postfix
    initialState = groundstate(N)
    initialStateDescription = "gs"
    
    # Phonon bare energies, i.e. trap frequencies
    να = να0_ul #assumes an atomic array of the type (ρa, 0, z)
    
    # Lamb-Dicke parameters
    ηα = ηα0 #assumes an atomic array of the type (ρa, 0, z)
    ηα = [0., 0., 0.]
    
    # Atomic dipole moment
    # d = []
    # d = chiralDipoleMoment(fiber, ρa)
    d = "chiral"
    
    if d == "chiral" && any([site[1] != ρa || site[2] != 0 for site in array]) throw(ArgumentError("d = 'dipole' assumes an array (ρa, 0, z)")) end
    
    # Whether to approximate the real part of the transverse part of the radiation Green's function 
    # using the corresponding part of the vacuum GF, as well as whether to scale the real part of the rad. GF
    # with the local radiation decay rates
    approx_Re_Grm_trans = true
    
    return (λ0=λ0, ω0=ω0, ρf0=ρf0, n0=n0, ρa0=ρa0, a0=a0, να0=να0, νR0=νR0, ηα0=ηα0,
            fiber=fiber,
            ω_specs=ω_specs, ρf_specs=ρf_specs, n_specs=n_specs,
            ω_range=ω_range, ρf_range=ρf_range, n_range=n_range,
            Δ_specs=Δ_specs, Δ_range=Δ_range,
            tspan=tspan, dtmax=dtmax,
            N=N, ρa=ρa, a=a, array=array, arrayDescription=arrayDescription,
            initialState=initialState, initialStateDescription=initialStateDescription,
            να=να, ηα=ηα,
            d=d,
            approx_Re_Grm_trans=approx_Re_Grm_trans)
end


function main()
    # Define system parameters
    SP = define_SP()
    # printSP(SP)
    
    
    # plot_propConst_inOutMom(SP)
    # plot_σBαTrajectories_σBαSS(SP)
    plot_transmission_vs_Δ(SP)
    
    return nothing
end


#================================================
    Perform calculations
================================================#



#================================================
    Generate figures
================================================#
function plot_propConst_inOutMom(SP)
    κ = scan_propConst(SP)
    
    # Plot propConst with a specific choice of ρf and n from the ranges
    ρf_ind = 1
    n_ind  = 1
    fig_propConst_vs_ω(SP.ω_range, κ[:, ρf_ind, n_ind], SP.ρf_range[ρf_ind], SP.n_range[n_ind])
    
    # Plot the momenta inside and outside the fiber
    h = in_momentum.(κ[:, ρf_ind, n_ind], SP.ω_range, SP.n_range[n_ind])
    q = out_momentum.(κ[:, ρf_ind, n_ind], SP.ω_range)
    fig_inout_momenta_vs_ω(SP.ω_range, h, q, SP.n_range[n_ind])
end


function plot_σBαTrajectories_σBαSS(SP)
    Δ = 0.0
    σ_SS, Bα_SS = calc_σBα_steadyState(SP, Δ)
    xTrajectories = timeEvolution(SP, Δ)
    
    times, σTrajectories, BαTrajectories = prep_times_σBαTrajectories(xTrajectories, SP.N)
    fig_σBαTrajectories_σBαSS(times, σTrajectories, BαTrajectories, σ_SS, Bα_SS)
end


function plot_transmission_vs_Δ(SP)
    σBα_scan = scan_σBα_steadyState(SP)
    t = calc_transmission.(Ref(SP), σBα_scan)
    
    fig_transmission_vs_Δ(SP.Δ_range, t)
end





println("\n -- Running main() -- \n")
@time main()




# TODO list:

# Consider using StaticArrays in some places?

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