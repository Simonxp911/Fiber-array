
#================================================
    Julia libraries
================================================#
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


#================================================
    Files
================================================#
include("calcs.jl")
include("prep.jl")
include("figs.jl")
include("physics.jl")
include("utility.jl")
include("save_load.jl")


#================================================
    Structures and constants
================================================#
struct Fiber
    radius::Real                           # radius of the fiber 
    refractive_index::Real                 # refractive index of the fiber material
    frequency::Real                        # frequency of the guided mode chosen as drive
    propagation_constant::Real             # propagation constant of that guided mode
    propagation_constant_derivative::Real  # derivative of propagation constant of that guided mode at its frequency
    inside_momentum::Real                  # momentum parameter pertaining to the inside of the fiber
    outside_momentum::Real                 # momentum parameter pertaining to the outside of the fiber
    s_parameter::Real                      # so-called s-parameter, calculated for later ease 
    normalization_constant::Real           # guided mode normalization constant
    postfix::String                        # postfix for saving/loading pertaining to the fiber
    
    function Fiber(radius::Real, refractive_index::Real, frequency::Real)
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
    
        return new(radius, refractive_index, frequency, κ, dκ, h, q, s, C, postfix)
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


struct SysPar
    ρf::Real                                        # Fiber radius 
    n::Real                                         # Index of refraction
    ω::Real                                         # Chosen frequency of guided mode for driving
    fiber::Fiber                                    # Fiber with the above specifications
    
    Δ_specs::Tuple{Real, Real, Int}                 # Specs for detuning for scanning time evolution, steady states, etc.
    Δ_range::AbstractRange                          # Range of detuning for scanning time evolution, steady states, etc.
    
    tspan::Tuple{Real, Real}                        # Time span (min and max value) for time evolution
    dtmax::Real                                     # Maximum allowed time step for time evolution
    initialState::Vector{<:Real}                    # Initial state for time evolution
    initialStateDescription::String                 # Description of initial state for postfix
    
    N::Int                                          # Number of atoms in array
    ρa::Union{Real, Nothing}                        # Radial coordinate of regular array of atoms (set to nothing if array is not regular)
    a::Union{Real, Nothing}                         # Lattice spacing of regular array of atoms (set to nothing if array is not regular)
    array::Vector{Vector{<:Real}}                   # Atomic array (by default initialized using the above specifications)
    arrayDescription::String                        # Description of the atomic array for postfix
    
    να::Vector{<:Real}                              # Trap frequencies, i.e. bare energies of phonons
    ηα::Vector{<:Real}                              # Lamb-Dicke parameters
    
    d::Union{Vector{<:Number}, String}              # Dipole moment of atoms
    incField_wlf::Vector{Tuple{<:Number, Int, Int}} # Vector of (weight, l, f) tuples for defining the incoming driving field
    
    # Whether to approximate the real, transverse part of the radiation Green's function 
    # using the corresponding part of the vacuum GF, as well as whether to scale the real part of the rad. GF
    # with the local radiation decay rates
    approx_Re_Grm_trans::Bool
    
    
    function SysPar(ρf::Real, n::Real, ω::Real,
                    Δ_specs::Tuple{Real, Real, Int},
                    tspan::Tuple{Real, Real}, dtmax::Real, initialState::Vector, initialStateDescription::String,
                    N::Int, ρa::Real, a::Real,
                    να::Vector, ηα::Vector,
                    d::Union{Vector, String}, incField_wlf::Vector, approx_Re_Grm_trans::Bool)
        
        fiber = Fiber(ρf, n, ω)
        Δ_range = range(Δ_specs...)
        array = get_array(N, ρa, a)
        arrayDescription = standardArrayDescription(N, ρa, a)

        return new(ρf, n, ω, fiber,
                   Δ_specs, Δ_range,
                   tspan, dtmax, initialState, initialStateDescription,
                   N, ρa, a, array, arrayDescription,
                   να, ηα,
                   d, incField_wlf, approx_Re_Grm_trans)
    end
    
    function SysPar(ρf::Real, n::Real, ω::Real,
                    Δ_specs::Tuple{Real, Real, Int},
                    tspan::Tuple{Real, Real}, dtmax::Real, initialState::Vector, initialStateDescription::String,
                    array::Vector{Vector}, arrayDescription::String,
                    να::Vector, ηα::Vector,
                    d::Union{Vector, String}, incField_wlf::Vector, approx_Re_Grm_trans::Bool)
        
        if d == "chiral" && any([site[1] != ρa || site[2] != 0 for site in array]) throw(ArgumentError("d = 'dipole' assumes an array (ρa, 0, z)")) end
        
        fiber = Fiber(ρf, n, ω)
        Δ_range = range(Δ_specs...)
        N  = length(array)
        ρa = nothing
        a  = nothing
        
        return new(ρf, n, ω, fiber,
                   Δ_specs, Δ_range,
                   tspan, dtmax, initialState, initialStateDescription,
                   N, ρa, a, array, arrayDescription,
                   να, ηα,
                   d, incField_wlf, approx_Re_Grm_trans)
    end
end


function Base.show(io::IO, SP::SysPar)
    println(io, "--- System Parameters ---")
    
    show(SP.fiber)
    println(io, "")
    
    println(io, "Specs for scan of time evolution and steady state")
    println(io, "Δ_specs: ", SP.Δ_specs)
    println(io, "")
    
    println(io, "Time span and maximum time step allowed in time evolution")
    println(io, "tspan: ", SP.tspan)
    println(io, "dtmax: ", SP.dtmax)
    println(io, "")
    
    println(io, "Description of the initial state")
    println(io, SP.initialStateDescription)
    println(io, "")
    
    println(io, "Description of the atomic array")
    println(io, SP.arrayDescription)
    println(io, "")
    
    println(io, "Trap frequencies, Lamb-Dicke parameters, and atomic dipole moment")
    println(io, "να: ", SP.να)
    println(io, "ηα: ", SP.ηα)
    if typeof(SP.d) == String println(io, "d: ", SP.d)
    else println(io, "d:  [", join(format_Complex_to_String.(SP.d), ", "), "]") end
    println(io, "")
    
    println(io, "Incoming field described in terms of weights, polarization indices, and direction indices")
    println(io, "incField_wlf: ", SP.incField_wlf)
    println(io, "")
    
    println(io, "Whether the real part, transverse part of the radiation GF has been approximated with corresponding part of the vacuum GF")
    println(io, "approx_Re_Grm_trans: ", SP.approx_Re_Grm_trans)
    println(io, "")
    
    println(io, "---  ---")
end


struct ω_ρf_n_ranges
    ω_specs::Tuple{Real, Real, Int}  # Specs for driving frequencies for scanning propagation constant etc.
    ρf_specs::Tuple{Real, Real, Int} # Specs for fiber radii for scanning propagation constant etc.
    n_specs::Tuple{Real, Real, Int}  # Specs for indices of refraction for scanning propagation constant etc.
    ω_range::AbstractRange           # Range of driving frequencies for scanning propagation constant etc.
    ρf_range::AbstractRange          # Range of fiber radii for scanning propagation constant etc.
    n_range::AbstractRange           # Range of indices of refraction for scanning propagation constant etc.
    
    function ω_ρf_n_ranges(ω_specs::Tuple{Real, Real, Int}, ρf_specs::Tuple{Real, Real, Int}, n_specs::Tuple{Real, Real, Int})
        ω_range  = range(ω_range...)
        ρf_range = range(ρf_specs...)
        n_range  = range(n_specs...)
        return new(ω_specs, ρf_specs, n_specs, ω_range, ρf_range, n_range)
    end
end


const ωa = 2π
const saveDir = "C:/Users/Simon/Forskning/Data/fiber_array_data/"


#================================================
    Regarding units
================================================#
# We have set μ0 = ϵ0 = ħ = 1 (and thus c = 1)
# Furthermore, we have set the driving strength Ω = d*Ein = 1 and the atomic decay rate γ = 1
# Most quantities we care about become independent of Ω (e.g. the transmission) 
# or are linear in it (we are focussing on the linear physics of the system),
# such that we can express them in terms of Ω (e.g. X/Ω)
# Likewise, γ enters as an overall scale of certain quantities (e.g. the collective energies)
# such that we can express them in terms of γ (e.g. X/γ)
# Finally, we will generally take the wavelength of the atomic transition to also be λ = 1
# This is again because it either cancels or quantities may be proportional to it
