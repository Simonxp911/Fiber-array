

# ================================================
#   Julia libraries
# ================================================
using OrdinaryDiffEq #solve differential equation for time-evolution or steady state
using NonlinearSolve #addition to OrdinaryDiffEq for nonlinear EoMs
using Integrals #for computing integrals
using LinearAlgebra #norm of vectors and other standard linear algebra
using JLD2 #saving and loading
using DelimitedFiles #read/write simple text data files
# using Plots; pythonplot() #plot using Python-Matplotlib as backend
# using CairoMakie #plotting, specialized for making png/pdf output
using Colors #for generating distinguishable colors
using LaTeXStrings #LaTeX formatting in string in plots
# using Polylogarithms #for calculating the linear array FT GF [POSSIBLY FATALLY BUGGED/DEPRECATED]
# using HomotopyContinuation #for finding fixed points (i.e. solving systems of polynomial equations)
using Bessels #Bessel functions for fiber equation and modes
# using StaticArrays #implements arrays with static size, which are faster for matrix manipulations of small matrices
using Printf #for formatting strings
using Random #for random permutations etc.
using Statistics #for calculating mean, standard deviation, etc.


# ================================================
#   Files
# ================================================
include("calcs.jl")
include("prep.jl")
include("figs.jl")
include("physics.jl")
include("utility.jl")
include("save_load.jl")


# ================================================
#   Structures and constants
# ================================================
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
        C  = norm_constant(refractive_index, κ, h, q, s, radius)
        postfix = get_postfix_Fiber(radius, refractive_index, frequency)
    
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
    
    ΔvariDependence::String                         # String identifier for the site-dependent detuning variation
    Δvari_args::Union{Tuple, Nothing}               # Arguments for the calculation of the site-dependent detuning variation
    ΔvariDescription::String                        # Description of the above for postfix
    
    tspan::Tuple{Real, Real}                        # Time span (min and max value) for time evolution
    dtmax::Real                                     # Maximum allowed time step for time evolution
    initialState::Vector{<:Real}                    # Initial state for time evolution
    initialStateDescription::String                 # Description of initial state for postfix
    
    arrayType::String                               # String defining the type of array used
    N_sites::Int                                    # Number of sites in array (which may or may not hold an atom)
    ρa::Union{Real, Nothing}                        # Radial coordinate of regular array of atoms (set to nothing if array is not regular)
    a::Union{Real, Nothing}                         # Lattice spacing of regular array of atoms (set to nothing if array is not regular)
    ff::Union{Real, Nothing}                        # Filling fraction of the array
    pos_unc::Union{Real, Vector, Nothing}           # Classical positional uncertainty
    n_inst::Int                                     # Number of instantiations of the array, for calculations involving random instantiations
    array::Vector{Vector}                           # Atomic array or list of different atomic array instantiations (by default initialized using the above specifications)
    arrayDescription::String                        # Description of the atomic array for postfix
    N::Int                                          # Number of atoms in array
    
    να::Vector{<:Real}                              # Trap frequencies, i.e. bare energies of phonons
    ηα::Vector{<:Real}                              # Lamb-Dicke parameters
    noPhonons::Bool                                 # Whether phonons are excluded or not from the calculations
    
    d::Union{Vector, String}                        # Dipole moment of atoms (one for each atom)
    dDescription::String                            # Description of the dipole moment for postfix
    incField_wlf::Vector{Tuple{<:Number, Int, Int}} # Vector of (weight, l, f) tuples for defining the incoming driving field
    
    save_Im_Grm_trans::Bool                         # Whether to save calculated Im_Grm_trans
    abstol_Im_Grm_trans::Real                       # Absolute tolerance in the calculations of Im_Grm_trans
    
    # Whether to approximate the transverse part of the radiation Green's function 
    # using the corresponding part of the vacuum GF, as well as whether to scale the real part of the rad. GF
    # with the local radiation decay rates
    # The the two booleans determine whether to approximate the real and imaginary part respectively.
    approx_Grm_trans::Tuple{Bool, Bool}
    
    save_steadyState::Bool                          # Whether to save calculated steady states
    save_timeEvol::Bool                             # Whether to save calculated time evolutions
    
    z_range::Union{AbstractRange, Nothing}          # Range of z values for calculating the radiation E-field
    x_range::Union{AbstractRange, Nothing}          # Range of x values for calculating the radiation E-field
    y_fix::Union{Real, Nothing}                                     # Fixed value of y for calculating the radiation E-field
    r_fields::Union{Matrix{Vector{<:Real}}, Nothing} # Range of positions for calculating the radiation E-field
    
    
    function SysPar(ρf::Real, n::Real, ω::Real,
                    Δ_specs::Tuple{Real, Real, Int},
                    ΔvariDependence::String, Δvari_args::Union{Tuple, Nothing}, ΔvariDescription::String,
                    tspan::Tuple{Real, Real}, dtmax::Real, initialState::Vector, initialStateDescription::String,
                    arrayType::String, N_sites::Int, ρa::Real, a::Real, ff::Real, pos_unc::Union{Real, Vector}, n_inst::Int, array::Vector, arrayDescription::String, N::Int,
                    να::Vector, ηα::Vector, noPhonons::Bool,
                    d::Union{Vector, String}, dDescription::String, incField_wlf::Vector, save_Im_Grm_trans::Bool, abstol_Im_Grm_trans::Real, approx_Grm_trans::Tuple,
                    save_steadyState::Bool, save_timeEvol::Bool, 
                    z_range::AbstractRange, x_range::AbstractRange, y_fix::Real)

        fiber = Fiber(ρf, n, ω)
        Δ_range = range(Δ_specs...)
        r_fields = [[x, y_fix, z] for z in z_range, x in x_range]
        
        return new(ρf, n, ω, fiber,
                   Δ_specs, Δ_range,
                   ΔvariDependence, Δvari_args, ΔvariDescription,
                   tspan, dtmax, initialState, initialStateDescription,
                   arrayType, N_sites, ρa, a, ff, pos_unc, n_inst, array, arrayDescription, N,
                   να, ηα, noPhonons,
                   d, dDescription, incField_wlf, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                   save_steadyState, save_timeEvol, 
                   z_range, x_range, y_fix, r_fields)
    end
    
    function SysPar(ρf::Real, n::Real, ω::Real,
                    Δ_specs::Tuple{Real, Real, Int},
                    ΔvariDependence::String, Δvari_args::Union{Tuple, Nothing}, ΔvariDescription::String,
                    tspan::Tuple{Real, Real}, dtmax::Real, initialState::Vector, initialStateDescription::String,
                    arrayType::String, N::Int, ρa::Real, a::Real,
                    να::Vector, ηα::Vector, noPhonons::Bool,
                    d::Union{Vector, String}, dDescription::String, incField_wlf::Vector, approx_Grm_trans::Tuple)
        
        fiber = Fiber(ρf, n, ω)
        Δ_range = range(Δ_specs...)
        ff = 1.0
        pos_unc = 0.0
        n_inst = 1
        array = get_array(arrayType, N, ρa, a, ff, pos_unc)
        arrayDescription = arrayDescript(arrayType, N_sites, ρa, a, ff, pos_unc)
        abstol_Im_Grm_trans = 1e-4
        save_Im_Grm_trans = true
        save_steadyState = true
        save_timeEvol = true
        z_range = nothing
        x_range = nothing
        y_fix   = nothing
        r_fields = nothing

        return new(ρf, n, ω, fiber,
                   Δ_specs, Δ_range,
                   ΔvariDependence, Δvari_args, ΔvariDescription,
                   tspan, dtmax, initialState, initialStateDescription,
                   arrayType, N, ρa, a, ff, pos_unc, n_inst, array, arrayDescription, N,
                   να, ηα, noPhonons,
                   d, dDescription, incField_wlf, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans,
                   save_steadyState, save_timeEvol, 
                   z_range, x_range, y_fix, r_fields)
    end
end


function Base.show(io::IO, SP::SysPar)
    println(io, "--- System Parameters ---")
    
    show(SP.fiber)
    println(io, "")
    
    println(io, "Specs for scan of time evolution and steady state")
    println(io, "Δ_specs: ", SP.Δ_specs)
    println(io, "")
    
    println(io, "Parameters of site-dependent detuning")
    println(io, "ΔvariDependence: ", SP.ΔvariDependence)
    println(io, "Δvari_args: ", SP.Δvari_args)
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
    println(io, "d: ", SP.dDescription)
    println(io, "")
    
    println(io, "Incoming field described in terms of weights, polarization indices, and direction indices")
    println(io, "incField_wlf: ", "[" * join(["(" * join([format_Complex_to_String(wlf[1]), wlf[2], wlf[3]], ", ") * ")" for wlf in SP.incField_wlf], ", ") * "]")
    println(io, "")
    
    println(io, "Whether to save the imaginary transverse part of the radiation Green's function")
    println(io, "save_Im_Grm_trans: ", SP.save_Im_Grm_trans)
    println(io, "")
    
    println(io, "Whether the real part, transverse part of the radiation GF has been approximated with corresponding part of the vacuum GF")
    println(io, "approx_Grm_trans: ", SP.approx_Grm_trans)
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


# ================================================
#   Regarding units
# ================================================
# We have set μ0 = ϵ0 = ħ = 1 (and thus c = 1)
# Furthermore, we have set the driving strength Ω = d*Ein = 1 and the atomic decay rate γ = 1
# Most quantities we care about become independent of Ω (e.g. the transmission) 
# or are linear in it (we are focussing on the linear physics of the system),
# such that we can express them in terms of Ω (e.g. X/Ω)
# Likewise, γ enters as an overall scale of certain quantities (e.g. the collective energies)
# such that we can express them in terms of γ (e.g. X/γ)
# Finally, we will generally take the wavelength of the atomic transition to also be λ = 1
# This is again because it either cancels or quantities may be proportional to it
