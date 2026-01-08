

# ================================================
#   Functions pertaining to atomic array itself
# ================================================
"""
Calculate a list of the atomic positions along the fiber,
possibly including imperfect filling fraction and (classical) positional uncertainty

Implemented arrayType's are\n
1Dchain: Standard 1D chain along fiber\n
doubleChain: Double 1D chain, positioned on opposite sides of the fiber\n
randomZ: Standard 1D chain along fiber, but with random position in z\n

For random arrays, the lattice spacing, a, defines the length along z of the array via L = (N - 1)a
"""
function get_array(arrayType, N_sites, ρa, a, ff=1, pos_unc=0)
    # Get the perfect array (no missing atoms, no position uncertainty)
    if arrayType == "1Dchain"
        array = [[ρa, 0, a*n] for n in 0:N_sites-1]
        
    elseif arrayType == "doubleChain"
        if !iseven(N_sites) throw(ArgumentError("The number of sites, N_sites = $N_sites, must be even for arrayType = doubleChain")) end
        arrayAbove = [[ ρa, 0, a*n] for n in 0:N_sites/2-1]
        arrayBelow = [[-ρa, 0, a*n] for n in 0:N_sites/2-1]
        array = vcat(arrayAbove, arrayBelow)
        
    elseif arrayType == "randomZ"
        array = [[ρa, 0, zn] for zn in a*(N_sites - 1)*sort(rand(N_sites))]
        
    else throw(ArgumentError("The arrayType = $arrayType has not been implemented in get_array")) 
    end
    
    if      ff != 1 array = remove_atoms_from_array(array, ff) end
    if pos_unc != 0 array = introduce_position_uncertainty_to_array_sites(array, pos_unc) end
    
    return array
end


"""
Returns the atomic array, or a list of different atomic array instantiations,
as well as the arrayDescription and the number of atoms 
"""
function get_array(arrayType, N_sites, ρa, a, ff, pos_unc, n_inst)
    arrayDescription = arrayDescript(arrayType, N_sites, ρa, a, ff, pos_unc)
    
    if n_inst == 1 
        array = get_array(arrayType, N_sites, ρa, a, ff, pos_unc)
        N = length(array)
        return array, arrayDescription, N
    else 
        array = [get_array(arrayType, N_sites, ρa, a, ff, pos_unc) for i in 1:n_inst] 
        N = length(array[1])
        return array, arrayDescription, N
    end
end


"""
Create imperfectly filled array according to the filling fraction ff
"""
function remove_atoms_from_array(array, ff)
    # Find total number of atoms and the number of atoms to be kept to match the desired filling fraction
    N = length(array)
    N_to_be_kept = Int(floor(N*ff))
    
    # Return N_to_be_kept of the original array sites
    return array[sort(randperm(N)[1:N_to_be_kept])]
end


"""
Create a (classically) disordered array according to a Gaussian distribution with width pos_unc
"""
function introduce_position_uncertainty_to_array_sites(array, pos_unc)
    # Generate normally-distributed random numbers for each coordinate of each atom
    N = length(array)
    # random_shift = pos_unc*randn.(fill(3, N))
    random_shift = broadcast(.*, [pos_unc], randn.(fill(3, N))) # works both for the cases of pos_unc=scalar and pos_unc=triplet
    
    # Return the randomly shifted array sites
    return array + random_shift
end


# ================================================
#   Functions pertaining to parameters of the fiber
# ================================================
"""
The momentum h related to the inside of the fiber, used in several fiber-related equations
"""
function in_momentum(κ, ω, n)
    return sqrt(n^2*ω^2 - κ^2)
end


"""
The momentum h related to the outside of the fiber, used in several fiber-related equations
"""
function out_momentum(κ, ω)
    # By taking the absolute value, the function is useful for both guided and radiation modes
    return sqrt(abs(κ^2 - ω^2))
end


"""
The s-parameter, used in several fiber-related equations
"""
function s_parameter(h, q, ρf)
    return (1/(h*ρf)^2 + 1/(q*ρf)^2)/(dbesselj(1, 1, h*ρf)/(h*ρf*besselj(1, h*ρf)) + dbesselk(1, 1, q*ρf)/(q*ρf*besselk(1, q*ρf)))
end


"""
The fiber guided mode normalization constant
"""
function norm_constant(n, κ, h, q, s, ρf)
    J0 = besselj(0, h*ρf)
    J1 = besselj(1, h*ρf)
    J2 = besselj(2, h*ρf)
    J3 = besselj(3, h*ρf)
    K0 = besselk(0, q*ρf)
    K1 = besselk(1, q*ρf)
    K2 = besselk(2, q*ρf)
    K3 = besselk(3, q*ρf)
    
    C_in_1 = (1 - s)^2*(J0^2 + J1^2) 
    C_in_2 = (1 + s)^2*(J2^2 - J1*J3) 
    C_in_3 = 2*(h/κ)^2*(J1^2 - J0*J2) 
    C_in = (n*q*K1/(h*J1))^2*(C_in_1 + C_in_2 + C_in_3)

    C_out_1 = (1 - s)^2*(K0^2 - K1^2)
    C_out_2 = (1 + s)^2*(K2^2 - K1*K3)
    C_out_3 = 2*(q/κ)^2*(K1^2 - K0*K2)
    C_out = -(C_out_1 + C_out_2 + C_out_3)

    return 1/sqrt(2π*ρf^2*(C_in + C_out))
end


"""
Implements the fiber eigenequation with all terms are moved to left, 
meaning the propagation constant is found as the roots of this.

Signature is chosen to conform with the requirements of NonlinearSolve.

params = ω, ρf, n
"""
function fiber_equation(x, params)
    ω, ρf, n = params
    
    h   = in_momentum(x, ω, n)
    q   = out_momentum(x, ω)
    hρf = h*ρf
    qρf = q*ρf
    
    J0  = dbesselj(0, 0, hρf)
    J1  = dbesselj(0, 1, hρf)
    K1  = dbesselk(0, 1, qρf)
    K1p = dbesselk(1, 1, qρf)
    
    A  = J0/(hρf*J1)
    B  = (n^2 + 1)/(2*n^2) * K1p/(qρf*K1)
    C  = -1/hρf^2
    D1 = ((n^2 - 1)/(2*n^2) * K1p/(qρf*K1))^2
    D2 = x^2/(n^2*ω^2) * (1/qρf^2 + 1/hρf^2)^2
    D  = sqrt(D1 + D2)
    return A + B + C + D
    
    
    # # For small arguments the Bessel functions are replaced by their asymptotic expressions,
    # # and the resulting products are reduced to avoid division by ~zero
    # tol = 1e-4
    # if hρf > tol
    #     J0 = besselj0(hρf)
    #     J1 = besselj1(hρf)
    #     A = hρf*qρf^2*J0/J1
    # else
    #     A = qρf^2*2
    # end
    
    # if qρf > tol
    #     K1  = besselk1(qρf)
    #     K1p = -(besselk0(qρf) + besselk(2, qρf))/2
    #     B  =  (n^2 + 1)/(2*n^2) * hρf^2*qρf*K1p/K1
    #     D1 = ((n^2 - 1)/(2*n^2) * hρf^2*qρf*K1p/K1)^2
    # else
    #     B  =  (n^2 + 1)/(2*n^2) * hρf^2*(-1)
    #     D1 = ((n^2 - 1)/(2*n^2) * hρf^2*(-1))^2
    # end
    
    # C = -qρf^2
    # D2 = x^2/(n^2*ω^2) * (qρf^2 + hρf^2)^2
    # D = sqrt(D1 + D2)
    # return A + B + C + D
end


"""
Calculates the dipole moment which yields a chiral fiber setup
"""
function chiralDipoleMoment(fiber, ρa)
    eρ, _, ez = guidedModeComps(fiber, ρa) 
    return 1im*[ez, 0, -eρ]/sqrt(abs2(eρ) + abs2(ez))    
end


"""
Returns a list of the dipole moments of each atom which yields a chiral fiber setup
for the case of 1Dchain or doubleChain arrays
"""
function chiralDipoleMoment(fiber, ρa, array::Vector{<:Vector})
    return [chiralDipoleMoment(fiber, ρa).*[sign(site[1]), 1, 1] for site in array]    
end


"""
Returns a list of the dipole moments of each atom which yields a chiral fiber setup
for the case of 1Dchain or doubleChain arrays, for the case of n_inst != 1
"""
function chiralDipoleMoment(fiber, ρa, array::Vector{<:Vector{<:Vector}})
    return [chiralDipoleMoment(fiber, ρa, array_inst) for array_inst in array]
end


"""
Returns a right-circular dipole moment vector
"""
function rightCircularDipoleMoment()
    return [1im, 0, 1]/sqrt(2)
end


"""
Returns a list of right-circular dipole moments 
"""
function rightCircularDipoleMoment(array::Vector{<:Vector})
    return [rightCircularDipoleMoment() for site in array]    
end


"""
Returns a list of right-circular dipole moments, for the case of n_inst != 1
"""
function rightCircularDipoleMoment(array::Vector{<:Vector{<:Vector}})
    return [rightCircularDipoleMoment(array_inst) for array_inst in array]
end


# ================================================
#   Functions pertaining to guided modes of the fiber
# ================================================
"""
Calculates the guided mode or its derivatives with
l: transverse polarization index
f: forward/backward index
r: field position
derivOrder: order of derivative
α: index indicating which derivative to take (α = 1,2,3 corresponds to x,y,z)
"""
function Egm(fiber, l, f, r, derivOrder=0, α=1)
    if l ∉ (-1, 1) || f ∉ (-1, 1) throw(ArgumentError("The indices l and f must be either -1 or +1")) end
    
    κ = fiber.propagation_constant
    ρ, ϕ, z = cylCoordinates(r)
    ρ_unit, ϕ_unit, z_unit = cylUnitVectors(r)
    
    # Put together the mode
    eρ, eϕ, ez = guidedModeComps(fiber, ρ)
    eμ = eρ*ρ_unit + l*eϕ*ϕ_unit + f*ez*z_unit
    Egm = eμ*exp(1im*l*ϕ)*exp(1im*f*κ*z)
    if derivOrder == 0 return Egm end
    
    # Calculate the ρ and ϕ derivatives if necessary
    if α ∈ (1, 2)
        if derivOrder >= 1
            # The ρ derivatives are found simply through the derivatives of eμ
            eρ_ρ, eϕ_ρ, ez_ρ = guidedModeComps(fiber, ρ, 1)
            eμ_ρ = eρ_ρ*ρ_unit + l*eϕ_ρ*ϕ_unit + f*ez_ρ*z_unit
            Egm_ρ = eμ_ρ*exp(1im*l*ϕ)*exp(1im*f*κ*z)
            
            # The ϕ derivatives can be found from the mode itself
            Egm_ϕ = 1im*l*Egm + (ρ_unit'*Egm)*ϕ_unit - (ϕ_unit'*Egm)*ρ_unit
        end
        if derivOrder >= 2
            eρ_ρρ, eϕ_ρρ, ez_ρρ = guidedModeComps(fiber, ρ, 2)
            eμ_ρρ = eρ_ρρ*ρ_unit + l*eϕ_ρρ*ϕ_unit + f*ez_ρρ*z_unit
            Egm_ρρ = eμ_ρρ*exp(1im*l*ϕ)*exp(1im*f*κ*z)
            
            # The ϕ derivatives can be found from the mode itself
            Egm_ϕϕ = -2*Egm + (z_unit'*Egm)*z_unit + 2im*l*((ρ_unit'*Egm)*ϕ_unit - (ϕ_unit'*Egm)*ρ_unit)
            
            # The mixed double derivative is found by combining the above relations
            Egm_ρϕ = 1im*l*Egm_ρ + (ρ_unit'*Egm_ρ)*ϕ_unit - (ϕ_unit'*Egm_ρ)*ρ_unit
        end
    end
    
    # The Cartezian derivatives are calculated from the cylindrical derivatives
    if derivOrder == 1
        if     α == 1 return cos(ϕ)*Egm_ρ - sin(ϕ)/ρ*Egm_ϕ
        elseif α == 2 return sin(ϕ)*Egm_ρ + cos(ϕ)/ρ*Egm_ϕ
        elseif α == 3 return 1im*f*κ*Egm
        end
    end
    if derivOrder == 2
        if     α == 1 return cos(ϕ)^2*Egm_ρρ + sin(ϕ)^2/ρ^2*Egm_ϕϕ - 2*cos(ϕ)*sin(ϕ)/ρ*Egm_ρϕ + sin(ϕ)^2/ρ*Egm_ρ + 2*cos(ϕ)*sin(ϕ)/ρ^2*Egm_ϕ
        elseif α == 2 return sin(ϕ)^2*Egm_ρρ + cos(ϕ)^2/ρ^2*Egm_ϕϕ + 2*cos(ϕ)*sin(ϕ)/ρ*Egm_ρϕ + cos(ϕ)^2/ρ*Egm_ρ - 2*cos(ϕ)*sin(ϕ)/ρ^2*Egm_ϕ
        elseif α == 3 return -κ^2*Egm
        end
    end
    if derivOrder > 2 throw(ArgumentError("Egm is not implemented for derivOrder > 2")) end
end


"""
Calculates the individual components of the guided mode or their derivatives.
"""
function guidedModeComps(fiber, ρ, derivOrder=0)
    # Set up some parameters
    κ, ρf = fiber.propagation_constant, fiber.radius
    h, q  = fiber.inside_momentum, fiber.outside_momentum
    s, C  = fiber.s_parameter, fiber.normalization_constant
    
    # Put together the components    
    if ρ < ρf
        eρ = 1im*C*q/h*dbesselk(0, 1, q*ρf)/dbesselj(0, 1, h*ρf)*h^derivOrder*((1 - s)*dbesselj(derivOrder, 0, h*ρ) - (1 + s)*dbesselj(derivOrder, 2, h*ρ))
        eϕ =    -C*q/h*dbesselk(0, 1, q*ρf)/dbesselj(0, 1, h*ρf)*h^derivOrder*((1 - s)*dbesselj(derivOrder, 0, h*ρ) + (1 + s)*dbesselj(derivOrder, 2, h*ρ))
        ez =   C*2*q/κ*dbesselk(0, 1, q*ρf)/dbesselj(0, 1, h*ρf)*h^derivOrder*dbesselj(derivOrder, 1, h*ρ)
    else
        eρ = 1im*C*q^derivOrder*((1 - s)*dbesselk(derivOrder, 0, q*ρ) + (1 + s)*dbesselk(derivOrder, 2, q*ρ))
        eϕ =    -C*q^derivOrder*((1 - s)*dbesselk(derivOrder, 0, q*ρ) - (1 + s)*dbesselk(derivOrder, 2, q*ρ))
        ez = C*q^derivOrder*2*q/κ*dbesselk(derivOrder, 1, q*ρ)
    end
    return eρ, eϕ, ez
end


"""
Calculates the guided mode Green's function or its derivatives 
(ignoring the step-function when it comes to derivatives, 
i.e. the step-function is never differentiated)
"""
function Ggm(fiber, r_field, r_source, derivOrder=(0, 0), α=1)
    if r_field[3] - r_source[3] == 0 && norm(r_field - r_source) != 0
        return 1im/(2*fiber.frequency)*fiber.propagation_constant_derivative*
                sum([heaviside(f*eps(r_field[3] - r_source[3]))*
                     Egm(fiber, l, f, r_field , derivOrder[1], α)*
                     Egm(fiber, l, f, r_source, derivOrder[2], α)' for l in (1, -1), f in (1, -1)])
    else
        return 1im/(2*fiber.frequency)*fiber.propagation_constant_derivative*
                sum([heaviside(f*(r_field[3] - r_source[3]))*
                     Egm(fiber, l, f, r_field , derivOrder[1], α)*
                     Egm(fiber, l, f, r_source, derivOrder[2], α)' for l in (1, -1), f in (1, -1)])
    end
end


# ================================================
#   Functions pertaining to radiation modes of the fiber
# ================================================
"""
Calculates the radiation modes or their derivatives 

Note that here, ω and κ can be any frequency and z-momentum of light
and are not necessarily related to the fiber-specific frequency and propagation constant
"""
function Erm(fiber, ω, κ, m, l, r, derivOrder=0, α=1)
    ρ, ϕ, z = cylCoordinates(r)
    ρ_unit, ϕ_unit, z_unit = cylUnitVectors(r)
    
    # Put together the mode
    eρ, eϕ, ez = radiationModeComps(fiber, ω, κ, m, l, ρ)
    eν = eρ*ρ_unit + eϕ*ϕ_unit + ez*z_unit
    Erm = eν*exp(1im*m*ϕ)*exp(1im*κ*z)
    if derivOrder == 0 return Erm end
    
    # Calculate the ρ and ϕ derivatives if necessary
    if α ∈ (1, 2)
        if derivOrder >= 1
            # The ρ derivatives are found simply through the derivatives of eν
            eρ_ρ, eϕ_ρ, ez_ρ = radiationModeComps(fiber, ω, κ, m, l, ρ, 1)
            eν_ρ = eρ_ρ*ρ_unit + eϕ_ρ*ϕ_unit + ez_ρ*z_unit
            Erm_ρ = eν_ρ*exp(1im*m*ϕ)*exp(1im*κ*z)
            
            # The ϕ derivatives can be found from the mode itself
            Erm_ϕ = 1im*m*Erm + (ρ_unit'*Erm)*ϕ_unit - (ϕ_unit'*Erm)*ρ_unit
        end
        if derivOrder >= 2
            eρ_ρρ, eϕ_ρρ, ez_ρρ = radiationModeComps(fiber, ω, κ, m, l, ρ, 2)
            eν_ρρ = eρ_ρρ*ρ_unit + eϕ_ρρ*ϕ_unit + ez_ρρ*z_unit
            Erm_ρρ = eν_ρρ*exp(1im*m*ϕ)*exp(1im*κ*z)
            
            # The ϕ derivatives can be found from the mode itself
            Erm_ϕϕ = -(m^2 + 1)*Erm + (z_unit'*Erm)*z_unit + 2im*m*((ρ_unit'*Erm)*ϕ_unit - (ϕ_unit'*Erm)*ρ_unit)
            
            # The mixed double derivative is found by combining the above relations
            Erm_ρϕ = 1im*m*Erm_ρ + (ρ_unit'*Erm_ρ)*ϕ_unit - (ϕ_unit'*Erm_ρ)*ρ_unit
        end
    end
    
    # The Cartezian derivatives are calculated from the cylindrical derivatives
    if derivOrder == 1
        if     α == 1 return cos(ϕ)*Erm_ρ - sin(ϕ)/ρ*Erm_ϕ
        elseif α == 2 return sin(ϕ)*Erm_ρ + cos(ϕ)/ρ*Erm_ϕ
        elseif α == 3 return 1im*κ*Erm
        end
    end
    if derivOrder == 2
        if     α == 1 return cos(ϕ)^2*Erm_ρρ + sin(ϕ)^2/ρ^2*Erm_ϕϕ - 2*cos(ϕ)*sin(ϕ)/ρ*Erm_ρϕ + sin(ϕ)^2/ρ*Erm_ρ + 2*cos(ϕ)*sin(ϕ)/ρ^2*Erm_ϕ
        elseif α == 2 return sin(ϕ)^2*Erm_ρρ + cos(ϕ)^2/ρ^2*Erm_ϕϕ + 2*cos(ϕ)*sin(ϕ)/ρ*Erm_ρϕ + cos(ϕ)^2/ρ*Erm_ρ - 2*cos(ϕ)*sin(ϕ)/ρ^2*Erm_ϕ
        elseif α == 3 return -κ^2*Erm
        end
    end
    if derivOrder > 2 throw(ArgumentError("Erm is not implemented for derivOrder > 2")) end
end


"""
Calculates the individual components of the radiation modes or their derivatives.

Uses the normalization convention of Kien, Rauschenbeutel 2017
"""
function radiationModeComps(fiber, ω, κ, m, l, ρ, derivOrder=0)
    ρf, n = fiber.radius, fiber.refractive_index
    
    # Set up some parameters
    h = in_momentum(κ, ω, n)
    q = out_momentum(κ, ω)
    
    V = [m*ω*κ/(ρf*h^2*q^2)*(1 - n^2)*dbesselj(0, m, h*ρf)*conj(dbesselh(0, m, j, q*ρf)) for j in 1:2]
    M = [  1/h*dbesselj(1, m, h*ρf)*conj(dbesselh(0, m, j, q*ρf)) - 1/q*dbesselj(0, m, h*ρf)*conj(dbesselh(1, m, j, q*ρf)) for j in 1:2]
    L = [n^2/h*dbesselj(1, m, h*ρf)*conj(dbesselh(0, m, j, q*ρf)) - 1/q*dbesselj(0, m, h*ρf)*conj(dbesselh(1, m, j, q*ρf)) for j in 1:2]
    
    η = sqrt((abs2(V[1]) + abs2(L[1]))/(abs2(V[1]) + abs2(M[1])))
    B = 1im*l*η
    
    C = [(-1)^j      *1im*π*q^2*ρf/4*(    L[j] + B*1im*V[j]) for j in 1:2]
    D = [(-1)^(j - 1)*1im*π*q^2*ρf/4*(1im*V[j] - B    *M[j]) for j in 1:2]
    
    N = 8*π*ω/q^2*(abs2(C[1]) + abs2(D[1]))
    A = 1/sqrt(N)
    
    # Put together the components
    if ρ < ρf 
        eρ = A*1im*h^(derivOrder - 2)*(κ*      h*dbesselj(1 + derivOrder, m, h*ρ) + B*ω*1im*m/ρ*dbesselj(derivOrder, m, h*ρ))
        eϕ = A*1im*h^(derivOrder - 2)*(κ*1im*m/ρ*dbesselj(derivOrder, m, h*ρ)     - B*ω*      h*dbesselj(1 + derivOrder, m, h*ρ))
        ez = A*h^derivOrder*dbesselj(derivOrder, m, h*ρ)
        
        if derivOrder == 1
            eρ += -A*1im/h^2*B*ω*1im*m/ρ^2*dbesselj(0, m, h*ρ)
            eϕ += -A*1im/h^2*  κ*1im*m/ρ^2*dbesselj(0, m, h*ρ)
        elseif derivOrder == 2
            eρ += -A*1im/h^2*(2*h*B*ω*1im*m/ρ^2*dbesselj(1, m, h*ρ) - 2*B*ω*1im*m/ρ^3*dbesselj(0, m, h*ρ))
            eϕ += -A*1im/h^2*(2*h*  κ*1im*m/ρ^2*dbesselj(1, m, h*ρ) - 2*  κ*1im*m/ρ^3*dbesselj(0, m, h*ρ))
        elseif derivOrder != 0 throw(ArgumentError("radiationModeComps is not implemented for derivOrder > 2"))
        end
    else
        eρ = A*1im*q^(derivOrder - 2)*sum([C[j]*κ*      q*dbesselh(1 + derivOrder, m, j, q*ρ) + D[j]*ω*1im*m/ρ*dbesselh(derivOrder, m, j, q*ρ) for j in 1:2])
        eϕ = A*1im*q^(derivOrder - 2)*sum([C[j]*κ*1im*m/ρ*dbesselh(derivOrder, m, j, q*ρ)     - D[j]*ω*      q*dbesselh(1 + derivOrder, m, j, q*ρ) for j in 1:2])
        ez = A*q^derivOrder          *sum([C[j]*dbesselh(derivOrder, m, j, q*ρ) for j in 1:2])
        
        if derivOrder == 1
            eρ += -A*1im/q^2*sum([D[j]*ω*1im*m/ρ^2*dbesselh(0, m, j, q*ρ) for j in 1:2])
            eϕ += -A*1im/q^2*sum([C[j]*κ*1im*m/ρ^2*dbesselh(0, m, j, q*ρ) for j in 1:2])
        elseif derivOrder == 2
            eρ += -A*1im/q^2*sum([2*q*D[j]*ω*1im*m/ρ^2*dbesselh(1, m, j, q*ρ) - 2*D[j]*ω*1im*m/ρ^3*dbesselh(0, m, j, q*ρ) for j in 1:2])
            eϕ += -A*1im/q^2*sum([2*q*C[j]*κ*1im*m/ρ^2*dbesselh(1, m, j, q*ρ) - 2*C[j]*κ*1im*m/ρ^3*dbesselh(0, m, j, q*ρ) for j in 1:2])
        elseif derivOrder != 0 throw(ArgumentError("radiationModeComps is not implemented for derivOrder > 2"))
        end
    end
    return eρ, eϕ, ez
end


"""
Calculates the radiation mode Green's function or its derivatives 
"""
function Grm(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1, save_Im_Grm_trans=true, abstol=1e-5, approx_Grm_trans=(false, false), interpolate_Im_Grm_trans=false, interpolation_Im_Grm_trans=nothing)
    # The Green's function is calculated in terms of the contributions: the longitudinal part, the imaginary transverse part, and the real transverse part
    
    # First, we calculate the longitudinal part
    G0_lo     = G0_long(ω, r_field, r_source, derivOrder, α)
    
    # Second, the real transverse part
    Re_Grm_tr = Re_Grm_trans(fiber, ω, r_field, r_source, derivOrder, α, approx_Grm_trans[1])

    # Then, the imaginary transverse part (there is no imaginary longitudinal part)
    Im_Grm_tr = Im_Grm_trans(fiber, ω, r_field, r_source, derivOrder, α, save_Im_Grm_trans, abstol, approx_Grm_trans[2], interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    
    return G0_lo + Re_Grm_tr + 1im*Im_Grm_tr
end


"""
Small wrapper for the calculation of the imaginary part of the transverse part of radiation mode 
Green's function or its derivatives that exploits the Onsager reciprocity to simpilify calculations
"""
function Im_Grm_trans(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1, save_Im_Grm_trans=true, abstol=1e-5, approx_Im_Grm_trans=false, interpolate_Im_Grm_trans=false, interpolation_Im_Grm_trans=nothing)
    if approx_Im_Grm_trans return imag(G0(ω, r_field, r_source, derivOrder, α)) end
    
    if interpolate_Im_Grm_trans
        if r_field[3] < r_source[3]
            return transpose(interpolation_Im_Grm_trans["$(reverse(derivOrder)), $α"](r_source[3] - r_field[3]))
        else
            return interpolation_Im_Grm_trans["$derivOrder, $α"](r_field[3] - r_source[3])
        end
    end
    
    if r_field[3] < r_source[3]
        return transpose(Im_Grm_trans_(fiber, ω, r_source, r_field, reverse(derivOrder), α, save_Im_Grm_trans, abstol))
    else
        return Im_Grm_trans_(fiber, ω, r_field, r_source, derivOrder, α, save_Im_Grm_trans, abstol)
    end
end


"""
Calculates the imaginary part of the transverse part of radiation mode Green's function or its derivatives 

Uses the normalization convention of Kien, Rauschenbeutel 2017
"""
function Im_Grm_trans_(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1, save_Im_Grm_trans=true, abstol=1e-5)
    # For distances greater than 100 times the wavelength, the coupling is effectively zero (the calculation below just yields random small values)
    # We should be careful with this cut-off if we ever do calculations that truly depend on long-distance interactions
    if norm(r_field - r_source) > 100
        return zeros(ComplexF64, 3, 3)
    else
        # The Green's function only depends on the difference in the z-coordinates, 
        # so we save the calculation according to that value, rather than the individual z-coordinates
        coords = ro.([r_field[1], r_field[2], r_source[1], r_source[2], r_field[3] - r_source[3]])
        postfix = get_postfix_Im_Grm_trans(ω, coords, derivOrder, α, abstol, fiber.postfix)
        filename = "IGrmt_" * postfix
        folder = "Im_Grm_trans/"
        
        if isfile(saveDir * folder * filename * ".txt") return load_as_txt(saveDir * folder, filename) end
        
        # Set up the integrand and the integral domain
        integrand(x, args) = Erm(fiber, ω, x, args..., r_field , derivOrder[1], α)*
                             Erm(fiber, ω, x, args..., r_source, derivOrder[2], α)'
        domain = (-ω + eps(1.0), ω - eps(1.0))
        
        # Perform the combined sum and integration
        Im_Grm_trans = zeros(3, 3)
        summand_m = ones(ComplexF64, 3, 3)
        m = 0
        while maximum(abs.(summand_m)) > abstol
            summand_m .= 0
            for l in (-1, 1)
                args = (m, l)
                prob = IntegralProblem(integrand, domain, args)
                integral = Integrals.solve(prob, HCubatureJL())
                summand_m += integral.u/(4*ω)
                if m != 0
                    args = (-m, l)
                    prob = IntegralProblem(integrand, domain, args)
                    integral = Integrals.solve(prob, HCubatureJL())
                    summand_m += integral.u/(4*ω)
                end
            end
            # summand_m is always real (after adding all combinations of l and (m, -m)), even though integral.u is not
            Im_Grm_trans += real(summand_m)
            m += 1
        end
        
        if save_Im_Grm_trans save_as_txt(Im_Grm_trans, saveDir * folder, filename) end
        return Im_Grm_trans
    end
end


"""
Calculates the Fourier transform in z of the imaginary part of the transverse part of 
radiation mode Green's function or its (x and y) derivatives 

Uses the normalization convention of Kien, Rauschenbeutel 2017
"""
function Im_Grm_trans_FT(fiber, a, ω, ρ_field, ρ_source, kz, derivOrder=(0, 0), α=1, save_Im_Grm_trans=true, abstol=1e-5)
    if α == 3 throw(ArgumentError("Im_Grm_trans_FT cannot calculate a derivative in z, as it is Fourier transformed in z")) end
    if a > 1/2 throw(ArgumentError("Im_Grm_trans_FT is only implemented for subwavelength lattice spacings a < 1/2")) end
    # For distances greater than 100 times the wavelength, the coupling is effectively zero (the calculation below just yields random small values)
    # We should be careful with this cut-off if we ever do calculations that truly depend on long-distance interactions
    # For a kz outside the light cone the result is zero, assuming a subwavelength array 
    if norm(ρ_field - ρ_source) > 100 || abs(kz) > ω
        return zeros(ComplexF64, 3, 3)
    else
        coords = ro.(vcat(ρ_field, ρ_source, kz))
        postfix = get_postfix_Im_Grm_trans(ω, coords, derivOrder, α, abstol, fiber.postfix)
        filename = "IGrmtFT_" * postfix
        folder = "Im_Grm_trans/"
        
        if isfile(saveDir * folder * filename * ".txt") return load_as_txt(saveDir * folder, filename) end
        
        # Set up the integrand and the integral domain
        summand(args) = Erm(fiber, ω, kz, args..., [ρ_field..., 0] , derivOrder[1], α)*
                        Erm(fiber, ω, kz, args..., [ρ_source..., 0], derivOrder[2], α)' 
        
        # Perform the combined sum and integration
        Im_Grm_trans_FT = zeros(ComplexF64, 3, 3)
        summand_m = ones(ComplexF64, 3, 3)
        m = 0
        while maximum(abs.(summand_m)) > abstol
            summand_m .= 0
            for l in (-1, 1)
                args = (m, l)
                summand_m += summand(args)
                if m != 0
                    args = (-m, l)
                    summand_m += summand(args)
                end
            end
            Im_Grm_trans_FT += summand_m * π/(2*ω*a)
            m += 1
        end
        
        # if save_Im_Grm_trans save_as_txt(Im_Grm_trans_FT, saveDir * folder, filename) end
        return Im_Grm_trans_FT
    end
    
    
    # Im_Grm_FT = Im_Grm_trans_FT.(Ref(SP.fiber), SP.a, ωa, Ref(ρ_field), Ref(ρ_source), kz_range)
    # Γkz = 2*real.(3*π/ωa*(Ref(d').*Im_Grm_FT.*Ref(d)))
    # display(GLMakie.Screen(), plot(kz_range, Γkz))
end


"""
Calculates the real part of the transverse part of radiation mode Green's function or its derivatives 
"""
function Re_Grm_trans(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1, approx_Re_Grm_trans=false)
    if approx_Re_Grm_trans
        return real(G0_trans(ω, r_field, r_source, derivOrder, α))
    else
        # Calculate from imaginary part using Kramers-Kronig relation
        throw(ArgumentError("The non-approximate calculation of real part of the transverse part of radiation mode Green's function or its derivatives in Re_Grm_trans has not been implemented"))
    end
end


"""
Calculates the real part of the radiation mode Green's function or its derivatives 
"""
function Re_Grm(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1, approx_Re_Grm_trans=false)
    # First, we calculate the longitudinal part
    G0_lo     = G0_long(ω, r_field, r_source, derivOrder, α)
    
    # Second, the real transverse part
    Re_Grm_tr = Re_Grm_trans(fiber, ω, r_field, r_source, derivOrder, α, approx_Re_Grm_trans)
    
    return G0_lo + Re_Grm_tr
end


# ================================================
#   Functions pertaining to vacuum Green's function
# ================================================
"""
Calculates the vacuum Green's function or its derivatives 

If the function is evaluated at the origin, r_field = r_source, where its real part is not defined,
it simply returns the imaginary part (corresponding to absorbing the divergent real part in the definition of ωa)
"""
function G0(ω, r_field, r_source, derivOrder=(0, 0), α=1)
    # Set up rvec and its length, as well as the coordinate with respect to which we are deriving and its corresponding unit vector
    rvec  = r_field - r_source
    r     = norm(rvec)
    αcoor = rvec[α]
    αhat  = zeros(3); αhat[α] = 1
    
    # Return only the imaginary part when evaluated at the origin
    if r == 0
        if sum(derivOrder) == 0 return 1im*ω/(6*π)*I(3) end
        if sum(derivOrder) == 1 return zeros(ComplexF64, 3, 3) end
        if sum(derivOrder) == 2 return 1im*ω^3/(15*π)*(-I(3) + 1/2*αhat*αhat') end
    end
    
    # Prepare rr dyad and its derivatives
    rr = drhatrhat(rvec, 0, α)
    if sum(derivOrder) >= 1 drr  = drhatrhat(rvec, 1, α) end
    if sum(derivOrder) >= 2 d2rr = drhatrhat(rvec, 2, α) end
    if sum(derivOrder) >= 3 throw(ArgumentError("G0 is not implemented for sum(derivOrder) > 2")) end
    
    # Calculate derivatives (including zeroth order)
    if sum(derivOrder) == 0
        G0 =                           (2/3*dbesselsphh(0, 0, 1, ω*r) - 1/3*dbesselsphh(0, 2, 1, ω*r))*I + dbesselsphh(0, 2, 1, ω*r)*rr
    elseif sum(derivOrder) == 1
        G0 = (              ω*αcoor/r*((2/3*dbesselsphh(1, 0, 1, ω*r) - 1/3*dbesselsphh(1, 2, 1, ω*r))*I + dbesselsphh(1, 2, 1, ω*r)*rr)
                                                                                                         + dbesselsphh(0, 2, 1, ω*r)*drr )
    elseif sum(derivOrder) == 2
        G0 = (          (ω*αcoor/r)^2*((2/3*dbesselsphh(2, 0, 1, ω*r) - 1/3*dbesselsphh(2, 2, 1, ω*r))*I + dbesselsphh(2, 2, 1, ω*r)*rr)
              + ω*(1 - αcoor^2/r^2)/r*((2/3*dbesselsphh(1, 0, 1, ω*r) - 1/3*dbesselsphh(1, 2, 1, ω*r))*I + dbesselsphh(1, 2, 1, ω*r)*rr)
                                                                                             + 2*ω*αcoor/r*dbesselsphh(1, 2, 1, ω*r)*drr
                                                                                                         + dbesselsphh(0, 2, 1, ω*r)*d2rr )
    end
    
    # Return result after appending a sign acoording to how many times the derivative was taken with respect to r_source
    # and multiplying by some global constants
    return (-1)^derivOrder[2]*1im*ω/(4*π)*G0
end


"""
Calculates the longitudinal part of the vacuum Green's function or its derivatives 

If the function is evaluated at the origin, r_field = r_source, where it is not defined,
it simply returns zero (corresponding to absorbing the divergence in the definition of ωa)

Note that this part of the Green's function is purely real
"""
function G0_long(ω, r_field, r_source, derivOrder=(0, 0), α=1)
    # Set up rvec and its length, as well as the coordinate with respect to which we are deriving
    rvec  = r_field - r_source
    r     = norm(rvec)
    αcoor = rvec[α]
    
    # Return zero when evaluated at the origin
    if r == 0 return zeros(3,3) end
    
    # Prepare rr dyad and its derivatives
    rr = drhatrhat(rvec, 0, α)
    if sum(derivOrder) >= 1 drr  = drhatrhat(rvec, 1, α) end
    if sum(derivOrder) >= 2 d2rr = drhatrhat(rvec, 2, α) end
    if sum(derivOrder) >= 3 throw(ArgumentError("G0_long is not implemented for sum(derivOrder) > 2")) end
    
    # Calculate derivatives (including zeroth order)
    if sum(derivOrder) == 0
        G0_long = (3*rr - I)/r^3
    elseif sum(derivOrder) == 1
        G0_long = (3*drr/r^3 - (3*rr - I)*3*αcoor/r^5)
    elseif sum(derivOrder) == 2
        G0_long = (3*d2rr/r^3 - 3*drr*6*αcoor/r^5 - (3*rr - I)*(3/r^5 - 15*αcoor^2/r^7))
    end
    
    # Return result after appending a sign acoording to how many times the derivative was taken with respect to r_source
    # and multiplying by some global constants
    return (-1)^derivOrder[2]*G0_long/(4*π*ω^2)
end


"""
Calculates the transverse part of the vacuum Green's function or its derivatives,
which can be be used to approximate the corresponding part of radiation Green's function

If the function is evaluated at the origin, r_field = r_source, where its real part is not defined,
it simply returns the imaginary part (corresponding to absorbing the divergent real part in the definition of ωa)
"""
function G0_trans(ω, r_field, r_source, derivOrder=(0, 0), α=1)
    return G0(ω, r_field, r_source, derivOrder, α) - G0_long(ω, r_field, r_source, derivOrder, α)
end


# ================================================
#   Functions pertaining to 1D Fourier transform (along z) of the vacuum Green's function
# ================================================
"""
Calculates the 1D Fourier transform of the vacuum Green's function (along z)
"""
function G0_1DFT(ω, a, kz)
    # Get the parallel and perpendicular components of the energy shifts and decay rates
    Re_G0_1DFT_para, Re_G0_1DFT_perp = Re_G0_1DFT(ω, a, kz)
    Im_G0_1DFT_para, Im_G0_1DFT_perp = Im_G0_1DFT(ω, a, kz)
    
    # Put together the result (note the sign on the real part, due to the defintion of the energy shift)
    return Diagonal([Re_G0_1DFT_perp + 1im*Im_G0_1DFT_perp, Re_G0_1DFT_perp + 1im*Im_G0_1DFT_perp, Re_G0_1DFT_para + 1im*Im_G0_1DFT_para])
end


"""
Calculates the parallel component and the perpendicular components of the real part 
of the 1D Fourier transform of the vacuum Green's function (along z)
"""
function Re_G0_1DFT(ω, a, kz)
    Re_G0_1DFT_para = real( 
                                      ( polylogarithm(3, exp(1im*(ω + kz)*a)) + polylogarithm(3, exp(1im*(ω - kz)*a)) )
                            - 1im*ω*a*( polylogarithm(2, exp(1im*(ω + kz)*a)) + polylogarithm(2, exp(1im*(ω - kz)*a)) )
                          )
    
    Re_G0_1DFT_perp = real( 
                            -         (  polylogarithm(3, exp(1im*(ω + kz)*a)) + polylogarithm(3, exp(1im*(ω - kz)*a)))
                            + 1im*ω*a*( polylogarithm(2, exp(1im*(ω + kz)*a)) + polylogarithm(2, exp(1im*(ω - kz)*a)) )
                            + ω^2*a^2*( polylogarithm(1, exp(1im*(ω + kz)*a)) + polylogarithm(1, exp(1im*(ω - kz)*a)) )
                          )
    
    Re_G0_1DFT_para *= 3/(2*ω^3*a^3)
    Re_G0_1DFT_perp *= 3/(4*ω^3*a^3)
    return Re_G0_1DFT_para, Re_G0_1DFT_perp
end


"""
Calculates the parallel component and the perpendicular components of the imaginary part 
of the 1D Fourier transform of the vacuum Green's function (along z)
"""
function Im_G0_1DFT(ω, a, kz)
    Im_G0_1DFT_para = 0.0
    Im_G0_1DFT_perp = 0.0
    q0 = 2*π/a
    
    if abs(kz) ≤ ω        
        Im_G0_1DFT_para += 1 - kz^2/ω^2
        Im_G0_1DFT_perp += 1 + kz^2/ω^2
    end
    m = 1
    while abs(kz + q0*m) ≤ ω        
        Im_G0_1DFT_para += 1 - (kz + q0*m)^2/ω^2
        Im_G0_1DFT_perp += 1 + (kz + q0*m)^2/ω^2
        m += 1
    end
    m = -1
    while abs(kz + q0*m) ≤ ω        
        Im_G0_1DFT_para += 1 - (kz + q0*m)^2/ω^2
        Im_G0_1DFT_perp += 1 + (kz + q0*m)^2/ω^2
        m -= 1
    end
    
    Im_G0_1DFT_para *= 3*π/(4*ω*a)
    Im_G0_1DFT_perp *= 3*π/(8*ω*a)
    return Im_G0_1DFT_para, Im_G0_1DFT_perp
end


# ================================================
#   Functions pertaining to the time evolution of the atomic and phononic degrees of freedom
# ================================================
"""
Implements the equations of motion for the atomic of freedom
to first order in the driving (for the case of no phonons)
"""
function EoMs!(dσdt, σ, Δ, Δvari, tildeΩ, tildeG)
    N = length(σ)
    for i in 1:N
            dσdt[i]  = -1im*( -Δ*σ[i] - tildeΩ[i] )
        for j in 1:N
            dσdt[i] += -1im*( -(Δvari[i, j] + tildeG[i, j])*σ[j] )
        end
    end
end


"""
Implements the equations of motion for the atomic and phononic degrees of freedom
to first order in the driving and to second order in the Lamb-Dicke parameter
"""
function EoMs!(dσdt, dBαdt, σ, Bα, Δ, Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2)
    N = length(σ)
    for i in 1:N
                     dσdt[i]        = -1im*( -Δ*σ[i] - tildeΩ[i] )
        for j in 1:N
                     dσdt[i]       += -1im*( -(Δvari[i, j] + tildeG[i, j])*σ[j] )
            for α in 1:3
                     dσdt[i]       += -1im*( -tildeGα1[α][i, j]*Bα[α][i, j] - tildeGα2[α][i, j]*Bα[α][j, j] )
                    dBαdt[α][i, j]  = -1im*( -Δ*Bα[α][i, j] - tildeGα2[α][j, i]*σ[i] )
                for k in 1:N
                    dBαdt[α][i, j] += -1im*( -(Δvari[k, j] + tildeFα[α][j, k])*Bα[α][i, k] )
                end
            end
        end
        for α in 1:3
                    dBαdt[α][i, i] += -1im*( -tildeΩα[α][i] )
            for k in 1:N
                    dBαdt[α][i, i] += -1im*( -tildeGα1[α][i, k]*σ[k] )
            end
            
        end
    end
end


"""
Implements the equations of motion for the atomic of freedom
to first order in the driving (for the case of no phonons)
including a third metastable level to facilitate EIT
""" 
function EoMs!(dσgedt, dσgsdt, σge, σgs, Δ, Δc, Δvari, tildeΩ, tildeΩc, tildeG)
    N = length(σge)
    for i in 1:N
            dσgedt[i]  = -1im*( -Δ*σge[i] - tildeΩ[i] - tildeΩc[i].*σgs[i] )
            dσgsdt[i]  = -1im*( -(Δ - Δc)*σgs[i] - conj(tildeΩc[i]).*σge[i] )
        for j in 1:N
            dσgedt[i] += -1im*( -(Δvari[i, j] + tildeG[i, j])*σge[j] )
            dσgsdt[i] += -1im*( -Δvari[i, j]*σgs[j] )
        end
    end
end


"""
Implements the equations of motion for the atomic and phononic degrees of freedom
to first order in the driving and to second order in the Lamb-Dicke parameter
including a third metastable level to facilitate EIT
""" 
function EoMs!(dσgedt, dσgsdt, dBαgedt, dBαgsdt, σge, σgs, Bαge, Bαgs, 
               Δ, Δc, να, Δvari, tildeΩ, tildeΩα, tildeΩc, tildeΩcα, tildeG, tildeFα, tildeGα1, tildeGα2)
    N = length(σge)
    for i in 1:N
                     dσgedt[i]        = -1im*( -Δ*σge[i] - tildeΩ[i] - tildeΩc[i].*σgs[i] )
                     dσgsdt[i]        = -1im*( -(Δ - Δc)*σgs[i] - conj(tildeΩc[i]).*σge[i] )
        for j in 1:N
                     dσgedt[i]       += -1im*( -(Δvari[i, j] + tildeG[i, j])*σge[j] )
                     dσgsdt[i]       += -1im*( -Δvari[i, j]*σgs[j] )
            for α in 1:3
                     dσgedt[i]       += -1im*( -tildeΩcα[α][i]*Bαgs[α][i, i] - tildeGα1[α][i, j]*Bαge[α][i, j] - tildeGα2[α][i, j]*Bαge[α][j, j] )
                    dBαgedt[α][i, j]  = -1im*( -Δ*Bαge[α][i, j] - tildeΩc[i]*Bαgs[α][i, j] - tildeGα2[α][j, i]*σge[i] )
                    dBαgsdt[α][i, j]  = -1im*( -(Δ - Δc - να[α])*Bαgs[α][i, j] - conj(tildeΩc[i])*Bαge[α][i, j] )
                for k in 1:N
                    dBαgedt[α][i, j] += -1im*( -(Δvari[k, j] + tildeFα[α][j, k])*Bαge[α][i, k] )
                end
            end
        end
        for α in 1:3
                    dBαgedt[α][i, i] += -1im*( -tildeΩα[α][i] - tildeΩcα[α][i]*σgs[i] )
                     dσgsdt[i]       += -1im*( -conj(tildeΩcα[α][i])*Bαge[α][i, i] )
                    dBαgsdt[α][i, i] += -1im*( -conj(tildeΩcα[α][i])*σge[i] )
            for k in 1:N
                    dBαgedt[α][i, i] += -1im*( -tildeGα1[α][i, k]*σge[k] )
            end
        end
    end
end


"""
Wraps the EoMs to conform with the requirements of NonlinearSolve (in the case of no phonons).

args = dσdt, σ, Δ, Δvari, tildeΩ, tildeG
"""
function EoMs_wrap_noPh(dxdt, x, args, t)
    # Unpack σ from x
    unpack_σFromx!(args[2], x)
    
    # Calculate and update dσdt
    EoMs!(args...)
        
    # Pack dσdt into dxdt
    pack_σIntox!(args[1], dxdt)
    
    return dxdt
end


"""
Wraps the EoMs to conform with the requirements of NonlinearSolve.

args = dσdt, dBαdt, σ, Bα, Δ, Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2
"""
function EoMs_wrap(dxdt, x, args, t)
    # Unpack σ, Bα from x
    unpack_σBαFromx!(args[3:4]..., x)
    
    # Calculate and update dσdt, dBαdt
    EoMs!(args...)
        
    # Pack dσdt and dBαdt into dxdt
    pack_σBαIntox!(args[1:2]..., dxdt)
    
    return dxdt
end


"""
Wraps the EoMs to conform with the requirements of NonlinearSolve (in the case of no phonons).
For the case of including a third metastable level to facilitate EIT.

args = dσgedt, dσgsdt, σge, σgs, Δ, Δc, Δvari, tildeΩ, tildeΩc, tildeG
"""
function EoMs_wrap_noPh_w3l(dxdt, x, args, t)
    # Unpack σge, σgs, Bαge, and Bαgs from x
    unpack_σgeσgsFromx!(args[3:4]..., x)
    
    # Calculate and update dσdt, dBαdt
    EoMs!(args...)
        
    # Pack dσgedt, dσgsdt, dBαgedt, and dBαgsdt into dxdt
    pack_σgeσgsIntox!(args[1:2]..., dxdt)
    
    return dxdt
end


"""
Wraps the EoMs to conform with the requirements of NonlinearSolve.
For the case of including a third metastable level to facilitate EIT.

args = dσgedt, dσgsdt, dBαgedt, dBαgsdt, σge, σgs, Bαge, Bαgs, 
       Δ, Δc, να, Δvari, tildeΩ, tildeΩα, tildeΩc, tildeΩcα, tildeG, tildeFα, tildeGα1, tildeGα2
"""
function EoMs_wrap_w3l(dxdt, x, args, t)
    # Unpack σge, σgs, Bαge, and Bαgs from x
    unpack_σgeσgsBαgeBαgsFromx!(args[5:8]..., x)
    
    # Calculate and update dσdt, dBαdt
    EoMs!(args...)
        
    # Pack dσgedt, dσgsdt, dBαgedt, and dBαgsdt into dxdt
    pack_σgeσgsBαgeBαgsIntox!(args[1:4]..., dxdt)
    
    return dxdt
end


"""
Implements the analytical solution for the steady state values of the atomic 
degrees of freedom to first order in the driving (in the case of no phonons)
"""
function steadyState(Δ, Δvari, tildeΩ, tildeG)
    return -(Δ*I + Δvari + tildeG)\tildeΩ
end


"""
Implements the analytical solution for the steady state values of the atomic and phononic 
degrees of freedom to first order in the driving and to second order in the Lamb-Dicke parameter.
"""
function steadyState(Δ, Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2)
    tildeFα_inv = inv.(Ref(Δ*I + Δvari) .+ tildeFα)
    
    # Calculate the coefficient matrices
    Cα = @. (  Diag(tildeGα1*tildeFα_inv)
             + tildeGα2*Diag(tildeFα_inv) )
    Dα = @. (  Diag(tildeGα1*tildeFα_inv)*tildeGα1
             + Diag(tildeGα1*tildeFα_inv*tildeGα2)
             + tildeGα2*Diag(tildeFα_inv)*tildeGα1
             + tildeGα2*Diag(tildeFα_inv*tildeGα2) )
    
    # Finally, we calculate the steady state values
    σ_SS  = -(Δ*I + Δvari + tildeG - sum(Dα))\(tildeΩ - sum(Cα.*tildeΩα))
    Bα_SS = -(Di.(tildeΩα .+ tildeGα1.*Ref(σ_SS)) .+ Ref(Di(σ_SS)).*transpose.(tildeGα2)).*transpose.(tildeFα_inv)
    return σ_SS, Bα_SS
end


"""
Implements the analytical solution for the steady state values of the atomic 
degrees of freedom to first order in the driving (in the case of no phonons)
For the case of including a third metastable level to facilitate EIT.
"""
function steadyState(Δ, Δc, Δvari, tildeΩ, tildeΩc, tildeG)
    # Prepare shifted parameter matrix
    tildeGp = tildeG - Di(abs2.(tildeΩc))/(Δ - Δc)
    
    # Finally, we calculate the steady state values
    σge_SS  = -(Δ*I + Δvari + tildeGp)\tildeΩ
    σgs_SS  = -conj(tildeΩc).*σge_SS./(Δ - Δc .+ di(Δvari))
    return σge_SS, σgs_SS
end


"""
Implements the analytical solution for the steady state values of the atomic and phononic 
degrees of freedom to first order in the driving and to second order in the Lamb-Dicke parameter.
For the case of including a third metastable level to facilitate EIT.
"""
function steadyState(Δ, Δc, να, Δvari, tildeΩ, tildeΩα, tildeΩc, tildeΩcα, tildeG, tildeFα, tildeGα1, tildeGα2)
    # Prepare shifted parameter matrices
    tildeGp  = tildeG - Di(abs2.(tildeΩc))/(Δ - Δc) - Di(sum(@. broadcast(abs2, tildeΩcα/(Δ - Δc - να))))
    tildeFαp = tildeFα .- Ref(Di(abs2.(tildeΩc)))./(Δ - Δc .- να)
    tildeGα2_1 = tildeGα2 .- Di.(broadcast.(.*, tildeΩcα, Ref(conj(tildeΩc))))./(Δ - Δc .- να) .- Di.(broadcast.(.*, conj(tildeΩcα), Ref(tildeΩc)))/(Δ - Δc)
    tildeGα2_2 = tildeGα2 .- Di.(broadcast.(.*, conj(tildeΩcα), Ref(tildeΩc)))./(Δ - Δc .- να) .- Di.(broadcast.(.*, tildeΩcα, Ref(conj(tildeΩc))))/(Δ - Δc)
    
    tildeFαp_inv = inv.(Ref(Δ*I + Δvari) .+ tildeFαp)
    
    # Calculate the coefficient matrices
    Cαpp = @. (  Diag(tildeGα1*tildeFαp_inv)
               + tildeGα2_1*Diag(tildeFαp_inv) )
    Dαpp = @. (  Diag(tildeGα1*tildeFαp_inv)*tildeGα1
               + Diag(tildeGα1*tildeFαp_inv*tildeGα2_2)
               + tildeGα2_1*Diag(tildeFαp_inv)*tildeGα1
               + tildeGα2_1*Diag(tildeFαp_inv*tildeGα2_2) )
    
    # Finally, we calculate the steady state values
    σge_SS  = -(Δ*I + Δvari + tildeGp - sum(Dαpp))\(tildeΩ - sum(Cαpp.*tildeΩα))
    Bαge_SS = -(Di.(tildeΩα .+ tildeGα1.*Ref(σge_SS)) .+ Ref(Di(σge_SS)).*transpose.(tildeGα2_2)).*transpose.(tildeFαp_inv)
    σgs_SS  = -(conj(tildeΩc).*σge_SS + sum(@. broadcast(*, conj(tildeΩcα), di(Bαge_SS))))./(Δ - Δc .+ di(Δvari))
    Bαgs_SS = -((Δ - Δc .- να).*Ref(I) .+ Ref(Δvari)).\(Ref(Di(conj(tildeΩc))).*Bαge_SS + Di.(broadcast.(*, conj(tildeΩcα), Ref(σge_SS))))
    return σge_SS, σgs_SS, Bαge_SS, Bαgs_SS
end


"""
Calculate the time evolution of the atomic coherences (in the case of no phonons)
using the eigenmodes approach
"""
function timeEvolution_eigenmodes_noPh(t, Δ, tildeΩ, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, initialState)
    # Prepare the initial state (in the real space basis and the eigenmode basis)
    σ0 = unpack_σFromx(initialState)
    
    # Prepare the initial state and the steady state in the eigenmode basis
    tildeσ0 = eigenModesMatrix_inv*σ0
    tildeσSS = -(eigenModesMatrix_inv*tildeΩ)./(Δ .+ eigenEnergies)
    
    # Put together the final result
    return eigenModesMatrix*( tildeσSS .+ (tildeσ0 - tildeσSS).*exp.(1im*(Δ .+ eigenEnergies)*t) )
end


"""
Calculate the time evolution of the atomic coherences and the atom-phonon correlations
using the eigenmodes approach
"""
function timeEvolution_eigenmodes(t, Δ, fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, initialState)
    # Prepare the initial state (in terms of the vectorized σBα)
    σBα0 = unpack_σBαFromx(initialState)
    σBαVec0 = pack_σBαIntoσBαVec(σBα0)
    
    # Prepare the initial state and the steady state in the eigenmode basis
    tildeσBαVec0  = eigenModesMatrix_inv*σBαVec0
    tildeσBαVecSS = -(eigenModesMatrix_inv*fullDrive)./(Δ .+ eigenEnergies)
    
    # Find the state at time t
    σBαVec_t = eigenModesMatrix*( tildeσBαVecSS .+ (tildeσBαVec0 - tildeσBαVecSS).*exp.(1im*(Δ .+ eigenEnergies)*t) )
    
    # Repack and return in the form of σ, Bα
    return unpack_σBαFromσBαVec(σBαVec_t)
end


"""
Prepare the groundstate in terms of the x-vector for time-evolution
"""
function groundstate(N, noPh=false, inc3l=false)
    if noPh
        if !inc3l return empty_xVector_noPh(N)
        else      return empty_xVector_noPh_w3l(N) end
    else
        if !inc3l return empty_xVector(N)
        else      return empty_xVector_w3l(N) end
    end
end


"""
Prepare a (spatially) Gaussian distribution along z of a single excitation
in terms of the x-vector for time-evolution.

The Gaussian is centered at zc, with a width given by w,
and has a momentum (i.e. phase) given by kz.

Whether the excitation is in the e- or the s-state can be determined
by setting whichState = 'e' or whichState = 's'.
"""
function GaussianState(N, array, kz, zc, w, whichState, noPh=false, inc3l=false)
    zs = [site[3] for site in array]
    σGauss = exp.(1im*kz*zs).*exp.(-(zs .- zc).^2/(2*w^2))
    σGauss /= norm(σGauss)
    
    if whichState == "e"
        σge = σGauss
        if inc3l σgs = empty_σVector(N) end
    elseif whichState == "s"
        if inc3l == false throw(ArgumentError("GaussianState with whichState = 's' prepares an excitation in the s-state and must have inc3l=true")) end
        σge = empty_σVector(N)
        σgs = σGauss
    end
    
    if noPh
        if !inc3l return pack_σIntox(σge)
        else      return pack_σgeσgsIntox(σge, σgs) end
    else
        Bαge, Bαgs = empty_BαVector(N), empty_BαVector(N)
        if !inc3l return pack_σBαIntox(σge, Bαge)
        else      return pack_σgeσgsBαgeBαgsIntox(σge, σgs, Bαge, Bαgs) end
    end
end


"""
Prepare a (spatially) Gaussian distribution along z of a single excitation in the s-state
in terms of the x-vector for time-evolution.

The Gaussian is centered in the middle of the array (with respect to the z-axis), with a width given by w
and has a momentum (i.e. phase) given by the fiber propagation constant.
"""
function Gaussian_sState(N, array, fiber, w, noPh=false, inc3l=false)
    if inc3l == false throw(ArgumentError("Gaussian_sState prepares an excitation in the s-state and must have inc3l=true")) end
    
    zs = [site[3] for site in array]
    zc = (maximum(zs) - minimum(zs))/2
    return GaussianState(N, array, fiber.propagation_constant, zc, w, "s", noPh, inc3l)
end


# ================================================
#   Functions pertaining to transport of light through the fiber
# ================================================
"""
Calculates the transmission through the guided mode (in the case of no phonons)
"""
function transmission(σ, tildeΩ, fiber)
    return 1 + 3im*π*fiber.propagation_constant_derivative/(2*ωa^2)*sum( conj(tildeΩ).*σ )
end


"""
Calculates the transmission through the guided mode
"""
function transmission(σ, Bα, tildeΩ, tildeΩα, fiber)
    return 1 + 3im*π*fiber.propagation_constant_derivative/(2*ωa^2)*sum( conj(tildeΩ).*σ + sum([conj(tildeΩα[α]).*di(Bα[α]) for α in 1:3]) )
end


"""
Calculate the transmission of light through the guided mode, using the eigenmodes approach 
(works both with or without including phonons)
"""
function transmission_eigenmodes(Δ, drive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, κ_prime)
    driveTimesModes = drive'*eigenModesMatrix
    modesTimesDrive = eigenModesMatrix_inv*drive
    return 1 - 3im*π*κ_prime/(2*ωa^2)*sum( [driveTimesModes[i]*modesTimesDrive[i]/(Δ + eigenEnergy) for (i, eigenEnergy) in enumerate(eigenEnergies)] )
end


"""
Calculate the transmission of light through the guided mode, using the eigenmodes approach 
(works both with or without including phonons)
"""
function transmission_eigenmodes_weights_resonances(Δ_range, drive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, κ_prime)
    driveTimesModes = drive'*eigenModesMatrix
    modesTimesDrive = eigenModesMatrix_inv*drive
    weights    = [3im*π*κ_prime/(2*ωa^2)*driveTimesModes[i]*modesTimesDrive[i] for i in eachindex(eigenEnergies)]
    resonances = [weights[i]./(Δ_range .+ eigenEnergy) for (i, eigenEnergy) in enumerate(eigenEnergies)]
    return weights, resonances 
end


"""
Calculates the transmission through the guided mode (in the case of no phonons)
for the case of independent decay
"""
function transmission_indepDecay(Δ, γ_gm, γ_rm, N)
    t1 = 1 - 1im*γ_gm/(Δ + 1im*(γ_gm + γ_rm)/2)
    return t1^N
end


"""
Calculates the transmission through the guided mode (in the case of no phonons)
for the case of independent decay
"""
function transmission_indepDecay(Δ, γ_gms, γ_rms)
    ts = @. 1 - 1im*γ_gms/(Δ + 1im*(γ_gms + γ_rms)/2)
    return prod(ts)
end


"""
Calculates the emission amplitude for the H-backward mode, corresponding to the reflection
assuming the atoms to be polarized in xz plane,
assuming the driving to be in the (l=1,f=1)+(l=-1,f=1) mode (H forward)
(in the case of no phonons)
"""
function reflection(σ, tildeΩ_refl, fiber)
    return 3im*π*fiber.propagation_constant_derivative/(2*ωa^2)*sum( conj(tildeΩ_refl).*σ )
end


"""
Calculates the emission amplitude for the H-backward mode, corresponding to the reflection
assuming the atoms to be polarized in xz plane,
assuming the driving to be in the (l=1,f=1)+(l=-1,f=1) mode (H forward)
"""
function reflection(σ, Bα, tildeΩ_refl, tildeΩα_refl, fiber)
    return 3im*π*fiber.propagation_constant_derivative/(2*ωa^2)*sum( conj(tildeΩ_refl).*σ + sum([conj(tildeΩα_refl[α]).*di(Bα[α]) for α in 1:3]) )
end


"""
Calculates the emission amplitudes in each of the four guided modes (assuming a 'single-mode' fiber),
assuming the driving to be in the (l=1,f=1)+(l=-1,f=1) mode (H forward)
(in the case of no phonons)
"""
function guidedEmissions(σ, fiber, d, array)
    emissions = []
    for incField_wlf in [[(1, 1,  1), (1, -1,  1)],
                         [(1, 1, -1), (1, -1, -1)],
                         [(1, 1,  1), (1,  1,  1)],
                         [(1, 1, -1), (1,  1, -1)]]
        tildeΩ = get_tildeΩs(fiber, d, incField_wlf, array, true)
        if incField_wlf == [(1, 1,  1), (1, -1,  1)]
            push!(emissions, 1 + 3im*π*fiber.propagation_constant_derivative/(2*ωa^2)*sum( conj(tildeΩ).*σ ))
        else
            push!(emissions,     3im*π*fiber.propagation_constant_derivative/(2*ωa^2)*sum( conj(tildeΩ).*σ ))
        end
    end
    return emissions
end


"""
Calculates the emission amplitudes in each of the four guided modes (assuming a 'single-mode' fiber),
assuming the driving to be in the (l=1,f=1)+(l=-1,f=1) mode (H-forward)
"""
function guidedEmissions(σ, Bα, fiber, d, ηα, array)
    emissions = []
    for incField_wlf in [[(1, 1,  1), (1, -1,  1)],
                         [(1, 1, -1), (1, -1, -1)],
                         [(1, 1,  1), (1,  1,  1)],
                         [(1, 1, -1), (1,  1, -1)]]
        tildeΩ, tildeΩα = get_tildeΩs(fiber, d, ηα, incField_wlf, array, true)
        if incField_wlf == [(1, 1,  1), (1, -1,  1)]
            push!(emissions, 1 + 3im*π*fiber.propagation_constant_derivative/(2*ωa^2)*sum( conj(tildeΩ).*σ + sum([conj(tildeΩα[α]).*di(Bα[α]) for α in 1:3]) ))
        else
            push!(emissions,     3im*π*fiber.propagation_constant_derivative/(2*ωa^2)*sum( conj(tildeΩ).*σ + sum([conj(tildeΩα[α]).*di(Bα[α]) for α in 1:3]) ))
        end
    end
    return emissions
end


"""
Calculates the (complex) index of refraction from the transmission amplitude through the fiber 
taking the atomic array to define a medium with a certain length 
"""
function refractiveIndex(t, arrayType, N_sites, a)
    if arrayType ∈ ("1Dchain", "randomZ")
        mediumLength = a*(N_sites - 1)
    elseif arrayType == "doubleChain"
        mediumLength = a*(N_sites/2 - 1)
    else 
        throw(ArgumentError("The arrayType = $arrayType has not been implemented in refractiveIndex")) 
    end
    
    return 1 + log(t)/(1im*ωa*mediumLength)
end


"""
Calculates the group velocity from the (potentially complex) index of refraction.
The velocity is given in units of the bare velocity (c = 1) and the frequency is assumed to be ωa.
"""
function groupVelocity(n, Δ_range)
    n_real = real.(n)
    dΔ     = Δ_range[2] - Δ_range[1]
    dndΔ   = zeros(length(Δ_range))
    
    dndΔ[1] = (n_real[2] - n_real[1])/dΔ
    for i in 2:length(Δ_range)-1
        dndΔ[i] = (n_real[i+1] - n_real[i-1])/(2*dΔ)
    end
    dndΔ[end] = (n_real[end] - n_real[end-1])/dΔ
    return 1 ./(n_real +  ωa*dndΔ)
end


# ================================================
#   Functions pertaining to the radiated light
# ================================================
"""
Calculate the radiated E-field, assuming no incoming radiation field
(in the case of no phonons)
"""
function radiation_Efield(σ, Grm_rrn, d)
    d_mag = sqrt(3π/ωa^3)
    return ωa^2*d_mag*sum( Grm_rrn.*d.*σ )
end


"""
Calculate the radiated E-field, assuming no incoming radiation field
"""
function radiation_Efield(σ, Bα, tildeGrm_rrn, tildeGα2rm_rrn, d)
    d_mag = sqrt(3π/ωa^3)
    return ωa^2*d_mag*sum( tildeGrm_rrn.*d.*σ + sum([tildeGα2rm_rrn[α].*d.*di(Bα[α]) for α in 1:3]) )
end


# ================================================
#   Functions pertaining to the eigenmodes of coupling matrix
# ================================================
"""
Find the eigenvalues and -vectors of a coupling matrix (i.e. the collective energies and modes)
"""
function spectrum(G)
    F = eigen(G)
    return F.values, collect.(eachcol(F.vectors))
end


"""
Find the eigenvalues and -vectors of a coupling matrix (i.e. the collective energies and modes), 
as well as the dominant k-vector (from a discrete Fourier transform) 
pertaining to each of those eigenvectors
"""
function spectrum_dominant_ks(G, a)
    eigenEnergies, eigenModes = spectrum(G)
    eigenModes_FT = discFourierTransform.(eigenModes, a)
    dominant_ks = [ks[argmax(abs.(eigenMode_FT))] for (ks, eigenMode_FT) in eigenModes_FT]
    return eigenEnergies, eigenModes, dominant_ks
end


"""
Find the eigenvalues and the matrices used to change basis to and from the corresponding eigenbasis
"""
function spectrum_basisMatrices(G)
    eigenEnergies, eigenModes = spectrum(G)
    return eigenEnergies, vectorOfCols2Matrix(eigenModes), inv(vectorOfCols2Matrix(eigenModes))
end


"""
Find the eigenvalues and the matrices used to change basis to and from the corresponding eigenbasis, 
as well as the dominant k-vector (from a discrete Fourier transform) 
pertaining to each of the eigenvectors
"""
function spectrum_dominant_ks_basisMatrices(G, a)
    eigenEnergies, eigenModes, dominant_ks = spectrum_dominant_ks(G, a)
    return eigenEnergies, dominant_ks, vectorOfCols2Matrix(eigenModes), inv(vectorOfCols2Matrix(eigenModes))
end


"""
Calculate collective energies from eigenenergies of a coupling matrix

Shallow function, but allows for a good abstraction
"""
function collEnergies(eigenEnergies)
    return -real(eigenEnergies), 2*imag(eigenEnergies)
end


"""
Calculate split collective energies from eigenmodes energies and split coupling matrix

Shallow function, but allows for a good abstraction
"""
function splitCollEnergies(eigenModesMatrix, Δvari, tildeG_gm, tildeG_rm)
    eigenEnergies_Δv   = [eigenmode'*Δvari*eigenmode for eigenmode in eachcol(eigenModesMatrix)]
    eigenEnergies_gm   = [eigenmode'*tildeG_gm*eigenmode for eigenmode in eachcol(eigenModesMatrix)]
    eigenEnergies_rm   = [eigenmode'*tildeG_rm*eigenmode for eigenmode in eachcol(eigenModesMatrix)]
    collΔ_Δv, collΓ_Δv = collEnergies(eigenEnergies_Δv)
    collΔ_gm, collΓ_gm = collEnergies(eigenEnergies_gm)
    collΔ_rm, collΓ_rm = collEnergies(eigenEnergies_rm)
    return collΔ_Δv, collΓ_Δv, collΔ_gm, collΓ_gm, collΔ_rm, collΓ_rm
end


"""
Calculate the weight ('population') in the excitation and excitation+phonon sectors respectively,
taking a vectorized version of the system state as input
"""
function statePopulations(σBαVec, noPhonons)
    if noPhonons return 1.0, 0.0 end
    
    N = Int((sqrt(12*length(σBαVec) + 1) - 1)/6) #if length(σBαVec) = N + 3N^2, then 12*length(x) + 1 = (6N + 1)^2, and N is equal to the following
    if norm(σBαVec) != 1
        σBαVec /= norm(σBαVec)
    end
    return norm(σBαVec[1:N])^2, norm(σBαVec[N+1:end])^2
end