

#================================================
    Functions pertaining to physics of the fiber and its modes
================================================#
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
function norm_constant(ω, κ, h, q, s, ρf)
    D_in_1 = (1 - s)*(1 + (1 - s)*(κ/h)^2)*(besselj(0, h*ρf)^2 + besselj(1, h*ρf)^2)
    D_in_2 = (1 + s)*(1 + (1 + s)*(κ/h)^2)*(besselj(2, h*ρf)^2 - besselj(1, h*ρf)*besselj(3, h*ρf))
    D_in = D_in_1 + D_in_2
    
    D_out_1 = (1 - s)*(1 - (1 - s)*(κ/q)^2)*(besselk(0, q*ρf)^2 - besselk(1, q*ρf)^2)
    D_out_2 = (1 + s)*(1 - (1 + s)*(κ/q)^2)*(besselk(2, q*ρf)^2 - besselk(1, q*ρf)*besselk(3, q*ρf))
    D_out = (besselj(1, h*ρf)/besselk(1, q*ρf))^2*(D_out_1 + D_out_2)
    
    return sqrt(4*ω/(π*ρf^2*κ*(D_in + D_out)))
end


"""
Implements the fiber eigenequation with all terms are moved to left, 
meaning the propagation constant is found as the roots of this.

Signature is chosen to conform with the requirements of NonlinearSolve.

params = ω, ρf, n
"""
function fiber_equation(x, params)
    ω, ρf, n = params
    
    h = in_momentum(x, ω, n)
    q = out_momentum(x, ω)
    hρf = h*ρf
    qρf = q*ρf
    
    J0 = besselj0(hρf)
    J1 = besselj1(hρf)
    K1 = besselk1(qρf)
    K1p = -(besselk0(qρf) + besselk(2, qρf))/2
    
    A = J0/(hρf*J1)
    B = (n^2 + 1)/(2*n^2) * K1p/(qρf*K1)
    C = -1/hρf^2
    D1 = ((n^2 - 1)/(2*n^2) * K1p/(qρf*K1))^2
    D2 = x^2/(n^2*ω^2) * (1/qρf^2 + 1/hρf^2)^2
    D = sqrt(D1 + D2)
    return A + B + C + D
    
    # TODO: Make this work?
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


# TODO: move this variable explanation somewhere else? these parameters are used in many places, and most functions dont have parameters explanation
"""
Calculates the guided mode or its derivatives with
l: transverse polarization index
f: forward/backward index
r: field position
derivOrder: order of derivative
α: index indicating which derivative to take (α = 1,2,3 corresponds to x,y,z)
"""
function Egm(fiber, l, f, r, derivOrder=0, α=1)
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
    q, s, C = fiber.outside_momentum, fiber.s_parameter, fiber.normalization_constant
    
    if ρ < ρf throw(ArgumentError("guidedModeComps is only implemented for field points outside the fiber")) end
    
    # Put together the components
    eρ = 1im*C*q^derivOrder*((1 - s)*dbesselk(derivOrder, 0, q*ρ) + (1 + s)*dbesselk(derivOrder, 2, q*ρ))
    eϕ =    -C*q^derivOrder*((1 - s)*dbesselk(derivOrder, 0, q*ρ) - (1 + s)*dbesselk(derivOrder, 2, q*ρ))
    ez = C*q^derivOrder*2*q/κ*dbesselk(derivOrder, 1, q*ρ)
    
    return eρ, eϕ, ez
end


"""
Calculates the guided mode Green's function or its derivatives 
(ignoring the step-function when it comes to derivatives, 
i.e. the step-function is never differentiated)
"""
function Ggm(fiber, r_field, r_source, derivOrder=(0, 0), α=1)
    return 1im/(2*fiber.frequency)*fiber.propagation_constant_derivative*
            sum([heaviside(f*(r_field[3] - r_source[3]))*
                 Egm(fiber, l, f, r_field , derivOrder[1], α)*
                 Egm(fiber, l, f, r_source, derivOrder[2], α)' for l in (1, -1), f in (1, -1)])
end


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
"""
function radiationModeComps(fiber, ω, κ, m, l, ρ, derivOrder=0)
    ρf, n = fiber.radius, fiber.refractive_index
    
    if ρ < ρf throw(ArgumentError("radiationModeComps is only implemented for field points outside the fiber")) end
    
    # Set up some parameters
    # TODO: these quantities do not depend on position, so perhaps they should be calculated further out/up, to avoid repetition (same comment goes for other quantities in the above mode-related functions)
    h = in_momentum(κ, ω, n)
    q = out_momentum(κ, ω)
    
    V = [m*ω*κ/(ρf*h^2*q^2)*(1 - n^2)*dbesselj(0, m, h*ρf)*conj(dbesselh(0, m, j, q*ρf)) for j in 1:2]
    M = [  1/h*dbesselj(1, m, h*ρf)*conj(dbesselh(0, m, j, q*ρf)) - 1/q*dbesselj(0, m, h*ρf)*conj(dbesselh(1, m, j, q*ρf)) for j in 1:2]
    L = [n^2/h*dbesselj(1, m, h*ρf)*conj(dbesselh(0, m, j, q*ρf)) - 1/q*dbesselj(0, m, h*ρf)*conj(dbesselh(1, m, j, q*ρf)) for j in 1:2]
    
    η = sqrt((abs2(V[1]) + abs2(L[1]))/(abs2(V[1]) + n^2*abs2(M[1])))
    B = 1im*l*η
    
    C = [(-1)^j      *1im*π*q^2*ρf/4*(    L[j] + B*1im*V[j]) for j in 1:2]
    D = [(-1)^(j - 1)*1im*π*q^2*ρf/4*(1im*V[j] + B    *M[j]) for j in 1:2]
    
    N = 8π*ω/q^2*(abs2(C[1]) + abs2(D[1]))
    A = 1/sqrt(N)
    
    # Put together the components
    eρ = A*1im*q^(derivOrder - 2)*sum([C[j]*κ*      q*dbesselh(1 + derivOrder, m, j, q*ρ) + D[j]*ω*1im*m/ρ*dbesselh(derivOrder, m, j, q*ρ) for j in 1:2])
    eϕ = A*1im*q^(derivOrder - 2)*sum([C[j]*κ*1im*m/ρ*dbesselh(derivOrder, m, j, q*ρ)     - D[j]*ω      *q*dbesselh(1 + derivOrder, m, j, q*ρ) for j in 1:2])
    ez = A*q^derivOrder          *sum([C[j]*dbesselh(derivOrder, m, j, q*ρ) for j in 1:2])
    
    if derivOrder == 1
        eρ += -A*1im/q^2*sum([D[j]*ω*1im*m/ρ^2*dbesselh(0, m, j, q*ρ) for j in 1:2])
        eϕ += -A*1im/q^2*sum([C[j]*κ*1im*m/ρ^2*dbesselh(0, m, j, q*ρ) for j in 1:2])
    elseif derivOrder == 2
        eρ += -A*1im/q^2*sum([2*q*D[j]*ω*1im*m/ρ^2*dbesselh(1, m, j, q*ρ) - 2*D[j]*ω*1im*m/ρ^3*dbesselh(0, m, j, q*ρ) for j in 1:2])
        eϕ += -A*1im/q^2*sum([2*q*C[j]*κ*1im*m/ρ^2*dbesselh(1, m, j, q*ρ) - 2*C[j]*κ*1im*m/ρ^3*dbesselh(0, m, j, q*ρ) for j in 1:2])
    elseif derivOrder != 0 throw(ArgumentError("radiationModeComps is not implemented for derivOrder > 2"))
    end
    
    return eρ, eϕ, ez
end


"""
Calculates the radiation mode Green's function or its derivatives 
"""
function Grm(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1, approx_Re_Grm_trans=false)
    # The Green's function is calculated in terms of the contributions: the longitudinal part, the imaginary transverse part, and the real transverse part
    
    # First, we calculate the longitudinal part
    G0_lo = G0_long(ω, r_field, r_source, derivOrder, α)
    
    # Second, the imaginary transverse part
    Im_Grm_tr = Im_Grm_trans(fiber, ω, r_field, r_source, derivOrder, α)
    
    # Finally, the real transverse part
    Re_Grm_tr = Re_Grm_trans(fiber, ω, r_field, r_source, derivOrder, α, approx_Re_Grm_trans)
    
    return G0_lo + 1im*Im_Grm_tr + Re_Grm_tr
end


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
        if sum(derivOrder) == 2 return 1im*ω^3/(8*π)*(-I(3) + 1/2*αhat*αhat') end
    end
    
    # Prepare rr dyad and its derivatives
    rr = drhatrhat(rvec, 0, α)
    if sum(derivOrder) >= 1 drr  = drhatrhat(rvec, 1, α) end
    if sum(derivOrder) >= 2 d2rr = drhatrhat(rvec, 2, α) end
    if sum(derivOrder) >= 3 throw(ArgumentError("G0 is not implemented for sum(derivOrder) > 2")) end
    
    # Calculate derivatives (including zeroth order)
    if sum(derivOrder) == 0
        G0 = 1im*ω/4*π*((2/3*dbesselh(0, 0, 1, ω*r) - 1/3*dbesselh(0, 2, 1, ω*r))*I + dbesselh(0, 2, 1, ω*r)*rr)
    elseif sum(derivOrder) == 1
        G0 =  1im*ω/4*π*ω*αcoor/r*((2/3*dbesselh(1, 0, 1, ω*r) - 1/3*dbesselh(1, 2, 1, ω*r))*I + dbesselh(1, 2, 1, ω*r)*rr)
            + 1im*ω/4*π*(dbesselh(1, 2, 1, ω*r)*drr)
    elseif sum(derivOrder) == 2
        G0 =  1im*ω/4*π*ω*(1 - αcoor^2/r^2)/r*((2/3*dbesselh(1, 0, 1, ω*r) - 1/3*dbesselh(1, 2, 1, ω*r))*I + dbesselh(1, 2, 1, ω*r)*rr)
            + 1im*ω/4*π*(ω*αcoor/r)^2*((2/3*dbesselh(2, 0, 1, ω*r) - 1/3*dbesselh(2, 2, 1, ω*r))*I + dbesselh(2, 2, 1, ω*r)*rr)
            + 1im*ω/4*π*(dbesselh(2, 2, 1, ω*r)*drr)
            + 1im*ω/4*π*(dbesselh(1, 2, 1, ω*r)*d2rr)
    end
    
    # Return result after appending a sign acoording to how many times the derivative was taken with respect to r_source
    return (-1)^derivOrder[2]*G0
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
        G0_long = (3*rr - I)/(4*π*ω^2*r^3)
    elseif sum(derivOrder) == 1
        G0_long = (3*drr/r^3 - (3*rr - I)*3*αcoor/r^5)/(4*π*ω^2)
    elseif sum(derivOrder) == 2
        G0_long = (3*d2rr/r^3 - 3*drr*6*αcoor/r^5 - (3*rr - I)*(3/r^5 - 15*αcoor^2/r^7))/(4*π*ω^2)
    end
    
    # Return result after appending a sign acoording to how many times the derivative was taken with respect to r_source
    return (-1)^derivOrder[2]*G0_long
end


"""
Calculates the real part of the transverse part of the vacuum Green's function or its derivatives,
which can be be used to approximate the corresponding part of radiation Green's function
"""
function Re_G0_trans(ω, r_field, r_source, derivOrder=(0, 0), α=1)
    return real(G0(ω, r_field, r_source, derivOrder, α)) - G0_long(ω, r_field, r_source, derivOrder, α)
end


"""
Small wrapper for the calculation of the imaginary part of the transverse part of radiation mode 
Green's function or its derivatives that explots the Onsager reciprocity to simpilify calculations
"""
function Im_Grm_trans(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1)
    if r_field[3] < r_source[3]
        return transpose(Im_Grm_trans_calc(fiber, ω, r_source, r_field, derivOrder, α))
    else
        return Im_Grm_trans_calc(fiber, ω, r_field, r_source, derivOrder, α)
    end
end


# TODO: only calculate needed components, and exploit that some entries appear to be zero depending on derivOrder and α,
# presumably because of the simple array structure
"""
Calculates the imaginary part of the transverse part of radiation mode Green's function or its derivatives 
"""
function Im_Grm_trans_calc(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1, overwrite_bool=false)
    # The Green's function only depends on the difference in the z-coordinates, 
    # so we save the calculation according to that value, rather than the individual z-coordinates
    coords = ro.([r_field[1], r_field[2], r_source[1], r_source[2], r_field[3] - r_source[3]])
    postfix = get_postfix(ω, coords, derivOrder, α, fiber.postfix)
    filename = "IGrmt_" * postfix
    folder = "Im_Grm_trans/"
    
    if isfile(saveDir * folder * filename * ".txt")
        if overwrite_bool 
            println("The imaginary part of the radiation Green's function for \n   $filename\nhas already been calculated.\n" *
                    "Recalculating and overwriting in 5 seconds...")
            sleep(5)
        else
            # println("Loading the imaginary part of the radiation Green's function")
            return load_as_txt(saveDir * folder, filename)
        end
    end
    
    # Set up the integrand, the integral domain, and the cut-off for the m-sum
    integrand(x, args) = Erm(fiber, ω, ω*cos(x), args..., r_field , derivOrder[1], α)*
                         Erm(fiber, ω, ω*cos(x), args..., r_source, derivOrder[2], α)'
    domain = (eps(1.0), π - eps(1.0))
    abstol = 1e-3
    
    # Perform the combined sum and integration
    Im_Grm_trans = zeros(3, 3)
    summand_m = ones(ComplexF64, 3, 3)
    m = 0
    while maximum(abs.(summand_m)) > abstol
        summand_m = zeros(ComplexF64, 3, 3)
        for l in (-1, 1)
            args = (m, l)
            prob = IntegralProblem(integrand, domain, args)
            integral = Integrals.solve(prob, HCubatureJL())
            summand_m += π/2*integral
            if m != 0
                args = (-m, l)
                prob = IntegralProblem(integrand, domain, args)
                integral = Integrals.solve(prob, HCubatureJL())
                summand_m += π/2*integral
            end
        end
        # summand_m is always real (after adding all combinations of l and (m, -m)), even though integral is not
        Im_Grm_trans += real(summand_m)
        m += 1
    end
    
    save_as_txt(Im_Grm_trans, saveDir * folder, filename)
    return Im_Grm_trans
end


"""
Calculates the real part of the transverse part of radiation mode Green's function or its derivatives 
"""
function Re_Grm_trans(fiber, ω, r_field, r_source, derivOrder=(0, 0), α=1, approx_Re_Grm_trans=false)
    if approx_Re_Grm_trans
        return Re_G0_trans(ω, r_field, r_source, derivOrder, α)
    else
        # Calculate from imaginary part using Kramers-Kronig relation
        throw(ArgumentError("The non-approximate calculation of real part of the transverse part of radiation mode Green's function or its derivatives in Re_Grm_trans has not been implemented"))
    end
end


"""
Calculates the dipole moment which yields a chiral fiber setup
"""
function chiralDipoleMoment(fiber, ρa)
    eρ, _, ez = guidedModeComps(fiber, ρa) 
    return 1im*[ez, 0, -eρ]/sqrt(abs2(eρ) + ez^2)    
end


#================================================
    Functions pertaining to atomic array
================================================#
"""
Calculate a list of the atomic positions along the fiber for 1D chain array,
possibly including imperfect filling fraction and (classical) positional uncertainty
"""
function get_array(N, ρa, a, ff=1, pos_unc=0)
    # Get the perfect array (no missing atoms, no randomness in position)
    array = [[ρa, 0, a*n] for n in 0:N-1]
    # array = [SVector{3}([ρa, 0, a*n]) for n in 1:N]
    
    if      ff != 1 array = remove_atoms_from_array(array, ff) end
    if pos_unc != 0 array = introduce_position_uncertainty_to_array_sites(array, pos_unc) end
    return array
end


"""
Create imperfectly filled array according to the filling fraction ff
"""
function remove_atoms_from_array(array, ff)
    # Find total number of atoms and the number of atoms to be kept to match the desired filling fraction
    N = length(array)
    N_to_be_kept = Int(floor(N*ff))
    
    # Return N_to_be_kept of the original array sites
    return array[randperm(N)[1:N_to_be_kept]]
end


"""
Create a (classically) disordered array according to a Gaussian distribution with width pos_unc
"""
function introduce_position_uncertainty_to_array_sites(array, pos_unc)
    # Generate normally-distributed random numbers for each coordinate of each atom
    N = length(array)
    random_shift = pos_unc*randn.(fill(3, N))
    
    # Return the randomly shifted array sites
    return array + random_shift
end


#================================================
    Functions pertaining to the time evolution of the atomic and phononic degrees of freedom
================================================#
"""
Implements the equations of motion for the atomic of freedom
to first order in the driving (for the case of no phonons)
"""
function EoMs!(dσdt, σ, tildeΩ, tildeG)
    dσdt  .= -1im*(
                   -tildeΩ - tildeG*σ
                  )
end


"""
Implements the equations of motion for the atomic and phononic degrees of freedom
to first order in the driving and to second order in the Lamb-Dicke parameter
"""
function EoMs!(dσdt, dBαdt, σ, Bα, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2)
    dσdt  .= -1im*(
                   -tildeΩ - tildeG*σ - sum(@. di(Bα*transpose(tildeGα1)) + tildeGα2*di(Bα))
                  )
    dBαdt .= -1im*(
                   -Bα.*transpose.(tildeFα) - Di.(tildeΩα + tildeGα1.*Ref(σ)) - Ref(Di(σ)).*transpose.(tildeGα2)
                  )
end


"""
Wraps the EoMs to conform with the requirements of NonlinearSolve (in the case of no phonons).

args = dσdt, σ, tildeΩ, tildeG
"""
function EoMs_wrap_noPh(dxdt, x, args, t)
    # Unpack args
    dσdt, σ = args[1:2]
    
    # Unpack σ from x
    unpack_σFromx!(σ, x)
    
    # Calculate and update dσdt
    # EoMs!(dσdt, σ, args[3:end]...)
    EoMs!(args...)
        
    # Pack dσdt into dxdt
    pack_σIntox!(dσdt, dxdt)
end


"""
Wraps the EoMs to conform with the requirements of NonlinearSolve.

args = dσdt, dBαdt, σ, Bα, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2
"""
function EoMs_wrap(dxdt, x, args, t)
    # Unpack args
    dσdt, dBαdt, σ, Bα = args[1:4]
    
    # Unpack σ, Bα from x
    unpack_σBαFromx!(σ, Bα, x)
    
    # Calculate and update dσdt, dBαdt
    # EoMs!(dσdt, dBαdt, σ, Bα, args[5:end]...)
    EoMs!(args...)
        
    # Pack dσdt and dBαdt into dxdt
    pack_σBαIntox!(dσdt, dBαdt, dxdt)
end


"""
Implements the analytical solution for the steady state values of the atomic 
degrees of freedom to first order in the driving (in the case of no phonons)
"""
function σ_steadyState(tildeΩ, tildeG)
    return -tildeG\tildeΩ
end


"""
Implements the analytical solution for the steady state values of the atomic and phononic 
degrees of freedom to first order in the driving and to second order in the Lamb-Dicke parameter.
"""
function σBα_steadyState(tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2)
    # TODO: use only zeroth order tildeFα? Reduce to second order in eta in other ways (i.e. calculate for eta=0 and eta≠0 to extract exact dependencies)? Time evolution anyways also finds results that are effectively higher order
    
    tildeFα_inv = inv.(tildeFα)
    
    # Calculate the coefficient matrices
    Cα = @. (  Diag(tildeGα1*tildeFα_inv)
             + tildeGα2*Diag(tildeFα_inv) )
    Dα = @. ( Diag(tildeGα1*tildeFα_inv)*tildeGα1
             + Diag(tildeGα1*tildeFα_inv*tildeGα2)
             + tildeGα2*Diag(tildeFα_inv)*tildeGα1
             + tildeGα2*Diag(tildeFα_inv*tildeGα2) )
    
    # Finally, we calculate the steady state values
    σ_SS  = -(tildeG - sum(Dα))\(tildeΩ - sum(Cα.*tildeΩα))
    Bα_SS = -(Di.(tildeΩα + tildeGα1.*Ref(σ_SS)) + Ref(Di(σ_SS)).*transpose.(tildeGα2)).*transpose(tildeFα_inv)
    return σ_SS, Bα_SS
end


"""
Prepare the groundstate in terms of the x-vector for time-evolution

Shallow function, but allows for a good abstraction
"""
function groundstate(N, noPh=false)
    if noPh
        return empty_xVector_noPh(N)
    else
        return empty_xVector(N)
    end
end


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