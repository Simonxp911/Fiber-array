

# ================================================
#   Functions pertaining to the propagation constant
# ================================================
"""
Find the roots of the fiber eigenequation.
"""
function calc_propConst(ω, ρf, n)
    xspan = (ω + eps(ω), n*ω - eps(n*ω))
    params = (ω, ρf, n)
    
    prob = IntervalNonlinearProblem(fiber_equation, xspan, params)
    sol = NonlinearSolve.solve(prob)
    return sol.u
end


function scan_propConst(ω_range, ρf_range, n_range)
    return calc_propConst.(reshape(ω_range,  (length(ω_range), 1, 1)),
                           reshape(ρf_range, (1, length(ρf_range), 1)),
                           reshape(n_range,  (1, 1, length(n_range))))
end


"""
Calculate the derivative of roots of the fiber eigenequation with respect to frequency
"""
function calc_propConstDerivative(ω, ρf, n, dω = 1e-9)
    κ_p = calc_propConst(ω + dω/2, ρf, n)
    κ_m = calc_propConst(ω - dω/2, ρf, n)
    return (κ_p - κ_m)/dω
end


# ================================================
#   Functions pertaining to calculating the driving and couplings of the system
# ================================================
function get_parameterMatrices(SP)
    return get_parameterMatrices(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn, SP.tildeG_flags, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
end


function get_parameterMatrices(noPhonons, ΔvariDependence, Δvari_args, fiber, d, να, ηα, incField_wlf, array, ΩDriveOn, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans, include3rdLevel, cDriveType, Δc, Ωc, cDriveArgs)
    Δvari = get_Δvari(ΔvariDependence, Δvari_args, array)
    if all(ηα .== 0)
        tildeΩ = get_tildeΩs(fiber, d, incField_wlf, array, ΩDriveOn)
        tildeG = get_tildeGs(fiber, d, array, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
        if include3rdLevel tildeΩc = get_tildeΩcs(cDriveType, Ωc, array, cDriveArgs) end
        
        if !include3rdLevel return Δvari, tildeΩ, tildeG
        else                return Δc, Δvari, tildeΩ, tildeΩc, tildeG end
        
    else
        tildeΩ, tildeΩα            = get_tildeΩs(fiber, d, ηα, incField_wlf, array, ΩDriveOn)
        tildeG, tildeGα1, tildeGα2 = get_tildeGs(fiber, d, ηα, array, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
        tildeFα                    = get_tildeFα(tildeG, να)
        if include3rdLevel tildeΩc, tildeΩcα = get_tildeΩcs(cDriveType, Ωc, ηα, array, cDriveArgs) end
        
        if noPhonons
            if !include3rdLevel return Δvari, tildeΩ, tildeG
            else                return Δc, Δvari, tildeΩ, tildeΩc, tildeG end
        else
            if !include3rdLevel return Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2
            else                return Δc, να, Δvari, tildeΩ, tildeΩα, tildeΩc, tildeΩcα, tildeG, tildeFα, tildeGα1, tildeGα2 end
        end
    end
end


function get_tildeΩs(fiber, d::String, incField_wlf, array, ΩDriveOn)
    if !ΩDriveOn return zeros(ComplexF64, size(array)) end
    
    if d == "chiral"
        # Set up coordinates and guided mode components (including their first and second order derivatives)
        ρa = abs(array[1][1])
        xn = [site[1] for site in array]
        zn = [site[3] for site in array]
        κ = fiber.propagation_constant
        eρ, eϕ, ez = guidedModeComps(fiber, ρa)
        eρ = -1im*eρ #remove overall imaginary unit for ease of expressions
        dNorm = sqrt(eρ^2 + ez^2)
        propPhase = exp.(1im*κ*zn)
        ϕPhase = sign.(xn)
        
        # Put together the driving 
        Ωn = sqrt(8)*eρ*ez/dNorm*propPhase.*ϕPhase
    else
        throw(ArgumentError("Dipole moment = '$d' was not recognized in get_tildeΩs"))
    end
    return Ωn
end


function get_tildeΩs(fiber, d::String, ηα, incField_wlf, array, ΩDriveOn)
    if !ΩDriveOn return zeros(ComplexF64, size(array)), [zeros(ComplexF64, size(array)) for α in 1:3]  end
    
    if d == "chiral"
        # Set up coordinates and guided mode components (including their first and second order derivatives)
        ρa = abs(array[1][1])
        xn = [site[1] for site in array]
        zn = [site[3] for site in array]
        κ = fiber.propagation_constant
        eρ   , eϕ   , ez    = guidedModeComps(fiber, ρa)
        eρ_ρ , eϕ_ρ , ez_ρ  = guidedModeComps(fiber, ρa, 1)
        eρ_ρρ, eϕ_ρρ, ez_ρρ = guidedModeComps(fiber, ρa, 2)
        eρ, eρ_ρ, eρ_ρρ = -1im.*[eρ, eρ_ρ, eρ_ρρ] #remove overall imaginary unit for ease of expressions
        dNorm = sqrt(eρ^2 + ez^2)
        propPhase = exp.(1im*κ*zn)
        ϕPhase = sign.(xn)
        
        # Put together the driving 
        Ωn   =  sqrt(8)*eρ*ez/dNorm*propPhase.*ϕPhase
        Ωnα  = [sqrt(2)*(ez*eρ_ρ + eρ*ez_ρ)/dNorm*propPhase,
                zeros(ComplexF64, size(propPhase)),
                1im*κ*Ωn]
        Ωnαα = [sqrt(2)*(ez*eρ_ρρ + eρ*ez_ρρ)/dNorm*propPhase,
                -sqrt(2)*ez*(3*eρ + 2*eϕ)/(dNorm*ρa^2)*propPhase + sqrt(2)*(ez*eρ_ρ + eρ*ez_ρ)/(dNorm*ρa)*propPhase,
                -κ^2*Ωn]
    else
        throw(ArgumentError("Dipole moment = '$d' was not recognized in get_tildeΩs"))
    end
    return Ωn + sum(@. ηα^2*Ωnαα)/(2*ωa^2), ηα.*Ωnα/ωa
end


function get_tildeΩs(fiber, d, incField_wlf, array, ΩDriveOn)
    if !ΩDriveOn return zeros(ComplexF64, size(array)) end
    
    # Calculate incoming field
    En = fill(zeros(ComplexF64, 3), size(array))
    for (w, l, f) in incField_wlf
        En += w*Egm.(Ref(fiber), l, f, array)
    end
    En /= sqrt(sum([abs2(w) for (w, l, f) in incField_wlf], init=0.0))
    
    # Put together the driving
    return adjoint.(d).*En
end


function get_tildeΩs(fiber, d, ηα, incField_wlf, array, ΩDriveOn)
    if !ΩDriveOn return zeros(ComplexF64, size(array)), [zeros(ComplexF64, size(array)) for α in 1:3]  end
    
    # Calculate incoming field and its second derivative
    En   = fill(zeros(ComplexF64, 3), size(array))
    Enα  = [deepcopy(En) for α in 1:3]
    Enαα = deepcopy(Enα)
    for (w, l, f) in incField_wlf
        En += w*Egm.(Ref(fiber), l, f, array)
        for α in 1:3
            Enα[α]  += w*Egm.(Ref(fiber), l, f, array, 1, α)
            Enαα[α] += w*Egm.(Ref(fiber), l, f, array, 2, α)
        end
    end
    normFactor = sqrt(sum([abs2(w) for (w, l, f) in incField_wlf], init=0.0))
    En   /= normFactor
    Enα  /= normFactor
    Enαα /= normFactor
    
    # Put together the driving
    Ωn   =  adjoint.(d).*En
    Ωnα  = [adjoint.(d).*Enα[α]  for α in 1:3]
    Ωnαα = [adjoint.(d).*Enαα[α] for α in 1:3]
    return Ωn + sum(@. ηα^2*Ωnαα)/(2*ωa^2), ηα.*Ωnα/ωa
end


function get_tildeΩcs(cDriveType, Ωc, array, cDriveArgs)
    if cDriveType == "constant"
        Ωcn = fill(Ωc, size(array))
        
    elseif cDriveType == "planeWave"
        kc = cDriveArgs.kc
        Ωcn = Ωc*exp.(1im*(Ref(transpose(kc)).*array))
        
    elseif cDriveType == "hyperbolic"
        N_sites, a = cDriveArgs.N_sites, cDriveArgs.a
        zs = [site[3] for site in array]
        Ωcn = @. Ωc*sqrt(N_sites/(N_sites - zs/a))

    end
    return Ωcn
end


function get_tildeΩcs(cDriveType, Ωc, ηα, array, cDriveArgs)
    if cDriveType == "constant"
        Ωcn   = fill(Ωc, size(array))
        Ωcnα  = [zeros(size(array)) for α in 1:3]
        Ωcnαα = deepcopy(Ωcnα)
        
    elseif cDriveType == "planeWave"
        kc = cDriveArgs.kc
        Ωcn   = Ωc*exp.(1im*(Ref(transpose(kc)).*array))
        Ωcnα  = [1im*kc[α]  *Ωcn for α in 1:3]
        Ωcnαα = [   -kc[α]^2*Ωcn for α in 1:3]
        
    elseif cDriveType == "hyperbolic"
        N_sites, a = cDriveArgs.N_sites, cDriveArgs.a
        zs = [site[3] for site in array]
        Ωcn = @. Ωc*sqrt(N_sites/(N_sites - zs/a))
        Ωcnα  = [zeros(size(array)),
                 zeros(size(array)),
                 @. Ωc*sqrt(N_sites)/(2*a*(N_sites - zs/a)^(3/2))]
        Ωcnαα = [zeros(size(array)),
                 zeros(size(array)),
                 @. Ωc*3*sqrt(N_sites)/(4*a^2*(N_sites - zs/a)^(5/2))]

    end
    return Ωcn + sum(@. ηα^2*Ωcnαα)/(2*ωa^2), ηα.*Ωcnα/ωa
end


function get_tildeGs(fiber, d::String, array, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    include_Ggm, include_Grm, include_Grm_offdiag = tildeG_flags
    
    if d == "chiral"
        if save_Im_Grm_trans
            ρa = abs(array[1][1])
            d = chiralDipoleMoment(fiber, ρa, array)
            return get_tildeGs(fiber, d, array, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
        else
            # Prepare some parameters and quantities
            ρa = abs(array[1][1])
            κ = fiber.propagation_constant
            κ_prime = fiber.propagation_constant_derivative
            N = length(array)
            eρ_gm, eϕ_gm, ez_gm = guidedModeComps(fiber, ρa)
            eρ_gm = -1im*eρ_gm #remove overall imaginary unit for ease of expressions
            
            # Prepare for the integral
            domain = (-ωa + eps(1.0), ωa - eps(1.0))
            function integrand(x, args) 
                eρ_rm, eϕ_rm, ez_rm = radiationModeComps(fiber, ωa, x, args[1], args[2], ρa)
                return abs2(-1im*ez_gm*eρ_rm + eρ_gm*ez_rm)/(eρ_gm^2 + ez_gm^2) * exp(1im*x*args[3])
            end
            
            # Start calculating
            Ggm_ = zeros(ComplexF64, N, N)
            Grm_ = deepcopy(Ggm_)
            for j in 1:N, i in 1:j
                z_rel = array[i][3] - array[j][3]
                
                # The guided mode GF (exploiting Onsager reciprocity)
                if include_Ggm
                    if z_rel >= 0
                        Ggm_[i, j] += 1im/(2*ωa)*κ_prime*heaviside(z_rel)*
                                    8*eρ_gm^2*ez_gm^2/(eρ_gm^2 + ez_gm^2)*exp(1im*κ*z_rel)
                    else
                        Ggm_[j, i] += 1im/(2*ωa)*κ_prime*heaviside(-z_rel)*
                                    8*eρ_gm^2*ez_gm^2/(eρ_gm^2 + ez_gm^2)*exp(-1im*κ*z_rel)
                    end
                end
                
                # The radiation mode GF
                if include_Grm && (include_Grm_offdiag || i == j)
                    # The real part
                    Re_Grm = 0.0im
                    if approx_Grm_trans[1]
                        if z_rel != 0
                            r = abs(z_rel)
                            Re_Grm = -ωa/(4*π)*(2/3*Bessels.sphericalbessely(0, ωa*r) - 1/3*Bessels.sphericalbessely(2, ωa*r) + Bessels.sphericalbessely(2, ωa*r)*eρ_gm^2/(eρ_gm^2 + ez_gm^2))
                        end
                    else
                        throw(ArgumentError("The non-approximate calculation of real part of the transverse part of radiation mode Green's function or its derivatives in get_tildeGs has not been implemented"))
                    end
                    
                    # The imaginary (transverse) part
                    Im_Grm_tr = 0.0im
                    if approx_Grm_trans[2]
                        if z_rel == 0
                            Im_Grm_tr = ωa/(6*π)
                        else
                            r = abs(z_rel)
                            Im_Grm_tr = ωa/(4*π)*((2/3*Bessels.sphericalbesselj(0, ωa*r) - 1/3*Bessels.sphericalbesselj(2, ωa*r)) + Bessels.sphericalbesselj(2, ωa*r)*eρ_gm^2/(eρ_gm^2 + ez_gm^2))
                        end
                    else
                        # Perform the combined sum and integration
                        summand_m = 2*abstol_Im_Grm_trans + 0.0im
                        m = 0
                        while abs(summand_m) > abstol_Im_Grm_trans
                            summand_m = 0.0im
                            for l in (-1, 1)
                                args = (m, l, z_rel)
                                prob = IntegralProblem(integrand, domain, args)
                                integral = Integrals.solve(prob, HCubatureJL())
                                summand_m += integral.u/(4*ωa)
                                if m != 0
                                    args = (-m, l, z_rel)
                                    prob = IntegralProblem(integrand, domain, args)
                                    integral = Integrals.solve(prob, HCubatureJL())
                                    summand_m += integral.u/(4*ωa)
                                end
                            end
                            Im_Grm_tr += summand_m
                            m += 1
                        end
                    end
                    Grm_[i, j] = Re_Grm + 1im*Im_Grm_tr
                    Grm_[j, i] = Re_Grm + 1im*conj(Im_Grm_tr) # Onsager reciprocity
                end
            end
            
            # Get the couplings by appropriately multiplying some constants
            Ggm_ *= 3*π/ωa
            Grm_ *= 3*π/ωa
            
            # Scale the real part of the radiation GF with the local radiation decay rates (if Re_Grm_trans is being approximated)
            if approx_Grm_trans[1] && include_Grm
                gammas = 2*diag(imag(Grm_))
                scaleFactors = sqrt.(gammas*gammas')
                Grm_ = real(Grm_).*scaleFactors + 1im*imag(Grm_)
            end
            
            # Put together the full Green's function and return it
            return Ggm_ + Grm_
        end
    else
        throw(ArgumentError("get_tildeGs(d::String) is only implemented for d = 'chiral'"))
    end
end


function get_tildeGs(fiber, d::String, ηα, array, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    if d == "chiral"
        ρa = abs(array[1][1])
        d = chiralDipoleMoment(fiber, ρa, array)
        return get_tildeGs(fiber, d, ηα, array, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    else
        throw(ArgumentError("get_tildeGs(d::String) is only implemented for d = 'chiral'"))
    end
end


function get_tildeGs(fiber, d, array, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa    
    include_Ggm, include_Grm, include_Grm_offdiag = tildeG_flags
    
    N = length(array)
    Ggm_ = fill(zeros(ComplexF64, 3, 3), N, N)
    Grm_ = deepcopy(Ggm_)
    # Calculate the guided and radiation mode Green's functions
    if include_Ggm
        for j in 1:N, i in 1:j
            Ggm_[i, j] = Ggm(fiber, array[i], array[j])
            Ggm_[j, i] = transpose(Ggm_[i, j]) #Onsager reciprocity
        end
    end
    if include_Grm
        for j in 1:N, i in 1:j
            if include_Grm_offdiag || i == j
                Grm_[i, j] = Grm(fiber, ωa, array[i], array[j], (0, 0), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
                Grm_[j, i] = transpose(Grm_[i, j]) #Onsager reciprocity
            end
        end
    end
    
    # Scale the real part of the radiation GF with the local radiation decay rates (if Re_Grm_trans is being approximated)
    if approx_Grm_trans[1] && include_Grm
        gammas = 2*3*π/ωa*(adjoint.(d).*diag(imag(Grm_)).*d)
        scaleFactors = sqrt.(gammas*gammas')
        Grm_ = real(Grm_).*scaleFactors + 1im*imag(Grm_)
    end
    
    # Put together the full Green's function
    G = Ggm_ + Grm_
    
    # Get the couplings by appropriately multiplying with the dipole moment and some constants
    return 3*π/ωa*(adjoint.(d).*G.*d)
end


function get_tildeGs(fiber, d, ηα, array, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    include_Ggm, include_Grm, include_Grm_offdiag = tildeG_flags
    
    N = length(array)
    Ggm_     = fill(zeros(ComplexF64, 3, 3), N, N)
    Ggm_α1   = [deepcopy(Ggm_) for α in 1:3]
    Ggm_α2   = deepcopy(Ggm_α1)
    Ggm_αα11 = deepcopy(Ggm_α1)
    Ggm_αα22 = deepcopy(Ggm_α1)
    Ggm_αα12 = deepcopy(Ggm_α1)
    Grm_     = deepcopy(Ggm_)
    Grm_α1   = deepcopy(Ggm_α1)
    Grm_α2   = deepcopy(Ggm_α1)
    Grm_αα11 = deepcopy(Ggm_α1)
    Grm_αα22 = deepcopy(Ggm_α1)
    Grm_αα12 = deepcopy(Ggm_α1)
    
    # Calculate the guided and radiation mode Green's function and the needed derivatives (notice that Ggm_αα12 and Grm_αα12 consist of vectors)    
    if include_Ggm
        for j in 1:N, i in 1:j
            Ggm_[i, j] = Ggm(fiber, array[i], array[j])
            Ggm_[j, i] = transpose(Ggm_[i, j]) #Onsager reciprocity
        end
    end
    if include_Grm
        for j in 1:N, i in 1:j
            if include_Grm_offdiag || i == j
                Grm_[i, j] = Grm(fiber, ωa, array[i], array[j], (0, 0), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
                Grm_[j, i] = transpose(Grm_[i, j]) #Onsager reciprocity
            end
        end 
    end
    
    if include_Ggm
        for j in 1:N, i in 1:N, α in 1:3
            Ggm_α1[α][i, j]   = Ggm(fiber, array[i], array[j], (1, 0), α)
            Ggm_αα11[α][i, j] = Ggm(fiber, array[i], array[j], (2, 0), α)
            if i == j
                Ggm_αα12[α][i, j] = Ggm(fiber, array[i], array[j], (1, 1), α)
            end
        end
    end
    if include_Grm
        for j in 1:N, i in 1:N, α in 1:3
            if include_Grm_offdiag || i == j
                Grm_α1[α][i, j]   = Grm(fiber, ωa, array[i], array[j], (1, 0), α, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
                Grm_αα11[α][i, j] = Grm(fiber, ωa, array[i], array[j], (2, 0), α, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
                if i == j
                    Grm_αα12[α][i, j] = Grm(fiber, ωa, array[i], array[j], (1, 1), α, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
                end
            end
        end
    end
    
    # Onsager reciprocity
    if include_Ggm
        for j in 1:N, i in 1:N, α in 1:3
            Ggm_α2[α][i, j]   = transpose(Ggm_α1[α][j, i])
            Ggm_αα22[α][i, j] = transpose(Ggm_αα11[α][j, i])
        end
    end
    if include_Grm
        for j in 1:N, i in 1:N, α in 1:3
            if include_Grm_offdiag || i == j
                Grm_α2[α][i, j]   = transpose(Grm_α1[α][j, i])
                Grm_αα22[α][i, j] = transpose(Grm_αα11[α][j, i])
            end
        end
    end
    
    # Scale the real part of the radiation GF with the local radiation decay rates (if Re_Grm_trans is being approximated)
    if approx_Grm_trans[1] && include_Grm
        gammas = 2*imag(3*π/ωa*(adjoint.(d).*diag(Grm_).*d))
        scaleFactors = sqrt.(gammas*gammas')
        Grm_     =  real(Grm_)       .*scaleFactors + 1im*imag(Grm_)
        Grm_α1   = [real(Grm_α1[α])  .*scaleFactors + 1im*imag(Grm_α1[α])   for α in 1:3]
        Grm_α2   = [real(Grm_α2[α])  .*scaleFactors + 1im*imag(Grm_α2[α])   for α in 1:3]
        Grm_αα11 = [real(Grm_αα11[α]).*scaleFactors + 1im*imag(Grm_αα11[α]) for α in 1:3]
        Grm_αα22 = [real(Grm_αα22[α]).*scaleFactors + 1im*imag(Grm_αα22[α]) for α in 1:3]
        # Grm_αα12 is purely imaginary
    end
    
    # Put together the full Green's function
    G      = Ggm_ + Grm_
    G_α1   = Ggm_α1 + Grm_α1
    G_α2   = Ggm_α2 + Grm_α2
    G_αα11 = Ggm_αα11 + Grm_αα11
    G_αα22 = Ggm_αα22 + Grm_αα22
    G_αα12 = Ggm_αα12 + Grm_αα12
    
    # Get the couplings by appropriately multiplying with the dipole moment and some constants
    Gnm     =  3*π/ωa*(adjoint.(d).*G        .*d)
    Gnmα1   = [3*π/ωa*(adjoint.(d).*G_α1[α]  .*d) for α in 1:3]
    Gnmα2   = [3*π/ωa*(adjoint.(d).*G_α2[α]  .*d) for α in 1:3]
    Gnmαα11 = [3*π/ωa*(adjoint.(d).*G_αα11[α].*d) for α in 1:3]
    Gnmαα22 = [3*π/ωa*(adjoint.(d).*G_αα22[α].*d) for α in 1:3]
    Gnnαα12 = [3*π/ωa*(adjoint.(d).*G_αα12[α].*d) for α in 1:3]
    
    # Put together tildeG
    return Gnm + sum(@. ηα^2*(Gnmαα11 + Gnmαα22 + 2*Gnnαα12))/(2*ωa^2), ηα.*Gnmα1/ωa, ηα.*Gnmα2/ωa
end


function get_tildeGs_split(SP)
    if all(SP.ηα .== 0)
        return get_tildeGs_split(SP.fiber, SP.d, SP.array, SP.tildeG_flags[3], SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)
    else 
        return get_tildeGs_split(SP.fiber, SP.d, SP.ηα, SP.array, SP.tildeG_flags[3], SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)
    end
end


function get_tildeGs_split(fiber, d, array, include_Grm_offdiag, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    tildeG_gm = get_tildeGs(fiber, d, array, (true, false, include_Grm_offdiag), save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    tildeG_rm = get_tildeGs(fiber, d, array, (false, true, include_Grm_offdiag), save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    return tildeG_gm, tildeG_rm
end


function get_tildeGs_split(fiber, d, ηα, array, include_Grm_offdiag, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    tildeG_gm, tildeGα1_gm, tildeGα2_gm = get_tildeGs(fiber, d, ηα, array, (true, false, include_Grm_offdiag), save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    tildeG_rm, tildeGα1_rm, tildeGα2_rm = get_tildeGs(fiber, d, ηα, array, (false, true, include_Grm_offdiag), save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    return tildeG_gm, tildeGα1_gm, tildeGα2_gm, tildeG_rm, tildeGα1_rm, tildeGα2_rm
end


function get_γs(SP)
    if SP.n_inst == 1
        site = SP.array[1]
    else
        site = SP.array[1][1]
    end
    
    if SP.d == "chiral"
        ρa = abs(site[1])
        d = chiralDipoleMoment(SP.fiber, ρa, [site])[1]
    else
        if SP.n_inst == 1
            d = SP.d[1]
        else
            d = SP.d[1][1]
        end
    end
    
    return get_γs(SP.fiber, d, SP.ηα, site, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans)
end


function get_γs(fiber, d, ηα, site, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    
    if all(ηα .== 0)
        # Calculate the guided and radiation mode Green's functions
        Ggm_ = Ggm(fiber, site, site)
        Grm_ = Grm(fiber, ωa, site, site, (0, 0), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
        
        # Get the decay rates by appropriately multiplying with the dipole moment and some constants (these are always real and we want to return them as real floats)
        γ_gm = real(2*3*π/ωa*(d'*imag(Ggm_)*d))
        γ_rm = real(2*3*π/ωa*(d'*imag(Grm_)*d))
    else
        # Calculate the guided and radiation mode Green's functions
        Ggm_     =  Ggm(fiber, site, site)
        Ggm_αα11 = [Ggm(fiber, site, site, (2, 0), α) for α in 1:3]
        Ggm_αα22 = [transpose(Ggm_αα11[α]) for α in 1:3]
        Ggm_αα12 = [Ggm(fiber, site, site, (1, 1), α) for α in 1:3]
        
        Grm_     =  Grm(fiber, ωa, site, site, (0, 0), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans)
        Grm_αα11 = [Grm(fiber, ωa, site, site, (2, 0), α, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans) for α in 1:3]
        Grm_αα22 = [transpose(Grm_αα11[α]) for α in 1:3]
        Grm_αα12 = [Grm(fiber, ωa, site, site, (1, 1), α, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans) for α in 1:3]
        
        # Get the decay rates by appropriately multiplying with the dipole moment and some constants (these are always real and we want to return them as real floats)
        γ_gm = real(2*3*π/ωa*(d'*imag( Ggm_ + sum(@. ηα^2*(Ggm_αα11 + Ggm_αα22 + 2*Ggm_αα12))/(2*ωa^2) )*d))
        γ_rm = real(2*3*π/ωa*(d'*imag( Grm_ + sum(@. ηα^2*(Grm_αα11 + Grm_αα22 + 2*Grm_αα12))/(2*ωa^2) )*d))
    end
    return γ_gm, γ_rm
end


function get_tildeFα(tildeG, να)
    return [tildeG - να[α]*I for α in 1:3]
end


function get_tildeG0(fiber, d::String, array)
    if d == "chiral"
        ρa = abs(array[1][1])
        d = chiralDipoleMoment(fiber, ρa, array)
        return get_tildeG0(fiber, d, array)
    else
        throw(ArgumentError("get_tildeG0(d::String) is only implemented for d = 'chiral'"))
    end
end


function get_tildeG0(fiber, d, array)
    N = length(array)
    
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    
    G0_ = fill(zeros(ComplexF64, 3, 3), N, N)
    for j in 1:N, i in 1:j
        # Calculate vacuum Green's functions
        G0_[i, j] = G0(ωa, array[i], array[j])
        G0_[j, i] = transpose(G0_[i, j]) #Onsager reciprocity
    end
    
    # Get the couplings by appropriately multiplying with the dipole moment and some constants
    return 3*π/ωa*(adjoint.(d).*G0_.*d)
end


function get_Grm_rrns(fiber, r_field, array, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans)
    return Grm.(Ref(fiber), ωa, Ref(r_field), array, Ref((0, 0)), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, Ref(approx_Grm_trans))
end


function get_Grm_rrns(fiber, ηα, r_field, array, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans)
    Grm_rrn     =  Grm.(Ref(fiber), ωa, Ref(r_field), array, Ref((0, 0)), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, Ref(approx_Grm_trans))
    Grm_rrnα2   = [Grm.(Ref(fiber), ωa, Ref(r_field), array, Ref((0, 1)), α, save_Im_Grm_trans, abstol_Im_Grm_trans, Ref(approx_Grm_trans)) for α in 1:3]
    Grm_rrnαα22 = [Grm.(Ref(fiber), ωa, Ref(r_field), array, Ref((0, 2)), α, save_Im_Grm_trans, abstol_Im_Grm_trans, Ref(approx_Grm_trans)) for α in 1:3]
    return Grm_rrn + sum(@. ηα^2*Grm_rrnαα22)/(2*ωa^2), ηα.*Grm_rrnα2/ωa
end


function get_Jgm(fiber, d, r1, r2)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa        
    return 3*π/ωa*d'*real(Ggm(fiber, r1, r2))*d
end


function get_Γgm(fiber, d, r1, r2)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    return 2*3*π/ωa*d'*imag(Ggm(fiber, r1, r2))*d
end


function get_Jrm(fiber, d, r1, r2, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    
    # Calculate the radiation mode Green's function
    Grm_ = Grm(fiber, ωa, r1, r2, (0, 0), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans)

    # Scale the real part of the radiation GF with the local radiation decay rates (if Re_Grm_trans is being approximated)
    if approx_Grm_trans[1]
        gamma1 = 2*3*π/ωa*d'*imag(Grm(fiber, ωa, r1, r1, (0, 0), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans))*d
        gamma2 = 2*3*π/ωa*d'*imag(Grm(fiber, ωa, r2, r2, (0, 0), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans))*d
        return 3*π/ωa*d'*real(Grm_)*d * sqrt(gamma1*gamma2)
    else
        return 3*π/ωa*d'*real(Grm_)*d
    end
end


function get_Γrm(fiber, d, r1, r2, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    return 2*3*π/ωa*d'*imag(Grm(fiber, ωa, r1, r2, (0, 0), 1, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans))*d
end


function get_Δvari(ΔvariDependence, Δvari_args, array)
    if ΔvariDependence == "flat"
        return Di(zeros(size(array)))
    
    elseif ΔvariDependence == "Gaussian"
        amp, edge_width = Δvari_args
        zn = [site[3] for site in array]
        Gaussian_edges(z) = begin
            if     minimum(zn) <= z <= minimum(zn) + edge_width return (1 - exp(log(0.01)*(z - (minimum(zn) + edge_width))^2/edge_width^2))/0.99
            elseif maximum(zn) - edge_width <= z <= maximum(zn) return (1 - exp(log(0.01)*(z - (maximum(zn) - edge_width))^2/edge_width^2))/0.99
            else return 0 end
        end
        return Di(amp*Gaussian_edges.(zn))
    
    elseif ΔvariDependence == "linear"
        amp, edge_width = Δvari_args
        zn = [site[3] for site in array]
        linear_edges(z) = begin
            if     minimum(zn) <= z <= minimum(zn) + edge_width return abs((z - (minimum(zn) + edge_width))/edge_width)
            elseif maximum(zn) - edge_width <= z <= maximum(zn) return abs((z - (maximum(zn) - edge_width))/edge_width)
            else return 0 end
        end
        return Di(amp*linear_edges.(zn))
    
    elseif ΔvariDependence == "parabolic"
        amp = Δvari_args[1]
        zn = [site[3] for site in array]
        zn_shifted = zn .- (maximum(zn) + minimum(zn))/2
        parabola = zn_shifted.^2
        parabola_normalized = parabola./maximum(parabola)
        return Di(amp*parabola_normalized)
    end
end


function get_fullDriveVector(SP)
    return get_fullDriveVector(SP.noPhonons, SP.fiber, SP.d, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn, SP.include3rdLevel)
end


function get_fullDriveVector(noPhonons, fiber, d, ηα, incField_wlf, array, ΩDriveOn, include3rdLevel)
    if all(ηα .== 0)
        tildeΩ = get_tildeΩs(fiber, d, incField_wlf, array, ΩDriveOn)
    else
        tildeΩ, tildeΩα = get_tildeΩs(fiber, d, ηα, incField_wlf, array, ΩDriveOn)
    end
    N = length(tildeΩ)
    
    if include3rdLevel O1 = zeros(ComplexF64, N) end
    
    # Construct the phononic drive
    if !noPhonons
        O2 = zeros(ComplexF64, N^2)
        D = [O2 for α in 1:3]
        for α in 1:3
            for i in 1:N
                D[α][i + (i - 1)*N] = tildeΩα[α][i]
            end
        end
    end
    
    # Put together and return the full drive vector
    if noPhonons
        if !include3rdLevel
            fullDrive = tildeΩ
        else
            fullDrive = [tildeΩ
                         O1]
        end
    else
        if !include3rdLevel
            fullDrive = [tildeΩ
                         D[1]
                         D[2]
                         D[3]]
        else
            fullDrive = [tildeΩ
                         D[1]
                         D[2]
                         D[3]
                         O1
                         O2
                         O2
                         O2]
        end
    end
    return fullDrive
end


function get_fullCouplingMatrix(SP)
    return get_fullCouplingMatrix(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn, SP.tildeG_flags, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
end


function get_fullCouplingMatrix(noPhonons, ΔvariDependence, Δvari_args, fiber, d, να, ηα, incField_wlf, array, ΩDriveOn, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans, include3rdLevel, cDriveType, Δc, Ωc, cDriveArgs)
    if noPhonons
        if !include3rdLevel Δvari, tildeΩ, tildeG                                                                   = get_parameterMatrices(noPhonons, ΔvariDependence, Δvari_args, fiber, d, να, ηα, incField_wlf, array, ΩDriveOn, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans, include3rdLevel, cDriveType, Δc, Ωc, cDriveArgs)
        else                Δc, Δvari, tildeΩ, tildeΩc, tildeG                                                      = get_parameterMatrices(noPhonons, ΔvariDependence, Δvari_args, fiber, d, να, ηα, incField_wlf, array, ΩDriveOn, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans, include3rdLevel, cDriveType, Δc, Ωc, cDriveArgs) end
    else
        if !include3rdLevel Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2                             = get_parameterMatrices(noPhonons, ΔvariDependence, Δvari_args, fiber, d, να, ηα, incField_wlf, array, ΩDriveOn, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans, include3rdLevel, cDriveType, Δc, Ωc, cDriveArgs)
        else                Δc, να, Δvari, tildeΩ, tildeΩα, tildeΩc, tildeΩcα, tildeG, tildeFα, tildeGα1, tildeGα2  = get_parameterMatrices(noPhonons, ΔvariDependence, Δvari_args, fiber, d, να, ηα, incField_wlf, array, ΩDriveOn, tildeG_flags, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans, interpolate_Im_Grm_trans, interpolation_Im_Grm_trans, include3rdLevel, cDriveType, Δc, Ωc, cDriveArgs) end
    end
    N = size(Δvari)[1] #change to get this directly as an argument..?
    
    # The excitation-excitation coupling
    G = tildeG + Δvari
    
    # The two-photon detuning (with Δ set to 0, as it merely represents an irrelevant shift of the diagonal)
    if include3rdLevel
        ΔI = -Δc*I(N)
        DiΩ = Di(tildeΩc)
    end
    
    # The excitation-phonon coupling matrices (without the third level)
    if !noPhonons
        O12 = zeros(ComplexF64, N, N^2)
        O22 = zeros(ComplexF64, N^2, N^2)
        N1α = [deepcopy(O12) for α in 1:3]
        N2α = deepcopy(N1α)
        M1α = [deepcopy(O12)' for α in 1:3]
        M2α = deepcopy(M1α)
        for α in 1:3
            for j in 1:N, i in 1:N
                N1α[α][i, i + (j - 1)*N] = tildeGα1[α][i, j]
                N2α[α][i, j + (j - 1)*N] = tildeGα2[α][i, j]
                M1α[α][i + (i - 1)*N, j] = tildeGα1[α][i, j]
                M2α[α][j + (i - 1)*N, j] = tildeGα2[α][i, j]
            end
        end
        Nα = N1α + N2α
        Mα = M1α + M2α
        Fα = [kron(tildeFα[α] + Δvari, I(N)) for α in 1:3]
        
        # The full excitation-phonon coupling matrix (without the third level)
        A = [G     Nα[1] Nα[2] Nα[3]
             Mα[1] Fα[1] O22   O22
             Mα[2] O22   Fα[2] O22
             Mα[3] O22   O22   Fα[3]]
    end
    
    # The detunings and coupling matrices when including the third level (with Δ set to 0, as it merely represents an irrelevant shift of the diagonal)
    if include3rdLevel && !noPhonons
        ΔIα = [(-Δc - να[α])*I(N^2) for α in 1:3]
        IDiΩ = kron(I(N), DiΩ)
        Pα = [deepcopy(O12) for α in 1:3]
        for α in 1:3
            for i in 1:N
                Pα[α][i, i + (i - 1)*N] = tildeΩcα[α][i]
            end
        end
        Qα = transpose.(Pα)
        
        # The full coupling matrices into and within the third-level sector
        B = [DiΩ   Pα[1] Pα[2] Pα[3]
             Qα[1] IDiΩ  O22   O22
             Qα[2] O22   IDiΩ  O22
             Qα[3] O22   O22   IDiΩ]
        
        C = [ΔI   O12    O12    O12
             O12' ΔIα[1] O22    O22
             O12' O22    ΔIα[2] O22
             O12' O22    O22    ΔIα[3]]
    end
    
    
    # Put together and return the full coupling matrix
    if noPhonons
        if !include3rdLevel
            fullCoupling = G
        else
            fullCoupling = [G    DiΩ
                            DiΩ' ΔI ]
        end
    else
        if !include3rdLevel
            fullCoupling = A
        else
            fullCoupling = [A  B
                            B' C]
        end
    end
    return fullCoupling
end


function get_fullΓrm(SP)
    tildeG_flags = (false, true, SP.tildeG_flags[3])
    fullCoupling_rm = get_fullCouplingMatrix(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn, tildeG_flags, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
    return 2*imag(fullCoupling_rm)
end


function get_fullΓrm_egSector(SP)
    tildeG_flags = (false, true, SP.tildeG_flags[3])
    include3rdLevel = false
    fullCoupling_rm_egSector = get_fullCouplingMatrix(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn, tildeG_flags, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
    return 2*imag(fullCoupling_rm_egSector)
end


# ================================================
#   Functions pertaining to the steady state of the atomic and phononic degrees of freedom
# ================================================
"""
Calculate the steady state values of atomic coherences and the atom-phonon correlations (σvar)
"""
function calc_steadyState(Δ, params, postfix, save_steadyState=true)
    filename = "SS_" * postfix
    folder = "steadyStates/"
    
    if isfile(saveDir * folder * filename * ".jld2") return load_as_jld2(saveDir * folder, filename) end
    
    # Calculate steady state 
    result = steadyState(Δ, params...)
    
    if save_steadyState save_as_jld2(result, saveDir * folder, filename) end
    return result
end


"""
Calculate the steady state values of atomic coherences and the atom-phonon correlations (σvar)
for parameters given by SP and a given detuning
"""
function calc_steadyState(SP, Δ)
    postfix = get_postfix_steadyState(Δ, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.noPhonons, SP.incField_wlf, SP.tildeG_flags, SP.arrayDescription, SP.fiber.postfix, SP.include3rdLevel, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
    params = get_parameterMatrices(SP)
    return calc_steadyState(Δ, params, postfix, SP.save_steadyState)
end


"""
Scan the steady state values of atomic coherences and the atom-phonon correlations (σvar) over the detuning
"""
function scan_steadyState(SP)
    postfixes = get_postfix_steadyState.(SP.Δ_range, SP.ΔvariDescription, SP.dDescription, Ref(SP.να), Ref(SP.ηα), SP.noPhonons, Ref(SP.incField_wlf), Ref(SP.tildeG_flags), SP.arrayDescription, SP.fiber.postfix, SP.include3rdLevel, SP.cDriveDescription, SP.Δc, SP.Ωc, Ref(SP.cDriveArgs))
    params = get_parameterMatrices(SP)
    # params = nothing
    return calc_steadyState.(SP.Δ_range, Ref(params), postfixes, SP.save_steadyState)
end


"""
Scan the steady state values of atomic coherences and the atom-phonon correlations (σvar) over the detuning
for a given array (and dipole moments)
"""
function scan_steadyState(SP, d, array)
    params = get_parameterMatrices(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, d, SP.να, SP.ηα, SP.incField_wlf, array, SP.ΩDriveOn, SP.tildeG_flags, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
    return calc_steadyState.(SP.Δ_range, Ref(params), "", false)
end


# ================================================
#   Functions pertaining to the time evolution of the atomic and phononic degrees of freedom
# ================================================
"""
Perform time evolution of the atomic and phononic degrees of freedom
"""
function calc_timeEvolution(Δ, params, N, initialState, tspan, dtmax, postfix, noPhonons, include3rdLevel, whichTimeEvolver, save_timeEvol=true)
    filename = "TE_" * postfix
    folder = "timeEvol/"
    
    if isfile(saveDir * folder * filename * ".txt") return load_as_txt(saveDir * folder, filename) end
    
    # Prepare initial state and args for time evolution
    initialState_x = pack_σvarIntox(initialState, N, noPhonons, include3rdLevel)
    args = empty_σvar(N, noPhonons, include3rdLevel), empty_σvar(N, noPhonons, include3rdLevel), Δ, params, N, noPhonons, include3rdLevel
    
    # Perform time evolution
    if whichTimeEvolver == "OrdinaryDiffEq"        
        prob = ODEProblem(EoMs_wrap, initialState_x, tspan, args)
        sol = OrdinaryDiffEq.solve(prob, Tsit5(), dtmax=dtmax)
        
    elseif whichTimeEvolver == "simple"
        timeEvol_args = (tspan=tspan, dtmax=dtmax)
        sol = timeEvol(EoMs_wrap, initialState_x, args, timeEvol_args, stepCondition_endOftspan, stepFunc_nothing, "all")
    end
    
    # Pack and save data
    formattedResult = zeros(length(sol.t), 1 + length(initialState_x))
    formattedResult[:, 1] .= sol.t
    for i in eachindex(sol.t)
        formattedResult[i, 2:end] .= sol.u[i]
    end
    if save_timeEvol save_as_txt(formattedResult, saveDir * folder, filename) end
    
    # Transform data in to σvar
    σvar_t = unpack_σvarFromx.(sol.u, N, noPhonons, include3rdLevel)
        
    return sol.t, σvar_t
end


"""
Perform time evolution for parameters given by SP
"""
function calc_timeEvolution(SP, Δ)
    postfix = get_postfix_timeEvolution(Δ, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.noPhonons, SP.incField_wlf, SP.tildeG_flags, SP.arrayDescription, SP.fiber.postfix, SP.initialStateDescription, SP.tspan, SP.dtmax, SP.include3rdLevel, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
    params = get_parameterMatrices(SP)
    return calc_timeEvolution(Δ, params, SP.N, SP.initialState, SP.tspan, SP.dtmax, postfix, SP.noPhonons, SP.include3rdLevel, SP.whichTimeEvolver, SP.save_timeEvol)
end
    

"""
Scan time evolutions over the detuning
"""
function scan_timeEvolution(SP)
    postfixes = get_postfix_timeEvolution.(SP.Δ_range, SP.ΔvariDescription, SP.dDescription, Ref(SP.να), Ref(SP.ηα), SP.noPhonons, Ref(SP.incField_wlf), Ref(SP.tildeG_flags), SP.arrayDescription, SP.fiber.postfix, SP.initialStateDescription, Ref(SP.tspan), SP.dtmax, SP.include3rdLevel, SP.cDriveDescription, SP.Δc, SP.Ωc, Ref(SP.cDriveArgs))
    params = get_parameterMatrices(SP)
    return calc_timeEvolution.(SP.Δ_range, Ref(params), SP.N, Ref(SP.initialState), Ref(SP.tspan), SP.dtmax, postfixes, SP.noPhonons, SP.include3rdLevel, SP.whichTimeEvolver, SP.save_timeEvol)
end


"""
Perform time evolution of the atomic, using the eigenmodes approach
"""
function calc_timeEvolution_eigenmodes(N, Δ, fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, initialState, tspan, dtmax, postfix, noPhonons, include3rdLevel, save_timeEvol=true)
    filename = "TE_eig_" * postfix
    folder = "timeEvol/"
    
    if isfile(saveDir * folder * filename * ".txt") return load_as_txt(saveDir * folder, filename) end
    
    # Calculate time evolution
    times = range(tspan..., Int(floor((tspan[2] - tspan[1])/dtmax)))
    σvar_t = timeEvolution_eigenmodes.(times, N, Δ, Ref(fullDrive), Ref(eigenEnergies), Ref(eigenModesMatrix), Ref(eigenModesMatrix_inv), Ref(initialState), noPhonons, include3rdLevel)
    x_t = pack_σvarIntox.(σvar_t, N, noPhonons, include3rdLevel)
    
    # Pack and save data
    formattedResult = vectorOfRows2Matrix([vcat(times[i], x_t[i]) for i in eachindex(times)])
    if save_timeEvol save_as_txt(formattedResult, saveDir * folder, filename) end
    
    return times, σvar_t
end


"""
Perform time evolution for parameters given by SP, using the eigenmodes approach
"""
function calc_timeEvolution_eigenmodes(SP, Δ)
    postfix = get_postfix_timeEvolution(Δ, SP.ΔvariDescription, SP.dDescription, SP.να, SP.ηα, SP.noPhonons, SP.incField_wlf, SP.tildeG_flags, SP.arrayDescription, SP.fiber.postfix, SP.initialStateDescription, SP.tspan, SP.dtmax, SP.include3rdLevel, SP.cDriveDescription, SP.Δc, SP.Ωc, SP.cDriveArgs)
    fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv = prepare_eigenmodesCalculation(SP)
    return calc_timeEvolution_eigenmodes(SP.N, Δ, fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, SP.initialState, SP.tspan, SP.dtmax, postfix, SP.noPhonons, SP.include3rdLevel, SP.save_timeEvol)
end


"""
Scan time evolutions over the detuning, using the eigenmodes approach
"""
function scan_timeEvolution_eigenmodes(SP)
    postfixes = get_postfix_timeEvolution.(SP.Δ_range, SP.ΔvariDescription, SP.dDescription, Ref(SP.να), Ref(SP.ηα), SP.noPhonons, Ref(SP.incField_wlf), Ref(SP.tildeG_flags), SP.arrayDescription, SP.fiber.postfix, SP.initialStateDescription, Ref(SP.tspan), SP.dtmax, SP.include3rdLevel, SP.cDriveDescription, SP.Δc, SP.Ωc, Ref(SP.cDriveArgs))
    fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv = prepare_eigenmodesCalculation(SP)
    return calc_timeEvolution_eigenmodes.(SP.N, SP.Δ_range, Ref(fullDrive), Ref(eigenEnergies), Ref(eigenModesMatrix), Ref(eigenModesMatrix_inv), Ref(SP.initialState), Ref(SP.tspan), SP.dtmax, postfixes, SP.noPhonons, SP.include3rdLevel, SP.save_timeEvol)
end


# ================================================
#   Functions pertaining to transmission through the fiber
# ================================================
"""
Calculate the transmission of light through the fiber in the chosen driving mode for parameters given by SP

The function assumes that σBα contains only σ if the Lamb-Dicke parameters are zero
"""
function calc_transmission(SP, σvar)
    if SP.noPhonons
        if SP.include3rdLevel σge = σvar[1] else σge = σvar end
        tildeΩ = get_tildeΩs(SP.fiber, SP.d, SP.incField_wlf, SP.array, SP.ΩDriveOn)
        return transmission(σge, tildeΩ, SP.fiber)
    else
        if SP.include3rdLevel σge, Bαge = σvar[[1, 3]] else σge, Bαge = σvar end
        tildeΩ, tildeΩα = get_tildeΩs(SP.fiber, SP.d, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn)
        return transmission(σge, Bαge, tildeΩ, tildeΩα, SP.fiber)
    end
end


"""
Calculate the transmission of light through the fiber in the chosen driving mode for parameters given by SP
for a given array (and dipole moments)

The function assumes that σBα contains only σ if the Lamb-Dicke parameters are zero
"""
function calc_transmission(SP, σvar, d, array)
    if SP.noPhonons
        if SP.include3rdLevel σge = σvar[1] else σge = σvar end
        tildeΩ = get_tildeΩs(SP.fiber, d, SP.incField_wlf, array, SP.ΩDriveOn)
        return transmission(σge, tildeΩ, SP.fiber)
    else
        if SP.include3rdLevel σge, Bαge = σvar[[1, 3]] else σge, Bαge = σvar end
        tildeΩ, tildeΩα = get_tildeΩs(SP.fiber, d, SP.ηα, SP.incField_wlf, array, SP.ΩDriveOn)
        return transmission(σge, Bαge, tildeΩ, tildeΩα, SP.fiber)
    end
end


"""
Calculate the transmission of light through the fiber in the chosen driving mode for parameters given by SP, 
using the eigenmodes approach
"""
function calc_transmission_eigenmodes(SP, Δ)
    fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv = prepare_eigenmodesCalculation(SP)
    return transmission_eigenmodes(Δ, fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, SP.fiber.propagation_constant_derivative)
end


"""
Scan the transmission of light through the fiber in the chosen driving mode over the detuning for parameters given by SP, 
using the eigenmodes approach
"""
function scan_transmission_eigenmodes(SP)
    fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv = prepare_eigenmodesCalculation(SP)
    return transmission_eigenmodes.(SP.Δ_range, Ref(fullDrive), Ref(eigenEnergies), Ref(eigenModesMatrix), Ref(eigenModesMatrix_inv), SP.fiber.propagation_constant_derivative)
end


"""
Calculate the transmission of light through the fiber in the chosen driving mode for parameters given by SP
for the case of independent decay
"""
function calc_transmission_indepDecay(SP, Δ)
    if SP.include3rdLevel throw(ArgumentError("calc_transmission_indepDecay assumes the third level (s) is excluded")) end
    
    if SP.noPhonons
        if SP.arrayType ∈ ("1Dchain", "randomZ")
            γ_gm, γ_rm = get_γs(SP)
            return transmission_indepDecay(Δ, γ_gm, γ_rm, SP.N)
        else
            γs = [get_γs(SP.fiber, SP.d, SP.ηα, site, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans) for site in SP.array]
            γ_gms = [γ_gm for (γ_gm, γ_rm) in γs]
            γ_rms = [γ_rm for (γ_gm, γ_rm) in γs]
            return transmission_indepDecay(Δ, γ_gms, γ_rms)
        end
    else
        if SP.arrayType ∈ ("1Dchain", "randomZ")
            Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2 = get_parameterMatrices(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array[1:1], SP.ΩDriveOn, (true, true, false), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
            σvar = calc_steadyState(Δ, (Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2), "", false)
            t1 = transmission(σvar..., tildeΩ, tildeΩα, SP.fiber)
            return t1^SP.N
        else
            Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2 = get_parameterMatrices(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn, (true, true, false), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
            ts = []
            for i in 1:SP.N
                ind = i + (i-1)*SP.N
                Δvari_i    = Δvari[ind:ind]
                tildeΩ_i   = tildeΩ[ind:ind]
                tildeΩα_i  = [tildeΩα[α][ind:ind] for α in 1:3]
                tildeG_i   = tildeG[ind:ind]
                tildeFα_i  = [tildeFα[α][ind:ind] for α in 1:3]
                tildeGα1_i = [tildeGα1[α][ind:ind] for α in 1:3]
                tildeGα2_i = [tildeGα2[α][ind:ind] for α in 1:3]
                σvar = calc_steadyState(Δ, (Δvari_i, tildeΩ_i, tildeΩα_i, tildeG_i, tildeFα_i, tildeGα1_i, tildeGα2_i), "", false)
                push!(ts, transmission(σvar..., tildeΩ_i, tildeΩα_i, SP.fiber))
            end
            return prod(ts)
        end
    end
end


"""
Scan the transmission of light through the fiber in the chosen driving mode for parameters given by SP
for the case of independent decay
"""
function scan_transmission_indepDecay(SP)
    if SP.include3rdLevel throw(ArgumentError("calc_transmission_indepDecay assumes the third level (s) is excluded")) end
    
    if SP.noPhonons
        if SP.arrayType ∈ ("1Dchain", "randomZ")
            γ_gm, γ_rm = get_γs(SP)
            return transmission_indepDecay.(SP.Δ_range, γ_gm, γ_rm, SP.N)
        else
            γs = [get_γs(SP.fiber, SP.d, SP.ηα, site, SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans) for site in SP.array]
            γ_gms = [γ_gm for (γ_gm, γ_rm) in γs]
            γ_rms = [γ_rm for (γ_gm, γ_rm) in γs]
            return transmission_indepDecay.(SP.Δ_range, Ref(γ_gms), Ref(γ_rms))
        end
    else
        # Analytical expressions could be found for this and a transmission_indepDecay could be defined,
        # but most likely it would not reduce complexity much (though it would maybe read a bit cleaner)
        if SP.arrayType ∈ ("1Dchain", "randomZ")
            if typeof(SP.d) == String d = SP.d else d = SP.d[1:1] end
            Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2 = get_parameterMatrices(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, d, SP.να, SP.ηα, SP.incField_wlf, SP.array[1:1], SP.ΩDriveOn, (true, true, false), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
            σvars = calc_steadyState.(SP.Δ_range, Ref((Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2)), "", false)
            t1 = [transmission(σvar..., tildeΩ, tildeΩα, SP.fiber) for σvar in σvars]
            return t1.^SP.N
        else
            Δvari, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2 = get_parameterMatrices(SP.noPhonons, SP.ΔvariDependence, SP.Δvari_args, SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.ΩDriveOn, (true, true, false), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, SP.approx_Grm_trans, SP.interpolate_Im_Grm_trans, SP.interpolation_Im_Grm_trans, SP.include3rdLevel, SP.cDriveType, SP.Δc, SP.Ωc, SP.cDriveArgs)
            ts = []
            for i in 1:SP.N
                ind = i + (i-1)*SP.N
                Δvari_i    =  reshape(Δvari[ind:ind], 1, 1)
                tildeΩ_i   =  tildeΩ[i:i]
                tildeΩα_i  = [tildeΩα[α][i:i] for α in 1:3]
                tildeG_i   =  reshape(tildeG[ind:ind], 1, 1)
                tildeFα_i  = [reshape(tildeFα[α][ind:ind], 1, 1) for α in 1:3]
                tildeGα1_i = [reshape(tildeGα1[α][ind:ind], 1, 1) for α in 1:3]
                tildeGα2_i = [reshape(tildeGα2[α][ind:ind], 1, 1) for α in 1:3]
                σvars = calc_steadyState.(SP.Δ_range, Ref((Δvari_i, tildeΩ_i, tildeΩα_i, tildeG_i, tildeFα_i, tildeGα1_i, tildeGα2_i)), "", false)
                push!(ts, [transmission(σvar..., tildeΩ_i, tildeΩα_i, SP.fiber) for σvar in σvars])
            end
            return reduce(.*, ts)
        end
    end
end


"""
Calculate the reflection of light through the fiber, assuming dipole moments in the xz plane for parameters given by SP
"""
function calc_reflection(SP, σvar)
    if typeof(SP.d) == String
        if SP.d == "chiral"
            d = chiralDipoleMoment(SP.fiber, SP.ρa, SP.array)
        else
            throw(ArgumentError("calc_reflection is not implemented for any String dipole moments other than 'chiral'"))
        end
    else
        d = SP.d
    end
    
    if SP.noPhonons
        if SP.include3rdLevel σge = σvar[1] else σge = σvar end
        tildeΩ_refl = get_tildeΩs(SP.fiber, d, [(1, 1, -1), (1, -1, -1)], SP.array, SP.ΩDriveOn)
        return reflection(σge, tildeΩ_refl, SP.fiber)
    else
        if SP.include3rdLevel σge, Bαge = σvar[[1, 3]] else σge, Bαge = σvar end
        tildeΩ_refl, tildeΩα_refl = get_tildeΩs(SP.fiber, d, SP.ηα, [(1, 1, -1), (1, -1, -1)], SP.array, SP.ΩDriveOn)
        return reflection(σge, Bαge, tildeΩ_refl, tildeΩα_refl, SP.fiber)
    end
end


# ================================================
#   Functions pertaining to the radiation E-field around the fiber
# ================================================
"""
Calculate the intensity of the radiated light 

The function assumes that σvar contains only σge
"""
function calc_radiation_Efield(σvar, fiber, d, r_field, array, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans::Tuple=(true, true))
    Grm_rrn = get_Grm_rrns(fiber, r_field, array, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans)
    return radiation_Efield(σvar, Grm_rrn, d)
end

"""
Calculate the intensity of the radiated light

The function assumes that σvar contains only σge and Bαge
"""
function calc_radiation_Efield(σvar, fiber, d, ηα, r_field, array, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans::Tuple=(true, true))
    tildeGrm_rrn, tildeGα2rm_rrn = get_Grm_rrns(fiber, ηα, r_field, array, save_Im_Grm_trans, abstol_Im_Grm_trans, approx_Grm_trans)
    return radiation_Efield(σvar..., tildeGrm_rrn, tildeGα2rm_rrn, d)
end


"""
Scan the intensity of the radiated light for parameters given by SP over the field points given in SP

The function assumes that σvar contains only σge or only σge and Bαge
"""
function scan_radiation_Efield(SP, σvar, approx_Grm_trans::Tuple=(true, true))
    if SP.include3rdLevel throw(ArgumentError("scan_radiation_Efield assumes the third level (s) is excluded")) end
    if SP.d == "chiral" d = chiralDipoleMoment(SP.fiber, SP.ρa, SP.array) else d = SP.d end
    
    if SP.noPhonons
        return calc_radiation_Efield.(Ref(σvar), Ref(SP.fiber), Ref(d), SP.r_fields, Ref(SP.array), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, Ref(approx_Grm_trans))
    else
        return calc_radiation_Efield.(Ref(σvar), Ref(SP.fiber), Ref(d), Ref(SP.ηα), SP.r_fields, Ref(SP.array), SP.save_Im_Grm_trans, SP.abstol_Im_Grm_trans, Ref(approx_Grm_trans))
    end
end


# ================================================
#   Functions pertaining to using the atomic array as a quantum memory
# ================================================
"""
Perform time-evolution specifically for the calculation of the memory retrieval error
"""
function calc_timeEvolution_forMemoryRetrievalError(SP, fullΓrm_egSector)
    # Prepare time evolution
    initialState_x = pack_σvarIntox(SP.initialState, SP.N, SP.noPhonons, SP.include3rdLevel)
    Δ = 0 # when the g-e drive is off, the dynamics is independent of the detuning 
    params = get_parameterMatrices(SP)
    args = empty_σvar(SP.N, SP.noPhonons, SP.include3rdLevel), empty_σvar(SP.N, SP.noPhonons, SP.include3rdLevel), Δ, params, SP.N, SP.noPhonons, SP.include3rdLevel
    timeEvol_args = (tspan=SP.tspan,
                     dtmax=SP.dtmax, 
                     stepFuncValLowerTol=SP.radDecayRateAndStateNorm_LowerTol)
    
    # During time evolution we will evaluate the instantaneous decay rate as well as the state norm (time evolution will end when both are small)
    stepFunc(t, xt, ΔxΔt, timeEvol_args) = (calc_radiativeDecayRate(SP, xt, fullΓrm_egSector), norm(xt)^2)
    
    # Perform time evolution
    sol = timeEvol(EoMs_wrap, initialState_x, args, timeEvol_args, stepCondition_stepFuncVal_isSmall, stepFunc, "timeAndStepFuncVal")
    
    # Extract the instantaneous decay rates
    radiativeDecayRates = [x[1] for x in sol.stepFuncVal]
    
    return sol.t, radiativeDecayRates
end


"""
Calculate memory retrieval error
"""
function calc_memoryRetrievalError(SP)
    fullΓrm_egSector = get_fullΓrm_egSector(SP) 
    times, radiativeDecayRates = calc_timeEvolution_forMemoryRetrievalError(SP, fullΓrm_egSector)
    return memoryRetrievalError(times, radiativeDecayRates)
end


"""
Calculate the radiative decay rate, using the eigenmodes approach
"""
function calc_memoryRetrievalError_eigenmodes(SP)
    fullΓrm = get_fullΓrm(SP)
    fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv = prepare_eigenmodesCalculation(SP)
    return memoryRetrievalError_eigenmodes(eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv, SP.initialState, fullΓrm, SP.noPhonons, SP.include3rdLevel)
end


# ================================================
#   Miscellaneous
# ================================================
function prepare_eigenmodesCalculation(SP)
    fullDrive = get_fullDriveVector(SP)
    fullCoupling = get_fullCouplingMatrix(SP)
    eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv = spectrum_basisMatrices(fullCoupling)
    return fullDrive, eigenEnergies, eigenModesMatrix, eigenModesMatrix_inv
end