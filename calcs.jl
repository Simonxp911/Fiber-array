

#================================================
    Functions pertaining to the propagation constant
================================================#
"""
Find the roots of the fiber eigenequation.

Using either NonlinearSolve or HomotopyContinuation to do so 
is determined by method = "NLsolve", "HTcont"
"""
function calc_propConst(ω, ρf, n, method="NLsolve")
    xspan = (ω + eps(ω), n*ω - eps(n*ω))
    params = (ω, ρf, n)
    
    if method == "NLsolve"
        prob = IntervalNonlinearProblem(fiber_equation, xspan, params)
        sol = NonlinearSolve.solve(prob)
        return sol.u
    elseif method == "HTcont"
        throw(ArgumentError("HomotopyContinuation approach is not implemented in calc_propConst..."))
    else
        throw(ArgumentError("method=$method was not recognized in calc_propConst..."))
    end
end


function scan_propConst(SP, overwrite_bool=false)
    postfix = get_postfix(SP.ω_specs, SP.ρf_specs, SP.n_specs)
    filename = "kappa_" * postfix
    
    if isfile(saveDir * filename * ".jld2")
        if overwrite_bool 
            println("The propagation constant for \n   $filename\nhas already been calculated.\n" *
                    "Recalculating and overwriting in 5 seconds...")
            sleep(5)
        else
            println("Loading propagation constant")
            return load_as_jld2(saveDir, filename)
        end
    end
    
    κ = calc_propConst.(reshape(SP.ω_range,  (length(SP.ω_range), 1, 1)),
                        reshape(SP.ρf_range, (1, length(SP.ρf_range), 1)),
                        reshape(SP.n_range,  (1, 1, length(SP.n_range))))
    
    save_as_jld2(κ, saveDir, filename)
    return κ
end


"""
Calculate the derivative of roots of the fiber eigenequation with respect to frequency
"""
function calc_propConstDerivative(ω, ρf, n, dω = 1e-9)
    κ_p = calc_propConst(ω + dω/2, ρf, n)
    κ_m = calc_propConst(ω - dω/2, ρf, n)
    return (κ_p - κ_m)/dω
end


#================================================
    Functions pertaining to calculating the driving and couplings of the system
================================================#
function get_parameterMatrices(fiber, d, να, ηα, incField_wlf, array, approx_Re_Grm_trans)
    # TODO: implement saving and loading of these quantities?
    
    if all(ηα .== 0)
        tildeΩ = get_tildeΩs(fiber, d, incField_wlf, array)
        tildeG = get_tildeGs(fiber, d, array, approx_Re_Grm_trans)
        return tildeΩ, tildeG
    else
        tildeΩ, tildeΩα            = get_tildeΩs(fiber, d, ηα, incField_wlf, array)
        tildeG, tildeGα1, tildeGα2 = get_tildeGs(fiber, d, ηα, array, approx_Re_Grm_trans)
        tildeFα                    = get_tildeFα(tildeG, να)
        return tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2
    end
end


function get_tildeΩs(fiber, d::String, incField_wlf, array)
    if d == "chiral"
        # Set up coordinates and guided mode components (including their first and second order derivatives)
        ρa = array[1][1]
        zn = [site[3] for site in array]
        κ = fiber.propagation_constant
        eρ, eϕ, ez = guidedModeComps(fiber, ρa)
        eρ = -1im*eρ #remove overall imaginary unit for ease of expressions
        dNorm = sqrt(eρ^2 + ez^2)
        propPhase = exp.(1im*κ*zn)
        
        # Put together the driving 
        Ωn = sqrt(8)*eρ*ez/dNorm*propPhase
    else
        throw(ArgumentError("Dipole moment = '{d}' was not recognized in get_tildeΩs"))
    end
    return Ωn
end


function get_tildeΩs(fiber, d::String, ηα, incField_wlf, array)
    if d == "chiral"
        # Set up coordinates and guided mode components (including their first and second order derivatives)
        ρa = array[1][1]
        zn = [site[3] for site in array]
        κ = fiber.propagation_constant
        eρ   , eϕ   , ez    = guidedModeComps(fiber, ρa)
        eρ_ρ , eϕ_ρ , ez_ρ  = guidedModeComps(fiber, ρa, 1)
        eρ_ρρ, eϕ_ρρ, ez_ρρ = guidedModeComps(fiber, ρa, 2)
        eρ, eρ_ρ, eρ_ρρ = -1im.*[eρ, eρ_ρ, eρ_ρρ] #remove overall imaginary unit for ease of expressions
        dNorm = sqrt(eρ^2 + ez^2)
        propPhase = exp.(1im*κ*zn)
        
        # Put together the driving 
        Ωn   =  sqrt(8)*eρ*ez/dNorm*propPhase
        Ωnα  = [sqrt(2)*(ez*eρ_ρ + eρ*ez_ρ)/dNorm*propPhase,
                zeros(ComplexF64, size(propPhase)),
                1im*κ*sqrt(8)*eρ*ez/dNorm*propPhase]
        Ωnαα = [sqrt(2)*(ez*eρ_ρρ + eρ*ez_ρρ)/dNorm*propPhase,
                -sqrt(2)*ez*(3*eρ + 2*eϕ)/(dNorm*ρa^2)*propPhase + sqrt(2)*(ez*eρ_ρ + eρ*ez_ρ)/(dNorm*ρa)*propPhase,
                -κ^2*Ωn]
    else
        throw(ArgumentError("Dipole moment = '{d}' was not recognized in get_tildeΩs"))
    end
    return Ωn + sum(@. ηα^2*Ωnαα)/(2*ωa^2), ηα.*Ωnα/ωa
end


function get_tildeΩs(fiber, d, incField_wlf, array)
    # Calculate incoming field
    En = fill(zeros(ComplexF64, 3), size(array))
    for (w, l, f) in incField_wlf
        En += w*Egm.(Ref(fiber), l, f, array)
    end
    En /= sqrt(sum([abs2(w) for (w, l, f) in incField_wlf]))
    
    # Put together the driving
    Ωn = Ref(d').*En
    return Ωn 
end


function get_tildeΩs(fiber, d, ηα, incField_wlf, array)
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
    normFactor = sqrt(sum([abs2(w) for (w, l, f) in incField_wlf]))
    En /= normFactor
    for α in 1:3
        Enα[α]  /= normFactor
        Enαα[α] /= normFactor
    end
    
    # Put together the driving
    Ωn   =  Ref(d').*En
    Ωnα  = [Ref(d').*Enα[α]  for α in 1:3]
    Ωnαα = [Ref(d').*Enαα[α] for α in 1:3]
    return Ωn + sum(@. ηα^2*Ωnαα)/(2*ωa^2), ηα.*Ωnα/ωa
end


function get_tildeGs(fiber, d::String, array, approx_Re_Grm_trans)
    # TODO: Implement non-lazy version of this? Presumably significantly faster when exploiting knowledge of which components etc. are actually needed, but also very messy...
    if d == "chiral"
        ρa = array[1][1]
        if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
        d = chiralDipoleMoment(fiber, ρa)
        return get_tildeGs(fiber, d, array, approx_Re_Grm_trans)
    else
        throw(ArgumentError("get_tildeGs(d::String) is only implemented for d = 'chiral'"))
    end
end


function get_tildeGs(fiber, d::String, ηα, array, approx_Re_Grm_trans)
    # TODO: Implement non-lazy version of this? Presumably significantly faster when exploiting knowledge of which components etc. are actually needed, but also very messy...
    if d == "chiral"
        ρa = array[1][1]
        if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
        d = chiralDipoleMoment(fiber, ρa)
        return get_tildeGs(fiber, d, ηα, array, approx_Re_Grm_trans)
    else
        throw(ArgumentError("get_tildeGs(d::String) is only implemented for d = 'chiral'"))
    end
end


function get_tildeGs(fiber, d, array, approx_Re_Grm_trans)
    N = length(array)
    
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    
    Ggm_ = fill(zeros(ComplexF64, 3, 3), N, N)
    Grm_ = deepcopy(Ggm_)
    for j in 1:N, i in 1:j
        # Calculate the guided and radiation mode Green's functions
        Ggm_[i, j] = Ggm(fiber, array[i], array[j])
        Grm_[i, j] = Grm(fiber, ωa, array[i], array[j], (0, 0), 1, approx_Re_Grm_trans)
        Ggm_[j, i] = transpose(Ggm_[i, j]) #Onsager reciprocity
        Grm_[j, i] = transpose(Grm_[i, j])
    end
    
    # Scale the real part of the radiation GF with the local radiation decay rates (if Re_Grm_trans is being approximated)
    if approx_Re_Grm_trans
        gammas = 2*3*π/ωa*(Ref(d').*diag(imag(Grm_)).*Ref(d))
        scaleFactors = sqrt.(gammas*gammas')
        Grm_ = real(Grm_).*scaleFactors + 1im*imag(Grm_)
    end
    
    # Put together the full Green's function
    G = Ggm_ + Grm_
    
    # Get the couplings by appropriately multiplying with the dipole moment and some constants
    Gnm = 3*π/ωa*(Ref(d').*G.*Ref(d))
    
    # Put together tildeG
    return Gnm
end


function get_tildeGs(fiber, d, ηα, array, approx_Re_Grm_trans)
    N = length(array)
    
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    
    Ggm_     = fill(zeros(ComplexF64, 3, 3), N, N)
    Ggm_α1   = [deepcopy(Ggm_) for α in 1:3]
    Ggm_α2   = deepcopy(Ggm_α1)
    Ggm_αα11 = deepcopy(Ggm_α1)
    Ggm_αα22 = deepcopy(Ggm_α1)
    Ggm_αα12 = [fill(zeros(ComplexF64, 3, 3), N) for α in 1:3]
    Grm_     = deepcopy(Ggm_)
    Grm_α1   = deepcopy(Ggm_α1)
    Grm_α2   = deepcopy(Ggm_α1)
    Grm_αα11 = deepcopy(Ggm_α1)
    Grm_αα22 = deepcopy(Ggm_α1)
    Grm_αα12 = deepcopy(Ggm_αα12)
    
    # Calculate the guided and radiation mode Green's function and the needed derivatives (notice that Ggm_αα12 and Grm_αα12 consist of vectors)    
    for j in 1:N, i in 1:j
        Ggm_[i, j] = Ggm(fiber, array[i], array[j])
        Grm_[i, j] = Grm(fiber, ωa, array[i], array[j], (0, 0), 1, approx_Re_Grm_trans)
        Ggm_[j, i] = transpose(Ggm_[i, j]) #Onsager reciprocity
        Grm_[j, i] = transpose(Grm_[i, j])
    end
    
    for j in 1:N, i in 1:N, α in 1:3
        Ggm_α1[α][i, j]   = Ggm(fiber, array[i], array[j], (1, 0), α)
        Ggm_αα11[α][i, j] = Ggm(fiber, array[i], array[j], (2, 0), α)
        
        Grm_α1[α][i, j]   = Grm(fiber, ωa, array[i], array[j], (1, 0), α, approx_Re_Grm_trans)
        Grm_αα11[α][i, j] = Grm(fiber, ωa, array[i], array[j], (2, 0), α, approx_Re_Grm_trans)            
        if i == j
            Ggm_αα12[α][i] = Ggm(fiber, array[i], array[i], (1, 1), α)
            Grm_αα12[α][i] = Grm(fiber, ωa, array[i], array[i], (1, 1), α, approx_Re_Grm_trans)
        end
    end
    
    # Onsager reciprocity
    for j in 1:N, i in 1:N, α in 1:3
        Ggm_α2[α][i, j]   = transpose(Ggm_α1[α][j, i])
        Ggm_αα22[α][i, j] = transpose(Ggm_αα11[α][j, i])
        Grm_α2[α][i, j]   = transpose(Grm_α1[α][j, i])
        Grm_αα22[α][i, j] = transpose(Grm_αα11[α][j, i])
    end
    
    # Scale the real part of the radiation GF with the local radiation decay rates (if Re_Grm_trans is being approximated)
    if approx_Re_Grm_trans
        gammas = 2*imag(3*π/ωa*(Ref(d').*diag(Grm_).*Ref(d)))
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
    Gnm     =     3*π/ωa*(Ref(d').*G        .*Ref(d))
    Gnmα1   = [   3*π/ωa*(Ref(d').*G_α1[α]  .*Ref(d))  for α in 1:3]
    Gnmα2   = [   3*π/ωa*(Ref(d').*G_α2[α]  .*Ref(d))  for α in 1:3]
    Gnmαα11 = [   3*π/ωa*(Ref(d').*G_αα11[α].*Ref(d))  for α in 1:3]
    Gnmαα22 = [   3*π/ωa*(Ref(d').*G_αα22[α].*Ref(d))  for α in 1:3]
    Gnnαα12 = [Di(3*π/ωa*(Ref(d').*G_αα12[α].*Ref(d))) for α in 1:3]
    
    # Put together tildeG
    return Gnm + sum(@. ηα^2*(Gnmαα11 + Gnmαα22 + 2*Gnnαα12))/(2*ωa^2), ηα.*Gnmα1/ωa, ηα.*Gnmα2/ωa
end


function get_tildeFα(tildeG, να)
    return [tildeG - να[α]*I for α in 1:3]
end


function get_Jgm(fiber, d, r1, r2)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa        
    return 3*π/ωa*d'*real(Ggm(fiber, r1, r2))*d
end


function get_Γgm(fiber, d, r1, r2)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    return 2*3*π/ωa*d'*imag(Ggm(fiber, r1, r2))*d
end


function get_Jrm(fiber, d, r1, r2, approx_Re_Grm_trans)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    
    # Calculate the radiation mode Green's function
    Grm_ = Grm(fiber, ωa, r1, r2, (0, 0), 1, approx_Re_Grm_trans)

    # Scale the real part of the radiation GF with the local radiation decay rates (if Re_Grm_trans is being approximated)
    if approx_Re_Grm_trans
        gamma1 = 2*3*π/ωa*d'*imag(Grm(fiber, ωa, r1, r1, (0, 0), 1, approx_Re_Grm_trans))*d
        gamma2 = 2*3*π/ωa*d'*imag(Grm(fiber, ωa, r2, r2, (0, 0), 1, approx_Re_Grm_trans))*d
        return 3*π/ωa*d'*real(Grm_)*d * sqrt(gamma1*gamma2)
    else
        return 3*π/ωa*d'*real(Grm_)*d
    end
end


function get_Γrm(fiber, d, r1, r2, approx_Re_Grm_trans)
    if fiber.frequency != ωa fiber = Fiber(fiber.radius, fiber.refractive_index, ωa) end #atoms always interact at frequency ω = ωa
    return 2*3*π/ωa*d'*imag(Grm(fiber, ωa, r1, r2, (0, 0), 1, approx_Re_Grm_trans))*d
end


#================================================
    Functions pertaining to the time evolution of the atomic and phononic degrees of freedom
================================================#
"""
Calculate the steady state values of atomic coherences σ and the atom-phonon correlations Bα
"""
function calc_σBα_steadyState(Δ, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2, postfix)
    filename = "sigmaBalpha_" * postfix
    folder = "steadyStates/"
    
    if isfile(saveDir * folder * filename * ".jld2") load_as_jld2(saveDir * folder, filename) end
    
    # Calculate steady state σBα
    result = σBα_steadyState(Δ, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2)
    
    save_as_jld2(result, saveDir * folder, filename)
    return result
end


"""
Calculate the steady state values of atomic coherences σ and the atom-phonon correlations Bα
(for the case of no phonons)
"""
function calc_σ_steadyState(Δ, tildeΩ, tildeG, postfix)
    filename = "sigma_" * postfix
    folder = "steadyStates/"
    
    if isfile(saveDir * folder * filename * ".jld2") load_as_jld2(saveDir * folder, filename) end
    
    # Calculate steady state σ
    result = σ_steadyState(Δ, tildeΩ, tildeG)
    
    save_as_jld2(result, saveDir * folder, filename)
    return result
end


"""
Calculate the steady state values of atomic coherences σ and the atom-phonon correlations Bα 
for parameters given by SP and a given detuning
"""
function calc_σBα_steadyState(SP, Δ)
    postfix = get_postfix(Δ, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.arrayDescription, SP.fiber.postfix)
    params = get_parameterMatrices(SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.approx_Re_Grm_trans)
    if all(SP.ηα .== 0) return calc_σ_steadyState(Δ, params..., postfix)
    else                return calc_σBα_steadyState(Δ, params..., postfix)
    end
end


"""
Scan the steady state values of atomic coherences σ and the atom-phonon correlations Bα over the detuning
"""
function scan_σBα_steadyState(SP)
    postfixes = get_postfix.(SP.Δ_range, Ref(SP.d), Ref(SP.να), Ref(SP.ηα), Ref(SP.incField_wlf), SP.arrayDescription, SP.fiber.postfix)
    params = get_parameterMatrices(SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.approx_Re_Grm_trans)
    if all(SP.ηα .== 0) return calc_σ_steadyState.(SP.Δ_range, Ref.(params)..., postfixes)
    else                return calc_σBα_steadyState.(SP.Δ_range, Ref.(params)..., postfixes)
    end
end


"""
Perform time evolution of the atomic and phononic degrees of freedom
"""
function timeEvolution(Δ, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2, N, initialState, tspan, dtmax, postfix)
    filename = "TE_" * postfix
    folder = "timeEvol/"
    
    if isfile(saveDir * folder * filename * ".txt") load_as_txt(saveDir * folder, filename) end
    
    args = empty_σVector(N), empty_BαVector(N), empty_σVector(N), empty_BαVector(N), 
           Δ, tildeΩ, tildeΩα, tildeG, tildeFα, tildeGα1, tildeGα2
    
    # Prepare the time evolution problem
    prob = ODEProblem(EoMs_wrap, initialState, tspan, args)
    
    # Perform the time evolution
    sol = OrdinaryDiffEq.solve(prob, Tsit5(), dtmax=dtmax)
    
    # Pack and save data
    formattedResult = vectorOfRows2Matrix([vcat(sol.t[i], sol.u[i]) for i in eachindex(sol.t)])
    save_as_txt(formattedResult, saveDir * folder, filename)
    
    return formattedResult
end


"""
Perform time evolution of the atomic and phononic degrees of freedom
(for the case of no phonons)
"""
function timeEvolution(Δ, tildeΩ, tildeG, N, initialState, tspan, dtmax, postfix)   
    filename = "TE_noPh_" * postfix
    folder = "timeEvol/"
    
    if isfile(saveDir * folder * filename * ".txt") load_as_txt(saveDir * folder, filename) end
    
    args = empty_σVector(N), empty_σVector(N), 
           Δ, tildeΩ, tildeG
    
    # Prepare the time evolution problem
    prob = ODEProblem(EoMs_wrap_noPh, initialState, tspan, args)
    
    # Perform the time evolution
    sol = OrdinaryDiffEq.solve(prob, Tsit5(), dtmax=dtmax)
    
    # Pack and save data
    formattedResult = vectorOfRows2Matrix([vcat(sol.t[i], sol.u[i]) for i in eachindex(sol.t)])
    save_as_txt(formattedResult, saveDir * folder, filename)
    
    return formattedResult
end


"""
Perform time evolution for parameters given by SP
"""
function timeEvolution(SP, Δ)
    postfix = get_postfix(Δ, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.arrayDescription, SP.fiber.postfix, SP.initialStateDescription, SP.tspan, SP.dtmax)
    params = get_parameterMatrices(SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.approx_Re_Grm_trans)
    return timeEvolution(Δ, params..., SP.N, SP.initialState, SP.tspan, SP.dtmax, postfix)
end
    

"""
Scan time evolutions over the detuning
"""
function scan_timeEvolution(SP)
    postfixes = get_postfix.(SP.Δ_range, Ref(SP.d), Ref(SP.να), Ref(SP.ηα), Ref(SP.incField_wlf), SP.arrayDescription, SP.fiber.postfix, SP.initialStateDescription, Ref(SP.tspan), SP.dtmax)
    params = get_parameterMatrices(SP.fiber, SP.d, SP.να, SP.ηα, SP.incField_wlf, SP.array, SP.approx_Re_Grm_trans)
    return timeEvolution(SP.Δ_range, Ref.(params)..., SP.N, Ref(SP.initialState), Ref(SP.tspan), SP.dtmax, postfixes)
end


#================================================
    Functions pertaining to transmission through the fiber
================================================#
"""
Calculate the transmission of light through the fiber (in the case of no phonons) 
in the same mode as used for the driving (note that there may be further transmitted
or reflected light in other modes)
"""
function calc_transmission(σ, fiber, d, incField_wlf, array)
    # postfix = get_postfix()
    # filename = "t_noPh_" * postfix
    
    # if isfile(saveDir * filename * ".txt")
    #     if overwrite_bool 
    #         println("The transmission for \n   $filename\nhas already been calculated.\n" *
    #                 "Recalculating and overwriting in 5 seconds...")
    #         sleep(5)
    #     else
    #         println("Loading transmission")
    #         return load_as_txt(saveDir, filename)
    #     end
    # end
    
    # Prepare some parameters and calculate transmission
    tildeΩ = get_tildeΩs(fiber, d, incField_wlf, array)
    t = transmission(σ, tildeΩ, fiber)
    
    # save_as_txt(t, saveDir, filename)
    return t
end


"""
Calculate the transmission of light through the fiber
in the same mode as used for the driving (note that there may be further transmitted
or reflected light in other modes)
"""
function calc_transmission(σBα, fiber, d, ηα, incField_wlf, array)
    # postfix = get_postfix()
    # filename = "t_" * postfix
    
    # if isfile(saveDir * filename * ".txt")
    #     if overwrite_bool 
    #         println("The transmission for \n   $filename\nhas already been calculated.\n" *
    #                 "Recalculating and overwriting in 5 seconds...")
    #         sleep(5)
    #     else
    #         println("Loading transmission")
    #         return load_as_txt(saveDir, filename)
    #     end
    # end
    
    # Prepare some parameters and calculate transmission
    tildeΩ, tildeΩα = get_tildeΩs(fiber, d, ηα, incField_wlf, array)
    t = transmission(σBα..., tildeΩ, tildeΩα, fiber)
    
    # save_as_txt(t, saveDir, filename)
    return t
end


"""
Calculate the transmission of light through the fiber in the chosen driving mode for parameters given by SP
    
The function assumes that σBα contains only σ if the Lamb-Dicke parameters are zero
"""
function calc_transmission(SP, σBα)
    if all(SP.ηα .== 0)
        return calc_transmission(σBα, SP.fiber, SP.d, SP.incField_wlf, SP.array)
    else
        return calc_transmission(σBα, SP.fiber, SP.d, SP.ηα, SP.incField_wlf, SP.array)
    end
end