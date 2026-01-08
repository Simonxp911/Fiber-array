

"""
EoMs: equations of motion (assumed to take the argument (dxdt_container, x, args, t))
x0: initial state
Δt: fixed time step
EoMs_args: arguments for EoMs
timeEvol_args: arguments for timeEvol, stepCondition, and stepFunc
stepCondition: condition to continue time evolution (evaluated at every step)
stepFunc: function that is evaluated at every step (and its value returned)
returnAll: if true, returns all times, states, and stepFuncVal (otherwise only the final step is returned)
"""
function timeEvol(EoMs, x0, EoMs_args, timeEvol_args, stepCondition=stepCondition_endOftspan, stepFunc=stepFunc_nothing, returnAll=true)
    # Prepare a function f with arguments t and x(t) that calculate the derivative at that time and with that state
    dxdt_container = zeros(length(x0))
    dxdt_function(t, x) = EoMs(dxdt_container, x, EoMs_args, t)
    
    # Prepare variables and containers
    t = timeEvol_args.tspan[1]
    Δt = timeEvol_args.dtmax
    xt = deepcopy(x0)
    ΔxΔt = zeros(length(x0))
    xt_intermediate = zeros(length(x0))
    
    # Prepare first value of ΔxΔt and stepFuncVal (to make stepCondition work and to have an initial stepFuncVal)
    RK4!(ΔxΔt, t, xt, Δt, dxdt_container, dxdt_function, xt_intermediate)
    stepFuncVal = stepFunc(t, xt, ΔxΔt, timeEvol_args)
    
    if returnAll
        t_all = [Float64(t)]
        x_all = [deepcopy(x0)]
        stepFuncVal_all = [stepFuncVal]
    end
    
    # A while loop is run to perform the time evolution using 4th order Runge-Kutta 
    while stepCondition(t, xt, ΔxΔt, stepFuncVal, timeEvol_args)
        # Calculate the finite-step state derivative 
        RK4!(ΔxΔt, t, xt, Δt, dxdt_container, dxdt_function, xt_intermediate)
    
        # Step the time and state forward
        t   += Δt
        xt .+= ΔxΔt.*Δt
        
        # Evaluate the stepFunc
        stepFuncVal = stepFunc(t, xt, ΔxΔt, timeEvol_args)
        
        if returnAll
            push!(t_all, t)
            push!(x_all, deepcopy(xt))
            push!(stepFuncVal_all, deepcopy(stepFuncVal))
        end
    end
    
    if returnAll return (t=t_all, u=x_all, stepFuncVal=stepFuncVal_all)
    else         return (t=[t]  , u=[xt] , stepFuncVal=[stepFuncVal]) 
    end
end


"""
4th order Runge-Kutta calculation of ΔxΔt (performed in-place)
"""
function RK4!(ΔxΔt, t, xt, Δt, dxdt_container, dxdt_function, xt_intermediate)
    # Calculate the finite-step state derivative ΔxΔt using 4th order Runge-Kutta
    # The code below implements the following steps in a way that is optimized to not initiate any new Arrays
    # k1 = dxdt_function(t, xt)
    # k2 = dxdt_function(t + Δt/2, xt + k1*Δt/2)
    # k3 = dxdt_function(t + Δt/2, xt + k2*Δt/2)
    # k4 = dxdt_function(t + Δt, xt + k3*Δt)
    # @. ΔxΔt = (k1 + 2*k2 + 2*k3 + k4)/6
    # The above only works if dxdt_function initiates a new dxdt_container every time it is called
    
    
    # Calculate the k variables
    dxdt_function(t, xt)
    @. ΔxΔt  = dxdt_container/6
    
    @. xt_intermediate = xt + dxdt_container*Δt/2
    dxdt_function(t + Δt/2, xt_intermediate)
    @. ΔxΔt += dxdt_container/3
    
    @. xt_intermediate = xt + dxdt_container*Δt/2
    dxdt_function(t + Δt/2, xt_intermediate)
    @. ΔxΔt += dxdt_container/3
    
    @. xt_intermediate = xt + dxdt_container*Δt
    dxdt_function(t + Δt, xt_intermediate)
    @. ΔxΔt += dxdt_container/6
end


"""
To stop time evolution at the end of tspan
"""
function stepCondition_endOftspan(t, xt, ΔxΔt, stepFuncVal, timeEvol_args)
    return t <= timeEvol_args.tspan[2]
end


"""
To stop time evolution at the steady state, as defined by 
timeEvol_args.stateDiffAbsTol and timeEvol_args.stateDiffRelTol
"""
function stepCondition_reachSS(t, xt, ΔxΔt, stepFuncVal, timeEvol_args)
    return (t <= timeEvol_args.tspan[2]) && (any(ΔxΔt .> timeEvol_args.stateDiffAbsTol) || any(ΔxΔt.*xt .> timeEvol_args.stateDiffRelTol))
end


"""
To stop time evolution when stepFuncVal is small, as defined by timeEvol_args.stepFuncValLowerTol,
and one-tenth of tspan has passed.
If stepFuncVal is an iterable (and timeEvol_args.stepFuncValLowerTol is too), all values are compared
and time evolution is stopped when all of them are small
"""
function stepCondition_stepFuncVal_isSmall(t, xt, ΔxΔt, stepFuncVal, timeEvol_args)
    return (t < timeEvol_args.tspan[1] + (timeEvol_args.tspan[2] - timeEvol_args.tspan[1])/10) || any(stepFuncVal .>= timeEvol_args.stepFuncValLowerTol)
end


"""
Default stepFunc, which simply returns nothing
"""
function stepFunc_nothing(t, xt, ΔxΔt, timeEvol_args)
    return nothing
end

