

function time_evol(EoMs, x0, tspan, Δt, args, evolToSS=false, returnAll=false, abstol=1e-5, reltol=1e-3)
    # evolToSS: if true, the evolution stops when the steady state is reached (as defined by abstol and reltol), and otherwise it stops at the end of tspan
    # returnAll: if true, returns all states found while evolving over tspan (otherwise only the final state is returned)
        # Notice that copying xt to the output, x_all, usually takes the same amount of time as actually calculating xt!
    
    # We assume EoMs takes the arguments (dxdt_container, x, args, t)
    # Prepare a function f with arguments t and x(t) that calculate the derivative at that time and with that state
    dxdt_container = zeros(length(x0))
    dxdt_function(t, x) = EoMs(dxdt_container, x, args, t)
     
    # Prepare a container for the state at all times if needed
    if returnAll
        t_all = [Float64(tspan[1])]
        x_all = [deepcopy(x0)]
    end
    
    # Prepare the condition for continuing the time-evolution
    function step_condition(t, xt, ΔxΔt)
        if evolToSS
            return (t <= tspan[2]) && (any(ΔxΔt .> abstol) || any(ΔxΔt.*xt .> reltol))
        else
            return t <= tspan[2]
        end
    end
    
    # A while loop is run to perform the time evolution using 4th order Runge-Kutta 
    # Initialize the time, state, and finite-step state derivative
    t = tspan[1]
    xt = deepcopy(x0)
    xt_intermediate = zeros(length(x0))
    ΔxΔt = ones(length(x0)) #use ones to make sure condition = true
    while step_condition(t, xt, ΔxΔt)
        # Calculate the finite-step state derivative 
        RK4!(ΔxΔt, t, xt, Δt, dxdt_container, dxdt_function, xt_intermediate)
    
        # Step the time and state forward
        t   += Δt
        xt .+= ΔxΔt.*Δt
        
        if returnAll
            push!(t_all, t)
            push!(x_all, deepcopy(xt))
        end
    end
    
    if returnAll
        return (t=t_all, u=x_all)
    else
        return (t=[t], u=[xt])
    end
end


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