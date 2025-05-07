
"""
Plot the propagation constant, κ, as a function of frequency, ω, 
for a specific value of the fiber radius, ρf, and the index of refraction, n.

It is assumed that the units of ρf and ω are nm and nm^-1 resp.
"""
function fig_propConst_vs_ω(ω_range, κ, ρf, n)
    # Start figure 
    fig = plot(reuse=false, size=(800, 600))

    # Plot the propagation constant, and the lines y = ω and y = nω
    plot!(ω_range, κ, label=L"$ y = \kappa(\omega) $", c=:black)
    plot!(ω_range, n * ω_range, label=L"$ y = n\omega $", c=:red)
    plot!(ω_range, ω_range, label=L"$ y = \omega $", c=:blue)

    # Finish figure
    plot!(ticks=:native)
    xlims!(extrema(ω_range))
    xlabel!(L"$ \omega $, [nm$^{-1}$]")
    ylabel!(L"[nm$^{-1}$]")
    title!("Propagation constant\n" *
           L"$ \rho_{f} = %$(round(ρf, sigdigits=3)) $nm, $ n = %$(round(n, sigdigits=3)) $")
    display(fig)
end


"""
Plot the inside/outside momenta of the fiber,
q = sqrt(κ^2 - ω^2), h = sqrt(n^2ω^2 - κ^2)

It is assumed that the units of ω nm^-1 resp.
"""
function fig_inout_momenta_vs_ω(ω_range, h, q, n)
    # Start figure 
    fig = plot(reuse=false, size=(800, 600))

    # Plot the propagation constant, and the lines y = ω and y = nω
    plot!(ω_range, h, label=L"$ h = \sqrt{n^2\omega^2 - \kappa^2} $", c=:blue)
    plot!(ω_range, q, label=L"$ q = \sqrt{\kappa^2 - \omega^2} $", c=:red)

    # Finish figure
    plot!(ticks=:native)
    xlims!(extrema(ω_range))
    xlabel!(L"$ \omega $, [nm$^{-1}$]")
    ylabel!(L"[nm$^{-1}$]")
    title!("Inside/outside momenta of the fiber (h and q) \n" *
           L"$ n = %$(round(n, sigdigits=3)) $")
    display(fig)
end


"""
Plot a coupling strength's magnitude and phase as a function of a relative spatial coordinate
"""
function fig_coupling_vs_x(x_range, coupling, x_label, y_label)
    # Start figure 
    fig = plot(reuse=false, size=(800, 600), layout=(1, 2))
    
    # Plot the propagation constant, and the lines y = ω and y = nω
    plot!(x_range, abs.(coupling) , label=L"$ |C| $" , c=:blue, subplot=1)
    plot!(x_range, angle.(coupling), label=L"arg$ (C) $", c=:red , subplot=2)

    # Finish figure
    plot!(ticks=:native)
    xlims!(extrema(x_range))
    xlabel!(x_label)
    ylabel!(y_label * ", magnitude", subplot=1)
    ylabel!(y_label * ", phase", subplot=1)
    display(fig)
end


"""
Plot time evolved atomic coherences σ
and the corresponding analytically calculated steady state values
(for the case of no phonons)
"""
function fig_σTrajectories_σSS(times, σTrajectories, σ_SS)
    colors = distinguishable_colors(length(σTrajectories), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    
    # Start figure 
    fig = plot(reuse=false, size=(800, 600))

    # Plot the trajectories
    for (i, traj) in enumerate(σTrajectories)
        plot!(times, real.(traj), label=false, c=colors[i], ls=:solid)
        plot!(times, imag.(traj), label=false, c=colors[i], ls=:dash )
        plot!([times[end] - (times[end] - times[1])/10, times[end]], real.(σ_SS[i])*ones(2),  label=false, c=colors[i], ls=:dashdot)
        plot!([times[end] - (times[end] - times[1])/10, times[end]], imag.(σ_SS[i])*ones(2),  label=false, c=colors[i], ls=:dashdot)
    end
    
    # Finish figure
    plot!(ticks=:native)
    xlims!(extrema(times))
    xlabel!(L"$ γt $")
    # ylabel!(L"")
    title!(L"$ σ $ trajectories")
    display(fig)
end


"""
Plot time evolved atomic coherences σ and the atom-phonon correlations Bα 
and the corresponding analytically calculated steady state values
"""
function fig_σBαTrajectories_σBαSS(times, σTrajectories, BαTrajectories, σ_SS, Bα_SS)
    colors = distinguishable_colors(length(σTrajectories) + 3*length(BαTrajectories[1]), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    
    # Start figure 
    fig = plot(reuse=false, size=(800, 600), layout=(1, 2))

    # Plot the trajectories
    for (i, traj) in enumerate(σTrajectories)
        plot!(times, real.(traj), label=false, c=colors[i], ls=:solid, subplot=1)
        plot!(times, imag.(traj), label=false, c=colors[i], ls=:dash , subplot=1)
        plot!([times[end] - (times[end] - times[1])/10, times[end]], real.(σ_SS[i])*ones(2),  label=false, c=colors[i], ls=:dashdot, subplot=1)
        plot!([times[end] - (times[end] - times[1])/10, times[end]], imag.(σ_SS[i])*ones(2),  label=false, c=colors[i], ls=:dashdot, subplot=1)
    end
    for α in 1:3
        clr_ind_offset = length(σTrajectories) + (α - 1)*length(BαTrajectories[1])
        for (i, traj) in enumerate(BαTrajectories[α])
            plot!(times, real.(traj), label=false, c=colors[clr_ind_offset + i], ls=:solid, subplot=2)
            plot!(times, imag.(traj), label=false, c=colors[clr_ind_offset + i], ls=:dash , subplot=2)
            plot!([times[end] - (times[end] - times[1])/10, times[end]], real.(Bα_SS[α][i])*ones(2), label=false, c=colors[clr_ind_offset + i], ls=:dashdot, subplot=2)
            plot!([times[end] - (times[end] - times[1])/10, times[end]], imag.(Bα_SS[α][i])*ones(2), label=false, c=colors[clr_ind_offset + i], ls=:dashdot, subplot=2)
        end
    end
    
    
    # Finish figure
    plot!(ticks=:native)
    xlims!(extrema(times))
    xlabel!(L"$ γt $")
    # ylabel!(L"")
    title!(L"$ σ $ trajectories", subplot=1)
    title!(L"$ B_α $ trajectories", subplot=2)
    display(fig)
end


"""
Plot magnitude and phase of transmission amplitude as a function of detuning
"""
function fig_transmission_vs_Δ(Δ_range, T, phase)
    # Start figure 
    fig = plot(reuse=false, size=(800, 600), layout=(1, 2))
    
    # Plot magnitude squared and the phase of the transmission 
    plot!(Δ_range, T    , label=L"$ |t|^2 $" , c=:blue, subplot=1)
    plot!(Δ_range, phase, label=L"arg$ (t) $", c=:red , subplot=2)

    # Finish figure
    plot!(ticks=:native)
    xlims!(extrema(Δ_range))
    # ylims!(ylims(fig[1])[1], 1, subplot=1)
    ylims!(0, 1, subplot=1)
    ylims!(-3.1415, 3.1415, subplot=2) #for some reason using π gives an error from Python
    xlabel!(L"$ \Delta/\gamma $")
    # ylabel!(L"")
    title!("Transmission coefficient", subplot=1)
    title!("Transmission phase", subplot=2)
    display(fig)
end


"""
Plot mean magnitude and phase of transmission amplitude as a function of detuning

with bands given by the standard deviation
"""
function fig_classDisorder_transmission_vs_Δ(Δ_range, T_means, T_stds, phase_means, phase_stds)
    # Start figure 
    fig = plot(reuse=false, size=(800, 600), layout=(1, 2))
    
    # Plot mean plus/minus the standard deviation of the magnitude squared and the phase of the transmission 
    plot!(Δ_range, T_means, ribbon=T_stds, fillalpha=0.35, label=L"$ |t|^2 $" , c=:blue, subplot=1)
    plot!(Δ_range, phase_means, ribbon=phase_stds, fillalpha=0.35, label=L"arg$ (t) $", c=:red , subplot=2)
    
    # Finish figure
    plot!(ticks=:native)
    xlims!(extrema(Δ_range))
    # ylims!(ylims(fig[1])[1], 1, subplot=1)
    ylims!(0, 1, subplot=1)
    ylims!(-3.1415, 3.1415, subplot=2) #for some reason using π gives an error from Python
    xlabel!(L"$ \Delta/\gamma $")
    # ylabel!(L"")
    title!("Transmission coefficient", subplot=1)
    title!("Transmission phase", subplot=2)
    display(fig)
end


"""
Plot the intensity (norm-squared) of the radiated E-field around the fiber
"""
function fig_E_radiation(z_range, x_range, intensity, ρf, array)
    # Start figure 
    fig = plot(reuse=false, size=(800, 600))
    
    # Plot the E-field intensity
    sat = 1e-2
    intensity[intensity .> sat] .= sat
    # contourf!(z_range, x_range, intensity', c=:viridis)
    heatmap!(z_range, x_range, intensity', c=:viridis)
    
    # Plot a representation of the fiber
    rectangle(w, h, x, y) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])
    plot!(rectangle(z_range[end] - z_range[1], 2*ρf, z_range[1], -ρf), lw=2, c=:gray, opacity=0.3, label=false)
    
    # Plot a representation of the atomic array
    scatter!([site[3] for site in array], [site[1] for site in array], mc=:black, markershape=:circle, ms=10, label=false)
    
    # Finish figure
    plot!(ticks=:native, aspect_ratio=:equal)
    xlims!(extrema(z_range))
    ylims!(extrema(x_range))
    xlabel!(L"$ z/λ $")
    ylabel!(L"$ x/λ $")
    title!(L"$ I/(\gamma/\lambda^3) $")
    display(fig)
end