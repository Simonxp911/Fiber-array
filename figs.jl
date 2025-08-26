

"""
Plot the propagation constant, κ, as a function of frequency, ω, 
for a specific value of the fiber radius, ρf, and the index of refraction, n.

It is assumed that the units of ρf and ω are nm and nm^-1 resp.
"""
function fig_propConst_vs_ω(ω_range, κ, ρf, n)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1], "Propagation constant\n" *
                    L"$ \rho_{f} = %$(round(ρf, sigdigits=3)) $nm, $ n = %$(round(n, sigdigits=3)) $")
    ax1 = Axis(fig[2, 1], limits=(extrema(ω_range), nothing), 
               xlabel=L"$ ω $, [nm$^{-1}$]", 
               ylabel=L"[nm$^{-1}$]")
    
    # Plot the propagation constant, and the lines y = ω and y = nω
    lines!(ax1, ω_range, κ          , label=L"$ y = κ(ω) $", color=:black)
    lines!(ax1, ω_range, n * ω_range, label=L"$ y = nω $"       , color=:red)
    lines!(ax1, ω_range, ω_range    , label=L"$ y = ω $"        , color=:blue)
    
    # Finish figure
    axislegend()
    display(GLMakie.Screen(), fig)
end


"""
Plot the inside/outside momenta of the fiber,
q = sqrt(κ^2 - ω^2), h = sqrt(n^2ω^2 - κ^2)

It is assumed that the units of ω nm^-1 resp.
"""
function fig_inout_momenta_vs_ω(ω_range, h, q, n)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1], "Inside/outside momenta of the fiber (h and q) \n" *
                    L"$ n = %$(round(n, sigdigits=3)) $")
    Axis(fig[2, 1], limits=(extrema(ω_range), nothing), 
                    xlabel=L"[nm$^{-1}$]", 
                    ylabel=L"$ ω $, [nm$^{-1}$]")
    
    # Plot the momenta
    lines!(ω_range, h, label=L"$ h = \sqrt{n^2ω^2 - \kappa^2} $", color=:blue)
    lines!(ω_range, q, label=L"$ q = \sqrt{\kappa^2 - ω^2} $", color=:red)
    
    # Finish figure
    axislegend()
    display(GLMakie.Screen(), fig)
end


"""
Plot a coupling strength's magnitude and phase as a function of a relative spatial coordinate
"""
function fig_coupling_vs_x(x_range, coupling, x_label, y_label)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make axes
    ax1 = Axis(fig[1, 1], limits=(extrema(x_range), nothing), 
               xlabel=x_label, 
               ylabel=latexstring(y_label * ", magnitude"))
    ax2 = Axis(fig[1, 2], limits=(extrema(x_range), nothing), 
               xlabel=x_label, 
               ylabel=latexstring(y_label * ", phase"))
    
    # Plot the coupling's magnitude and phase
    lines!(ax1, x_range, abs.(coupling) , label=L"$ |C| $" , color=:blue)
    lines!(ax2, x_range, angle.(coupling), label=L"arg$ (C) $", color=:red)
    
    # Finish figure
    axislegend.([ax1, ax2])
    display(GLMakie.Screen(), fig)
end


"""
Plot the atomic array and the fiber in 3D
"""
function fig_arrayIn3D(array, x_range, z_range, ρf)
    # Start figure 
    fig = Figure(size=(900, 600))
    
    # Make title and axis
    Label(fig[1, 1], L"The atomic array and fiber$$", tellwidth=false)
    zWidth  = maximum(z_range) - minimum(z_range)
    xHeight = maximum(x_range) - minimum(x_range)
    Axis3(fig[2, 1], limits=(extrema(z_range), extrema(x_range), extrema(x_range)), 
                     yreversed = true,
                     xlabel=L"$ z/λ_{a} $", 
                     ylabel=L"$ y/λ_{a} $", 
                     zlabel=L"$ x/λ_{a} $", 
                     aspect=(zWidth, xHeight, xHeight)./maximum((zWidth, xHeight)))
    
    # Plot the atoms
    radius = 0.1
    θs = range(0, π, 20)
    φs = range(0, 2π, 20)
    xSph = radius.*[cos(φ)*sin(θ) for θ in θs, φ in φs]
    ySph = radius.*[sin(φ)*sin(θ) for θ in θs, φ in φs]
    zSph = radius.*[cos(θ) for θ in θs, φ in φs]
    for site in array
        surface!(zSph .+ site[3], ySph .+ site[2], xSph .+ site[1], colormap=:grays)
    end
    
    # Plot the fiber
    xCyl = ρf.*[cos(φ) for z in z_range, φ in φs]
    yCyl = ρf.*[sin(φ) for z in z_range, φ in φs]
    zCyl = [z for z in z_range, φ in φs]
    surface!(zCyl, yCyl, xCyl, colormap=(:grays, 0.3))
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot time evolved atomic coherences σ
and the corresponding analytically calculated steady state values
(for the case of no phonons)
"""
function fig_σTrajectories_σSS(times, σTrajectories, σ_SS)
    colors = distinguishable_colors(length(σTrajectories), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1], L"$ σ $ trajectories", tellwidth=false)
    Axis(fig[2, 1], limits=(extrema(times), nothing), 
                    xlabel=L"$ γ_{a}t $")
    
    # Plot the trajectories
    for (i, traj) in enumerate(σTrajectories)
        lines!(times, real.(traj), color=colors[i], linestyle=:solid)
        lines!(times, imag.(traj), color=colors[i], linestyle=:dash )
        lines!([times[end] - (times[end] - times[1])/10, times[end]], real.(σ_SS[i])*ones(2), color=colors[i], linestyle=:dashdot)
        lines!([times[end] - (times[end] - times[1])/10, times[end]], imag.(σ_SS[i])*ones(2), color=colors[i], linestyle=:dashdot)
    end
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot time evolved atomic coherences σ and the atom-phonon correlations Bα 
and the corresponding analytically calculated steady state values
"""
function fig_σBαTrajectories_σBαSS(times, σTrajectories, BαTrajectories, σ_SS, Bα_SS)
    colors = distinguishable_colors(length(σTrajectories) + 3*length(BαTrajectories[1]), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1], L"$ σ $ trajectories", tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(times), nothing), 
               xlabel=L"$ γ_{a}t $")
    Label(fig[1, 2], L"$ B_α $ trajectories", tellwidth=false)
    ax2 = Axis(fig[2, 2], limits=(extrema(times), nothing), 
               xlabel=L"$ γ_{a}t $")
    
    # Plot the trajectories
    for (i, traj) in enumerate(σTrajectories)
        lines!(ax1, times, real.(traj), color=colors[i], linestyle=:solid)
        lines!(ax1, times, imag.(traj), color=colors[i], linestyle=:dash)
        lines!(ax1, [times[end] - (times[end] - times[1])/10, times[end]], real.(σ_SS[i])*ones(2),  color=colors[i], linestyle=:dashdot)
        lines!(ax1, [times[end] - (times[end] - times[1])/10, times[end]], imag.(σ_SS[i])*ones(2),  color=colors[i], linestyle=:dashdot)
    end
    for α in 1:3
        clr_ind_offset = length(σTrajectories) + (α - 1)*length(BαTrajectories[1])
        for (i, traj) in enumerate(BαTrajectories[α])
            lines!(ax2, times, real.(traj), color=colors[clr_ind_offset + i], linestyle=:solid)
            lines!(ax2, times, imag.(traj), color=colors[clr_ind_offset + i], linestyle=:dash)
            lines!(ax2, [times[end] - (times[end] - times[1])/10, times[end]], real.(Bα_SS[α][i])*ones(2), color=colors[clr_ind_offset + i], linestyle=:dashdot)
            lines!(ax2, [times[end] - (times[end] - times[1])/10, times[end]], imag.(Bα_SS[α][i])*ones(2), color=colors[clr_ind_offset + i], linestyle=:dashdot)
        end
    end
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot magnitude and phase of transmission amplitude as a function of detuning
"""
function fig_transmission_vs_Δ(Δ_range, T, phase, titl)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    Label(fig[2, 1], L"Transmission coefficient, $ |t|^2 $", tellwidth=false)
    ax1 = Axis(fig[3, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/γ_{a} $")
    Label(fig[2, 2], L"Transmission phase, arg$ (t) $", tellwidth=false)
    if all(-π .<= phase .<= π)
        ax2 = Axis(fig[3, 2], limits=(extrema(Δ_range)..., -π, π),
                yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
                xlabel=L"$ \Delta/γ_{a} $")
    else
        ax2 = Axis(fig[3, 2], limits=(extrema(Δ_range)..., nothing, nothing),
                xlabel=L"$ \Delta/γ_{a} $")
    end
               
    # Plot magnitude squared and the phase of the transmission 
    lines!(ax1, Δ_range, T    , color=:blue)
    lines!(ax2, Δ_range, phase, color=:red)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot magnitude, phase, unfolded phase, unfolded phase divided by N, and slope of phase
of transmission amplitude as a function of detuning
"""
function fig_transmission_vs_Δ_phaseDetails(Δ_range, T, phase, unwrappedPhase, phasePerAtom, phaseSlope, titl)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    Label(fig[2, 1], L"Transmission coefficient, $ |t|^2 $", tellwidth=false)
    ax1 = Axis(fig[3, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/γ_{a} $")
    Label(fig[2, 2], L"Transmission phase, arg$ (t) $", tellwidth=false)
    ax2 = Axis(fig[3, 2], limits=(extrema(Δ_range)..., -π, π),
            yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
            xlabel=L"$ \Delta/γ_{a} $")
    
    Label(fig[4, 1], L"Unwrapped phase, arg$ (t) $", tellwidth=false)
    ax3 = Axis(fig[5, 1], limits=(extrema(Δ_range)..., nothing, nothing),
            xlabel=L"$ \Delta/γ_{a} $")
    # ax4 = Axis(fig[5, 1], limits=(extrema(Δ_range)..., nothing, nothing),
    #         ylabel=L"Phase per atom, arg$ (t)/N $",
    #         yaxisposition=:right)
    # hidespines!(ax4)
    # hidexdecorations!(ax4)
    Label(fig[4, 2], L"Phase slope, arg$ '(t) $", tellwidth=false)
    ax5 = Axis(fig[5, 2], limits=(extrema(Δ_range)..., nothing, nothing),
            xlabel=L"$ \Delta/γ_{a} $")
               
    # Plot magnitude squared and the phase of the transmission 
    lines!(ax1, Δ_range, T    , color=:blue)
    lines!(ax2, Δ_range, phase, color=:red)
    
    # Plot unwrapped phase, phase per atom, and slope of phase
    lines!(ax3, Δ_range, unwrappedPhase, color=:red)
    # lines!(ax4, Δ_range, phasePerAtom, color=:purple)
    lines!(ax5, Δ_range[1:end-1] .+ diff(Δ_range)./2, phaseSlope, color=:black)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot magnitude and phase of the transmission and reflection amplitudes as a function of detuning
"""
function fig_transmissionAndReflection_vs_Δ(Δ_range, T, tPhase, R, rPhase, titl)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    Label(fig[2, 1], L"Transmission coefficient, $ |t|^2 $", tellwidth=false)
    ax1 = Axis(fig[3, 1], limits=(extrema(Δ_range)..., 0, 1))
    Label(fig[2, 2], L"Transmission phase, arg$ (t) $", tellwidth=false)
    ax2 = Axis(fig[3, 2], limits=(extrema(Δ_range)..., -π, π),
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]))
    Label(fig[4, 1], L"Reflection coefficient, $ |r|^2 $", tellwidth=false)
    ax3 = Axis(fig[5, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/γ_{a} $")
    Label(fig[4, 2], L"Reflection phase, arg$ (r) $", tellwidth=false)
    ax4 = Axis(fig[5, 2], limits=(extrema(Δ_range)..., -π, π),
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=L"$ \Delta/γ_{a} $")
               
    # Plot magnitudes squared and the phases
    lines!(ax1, Δ_range, T     , color=:blue)
    lines!(ax2, Δ_range, tPhase, color=:red)
    lines!(ax3, Δ_range, R     , color=:blue)
    lines!(ax4, Δ_range, rPhase, color=:red)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot mean magnitude and phase of transmission amplitude as a function of detuning
with bands given by the standard deviation
"""
function fig_imperfectArray_transmission_vs_Δ(Δ_range, T_means, T_stds, phase_means, phase_stds, titl)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    Label(fig[2, 1], L"Transmission coefficient, $ |t|^2 $", tellwidth=false)
    ax1 = Axis(fig[3, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/γ_{a} $")
    Label(fig[2, 2], L"Transmission phase, arg$ (t) $", tellwidth=false)
    ax2 = Axis(fig[3, 2], limits=(extrema(Δ_range)..., -π, π), 
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=L"$ \Delta/γ_{a} $")
               
    # Plot magnitude squared and the phase of the transmission with bands for standard deviations
    lines!(ax1, Δ_range, T_means , color=:blue)
    band!( ax1, Δ_range, T_means + T_stds, T_means - T_stds , color=(:blue, 0.35))
    lines!(ax2, Δ_range, phase_means, color=:red)
    band!( ax2, Δ_range, phase_means + phase_stds, phase_means - phase_stds , color=(:red, 0.35))
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot mean magnitude and phase of transmission amplitude as a function of detuning
with bands given by the standard deviation
for several values of ff and ηα
"""
function fig_compareImperfectArray_transmission_vs_Δ(Δ_range, T_meanss, T_stdss, phase_meanss, phase_stdss, labels, titl)
    colors = distinguishable_colors(length(labels), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    Label(fig[2, 1], L"Transmission coefficient, $ |t|^2 $", tellwidth=false)
    ax1 = Axis(fig[3, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/γ_{a} $")
    Label(fig[2, 2], L"Transmission phase, arg$ (t) $", tellwidth=false)
    ax2 = Axis(fig[3, 2], limits=(extrema(Δ_range)..., -π, π), 
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=L"$ \Delta/γ_{a} $")
    
    # Plot magnitude squared and the phase of the transmission with bands for standard deviations
    for (i, label) in enumerate(labels)
        lines!(ax1, Δ_range, T_meanss[i], label=label , color=colors[i])
        band!( ax1, Δ_range, T_meanss[i] + T_stdss[i], T_meanss[i] - T_stdss[i] , color=(colors[i], 0.35))
        lines!(ax2, Δ_range, phase_meanss[i], color=colors[i])
        band!( ax2, Δ_range, phase_meanss[i] + phase_stdss[i], phase_meanss[i] - phase_stdss[i] , color=(colors[i], 0.35))
    end
    
    # Finish figure
    axislegend(ax1, position=:lb)
    display(GLMakie.Screen(), fig)
end


"""
Plot mean magnitude and phase of transmission amplitude for a specific detuning
as a function of any quantity x
with error bars given by the standard deviation
"""
function fig_compareImperfectArray_transmission_vs_X(xs, x_label, T_means, T_stds, phase_means, phase_stds, T_indepDecays, phase_indepDecays, titl)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    Label(fig[2, 1], L"Transmission coefficient, $ |t|^2 $", tellwidth=false)
    # ax1 = Axis(fig[3, 1], limits=(nothing, nothing, 0, 1), 
    #            xlabel=x_label)
    ax1 = Axis(fig[3, 1], xlabel=x_label, yscale=log10)
    Label(fig[2, 2], L"Transmission phase, arg$ (t) $", tellwidth=false)
    ax2 = Axis(fig[3, 2], limits=(nothing, nothing, -π, π), 
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=x_label)
    
    # Plot magnitude squared and the phase of the transmission with bands for standard deviations
    errorbars!(ax1, xs, T_means, T_stds, color=:black, whiskerwidth=10, label=L"Incl. inter.$$")
    scatter!(ax1, xs, T_means, color=:black, markersize=6)
    errorbars!(ax2, xs, phase_means, phase_stds, color=:black, whiskerwidth=10)
    scatter!(ax2, xs, phase_means, color=:black, markersize=6)
    
    # Compare with independent decay case
    lines!(ax1, xs, T_indepDecays, color=:red, label=L"Indep. decay$$")
    lines!(ax2, xs, phase_indepDecays, color=:red)
    
    # Finish figure
    axislegend(ax1, position=:rt)
    display(GLMakie.Screen(), fig)
end


"""
Plot mean magnitude and phase of transmission amplitude for a specific detuning and multiple ff
as a function of N
with error bars given by the standard deviation
"""
function fig_compareImperfectArray_transmission_vs_N(Ns, T_means, T_stds, phase_means, phase_stds, T_indepDecays, phase_indepDecays, T_fits, phase_fits, β_indepDecay, labels, titl)
    # Prepare colors
    colors = distinguishable_colors(length(labels), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Start figure 
    fig = Figure(size=(1200, 700))
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    Label(fig[2, 1], L"Transmission coefficient, $ |t|^2 $", tellwidth=false)
    ax1 = Axis(fig[3, 1], xlabel=L"$ N $", yscale=log10)
    Label(fig[2, 2], L"Transmission phase, arg$ (t) $", tellwidth=false)
    ax2 = Axis(fig[3, 2], xlabel=L"$ N $")
    
    # Compare with independent decay case
    lines!(ax1, Ns, T_indepDecays, color=:black, label=L"Indep. decay, $ β = %$(round(β_indepDecay, digits=3)) $")
    lines!(ax2, Ns, phase_indepDecays, color=:black)
    
    # Plot magnitude squared and the phase of the transmission with bands for standard deviations
    for (i, label) in enumerate(labels)
        # Sometimes the errorbar reaches into the negative domain, which the logplot cannot handle, so instead the lower error is set at T*0.001
        lowerrors = T_stds[i, :]
        for (j, error) in enumerate(lowerrors)
            if T_means[i, j] - error < 0
                lowerrors[j] = T_means[i, j]*0.999
            end
        end
        
        errorbars!(ax1, Ns, T_means[i, :], lowerrors, T_stds[i, :], color=colors[i], whiskerwidth=10, label=label)
        scatter!(ax1, Ns, T_means[i, :], color=colors[i], markersize=6)
        lines!(ax1, Ns, T_fits[i, :], color=colors[i])
        errorbars!(ax2, Ns, phase_means[i, :], phase_stds[i, :], color=colors[i], whiskerwidth=10)
        scatter!(ax2, Ns, phase_means[i, :], color=colors[i], markersize=6)
        lines!(ax2, Ns, phase_fits[i, :], color=colors[i])
    end
    
    # Finish figure
    axislegend(ax1, position=:lb)
    display(GLMakie.Screen(), fig)
end


"""
Plot the effective β-factor as a function of detuning for multiple ff
"""
function fig_βfactor_vs_Δ(Δ_range, βs, β_indepDecay, labels, titl)
    # Prepare colors
    colors = distinguishable_colors(length(labels), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1], titl, tellwidth=false)
    ax1 = Axis(fig[2, 1], xlabel=L"$ Δ/γ_{a} $", ylabel=L"$ β $-factor")
    
    # Compare with independent decay case
    hlines!(ax1, β_indepDecay, color=:black, label=L"Indep. decay$$")
    
    # Plot magnitude squared and the phase of the transmission with bands for standard deviations
    for (i, label) in enumerate(labels)
        lines!(ax1, Δ_range, βs[i, :], color=colors[i], label=label)
    end
    
    # Finish figure
    axislegend(ax1, position=:lb)
    display(GLMakie.Screen(), fig)
end


"""
Plot the effective decay rates as a function of detuning for multiple ff
"""
function fig_effectiveβΔ_vs_Δ(Δ_range, β_effs, Δ_effs, β_indepDecay, labels, titl)
    # Prepare colors
    colors = distinguishable_colors(length(labels), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    Label(fig[2, 1], L"$ β_{eff} $", tellwidth=false)
    ax1 = Axis(fig[3, 1], xlabel=L"$ Δ/γ_{a} $")
    Label(fig[2, 2], L"$ [Δ/(γ_{gm} + γ_{rm})]_{eff} $", tellwidth=false)
    ax2 = Axis(fig[3, 2], xlabel=L"$ Δ/γ_{a} $")
    
    # Compare with independent decay case
    hlines!(ax1, β_indepDecay, color=:black, label=L"Indep. decay$$")
    lines!(ax2, Δ_range, Δ_range, color=:black)
    
    # Plot magnitude squared and the phase of the transmission with bands for standard deviations
    for (i, label) in enumerate(labels)
        lines!(ax1, Δ_range, β_effs[i, :], color=colors[i], label=label)
        lines!(ax2, Δ_range, Δ_effs[i, :], color=colors[i])
    end
    
    # Finish figure
    axislegend(ax1, position=:lt)
    display(GLMakie.Screen(), fig)
end


"""
Plot the intensity (norm-squared) of the radiated E-field around the fiber
"""
function fig_radiation_Efield(z_range, x_range, intensity, ρf, array)
    # Start figure 
    zWidth  = maximum(z_range) - minimum(z_range)
    xHeight = maximum(x_range) - minimum(x_range)
    fig = Figure(size=(zWidth/xHeight*300, 300))
    
    # Make title and axis
    Label(fig[1, 1], L"$ I/(γ_{a}/λ_{a}^3) $", tellwidth=false)
    Axis(fig[2, 1], limits=(extrema(z_range), extrema(x_range)), 
                    xlabel=L"$ z/λ_{a} $", 
                    ylabel=L"$ x/λ_{a} $", 
                    aspect=DataAspect())
    
    # Plot the E-field intensity
    # sat = 1e1
    # intensity[intensity .> sat] .= sat
    # contourf!(z_range, x_range, intensity, colormap =:viridis)
    hm = heatmap!(z_range, x_range, intensity, colormap =:viridis)
    Colorbar(fig[2, 2], hm)
    
    # Plot a representation of the fiber
    poly!([(z_range[1], -ρf), (z_range[end], -ρf), (z_range[end], ρf), (z_range[1], ρf)], color=(:gray, 0.3))
    
    # Plot a representation of the atomic array
    scatter!([site[3] for site in array], [site[1] for site in array], color=:black, marker=:circle, markersize=10)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot a (complex) function (represented by an N-vector v) on a 1D chain,
as well as its discrete Fourier transform.
"""
function fig_fOnChain(rs, v, ks, vFT)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 2], L"Real space$$", tellwidth=false)
    Label(fig[1, 3], L"Momentum space$$", tellwidth=false)
    Label(fig[2, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[3, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[2, 2])
    ax2 = Axis(fig[2, 3])
    ax3 = Axis(fig[3, 2], xlabel=L"$ z/λ_{a} $")
    ax4 = Axis(fig[3, 3], xlabel=L"$ λ_{a}k_z $")
    
    # Plot the real space function
    lines!(ax1, rs, abs.(v)  , color=:blue)
    lines!(ax3, rs, angle.(v), color=:blue)
    
    # Plot the k-space function
    lines!(ax2, ks, abs.(vFT)  , color=:red)
    lines!(ax4, ks, angle.(vFT), color=:red)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot a (complex) function (represented by an (N, N)-matrix M) on a 2D square lattice,
as well as its discrete Fourier transform.
"""
function fig_fOnSquare(rs, M, ks, MFT)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 2], L"Real space$$", tellwidth=false)
    Label(fig[1, 4], L"Momentum space$$", tellwidth=false)
    Label(fig[2, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[3, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[2, 2])
    ax2 = Axis(fig[2, 4])
    ax3 = Axis(fig[3, 2], xlabel=L"$ z/λ_{a} $")
    ax4 = Axis(fig[3, 4], xlabel=L"$ λ_{a}k_z $")
    
    # Plot the real space function
    hm1 = heatmap!(ax1, rs, rs, abs.(M), colormap=:viridis)
    hm3 = heatmap!(ax3, rs, rs, angle.(M), colormap=:RdBu)
    Colorbar(fig[2, 3], hm1)
    Colorbar(fig[3, 3], hm3)
    
    # Plot the k-space function
    hm2 = heatmap!(ax2, ks, ks, abs.(MFT), colormap=:viridis)
    hm4 = heatmap!(ax4, ks, ks, angle.(MFT), colormap=:RdBu)
    Colorbar(fig[2, 5], hm2)
    Colorbar(fig[3, 5], hm4)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot the eigenvectors of a coupling matrix in real space, as well as their Fourier transform.
Furthermore, plot the emission pattern of that mode.
"""
function fig_GnmEigenModes(zs, eigen_σ, ks, eigen_σ_FT, z_range, x_range, intensity, ρf, array, collΔ, collΓ_gm, collΓ_rm, κ, titl)
    # Start figure 
    fig = Figure(size=(800, 700))
    
    # Make title and axis
    Label(fig[1, 1:3], titl, tellwidth=false)
    Label(fig[2, 2:3], L"$ Δ_{i}/γ_{a} $, $ Γ_{i, gm}/γ_{a} $, $ Γ_{i, rm}/γ_{a} $ = %$(round(collΔ, digits=3)), %$(round(collΓ_gm, digits=3)), %$(round(collΓ_rm, digits=3)), ", tellwidth=false)
    Label(fig[3, 2], L"Real space$$", tellwidth=false)
    Label(fig[3, 3], L"Momentum space$$", tellwidth=false)
    Label(fig[4, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[5, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[4, 2])
    ax2 = Axis(fig[4, 3])
    ax3 = Axis(fig[5, 2], xlabel=L"$ z/λ_{a} $")
    ax4 = Axis(fig[5, 3], xlabel=L"$ λ_{a}k_z $")
    ax5 = Axis(fig[6, 2:3], aspect=DataAspect(), xlabel=L"$ z/λ_{a} $", ylabel=L"$ x/λ_{a} $")
    
    # Plot real space 
    lines!(ax1, zs, abs.(eigen_σ)  , color=:blue)
    lines!(ax3, zs, angle.(eigen_σ), color=:blue)
    
    # Plot k-space 
    lines!(ax2, ks, abs.(eigen_σ_FT)  , color=:red)
    lines!(ax4, ks, angle.(eigen_σ_FT), color=:red)
    
    # Mark the light cone and position of propagation constant
    vlines!(ax2, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
    vlines!(ax2, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    axislegend(ax2, position=:ct)
    
    # Plot the E-field intensity
    hm = heatmap!(ax5, z_range, x_range, intensity, colormap=:viridis)
    Colorbar(fig[6, 4], hm)
    
    # Plot a representation of the fiber
    poly!(ax5, [(z_range[1], -ρf), (z_range[end], -ρf), (z_range[end], ρf), (z_range[1], ρf)], color=(:gray, 0.3))
    
    # Plot a representation of the atomic array
    scatter!(ax5, [site[3] for site in array], [site[1] for site in array], color=:black, marker=:circle, markersize=10)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot the eigenvectors of a coupling matrix in real space, as well as their Fourier transform.
Furthermore, plot the emission pattern of that mode.

For the case of including phonons.
"""
function fig_GnmEigenModes(zs, eigen_σ, eigen_diagBα, ks, eigen_σ_FT, eigen_diagBα_FT, z_range, x_range, intensity, ρf, array, eigval, κ, titl)
    # Prepare colors
    colors = distinguishable_colors(8, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Start figure 
    fig = Figure(size=(800, 700))
    
    # Make title and axis
    Label(fig[1, 1:3], titl, tellwidth=false)
    Label(fig[2, 2:3], latexstring(L"Eigval. $ = " * format_Complex_to_String(eigval) * L"$"), tellwidth=false)
    Label(fig[3, 2], L"Real space$$", tellwidth=false)
    Label(fig[3, 3], L"Momentum space$$", tellwidth=false)
    Label(fig[4, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[5, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[4, 2])
    ax2 = Axis(fig[4, 3])
    ax3 = Axis(fig[5, 2], xlabel=L"$ z/λ_{a} $")
    ax4 = Axis(fig[5, 3], xlabel=L"$ λ_{a}k_z $")
    ax5 = Axis(fig[6, 2:3], aspect=DataAspect(), xlabel=L"$ z/λ_{a} $", ylabel=L"$ x/λ_{a} $")
    
    # Plot real space 
    lines!(ax1, zs, abs.(eigen_σ), color=colors[1], label=L"$ ⟨σ_{n}⟩ $")
    for α in 1:3 lines!(ax1, zs, abs.(eigen_diagBα[α]), linestyle=:dash, color=colors[1 + α], label=L"$ ⟨b_{%$(α)n}σ_{n}⟩ $") end
    lines!(ax3, zs, angle.(eigen_σ), color=colors[1])
    for α in 1:3 lines!(ax3, zs, angle.(eigen_diagBα[α]), linestyle=:dash, color=colors[1 + α]) end
    
    # Plot  k-space 
    lines!(ax2, ks, abs.(eigen_σ_FT)  , color=colors[5])
    for α in 1:3 lines!(ax2, ks, abs.(eigen_diagBα_FT[α]), linestyle=:dash, color=colors[5 + α]) end
    lines!(ax4, ks, angle.(eigen_σ_FT), color=colors[5], label=L"FT$ [⟨σ_{n}⟩](k) $")
    for α in 1:3 lines!(ax4, ks, angle.(eigen_diagBα_FT[α]), linestyle=:dash, color=colors[5 + α], label=L"FT$ [⟨b_{%$(α)n}σ_{n}⟩](k) $") end
    
    # Mark the light cone and position of propagation constant
    vlines!(ax2, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
    vlines!(ax2, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    
    # Plot the E-field intensity
    hm = heatmap!(ax5, z_range, x_range, intensity, colormap=:viridis)
    Colorbar(fig[6, 4], hm)
    
    # Plot a representation of the fiber
    poly!(ax5, [(z_range[1], -ρf), (z_range[end], -ρf), (z_range[end], ρf), (z_range[1], ρf)], color=(:gray, 0.3))
    
    # Plot a representation of the atomic array
    scatter!(ax5, [site[3] for site in array], [site[1] for site in array], color=:black, marker=:circle, markersize=10)
    
    # Finish figure
    axislegend(ax1)
    axislegend(ax2)
    axislegend(ax4)
    display(GLMakie.Screen(), fig)
end


"""
Plot the collective energies of a coupling matrix as a function of the dominant k
in their discrete Fourier transform (i.e. plot the band structure of the coupling matrix)
"""
function fig_eigenEnergies_vs_k(dominant_ks, collΔ, collΓ, weights_abs, κ, titl)
    # Start figure 
    fig = Figure(size=(800, 900))
    
    # Make title and axis
    Label(fig[1, 1], titl, tellwidth=false)
    Label(fig[end+1, 1], L"Collective energies, $ Δ_\text{coll}/γ_{a} $", tellwidth=false)
    ax1 = Axis(fig[end+1, 1])
    Label(fig[end+1, 1], L"Collective decay rates, $ Γ_\text{coll}/γ_{a} $", tellwidth=false)
    ax2 = Axis(fig[end+1, 1]) #, yscale=log10)
    Label(fig[end+1, 1], L"Resonance weights, $ |w| $", tellwidth=false)
    ax3 = Axis(fig[end+1, 1], 
               xlabel=L"$ λ_{a}k_z $",
               yscale=log10)
    
    # Plot
    scatter!(ax1, dominant_ks, collΔ, color=:blue)
    scatter!(ax2, dominant_ks, collΓ, color=:red)
    scatter!(ax3, dominant_ks, weights_abs, color=:black)
    
    # Mark the light cone and position of propagation constant
    for ax in [ax1, ax2, ax3]
        vlines!(ax, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
        vlines!(ax, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    end
    
    # Finish figure
    axislegend(ax1, position=:ct)
    display(GLMakie.Screen(), fig)
end


"""
Plot and compare the collective energies of a number of coupling matrices as a function of the dominant k
in their discrete Fourier transform (i.e. plot the band structure of the coupling matrices)
"""
function fig_compareEigenEnergies_vs_k(dominant_kss, collΔs, collΓs, kz_range, collΔ_inf, collΓ_inf, Ns, κ, titl)
    # Prepare colors
    colors = distinguishable_colors(length(Ns), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
    
    # Start figure 
    fig = Figure(size=(800, 900))
    
    # Make title and axis
    Label(fig[1, 1], titl, tellwidth=false)
    Label(fig[end+1, 1], L"Collective energies, $ Δ_\text{coll}/γ_{a} $", tellwidth=false)
    ax1 = Axis(fig[end+1, 1])
    Label(fig[end+1, 1], L"Collective decay rates, $ Γ_\text{coll}/γ_{a} $", tellwidth=false)
    ax2 = Axis(fig[end+1, 1],
               xlabel=L"$ λ_{a}k_z $")
    
    # Plot the infinite case energies
    lines!(ax1, kz_range, collΔ_inf, color=:black)
    lines!(ax2, kz_range, collΓ_inf, color=:black)
    
    # Plot the finite case energies
    for (i, (dominant_ks, collΔ, collΓ, N)) in enumerate(zip(dominant_kss, collΔs, collΓs, Ns))
        scatter!(ax1, dominant_ks, collΔ, color=colors[i], label=L"N = %$(N)")
        scatter!(ax2, dominant_ks, collΓ, color=colors[i])
    end
    
    # Mark the light cone and position of propagation constant
    vlines!(ax1, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1)
    vlines!(ax1, [κ], color=:red, linestyle=:dash, linewidth=1)
    vlines!(ax2, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
    vlines!(ax2, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    
    # Finish figure
    axislegend.([ax1, ax2], position=:lt)
    display(GLMakie.Screen(), fig)
end


"""
Plot magnitude and phase of transmission amplitude as a function of detuning
and mark the position of the eigvals of Gnm
"""
function fig_loss_withGnmeigenEnergies(Δ_range, L, resonances_abs, collΔ, collΓ, weights_abs, titl)
    # Start figure 
    fig = Figure(size=(800, 900))
    
    # Plot the loss
    Label(fig[1, 1], titl, tellwidth=false)
    Label(fig[end+1, 1], L"Loss coefficient, individual resonances superimposed $$", tellwidth=false)
    ax1 = Axis(fig[end+1, 1], limits=(extrema(Δ_range)..., 0, nothing), 
    # ax1 = Axis(fig[end+1, 1], limits=(extrema(Δ_range)..., nothing, nothing), 
               xlabel=L"$ Δ/γ_{a} $",
               ylabel=L"$ 1 - |t|^2 $")
               lines!(ax1, Δ_range, L, color=:blue, label=L"Loss, $ 1 - |t|^2 $")
    
    # Plot the resonances
    for resonance in resonances_abs
        lines!(ax1, Δ_range, resonance, color=:skyblue, label=L"Resonances, $ \left|\frac{γ_{a}w}{(Δ - Δ_\text{coll} + iΓ_\text{coll}/2)}\right| $")
    end
    axislegend(position=:lt, unique=true)
            
    # Mark the position of the resonances
    vlines!(ax1, collΔ, linewidth=0.2, color=:purple, label=false)
    
    # Plot the Gnm eigenmode decay rates as a function of their energy
    Label(fig[end+1, 1], L"Collective decay rates, $ Γ_\text{coll}/γ_{a} $", tellwidth=false)
    ax2 = Axis(fig[end+1, 1], limits=(extrema(Δ_range), nothing), 
               xlabel=L"$ Δ_\text{coll}/γ_{a} $",
               ylabel=L"$ Γ_\text{coll}(Δ_\text{coll})/γ_{a} $", 
               yscale=log10)
    scatter!(ax2, collΔ, collΓ, color=:purple)
    
    # Plot the weights of the resonances
    Label(fig[end+1, 1], L"Resonance weights, $ |w| $", tellwidth=false)
    ax3 = Axis(fig[end+1, 1], limits=(extrema(Δ_range), nothing), 
               xlabel=L"$ Δ_\text{coll}/γ_{a} $",
               ylabel=L"$ |w(Δ_\text{coll})| $",
               yscale=log10)
    scatter!(ax3, collΔ, weights_abs, color=:black)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


"""
Plot the magnitude and phase of the coherences of a state, as well as it discrete
Fourier transform and its emission pattern
"""
function fig_state(rs, v, ks, vFT, dom_ks, collΓ_gm, collΓ_rm, v_eigenModes, z_range, x_range, intensity, ρf, array, κ, titl)
    # Start figure 
    fig = Figure(size=(1600, 900))
    
    # Make title and axis
    Label(fig[1, 1:3], titl, tellwidth=false)
    Label(fig[2, 2], L"Real space$$", tellwidth=false)
    Label(fig[2, 3], L"Momentum space$$", tellwidth=false)
    Label(fig[2, 4], L"In terms of eigenmodes, vs. dominant $ k_{z} $", tellwidth=false)
    Label(fig[2, 5], L"In terms of eigenmodes, vs. $ \tilde{Γ}_{i, gm} $", tellwidth=false)
    Label(fig[2, 6], L"In terms of eigenmodes, vs. $ \tilde{Γ}_{i, rm} $", tellwidth=false)
    Label(fig[3, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[4, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[3, 2])
    ax2 = Axis(fig[3, 3])
    ax3 = Axis(fig[4, 2], xlabel=L"$ z/λ_{a} $")
    ax4 = Axis(fig[4, 3], xlabel=L"$ λ_{a}k_z $")
    ax5 = Axis(fig[5, 2:3][1, 1], aspect=DataAspect(), xlabel=L"$ z/λ_{a} $", ylabel=L"$ x/λ_{a} $")
    ax6 = Axis(fig[3, 4])
    ax7 = Axis(fig[4, 4], xlabel=L"dominant $ λ_{a}k_{z} $")
    ax8 = Axis(fig[3, 5])
    ax9 = Axis(fig[4, 5], xlabel=L"$ \tilde{Γ}_{i. gm}/γ_{a} $")
    ax10 = Axis(fig[3, 6])
    ax11 = Axis(fig[4, 6], xlabel=L"$ \tilde{Γ}_{i, rm}/γ_{a} $")
    
    # Plot the real space function
    lines!(ax1, rs, abs.(v)  , color=:blue)
    lines!(ax3, rs, angle.(v), color=:blue)
    
    # Plot the k-space function
    lines!(ax2, ks, abs.(vFT)  , color=:red)
    lines!(ax4, ks, angle.(vFT), color=:red)
    
    # Mark the light cone and position of propagation constant
    vlines!(ax2, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
    vlines!(ax2, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    axislegend(ax2, position=:ct)
    
    # Plot the E-field intensity
    hm = heatmap!(ax5, z_range, x_range, intensity, colormap=:viridis)#, colorrange = (0, 0.05))
    Colorbar(fig[5, 2:3][1, 2], hm)
    
    # Plot a representation of the fiber
    poly!(ax5, [(z_range[1], -ρf), (z_range[end], -ρf), (z_range[end], ρf), (z_range[1], ρf)], color=(:gray, 0.3))
    
    # Plot a representation of the atomic array
    scatter!(ax5, [site[3] for site in array], [site[1] for site in array], color=:black, marker=:circle, markersize=10)
    
    # Plot the eigenmodes function vs dominant kz
    sortInds = sortperm(dom_ks)
    lines!(ax6, dom_ks[sortInds], abs.(v_eigenModes)[sortInds]  , color=:purple)
    lines!(ax7, dom_ks[sortInds], angle.(v_eigenModes)[sortInds], color=:purple)
    vlines!(ax6, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1)
    vlines!(ax6, [κ], color=:red, linestyle=:dash, linewidth=1)
    
    # Plot the eigenmodes function vs guided decay rate
    sortInds = sortperm(collΓ_gm)
    lines!(ax8, collΓ_gm[sortInds], abs.(v_eigenModes)[sortInds]  , color=:black)
    lines!(ax9, collΓ_gm[sortInds], angle.(v_eigenModes)[sortInds], color=:black)
    
    # Plot the eigenmodes function vs radiation decay rate
    sortInds = sortperm(collΓ_rm)
    lines!(ax10, collΓ_rm[sortInds], abs.(v_eigenModes)[sortInds]  , color=:black)
    lines!(ax11, collΓ_rm[sortInds], angle.(v_eigenModes)[sortInds], color=:black)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


# save("C:\\Users\\Simon\\Downloads\\FILENAME.png", fig)
