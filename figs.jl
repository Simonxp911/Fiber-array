
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
               xlabel=L"$ \omega $, [nm$^{-1}$]", 
               ylabel=L"[nm$^{-1}$]")
    
    # Plot the propagation constant, and the lines y = ω and y = nω
    lines!(ax1, ω_range, κ          , label=L"$ y = κ(\omega) $", color=:black)
    lines!(ax1, ω_range, n * ω_range, label=L"$ y = n\omega $"       , color=:red)
    lines!(ax1, ω_range, ω_range    , label=L"$ y = \omega $"        , color=:blue)
    
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
                    ylabel=L"$ \omega $, [nm$^{-1}$]")
    
    # Plot the momenta
    lines!(ω_range, h, label=L"$ h = \sqrt{n^2\omega^2 - \kappa^2} $", color=:blue)
    lines!(ω_range, q, label=L"$ q = \sqrt{\kappa^2 - \omega^2} $", color=:red)
    
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
    
    # Make title and axis
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
                    xlabel=L"$ γt $")
    
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
               xlabel=L"$ γt $")
    Label(fig[1, 2], L"$ B_α $ trajectories", tellwidth=false)
    ax2 = Axis(fig[2, 2], limits=(extrema(times), nothing), 
               xlabel=L"$ γt $")
    
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
    Label(fig[1, :], titl, tellwidth=false)
    Label(fig[2, 1], "Transmission coefficient", tellwidth=false)
    ax1 = Axis(fig[3, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/\gamma $")
    Label(fig[2, 2], "Transmission phase", tellwidth=false)
    ax2 = Axis(fig[3, 2], limits=(extrema(Δ_range)..., -π, π),
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=L"$ \Delta/\gamma $")
               
    # Plot magnitude squared and the phase of the transmission 
    lines!(ax1, Δ_range, T    , label=L"$ |t|^2 $" , color=:blue)
    lines!(ax2, Δ_range, phase, label=L"arg$ (t) $", color=:red)
    
    # Finish figure
    axislegend.([ax1, ax2])
    display(GLMakie.Screen(), fig)
end


"""
Plot mean magnitude and phase of transmission amplitude as a function of detuning

with bands given by the standard deviation
"""
function fig_classDisorder_transmission_vs_Δ(Δ_range, T_means, T_stds, phase_means, phase_stds)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1], "Transmission coefficient", tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/\gamma $")
    Label(fig[1, 2], "Transmission phase", tellwidth=false)
    ax2 = Axis(fig[2, 2], limits=(extrema(Δ_range)..., -π, π), 
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=L"$ \Delta/\gamma $")
               
    # Plot magnitude squared and the phase of the transmission with bands for standard deviations
    lines!(ax1, Δ_range, T_means, label=L"$ |t|^2 $" , color=:blue)
    band!( ax1, Δ_range, T_means + T_stds, T_means - T_stds , color=(:blue, 0.35))
    lines!(ax2, Δ_range, phase_means, label=L"arg$ (t) $", color=:red)
    band!( ax2, Δ_range, phase_means + phase_stds, phase_means - phase_stds , color=(:red, 0.35))
    
    # Finish figure
    axislegend.([ax1, ax2])
    display(GLMakie.Screen(), fig)
end


"""
Plot the intensity (norm-squared) of the radiated E-field around the fiber
"""
function fig_radiation_Efield(z_range, x_range, intensity, ρf, array)
    # Start figure 
    fig = Figure(size=(800, 600))
    
    # Make title and axis
    Label(fig[1, 1], L"$ I/(\gamma/\lambda^3) $", tellwidth=false)
    Axis(fig[2, 1], limits=(extrema(z_range), extrema(x_range)), 
                    xlabel=L"$ z/\lambda $", 
                    ylabel=L"$ x/\lambda $", 
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
    Label(fig[1, 2], "Real space")
    Label(fig[1, 3], "Momentum space")
    Label(fig[2, 1], "Magnitude", rotation = pi/2)
    Label(fig[3, 1], "Phase", rotation = pi/2)
    ax1 = Axis(fig[2, 2])
    ax2 = Axis(fig[2, 3])
    ax3 = Axis(fig[3, 2], xlabel=L"$ z/λ $")
    ax4 = Axis(fig[3, 3], xlabel=L"$ λk_z $")
    
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
Plot the eigenvectors of a coupling matrix in real space, as well as its Fourier transform.
Furthermore, plot the emission pattern of that mode.
"""
function fig_GnmEigenModes(rs, v, ks, vFT, z_range, x_range, intensity, ρf, array, eigval, κ, titl)
    # Start figure 
    fig = Figure(size=(800, 700))
    
    # Make title and axis
    Label(fig[1, :], titl, tellwidth=false)
    Label(fig[2, 2:3], latexstring(L"Eigval. $ = " * format_Complex_to_String(eigval) * L"$"), tellwidth=false)
    Label(fig[3, 2], L"Real space$$", tellwidth=false)
    Label(fig[3, 3], L"Momentum space$$", tellwidth=false)
    Label(fig[4, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[5, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[4, 2])
    ax2 = Axis(fig[4, 3])
    ax3 = Axis(fig[5, 2], xlabel=L"$ z/λ $")
    ax4 = Axis(fig[5, 3], xlabel=L"$ λk_z $")
    ax5 = Axis(fig[6, 2:3], aspect=DataAspect(), xlabel=L"$ z/λ $", ylabel=L"$ x/λ $")
    
    # Plot the real space function
    lines!(ax1, rs, abs.(v)  , color=:blue)
    lines!(ax3, rs, angle.(v), color=:blue)
    
    # Plot the k-space function
    lines!(ax2, ks, abs.(vFT)  , color=:red)
    lines!(ax4, ks, angle.(vFT), color=:red)
    
    # Mark the light cone and position of propagation constant
    vlines!(ax2, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
    vlines!(ax2, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    axislegend(ax2)
    
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
Plot the collective energies of a coupling matrix as a function of the dominant k
in their discrete Fourier transform (i.e. plot the band structure of the coupling matrix)
"""
function fig_eigenEnergies_vs_k(dominant_ks, collΔ, collΓ, weights_abs, κ, titl)
    # Start figure 
    fig = Figure(size=(800, 900))
    
    # Make title and axis
    Label(fig[1, :], titl, tellwidth=false)
    Label(fig[end+1, 1], L"Collective energies, $ Δ_\text{coll}/γ $", tellwidth=false)
    ax1 = Axis(fig[end+1, 1])
    Label(fig[end+1, 1], L"Collective decay rates, $ Γ_\text{coll}/γ $", tellwidth=false)
    ax2 = Axis(fig[end+1, 1])
    Label(fig[end+1, 1], L"Resonance weights, $ |w| $", tellwidth=false)
    ax3 = Axis(fig[end+1, 1], 
               xlabel=L"$ λk_z $",
               yscale=log10)
    
    # Plot the real space function
    scatter!(ax1, dominant_ks, collΔ, color=:blue)
    scatter!(ax2, dominant_ks, collΓ, color=:red)
    scatter!(ax3, dominant_ks, weights_abs, color=:black)
    
    # Mark the light cone and position of propagation constant
    for ax in [ax1, ax2, ax3]
        vlines!(ax, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
        vlines!(ax, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    end
    
    # Finish figure
    axislegend(ax1)
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
    Label(fig[1, :], titl, tellwidth=false)
    Label(fig[end+1, 1], L"Loss coefficient, individual resonances superimposed $$", tellwidth=false)
    ax1 = Axis(fig[end+1, 1], limits=(extrema(Δ_range)..., 0, nothing), 
               xlabel=L"$ Δ/γ $",
               ylabel=L"$ 1 - |t|^2 $")
               lines!(ax1, Δ_range, L, color=:blue, label=L"Loss, $ 1 - |t|^2 $")
    
    # Plot the resonances
    for resonance in resonances_abs
        lines!(ax1, Δ_range, resonance, color=:skyblue, label=L"Resonances, $ \left|\frac{γw}{(Δ - Δ_\text{coll} + iΓ_\text{coll}/2)}\right| $")
    end
    axislegend(position=:lt, unique=true)
            
    # Mark the position of the resonances
    vlines!(ax1, collΔ, linewidth=0.2, color=:purple, label=false)
    
    # Plot the Gnm eigenmode decay rates as a function of their energy
    Label(fig[end+1, 1], L"Collective decay rates, $ Γ_\text{coll}/γ $", tellwidth=false)
    ax2 = Axis(fig[end+1, 1], limits=(extrema(Δ_range), nothing), 
               xlabel=L"$ Δ_\text{coll}/γ $",
               ylabel=L"$ Γ_\text{coll}(Δ_\text{coll})/γ $", 
               yscale=log10)
    scatter!(ax2, collΔ, collΓ, color=:purple)
    
    # Plot the weights of the resonances
    Label(fig[end+1, 1], L"Resonance weights, $ |w| $", tellwidth=false)
    ax3 = Axis(fig[end+1, 1], limits=(extrema(Δ_range), nothing), 
               xlabel=L"$ Δ_\text{coll}/γ $",
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
function fig_state(rs, v, ks, vFT, z_range, x_range, intensity, ρf, array, κ, titl)
    # Start figure 
    fig = Figure(size=(800, 900))
    
    # Make title and axis
    Label(fig[1, :], titl, tellwidth=false)
    Label(fig[2, 2], L"Real space$$", tellwidth=false)
    Label(fig[2, 3], L"Momentum space$$", tellwidth=false)
    Label(fig[3, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[4, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[3, 2])
    ax2 = Axis(fig[3, 3])
    ax3 = Axis(fig[4, 2], xlabel=L"$ z/λ $")
    ax4 = Axis(fig[4, 3], xlabel=L"$ λk_z $")
    ax5 = Axis(fig[5, 2:3], aspect=DataAspect(), xlabel=L"$ z/λ $", ylabel=L"$ x/λ $")
    
    # Plot the real space function
    lines!(ax1, rs, abs.(v)  , color=:blue)
    lines!(ax3, rs, angle.(v), color=:blue)
    
    # Plot the k-space function
    lines!(ax2, ks, abs.(vFT)  , color=:red)
    lines!(ax4, ks, angle.(vFT), color=:red)
    
    # Mark the light cone and position of propagation constant
    vlines!(ax2, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
    vlines!(ax2, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    axislegend(ax2)
    
    # Plot the E-field intensity
    hm = heatmap!(ax5, z_range, x_range, intensity, colormap=:viridis)
    Colorbar(fig[5, 4], hm)
    
    # Plot a representation of the fiber
    poly!(ax5, [(z_range[1], -ρf), (z_range[end], -ρf), (z_range[end], ρf), (z_range[1], ρf)], color=(:gray, 0.3))
    
    # Plot a representation of the atomic array
    scatter!(ax5, [site[3] for site in array], [site[1] for site in array], color=:black, marker=:circle, markersize=10)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end



