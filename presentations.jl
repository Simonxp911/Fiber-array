

"""
Presentation version of fig_transmission_vs_Δ, 
showing magnitude and phase of transmission amplitude as a function of detuning
"""
function fig_presentation_transmission_vs_Δ(Δ_range, T, phase, titl)
    # Start figure 
    fig = Figure(size=(700, 300), fontsize=24)
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/γ_{a} $",
               ylabel=L"$ |t|^2 $")
    ax2 = Axis(fig[2, 2], limits=(extrema(Δ_range)..., -π, π),
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=L"$ \Delta/γ_{a} $",
               ylabel=L"arg$ (t) $")
               
    # Plot magnitude squared and the phase of the transmission 
    lines!(ax1, Δ_range, T    , color=:blue)
    lines!(ax2, Δ_range, phase, color=:red)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation version of fig_radiation_Efield, 
showing the intensity (norm-squared) of the radiated E-field around the fiber
"""
function fig_presentation_radiation_Efield(z_range, x_range, intensity, ρf, array)
    # Start figure 
    zWidth  = maximum(z_range) - minimum(z_range)
    xHeight = maximum(x_range) - minimum(x_range)
    fig = Figure(size=(1000, xHeight/zWidth*1000), fontsize=24)
    
    # Make title and axis
    Axis(fig[1, 1], limits=(extrema(z_range), extrema(x_range)), 
                    xlabel=L"$ z/λ_{a} $", 
                    ylabel=L"$ x/λ_{a} $")
    colsize!(fig.layout, 1, Aspect(1, zWidth/xHeight))
    
    # Plot the E-field intensity
    hm = heatmap!(z_range, x_range, intensity, colormap =:viridis)
    Colorbar(fig[1, 2], hm, label=L"$ I/(γ_{a}/λ_{a}^3) $")
    
    # Plot a representation of the fiber
    poly!([(z_range[1], -ρf), (z_range[end], -ρf), (z_range[end], ρf), (z_range[1], ρf)], color=(:gray, 0.5))
    
    # Plot a representation of the atomic array
    scatter!([site[3] for site in array], [site[1] for site in array], color=:black, marker=:circle, markersize=10)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation version of fig_eigenEnergies_vs_k,
showing the collective energies of a coupling matrix as a function of the dominant k
in their discrete Fourier transform (i.e. plot the band structure of the coupling matrix)
"""
function fig_presentation_eigenEnergies_vs_k(dominant_ks, collΔ, collΓ, weights_abs, κ, titl)
    # Start figure 
    fig = Figure(size=(700, 800), fontsize=24)
    
    # Make title and axis
    Label(fig[1, 1], titl, tellwidth=false)
    ax1 = Axis(fig[end+1, 1],
               ylabel=L"$ \tilde{Δ}_{i}/γ_{a} $")
    ax2 = Axis(fig[end+1, 1],
               ylabel=L"$ \tilde{Γ}_{i}/γ_{a} $",
               yscale=log10, limits=(nothing, nothing, 1e-6, 5))
    ax3 = Axis(fig[end+1, 1], 
               xlabel=L"$ λ_{a}k_z $",
               ylabel=L"$ |w_{i}|Γ_{\text{gm}}/γ_{a} $",
               yscale=log10)
    
    # Plot
    scatter!(ax1, dominant_ks, collΔ, color=:blue)
    scatter!(ax2, dominant_ks, collΓ, color=:red)
    scatter!(ax3, dominant_ks, weights_abs, color=:black)
    
    # Mark the light cone and position of propagation constant
    for ax in [ax1, ax2, ax3]
        vlines!(ax, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
        vlines!(ax, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"Prop. const. $$")
    end
    
    # Finish figure
    axislegend(ax1, position=:ct)
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation version of fig_loss_withGnmeigenEnergies,
showing magnitude and phase of transmission amplitude as a function of detuning
and mark the position of the eigvals of Gnm
"""
function fig_presentation_loss_withGnmeigenEnergies(Δ_range, L, resonances_abs, collΔ, collΓ, weights_abs, titl)
    # Start figure 
    fig = Figure(size=(600, 1000), fontsize=24)
    
    # Plot the loss
    Label(fig[0, 1], titl, tellwidth=false)
    ax1 = Axis(fig[1, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ Δ/γ_{a} $",
               ylabel=L"$ 1 - |t|^2 $")
               lines!(ax1, Δ_range, L, color=:blue, label=L"Loss, $ 1 - |t|^2 $")
    
    # Plot the resonances
    for resonance in resonances_abs
        lines!(ax1, Δ_range, resonance, color=:skyblue, label=L"Peaks, $ \left|\frac{Γ_{\text{gm}}w_{i}}{(Δ - \tilde{Δ}_{i} + i\tilde{Γ}_{i}/2)}\right| $")
    end
    # Legend(fig[0, 1], ax1, unique=true, tellwidth=false)
    # axislegend(position=:lt, unique=true, labelsize=16)
            
    # Mark the position of the resonances
    vlines!(ax1, collΔ, linewidth=0.2, color=:purple, label=false)
    
    # Plot the Gnm eigenmode decay rates as a function of their energy
    ax2 = Axis(fig[end+1, 1], limits=(extrema(Δ_range), nothing), 
               xlabel=L"$ \tilde{Δ}_{i}/γ_{a} $",
               ylabel=L"$ \tilde{Γ}_{i}/γ_{a} $", 
               yscale=log10)
    scatter!(ax2, collΔ, collΓ, color=:red)
    
    # Plot the weights of the resonances
    ax3 = Axis(fig[end+1, 1], limits=(extrema(Δ_range), nothing), 
               xlabel=L"$ \tilde{Δ}_{i}/γ_{a} $",
               ylabel=L"$ |w_{i}|Γ_{\text{gm}}/γ_{a} $",
               yscale=log10)
    scatter!(ax3, collΔ, weights_abs, color=:black)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation figure,
showing the detuning variation across the array
"""
function fig_presentation_Δvari(zs, Δvari)
    # Start figure 
    fig = Figure(size=(700, 200), fontsize=24)
    
    # Make axis
    Axis(fig[1, 1], 
         xlabel=L"$ z_{n}/λ_{a} $",
         ylabel=L"$ Δ $-variation")
               
    # Plot Δvari
    lines!(zs, Δvari, color=:black)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation version of fig_state,
showing the magnitude and phase of the coherences of a state, as well as it discrete
Fourier transform and its emission pattern
"""
function fig_presentation_state(rs, v, ks, vFT, z_range, x_range, intensity, ρf, array, κ, titl)
    # Start figure 
    fig = Figure(size=(700, 600), fontsize=24)
    
    # Make title and axis
    Label(fig[1, 2], L"Real space$$", tellwidth=false)
    Label(fig[1, 3], L"Momentum space$$", tellwidth=false)
    Label(fig[2, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[3, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[2, 2])
    ax2 = Axis(fig[2, 3])
    ax3 = Axis(fig[3, 2], xlabel=L"$ z/λ_{a} $")
    ax4 = Axis(fig[3, 3], xlabel=L"$ λ_{a}k_z $")
    ax5 = Axis(fig[4, 2:3], aspect=DataAspect(), xlabel=L"$ z/λ_{a} $", ylabel=L"$ x/λ_{a} $")
    
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
    hm = heatmap!(ax5, z_range, x_range, intensity, colormap=:viridis)
    Colorbar(fig[4, 4], hm)
    
    # Plot a representation of the fiber
    poly!(ax5, [(z_range[1], -ρf), (z_range[end], -ρf), (z_range[end], ρf), (z_range[1], ρf)], color=(:gray, 0.3))
    
    # Plot a representation of the atomic array
    scatter!(ax5, [site[3] for site in array], [site[1] for site in array], color=:black, marker=:circle, markersize=10)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation of the phonon levels of the atomic array along the fiber
"""
function fig_presentation_arrayIn3D_and_harmonicTrap()
    ρf = 1
    ρa = 2
    a  = 4
    array = [[ρa, 0, n*a] for n in -1:1]
    x_range = range(-ρa, 3*ρa, 100)
    y_range = range(-2*ρa, 2*ρa, 100)
    z_range = range(-a - 1, a + 1, 100)
    
    # Start figure 
    fig = Figure(size=(900, 600), fontsize=24)
    
    # Make title and axis
    zWidth  = maximum(z_range) - minimum(z_range)
    xHeight = maximum(x_range) - minimum(x_range)
    yDepth  = maximum(y_range) - minimum(y_range)
    Axis3(fig[1, 1], limits=(extrema(z_range), extrema(y_range), extrema(x_range)), 
                     yreversed = true,
                     xlabel=L"$ z/λ_{a} $", 
                     ylabel=L"$ y/λ_{a} $", 
                     zlabel=L"$ x/λ_{a} $", 
                     aspect=(zWidth, yDepth, xHeight)./maximum((zWidth, yDepth, xHeight)))
    
    # Plot the atoms
    radius = 0.25
    θs = range(0, π, 20)
    φs = range(0, 2π, 20)
    xSph = radius.*[cos(φ)*sin(θ) for θ in θs, φ in φs]
    ySph = radius.*[sin(φ)*sin(θ) for θ in θs, φ in φs]
    zSph = radius.*[cos(θ) for θ in θs, φ in φs]
    for site in array
        surface!(zSph .+ site[3], ySph .+ site[2], xSph .+ site[1], colormap=:grays)
    end
    
    # Plot the harmonic trap
    x_trap = @. z_range^2 + ρa - 2*radius
    y_trap = zeros(size(z_range))
    lines!(z_range, y_trap, x_trap, color=(:black, 0.9), linewidth=3)
    for shift in 0:3
        lines!([-sqrt(2*radius + shift), sqrt(2*radius + shift)], [0, 0], [ρa + shift, ρa + shift], color=(:black, 0.9), linewidth=2.5)
        text!(-sqrt(2*radius + shift) - 1.5, 0, ρa + shift, text=L"$ | %$(shift) \rangle $", clip_planes=Plane3f[])
    end
    
    # Plot inner state labels
    text!(4, 0, 3.5, text=L"$ | g \rangle $", clip_planes=Plane3f[])
    text!(4, 0, 4.5, text=L"$ | e \rangle $", clip_planes=Plane3f[])
    
    # Plot the fiber
    xCyl = ρf.*[cos(φ) for z in z_range, φ in φs]
    yCyl = ρf.*[sin(φ) for z in z_range, φ in φs]
    zCyl = [z for z in z_range, φ in φs]
    surface!(zCyl, yCyl, xCyl, colormap=(:grays, 0.3))
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation version of fig_compareImperfectArray_transmission_vs_Δ,
showing mean magnitude and phase of transmission amplitude as a function of detuning
with bands given by the standard deviation
for several values of ff and ηα
"""
function fig_presentation_compareImperfectArray_transmission_vs_Δ(Δ_range, T_meanss, T_stdss, phase_meanss, phase_stdss, labels, titl)
    colors = distinguishable_colors(length(labels), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    # colors = distinguishable_colors(5, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    # colors = [colors[5]]
    
    # Start figure 
    fig = Figure(size=(900, 300), fontsize=24)
    # fig = Figure(size=(500, 300), fontsize=12)
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    # Label(fig[1, 1], titl, tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(Δ_range)..., 0, 1), 
               xlabel=L"$ \Delta/γ_{a} $",
               ylabel=L"$ |t|^2 $")
    ax2 = Axis(fig[2, 2], limits=(extrema(Δ_range)..., -π, π),
               yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=L"$ \Delta/γ_{a} $",
               ylabel=L"arg$ (t) $")
    
    # Plot magnitude squared and the phase of the transmission with bands for standard deviations
    for (i, label) in enumerate(labels)
        # lines!(ax1, Δ_range, T_meanss[i], color=colors[i])
        lines!(ax1, Δ_range, T_meanss[i], color=colors[i], label=L"$ %$(label) $")
        band!( ax1, Δ_range, T_meanss[i] + T_stdss[i], T_meanss[i] - T_stdss[i], color=(colors[i], 0.35))
        lines!(ax2, Δ_range, phase_meanss[i], color=colors[i], label=label)
        # lines!(ax2, Δ_range, phase_meanss[i], color=colors[i], label=L"$ %$(label) $")
        band!( ax2, Δ_range, phase_meanss[i] + phase_stdss[i], phase_meanss[i] - phase_stdss[i], color=(colors[i], 0.35))
    end
    
    # Finish figure
    # axislegend(ax2, position=:lb, titlesize=16, labelsize=16)
    # axislegend(ax2, L"$ ff = \dots $", position=:lb, titlesize=16, labelsize=16)
    # axislegend(ax2, L"$ ff = \dots $", position=:lb, titlesize=16, labelsize=16, nbanks=2)
    # axislegend(ax2, L"$ r = \dots $", position=:rt, titlesize=16, labelsize=16)
    # Legend(fig[2, 2], ax1, L"$ ff = \dots $", titlesize=12, labelsize=12, nbanks=2)
    Legend(fig[2, 3], ax1, L"$ ff = \dots $", titlesize=12, labelsize=12, nbanks=2)
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation version of fig_compareImperfectArray_transmission_vs_ffOrηα,
showing mean magnitude and phase of transmission amplitude as a function for a specific detuning
as a function of ff or ηα
with error bars given by the standard deviation
"""
function fig_presentation_compareImperfectArray_transmission_vs_ffOrηα(xs, x_label, T_means, T_stds, phase_means, phase_stds, T_noRadInt, phase_noRadInt, titl)
    # Start figure 
    fig = Figure(size=(700, 300), fontsize=24)
    
    # Make title and axis
    Label(fig[1, 1:2], titl, tellwidth=false)
    # ax1 = Axis(fig[2, 1], limits=(0.0, 0.8, 0, 1), 
    ax1 = Axis(fig[2, 1], limits=(0.0, 0.8, nothing, nothing), 
               xlabel=x_label,
               ylabel=L"$ |t|^2 $")
    # ax2 = Axis(fig[2, 2], limits=(0.0, 0.8, -π, π),
    ax2 = Axis(fig[2, 2], limits=(0.0, 0.8, nothing, nothing),
            #    yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
               xlabel=x_label,
               ylabel=L"arg$ (t) $")
    
    # Add line or scatter for independent decay calculation
    t2 = 0.4507100477120275 #N=100, Δ=-1
    argt = 0.797941863152475
    hlines!(ax1, t2, color=:black, label=L"Ind. decay $$")
    hlines!(ax2, argt, color=:black, label=L"Ind. decay $$")
    # t2s = [0.889, 0.791, 0.704, 0.626, 0.557, 0.495, 0.392, 0.310, 0.245, 0.194]
    # argts = [0.237, 0.469, 0.704, 0.939, 1.173, 1.408, 1.878, 2.346, 2.816, -2.999]
    # scatter!(ax1, xs, T_noRadInt, color=:black, label=L"Ind. decay $$")
    # scatter!(ax2, xs, phase_noRadInt, color=:black, label=L"Ind. decay $$")
    
    # Plot magnitude squared and the phase of the transmission 
    x_n = length(xs)
    for (color, label, T_mean, T_std, phase_mean, phase_std) in zip([:blue, :red], 
                                                                    [L"Ordered$$", L"Random$$"],
                                                                    [T_means[1:x_n], T_means[x_n+1:end]], 
                                                                    [T_stds[1:x_n], T_stds[x_n+1:end]], 
                                                                    [phase_means[1:x_n], phase_means[x_n+1:end]], 
                                                                    [phase_stds[1:x_n], phase_stds[x_n+1:end]])
        errorbars!(ax1, xs, T_mean, T_std, color=color, whiskerwidth=10, label=label)
        errorbars!(ax2, xs, phase_mean, phase_std, color=color, whiskerwidth=10, label=label)
    end
    
    # Finish figure
    # axislegend(ax1, position=:lt, labelsize=16)
    axislegend(ax2, position=:lb, labelsize=16)
    display(GLMakie.Screen(), fig)
    return fig
end


"""
Presentation version of fig_GnmEigenModes, 
showing the eigenvectors of a coupling matrix in real space, as well as their Fourier transform.
Furthermore, plot the emission pattern of that mode.
"""
function fig_presentation_GnmEigenModes(zs, eigen_σ, ks, eigen_σ_FT, z_range, x_range, intensity, ρf, array, eigval, κ, titl)
    # Start figure 
    fig = Figure(size=(800, 700), fontsize=24)
    
    # Make title and axis
    Label(fig[1, 1:3], titl, tellwidth=false)
    Label(fig[2, 2], L"Real space$$", tellwidth=false)
    Label(fig[2, 3], L"Momentum space$$", tellwidth=false)
    Label(fig[3, 1], L"Magnitude$$", rotation = pi/2, tellheight=false)
    Label(fig[4, 1], L"Phase$$", rotation = pi/2, tellheight=false)
    ax1 = Axis(fig[3, 2])
    ax2 = Axis(fig[3, 3])
    ax3 = Axis(fig[4, 2], limits=(nothing, nothing, -π, π),
                          yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
                          xlabel=L"$ z/λ_{a} $")
    ax4 = Axis(fig[4, 3], limits=(nothing, nothing, -π, π),
                          yticks=([-π, -π/2, 0, π/2, π], [L"$ -π $", L"$ -π/2 $", L"$ 0 $", L"$ π/2 $", L"$ π $"]),
                          xlabel=L"$ λ_{a}k_z $")
    ax5 = Axis(fig[5, 2:3][1, 1], aspect=DataAspect(), xlabel=L"$ z/λ_{a} $", ylabel=L"$ x/λ_{a} $")
    
    # Plot real space 
    lines!(ax1, zs, abs.(eigen_σ)  , color=:blue)
    lines!(ax3, zs, angle.(eigen_σ), color=:blue)
    
    # Plot k-space 
    lines!(ax2, ks, abs.(eigen_σ_FT)  , color=:red)
    lines!(ax4, ks, angle.(eigen_σ_FT), color=:red)
    
    # Mark the light cone and position of propagation constant
    vlines!(ax2, [-ωa, ωa], color=:black, linestyle=:dash, linewidth=1, label=L"Light cone$$")
    vlines!(ax2, [κ], color=:red, linestyle=:dash, linewidth=1, label=L"$ \kappa $")
    axislegend(ax2, position=:ct, labelsize=16)
    
    # Plot the E-field intensity
    hm = heatmap!(ax5, z_range, x_range, intensity, colormap=:viridis, colorrange = (0, 0.2))
    Colorbar(fig[5, 2:3][1, 2], hm)
    
    # Plot a representation of the fiber
    poly!(ax5, [(z_range[1], -ρf), (z_range[end], -ρf), (z_range[end], ρf), (z_range[1], ρf)], color=(:gray, 0.3))
    
    # Plot a representation of the atomic array
    scatter!(ax5, [site[3] for site in array], [site[1] for site in array], color=:black, marker=:circle, markersize=10)
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    return fig
end