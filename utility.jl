
"""
Calculate the polylogarithm of order s evaluated at z.

The precision of the calculation is determined by abstol.
"""
function polylogarithm(z, s, abstol=1e-3)
    sum = 0.0 + 0.0im
    n = 1
    summand = 10.0*abstol + 0.0im
    while abs(summand) > abstol
        summand = z^n/n^s
        sum += summand
        n += 1
    end
    return sum
end


"""
The Heaviside step function evaluated at x.

Returns 1/2 at x = 0.
"""
function heaviside(x)
    if x > 0  return 1
    elseif x == 0  return 1/2
    else return 0
    end
end


"""
Take a matrix M and return its diagonal as a vector
"""
function di(M)
    return diag(M)
end


"""
Take a vector v and return a matrix whose diagonal is v, and all other entries are zero
"""
function Di(v)
    return Diagonal(v)
end


"""
Take a matrix M and and return a matrix with the same diagonal, but all other entries are zero
"""
function Diag(M)
    return Diagonal(diag(M))
end


"""
Returns the length of the x-vector for a given N
"""
function lengthx_fromN(N)
    return 2*(N + 3*N^2)
end


"""
Returns the value of N for a given x-vector
"""
function N_fromx(x)
     #if length(x) = 2(N + 3N^2), then 6*length(x) + 1 = (6N + 1)^2, and N is equal to the following
    return Int((sqrt(6*length(x) + 1) - 1)/6)
end


"""
Prepare an empty x-vector (i.e. all entries are zero)
"""
function empty_xVector(N)
    return zeros(lengthx_fromN(N))
end


"""
Prepare σBα vectors (i.e. all entries are zero)
"""
function empty_σBαVectors(N)
    return zeros(ComplexF64, N), [zeros(ComplexF64, N, N) for α in 1:3]
end


"""
Pack the σ and Bα entries into the x vector to facilitate NonlinearSolve

Assumes that σ is an N_vector, Bα is a 3-vector consisting of NxN-matrices, and x is a 2(N + 3N^2)-vector
"""
function pack_σBαIntox!(σ, Bα, x)
    N = length(σ)
    x[1     : N]   .= real.(σ)
    x[1 + N : 2*N] .= imag.(σ)
    for α in 1:3
        x[1 + 2*N +         (α - 1)*N^2 : 2*N         + α*N^2] .= real.(flatten(Bα[α]))
        x[1 + 2*N + 3*N^2 + (α - 1)*N^2 : 2*N + 3*N^2 + α*N^2] .= imag.(flatten(Bα[α]))
    end
end


function pack_σBαIntox(σ, Bα)
    N = length(σ)
    x = empty_xVector(N)
    pack_σBαIntox!(σ, Bα, x)
    return x
end


"""
Unpack the σ and Bα entries from the x vector

Assumes that σ is an N_vector, Bα is a 3-vector consisting of NxN-matrices, and x is a 2(N + 3N^2)-vector
"""
function unpack_σBαFromx!(σ, Bα, x)
    N = length(σ)
    σ .= x[1 : N] + 1im*x[1 + N : 2*N]
    for α in 1:3
        Bα[α] .= reshape(x[1 + 2*N + (α - 1)*N^2 : 2*N + α*N^2] + 1im*x[1 + 2*N + 3*N^2 + (α - 1)*N^2 : 2*N + 3*N^2 + α*N^2], (N, N))
    end
end


function unpack_σBαFromx(x)
    N = N_fromx(x)
    σ, Bα = empty_σBαVectors(N)
    unpack_σBαFromx!(σ, Bα, x)
    return σ, Bα
end


"""
Return a vector that consists of the concatenated columns of A
"""
function flatten(A)
    return reduce(vcat, A)
end


"""
Take a vector in Cartesian coordinates (x, y, z) and return it in cylindrical coordinates (ρ, ϕ, z)
"""
function cylCoordinates(r)
    return sqrt(r[1]^2 + r[2]^2), atan(r[2]/r[1]), r[3]
end


"""
Take a position vector in Cartesian coordinates (x, y, z) and return the cylindrical unit vectors at that point (ρ_unit, ϕ_unit, z_unit)
"""
function cylUnitVectors(r)
    ρ_unit = [r[1], r[2], 0]/sqrt(r[1]^2 + r[2]^2)
    ϕ_unit = [-r[2], r[1], 0]/sqrt(r[1]^2 + r[2]^2)
    z_unit = [0, 0, 1]
    # ρ_unit = SVector{3}([r[1], r[2], 0]/sqrt(r[1]^2 + r[2]^2))
    # ϕ_unit = SVector{3}([-r[2], r[1], 0]/sqrt(r[1]^2 + r[2]^2))
    # z_unit = SVector{3}([0, 0, 1])
    return ρ_unit, ϕ_unit, z_unit
end


"""
Derivative of the Bessel functions of the first kind, J, of order n evaluated at x
"""
function dbesselj(derivOrder, n, x)
    return sum([(-1)^i*binomial(derivOrder, i)*besselj(n - (derivOrder - 2*i), x) for i in 0:derivOrder])/2^derivOrder
end


"""
Derivative of the Hankel functions H of the j'th kind of order n evaluated at x
"""
function dbesselh(derivOrder, n, j, x)
    return sum([(-1)^i*binomial(derivOrder, i)*besselh(n - (derivOrder - 2*i), j, x) for i in 0:derivOrder])/2^derivOrder
end


"""
Derivative of the modified Bessel functions of the second kind, K, of order n evaluated at x
"""
function dbesselk(derivOrder, n, x)
    return (-1)^derivOrder*sum([binomial(derivOrder, i)*besselk(n - (derivOrder - 2*i), x) for i in 0:derivOrder])/2^derivOrder
end


"""
Derivatives of the dyad rhat*rhat'
"""
function drhatrhat(rvec, derivOrder, α)
    # Set up length and direction of rvec, as well as the coordinate with respect to which we are deriving and its corresponding unit vector
    r     = norm(rvec)
    rhat  = rvec/r
    αcoor = rvec[α]
    αhat  = zeros(3); αhat[α] = 1
    
    # Calculate derivatives of rvec (starting with the zeroth order derivative)
    drvec = [rvec]
    if derivOrder >= 1 push!(drvec, (αhat - αcoor/r*rhat)/r) end
    if derivOrder >= 2 push!(drvec, ((3*αcoor^2/r^2 - 1)*rhat - 2*αcoor/r*αhat)/r^2) end
    if derivOrder >= 3 throw(ArgumentError("drr is not implemented for derivOrder > 2")) end
    
    # Put together the derivative of rvec*rvec'
    return sum([binomial(derivOrder, i)*drvec[i+1]*drvec[derivOrder - i + 1]' for i in 0:derivOrder])
end


"""
Round off to a standard number (4) of significant digits
"""
function ro(x, sigdigits=4)
    return round(x, sigdigits=sigdigits)
end


"""
arrayDescription for the usual 1D chain along the fiber, as given by get_array()
"""
function standardArrayDescription(N, ρa, a)
    return "stA_N_$(N)_rhoa_$(ρa)_a_$(a)"
end


"""
A nice way to format complex numbers in strings
"""
function format_Complex_to_String(z)
    # Return a zero if z = 0
    if z == 0 return "0" end
    
    # Format the real and imaginary part
    r_string = real(z) != 0 ? @sprintf("%.3f", real(z))        : ""
    i_string = imag(z) != 0 ? @sprintf("%.3f", imag(z)) * "im" : ""
    
    # If the imaginary part is negative is has a sign, otherwise we add a plus, but only if the real part is nonzero
    if real(z) != 0 && imag(z) > 0 i_string = "+" * i_string end
    
    return r_string * i_string
end


"""
Returns a matrix whose rows are given by the entries of the given vector
"""
function vectorOfRows2Matrix(x)
    return Matrix(reduce(hcat, x)')
end


"""
Print out the details of the System Parameters
"""
function printSP(SP)
    println("--- System Parameters ---")
    
    println("The bare parameters from the article 'Magic-wavelength nanofiber-based two-color dipole trap with sub-λ/2 spacing'")
    for key in [:λ0, :ω0, :ρf0, :n0, :ρa0, :a0, :να0, :νR0, :ηα0]
        println(key, ": ", SP[key])
    end
    println("")
    
    show(SP.fiber)
    println("")
    
    println("Specs for scan of guided mode propagation constant")
    for key in [:ω_specs, :ρf_specs, :n_specs]
        println(key, ": ", SP[key])
    end
    println("")
    
    println("Specs for scan of time evolution and steady state")
    println("Δ_specs: ", SP.Δ_specs)
    println("")
    
    println("Time spand and maximum time step allowed in time evolution")
    println("tspan: ", SP.tspan)
    println("dtmax: ", SP.dtmax)
    println("")
    
    println("Description of the atomic array")
    println(SP.arrayDescription)
    println("")
    
    println("Description of the initial state")
    println(SP.initialStateDescription)
    println("")
    
    
    println("Trap frequencies, Lamb-Dicke parameters, and atomic dipole moment")
    println("να: ", SP.να)
    println("ηα: ", SP.ηα)
    if typeof(SP.d) == String println("d: ", SP.d)
    else println("d:  [", join(format_Complex_to_String.(SP.d), ", "), "]") end
    println("")
    
    println("---  ---")
end