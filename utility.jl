
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
Prepare an empty x-vector (i.e. all entries are zero), for the case of no phonons
"""
function empty_xVector_noPh(N)
    return zeros(2*N)
end


"""
Prepare an empty x-vector (i.e. all entries are zero)
"""
function empty_xVector(N)
    return zeros(2*(N + 3*N^2))
end


"""
Prepare an empty σ vector (i.e. all entries are zero)
"""
function empty_σVector(N)
    return zeros(ComplexF64, N)
end


"""
Prepare an empty Bα vector (i.e. all entries are zero)
"""
function empty_BαVector(N)
    return [zeros(ComplexF64, N, N) for α in 1:3]
end


"""
Pack the σ entries into the x vector to facilitate NonlinearSolve

Assumes that σ is an N-vector, and x is a 2N-vector
"""
function pack_σIntox!(σ, x)
    N = length(σ)
    x[1     : N]   .= real.(σ)
    x[1 + N : 2*N] .= imag.(σ)
end


function pack_σIntox(σ)
    N = length(σ)
    x = empty_xVector_noPh(N)
    pack_σIntox!(σ, x)
    return x
end


"""
Pack the σ and Bα entries into the x vector to facilitate NonlinearSolve

Assumes that σ is an N-vector, Bα is a 3-vector consisting of NxN-matrices, and x is a 2(N + 3N^2)-vector
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
Unpack the σ entries from the x vector

Assumes that σ is an N-vector, and x is a 2N-vector
"""
function unpack_σFromx!(σ, x)
    N = length(σ)
    σ .= x[1 : N] + 1im*x[1 + N : 2*N]
end


function unpack_σFromx(x)
    N = Int(length(x)/2) #if length(x) = 2N
    σ = empty_σVector(N)
    unpack_σFromx!(σ, x)
    return σ
end


"""
Unpack the σ and Bα entries from the x vector

Assumes that σ is an N-vector, Bα is a 3-vector consisting of NxN-matrices, and x is a 2(N + 3N^2)-vector
"""
function unpack_σBαFromx!(σ, Bα, x)
    N = length(σ)
    σ .= x[1 : N] + 1im*x[1 + N : 2*N]
    for α in 1:3
        Bα[α] .= reshape(x[1 + 2*N + (α - 1)*N^2 : 2*N + α*N^2] + 1im*x[1 + 2*N + 3*N^2 + (α - 1)*N^2 : 2*N + 3*N^2 + α*N^2], (N, N))
    end
end


function unpack_σBαFromx(x)
    N = Int((sqrt(6*length(x) + 1) - 1)/6) #if length(x) = 2(N + 3N^2), then 6*length(x) + 1 = (6N + 1)^2, and N is equal to the following
    σ, Bα = empty_σVector(N), empty_BαVector(N)
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
Remove singleton dimensions from an array
"""
function squeeze(A::AbstractArray)
    singleton_dims = tuple((d for d in 1:ndims(A) if size(A, d) == 1)...)
    return dropdims(A, dims=singleton_dims)
end


"""
Take a vector in Cartesian coordinates (x, y, z) and return it in cylindrical coordinates (ρ, ϕ, z)
"""
function cylCoordinates(r)
    return sqrt(r[1]^2 + r[2]^2), atan(r[2], r[1]), r[3]
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
Derivative of the modified Bessel functions of the second kind, K, of order n evaluated at x
"""
function dbesselk(derivOrder, n, x)
    return (-1)^derivOrder*sum([binomial(derivOrder, i)*besselk(n - (derivOrder - 2*i), x) for i in 0:derivOrder])/2^derivOrder
end


"""
Derivative of the Hankel functions H of the j'th kind of order n evaluated at x
"""
function dbesselh(derivOrder, n, j, x)
    return sum([(-1)^i*binomial(derivOrder, i)*besselh(n - (derivOrder - 2*i), j, x) for i in 0:derivOrder])/2^derivOrder
end


"""
Derivative of the spherical Hankel functions h of the j'th kind of order n evaluated at x
"""
function dbesselsphh(derivOrder, n, j, x)
    if !(derivOrder isa Integer && derivOrder >= 0)
        throw(ArgumentError("The order of differentiation derivOrder (= $derivOrder) take a non-negative integer value in dbesselsphh"))
    end
    
    if derivOrder == 0
        if j == 1
            return sphericalbesselj(n, x) + 1im*sphericalbessely(n, x)
        elseif j == 2
            return sphericalbesselj(n, x) - 1im*sphericalbessely(n, x)
        else
            throw(ArgumentError("The 'kind'-index j (= $j) must take values 1 or 2 in dbesselsphh"))
        end
    else
        return dbesselsphh(derivOrder - 1, n - 1, j, x) - (n + 1)*sum([binomial(derivOrder - 1, i)*(-1)^i*factorial(i)/x^(i + 1)*dbesselsphh(derivOrder - 1 - i, n, j, x) for i in 0:derivOrder - 1])
    end
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
    
    # Calculate derivatives of rhat (starting with the zeroth order derivative)
    drhat = [rhat]
    if derivOrder >= 1 push!(drhat, (αhat - αcoor/r*rhat)/r) end
    if derivOrder >= 2 push!(drhat, ((3*αcoor^2/r^2 - 1)*rhat - 2*αcoor/r*αhat)/r^2) end
    if derivOrder >= 3 throw(ArgumentError("drr is not implemented for derivOrder > 2")) end
    
    # Put together the derivative of rhat*rhat'
    return sum([binomial(derivOrder, i)*drhat[i+1]*drhat[derivOrder - i + 1]' for i in 0:derivOrder])
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
function standardArrayDescription(N, ρa, a, ff, pos_unc)
    if pos_unc isa Number
        return "stA_N_$(N)_rhoa_$(ro(ρa))_a_$(ro(a))_ff_$(ro(ff))_pu_$(ro(pos_unc))"
    else
        return "stA_N_$(N)_rhoa_$(ro(ρa))_a_$(ro(a))_ff_$(ro(ff))_pu_$(join(ro.(pos_unc), ","))"
    end
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
    return Matrix(transpose(reduce(hcat, x)))
end


"""
Returns a matrix whose columns are given by the entries of the given vector
"""
function vectorOfCols2Matrix(x)
    return reduce(hcat, x)
end


"""
Find the eigenvalues and eigenvectors of a matrix. Return these as two sorted vectors.
"""
function eigbasis(A)
    F = eigen(A)
    return F.values, eachcol(F.vectors)
end


"""
Calculate the discrete Fourier transform of a function (represented by an N-vector v) 
on a 1D chain with lattice spacing a and length N = length(v), assumed to be on the form (0:N-1)*a. 
    
Returns an N-vector of the relevant k-values and the transformed function as an N-vector.
"""
function discFourierTransform(v, a, cont_k=false, k_n=300)
    N = length(v)
    if cont_k
        ks = range(0, 2π/a, k_n) .- π/a
    else
        ks = 2π/(a*N)*(0:N-1) .- π/a
    end
    vFT = discFourierTransform.(Ref(v), a, ks)
    return ks, vFT
end


"""
Calculate the discrete Fourier transform of a function (represented by an N-vector v) 
on a 1D chain with lattice spacing a and length N = length(v), assumed to be on the form (0:N-1)*a. 
    
Returns the transformed function as an N-vector.
"""
function discFourierTransform(v, a, k)
    N = length(v)
    return sum( v.*exp.(-1im*k*(0:N-1)*a) )/sqrt(N)
end
