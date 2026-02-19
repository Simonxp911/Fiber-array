

# ================================================
#   Mathematical functions
# ================================================
"""
Calculate the polylogarithm of order s evaluated at z.

The precision of the calculation is determined by abstol.
"""
function polylogarithm(s, z, abstol=1e-4)
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
        return (n*dbesselsphh(derivOrder - 1, n - 1, j, x) - (n + 1)*dbesselsphh(derivOrder - 1, n + 1, j, x))/(2*n + 1)
    end
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
    return ρ_unit, ϕ_unit, z_unit
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


# ================================================
#   Functions pertaining to array manipulation
# ================================================
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
Take a vector of phases that have been calculated modulo 2π and unwrap these. 
There are different approaches: 'minimizeDiff' minimizes differences between
data points (the default); 'forcePositiveSlope' ensures increasing phases; 
'forceNegativeSlope' ensures decreasing phases.
"""
function unwrapPhase(phases, approach="minimizeDiff")
    if approach ∉ ("minimizeDiff", "forcePositiveSlope", "forceNegativeSlope")
        throw(ArgumentError("approach = $approach has not been implemented in unwrapPhase"))
    end
    
    jumps = diff(phases)
    phases_unwrap = deepcopy(phases)
    for (i, jump) in enumerate(jumps)
        if (approach == "minimizeDiff"       && abs(jump) > π) ||
           (approach == "forcePositiveSlope" && jump < 0) ||
           (approach == "forceNegativeSlope" && jump > 0)
            phases_unwrap[i+1:end] .-= sign(jump)*2π
        end
    end
    return phases_unwrap
end


"""
Take a vector of unwrapped (i.e. close to continous) phases and wraps them 
i.e. returns them modulo 2π in the interval -π to π
"""
function wrapPhase(phases)
    return @. mod2pi(phases + π) - π
end


# ================================================
#   Functions pertaining to how the dynamical variables are stored
# ================================================
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
Prepare an empty x-vector (i.e. all entries are zero)
"""
function empty_xVector_noPh_w3l(N)
    return zeros(4*N)
end


"""
Prepare an empty x-vector (i.e. all entries are zero)
"""
function empty_xVector_w3l(N)
    return zeros(4*(N + 3*N^2))
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
    @. x[1     : N]   = real(σ)
    @. x[1 + N : 2*N] = imag(σ)
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
    @. x[1     : N]   = real(σ)
    @. x[1 + N : 2*N] = imag(σ)
    for α in 1:3, i in 1:N, j in 1:N
        x[i + (j - 1)*N + 2*N +         (α - 1)*N^2] = real(Bα[α][i, j])
        x[i + (j - 1)*N + 2*N + 3*N^2 + (α - 1)*N^2] = imag(Bα[α][i, j])
    end
end


function pack_σBαIntox(σ, Bα)
    N = length(σ)
    x = empty_xVector(N)
    pack_σBαIntox!(σ, Bα, x)
    return x
end


function pack_σBαIntox!(σBα, x)
    return pack_σBαIntox!(σBα[1], σBα[2], x)
end


function pack_σBαIntox(σBα)
    return pack_σBαIntox(σBα[1], σBα[2])
end


"""
Pack the σge and σgs entries into the x vector to facilitate NonlinearSolve

Assumes that σge and σgs are N-vectors, and x is a 4N-vector
"""
function pack_σgeσgsIntox!(σge, σgs, x)
    N = length(σge)
    @. x[1      :   N] = real(σge)
    @. x[1 +   N: 2*N] = imag(σge)
    @. x[1 + 2*N: 3*N] = real(σgs)
    @. x[1 + 3*N: 4*N] = imag(σgs)
end


function pack_σgeσgsIntox(σge, σgs)
    N = length(σge)
    x = empty_xVector_noPh_w3l(N)
    pack_σgeσgsIntox!(σge, σgs, x)
    return x
end


"""
Pack the σge, σgs, Bαge, and Bαgs entries into the x vector to facilitate NonlinearSolve

Assumes that σge and σgs are N-vectors, Bαge and Bαgs are 3-vectors consisting of NxN-matrices, and x is a 4(N + 3N^2)-vector
"""
function pack_σgeσgsBαgeBαgsIntox!(σge, σgs, Bαge, Bαgs, x)
    N = length(σge)
    @. x[1      :   N] = real(σge)
    @. x[1 +   N: 2*N] = imag(σge)
    @. x[1 + 2*N: 3*N] = real(σgs)
    @. x[1 + 3*N: 4*N] = imag(σgs)
    for α in 1:3, i in 1:N, j in 1:N
        x[i + (j - 1)*N + 4*N         + (α - 1)*N^2] = real(Bαge[α][i, j])
        x[i + (j - 1)*N + 4*N + 3*N^2 + (α - 1)*N^2] = imag(Bαge[α][i, j])
        x[i + (j - 1)*N + 4*N + 6*N^2 + (α - 1)*N^2] = real(Bαgs[α][i, j])
        x[i + (j - 1)*N + 4*N + 9*N^2 + (α - 1)*N^2] = imag(Bαgs[α][i, j])
    end
end


function pack_σgeσgsBαgeBαgsIntox(σge, σgs, Bαge, Bαgs)
    N = length(σge)
    x = empty_xVector_w3l(N)
    pack_σgeσgsBαgeBαgsIntox!(σge, σgs, Bαge, Bαgs, x)
    return x
end


"""
Unpack the σ entries from the x vector

Assumes that σ is an N-vector, and x is a 2N-vector
"""
function unpack_σFromx!(σ, x)
    N = length(σ)
    for i in 1:N
        σ[i] = x[i] + 1im*x[i + N]
    end
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
    for i in 1:N
        σ[i] = x[i] + 1im*x[i + N]
        for α in 1:3, j in 1:N
            Bα[α][i, j] = x[i + (j - 1)*N + 2*N + (α - 1)*N^2] + 1im*x[i + (j - 1)*N + 2*N + 3*N^2 + (α - 1)*N^2]
        end
    end
end


function unpack_σBαFromx(x)
    N = Int((sqrt(6*length(x) + 1) - 1)/6) #if length(x) = 2(N + 3N^2), then 6*length(x) + 1 = (6N + 1)^2, and N is equal to the following
    σ, Bα = empty_σVector(N), empty_BαVector(N)
    unpack_σBαFromx!(σ, Bα, x)
    return σ, Bα
end


"""
Unpack the σge, σgs entries from the x vector

Assumes that σge and σgs are N-vectors, and x is a 4N-vector
"""
function unpack_σgeσgsFromx!(σge, σgs, x)
    N = length(σge)
    for i in 1:N
        σge[i] = x[i]       + 1im*x[i +   N]
        σgs[i] = x[i + 2*N] + 1im*x[i + 3*N]
    end
end


function unpack_σgeσgsFromx(x)
    N = Int(length(x)/4) #if length(x) = 4N
    σge, σgs = empty_σVector(N), empty_σVector(N)
    unpack_σgeσgsFromx!(σge, σgs, x)
    return σge, σgs
end


"""
Unpack the σge, σgs, Bαge, and Bαgs entries from the x vector

Assumes that σge and σgs are N-vectors, Bαge and Bαgs are 3-vectors consisting of NxN-matrices, and x is a 4(N + 3N^2)-vector
"""
function unpack_σgeσgsBαgeBαgsFromx!(σge, σgs, Bαge, Bαgs, x)
    N = length(σge)
    for i in 1:N
        σge[i] = x[i]       + 1im*x[i   + N]
        σgs[i] = x[i + 2*N] + 1im*x[i + 3*N]
        for α in 1:3, j in 1:N
            Bαge[α][i, j] = x[i + (j - 1)*N + 4*N         + (α - 1)*N^2] + 1im*x[i + (j - 1)*N + 4*N + 3*N^2 + (α - 1)*N^2]
            Bαgs[α][i, j] = x[i + (j - 1)*N + 4*N + 6*N^2 + (α - 1)*N^2] + 1im*x[i + (j - 1)*N + 4*N + 9*N^2 + (α - 1)*N^2]
        end
    end
end


function unpack_σgeσgsBαgeBαgsFromx(x)
    N = Int((sqrt(3*length(x) + 1) - 1)/6) #if length(x) = 4(N + 3N^2), then 3*length(x) + 1 = (6N + 1)^2, and N is equal to the following
    σge, σgs, Bαge, Bαgs = empty_σVector(N), empty_σVector(N), empty_BαVector(N), empty_BαVector(N)
    unpack_σgeσgsBαgeBαgsFromx!(σge, σgs, Bαge, Bαgs, x)
    return σge, σgs, Bαge, Bαgs
end


"""
Pack the σ and Bα entries into the vectorized σBα
"""
function pack_σBαIntoσBαVec(σ, Bα)
    return vcat(σ, vec.(Bα)...)
end


function pack_σBαIntoσBαVec(σBα)
    return vcat(σBα[1], vec.(σBα[2])...)
end


"""
Unpack the σ and Bα entries from the vectorized σBα
"""
function unpack_σBαFromσBαVec(σBαVec)
    N = Int(round(sqrt(length(σBαVec)/3 + 1/36) - 1/6))
    σ = σBαVec[1:N]
    Bα = [reshape(σBαVec[N + (α - 1)*N^2 + 1:N + α*N^2], (N, N)) for α in 1:3]
    return σ, Bα
end


# ================================================
#   Functions pertaining to string labels and descriptions
# ================================================
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
Round off to a standard number (4) of significant digits
"""
function ro(x, sigdigits=4)
    return round(x, sigdigits=sigdigits)
end


"""
arrayDescription for the arrays given by get_array()
"""
function arrayDescript(arrayType, N_sites, ρa, a, ff, pos_unc)
    if pos_unc isa Number
        if pos_unc == 0 return arrayType * "_Nsit_$(N_sites)_rhoa_$(ro(ρa))_a_$(ro(a))_ff_$(ro(ff))"
        else            return arrayType * "_Nsit_$(N_sites)_rhoa_$(ro(ρa))_a_$(ro(a))_ff_$(ro(ff))_pu_$(ro(pos_unc))" end
    else
        return arrayType * "_Nsit_$(N_sites)_rhoa_$(ro(ρa))_a_$(ro(a))_ff_$(ro(ff))_pu_$(join(ro.(pos_unc), ","))"
    end
end


"""
ΔvariDescription for the Δ-variation/potential over the array
"""
function ΔvariDescript(ΔvariDependence, Δvari_args)
    if ΔvariDependence == "flat"
        return ΔvariDependence
    elseif ΔvariDependence ∈ ("Gaussian", "linear", "parabolic")
        return ΔvariDependence * "_" * join(ro.(Δvari_args), ",")
    else 
        throw(ArgumentError("ΔvariDescription has not been implemented for ΔvariDependence = " * ΔvariDependence))
    end
end


"""
Format a number as a string with an explicit sign
"""
function formatNumberWithSign(x)
    if x < 0
        return "$x"
    else
        return "+$x"
    end
end


# ================================================
#   Functions pertaining to Fourier transformation
# ================================================
"""
Calculate the discrete Fourier transform of a function (represented by an N-vector v) 
on a 1D chain with lattice spacing a and length N, assumed to be on the form (0:N-1)*a. 
    
Returns an N-vector of the relevant k-values and the transformed function as an N-vector.
"""
function discFourierTransform(v::Vector, a::Real, cont_k::Bool=false, k_n::Int=300)
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
Calculate the discrete Fourier transform of a function (represented by an (N, N)-matrix M) 
on a 2D square lattice with lattice spacing a and size (N, N), assumed to be on the form (0:N-1)*a ⊗ (0:N-1)*a. 
    
Returns an N-vector of the relevant k-values and the transformed function as an (N, N)-matrix.
"""
function discFourierTransform(M::Matrix, a::Real, cont_k::Bool=false, k_n::Int=300)
    N = size(M)[1]
    if cont_k
        ks = range(0, 2π/a, k_n) .- π/a
    else
        ks = 2π/(a*N)*(0:N-1) .- π/a
    end
    MFT = discFourierTransform.(Ref(M), a, ks, ks')
    return ks, MFT
end


"""
Calculate the discrete Fourier transform of a function (represented by an N-vector v) 
on a 1D chain with lattice spacing a and length N, assumed to be on the form (0:N-1)*a. 
    
Returns the transformed function as an N-vector.
"""
function discFourierTransform(v::Vector, a::Real, k::Real)
    N = length(v)
    return sum( v.*exp.(-1im*k*(0:N-1)*a) )/sqrt(N)
end


"""
Calculate the discrete Fourier transform of a function (represented by an (N, N)-matrix M) 
on a 2D square lattice with lattice spacing a and size (N, N), assumed to be on the form (0:N-1)*a ⊗ (0:N-1)*a. 
    
Returns the transformed function as an (N, N)-matrix.
"""
function discFourierTransform(M::Matrix, a::Real, kx::Real, ky::Real)
    N = size(M)[1]
    return sum( M.*(exp.(-1im*kx*(0:N-1)*a)*exp.(-1im*ky*(0:N-1)'*a)) )/N
end


# ================================================
#   Functions pertaining to interpolation
# ================================================
"""
Interpolate a 1D function, F(z), (whose function values may be scalar or array) from a set of 
equally spaced points Fs_known = F.(zs_known) onto a set of points zs_target
that are contained within the interval covered by zs_known
"""
function interpolation1D_atTargetValues(zs_known, Fs_known, zs_target)
    if !all(diff(zs_known) .≈ zs_known[2] - zs_known[1]) throw(ArgumentError("interpolateCubic_1D assumes the known points of the function to be on a regular 1D grid")) end
    if !all(minimum(zs_known) .< zs_target .< maximum(zs_known)) throw(ArgumentError("interpolateCubic_1D assumes the zs_target to be within the interval covered by zs_known")) end
    
    itp = interpolation1D_asFunction(zs_known, Fs_known)
    return evalNestedFunc.(itp, zs_target)
end


"""
Interpolate a 1D function, F(z), (whose function values may be scalar or array) from a set of 
equally spaced points Fs_known = F.(zs_known) and return the interpolation function 
(which can then be evaluated as needed)
"""
function interpolation1D_asFunction(zs_known, Fs_known)
    if !all(diff(zs_known) .≈ zs_known[2] - zs_known[1]) throw(ArgumentError("interpolateCubic_1D assumes the known points of the function to be on a regular 1D grid")) end
    
    itp = interpolation1D_baseInterpolation(Fs_known) 
    return z -> evalNestedFunc(itp, interpolateIndex(zs_known, z))
end


"""
Get the interpolation object for a set of equally spaced points, Fs_known,
preserving whatever array structure Fs_known may have
"""
function interpolation1D_baseInterpolation(Fs_known)
    if typeof(Fs_known[1]) <: AbstractArray
        itp = Array{Any}(undef, size(Fs_known[1])) 
        for i in eachindex(Fs_known[1])
            Fs_known_i = [F[i] for F in Fs_known]
            itp[i] = interpolation1D_baseInterpolation(Fs_known_i)
        end
    else
        itp = interpolate(Fs_known, BSpline(Cubic(Interpolations.Flat(OnGrid()))))
    end
    return itp
end


"""
Assuming zs_known to be a list of sorted zs, and z_target to be within 
the interval covered by zs_known, return the 
"""
function interpolateIndex(zs_known, z_target)
    if !issorted(zs_known) throw(ArgumentError("interpolateIndex assumes zs_known to be sorted")) end
    if !(minimum(zs_known) <= z_target <= maximum(zs_known)) throw(ArgumentError("interpolateIndex assumes z_target to be within the interval covered by zs_known. Possibly you are trying to evaluate an interpolation function outside the originally given interval.")) end
    
    nearestIndices = sortperm(abs.(zs_known .- z_target))
    leftIndex, rightIndex = sort(nearestIndices[1:2])
    return leftIndex + (z_target - zs_known[leftIndex])/(zs_known[rightIndex] - zs_known[leftIndex])
end


"""
For f a function or f an array of functions, or array of arrays of functions, etc,
evalute that or those functions at x, keeping the array structure
"""
function evalNestedFunc(f, x)
    result = 0.0
    try
        result = f(x)
    catch
        if typeof(f) <: AbstractArray
            result = Array{Any}(undef, size(f)) 
            for i in eachindex(f)
                result[i] = evalNestedFunc(f[i], x)
            end
        else
            throw(ArgumentError("evalNestedFunc assumes f to be either a function or a nested array of functions"))
        end
    end
    return result
end


"""
Construct the interpolation of the Im_Grm_trans. Returns a dcitionary that takes a string
'derivOrder, α' as argument and returns a function of z.
"""
function interpolation1D_Im_Grm_trans(fiber, N_sites, ρa, a, ηα)
    array = get_array("1Dchain", N_sites, ρa, a)
    N = length(array)
    zs = [site[3] for site in array]
    
    derivOrder_α = [((0, 0), 1)]
    if !all(ηα .== 0)
        derivOrders = [(1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]
        αs = 1:3
        derivOrder_α = vcat(derivOrder_α, flatten([(derivOrder, α) for derivOrder in derivOrders, α in αs]))
    end
    
    itp = Dict()
    for (derivOrder, α) in derivOrder_α
        if derivOrder == (1, 1)
            Im_Grm = Im_Grm_trans(fiber, ωa, array[1], array[1], derivOrder, α)
            itp["$derivOrder, $α"] = f(z) = z == 0.0 ? Im_Grm : throw(ArgumentError("The derivOrder = (1, 1) interpolation1D_Im_Grm_trans can only be evaluated at z = 0.0"))
        else
            Im_Grm = fill(zeros(ComplexF64, 3, 3), N)
            for i in 1:N
                Im_Grm[i] = Im_Grm_trans(fiber, ωa, array[i], array[1], derivOrder, α)
            end
            itp["$derivOrder, $α"] = interpolation1D_asFunction(zs, Im_Grm)
        end
    end
    
    return itp
end


# ================================================
#   Functions pertaining to fitting
# ================================================
"""
Given vectors xdata and ydata (i.e. matching points of independent and dependent variables respectively),
and a model function with signature model(x, p), where p is a collection of parameters, returns the optimal
set of parameters pmin that minimize the sum of absolute squared differences between the data and the model,
with p0 as an initial guess. A vector ydataDeviations can be given for the uncertainties of the ydata, which
will be used as weigths in the minimization. The parameters can be constrained by setting lowerConstr and 
upperConstr.
"""
function fitComplexData(xdata, ydata, model, p0; ydataDeviations=ones(size(ydata)))
    sumOfSquares(p) = sum(abs2.( (ydata .- model.(xdata, Ref(p)))./ydataDeviations ))
    res = Optim.optimize(sumOfSquares, p0)
    return Optim.minimizer(res)
end
# function fitComplexData(xdata, ydata, model, p0; ydataDeviations=ones(size(ydata)), lowerConstr=fill(-Inf, length(p0)), upperConstr=fill(Inf, length(p0)))
#     sumOfSquares(p) = sum(abs2.( (ydata .- model.(xdata, Ref(p)))./ydataDeviations ))
#     # res = optimize(sumOfSquares, lowerConstr, upperConstr, p0)
#     # res = optimize(sumOfSquares, lowerConstr, upperConstr, p0, BFGS(linesearch=LineSearches.BackTracking()))
#     res = optimize(sumOfSquares, lowerConstr, upperConstr, p0, LBFGS())
#     # res = optimize(sumOfSquares, p0)
#     return Optim.minimizer(res)
# end