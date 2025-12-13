import Pkg; Pkg.add("HomotopyContinuation")
import Pkg; Pkg.add("LinearAlgebra")
import Pkg; Pkg.add("Distributions")
import Pkg; Pkg.add("JLD2")
import Pkg; Pkg.add("Arblib")
import Pkg; Pkg.add("GenericLinearAlgebra")
import Pkg; Pkg.add("KrylovKit")
using HomotopyContinuation
using LinearAlgebra
using Distributions
using JLD2
using Arblib
using GenericLinearAlgebra
using KrylovKit

using Base.Threads
@show(Threads.nthreads())




r = 
9

as = 8

bs = 
4

cs = 
4







@var a[1:r,
1:(as-1)] b[1:r,
1:(bs-1)] c[1:r,
1:cs]

T = 
sum( kron(c[i,:], [b[i,:];1], [a[i,:];1])
for i 
in 1:r)





N = as*bs*cs

dim = r*(as-1+bs-1+cs)

cdim = N - dim
@show (r,as,bs,cs,cdim)
flush(stdout)


@var A[1:N,
1:cdim]

@var B[1:N]

@var t[1:cdim]



vars = [vec([a b c]); t]

params = 
vec([A B])



L = A * t + B

fiber_codim = r*(as-1+bs-1+cs)-(as*bs*cs-cdim)
ll = rand(Uniform(0,1), fiber_codim, r*(as-1+bs-1+cs))


F = 
System(T - L, variables 
= vars, parameters = params)

TF = 
System(T, variables = 
vec([a b c]))



source_point = 
rand(Uniform(0,1), length(vars)
- 
length(t))

x0 = 
vec(source_point)
fiber_slices = ll * vec([a b c]) - ll * vec(source_point)
F = 
System([T - L; fiber_slices], variables 
= vars, parameters = params)



in_sys = 
InterpretedSystem(F) # this takes not so long anymore



@load "monodromy_results.jld2" sols
@load "p0_param.jld2" p0


import HomotopyContinuation: CertificationCache, CertificationParameters, 
                             certify_solution

cache = CertificationCache(in_sys)
cert_params = CertificationParameters(p0; prec=256)

refined_certs = map(((i, sol),) -> begin
    print("\rRefining solution $i / $(length(sols)))    ")
    flush(stdout)
    
    map(j -> midpoint(j), certified_solution_interval_after_krawczyk(
        HomotopyContinuation.extended_prec_certify_solution(
            in_sys,
            sol,
            sol,
            cert_params,
            cache,
            i,
            false;
            max_precision=256,
            extended_certificate=true
        )))
end, enumerate(sols));

@show("refining done!")
@load "unique_indices.jld2" unique_indices
unique_preimages = refined_certs[unique_indices];


function to_bigfloat_certs(refined_certs::Vector{<:AbstractMatrix{Complex{Arf}}})
    n = length(refined_certs)
    out = Vector{Matrix{Complex{BigFloat}}}(undef, n)

    for idx in 1:n
        M = refined_certs[idx]
        m, k = size(M)
        A = Matrix{Complex{BigFloat}}(undef, m, k)

        for i in 1:m, j in 1:k
            z = M[i, j]  # Complex{Arf}
            A[i, j] = Complex{BigFloat}(BigFloat(real(z)), BigFloat(imag(z)))
        end

        out[idx] = A
    end

    return out
end

unique_preimages = to_bigfloat_certs(unique_preimages)

t_coords = map(i -> i[end-1:end], unique_preimages);#extract_t_coords(sols, r, as, bs, cs, cdim)
homogeneous_coords = map(t -> [t; 1], t_coords);




function precompute_powers(coords, max_degree)
    n_points = length(coords)
    n_vars = length(coords[1])
    max_power_idx = ceil(Int, log2(max_degree)) + 1
    
    powers = Vector{Vector{Vector{Complex{BigFloat}}}}(undef, n_points)
    
    println("Precomputing powers...")
    for i in 1:n_points
        if i % 1000 == 0
            println("Point $i / $n_points")
        end
        
        pt = coords[i]
        powers[i] = Vector{Vector{Complex{BigFloat}}}(undef, n_vars)
        
        for j in 1:n_vars
            powers[i][j] = Vector{Complex{BigFloat}}(undef, max_power_idx)
            powers[i][j][1] = pt[j]  # x^1
            
            for k in 2:max_power_idx
                powers[i][j][k] = powers[i][j][k-1]^2  # x^(2^(k-1)) → x^(2^k)
            end
        end
    end
    
    return powers
end

function eval_monomial_fast(powers_pt, exponents)
    result = one(Complex{BigFloat})
    
    for (var_idx, exp) in enumerate(exponents)
        if exp == 0
            continue
        elseif exp == 1
            result *= powers_pt[var_idx][1]
        else
            var_contrib = one(Complex{BigFloat})
            e = exp
            power_idx = 1
            
            while e > 0
                if (e & 1) == 1  # if bit is set
                    var_contrib *= powers_pt[var_idx][power_idx]
                end
                e >>= 1
                power_idx += 1
            end
            
            result *= var_contrib
        end
    end
    
    return result
end

function generate_monomials(degree, n_vars)
    """
    Generate all multi-indices (i₁, ..., iₙ) where i₁ + ... + iₙ = degree
    """
    monomials = []
    
    function generate_recursive(current, remaining_degree, remaining_vars)
        if remaining_vars == 1
            push!(monomials, [current; remaining_degree])
            return
        end
        
        for i in 0:remaining_degree
            generate_recursive([current; i], remaining_degree - i, remaining_vars - 1)
        end
    end
    
    generate_recursive(Int[], degree, n_vars)
    return monomials
end


function build_interpolation_matrix_fast(homog_coords, degree)
    n_points = length(homog_coords)
    n_vars = length(homog_coords[1])
    n_monomials = binomial(degree + n_vars - 1, n_vars - 1)
    
    powers = precompute_powers(homog_coords, degree)
    
    monomials = generate_monomials(degree, n_vars)
    
    M = zeros(Complex{BigFloat}, n_points, n_monomials)
    
    println("Building matrix with fast evaluation...")
    for row in 1:n_points
        if row % 1000 == 0
            println("Row $row / $n_points")
        end
        
        for (col, exponents) in enumerate(monomials)
            M[row, col] = eval_monomial_fast(powers[row], exponents)
        end
        
        row_norm = maximum(abs(M[row, col]) for col in 1:n_monomials)
        if row_norm > 1e-30
            M[row, :] ./= row_norm
        end
    end
    
    return M
end
flush(stdout)
@show("building matrix is about to begin")
flush(stdout)

M = build_interpolation_matrix_fast(homogeneous_coords, 76);

F = qr(M, ColumnNorm())
d = abs.(diag(F.R))
d_min = minimum(d)

println("QR-based min diag ≈ ", d_min)

if d_min > big"1e-40"
    println("Numerically nonsingular ✔ ")
else
    println("Near-singular or singular-ish ✘ ")
end
