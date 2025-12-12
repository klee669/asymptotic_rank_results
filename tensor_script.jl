import Pkg; Pkg.add("HomotopyContinuation")
import Pkg; Pkg.add("LinearAlgebra")
import Pkg; Pkg.add("Distributions")
using HomotopyContinuation
using LinearAlgebra
using Distributions






r = 
17

as = 14

bs = 
7

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

A0 = 
rand(Uniform(0,1), N, cdim)

t0 = 
rand(Uniform(0,1), cdim)

B0 = 
TF(x0) - A0 * t0

p0 = 
vec([A0 B0])

norm(F([x0; t0], p0))
# test

@show ("sys has been made")
flush(stdout)


in_sys = 
InterpretedSystem(F) # this takes not so long anymore
@show ("int_sys has been made")
flush(stdout)



res = 
monodromy_solve(in_sys, [x0; t0], p0; target_solutions_count 
= 10000)
flush(stdout)


sols = solutions(res)

function shape_solutions(sols)
    sol_list = [];
    for sol in sols
        locator=1;
        reshaped_sol = [];
            sol_a = sol[locator:locator+r*(as-1)-1];
            locator = locator + r*(as-1);
            sol_b = sol[locator:locator+r*(bs-1)-1];
            locator = locator + r*(bs-1);
            sol_c = sol[locator:locator+r*cs-1];
            reshaped_sol = vec([sol_a; sol_b; sol_c])
        push!(sol_list, reshaped_sol)
    end
    return sol_list
end
ssols = shape_solutions(sols)
images = map(i -> evaluate(T, vec([a b c]) => i), ssols)

function group_similar_vectors(vectors; tol=1e-8)
    groups = Dict()
    for (i, v) in enumerate(vectors)
        found = false
        for (key, _) in groups
            if maximum(abs.(v - key)) < tol
                groups[key] = push!(groups[key], i)
                found = true
                break
            end
        end
	if !found
            groups[v] = [i]
        end
    end
    return groups
end

groups = group_similar_vectors(images)

@show groups
@show length(groups)


