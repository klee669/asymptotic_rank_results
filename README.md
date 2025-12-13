This repository contains the results and Julia scripts for the paper  
*“Asymptotic rank bounds: a numerical census.”*

The `.jl` files can be run line by line after setting the desired values of
`r`, `as`, `bs`, and `cs`.

`tensor_script.jl`: for codim 1 cases

`interpolation_tensor.jl`: for interpolation for codim > 1 cases.

**Caution.**  
If the expected codimension is not satisfied automatically, the value of
`cdim` must be specified manually.
