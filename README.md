# Rigorous numerics using Zernike Polynomials and existence proofs for PDEs on the disk



Table of contents:


* [Introduction](#introduction)
* [Rigorous numerics using Zernike Polynomials](#Rigorous-numerics-using-Zernike-Polynomials)
   * [Application to PDEs on the disk](#application-to-PDE-on-the-disk)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Julia code is a complement to the article 

#### [[1]](https://arxiv.org/abs/2411.18361) : "Validated matrix multiplication transform for orthogonal polynomials with applications to computer-assisted proofs for PDEs", M. Cadiot, J. Jaquette, J-P. Lessard and A. Takayasu, [ArXiv Link](https://arxiv.org/abs/2411.18361).

It provides the necessary rigorous computations of the bounds presented along the paper. The computations are performed using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). 

# Rigorous numerics using Zernike Polynomials



## Application to PDEs on the disk

The attached code allows to establish existence proofs of solutions to the following elliptic PDEs with Dirichlet boundary condition on the unit disk.

$$\Delta v (z) +  \overline{z}^m v^2(z) =0$$

$$\Delta v(z) +  {z}^{-1} v^2(z) =0$$

The code "proof_z_to_the_m.jl" treats the first equation for  $m = 0,1,2,20$ and the code "proof_1_over_r.jl" treats the second one. In each case we provide an approximate solution and prove the existence of a true solution in a vicinity of the approximate one.  The proof is computer-assisted and relies on the analysis derived in Section 5 of [[1]](https://arxiv.org/abs/2411.18361).

 
 # Utilisation and References
 
 The code is build using the following packages :

 - [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
 - [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
 - [JLD2](https://github.com/JuliaIO/JLD2.jl)
 
 
 # License and Citation
 
This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
  
If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{ZernikePolynomials.jl,
  author = {Matthieu Cadiot},
  title  = {ZernikePolynomials.jl},
  url    = {https://github.com/matthieucadiot/ZernikePolynomials.jl},
  note = {\url{ https://github.com/matthieucadiot/ZernikePolynomials.jl},
  year   = {2024},
}
```


# Contact

You can contact me at :

matthieu.cadiot@mail.mcgill.ca
