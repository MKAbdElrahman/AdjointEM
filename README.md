# AdjointEM
There  several goals for this project:
1. Give an introduction to basics of the most common computational methods 
in electromagnetics including the finite difference method, the finite element method and the method of moments.
2. Implement simple one dimensional example for each computational technique.
3. Understand constrained optimization  or inverse design of photonic devices.  
4. Explore the different approaches for inverse design and understand thier strengths and weaknesses. 
5. Understand the typical workflow when carrying inverse design, including
    - Choosing a set of design objectives.
    - Defining a set of  constraints (fabrication, material,footprint, bandwidth, radiation losses,directivity,   ..).
    - Understand the factors that help in discovering final robust  devices.
    - Choosing a suitable  filtering and projection technique for feature size control.

5- Once, we get the basics in 1D, we will move to 2D and 3D and add more material models (anisotropic, dispersive, nonlinear, .. ) 

At the end , we should be able to model and inverse design devices like waveguides bends, tapers, polarization rotators, wavelength multiplexers, spatial mode multiplexers,  nonlinear switches and  more!. We should be able to inverse design any photonic device :). 
## Under Construction 
### Handouts 
- This handout [handout](https://github.com/MKAbdElrahman/AdjointEM/blob/master/handouts/Approximate%20Finite%20Difference%20and%20%20Exact%20Automatic%20Deravtives/fd_auto_diff.pdf) introuduces the different appraoches used for calculating gradinets  or sensitivites of design objectives including the numerical finite difference and exact automatic differentiation.
- This [handout](https://github.com/MKAbdElrahman/AdjointEM/blob/master/handouts/AdjointEM%20Methods/Adjoint_EM.pdf) derives the adjoint method for  Maxwell's equation in time domain. 
