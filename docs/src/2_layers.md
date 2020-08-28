# Nondimensional 2-Layer QG without Coriolis
These equations are from [The vortex gas scaling regime of baroclinic turbulence](https://www.pnas.org/content/117/9/4491) by Basile Gallet and Raffaele Ferrari.

The equations are

```math
\begin{aligned}
\partial_t q_1 &= - J(\psi_1 , q_1) - \nu \Delta^4 q_1 + Q \sin( y / L)
\\
\partial_t q_2 &= - J(\psi_2 , q_2) - \nu \Delta^4 q_2 - 2 \kappa \nabla^2 \psi_2 - Q \sin( y / L)
\\
 q_1 &= \nabla^2 \psi_1 + \frac{1}{2 \lambda} (\psi_2 - \psi_1)
 \\
 q_2 &= \nabla^2 \psi_2 + \frac{1}{2 \lambda} (\psi_1 - \psi_2)
\end{aligned}
```
with 
```math
\begin{aligned}
J(a,b) &= (\partial_x a) (\partial_y b) - (\partial_x b) (\partial_y a ) 
\\ 
L = 1
\end{aligned}
```

