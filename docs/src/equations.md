# Equations

Here we will describe the system of equations that we will be solve

## N-Layer QG

The equations are

```math
\begin{aligned}
\partial_t q_n &= - R J(\psi_n , q_n) + R G_n'
\\
 q_1 &= R \nabla^2 \psi_1 + y - \frac{R L_y^2}{\delta_1 L_T^2} \left\{ \frac{\psi_1 - \psi_2}{\Delta \sigma_{21}}\right\}
  \\
 q_n &= R \nabla^2 \psi_n + y - \frac{R L_y^2}{\delta_n L_T^2} \left\{ \frac{\psi_n - \psi_{n-1}}{\Delta \sigma_{n,n-1}} +
\frac{\psi_{n+1} - \psi_{n}}{\Delta \sigma_{n+1,n}}
  \right\}
 \\
 q_N &= R \nabla^2 \psi_N + y - \frac{R L_y^2}{\delta_N L_T^2} \left\{ \frac{\psi_N - \psi_{N-1}}{\Delta \sigma_{N, N-1}}\right\}
 + \frac{f_0 d}{\beta L_y \delta_N}
\end{aligned}
```math
with
```math
\begin{aligned}
G_1' &= \frac{1}{\delta_1} \frac{U_s}{U_c} (w_0 - W_{21}) - \frac{\nu}{\beta L_y^5} \nabla^6 \psi_1
\\
G_n' &= \frac{n}{\delta_n} \frac{U_s}{U_c} (W_{n,n-1} - W_{n+1,n}) - \frac{\nu}{\beta L_y^5} \nabla^6 \psi_n
\\
G_N' &= \frac{N}{\delta_N} \frac{U_s}{U_c} (W_{N,N-1}) - \frac{\nu}{\beta L_y^5} \nabla^6 \psi_N - \frac{\epsilon}{\beta L_y} \nabla^2 \psi_N
\end{aligned}
```
