# Networked-Control-QuadTank


ALESSANDRO PEVERALI,
ANGELO MORONCELLI 
 
Overall, the linearized system dynamics (around a given equilibrium point) can be described by the following model:
𝑥̇ = 𝐴𝑥 + 𝐵1𝑢1 + 𝐵2𝑢2 𝑦𝑖 =𝐶𝑖𝑥,𝑖=1,2

Problem:

1. Generate the system matrices (both continuous-time and discrete-time, the latter with a sampling time of choice, selected according to the system dynamics). 
Perform the following analysis:
    a. Compute the eigenvalues and the spectral abscissa of the (continuous-time) system. Is it open-loop asymptotically stable?
    b. Compute the eigenvalues and the spectral radius of the (discrete-time) system. Is it open-loop asymptotically stable?

2. For different control structures (i.e., centralized, decentralized, and different distributed schemes) perform the following actions
    a. Compute the continuous-time fixed modes
    b. Compute the discrete-time fixed modes
    c. Compute, if possible, the CONTINUOUS-TIME control gains using LMIs to achieve the
       desired performances. Apply, for better comparison, different criteria for computing
       the control laws.
    d. Compute, if possible, the DISCRETE-TIME control gains using LMIs to achieve the
       desired performances. Apply, for better comparison, different criteria for computing
       the control laws.
    e. Analyze the properties of the so-obtained closed-loop systems (e.g., stability,
       eigenvalues) and compute the closed-loop system trajectories (generated both in continuous-time and in discrete-time) of the system variables, starting from a common random initial condition.


