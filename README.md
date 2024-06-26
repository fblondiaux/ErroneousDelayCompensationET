# Postural Perturbation Simulation Code

This MATLAB code generates simulations for a postural perturbation task, investigating the impact of delay estimation errors on motor control.

## Overview

This code simulates a postural perturbation task with varying perturbation magnitudes and conditions. The simulations focus on the accurate or underestimation of the delay in the state estimator. The results from this code are used to generate Figure 4 and 5 in the associated paper.

## Simulations

The simulations are performed using the LQG (Linear Quadratic Gaussian) framework, integrating state and motor commands over a specified delay interval.

### Main Script: `runSim.m`

- Defines simulation parameters.
- Runs simulations for both "Healthy Controls" (HC) and "Essential Tremor" (ET) conditions.
- Generates arm angle, angular velocity, control signal, and Power Spectral Density (PSD) plots.

### Integration Function: `Libs/integrate.m`

- Integrates state and motor commands over the delay interval.
- Takes into account proportional errors in delay parameters.
- Can include noise in the delay estimation.

### Simulation Function: `Libs/simulation.m`

- Simulates a single movement using the LQG framework.
- Accounts for delay errors, perturbation magnitude, and scaling factors for matrix estimations.

### Remaining Files

- Remaining files are functions related to plotting or to small variations of the simulation parameters.

## Quickstart
Run the script 'runSim.m' to obtain the basic simulations and generate figures.

```
runSim; 
```

> 🚧 Success
> 
> Our implementation has been tested on MATLAB >= R2023a

Sample of output:
![runSim](https://github.com/fblondiaux/ErroneousDelayCompensationET/blob/main/Figures/runSim.png)

## Matlab toolboxes

    Statistics and Machine Learning Toolbox

## License

This code is provided under the GPL-3.0 license. See the LICENSE.md file for details.

## Reference
This code is a part of research article:  
> 🥈 Blondiaux, F., Colmant, L., Lebrun, L., Hanseeuw, B., & Crevecoeur, F. (2024). Erroneous compensation for long-latency feedback delays as origin of Essential Tremor. bioRxiv, 2024-01.

For any questions or issues, please refer to the corresponding author of the publication.

