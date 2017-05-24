# Fenton-Karma model of cardiac action potential

The Fenton-Karma model is a simplified, three-variable model of cardiac action potential orignally described in [this paper](http://aip.scitation.org/doi/10.1063/1.166311).

## How to cite

Ashikaga H. https://github.com/ashikagah/Fenton-Karma/

## Installation
Clone the github repository.
```
$ git clone https://github.com/ashikagah/Fenton-Karma
```

## Usage

In MATLAB command window, 

```
>> cd Fenton-Karma
>> demo
```
The demo script will:

1. Generate spiral wave breakup

The function `fk_breakup.m` is a MATLAB implementation of the Fenton-Karma model in two dimensions (2-D). The function will an s1-s2  stimulation sequence to induce one spiral wave that breaks up into multiple spiral waves in a 2-D lattice. It will save the time series of the excitation variable (_ts_) in a file `orig.mat` and a movie file `orig_movie.avi`.

2. Convert to phase map

The function `phase_map.m` will map the excitation variable (_ts_) to phase [-pi, pi] using Hilbert transform. It will save the time series of the phase (_p_) in a file `phase.mat` and a movie file `phase_movie.avi`.

3. Identify phase singularities

The function `phase_singularity.m` will detect the phase singularities (= rotors) using topological charge. It will save the time series of the phase singularities (_ps_) in a file `singularity.mat` and a movie file `singularity_movie.avi`.

## Variables
This model involves three normalized state variables:
- **_V_**: transmembrane voltage
- _**u**_: gate variable for inactivation of a fast inward current (_I~fi~_) after depolarization and its reactivation after repolarization
- **_w_**: gate variable for inactivation and reactivation of a slow inward current (_I~si~_)

## Spatial domain
- Matrix size: 500 x 500
- Grid spacing: 0.025 cm
- Grid size: 12.5 x 12.5 cm

## Numerical solution
Model equations are solved using a finite difference method for spatial derivatives and explicit Euler integration for time derivatives. Newman boundary conditions are assumed. 

## Licence
MIT

## References
1. Fenton FH and Karma A AD. Vortex dynamics in three-dimensional continuous myocardium with fiber rotation: Filament instability and fibrillationa. Chaos 8: 20-47, 1998
2. Bray MA and Wikswo JP. Use of topological charge to determine filament location and dynamics in a numerical model of scroll wave activity. _IEEE Trans Biomed Eng_ 49: 1086-93, 2002
3. Hammer P. Spiral waves in monodomain reaction-diffusion model. [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/22492-spiral-waves-in-monodomain-reaction-diffusion-model)
4. Spottiswoode BS. 2D phase unwrapping algorithms. [MATLAB File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/22504-2d-phase-unwrapping-algorithms?focused=5111677&tab=function)
