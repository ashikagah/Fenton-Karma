# Fenton-Karma model of cardiac action potential

The Fenton-Karma model is a simplified, three-variable model of cardiac action potential orignally described in [this paper](http://aip.scitation.org/doi/10.1063/1.166311).

## Installation
Clone the github repository.
```
$ git clone https://github.com/ashikagah/Fenton-Karma
```

## Usage
The script `FK_2D.m` is a MATLAB implementation of the Fenton-Karma model in two dimensions. Simply select the whole script by `Ctrl+A`, and right-click to select `Evaluate Selection`. The script will create a movie file `FK2D_00_05.avi`.

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
Model equations are solved using a finite difference method for spatial derivatives and explicit Euler integration for time derivatives. Newman boundary conditions are assumed. Model parameters are taken from the parameter set 6 in `2002_FentonFH_Chaos.pdf`. 

## Licence
MIT
