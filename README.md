# PlanckFunctions

[![Build Status](https://github.com/Manarom/PlanckFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Manarom/PlanckFunctions.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/PlanckFunctions.jl)

## General Description
PlanckFunctions.jl is a high-performance Julia package designed to calculate blackbody thermal emission spectra and their analytical derivatives.

### Key Features
* **Spectral Intensity**: Evaluate the precise blackbody spectral intensity (Planck function).
* **Analytical Derivatives**: Compute exact first and second derivatives with respect to wavelength and temperature.
* **Fast Band Integration**: Perform rapid, customized integration over specified wavelength regions.
* **Averaged Coefficients**: Calculate Rosseland- and Planck-averaged attenuation coefficients.
* **Ecosystem Ready**: Integration with `ChainRulesCore.jl` (Automatic Differentiation) 

## Advanced Infrastructure Capabilities

### 1. Automatic Differentiation (AD) Compatibility
The package defines custom `frule` and `rrule` expressions for main functions 

```julia
using Zygote, PlanckFunctions

Zygote.gradient(ibb, 2.0 , 1500.0)
Zygote.gradient(t->band_power(t , ; λₗ = 1.2 , λᵣ =4.0), 1500.0) 
∇ₜband_power(1500.0 , ; λₗ = 1.2 , λᵣ =4.0) # the same analytically

```
### 2. Physical Uncertainty Propagation
Uncertainty propagation using  `Measurements.jl`:

```julia
using Measurements, PlanckFunctions

ibb(1.5, 1500.0 ± 2.5) # Automatically returns an intensity value with precise error bounds
spectral_band_ratio((1.0,2.0) , (2.0,3.0) , 1374.5 ± 0.5 )
```
### 3. Symbolic representation

If you load `Symbolics.jl`, `PlanckFunctions` automatically activates an extension allowing you to convert functions into symbolic expressions:

```julia
using PlanckFunctions
using Symbolics

ibb_sym = symbolize(ibb) # returns symbolic representation of Planck function 
```
## Documentation
Full documentation, including advanced examples and API details, is available at the [Documentation Hub](https://manarom.github.io/PlanckFunctions.jl/).

## Contact & Contribution

To report bugs, suggest features, or contribute code, please submit an issue or pull request through the [GitHub Repository](https://github.com/Manarom/PlanckFunctions.jl).

