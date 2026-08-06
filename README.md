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
* **Ecosystem Ready**: Full out-of-the-box integration with `ChainRulesCore.jl` (Automatic Differentiation) and `ModelingToolkit.jl`.

## Advanced Infrastructure Capabilities

### 1. Automatic Differentiation (AD) Compatibility
The package defines custom `frule` and `rrule` expressions. It hooks directly into AD engines like `Zygote.jl`, `Enzyme.jl`, and `ForwardDiff.jl`, making them completely allocation-free and robust against array mutation errors:

```julia
using Zygote, PlanckFunctions

# Zygote automatically intercepts our custom analytical pullback for maximum performance
grad = Zygote.gradient(T -> ibb(1.5, T), 1500.0)
```

### 2. Physical Uncertainty Propagation
Because the core math is built using abstract type constraints, you can seamlessly track experimental error boundaries through your thermal equations using `Measurements.jl`:

```julia
using Measurements, PlanckFunctions

T = 1500.0 ± 2.5 # Temperature with error bounds
intensity = ibb(1.5, T) # Automatically returns an intensity value with precise error bounds
```

## Documentation
Full documentation, including advanced examples and API details, is available at the [Documentation Hub](https://manarom.github.io/PlanckFunctions.jl/).

## Contact & Contribution

To report bugs, suggest features, or contribute code, please submit an issue or pull request through the [GitHub Repository](https://github.com/Manarom/PlanckFunctions.jl).

