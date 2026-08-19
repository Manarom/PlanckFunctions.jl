# PlanckFunctions.jl

[![Build Status](https://github.com/Manarom/PlanckFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Manarom/PlanckFunctions.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/PlanckFunctions.jl)

`PlanckFunctions.jl` is a high-performance, type-stable Julia package designed for computing blackbody thermal emission spectra, integrated band power, and their exact analytical derivatives. 

Engineered specifically for optical physics, radiative heat transfer, and radiation pyrometry, the package bypasses numerical differentiation pitfalls by leveraging exact mathematical derivatives combined with robust numerical integration schemes.

---

## Key Features

* **Exact Analytical Derivatives**: Fast computation of first and second derivatives with respect to temperature ($T$) and wavelength ($\lambda$).
* **Combined Evaluation (`Dₜ` Functions)**: Evaluate a base function and its first two temperature derivatives simultaneously in a single pass, cutting execution time for root-finding algorithms (e.g., Halley's method).
* **Functional Differentiating Operators**: Use standard notation (`∇ₜ`, `∇²ₜ`, `∇ₗ`, `∇²ₗ`) to map a base radiation function directly to its respective analytical derivative function.
* **Domain-Specific Language (`∫ₗ` Operator)**: Express total, band, and discrete experimental data integrals using standard mathematical syntax.
* **Averaged Spectral Coefficients**: Built-in support for Rosseland- and Planck-averaged attenuation coefficients and total emittance.

---

## Physical Units

To ensure consistency across physical domains, the package strictly adheres to the following unit convention:
* **Wavelength (``\lambda``)**: Microns (`μm`)
* **Temperature (``T``)**: Kelvins (`K`)
* **Spectral Intensity (``I_{\lambda}``)**: Watts per square meter per steradian per micron (`W/m²⋅sr⋅μm`)
* **Total / Band Power (``P``)**: Watts per square meter per steradian (``W/m²⋅sr``)

---

## Core Usage & Examples

### 1. Basic Evaluation & Exact Derivatives
Evaluate the Planck blackbody distribution (`ibb`) and its respective temperature gradients:

```julia
using PlanckFunctions

λ = 1.5   # microns
T = 1500.0 # Kelvins

# Individual evaluations
intensity = ibb(λ, T)
dI_dT     = ∇ₜibb(λ, T)
d2I_dT2   = ∇²ₜibb(λ, T)

# Combined high-performance evaluation (returns a NTuple{3, Float64})
val, dval, d2val = Dₜibb(λ, T)
```

### 2. Functional Differentiating Operators (`∇ₜ`, `∇²ₜ`, `∇ₗ`, `∇²ₗ`)
Instead of manually typing specific derivative function names, you can use high-level differentiating operators directly as functional mappings.

```julia
using PlanckFunctions

# Temperature derivatives mapped dynamically across the core API
# (Supports: ibb, power, band_power, spectral_ratio, spectral_band_ratio, etc.)
(∇²ₜ(ibb) == ∇²ₜibb)     # true
(∇ₜ(∇ₜibb) == ∇²ₜibb)    # true

# Practical usage on a band power function
dPower_dT = ∇ₜ(band_power)
dPower_dT(1500.0; λₗ = 1.2, λᵣ = 4.0)
```

### 3. High-Level Integration via the `∫ₗ` Operator
The package provides a mathematical syntax sugar operator `∫ₗ` that constructs highly optimized, zero-allocation callable functors based on your integration domain.

```julia
using PlanckFunctions

# Unary Operator: Analytical full-spectrum integration (Stefan-Boltzmann law: σ * T^4)
total_power = ∫ₗ(ibb)
total_power(1500.0)

# Band Integration: Continuous integration over a specific spectral band [λ₁, λ₂]
band_integrator = ∫ₗ(ibb, 1.2, 4.0)
band_integrator(1500.0) # Returns integrated power within 1.2 - 4.0 μm

# Tabular Integration: Integrates discrete experimental data (e.g., spectral emissivity α) 
# against the Planck distribution over a wavelength grid λ
λ_grid = collect(range(0.5, 3.0, length=100))
α_data = rand(100) # Experimental spectral data

discrete_integrator = ∫ₗ(ibb, λ_grid, α_data)
discrete_integrator(1500.0)
```

---

## Infrastructure Capabilities

### 1. Automatic Differentiation (AD) Compatibility
`PlanckFunctions.jl` defines custom `frule` and `rrule` expressions via `ChainRulesCore.jl`. This allows you to seamlessly plug the package into any modern Julia AD framework (e.g., Zygote) without encountering numerical instabilities or performance degradation during gradient calculation.

```julia
using Zygote, PlanckFunctions

# Compute gradient of spectral intensity
Zygote.gradient(ibb, 2.0, 1500.0)

# AD of band power exactly matches the analytical temperature derivative function
Zygote.gradient(t -> band_power(t; λₗ = 1.2, λᵣ = 4.0), 1500.0)
```

### 2. Physical Uncertainty Propagation
With built-in support for `Measurements.jl`, you can feed values with experimental uncertainties directly into the functions. The package will propagate the errors mathematically using linear error propagation:

```julia
using Measurements, PlanckFunctions

# Automatically returns an intensity value with precise error bounds
I_measured = ibb(1.5, 1500.0 ± 2.5) 

# Propagate error through a multi-band spectral ratio calculation
ratio = spectral_band_ratio((1.0, 2.0), (2.0, 3.0), 1374.5 ± 0.5)
```

### 3. Symbolic Representation & Code Generation
When `Symbolics.jl` is loaded, `PlanckFunctions.jl` automatically activates a package extension enabling you to convert complex radiation equations into pure symbolic expressions for analytical manipulation or specialized code generation:

```julia
using Symbolics, PlanckFunctions
Symbolics.@variables λ T λ1 λ2
d2I_sym = symbolize(∇²ₜibb , λ, T) # Returns the symbolic representation of the Planck function second derivative
# Easily differentiate or simplify symbolic equations
d3I_sym = Symbolics.derivative(d2I_sym , T) # returns the third derivaitve 
d3I_numeric = build_function(d3I_sym, λ, T; expression=false) 
val = d3I_numeric(2.0, 1500.0) 
```

### 4. Continuous Integration & Leibniz Differentiation (QuadGK Extension)
When `QuadGK.jl` and `StaticArrays.jl` are loaded, an extension is automatically triggered to support high-accuracy continuous spectrum integration. 

```julia
using PlanckFunctions
using QuadGK, StaticArrays # Automatically activates the extension

# Define an analytical temperature-dependent emissivity property
ε(λ, T)   = 0.7 + 0.05 * λ + 1.2e-4 * T
dε(λ, T)  = 1.2e-4
ddε(λ, T) = 0.0

q = AnalyticalSpectralQuantity(ε, dε, ddε)

# Evaluate value, 1st and 2nd derivatives simultaneously on a band interval
# Ideal for zero-allocation root-finding solvers (e.g., Halley's method)
val, d1, d2 = Dₜplanck_weighted(q, 1.0, 1.6, 1450.0)

# Also supports multi-band ratios for ratio pyrometry
ratio_tuple = Dₜplanck_weighted_ratio(q, (1.0, 1.2), (1.4, 1.6), 1450.0)
```

---
## Documentation

Full documentation, including mathematical derivations, pyrometry-specific examples, performance benchmarks, and detailed API documentation, is available at the [Documentation Hub](https://manarom.github.io/PlanckFunctions.jl/).

## Contact & Contribution

`PlanckFunctions.jl` is an open-source project. If you encounter bugs, want to suggest new features (such as supporting multi-channel pyrometer calibrations), or wish to contribute code modifications:
1. Open an issue on the [GitHub Issues tracker](https://github.com).
2. Submit a Pull Request through the [GitHub Repository](https://github.com/Manarom/PlanckFunctions.jl).

We welcome contributions regarding more spectral averaging formulas, fast approximations, and integration benchmarks.
