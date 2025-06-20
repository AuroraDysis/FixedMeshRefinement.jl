# FixedMeshRefinement.jl: one-dimensional fixed mesh refinement

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AuroraDysis.github.io/FixedMeshRefinement.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AuroraDysis.github.io/FixedMeshRefinement.jl/dev/)
[![Build Status](https://github.com/AuroraDysis/FixedMeshRefinement.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AuroraDysis/FixedMeshRefinement.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AuroraDysis/FixedMeshRefinement.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AuroraDysis/FixedMeshRefinement.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/B/FixedMeshRefinement.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/B/FixedMeshRefinement.html)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# Examples

- [Wave](examples/wave/main.jl): Wave equation with reflection boundary condition.

```bash
cd examples/wave

# 3 levels, with default subcycling, patch is at the left of the domain
julia main.jl ./config/wave_3levels.toml

# 3 levels, with Mongwane's subcycling, patch is at the left of the domain
julia main.jl ./config/wave_3levels_mongwane.toml

# 3 levels, with default subcycling, patch is at the middle of the domain
julia main.jl ./config/wave_3levels_mid_patch.toml

# 3 levels, with Mongwane's subcycling, patch is at the middle of the domain
julia main.jl ./config/wave_3levels_mid_patch_mongwane.toml
```

# Acknowledgments

This package was originally adapted from [lwJi/Infino.jl](https://github.com/lwJi/Infino.jl), created by [lwJi](https://github.com/lwJi). Other useful resources are [JLRipley314/one-dim-amr: A simple one-dimensional AMR code](https://github.com/JLRipley314/one-dim-amr) and [AMReX-Codes/amrex: AMReX: Software Framework for Block Structured AMR](https://github.com/AMReX-Codes/amrex).

# References

- [[2503.09629] GPU-accelerated Subcycling Time Integration with the Einstein Toolkit](https://arxiv.org/abs/2503.09629)
