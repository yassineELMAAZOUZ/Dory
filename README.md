
# Dory.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://a-kulkarn.github.io/Dory/dev)

Package to extend the functionality of Nemo/Hecke. Notable additions include:

## Basic utilities
- Allows Julia broadcasting for AbstractAlgebra matrices.
- Convenient constructors for AbstractAlgebra matrices.
- Indexing functions for AbstractAlgebra matrices.

## padic linear algebra:
- padic qr-factorization.
- padic singular value decomposition.
- padically stable solving of linear systems.
- padically stable hessenburg form.
- eigenvector solver (power and inverse iterations). [Only implemented for matrices defined over Qp]
- block schur form. [Only implemented for matrices defined over Qp]
