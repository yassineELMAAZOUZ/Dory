# Dory.jl Documentation 


Dory is a package to extend the functionality of Hecke.jl.

The main features are:
- Convenience functions, such as various matrix get/set index methods for matrices supporting various types.
- Generic eigenvector methods.
- p-adic Linear algebra.

## Broadcasting

Broadcasting is enabled for the `AbstractAlgebra.Generic.MatElem{T}` type. the syntax is `your_function.(A)`.


The return type is either of type `AbstractAlgebra.Generic.MatElem{T}` or an array if the output type of `your_function` is not a subtype of `NCRingElem`.

## pAdic linear algebra.

```@docs
padic_qr
```


```@docs
rectangular_solve
```

```@docs
eigen
```

```@docs
eigvecs
```

```@docs
eigspaces
```

```@docs
singular_values
```