title: test_updclie


# Description

Test of climatology routine `updclie.F90`. Same as basic forecast test except
with:
```fortran
  LMCCEC=true,
  LMCCIEC=false,
```
to cause the climate fields to vary in time. Climate fields are updated without
time interpolation.
