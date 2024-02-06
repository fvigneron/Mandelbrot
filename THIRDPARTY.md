# Attribution notice

## Included libraries

- This software includes the [MPFR](https://www.mpfr.org) library to handle arbitrary high-precision arithmetic,
which is distributed under [GNU LGPL-3](https://www.gnu.org/licenses/lgpl-3.0.html).
- Though we don't include it explicitely, MPFR includes in turn the [GMP](https://gmplib.org) library,
which can either be distributed under [GNU LGPL-3](https://www.gnu.org/licenses/lgpl-3.0.html)
or [GNU GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html).

## Alternative libraries

Note that the applications `hypQuick` and `misQuick` use hardware arithmetic (FP80 numbers).

We **welcome** the idea of using alternative libraries for high-precision arithmetic, in particular if your
favorite library works on GPU or in a mixed CPU/GPU environment. If you are interested, please contact me
(<francois.vigneron@univ-reims.fr>) so we can discuss evolving `Mandelbrot` together :wink:.
