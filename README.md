# Residue Computation in Maxima

This project reworks the Maxima CAS code for computing residues—that is, the coefficient of the reciprocal term in a Laurent series. In future updates, we may revise parts of the definite integration code that depend on residue calculations.

Internally, the new user-level function `residue` uses an extensible dispatch system that selects among several well-known methods, including Taylor series expansion, power series techniques, and asymptotic Taylor methods. 

## What is Included
The repository includes:
- A new regression test file
- Updated source code for residue computation

## Usage

To use the package in Maxima, start by loading the file. To do this, copy the file "residue.lisp" to a folder
that Maxima can find and enter the command ``load("residue.lisp")$``.  After loading the file, to compute the residue of `expr` with respect to `var` at `point`, enter ``residue(expr, var, point)``.

## Examples

The package optionally prints information about which method was successful. 
Below are sample computations with informational messages enabled. 
~~~
(%i1) load("residue.lisp")$

(%i2) residue(1/(x^4+x+1)^2,x,a);
               4
Is 0 equal to a  + a + 1?

y;
Residue method: residue-rational succeeded.
                                      2
                                  36 a  + 48 a
(%o2)                   - ────────────────────────────
                              3        2
                          80 a  + 352 a  + 512 a + 255

(%i3) residue(exp(-1/x),x,0);
Residue method: residue-by-taylor-asym succeeded.
(%o3)                                 - 1

(%i4) residue(x*exp(-1/x),x,0);
Residue method: residue-by-taylor-asym succeeded.
                                       1
(%o4)                                  ─
                                       2

(%i5) residue(sin(x)/x^n,x,0);
   n - 2
Is ───── an integer?
     2

yes;
        n - 2
Is 0 <= ─────?
          2

yes;
Residue method: residue-by-powerseries succeeded.
                                       n - 2
                                       ─────
                                         2
                                  (- 1)
(%o5)                             ──────────
                                   (n - 1)!
(%i6) residue(tan(x)/x^n,x,0);
   n
Is ─ an integer?
   2

yes;
        n
Is 0 <= ─?
        2

yes;
Residue method: residue-by-powerseries succeeded.
                            n/2 - 1   n       n
                       (- 1)        (2  - 1) 2  bern(n)
(%o6)                  ────────────────────────────────
                                      n!
(%i7)

(%i8) residue(exp(%i*x*z)/(z^3+z+1),z,a);
    3
Is a  + a + 1 equal to 0?

y;
Residue method: residue-by-taylor succeeded.
                                     %i a x
                                   %e
(%o8)                              ────────
                                      2
                                   3 a  + 1
~~~

## Related Work
See also Guo Yicong’s [Maxima residue package](https://github.com/guoyicong/Maxima_residue/tree/master).

 
## License

This project is released under the [GNU General Public License v2.0](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).


## Contributing

We welcome contributions. You can add new residue methods, report bugs (either via the Maxima mailing list or the GitHub issue tracker), suggest improvements or better algorithms, or extend integration support.




