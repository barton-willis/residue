# Residue Computation in Maxima

This project reworks the Maxima CAS code for computing residues. In future updates, we may revise parts of the definite integration code that depend on residue calculations.

We've introduced a new user-level function `residue` which dispatches various methods, including several new ones. This system is designed to be readily extensible. We've also updated some of the older internal functions to use this new method dispatch system.

So far, the repository includes:
- A new regression test file
- Updated source code for residue computation

## Examples and Debugging

For debugging, the dispatch system optionally prints a message indicating which method succeeded. For these examples, we have turned on informational messages.

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
