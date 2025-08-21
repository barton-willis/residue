Residue Computation in Maxima

This project reworks the Maxima CAS code for computing residues. In future updates, we may revise parts of the definite integration code that depend on residue calculations.

We've introduced a new user-level function `residue` which dispatches various methods, including several new ones. We've also updated some of the older internal functions to use this dispatch system.

So far, the repository includes:
- A new regression test file
- Updated source code for residue computation

Here are a few examples. For debugging the dispatch method optionally prints a message about which method is successful. 
~~~
(%i1) load("residu.lisp")$

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
~~~
