Residue Computation in Maxima

This project reworks the Maxima CAS code for computing residues. In future updates, we plan to revisit parts of the definite integration code that depend on residue calculations.

We've introduced a new user-level function `residue` which dispatches various methods, including several new ones. We've also updated some of the older internal functions to use this dispatch system.

So far, the repository includes:
- A new regression test file
- Updated source code for residue computation


