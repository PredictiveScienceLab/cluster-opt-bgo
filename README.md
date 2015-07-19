# cluster-opt-bgo
Bayesian Global Optimization for Minimum Energy Cluster Identification
======================================================================

This package contains examples of application of Bayesian Global Optimization
(BGO) to the identification of minimum energy clusters. The purpose of the code
is educational and, I presume, it is quite easy to break it.

Potential Energy of Arbitrary Clusters
--------------------------------------

We use the [Atomistic Simulations Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)
 module for both the representation of clusters and the computation of their
energy.
Despite the fact that ASE can serve as an interface to a wide variety of ab
initio calculators, we opt for the use of the
[Effective Medium Theory (EMT)](https://wiki.fysik.dtu.dk/ase/ase/calculators/emt.html#module-ase.calculators.emt)
calculator,
since it is built in ASE and, therefore, it does not require external, 
potentially proprietary, software.
As a result, **do not expect to get the real minimum energies and/or structures**.
