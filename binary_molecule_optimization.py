"""
Discover the minimum energy configuration of a cluster of atoms.

Author:
    Ilias Bilionis

Date:
    7/18/2015

"""


import numpy as np
from ase import *
from ase.calculators.emt import EMT 
from ase.optimize import BFGS
from geometry import *
from plots import *
import GPy
import pydes 
import sys


def eval_en(x, mol):
    """
    Evaluate the energy of an atom.
    """
    mol.set_positions(x)
    return [mol.get_potential_energy()]


if __name__ == '__main__':
    # Fix the random seed in order to ensure reproducibility of results:
    np.random.seed(13456)
    # The molecule type (user)
    chem_formula = 'O2'
    # Number of initial data pool (user)
    n_init = 2
    # Number of candidate test points (user)
    n_design = 1000
    # Maximum iterations for BGO (user)
    max_it = 100
    # Tolerance for BGO (user)
    tol = 1e-3
    # Minimum and maximum distance between atoms (user)
    lower_dist = .6
    upper_dist = 2.
    # Start the algorithm
    print 'CLUSTER OPTIMIZATION USING BGO'.center(80)
    print '=' * 80
    print '{0:20s}: {1:8s}'.format('Chemical formula', chem_formula)
    print '{0:20s}: {1:d}'.format('Init. pool size', n_init)
    print '{0:20s}: {1:d}'.format('Design. pool size', n_design)
    print '{0:20s}: {1:d}'.format('Max BGO iter.', max_it)
    print '{0:20s}: {1:e}'.format('Tol. of BGO:', tol)
    print '=' * 80
    print '> starting computations'
    # A representation of the molecule in ase
    molecule = Atoms(chem_formula)
    # Assign an energy calculator to the molecule
    molecule.set_calculator(EMT())
    # Lower and upper bounds to distance matrix
    num_atoms = molecule.get_number_of_atoms()
    L = np.ones((num_atoms, num_atoms)) * lower_dist
    np.fill_diagonal(L, 0.)
    U = np.ones((num_atoms, num_atoms)) * upper_dist
    np.fill_diagonal(U, 0.)
    # The initial pool of coordinates in Cartesian form
    R_init, X_init = usample_many(L, U, n_init)
    # The energies of all the atoms in the initial pool
    E_init = np.array([eval_en(x, molecule) for x in X_init])
    # The design pool
    R_design, X_design = usample_many(L, U, n_design)
    # Let's build a Gaussian process that represents this.
    bgo = pydes.GlobalOptimizer(R_init, R_design, eval_en, args=(molecule,),
                                X_masked_init=X_init, X_masked_design=X_design)
    sys.stdout.write('> evaluating energies of all atoms in design pool... ')
    sys.stdout.flush()
    # To find the best achievable minimum:
    E_design = np.array([eval_en(x, molecule) for x in X_design])
    min_energy_over_design_pool = np.min(E_design)  # We want to see how fast BGO is
                                                    # going to find this
    sys.stdout.write('Done!\n')
    print '> Best achievable energy: {0:2.3f}'.format(min_energy_over_design_pool)
    print '> starting BGO'
    try:
        bgo.optimize(max_it=max_it, tol=tol,
                     callback_func=plot_1d_callback,
                     callback_func_args=(bgo,),
                     add_at_least=1)
    except KeyboardInterrupt:
        print '> keyboard interruption'
    X_best = bgo.best_masked_design
    molecule.set_positions(X_best)
    print '> plotting the results'
    dyn = BFGS(molecule)
    dyn.run()
    make_plots(bgo, molecule)