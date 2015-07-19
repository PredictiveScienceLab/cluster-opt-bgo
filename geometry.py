"""
Some functions that deal with the geometry of clusters.

Author:
    Ilias Bilionis

Date:
    7/18/2015

"""


__all__ = ['sample_coordinates', 'sample_many_coordinates',
           'get_distance_matrices', 'usample', 'floyd', 'mme',
           'usample_many']


import numpy as np
import scipy.spatial as spt 
import sys


def _sample_one_more(X, box, r):
    """
    Sample one more atom.
    """
    if X.shape[0] == 0:
        return box[:, 0] + (box[:, 1] - box[:, 0]) * np.random.rand(1, 3)
    while True:
        x = box[:, 0] + (box[:, 1] - box[:, 0]) * np.random.rand(1, 3)
        d = spt.distance_matrix(X, x)
        if (d > 2. * r).all():
            return x


def sample_coordinates(n, box=[(0, 1), (0, 1), (0, 1)], r=0.1):
    """
    Sample the coordinates of a cluster.

    :param n:   The number of atoms in the molecule.
    :param box: A box in which the molecule is confined.
    :param r:   The radius of each atom.
    """
    box = np.array(box)
    X = np.ndarray((n, 3))
    for i in xrange(n):
        X[i, :] = _sample_one_more(X[:i, :], box, r).flatten()
    return X


def sample_many_coordinates(s, n, box=None, r=.3):
    if box is None:
        #box = np.array([(0, 2 * n * 2 * r) * 3])
        box = np.array([(0, 2) * 3])
    X = np.ndarray((s, n, 3))
    for i in xrange(s):
        X[i, :, :] = sample_coordinates(n, box, r)
    return X


def get_distance_matrices(X):
    """
    Write
    """
    return np.array([spt.distance.pdist(x) for x in X])


def usample(L, U):
    """
    Samples a distance matrix from the configuration space C(L, U).
    """
    n = L.shape[0]
    D = np.zeros(L.shape)
    while True:
        L_in = L.copy()
        U_in = U.copy()
        for i in xrange(n - 1):
            for j in xrange(i + 1, n):
                D[i, j] = L_in[i, j] + np.random.rand() * (U_in[i, j] - L_in[i, j])
                D[j, i] = D[i, j]
                U_in[i, j] = D[i, j]
                U_in[j, i] = D[i, j]
                L_in[i, j] = D[i, j]
                L_in[j, i] = L_in[i, j]
                floyd(L_in, U_in)
        D, X = mme(D)
        if (L <= D).all() and (D <= U).all():
            break
    return D, X


def usample_many(L, U, num_samples):
    """
    Sample many distance matrices and the corresponding coordinates.
    """
    n = L.shape[0]
    X = np.ndarray((num_samples, n, 3))
    D = np.ndarray((num_samples, n * (n - 1) / 2))
    t = len(str(num_samples))
    for i in xrange(num_samples):
        sys.stdout.write('> sampling {0:s} of {1:s}\r'.format(str(i + 1).zfill(t),
                                                              str(num_samples)))
        sys.stdout.flush()
        d, x = usample(L, U)
        d = spt.distance.squareform(d)
        D[i, :] = d 
        X[i, :, :] = x
    sys.stdout.write('\n')
    return D, X


def floyd(L, U):
    """
    Update lower and upper bounds of the distance matrix by enforcing the triangle
    inequality.
    """
    n = L.shape[0]
    for k in xrange(n):
        for i in xrange(n - 1):
            for j in xrange(i + 1, n):
                if U[i, j] > U[i, k] + U[k, j]:
                    U[i, j] = U[i, k] + U[k, j]
                    U[j, i] = U[i, j]
                if L[i, j] < L[i, k] - U[k, j]:
                    L[i, j] = L[i, k] - U[k, j]
                    L[j, i] = L[i, j]
                if L[i, j] < L[j, k] - U[k, i]:
                    L[i, j] = L[j, k] - U[k, i]
                    L[j, i] = L[i, j]
                if L[i, j] > U[i, j]:
                    raise ValueError('Bad Bounds')


def mme(D):
    """
    Projects D to the closest distance matrix.
    """
    n = D.shape[0]
    W = np.ndarray((n, n))
    d_cm = np.ndarray((n,))
    for i in xrange(n):
        d_cm[i] = np.sum(D[i, :] ** 2) / n 
        for j in xrange(n):
            for k in xrange(j + 1, n):
                d_cm[i] -= D[j, k] ** 2 / n ** 2
    for i in xrange(n):
        for j in xrange(n):
            W[i, j] = .5 * (d_cm[i] + d_cm[j] - D[i, j] ** 2)
    lam, w = np.linalg.eig(W)
    lam = lam[::-1][:3]
    w = w[:, ::-1][:, :3]
    X = np.ndarray((n, 3))
    for i in xrange(min(3, n)):
        X[:, i] = np.sqrt(lam[i]) * w[:, i]
    return spt.distance.squareform(spt.distance.pdist(X)), X