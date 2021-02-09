from math import sqrt
from numpy import linalg
import networkx as nx
from scipy.spatial import Delaunay, ConvexHull
import re
import itertools
import numpy as np
from scipy.spatial import ConvexHull




def unit_vector(x):
    return x/np.linalg.norm(x)


def unit_normal(a, b, c):
    x = np.linalg.det([[1, a[1], a[2]],
                       [1, b[1], b[2]],
                       [1, c[1], c[2]]])
    y = np.linalg.det([[a[0], 1, a[2]],
                       [b[0], 1, b[2]],
                       [c[0], 1, c[2]]])
    z = np.linalg.det([[a[0], a[1], 1],
                       [b[0], b[1], 1],
                       [c[0], c[1], 1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

# area of polygon poly


def poly_area(poly):
    if len(poly) < 3:  # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)




def fit_hypersphere(data, method="Hyper"):
    """returns a hypersphere of the same dimension as the 
        collection of input tuples
                (radius, (center))


    """
    num_points = len(data)

    if num_points == 0:
        return (0, None)
    if num_points == 1:
        return (0, data[0])
    dimen = len(data[0])        # dimensionality of hypersphere

    if num_points < dimen+1:
        raise ValueError(
            "Error: fit_hypersphere needs at least {} points to fit {}-dimensional sphere, but only given {}".format(dimen+1, dimen, num_points))

    # central dimen columns of matrix  (data - centroid)
    central = np.matrix(data, dtype=float)      # copy the data
    centroid = np.mean(central, axis=0)
    for row in central:
        row -= centroid

    square_mag = [sum(a*a for a in row.flat) for row in central]
    square_mag = np.matrix(square_mag).transpose()

    if method == "Taubin":

        mean_square = square_mag.mean()
        data_Z = np.bmat(
            [[(square_mag-mean_square)/(2*sqrt(mean_square)), central]])

        u, s, v = linalg.svd(data_Z, full_matrices=False)
        param_vect = v[-1, :]
        # convert from (dimen+1) x 1 matrix to list
        params = [x for x in np.asarray(param_vect)[0]]
        params[0] /= 2*sqrt(mean_square)
        params.append(-mean_square*params[0])
        params = np.array(params)

    else:

        data_Z = np.bmat([[square_mag, central, np.ones((num_points, 1))]])

        u, s, v = linalg.svd(data_Z, full_matrices=False)

        if s[-1]/s[0] < 1e-12:

            param_vect = v[-1, :]
            params = np.asarray(param_vect)[0]

        else:
            Y = v.H*np.diag(s)*v  # v.H gives adjoint
            Y_inv = v.H*np.diag([1./x for x in s])*v

            Ninv = np.asmatrix(np.identity(dimen+2, dtype=float))
            if method == "Hyper":
                Ninv[0, 0] = 0
                Ninv[0, -1] = 0.5
                Ninv[-1, 0] = 0.5
                Ninv[-1, -1] = -2*square_mag.mean()
            elif method == "Pratt":
                Ninv[0, 0] = 0
                Ninv[0, -1] = -0.5
                Ninv[-1, 0] = -0.5
                Ninv[-1, -1] = 0
            else:
                raise ValueError(
                    "Error: unknown method: {} should be 'Hyper', 'Pratt', or 'Taubin'")

            matrix_for_eigen = Y*Ninv*Y

            eigen_vals, eigen_vects = linalg.eigh(matrix_for_eigen)

            positives = [x for x in eigen_vals if x > 0]
            if len(positives)+1 != len(eigen_vals):
                # raise ValueError("Error: for method {} exactly one eigenvalue should be negative: {}".format(method,eigen_vals))
                print("Warning: for method {} exactly one eigenvalue should be negative: {}".format(
                    method, eigen_vals), file=stderr)
            smallest_positive = min(positives)

            A_colvect = eigen_vects[:, list(
                eigen_vals).index(smallest_positive)]

            param_vect = (Y_inv*A_colvect).transpose()

            # convert from (dimen+2) x 1 matrix to array of (dimen+2)
            params = np.asarray(param_vect)[0]



    radius = 0.5 * \
        sqrt(sum(a*a for a in params[1:-1]) - 4 *
             params[0]*params[-1])/abs(params[0])
    center = -0.5*params[1:-1]/params[0]

    center += np.asarray(centroid)[0]

    return (radius, center)