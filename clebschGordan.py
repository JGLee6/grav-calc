# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 10:46:34 2020

@author: John
"""
import numpy as np
import scipy.special as sp

# XXX: Try to implement the Wei method using binomial coefficients (sp.comb)


def cgCoeff(j1, j2, m1, m2, J, M):
    """
    Calculate a specific Clebsch-Gordan coefficient based on the formula given
    on Wikipedia, < j1 j2; m1 m2 | J M >. This is useful for the calculation of
    translated multipole coefficients. At some point it may be useful to use a
    recursive calculation to speed this up.

    Inputs
    ------
    j1 : int
    j2 : int
    m1 : int
    m2 : int
    J : int
    M : int

    Returns
    -------
    CG : float
        Clebsch-Gordan coefficient < j1 j2; m1 m2 | J M >.

    Reference
    ---------
    *https://en.wikipedia.org/wiki/Table_of_Clebsch%E2%80%93Gordan_coefficients
    """
    # Make sure the M = m1 + m2
    if (j1 < abs(m1)) or (j2 < abs(m2)) or (J < abs(M)):
        print('Invalid (j1, m1), (j2, m2), or (J, M) pair')
        return 0
    # Make sure the combined angular momentum makes sense, J <= j1+j2
    if (J > (j1+j2)) or (J < abs(j1-j2)):
        print('Invalid combined angular momentum: J > j1+j2 or J < |j1-j2|')
        return 0
    delta = (M == (m1 + m2))
    # If we picked a negative M, then use symmetry property
    if M < 0:
        delta *= (-1)**(J-j1-j2)
        M = -M
        m1 = -m1
        m2 = -m2

    kmin = max([0, J-j1-j2, m1-j1, j1+m1-J])
    kmax = min([J-j2+m1, j2+m2])
    fac1 = sp.gammaln(J+j1-j2+1) + sp.gammaln(J-j1+j2+1)
    fac1 += sp.gammaln(j1+j2-J+1) - sp.gammaln(j1+j2+J+2)
    fac1 = np.sqrt((2*J+1)*np.exp(fac1))
    fac2 = sp.gammaln(J+M+1) + sp.gammaln(J-M+1)
    fac2 += sp.gammaln(j1-m1+1) + sp.gammaln(j1+m1+1)
    fac2 += sp.gammaln(j2-m2+1) + sp.gammaln(j2+m2+1)
    fac2 = np.sqrt(np.exp(fac2))
    sumCG = 0
    for k in range(kmin, kmax+1):
        denFac = sp.gammaln(k+1) + sp.gammaln(j1+j2-J-k+1)
        denFac += sp.gammaln(j1-m1-k+1) + sp.gammaln(j2+m2-k+1)
        denFac += sp.gammaln(J-j2+m1+k+1) + sp.gammaln(J-j1-m2+k+1)
        denFac = np.exp(denFac)
        sumCG += (-1)**k/denFac
    CG = delta*fac1*fac2*sumCG
    return CG
