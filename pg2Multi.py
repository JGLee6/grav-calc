# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 15:01:34 2016

@author: cenpa
"""
import numpy as np
import sympy.physics.quantum.cg as cg
import scipy.special as sp

BIG_G = 6.67428e-11


def id2lm(idx):
    """
    Changes between specified index to (l,m) pair
    """
    ql = int(np.sqrt(2*idx))-1
    qm = idx-(ql*(ql+1))//2
    return ql, qm


def qmoment(l, m, massArray):
    """
    Computes the small q(l, m) inner multipole moment of a point mass array by
    evaluating the regular solid harmonic at each point-mass position.

    Inputs
    ------
    l : int
        Multipole moment order
    m : int
        Multipole moment order, m < l
    massArray : ndarray
        Nx4 array of point masses [m, x, y, z]

    Returns
    -------
    qlm : complex
        Complex-valued inner multipole moment
    """
    r = np.sqrt(massArray[:, 1]**2 + massArray[:, 2]**2 + massArray[:, 3]**2)
    theta = np.arccos(massArray[:, 3]/r)
    phi = np.arctan2(massArray[:, 2], massArray[:, 1]) % (2*np.pi)
    qlm = massArray[:, 0]*r**l*np.conj(sp.sph_harm(m, l, phi, theta))
    qlm = np.sum(qlm)
    return qlm


def Qmomentb(l, m, massArray):
    """
    Computes the large Q(l, m) outer multipole moment of a point mass array by
    evaluating the irregular solid harmonic at each point-mass position.

    Inputs
    ------
    l : int
        Multipole moment order
    m : int
        Multipole moment order, m <= l
    massArray : ndarray
        Nx4 array of point masses [m, x, y, z]

    Returns
    -------
    qlm : complex
        Complex-valued inner multipole moment
    """
    r = np.sqrt(massArray[:, 1]**2 + massArray[:, 2]**2 + massArray[:, 3]**2)
    if (r == 0).any():
        print('Outer multipole moments cannot be evaluated at the origin.')
        return 0
    theta = np.arccos(massArray[:, 3]/r)
    phi = np.arctan2(massArray[:, 2], massArray[:, 1]) % (2*np.pi)
    qlm = massArray[:, 0]*sp.sph_harm(m, l, phi, theta)/r**(l+1)
    qlm = np.sum(qlm)
    return qlm


def jmoment(l, m, massArray):
    r = np.sqrt(massArray[:, 1]**2 + massArray[:, 2]**2 + massArray[:, 3]**2)
    jvals = np.zeros(len(r))
    theta = np.arccos(massArray[:, 3]/r)
    phi = np.arctan2(massArray[:, 2], massArray[:, 1]) % (2*np.pi)
    for k in range(len(r)):
        jvals[k] = sp.sph_jn(l, r[k])[0][-1]
    ylm = massArray[:, 0]*jvals*np.conj(sp.sph_harm(m, l, phi, theta))
    ylm = np.sum(ylm)
    return ylm


def qmoments(l, massArray):
    """
    Computes all q(l, m) inner multipole moments of a point mass array up to a
    given maximum order, l. It does so by evaluating the regular solid harmonic
    at each point-mass position.

    Inputs
    ------
    l : int
        Maximum multipole moment order
    massArray : ndarray
        Nx4 array of point masses [m, x, y, z]

    Returns
    -------
    qlms : ndarry, complex
        Complex-valued inner multipole moments up to order l.
    """
    nlm = (l+1)*(l+2)//2
    qlms = np.zeros(nlm, dtype='complex')
    qlmsb = np.zeros([l+1, 2*l+1], dtype='complex')
    ctr = 0
    r = np.sqrt(massArray[:, 1]**2 + massArray[:, 2]**2 + massArray[:, 3]**2)
    theta = np.arccos(massArray[:, 3]/r)
    phi = np.arctan2(massArray[:, 2], massArray[:, 1]) % (2*np.pi)
    for n in range(l+1):
        rl = r**n
        for m in range(n+1):
            qlm = massArray[:, 0]*rl*np.conj(sp.sph_harm(m, n, phi, theta))
            qlms[ctr] = np.sum(qlm)
            qlmsb[n, l+m] = np.sum(qlm)
            print(n, m)
            ctr += 1
    qlmsb += np.conj(np.fliplr(qlmsb))
    qlmsb[:, l] /= 2
    return qlms, qlmsb


def Qmomentsb(l, massArray):
    """
    Computes all Q(l, m) outer multipole moments of a point mass array up to a
    a given maximum order, l. It does so by evaluating the irregular solid
    harmonic at each point-mass position.

    Inputs
    ------
    l : int
        Maximum multipole moment order
    massArray : ndarray
        Nx4 array of point masses [m, x, y, z]

    Returns
    -------
    qlms : ndarry, complex
        Complex-valued outer multipole moments up to order l.
    """
    nlm = (l+1)*(l+2)//2
    Qlms = np.zeros(nlm, dtype='complex')
    Qlmsb = np.zeros([l+1, 2*l+1], dtype='complex')
    ctr = 0
    r = np.sqrt(massArray[:, 1]**2 + massArray[:, 2]**2 + massArray[:, 3]**2)
    if (r == 0).any():
        print('Outer multipole moments cannot be evaluated at the origin.')
        return Qlms
    theta = np.arccos(massArray[:, 3]/r)
    phi = np.arctan2(massArray[:, 2], massArray[:, 1]) % (2*np.pi)
    for n in range(l+1):
        rl1 = r**(n+1)
        for m in range(n+1):
            Qlm = massArray[:, 0]*sp.sph_harm(m, n, phi, theta)/rl1
            Qlms[ctr] = np.sum(Qlm)
            Qlmsb[n, l+m] = np.sum(Qlm)
            print(n, m)
            ctr += 1
    Qlmsb += np.conj(np.fliplr(Qlmsb))
    Qlmsb[:, l] /= 2
    return Qlms, Qlmsb


def jmoments(l, massArray):
    r = np.sqrt(massArray[:, 1]**2 + massArray[:, 2]**2 + massArray[:, 3]**2)
    theta = np.arccos(massArray[:, 3]/r)
    phi = np.arctan2(massArray[:, 2], massArray[:, 1]) % (2*np.pi)

    nlm = (l+1)*(l+2)//2
    nP = len(r)
    ylms = np.zeros([nlm, nP], dtype='complex')
    jvals = np.zeros([nP, l+1])
    for k in range(nP):
        jvals[k] = sp.sph_jn(l, r[k])[0]

    ctr = 0
    for n in range(l+1):
        for m in range(n+1):
            print(l, m)
            ylms[ctr] = massArray[:, 0]*jvals[:, n]*np.conj(sp.sph_harm(m, n, phi, theta))
            ctr += 1

    ylms = np.sum(ylms, 1)
    return ylms


def torque_lm(qlm, Qlm, L=None):
    r"""
    Returns all gravitational torque_lm moments up to l=10 computed from sensor
    and source multipole moments. It assumes the sensor (interior) moments sit
    in a rotating frame of a turntable so that

    .. math::
        \bar{q_{lm}} = q_{lm}e^{-im\phi_{TT}}
    Then the torque is given by

    .. math::
        \tau = -4\pi i G \sum_{l=0}^{\infty}\frac{1}{2l+1}
        \sum_{m=-l}^{l}m\ q_{lm}Q_{lm}e^{-im\phi_{TT}}

    .. math::
        = 4\pi i G \sum_{l=0}^{\infty}\frac{1}{2l+1}\sum_{m=0}^{l}m\
        (q*_{lm}Q*_{lm}e^{im\phi_{TT}} - q_{lm}Q_{lm}e^{-im\phi_{TT}})

    Since the indices l and m are identical, we may simply do an element-wise
    multiplication and sum along rows.

    Inputs
    ------
    qlm : ndarray
        10x20 array of sensor (interior) lowest order multipole moments. The
        data should be lower triangular, with l denoting row and m/2 denoting
        real columns and m/2+1 denoting imaginary columns.
    Qlm : ndarray
        10x20 array of source (exterior) lowest order multipole moments. The
        data should be lower triangular, with l denoting row and m/2 denoting
        real columns and m/2+1 denoting imaginary columns.

    Returns
    -------
    nlm : ndarray
        10x20 array of torque multipole moments. The data should be lower
        triangular, with l denoting row and m/2 denoting cosine(m*phi) columns
        and m/2+1 denoting sine(m*phi) columns.
    """
    lqlm = len(qlm)
    lQlm = len(Qlm)
    minL = min([lqlm, lQlm])-1

    ls = np.arange(minL+1)
    lfac = 1/(2*ls+1)
    ms = np.arange(-minL, minL+1)
    nlm = 4*np.pi*BIG_G*1j*np.outer(lfac, ms)*qlm*Qlm

    nm = np.sum(nlm, 0)
    nc = nm[minL:] + nm[minL::-1]
    ns = nm[minL:] - nm[minL::-1]

    return nlm, nc, ns


def translate_qlm(qlm, rPrime):
    r"""
    Takes in a 10x20 q_lm interior set of moments up to l=10 and returns the
    10x20 q_LM interior moments up to l=10 translated by the vector rPrime.

    Inputs
    ------
    qlm : ndarray
        10x20 array of sensor (interior) lowest order multipole moments. The
        data should be lower triangular, with l denoting row and m/2 denoting
        real columns and m/2+1 denoting imaginary columns.
    rPrime : ndarray
        [x,y,z] translation from the origin where q_lm components were computed

    Returns
    -------
    qLM : ndarray
        10x20 array of sensor (interior) lowest order multipole moments. The
        data should be lower triangular, with l denoting row and m/2 denoting
        real columns and m/2+1 denoting imaginary columns. Moments are for the
        translated moments by rPrime.
    """
    rP = np.sqrt(rPrime[0]**2+rPrime[1]**2+rPrime[2]**2)
    thetaP = np.arctan2(rPrime[1], rPrime[0])
    phiP = np.arccos(rPrime[2]/rP)

    # Conjugate spherical harmonics
    ylmS = np.zeros([10, 20])
    qLM = np.zeros([10, 20])

    m = np.arange(10)
    for l in range(10):
        sphHarmL = sp.sph_harm(m[:l], l, thetaP, phiP)
        ylmS[l, :2*l:2] = sphHarmL.real
        ylmS[l, 1:2*l+1:2] = -sphHarmL.imag

    for L in range(10):
        for M in range(L):
            for l in range(L):
                lP = L-l
                mP = np.arange(lP)
                fac = np.sqrt(4.*np.pi*
                        np.exp(sp.gammaln(4*L+2)-sp.gammaln(2*lP+2)-sp.gammaln(2*l+2)))
                for m in range(l):
                    cgCoeff = cg.clebsch_gordan(lP, mP, l, m, L, M)
                    qLM[L, M] += fac*(rP**lP)*qlm[l, 2*m]*np.dot(cgCoeff, ylmS[lP, :2*lP:2])
                    qLM[L, M] += fac*(rP**lP)*qlm[l, 2*m+1]*np.dot(cgCoeff, ylmS[lP, 1:2*lP+1:2])

    return qLM
