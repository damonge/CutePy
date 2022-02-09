import cutepy as cute
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre


def test_correlate_scaled(exponent):
    nside = 256
    npix = hp.nside2npix(nside)
    th, _ = hp.pix2ang(nside, np.arange(npix))
    f, hm, bins = cute.correlate_scaled(np.array([0.0]),
                                        np.array([0.0]),
                                        np.array([1.0]),
                                        np.array([1.0]),
                                        th**exponent,
                                        np.ones(npix),
                                        0.1, 20)
    thm = 0.5*(bins[1:] + bins[:-1])
    plt.plot(thm, f, 'k.')
    plt.plot(thm, thm**exponent, 'r-')


test_correlate_scaled(1.)
test_correlate_scaled(1.5)
test_correlate_scaled(2.)
plt.show()


def test_2pcf(exponent):
    nside = 64
    npix = hp.nside2npix(nside)
    ls = np.arange(1000)
    cls = (nside/(ls+10.))**exponent

    # Estimate from 10 sims
    fs = []
    for i in range(10):
        mp = hp.synfast(cls, nside)
        f, hm, bins = cute.correlate_2pcf(mp, np.ones(npix), 20., 10,
                                          scale=0., smooth_kind='Gaussian')
        fs.append(f)
    fmean = np.mean(fs, axis=0)
    fstd = np.std(fs, axis=0)

    # Brute-force theory calculation
    thm = 0.5*(bins[1:] + bins[:-1])
    cth = np.cos(np.radians(thm))
    fpred = np.zeros_like(f)
    for l, cl in zip(ls, cls):
        leg = legendre(int(l))
        fpred += (2*l+1)*cl*leg(cth)/(4*np.pi)
    plt.errorbar(thm, fmean, yerr=fstd, fmt='k.')
    plt.plot(thm, fpred, 'r-')

test_2pcf(1.)
plt.show()
