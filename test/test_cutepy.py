import cutepy as cute
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


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
