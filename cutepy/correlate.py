from cutepy import cutelib as lib
import healpy as hp
import numpy as np


def correlate_scaled(theta, phi, weight, scale,
                     field, mask, x_max, nx):
    pos = (np.array([np.cos(theta), phi,
                     weight, scale]).T).flatten()
    nside = hp.npix2nside(len(field))
    res = lib.cute_correlation_wrap(pos, field, mask,
                                    x_max, nx, nside, 2*nx)
    hf, hm = res.reshape([2, nx])
    f = np.zeros_like(hf)
    f[hm > 0] = hf[hm > 0] / hm[hm > 0]
    bins = np.linspace(0., x_max, nx+1)
    return f, hm, bins


def correlate_scaled_2D(theta, phi, weight, scale,
                        field, mask, x_max, nx, na):
    pos = (np.array([np.cos(theta), phi,
                     weight, scale]).T).flatten()
    nside = hp.npix2nside(len(field))
    res = lib.cute_correlation_2D_wrap(pos, field, mask,
                                       x_max, nx, na, nside,
                                       2*nx*na)
    hf, hm = res.reshape([2, nx*na])
    f = np.zeros_like(hf)
    f[hm > 0] = hf[hm > 0] / hm[hm > 0]
    f = f.reshape([nx, na])
    hm = hm.reshape([nx, na])
    bins_x = np.linspace(0., x_max, nx+1)
    bins_a = np.linspace(0., 2*np.pi, na+1)
    return f, hm, bins_x, bins_a
