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


def _window_ell(ells, theta, smooth_kind):
    if smooth_kind == 'TopHat':
        win = np.ones_like(ells)
        win[ells > np.pi/theta] = 0.
    elif smooth_kind == 'Gaussian':
        win = np.exp(-ells*(ells+1)*theta**2)
    else:
        raise ValueError(f"Unknown smoothing type {smooth_kind}")
    return win


def correlate_line(field, mask, theta_max_deg, n_theta,
                   w_phase=True, scale='per_bin',
                   smooth_kind='TopHat', use_C=True):
    # Compute alms
    alm = hp.map2alm(field)

    # Compute phase alms if needed
    if w_phase:
        abs_alm = np.abs(alm)
        good = abs_alm > 0
        alm[good] = alm[good] / abs_alm[good]

    # Array of angles in radians
    theta_max = np.radians(theta_max_deg)

    # Back to map
    npix = len(field)
    nside = hp.npix2nside(npix)
    ells = np.arange(3*nside)
    if scale=='per_bin':
        per_bin = 1
        maps = np.zeros([n_theta, npix])
        thetas = np.linspace(0, theta_max, n_theta+1)
        # Mid-point
        thetas = 0.5*(thetas[1:]+thetas[:-1])
        for ith, th in enumerate(thetas):
            fl = _window_ell(ells, th, smooth_kind)
            maps[ith] = hp.alm2map(hp.almxfl(alm, fl), nside)
    else:
        per_bin = 0
        th = np.radians(scale)
        fl = _window_ell(ells, th, smooth_kind)
        maps = np.array([hp.alm2map(hp.almxfl(alm, fl), nside)])

    # Compute histograms
    if use_C:
        res = lib.cute_line_correlation_wrap(maps.flatten(), mask,
                                             theta_max, n_theta,
                                             nside, per_bin, 2*n_theta)
        hf, hm = res.reshape([2, n_theta])
    else:
        vecs = np.array(hp.pix2vec(nside, np.arange(npix))).T
        hf = np.zeros(n_theta)
        hm = np.zeros(n_theta)
        for i1 in range(npix):
            if mask[i1] <= 0:
                continue
            v1 = vecs[i1]
            listpix = hp.query_disc(nside, v1, theta_max*2.2)
            for i3 in listpix:
                if i3 <= i1:
                    continue
                if mask[i3] <= 0:
                    continue
                v3 = vecs[i3]
                v2 = 0.5*(v1+v3)
                i2 = hp.vec2pix(nside, v2[0], v2[1], v2[2])
                if mask[i2] <= 0:
                    continue
                cth = np.dot(v3,v1)
                theta = np.arccos(cth)*0.5
                i_theta = int(n_theta*theta/theta_max)
                if i_theta >= n_theta:
                    continue
                if per_bin:
                    mp = maps[i_theta]
                else:
                    mp = maps[0]
                hf[i_theta] += mp[i1]*mp[i2]*mp[i3]
                hm[i_theta] += 1

    # Combine into LCF
    f = np.zeros_like(hf)
    f[hm > 0] = hf[hm > 0] / hm[hm > 0]
    bins = np.linspace(0., theta_max_deg, n_theta+1)
    return f, hm, bins
