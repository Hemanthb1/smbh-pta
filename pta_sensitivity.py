"""Continuous-wave (CW) sensitivity of the NANOGrav 12.5-yr data set with hasasia.

Approximations (state these when quoting numbers):
  * TOAs are epoch-averaged per receiver/backend (1-day groups), with ECORR
    added in quadrature -- the usual narrowband reduction.
  * The timing model is approximated by a quadratic spin-down, astrometry
    (position, proper motion, parallax) and one offset per backend; DMX is not
    modelled, so the curve is somewhat optimistic.
  * Intrinsic red noise from the public single-pulsar noise dictionary.
    Optionally the measured common-spectrum process / GWB is added as noise.
  * Optional "future" projection: each pulsar's last 3 yr of epochs are
    repeated forward in time with the same cadence and errors.
"""
import glob
import json
import os
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import hasasia.sensitivity as hsen
import hasasia.skymap as hsky

YR = 365.25 * 86400.0


def _read_par(path):
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            out[p[0]] = p[1]
    if "RAJ" in out:
        c = SkyCoord(out["RAJ"], out["DECJ"], unit=(u.hourangle, u.deg))
    else:
        lam = float(out.get("LAMBDA", out.get("ELONG")))
        bet = float(out.get("BETA", out.get("ELAT")))
        c = SkyCoord(lam * u.deg, bet * u.deg, frame="barycentrictrueecliptic")
    c = c.icrs
    return c.ra.rad, c.dec.rad


def _read_tim(path):
    """Return MJD, error [s], backend flag for uncommented TOA lines."""
    mjd, err, be = [], [], []
    for line in open(path):
        if line.startswith(("C ", "#", "FORMAT", "MODE", "INCLUDE")) or not line.strip():
            continue
        p = line.split()
        if len(p) < 5:
            continue
        try:
            m, e = float(p[2]), float(p[3])
        except ValueError:
            continue
        flag = p[p.index("-f") + 1] if "-f" in p else "default"
        mjd.append(m); err.append(e * 1e-6); be.append(flag)
    return np.array(mjd), np.array(err), np.array(be)


def _epoch_average(name, mjd, err, be, noise):
    """Average TOAs within 1 day per backend; return t[s], sigma[s], backend."""
    T, S, B = [], [], []
    for b in np.unique(be):
        m = be == b
        efac = noise.get(f"{name}_{b}_efac", 1.0)
        equad = 10 ** noise.get(f"{name}_{b}_log10_tnequad", -10.0)
        ecorr = 10 ** noise.get(f"{name}_{b}_log10_ecorr", -10.0)
        t = mjd[m]; s2 = (efac * err[m]) ** 2 + equad ** 2
        o = np.argsort(t); t, s2 = t[o], s2[o]
        start = 0
        for i in range(1, len(t) + 1):
            if i == len(t) or t[i] - t[start] > 1.0:
                w = 1.0 / s2[start:i]
                T.append(np.sum(w * t[start:i]) / w.sum())
                S.append(np.sqrt(1.0 / w.sum() + ecorr ** 2))
                B.append(b)
                start = i
    o = np.argsort(T)
    return np.array(T)[o] * 86400.0, np.array(S)[o], np.array(B)[o]


def _extend(t, s, b, extra_years):
    """Repeat the last 3 yr of epochs forward to mimic continued observing."""
    if extra_years <= 0:
        return t, s, b
    block = t >= t.max() - 3 * YR
    tt, ss, bb = [t], [s], [b]
    shift, k = 0.0, 0
    while shift < extra_years * YR:
        k += 1
        shift = k * 3 * YR
        keep = t[block] + shift <= t.max() + extra_years * YR
        tt.append(t[block][keep] + shift); ss.append(s[block][keep]); bb.append(b[block][keep])
    t2, s2, b2 = map(np.concatenate, (tt, ss, bb))
    o = np.argsort(t2)
    return t2[o], s2[o], b2[o]


def load_pulsars(data_dir, extra_years=0.0, A_gwb=None, gamma_gwb=13 / 3.0,
                 min_span_yr=3.0, verbose=True):
    noise = json.load(open(os.path.join(data_dir, "channelized_12p5yr_v3_full_noisedict.json")))
    tims = sorted(glob.glob(os.path.join(data_dir, "tim", "*.tim")))
    rn_freqs = None
    specs = []
    for tim in tims:
        name = os.path.basename(tim).split("_")[0]
        par = sorted(glob.glob(os.path.join(data_dir, "par", name + "_*.par")))[0]
        ra, dec = _read_par(par)
        mjd, err, be = _read_tim(tim)
        t, s, b = _epoch_average(name, mjd, err, be, noise)
        t, s, b = _extend(t, s, b, extra_years)
        span = (t.max() - t.min()) / YR
        if span < min_span_yr:
            continue
        # design matrix: quadratic + astrometry + per-backend offsets
        M = hsen.create_design_matrix(t, RADEC=True, PROPER=True, PX=True)
        backends = np.unique(b)
        if len(backends) > 1:
            jumps = np.array([(b == bk).astype(float) for bk in backends[1:]]).T
            M = np.hstack([M, jumps])
        # noise covariance: white + intrinsic red (+ optional GWB as noise)
        Tspan = t.max() - t.min()
        rn_freqs = np.logspace(np.log10(1 / (5 * Tspan)), np.log10(2e-7), 400)
        psd = np.zeros_like(rn_freqs)
        kA, kg = f"{name}_red_noise_log10_A", f"{name}_red_noise_gamma"
        if kA in noise:
            psd += hsen.red_noise_powerlaw(A=10 ** noise[kA], gamma=noise[kg], freqs=rn_freqs)
        if A_gwb:
            psd += hsen.red_noise_powerlaw(A=A_gwb, gamma=gamma_gwb, freqs=rn_freqs)
        N = np.diag(s ** 2) + hsen.corr_from_psd(freqs=rn_freqs, psd=psd, toas=t)
        psr = hsen.Pulsar(toas=t, toaerrs=s, phi=ra, theta=np.pi / 2 - dec,
                          designmatrix=M, N=N)
        psr.name = name
        specs.append(psr)
        if verbose:
            print(f"  {name}: {len(t):5d} epochs, {span:5.1f} yr")
    return specs


def spectra(psrs, freqs):
    out = []
    for p in psrs:
        sp = hsen.Spectrum(p, freqs=freqs)
        sp.NcalInv  # compute once
        out.append(sp)
    return out


def sky_threshold(specs, ra, dec, snr=1.0):
    """h0 needed for a given SNR at each source position, vs frequency.

    Returns (freqs, h0_thresh[n_freq, n_src], sky_avg[n_freq])."""
    ss = hsky.SkySensitivity(specs, np.pi / 2 - np.atleast_1d(dec), np.atleast_1d(ra))
    h_src = ss.h_thresh(SNR=snr)                       # (n_freq, n_src)
    dc = hsen.DeterSensitivityCurve(specs, pulsar_term=False)
    h_avg = snr * np.sqrt(dc.S_eff / dc.Tspan)
    return ss.freqs, h_src, h_avg
