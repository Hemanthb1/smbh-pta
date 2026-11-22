"""Which periodic-quasar SMBHB candidates are within reach of a PTA?

For each candidate (Graham et al. 2015, CRTS):
  f_GW = 2 / P_obs            (optical period taken as the observed orbital period)
  M_c,z = (1+z) M q^{3/5} / (1+q)^{6/5},  M = virial BH mass taken as total mass
  h0   = hasasia's circular-binary amplitude with M_c,z, f_GW and D_L
Strain is shown for q = 0.5 (point) and q = 0.05-1 (error bar).
Sensitivity: NANOGrav 12.5-yr (real TOAs + noise models, see pta_sensitivity.py),
evaluated at each candidate's sky position, plus a +10-yr projection.

Usage: python run_pta_reach.py --data /path/to/12p5yr_stochastic_analysis/tutorials/data
"""
import argparse
import json
import os
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
import hasasia.skymap as hsky
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pta_sensitivity as ps

DAY = 86400.0
SNR_DET = 3.0          # nominal detection threshold
A_GWB = 1.92e-15       # NANOGrav 12.5-yr common-process amplitude (used as extra noise)
# Published NG12.5 CW search (Arzoumanian et al. 2023, ApJL 951, L28):
# sky-averaged 95% upper limit h0 < 6.82e-15 at f_GW = 7.65 nHz (CW + CRN model)
H95_PUB, F95_PUB = 6.82e-15, 7.65e-9


def chirp_mass(M, q):
    return M * q ** 0.6 / (1 + q) ** 1.2


def load_candidates(path):
    df = pd.read_csv(path)
    c = SkyCoord(df.ra, df.dec, unit=(u.hourangle, u.deg))
    df["ra_rad"], df["dec_rad"] = c.ra.rad, c.dec.rad
    df["f_gw"] = 2.0 / (df.period_obs_days * DAY)
    df["dL_mpc"] = Planck18.luminosity_distance(df.z.values).to(u.Mpc).value
    M = 10 ** df.logM.values
    for q in (0.05, 0.5, 1.0):
        Mcz = (1 + df.z.values) * chirp_mass(M, q)
        df[f"h0_q{q}"] = hsky.h0_circ(Mcz, df.dL_mpc.values, df.f_gw.values).decompose().value
    return df


def evaluate(df, freqs, h_src, label):
    """SNR of each candidate (q=0.5 and q=1) given per-source threshold at SNR=1."""
    snr05, snr1 = [], []
    for i, f in enumerate(df.f_gw.values):
        h1 = np.interp(np.log(f), np.log(freqs), h_src[:, i])
        snr05.append(df["h0_q0.5"].values[i] / h1)
        snr1.append(df["h0_q1.0"].values[i] / h1)
    df[f"snr_q0.5_{label}"] = snr05
    df[f"snr_q1_{label}"] = snr1
    df[f"h_needed_snr{SNR_DET:g}_{label}"] = SNR_DET * df["h0_q0.5"].values / np.array(snr05)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--cands", default="data/graham2015_candidates.csv")
    ap.add_argument("--future", type=float, default=10.0, help="extra years for projection")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    df = load_candidates(args.cands)
    print(f"{len(df)} candidates with masses")

    curves = {}
    setups = {
        "NG12.5": dict(extra_years=0.0, A_gwb=None),
        "NG12.5+GWB": dict(extra_years=0.0, A_gwb=A_GWB),
        f"+{args.future:g}yr+GWB": dict(extra_years=args.future, A_gwb=A_GWB),
    }
    for label, kw in setups.items():
        print(f"\nBuilding {label}")
        psrs = ps.load_pulsars(args.data, verbose=(label == "NG12.5"), **kw)
        Tspan = max(p.toas.max() for p in psrs) - min(p.toas.min() for p in psrs)
        freqs = np.logspace(np.log10(1 / Tspan), np.log10(1e-7), 300)
        specs = ps.spectra(psrs, freqs)
        f, h_src, h_avg = ps.sky_threshold(specs, df.ra_rad.values, df.dec_rad.values, snr=1.0)
        evaluate(df, f, h_src, label)
        curves[label] = (f, h_avg)
        print(f"  {len(psrs)} pulsars, Tspan = {Tspan / ps.YR:.1f} yr")

    # calibrate the approximate curve against the published 95% upper limit:
    # k = SNR at which our sky-averaged curve equals the published h95
    f, h = curves["NG12.5+GWB"]
    k = H95_PUB / np.interp(np.log(F95_PUB), np.log(f), h)
    df["h0_over_h95eq_q0.5"] = df["snr_q0.5_NG12.5+GWB"] / k
    df["h0_over_h95eq_q1"] = df["snr_q1_NG12.5+GWB"] / k
    print(f"\nCalibration: published h95 at {F95_PUB*1e9:.2f} nHz corresponds to SNR = {k:.2f} on our curve")
    df.to_csv(f"{args.out}/candidate_snr.csv", index=False)

    # ---------- figure ----------
    fig, ax = plt.subplots(figsize=(8.5, 6))
    styles = {"NG12.5": ("k", "-"), "NG12.5+GWB": ("C3", "-"),
              f"+{args.future:g}yr+GWB": ("C0", "--")}
    for label, (f, h) in curves.items():
        col, ls = styles[label]
        ax.loglog(f, SNR_DET * h, color=col, ls=ls, lw=1.6,
                  label=f"{label} sky-avg., SNR={SNR_DET:g}")
    ax.plot(F95_PUB, H95_PUB, "k*", ms=12, label="NG12.5 published 95% UL (sky-avg.)")
    yerr = [df["h0_q0.5"] - df["h0_q0.05"], df["h0_q1.0"] - df["h0_q0.5"]]
    ax.errorbar(df.f_gw, df["h0_q0.5"], yerr=yerr, fmt="o", ms=3.5, color="C2",
                ecolor="C2", alpha=0.6, lw=0.8, label="Graham+15 candidates (q=0.5; bar q=0.05–1)")
    ax.set_xlabel("GW frequency $f_{\\rm GW}=2/P_{\\rm obs}$ [Hz]")
    ax.set_ylabel("strain amplitude $h_0$")
    ax.set_xlim(2e-9, 1e-7)
    ax.set_ylim(1e-19, 1e-12)
    for key, lab in [("J164452", "SDSS J1644+4307"), ("PG 1302", "PG 1302$-$102")]:
        r = df[df.name.str.contains(key, regex=False)].iloc[0]
        ax.annotate(lab, (r.f_gw, r["h0_q0.5"]), xytext=(8, -14), textcoords="offset points",
                    fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.text(3.17e-8, 2e-13, " 1/yr (astrometry fit)", fontsize=7, color="0.4")
    ax.set_title("Periodic-quasar SMBHB candidates vs NANOGrav CW sensitivity")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(f"{args.out}/strain_vs_sensitivity.png", dpi=150)

    # ---------- summary ----------
    top = df.sort_values("snr_q1_NG12.5+GWB", ascending=False).head(10)
    cols = ["name", "z", "logM", "period_obs_days", "h0_q0.5", "snr_q0.5_NG12.5+GWB",
            "snr_q1_NG12.5+GWB", "snr_q1_+%gyr+GWB" % args.future, "h0_over_h95eq_q1"]
    summary = {k: int((df[f"snr_q1_{k}"] >= SNR_DET).sum()) for k in setups}
    summary_q05 = {k: int((df[f"snr_q0.5_{k}"] >= SNR_DET).sum()) for k in setups}
    print("\nTop 10 by SNR (q=1, NG12.5 with GWB noise):")
    print(top[cols].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"\nN with SNR>={SNR_DET:g} (q=1):", summary)
    print(f"N with SNR>={SNR_DET:g} (q=0.5):", summary_q05)
    json.dump({"n_candidates": len(df), "snr_threshold": SNR_DET,
               "n_detectable_q1": summary, "n_detectable_q05": summary_q05,
               "max_snr": {kk: float(df[f"snr_q1_{kk}"].max()) for kk in setups},
               "calibration_snr_at_published_h95": float(k)},
              open(f"{args.out}/summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
