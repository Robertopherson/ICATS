"""Finite-ensemble checks against the configured preparation model."""
import json
from pathlib import Path

import numpy as np
from scipy.special import gammainc

from .constants import kboltz


def mean_statistics(values, reference=None):
    x = np.asarray(values, dtype=float)
    if not len(x):
        return {"n": 0}
    result = {"n": len(x), "mean": float(x.mean()),
              "standard_error": float(x.std(ddof=1)/np.sqrt(len(x))) if len(x) > 1 else None}
    if reference is not None:
        result["reference_mean"] = float(reference)
        result["mean_error"] = float(x.mean()-reference)
    return result


def distribution_statistics(counts, probabilities):
    counts = np.asarray(counts, dtype=float)
    p = np.asarray(probabilities, dtype=float)
    size = max(len(counts), len(p))
    counts = np.pad(counts, (0, size-len(counts)))
    p = np.pad(p, (0, size-len(p)))
    if np.any(counts < 0) or np.any(p < 0) or not np.isclose(p.sum(), 1):
        raise ValueError("Expected nonnegative counts and normalized probabilities")
    n = int(counts.sum())
    if n == 0:
        return {"n": 0}
    observed = counts/n
    return {"n": n, "total_variation": float(.5*np.abs(observed-p).sum()),
            "noise_scale": float(.5*np.sqrt(2/np.pi/n)*np.sqrt(p*(1-p)).sum()),
            "rare_reference_mass": float(p[n*p < 5].sum()),
            "outside_reference_count": int(counts[p == 0].sum()),
            "counts": counts.astype(int).tolist(), "reference": p.tolist()}


def discrete_statistics(values, probabilities):
    values = np.asarray(values, dtype=int)
    if np.any(values < 0):
        raise ValueError("State indices must be nonnegative")
    result = distribution_statistics(np.bincount(values), probabilities)
    result["moments"] = mean_statistics(values, np.dot(np.arange(len(probabilities)), probabilities))
    if len(values):
        labels = np.arange(len(probabilities))
        variance = np.dot(probabilities, (labels-result["moments"]["reference_mean"])**2)
        result["moments"]["reference_standard_error"] = float(np.sqrt(variance/len(values)))
    return result


def uniform_cdf_statistics(cdf_values, bins=20):
    """Use fixed equal-probability bins, keeping any out-of-support draws."""
    u = np.asarray(cdf_values, dtype=float)
    if not np.all(np.isfinite(u)):
        raise ValueError("Nonfinite sample CDF")
    counts = np.histogram(u, np.linspace(0, 1, bins+1))[0]
    counts = np.r_[counts, np.count_nonzero((u < 0) | (u > 1))]
    return distribution_statistics(counts, np.r_[np.full(bins, 1/bins), 0])


def quantum_j_probabilities(pa, pb, orbital_probabilities):
    """Unpolarized distinguishable-partner channel mixture, both parities."""
    pa, pb, pl = (np.asarray(p, dtype=float) for p in (pa, pb, orbital_probabilities))
    for p in (pa, pb, pl):
        if p.ndim != 1 or np.any(p < 0) or not np.isclose(p.sum(), 1):
            raise ValueError("Channel populations must be normalized")
    coef = np.zeros(len(pa)+len(pb)-1)
    for a, p_a in enumerate(pa):
        for b, p_b in enumerate(pb):
            coef[abs(a-b):a+b+1] += p_a*p_b/((2*a+1)*(2*b+1))
    # Difference-array accumulation avoids explicitly looping over every J.
    delta = np.zeros(len(coef)+len(pl))
    for s, c in enumerate(coef):
        for l, p_l in enumerate(pl):
            weight = c*p_l/(2*l+1)
            delta[abs(l-s)] += weight
            delta[l+s+1] -= weight
    j = np.arange(len(delta)-1)
    result = np.maximum(0, np.cumsum(delta)[:-1])*(2*j+1)
    if not np.isclose(result.sum(), 1, atol=1e-10):
        raise ValueError("Quantum channel count did not normalize")
    return result/result.sum()


def orbital_channel_probabilities(cap):
    """Map geometric magnitude cells [L,L+1) to quantum orbital channels."""
    if cap <= 0:
        raise ValueError("Positive orbital label cap required")
    r2 = cap*(cap+1)
    l = np.arange(int(np.ceil(np.sqrt(r2))))
    return np.maximum(0, np.minimum((l+1)**2, r2)-l**2)/r2


def rotor_energy_mean(mol, probabilities, temperature):
    if mol.sp.na == 1 or temperature <= 0:
        return 0.
    result = 0.
    for j, probability in enumerate(probabilities):
        if "asym" in mol.sp.rsym:
            levels = np.asarray(mol.sp.RR[j][1])
        elif mol.sp.rsym in ("linear", "spherical"):
            levels = np.zeros(1)
        else:
            levels = mol.sp.Jzc*np.arange(-j, j+1)**2
        conditional = np.exp(-(levels-levels.min())/(kboltz*temperature))
        conditional /= conditional.sum()
        result += probability*(mol.sp.J2c*j*(j+1)+np.dot(conditional, levels))
    return float(result)


def write_ensemble_statistics(sc, workers, directory):
    """Write a final report; never change the sampling or acceptance."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    report = {"notes": [
        "References describe the configured sampling model, not exact quantum dynamics.",
        "Standard errors assume independent draws; no universal convergence pass/fail.",
        "TV noise is a large-count estimate; rare_reference_mass flags sparsely sampled states.",
        "A zero observed SE does not establish precision for unobserved rare states; compare reference SE.",
        "Small marginal errors do not establish correct correlations or cross sections.",
        "Leading-Wigner phase-space entries pool modes at fixed excitation number."],
        "checks": {}, "skipped": {}}
    checks, skipped = report["checks"], report["skipped"]
    data, ip = sc.sdat, sc.ip
    if ip.usewang:
        skipped["unconditional references"] = "WL changes accepted populations; unweighted canonical references are not acceptance-conditioned."
    else:
        populations = []
        for mi, mol in enumerate(sc.mol):
            key = f"m{mi}"
            dist = workers[0].mol[mi].dist
            p = np.asarray(dist.get("rotJ", {}).get("probabilities", [1.0]))
            populations.append(p)
            checks[f"{key} rotational J"] = discrete_statistics(data["rot"][key]["J"], p)
            checks[f"{key} rotational eigenenergy (hartree)"] = {
                "moments": mean_statistics(data["rot"][key]["energy"],
                                           rotor_energy_mean(mol, p, ip.Trot))}
            vib = data["vib"][key]
            if ip.vib_mode == "rigid" or not mol.sp.nm:
                skipped[f"{key} vibration"] = "Rigid preparation or no vibrational modes."
                continue
            maxv = dist["nvib"]["MaxV"]
            n = np.arange(maxv+1)
            for mode, values in enumerate(vib["vi"]):
                w = mol.sp.w[mol.sp.ntr+mode]
                p_v = np.exp(-w*n/(kboltz*ip.Tvib)) if ip.Tvib > 0 else (n == 0).astype(float)
                p_v /= p_v.sum()
                result = discrete_statistics(values, p_v)
                result["omitted_thermal_mass"] = float(np.exp(-w*(maxv+1)/(kboltz*ip.Tvib))) if ip.Tvib > 0 else 0.
                checks[f"{key} mode {mode} excitation"] = result
            if mol.ip.nfreeze:
                skipped[f"{key} phase space"] = "Pooled records contain frozen modes."
                continue
            for state, (q, p_v) in enumerate(zip(vib["Q"], vib["P"])):
                q, p_v = np.asarray(q), np.asarray(p_v)
                if not len(q):
                    continue
                r2 = q*q+p_v*p_v
                result = uniform_cdf_statistics(gammainc(2*state+1, r2))
                result["moments"] = mean_statistics(r2/2, state+.5)
                result["Q_moments"] = mean_statistics(q, 0)
                result["P_moments"] = mean_statistics(p_v, 0)
                checks[f"{key} v={state} leading-Wigner energy"] = result
                checks[f"{key} v={state} phase"] = uniform_cdf_statistics(
                    (np.arctan2(p_v, q)+np.pi)/(2*np.pi))
        cap = float(ip.MaxL)
        if ip.FixedB is None and cap > 0:
            l = np.asarray(data["orb"]["iL"])
            geometric = getattr(ip, "orbital_sampling", "geometric") == "geometric" and ip.isotropic
            cdf = l*(l+1)/(cap*(cap+1)) if geometric else l/cap
            checks["orbital L"] = uniform_cdf_statistics(cdf)
            if geometric and all(m.sp.na == 1 or m.ip.ordist not in ("pdf", "fixed") for m in sc.mol):
                target = quantum_j_probabilities(*populations, orbital_channel_probabilities(cap))
                j = np.asarray(data["orb"]["iJ"])
                bins = np.floor(np.sqrt(j*(j+1))).astype(int)
                checks["total J quantum comparison"] = distribution_statistics(np.bincount(bins), target)
                checks["total J quantum comparison"]["reference_kind"] = "Quantum channel count; classical magnitude bins [J,J+1); finite-J approximation."
                result = checks["total J quantum comparison"]
                np.savetxt(directory / "ensemble_J.csv",
                           np.column_stack((np.arange(len(result["reference"])),
                                            result["reference"], np.array(result["counts"])/result["n"])),
                           delimiter=",", header="J,quantum,sampled", comments="")
            else:
                skipped["total J quantum comparison"] = "Requires geometric orbital sampling and unpolarized molecular preparation."
        else:
            skipped["orbital and total J"] = "Fixed-b or unsupported orbital preparation."
    checks["cos(L,Jab)"] = {"moments": mean_statistics(data["orb"]["cosLJab_thet"])}
    (directory / "ensemble_statistics.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    lines = ["ENSEMBLE STATISTICS", "", *report["notes"], "",
             "TV = total variation; noise = approximate finite-sample TV scale.",
             "Rare mass = reference probability in bins with expected count < 5.", ""]
    for name, result in checks.items():
        lines.append(name)
        if "total_variation" in result:
            lines.append(f"  N={result['n']}  TV={result['total_variation']:.5g}  noise~{result['noise_scale']:.5g}  rare mass={result['rare_reference_mass']:.5g}")
        for key in ("moments", "Q_moments", "P_moments"):
            m = result.get(key, {})
            if m.get("n", 0):
                se = "unavailable" if m["standard_error"] is None else f"{m['standard_error']:.4g}"
                line = f"  {key}: mean={m['mean']:.7g}, SE={se}"
                if "reference_mean" in m:
                    line += f", reference={m['reference_mean']:.7g}"
                if "reference_standard_error" in m:
                    line += f", reference SE={m['reference_standard_error']:.4g}"
                lines.append(line)
    lines += ["", "Not evaluated:"] + [f"  {k}: {v}" for k, v in skipped.items()]
    path = directory / "ensemble_statistics.txt"
    path.write_text("\n".join(lines)+"\n")
    return str(path)
