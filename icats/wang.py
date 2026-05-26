"""Wang-Landau stored-umbrella helpers."""

from __future__ import annotations

import pickle
from typing import Any, Dict, Optional, Tuple

import numpy as np


def metadata_from_input(ip) -> Dict[str, Any]:
    return {
        "format": "icats-wang-umbrella",
        "version": 1,
        "maxj": int(ip.MaxJ),
        "maxl": int(ip.MaxL),
        "peak_jab": int(getattr(ip, "PeakJab", 0)),
        "wlmode": str(ip.wlmode),
        "wl_ff": float(ip.wl_ff),
        "wl_nstep_mult": int(ip.wl_nstep_mult),
        "wl_flatness": float(ip.wl_flatness),
        "wl_wn_factor": float(ip.wl_wn_factor),
        "wl_wn": None if ip.wl_wn is None else int(ip.wl_wn),
        "wl_tol": float(ip.wl_tol),
    }


def save(path: str, uu, iwld, td, metadata: Dict[str, Any]) -> None:
    payload = {
        "metadata": metadata,
        "uu": np.asarray(uu, dtype=float),
        "iwld": np.asarray(iwld, dtype=float),
        "td": np.asarray(td, dtype=float),
    }
    with open(path, "wb") as f:
        pickle.dump(payload, f)


def load(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[Dict[str, Any]]]:
    with open(path, "rb") as f:
        stored = pickle.load(f)

    metadata = None
    if isinstance(stored, dict):
        try:
            uu = stored["uu"]
            iwld = stored["iwld"]
            td = stored["td"]
        except KeyError as exc:
            raise ValueError(
                "Invalid Wang-Landau umbrella file: missing key "
                + str(exc)
                + " in "
                + path
            ) from exc
        metadata = stored.get("metadata")
    else:
        try:
            uu, iwld, td = stored
        except Exception as exc:
            raise ValueError("Invalid Wang-Landau umbrella file format: " + path) from exc

    return (
        np.asarray(uu, dtype=float),
        np.asarray(iwld, dtype=float),
        np.asarray(td, dtype=float),
        metadata,
    )


def validate(path: str, uu, iwld, td, metadata, expected_metadata: Dict[str, Any]) -> str:
    """Validate a stored umbrella.

    Returns a warning string for old metadata-free files, or an empty string.
    Raises ValueError when reuse would be unsafe.
    """
    maxj = int(expected_metadata["maxj"])
    if len(td) < maxj:
        raise ValueError(
            "Existing Wang-Landau umbrella is too short for this input: "
            + path
            + "\nRequested maxJ = "
            + str(maxj)
            + ", but stored td has only "
            + str(len(td))
            + " bins.\nMove or rename the existing wang.pkl, then rerun so a new "
            + "umbrella is generated. The existing file was not overwritten."
        )
    if len(td) == 0 or not np.all(np.isfinite(td)):
        raise ValueError("Invalid Wang-Landau umbrella weights in: " + path)

    if metadata is None:
        return (
            "WARNING: Wang-Landau umbrella has no metadata; "
            "only length/finite-value checks were possible.\n"
        )

    checks = [
        ("maxj", int),
        ("maxl", int),
        ("wlmode", str),
        ("wl_nstep_mult", int),
        ("wl_wn", lambda x: None if x is None else int(x)),
    ]
    mismatches = []
    for key, conv in checks:
        if conv(metadata.get(key)) != conv(expected_metadata.get(key)):
            mismatches.append(
                key
                + " stored="
                + str(metadata.get(key))
                + " requested="
                + str(expected_metadata.get(key))
            )

    for key in ("wl_ff", "wl_flatness", "wl_wn_factor", "wl_tol"):
        if key not in metadata or not np.isclose(float(metadata.get(key)), float(expected_metadata.get(key)), rtol=1e-10, atol=1e-12):
            mismatches.append(
                key
                + " stored="
                + str(metadata.get(key))
                + " requested="
                + str(expected_metadata.get(key))
            )

    if mismatches:
        raise ValueError(
            "Existing Wang-Landau umbrella metadata does not match this input: "
            + path
            + "\n"
            + "\n".join(mismatches)
            + "\nMove or rename the existing wang.pkl, then rerun so a new "
            + "umbrella is generated. The existing file was not overwritten."
        )

    return ""


def load_validated(path: str, expected_metadata: Dict[str, Any]):
    uu, iwld, td, metadata = load(path)
    warning = validate(path, uu, iwld, td, metadata, expected_metadata)
    return uu, iwld, td, warning
