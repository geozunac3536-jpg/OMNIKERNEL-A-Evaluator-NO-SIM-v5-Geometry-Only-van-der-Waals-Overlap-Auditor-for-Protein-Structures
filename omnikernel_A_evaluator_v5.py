# -*- coding: utf-8 -*-
"""
OMNIKERNEL-A — Evaluator puro (NO SIM) v5
- Σ: contactos CA-CA no vecinos (estructura/coherencia)
- φ: fricción material por overlap de radios de van der Waals (VdW)
     overlap = (r_i + r_j)*vdw_scale - d_ij
     contribuye si overlap > overlap_threshold
- Filtro químico: excluye pares con |Δresid| < steric_resid_sep_min
- HETATM opcional (include_hetatm)
- Auditoría: pares testeados, pares con overlap, suma y máximo overlap
"""

import json, hashlib, argparse
from pathlib import Path
import numpy as np


# Radios VdW típicos (Å). Suficientes para auditoría geométrica.
# (No pretende ser tabla completa de química cuántica; es un modelo material estable y trazable.)
VDW = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "SE": 1.90,
    # Metales/iones comunes (aprox. útiles para “contactos raros”)
    "MG": 1.73,
    "NA": 2.27,
    "K": 2.75,
    "CA": 2.31,
    "ZN": 1.39,
    "FE": 1.94,
    "CU": 1.40,
    "MN": 1.79,
    "CO": 2.00,
    "NI": 1.63
}
VDW_DEFAULT = 1.70  # fallback (C)


def sha256_of_dict(d: dict) -> str:
    b = json.dumps(d, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(b).hexdigest()[:16]


def infer_element(atom_name: str) -> str:
    """
    Infere símbolo químico desde el nombre PDB de átomo.
    Regla pragmática:
      - si empieza con letra, usa 1 o 2 letras
      - normaliza halógenos y metales comunes (CL, BR, ZN, FE, MG, NA, CA, etc.)
    """
    s = atom_name.strip().upper()
    if not s:
        return "C"

    # PDB: "CA" puede ser carbono alfa, pero elemento real es C.
    # Elemento se infiere mejor por primera letra si es atom de proteína:
    # "CA" (carbon alpha) -> C
    # Para iones HETATM, "CA" suele ser calcio (Ca). Diferenciación perfecta requiere columna ELEMENT,
    # pero aquí mantenemos modelo determinista:
    # - si atom_name es exactamente "CA" y no es backbone context, podría ser Ca.
    #   Sin columna ELEMENT, tomamos criterio conservador: tratar "CA" como C por defecto.
    #   (Si quieres priorizar ion Ca, cambia abajo.)
    if s == "CA":
        return "C"

    # Halógenos/metales de 2 letras si corresponde
    two = s[:2]
    if two in {"CL", "BR", "ZN", "FE", "MG", "NA", "MN", "CO", "NI", "CU", "SE"}:
        return two

    # Iodo puede venir como " I" o "I1"
    if s.startswith("I"):
        return "I"

    # Si empieza con letra, toma 1 letra
    ch0 = s[0]
    if ch0.isalpha():
        return ch0

    return "C"


def read_pdb_atoms(path: Path, chain=None, include_hetatm=True):
    """
    Lee átomos pesados desde PDB (ATOM y opcional HETATM), excluye H.
    Devuelve:
      coords (N,3), resid (N,), elem (N,), atom_name (N,)
    """
    coords = []
    resid = []
    elem = []
    aname = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[0:6].strip()
            if rec == "ATOM":
                pass
            elif rec == "HETATM" and include_hetatm:
                pass
            else:
                continue

            # cadena
            ch = line[21].strip()
            if chain is not None and ch != chain:
                continue

            atom = line[12:16].strip()
            # Excluir H por nombre
            if atom.upper().startswith("H"):
                continue

            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                r = int(line[22:26])
            except Exception:
                continue

            e = infer_element(atom)
            # Excluir H por elemento inferido
            if e == "H":
                continue

            coords.append((x, y, z))
            resid.append(r)
            elem.append(e)
            aname.append(atom)

    if not coords:
        return (np.zeros((0,3), float),
                np.zeros((0,), int),
                np.array([], dtype=object),
                np.array([], dtype=object))

    return (np.array(coords, float),
            np.array(resid, int),
            np.array(elem, dtype=object),
            np.array(aname, dtype=object))


def read_pdb_CA(path: Path, chain=None):
    coords = []
    resid = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            ch = line[21].strip()
            if chain is not None and ch != chain:
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                r = int(line[22:26])
            except Exception:
                continue
            coords.append((x, y, z))
            resid.append(r)
    if not coords:
        return np.zeros((0,3), float), np.zeros((0,), int)
    return np.array(coords, float), np.array(resid, int)


def compute_sigma_CA(ca_coords: np.ndarray, sep_min: int, cutoff: float, weight: float):
    N = ca_coords.shape[0]
    if N < 4:
        return 0.0
    D = np.linalg.norm(ca_coords[:, None, :] - ca_coords[None, :, :], axis=2) + 1e-12
    sigma = 0.0
    for i in range(N):
        for j in range(i + sep_min, N):
            d = float(D[i, j])
            if d < cutoff:
                sigma += weight / (d * d + 0.2)
    return float(sigma)


def compute_phi_vdw(coords: np.ndarray,
                    resid: np.ndarray,
                    elem: np.ndarray,
                    resid_sep_min: int,
                    vdw_scale: float,
                    overlap_thr: float,
                    phi_weight: float):
    """
    φ por overlap VdW:
      overlap = (r_i + r_j)*vdw_scale - d
      contribuye si overlap > overlap_thr
      contribución: phi_weight * overlap/(d+0.05)  (suave, no explosiva)
    Devuelve: phi, overlap_pairs, tested_pairs, max_overlap
    """
    M = coords.shape[0]
    if M < 2:
        return 0.0, 0, 0, 0.0

    # radios por átomo
    r = np.array([VDW.get(str(e), VDW_DEFAULT) for e in elem], float)

    D = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2) + 1e-12

    phi = 0.0
    overlap_pairs = 0
    tested = 0
    max_ov = 0.0

    for i in range(M):
        ri = int(resid[i])
        for j in range(i + 1, M):
            rj = int(resid[j])
            if abs(ri - rj) < resid_sep_min:
                continue
            tested += 1
            d = float(D[i, j])
            ov = (r[i] + r[j]) * vdw_scale - d
            if ov > overlap_thr:
                overlap_pairs += 1
                if ov > max_ov:
                    max_ov = ov
                phi += phi_weight * (ov / (d + 0.05))

    return float(phi), int(overlap_pairs), int(tested), float(max_ov)


def entropy_proxy_CA(ca_coords: np.ndarray) -> float:
    if ca_coords.shape[0] == 0:
        return 0.0
    cm = ca_coords.mean(axis=0, keepdims=True)
    rg = float(np.sqrt(np.mean(np.sum((ca_coords - cm) ** 2, axis=1))))
    return float(np.log(rg + 1.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--brain", required=True)
    ap.add_argument("--pdbA", required=True)
    ap.add_argument("--pdbB", default=None)
    ap.add_argument("--out", default="out")
    args = ap.parse_args()

    brain = json.loads(Path(args.brain).read_text(encoding="utf-8"))
    cfg_hash = sha256_of_dict(brain)
    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)

    chain = brain["inputs"].get("chain", None)
    include_hetatm = bool(brain["inputs"].get("include_hetatm", True))

    m = brain["metrics"]
    sep_min = int(m["sequence_separation_min"])
    contact_cutoff = float(m["contact_cutoff_Ang"])
    sigma_w = float(m["sigma_weight_base"])

    resid_sep_min = int(m.get("steric_resid_sep_min", 2))
    vdw_scale = float(m.get("vdw_scale", 1.0))
    overlap_thr = float(m.get("overlap_threshold_Ang", 0.0))
    phi_w = float(m.get("phi_weight_base", 1.0))

    report_norm = bool(m.get("report_normalized", True))

    # --- A ---
    caA, _ = read_pdb_CA(Path(args.pdbA), chain=chain)
    heavyA, residA, elemA, anameA = read_pdb_atoms(Path(args.pdbA), chain=chain, include_hetatm=include_hetatm)

    sigmaA = compute_sigma_CA(caA, sep_min, contact_cutoff, sigma_w)
    phiA, ovPairsA, testedA, maxOvA = compute_phi_vdw(
        heavyA, residA, elemA, resid_sep_min, vdw_scale, overlap_thr, phi_w
    )
    HA = entropy_proxy_CA(caA)

    report = {
        "schema": "tcds.omnikernel.evaluator.report.v5",
        "config_hash": cfg_hash,
        "mode": "evaluator_only",
        "inputs": {
            "pdbA": Path(args.pdbA).name,
            "pdbB": Path(args.pdbB).name if args.pdbB else None,
            "chain": chain,
            "include_hetatm": include_hetatm
        },
        "params": {
            "contact_cutoff_Ang": contact_cutoff,
            "sequence_separation_min": sep_min,
            "steric_resid_sep_min": resid_sep_min,
            "vdw_scale": vdw_scale,
            "overlap_threshold_Ang": overlap_thr,
            "phi_weight_base": phi_w
        },
        "A": {
            "N_CA": int(caA.shape[0]),
            "N_atoms_heavy": int(heavyA.shape[0]),
            "Sigma": float(sigmaA),
            "Phi": float(phiA),
            "Phi_overlap_pairs": int(ovPairsA),
            "Phi_pairs_tested": int(testedA),
            "Phi_max_overlap_Ang": float(maxOvA),
            "H_proxy": float(HA)
        },
        "B": None,
        "DeltaS": None,
        "normalized": None,
        "validation": {"balance_ok": None, "e_veto_ok": None, "verdict": None, "notes": []}
    }

    if report_norm:
        nca = max(report["A"]["N_CA"], 1)
        nh = max(report["A"]["N_atoms_heavy"], 1)
        report["normalized"] = {
            "A": {
                "Sigma_per_CA": report["A"]["Sigma"] / nca,
                "Phi_per_heavy_atom": report["A"]["Phi"] / nh,
                "Overlap_pairs_per_heavy_atom": report["A"]["Phi_overlap_pairs"] / nh
            }
        }

    # --- B (opcional, para E-Veto conmensurable) ---
    if args.pdbB:
        caB, _ = read_pdb_CA(Path(args.pdbB), chain=chain)
        heavyB, residB, elemB, anameB = read_pdb_atoms(Path(args.pdbB), chain=chain, include_hetatm=include_hetatm)

        sigmaB = compute_sigma_CA(caB, sep_min, contact_cutoff, sigma_w)
        phiB, ovPairsB, testedB, maxOvB = compute_phi_vdw(
            heavyB, residB, elemB, resid_sep_min, vdw_scale, overlap_thr, phi_w
        )
        HB = entropy_proxy_CA(caB)
        dS = float(HB - HA)

        report["B"] = {
            "N_CA": int(caB.shape[0]),
            "N_atoms_heavy": int(heavyB.shape[0]),
            "Sigma": float(sigmaB),
            "Phi": float(phiB),
            "Phi_overlap_pairs": int(ovPairsB),
            "Phi_pairs_tested": int(testedB),
            "Phi_max_overlap_Ang": float(maxOvB),
            "H_proxy": float(HB)
        }
        report["DeltaS"] = dS

        if report_norm:
            ncaB = max(report["B"]["N_CA"], 1)
            nhB = max(report["B"]["N_atoms_heavy"], 1)
            report["normalized"]["B"] = {
                "Sigma_per_CA": report["B"]["Sigma"] / ncaB,
                "Phi_per_heavy_atom": report["B"]["Phi"] / nhB,
                "Overlap_pairs_per_heavy_atom": report["B"]["Phi_overlap_pairs"] / nhB
            }

    # --- Validación ---
    v = brain["validation"]
    Q = float(v.get("Q_eff_constant", 1.0))
    use_balance = bool(v.get("use_balance_law", True))
    balance_ok = True if not use_balance else (Q * report["A"]["Sigma"] > report["A"]["Phi"])
    report["validation"]["balance_ok"] = bool(balance_ok)

    # E-Veto solo si hay ΔS y es conmensurable por N_CA
    e = v.get("e_veto", {"enabled": False})
    if bool(e.get("enabled", False)) and report["DeltaS"] is not None:
        if bool(e.get("requires_conmensurable", True)) and (report["B"]["N_CA"] != report["A"]["N_CA"]):
            report["validation"]["e_veto_ok"] = None
            report["validation"]["notes"].append("E_VETO_SKIPPED: non-conmensurable (N_CA differs)")
        else:
            thr = float(e.get("deltaS_threshold", -0.2))
            report["validation"]["e_veto_ok"] = bool(report["DeltaS"] <= thr)
    else:
        report["validation"]["e_veto_ok"] = None

    verdict = "PASS" if report["validation"]["balance_ok"] else "FAIL_BALANCE"
    if report["DeltaS"] is not None and report["validation"]["e_veto_ok"] is False:
        verdict = "FAIL_E_VETO"
    report["validation"]["verdict"] = verdict

    out_path = out_dir / f"{brain['logging'].get('out_prefix','omnikernel_A_eval_v5')}_{cfg_hash}.json"
    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("✅ Evaluación completa (NO SIM, v5). Reporte:", out_path)


if __name__ == "__main__":
    main()
