# OMNIKERNEL-A Evaluator (NO-SIM) v5 — VdW Overlap Edition

**DOI:** https://doi.org/10.5281/zenodo.18194092  
**Tipo:** Software (Zenodo)  
**Licencia:** Apache License 2.0

## Qué es
OMNIKERNEL-A (NO-SIM) v5 es un **evaluador geométrico** para auditar **estructuras reales** de proteínas en formato **PDB**.

- **No** ejecuta dinámica molecular.
- **No** usa entrenamiento estadístico.
- **No** genera datos simulados.
- Trabaja únicamente con **coordenadas atómicas** presentes en el archivo PDB (ATOM y, opcionalmente, HETATM).

## Métricas (lectura directa)
- **Σ (Sigma):** coherencia/topología de contactos **Cα–Cα** (no vecinos) con un cutoff (Å).
- **φ (Phi):** fricción material por **overlap VdW** entre átomos pesados:
  - overlap = (r_i + r_j)·vdw_scale − d_ij
  - aporta a φ si overlap > overlap_threshold
- **H_proxy:** proxy geométrico (log(Rg+1)) para trazabilidad, *no* entropía termodinámica.
- **Overlaps:** conteo de pares con overlap (y máximo overlap en Å).

## Ley de balance (auditoría)
Se evalúa el criterio:
**Q_eff · Σ > φ**

En esta versión, **Q_eff** es constante configurable (por defecto 1.0) y sirve como escala de comparación.

## Cómo correr (Colab o local)
Ejemplo con un PDB real:

```bash
python omnikernel_A_evaluator_v5.py --brain brain_evaluator_v5_vdw.json --pdbA 1CRN.pdb --out out
```

Comparación A vs B (solo si son conmensurables y deseas reportar ambos estados):

```bash
python omnikernel_A_evaluator_v5.py --brain brain_evaluator_v5_vdw.json --pdbA A.pdb --pdbB B.pdb --out out
```

Salida:
- `out/omnikernel_A_eval_v5_<config_hash>.json`

## Viewer HTML
Incluye dos visores:
- `omnikernel_protein_viewer_nosim_v5.html` (HUD + controles)
- `omnikernel_protein_metrics_only.html` (solo HUD/métricas)

> Nota: si el visor intenta descargar un PDB desde RCSB, requiere internet en el entorno donde abras el HTML.

## Autoría
**Genaro Carrasco Ozuna**  
ORCID: https://orcid.org/0009-0005-6358-9910
