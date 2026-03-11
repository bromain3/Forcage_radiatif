"""
main.py — Simulation de transfert radiatif atmosphérique (PyTorch)

Physique implémentée :
  - Loi de Planck (corps noir)
  - Profil d'absorption Lorentzien (élargissement par pression, données HITRAN)
  - Équation de transfert radiatif 1D plan-parallèle
  - Modèle atmosphérique US Standard 1976 avec profils de concentration
    variables par espèce (définis dans atmosphere_pur_python.py)

Architecture :
  solver.py      → optique pure (Planck, alpha multi-espèces) — PyTorch
  integrator.py  → schéma numérique RTE + quadratures spectrales — PyTorch
  dataHandler.py → lecture HITRAN + getTau vectorisé — PyTorch
  main.py        → paramètres, orchestration, visualisation

Choix du device / dtype :
  Modifier DEVICE et DTYPE ci-dessous. Le choix se propage automatiquement
  à tous les modules (Solver, HITRANData, RTEIntegrator).
  Sur GPU, aucun transfert CPU↔GPU n'a lieu pendant la propagation RTE.
"""

import torch
import numpy as np
import matplotlib.pyplot as plt

import atmosphere_pur_python as model_atmo
import dataHandler as dh
from solver import Solver
from integrator import Integrator, RTEIntegrator

# ============================================================
# PARAMÈTRES DEVICE / DTYPE  ← choisir ici
# ============================================================

DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# DEVICE = torch.device('cpu')
# DEVICE = torch.device('mps')   # Apple Silicon

DTYPE  = torch.float64
# DTYPE = torch.float32          # plus rapide, moins précis

# ============================================================
# PARAMÈTRES NUMÉRIQUES  ← ajuster ici
# ============================================================

# --- Grille spectrale ---
NU_MIN = 200.0
NU_MAX = 1600.0
N_NU   = 10000

# --- Discrétisation verticale ---
DZ_M   = 100.0    # pas vertical [m]
Z_MAX  = 15_000.0 # altitude maximale [m]

# --- Espèces radiativement actives ---
SPECIES = ['CO2']

# --- Sortie ---
OUTPUT_FILE = 'spectre_radiatif.dat'

# ============================================================
# INITIALISATION
# ============================================================

print("=== Transfert radiatif atmosphérique (PyTorch) ===\n")
print(f"  Device  : {DEVICE}")
print(f"  Dtype   : {DTYPE}")
print(f"  Espèces : {SPECIES}")
print(f"  Grille  : {N_NU} points, nu ∈ [{NU_MIN}, {NU_MAX}] cm-1")
print(f"  Vertical: dz = {DZ_M} m, z_max = {Z_MAX} m "
      f"({int(Z_MAX / DZ_M)} couches)\n")

atm          = model_atmo.StandardAtm()
data_handler = dh.DataHandler()

# Instanciation du solveur (charge et convertit les données HITRAN sur device)
print("[1/3] Chargement des données HITRAN...")
sol = Solver(SPECIES, atm=atm, data_handler=data_handler,
             device=DEVICE, dtype=DTYPE)

# Grille spectrale créée directement sur le device — pas de transfert CPU→GPU
nu = sol.make_nu(NU_MIN, NU_MAX, N_NU)

# Instanciation de l'intégrateur RTE
rte = RTEIntegrator(sol)

# ============================================================
# PROPAGATION RADIATIVE
# ============================================================

print(f"\n[2/3] Propagation RTE de z=0 à z={Z_MAX:.0f} m par pas de {DZ_M} m...")
I_surface, I_top, z_steps = rte.propagate(nu, DZ_M, Z_MAX, verbose=True)

# ============================================================
# INTÉGRALES SPECTRALES
# ============================================================

print("\n[3/3] Intégrales spectrales...")
integrator = Integrator()

T_surf       = sol.surface_temperature()
flux_surface = torch.pi * integrator.integrateTrapezoidal(nu, I_surface)
flux_top     = torch.pi * integrator.integrateTrapezoidal(nu, I_top)
flux_SB      = 5.670374419e-8 * T_surf**4

print(f"  T_surface           = {T_surf:.2f} K")
print(f"  Flux surface (π∫B)  = {flux_surface:.2f} W/m²   (Stefan-Boltzmann = {flux_SB:.2f} W/m²)")
print(f"  Flux sommet  (π∫I)  = {flux_top:.2f} W/m²")
print(f"  Réduction           = {100*(1 - flux_top/flux_surface):.1f}%  "
      f"(énergie retenue par l'atmosphère)")

# ============================================================
# SAUVEGARDE  (rapatrier les tenseurs sur CPU pour numpy/savetxt)
# ============================================================

nu_np        = nu.cpu().numpy()
I_surface_np = I_surface.cpu().numpy()
I_top_np     = I_top.cpu().numpy()

header = (f"Transfert radiatif — espèces: {SPECIES}, dz={DZ_M}m, z_max={Z_MAX}m\n"
          f"nu(cm-1)   I_surface(W/m2/sr/cm-1)   I_top(W/m2/sr/cm-1)")
np.savetxt(OUTPUT_FILE,
           np.column_stack([nu_np, I_surface_np, I_top_np]),
           header=header, fmt='%.4f  %.6e  %.6e')
print(f"\n  Données sauvegardées : '{OUTPUT_FILE}'")

# ============================================================
# VISUALISATION
# ============================================================

fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
fig.suptitle(f"Transfert radiatif — espèces : {SPECIES} — dz = {DZ_M} m, z_max = {Z_MAX/1000:.0f} km",
             fontsize=12)

axes[0].plot(nu_np, I_surface_np, color='tomato',    lw=1.0, label=f'Surface — corps noir {T_surf:.1f} K')
axes[0].plot(nu_np, I_top_np,     color='steelblue', lw=1.0, label=f'Sommet atmosphère ({Z_MAX/1000:.0f} km)')
axes[0].set_ylabel('Luminance [W/m²/sr/cm⁻¹]')
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)
axes[0].set_title("Spectre radiatif")

delta = I_surface_np - I_top_np
axes[1].fill_between(nu_np, delta, color='mediumseagreen', alpha=0.5, label='Surface − Sommet')
axes[1].plot(nu_np, delta, color='green', lw=0.7)
axes[1].set_ylabel('ΔLuminance [W/m²/sr/cm⁻¹]')
axes[1].set_xlabel('Nombre d\'onde [cm⁻¹]')
axes[1].set_title("Flux retenu par l'atmosphère (effet de serre)")
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('spectre_radiatif.png', dpi=150)
print("  Graphique sauvegardé : 'spectre_radiatif.png'")
plt.show()