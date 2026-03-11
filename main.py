"""
main.py — Simulation de transfert radiatif atmosphérique

Physique implémentée :
  - Loi de Planck (corps noir)
  - Profil d'absorption Lorentzien (élargissement par pression, données HITRAN)
  - Équation de transfert radiatif 1D plan-parallèle
  - Modèle atmosphérique US Standard 1976 avec profils de concentration
    variables par espèce (définis dans atmosphere_pur_python.py)

Architecture :
  solver.py      → optique pure (Planck, alpha multi-espèces)
  integrator.py  → schéma numérique RTE + quadratures spectrales
  main.py        → paramètres, orchestration, visualisation

Compilation équivalente C :
  gcc -O2 -o atmosphere main.c -lm && ./atmosphere data_CO2.par
"""

import numpy as np
import matplotlib.pyplot as plt

import atmosphere_pur_python as model_atmo
import dataHandler as dh
from solver import Solver
from integrator import Integrator, RTEIntegrator

# ============================================================
# PARAMÈTRES NUMÉRIQUES  ← ajuster ici
# ============================================================

# --- Grille spectrale ---
NU_MIN = 200.0     # nombre d'onde minimal   [cm-1]
NU_MAX = 1600.0    # nombre d'onde maximal   [cm-1]
N_NU   = 1000      # nombre de points spectraux
#   ↑ augmenter pour résoudre les raies fines (ex. 10 000)
#   ↓ diminuer pour un test rapide

# --- Discrétisation verticale ---
DZ_M   = 100.0      # pas vertical             [m]
#   Recommandation : 10 m (précis) à 100 m (rapide)
#   Le schéma exact est stable pour tout DZ — pas de contrainte CFL.
#   En revanche, un DZ trop grand rate les gradients de T et de concentration.

Z_MAX  = 15_000.0  # altitude maximale        [m]
#   15 000 m couvre la troposphère (essentiel pour l'effet de serre CO2/H2O)
#   86 000 m pour la mésosphère complète (très lent)

# --- Espèces radiativement actives ---
# Chaque espèce doit avoir :
#   1. Un fichier  data/bandes_{espèce}.par   (données HITRAN)
#   2. Une méthode C{espèce}(h) dans atmosphere_pur_python.StandardAtm
#      retournant la fraction molaire à l'altitude h
#   3. Une entrée dans le dict speciesConcentration de getAtmProperties()
#
# Espèces disponibles selon atmosphere_pur_python.py :
#   'CO2'  →  420 ppm constant
#   'H2O'  →  1 % constant (profil exponentiel commenté dans atm)
#   'CH4'  →  1.9 ppm constant
#   'N2O'  →  0.33 ppm constant
#
# Pour ajouter une espèce : décommenter/modifier le profil dans atm,
# ajouter la clé dans speciesConcentration, et télécharger le .par HITRAN.

#SPECIES = ['CO2']
SPECIES = ['CO2', 'H2O']          # ajouter H2O (effet dominant dans l'IR)
#SPECIES = ['CO2', 'H2O', 'CH4', 'N2O']  # simulation complète

# --- Sortie ---
OUTPUT_FILE = 'spectre_radiatif.dat'

# ============================================================
# INITIALISATION
# ============================================================

print("=== Transfert radiatif atmosphérique ===\n")
print(f"  Espèces : {SPECIES}")
print(f"  Grille  : {N_NU} points, nu ∈ [{NU_MIN}, {NU_MAX}] cm-1")
print(f"  Vertical: dz = {DZ_M} m, z_max = {Z_MAX} m "
      f"({int(Z_MAX / DZ_M)} couches)\n")

atm          = model_atmo.StandardAtm()
data_handler = dh.DataHandler()

nu = np.linspace(NU_MIN, NU_MAX, N_NU)

# Instanciation du solveur optique (chargement HITRAN)
print("[1/3] Chargement des données HITRAN...")
sol = Solver(SPECIES, atm=atm, data_handler=data_handler)

# Instanciation de l'intégrateur RTE
rte = RTEIntegrator(sol)

# ============================================================
# PROPAGATION RADIATIVE
# ============================================================

print(f"\n[2/3] Propagation RTE de z=0 à z={Z_MAX:.0f} m par pas de {DZ_M} m...")
I_surface, I_top, z_steps = rte.propagate(nu, DZ_M, Z_MAX, verbose=True)

# ============================================================
# INTÉGRALES SPECTRALES (validation / analyse)
# ============================================================

print("\n[3/3] Intégrales spectrales...")
integrator = Integrator()

# Flux total émis (Stefan-Boltzmann attendu si corps noir pur)
T_surf       = sol.surface_temperature()
flux_surface = np.pi * integrator.integrateTrapezoidal(nu, I_surface)
flux_top     = np.pi * integrator.integrateTrapezoidal(nu, I_top)
flux_SB      = 5.670374419e-8 * T_surf**4   # W/m² (Stefan-Boltzmann)

print(f"  T_surface           = {T_surf:.2f} K")
print(f"  Flux surface (π∫B)  = {flux_surface:.2f} W/m²   (Stefan-Boltzmann = {flux_SB:.2f} W/m²)")
print(f"  Flux sommet  (π∫I)  = {flux_top:.2f} W/m²")
print(f"  Réduction           = {100*(1 - flux_top/flux_surface):.1f}%  "
      f"(énergie retenue par l'atmosphère)")

# ============================================================
# SAUVEGARDE
# ============================================================

header = (f"Transfert radiatif — espèces: {SPECIES}, dz={DZ_M}m, z_max={Z_MAX}m\n"
          f"nu(cm-1)   I_surface(W/m2/sr/cm-1)   I_top(W/m2/sr/cm-1)")
np.savetxt(OUTPUT_FILE,
           np.column_stack([nu, I_surface, I_top]),
           header=header, fmt='%.4f  %.6e  %.6e')
print(f"\n  Données sauvegardées : '{OUTPUT_FILE}'")

# ============================================================
# VISUALISATION
# ============================================================

fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
fig.suptitle(f"Transfert radiatif — espèces : {SPECIES} — dz = {DZ_M} m, z_max = {Z_MAX/1000:.0f} km",
             fontsize=12)

# Spectre surface vs sommet
axes[0].plot(nu, I_surface, color='tomato',    lw=1.0, label=f'Surface — corps noir {T_surf:.1f} K')
axes[0].plot(nu, I_top,     color='steelblue', lw=1.0, label=f'Sommet atmosphère ({Z_MAX/1000:.0f} km)')
axes[0].set_ylabel('Luminance [W/m²/sr/cm⁻¹]')
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)
axes[0].set_title("Spectre radiatif")

# Flux absorbé par l'atmosphère (effet de serre)
delta = I_surface - I_top
axes[1].fill_between(nu, delta, color='mediumseagreen', alpha=0.5, label='Surface − Sommet')
axes[1].plot(nu, delta, color='green', lw=0.7)
axes[1].set_ylabel('ΔLuminance [W/m²/sr/cm⁻¹]')
axes[1].set_xlabel('Nombre d\'onde [cm⁻¹]')
axes[1].set_title("Flux retenu par l'atmosphère (effet de serre)")
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('spectre_radiatif.png', dpi=150)
print("  Graphique sauvegardé : 'spectre_radiatif.png'")
plt.show()