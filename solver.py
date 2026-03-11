"""
solver.py — Calculs optiques purs (PyTorch)

Ce module est responsable uniquement de la physique spectroscopique :
  - Loi de Planck
  - Coefficient d'absorption total alpha(nu, z) pour N espèces

Toutes les opérations numériques utilisent PyTorch.
Le device et le dtype se propagent à tous les éléments internes :
  - le tenseur nu d'entrée doit déjà être sur le bon device
  - les tenseurs HITRAN sont convertis lors de l'initialisation
  - aucun transfert CPU↔GPU n'a lieu dans alpha() ou planck()

Unités cohérentes avec dataHandler.getTau :
  - nu      [cm-1]
  - alpha   [cm-1]           (coefficient d'absorption volumique)
  - planck  [W/m²/sr/cm-1]  (luminance corps noir par nombre d'onde)
"""

import torch
import atmosphere_pur_python as model_atmo
import dataHandler as dh

# ============================================================
# CONSTANTES PHYSIQUES
# ============================================================
H_PLANCK = 6.62607015e-34   # constante de Planck  [J.s]
C_LIGHT  = 299_792_458.0    # vitesse de la lumière [m/s]
KB       = 1.380649e-23     # constante de Boltzmann [J/K]
ONE_ATM  = 101_325.0        # 1 atmosphère [Pa]


class Solver:
    """
    Classe de calcul optique — PyTorch.

    Paramètres
    ----------
    species      : list[str]        espèces radiativement actives
    atm          : StandardAtm      modèle atmosphérique (Python pur, inchangé)
    data_handler : DataHandler      gestionnaire HITRAN
    device       : torch.device     device cible (cpu / cuda / mps …)
    dtype        : torch.dtype      précision flottante (float32 / float64)
    """

    def __init__(self, species, atm=None, data_handler=None,
                 device=None, dtype=None):
        self.species      = species
        self.atm          = atm          if atm          is not None else model_atmo.StandardAtm()
        self.data_handler = data_handler if data_handler is not None else dh.DataHandler()
        self.device       = device if device is not None else torch.device('cpu')
        self.dtype        = dtype  if dtype  is not None else torch.float64

        # Chargement + conversion torch des données HITRAN
        self.species_data = {}
        for sp in self.species:
            raw  = self.data_handler.readData(sp)
            iso1 = raw.getIsotopologueByID(1)
            iso1.to_torch(device=self.device, dtype=self.dtype)   # ← conversion ici
            self.species_data[sp] = iso1
            print(f"  [Solver] {sp} : {len(iso1)} raies chargées "
                  f"(isotopologue 1, device={self.device}, dtype={self.dtype})")

    # ----------------------------------------------------------
    # LOI DE PLANCK
    # ----------------------------------------------------------

    def planck(self, nu: torch.Tensor, T: float) -> torch.Tensor:
        """
        Luminance spectrale du corps noir.

        B(nu, T) = 2 h c² × 1e8 × nu³ / (exp(x) − 1)
        avec  x = h c × 100 × nu / (kb T)

        Paramètres
        ----------
        nu : torch.Tensor [cm-1]  — doit être sur self.device
        T  : float [K]

        Retour
        ------
        B  : torch.Tensor, même shape que nu  [W/m²/sr/cm-1]
        """
        x  = H_PLANCK * C_LIGHT * 100.0 * nu / (KB * T)
        # Garde-fou overflow : on masque les nu nuls et les x trop grands
        safe = (nu > 0.0) & (x < 700.0)
        B    = torch.zeros_like(nu)
        xs   = x[safe]
        B[safe] = (2.0 * H_PLANCK * C_LIGHT**2 * 1e8
                   * nu[safe]**3 / (torch.exp(xs) - 1.0))
        return B

    # ----------------------------------------------------------
    # COEFFICIENT D'ABSORPTION TOTAL
    # ----------------------------------------------------------

    def alpha(self, nu: torch.Tensor, z: float) -> torch.Tensor:
        """
        Coefficient d'absorption volumique total alpha(nu) à l'altitude z.

        Les T, P, concentrations sont récupérés depuis le modèle atmosphérique
        (Python pur, scalaires), puis getTau() est appelé pour chaque espèce —
        tout se passe sur self.device, pas de transfert CPU↔GPU.

        Paramètres
        ----------
        nu : torch.Tensor [cm-1]  — doit être sur self.device
        z  : float, altitude géométrique [m]

        Retour
        ------
        alpha_total : torch.Tensor, même shape que nu  [cm-1]
        """
        T, P_pa, _, species_conc = self.atm.getAtmProperties(z)
        P_atm = P_pa / ONE_ATM

        alpha_total = torch.zeros(len(nu), dtype=self.dtype, device=self.device)

        for sp, data in self.species_data.items():
            x_sp = species_conc.get(sp, 0.0)
            if x_sp == 0.0:
                continue
            alpha_total += data.getTau(nu, T, P_atm, x_sp)

        return alpha_total

    # ----------------------------------------------------------
    # UTILITAIRES
    # ----------------------------------------------------------

    def surface_temperature(self) -> float:
        """Retourne la température de surface [K]."""
        T, _, _, _ = self.atm.getAtmProperties(0.0)
        return T

    def make_nu(self, nu_min: float, nu_max: float, n_nu: int) -> torch.Tensor:
        """
        Crée la grille spectrale directement sur self.device / self.dtype.
        Permet d'éviter tout transfert CPU→GPU dans main.py.

        Retour
        ------
        nu : torch.Tensor, shape (n_nu,)  [cm-1]  sur self.device
        """
        return torch.linspace(nu_min, nu_max, n_nu,
                              dtype=self.dtype, device=self.device)