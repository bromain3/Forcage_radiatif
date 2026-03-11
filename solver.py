"""
solver.py — Calculs optiques purs

Ce module est responsable uniquement de la physique spectroscopique :
  - Loi de Planck
  - Coefficient d'absorption total alpha(nu, z) pour N espèces

Il ne contient aucune logique d'intégration numérique.
La séparation est volontaire : changer le schéma d'intégration ne touche
pas ce fichier, et vice-versa.

Unités cohérentes avec dataHandler.getTau :
  - nu      [cm-1]
  - alpha   [cm-1]           (coefficient d'absorption volumique)
  - planck  [W/m²/sr/cm-1]  (luminance corps noir par nombre d'onde)
"""

import numpy as np
import atmosphere_pur_python as model_atmo
import dataHandler as dh

# ============================================================
# CONSTANTES PHYSIQUES
# ============================================================
H_PLANCK = 6.62607015e-34   # constante de Planck        [J.s]
C_LIGHT  = 299_792_458.0    # vitesse de la lumière       [m/s]
KB       = 1.380649e-23     # constante de Boltzmann      [J/K]
ONE_ATM  = 101_325.0        # 1 atmosphère                [Pa]


class Solver:
    """
    Classe de calcul optique.

    Responsabilités :
      - Calcul de la luminance de Planck B(nu, T)
      - Calcul du coefficient d'absorption volumique alpha(nu, z)
        pour une liste arbitraire d'espèces dont les profils de
        concentration sont fournis par le modèle d'atmosphère.

    Le calcul d'alpha délègue à HITRANData.getTau() de dataHandler.py
    pour le profil Lorentzien. getTau retourne déjà alpha [cm-1]
    en tenant compte de la densité numérique locale de chaque espèce.
    """

    def __init__(self, species, atm=None, data_handler=None):
        """
        Paramètres
        ----------
        species : list[str]
            Liste des espèces à prendre en compte, ex. ['CO2', 'H2O', 'CH4'].
            Chaque espèce doit avoir un fichier data/bandes_{espèce}.par
            et être référencée dans atmosphere_pur_python.getAtmProperties().
        atm : StandardAtm, optionnel
            Instance du modèle atmosphérique (permet d'injecter un modèle
            personnalisé avec des profils de concentration différents).
        data_handler : DataHandler, optionnel
            Instance du gestionnaire de données HITRAN.
        """
        self.species      = species
        self.atm          = atm          if atm          is not None else model_atmo.StandardAtm()
        self.data_handler = data_handler if data_handler is not None else dh.DataHandler()

        # Chargement des données HITRAN : isotopologue principal (id=1) pour chaque espèce
        self.species_data = {}
        for sp in self.species:
            raw = self.data_handler.readData(sp)
            self.species_data[sp] = raw.getIsotopologueByID(1)
            print(f"  [Solver] {sp} : {len(self.species_data[sp])} raies chargées (isotopologue 1)")

    # ----------------------------------------------------------
    # LOI DE PLANCK
    # ----------------------------------------------------------

    def planck(self, nu, T):
        """
        Luminance spectrale du corps noir en nombre d'onde.

        B(nu, T) = 2 h c² × 100^4 × nu³ / (exp(x) - 1)
        avec  x = h c (100 nu) / (kb T)

        Le facteur 100^4 provient de la conversion de l'unité spectrale :
          dnu [cm-1] = dnu_SI [m-1] / 100
          => B [W/m²/sr/cm-1] = B_SI [W/m²/sr/m-1] × 100
        => préfacteur global = 2 h c² × 100^3 × 100 = 2 h c² × 1e8

        Paramètres
        ----------
        nu : float ou np.array  [cm-1]
        T  : float              [K]

        Retour
        ------
        B  : même shape que nu  [W/m²/sr/cm-1]
        """
        nu = np.asarray(nu, dtype=float)
        x  = H_PLANCK * C_LIGHT * 100.0 * nu / (KB * T)
        B  = np.zeros_like(nu)
        ok = (nu > 0.0) & (x < 700.0)          # garde-fou overflow
        B[ok] = (2.0 * H_PLANCK * C_LIGHT**2 * 1e8
                 * nu[ok]**3 / (np.exp(x[ok]) - 1.0))
        return B

    # ----------------------------------------------------------
    # COEFFICIENT D'ABSORPTION TOTAL
    # ----------------------------------------------------------

    def alpha(self, nu, z):
        """
        Coefficient d'absorption volumique total alpha(nu) à l'altitude z.

        Pour chaque espèce configurée, on récupère :
          - la fraction molaire x_sp(z) depuis le modèle atmosphérique
            (profil potentiellement variable avec l'altitude)
          - T(z) et P(z) depuis le modèle atmosphérique
        puis on appelle getTau(nu, T, P_atm, x_sp) de dataHandler,
        qui retourne alpha_sp [cm-1] via le profil Lorentzien.

        La somme sur les espèces donne l'absorption totale.

        Paramètres
        ----------
        nu : np.array  [cm-1]
        z  : float, altitude géométrique  [m]

        Retour
        ------
        alpha_total : np.array, même shape que nu  [cm-1]
        """
        T, P_pa, _, species_conc = self.atm.getAtmProperties(z)
        P_atm = P_pa / ONE_ATM

        alpha_total = np.zeros(len(nu), dtype=float)

        for sp, data in self.species_data.items():
            x_sp = species_conc.get(sp, 0.0)
            if x_sp == 0.0:
                continue
            alpha_sp = data.getTau(nu, T, P_atm, x_sp)
            alpha_total += alpha_sp

        return alpha_total

    def surface_temperature(self):
        """Retourne la température de surface [K] selon le modèle atmosphérique."""
        T, _, _, _ = self.atm.getAtmProperties(0.0)
        return T