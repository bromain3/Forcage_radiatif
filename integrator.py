"""
integrator.py — Intégration numérique

Ce module est responsable uniquement des schémas d'intégration :
  - Intégration scalaire (quadrature) pour valider des intégrales spectrales
  - Intégration de l'équation de transfert radiatif (RTE)

Il est complètement indépendant de la physique : il reçoit des fonctions
ou des tableaux en entrée et retourne des valeurs numériques.

Équation de transfert radiatif (1D, plan-parallèle, propagation verticale) :
  dI/dz = -alpha(z) * [I(z) - B(nu, T(z))]

Schéma d'intégration par couche homogène :
  tau   = alpha(z_mid) * dz_cm         [sans dimension]
  I_out = I_in * exp(-tau) + B(T_mid) * (1 - exp(-tau))

Ce schéma est exact pour une couche d'absorption constante. Il est
inconditionnellement stable (contrairement à Euler explicite) et garantit
physiquement 0 <= I_out <= max(I_in, B) pour tout tau >= 0.
"""

import numpy as np


class Integrator:
    """
    Intégrateurs scalaires (quadrature 1D sur tableaux).

    Ces méthodes opèrent sur des tableaux X, Y représentant une fonction
    échantillonnée — elles n'ont pas accès à la fonction analytique.
    """

    def integrateSquare(self, X, Y):
        """Méthode des rectangles (ordre 1). Borne gauche."""
        dX = np.diff(X)
        return float(np.sum(dX * Y[:-1]))

    def integrateTrapezoidal(self, X, Y):
        """Méthode des trapèzes (ordre 2)."""
        dX = np.diff(X)
        midY = 0.5 * (Y[:-1] + Y[1:])
        return float(np.sum(dX * midY))

    def integrateModifiedSimpson(self, X, Y):
        """
        Règle de Simpson par intervalles avec point milieu interpolé (ordre 2).
        Équivalent aux trapèzes pour une interpolation linéaire entre points.
        Pour un gain réel d'ordre, il faudrait passer la fonction analytique.
        """
        dX = np.diff(X)
        midY = 0.5 * (Y[:-1] + Y[1:])
        return float(np.sum(dX / 6.0 * (Y[:-1] + 4.0 * midY + Y[1:])))

    def integrateMonteCarlo(self, X, Y, niter=10_000):
        """
        Intégration Monte-Carlo par interpolation linéaire sur [X[0], X[-1]].
        Utile comme référence stochastique ou pour des intégrandes irréguliers.

        Paramètres
        ----------
        niter : int, nombre de tirages aléatoires.
        """
        x_rand = np.random.uniform(X[0], X[-1], niter)
        y_rand = np.interp(x_rand, X, Y)
        return float((X[-1] - X[0]) * np.mean(y_rand))


class RTEIntegrator:
    """
    Intégrateur de l'équation de transfert radiatif (RTE).

    Propage un spectre I(nu) depuis z=0 jusqu'à z_max en appliquant
    le schéma exact couche homogène à chaque pas dz.

    Le schéma est physiquement motivé :
      - Si tau >> 1 (couche opaque) : I_out -> B(T)  (émission locale)
      - Si tau << 1 (couche transparente) : I_out -> I_in  (transmission)
    Il est donc correct dans les deux limites sans ajustement numérique.

    Paramètres
    ----------
    solver : instance de Solver (optique)
        Fournit solver.alpha(nu, z) et solver.planck(nu, T)
    """

    def __init__(self, solver):
        self.solver = solver

    def propagate(self, nu, dz_m, z_max, verbose=True):
        """
        Propage le rayonnement de z=0 à z_max.

        Paramètres
        ----------
        nu     : np.array  [cm-1]   grille spectrale
        dz_m   : float     [m]      pas vertical
        z_max  : float     [m]      altitude maximale
        verbose: bool               affiche la progression

        Retour
        ------
        I_surface : np.array  [W/m²/sr/cm-1]  spectre corps noir à z=0
        I_top     : np.array  [W/m²/sr/cm-1]  spectre au sommet z_max
        z_steps   : np.array  [m]              altitudes des couches
        """
        dz_cm   = dz_m * 100.0          # conversion m -> cm (cohérence avec alpha [cm-1])
        n_steps = int(z_max / dz_m)
        z_steps = np.arange(n_steps) * dz_m

        # Condition initiale : corps noir à la surface
        T_surf    = self.solver.surface_temperature()
        I_surface = self.solver.planck(nu, T_surf)
        I         = I_surface.copy()

        progress_interval = max(1, n_steps // 20)

        for step in range(n_steps):
            # Altitude du centre de la couche (convention milieu)
            z_mid = (step + 0.5) * dz_m

            # Propriétés optiques et thermiques au centre de la couche
            T_mid, P_pa, _, _ = self.solver.atm.getAtmProperties(z_mid)
            alpha             = self.solver.alpha(nu, z_mid)

            # Épaisseur optique de la couche
            tau   = alpha * dz_cm        # [sans dimension]

            # Luminance d'équilibre local
            B     = self.solver.planck(nu, T_mid)

            # Schéma exact couche homogène
            e_tau = np.exp(-tau)
            I     = I * e_tau + B * (1.0 - e_tau)

            if verbose and step % progress_interval == 0:
                pct = 100.0 * step / n_steps
                print(f"    z = {step * dz_m:7.1f} m  ({pct:5.1f}%)")

        if verbose:
            print(f"    z = {z_max:7.1f} m  (100.0%)")

        return I_surface, I, z_steps