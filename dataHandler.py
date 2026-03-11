import requests
import torch
import data.corespondance_Moecule_ID as corr

kb = 1.380_649e-23
c  = 299_792_458
Na = 6.02214129e23


class HITRANStructure:
    """
    Stores one line of data from the HITRAN API.
    Unchanged from the original — pure Python data container.
    """
    def __init__(self, molecule_id, local_iso_id, nu, sw, a, gamma_air, gamma_self,
                 elower, n_air, delta_air, global_upper_quanta, global_lower_quanta,
                 local_upper_quanta, local_lower_quanta, ierr, iref,
                 line_mixing_flag, gp, gpp):
        self.molecule_id          = molecule_id
        self.local_iso_id         = local_iso_id
        self.nu                   = nu
        self.sw                   = sw
        self.a                    = a
        self.gamma_air            = gamma_air
        self.gamma_self           = gamma_self
        self.elower               = elower
        self.n_air                = n_air
        self.delta_air            = delta_air
        self.global_upper_quanta  = global_upper_quanta
        self.global_lower_quanta  = global_lower_quanta
        self.local_upper_quanta   = local_upper_quanta
        self.local_lower_quanta   = local_lower_quanta
        self.ierr                 = ierr
        self.iref                 = iref
        self.line_mixing_flag     = line_mixing_flag
        self.gp                   = gp
        self.gpp                  = gpp

    def __str__(self):
        res = ''
        for k in self.__dict__.keys():
            res += k + ' = ' + str(self.__dict__[k]) + '\n'
        return res


class HITRANData(list):
    """
    High-level class for handling HITRAN data (list of HITRANStructure).

    After calling to_torch(), spectroscopic arrays are stored as torch tensors
    on the chosen device/dtype so that getTau() never touches CPU when on GPU.
    """

    def __init__(self, specie, isIsotpologueData=False):
        list.__init__(self)
        self.isIsotpologueData = isIsotpologueData
        self.specie            = specie

        # torch tensors — populated by to_torch()
        self._gamma_air  = None
        self._gamma_self = None
        self._n_air      = None
        self._delta_air  = None
        self._nui        = None
        self._Sw         = None

    # ------------------------------------------------------------------
    # Filtering helpers (unchanged API)
    # ------------------------------------------------------------------

    def getIsotopologueByID(self, id):
        res = HITRANData(self.specie, isIsotpologueData=True)
        for line in self:
            if line.local_iso_id == id:
                res.append(line)
        return res

    def getData(self, data):
        return [getattr(line, data) for line in self]

    # ------------------------------------------------------------------
    # Torch conversion
    # ------------------------------------------------------------------

    def to_torch(self, device=None, dtype=None):
        """
        Convert spectroscopic arrays to torch tensors and move them to
        *device* with *dtype*.  Must be called before getTau().

        Parameters
        ----------
        device : torch.device or str, optional  (default: cpu)
        dtype  : torch.dtype, optional          (default: torch.float64)

        Returns
        -------
        self  (for chaining)
        """
        if not self.isIsotpologueData:
            raise ValueError('to_torch() can only be called on a single-isotopologue HITRANData.')

        device = device if device is not None else torch.device('cpu')
        dtype  = dtype  if dtype  is not None else torch.float64

        def _t(lst):
            return torch.tensor(lst, dtype=dtype, device=device)

        self._gamma_air  = _t(self.getData('gamma_air'))
        self._gamma_self = _t(self.getData('gamma_self'))
        self._n_air      = _t(self.getData('n_air'))
        self._delta_air  = _t(self.getData('delta_air'))
        self._nui        = _t(self.getData('nu'))
        self._Sw         = _t(self.getData('sw'))

        self._device = device
        self._dtype  = dtype
        return self

    # ------------------------------------------------------------------
    # getTau — fully vectorised, no Python loops
    # ------------------------------------------------------------------

    def getTau(self, nu: torch.Tensor, T: float, P: float, specieFraction: float) -> torch.Tensor:
        """
        Computes optical depth per cm for a batch of wavenumbers.

        All operations are performed on the tensor device — no CPU↔GPU
        transfers occur when to_torch() was called with a GPU device.

        Parameters
        ----------
        nu             : torch.Tensor, shape (N,)  [cm-1]   — already on device
        T              : float  [K]
        P              : float  [atm]
        specieFraction : float  mol fraction of the species

        Returns
        -------
        tau : torch.Tensor, shape (N,)  [cm-1]
        """
        if self._nui is None:
            raise RuntimeError('Call to_torch() before getTau().')

        T0 = 296.0   # reference temperature [K]

        # --- Lorentzian half-width (n_lines,) ---
        gamma = ((T0 / T) ** self._n_air
                 * (self._gamma_air  * P * (1.0 - specieFraction)
                    + self._gamma_self * P * specieFraction))            # (n_lines,)

        # --- Shifted line centres (n_lines,) ---
        nu_shifted = self._nui + self._delta_air * P                     # (n_lines,)

        # --- Lorentzian profile  (N_nu, n_lines) — broadcast ---
        # nu          : (N, 1)
        # nu_shifted  : (1, n_lines)   via unsqueeze
        delta_nu = nu.unsqueeze(1) - nu_shifted.unsqueeze(0)            # (N, n_lines)
        fL = (1.0 / torch.pi) * gamma / (gamma ** 2 + delta_nu ** 2)   # (N, n_lines)

        # --- Number density factor ---
        # n [molecules/cm3] = P_pa / (kb * T)  ×  (1e-2)^3
        P_pa  = P * 101_325.0
        n_mol = (P_pa / (kb * T)) * (1e-2) ** 3                         # scalar

        # --- Sum over lines ---
        tau = (fL @ self._Sw) * n_mol * specieFraction                  # (N,)

        return tau

    # ------------------------------------------------------------------
    # Visualisation helper (unchanged API)
    # ------------------------------------------------------------------

    def plotData(self, y='linear'):
        import matplotlib.pyplot as plt
        if not self.isIsotpologueData:
            raise ValueError('Impossible to plot multi-isotopologue data.')
        nu = self.getData('nu')
        Sw = self.getData('sw')
        plt.plot(nu, Sw)
        plt.yscale(y)
        plt.title(f'Absorption coefficient of {self.specie}')
        plt.xlabel('Frequency [cm-1]')
        plt.ylabel('Absorption coefficient [cm-1 / (molecule * cm-2)]')


class DataHandler:
    """
    Handles all operations related to stored data: downloading and reading.
    Unchanged public API.
    """

    API_FORMAT = 'https://hitran.org/lbl/api?iso_ids_list={}&numin={}&numax={}'

    def download(self, specie, nu_min, nu_max):
        ids = []
        for val in corr.ISO.values():
            if val[-1] == specie:
                ids.append(val[0])
        if len(ids) == 0:
            raise ValueError(f'No data found for specie {specie}')
        iso_ids_list = ','.join(str(i) for i in ids)
        data = requests.get(DataHandler.API_FORMAT.format(iso_ids_list, nu_min, nu_max)).text
        with open(f'data/bandes_{specie}.par', 'w') as file:
            file.write(data)

    def readData(self, specie):
        res = HITRANData(specie)
        with open(f'data/bandes_{specie}.par', 'r') as file:
            data = file.readlines()
        for line in data:
            res.append(HITRANStructure(
                molecule_id           = int(line[0:2],  36),
                local_iso_id          = int(line[2:3],  36),
                nu                    = float(line[3:15]),
                sw                    = float(line[15:25]),
                a                     = float(line[25:35]),
                gamma_air             = float(line[35:40]),
                gamma_self            = float(line[40:45]),
                elower                = float(line[45:55]),
                n_air                 = float(line[55:59]),
                delta_air             = float(line[59:67]),
                global_upper_quanta   = line[67:82],
                global_lower_quanta   = line[82:97],
                local_upper_quanta    = line[97:112],
                local_lower_quanta    = line[112:127],
                ierr                  = line[127:128],
                iref                  = line[128:130],
                line_mixing_flag      = line[130:131],
                gp                    = line[131:138],
                gpp                   = line[138:145],
            ))
        return res


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    a = DataHandler()

    #a.download('CO2', 200, 1600)
    #a.download('H2O', 200, 1600)
    #a.download('CH4', 200, 1600)

    CO2 = a.readData('CO2')
    CO2_iso_1 = CO2.getIsotopologueByID(1)
    for l in CO2_iso_1:
        print(l.local_iso_id)
    CO2_iso_1.plotData('log')
    plt.show()

    H2O = a.readData('H2O')
    H2O_iso_1 = H2O.getIsotopologueByID(1)
    for l in H2O_iso_1:
        print(l.local_iso_id)
    H2O_iso_1.plotData('log')
    plt.show()

    H2O = a.readData('CH4')
    H2O_iso_1 = H2O.getIsotopologueByID(1)
    for l in H2O_iso_1:
        print(l.local_iso_id)
    H2O_iso_1.plotData('log')
    plt.show()

    nu = np.linspace(200,1600,10000)
    L = CO2_iso_1.getTau(nu, 300, 1, 420e-6)  # cm-1
    plt.plot(nu, L)
    plt.yscale('log')
    plt.show()