import requests  # in C, an equivalent in libcurl
import matplotlib.pyplot as plt  # used to plot data. Plotting data can be a bit harder in C and should be replaced by saving data at first
import data.corespondance_Moecule_ID as corr
import numpy as np
import scipy.signal as signal


kb = 1.380_649e-23
c = 299_792_458
Na = 6.02214129e23
class HITRANStructure:
    """
    Class storing one line of data provided by the HITRAN API
    """
    def __init__(self, molecule_id, local_iso_id, nu, sw, a, gamma_air, gamma_self, elower, n_air, delta_air, global_upper_quanta, global_lower_quanta, local_upper_quanta, local_lower_quanta, ierr, iref, line_mixing_flag, gp, gpp):
        """
        Constructor of the class. Each input corresponds to one line of data provided by the HITRAN API and is stored in the corresponding instance variable
        :param Download HITRAN data from the website and look at the README for more insight on the input data
        """
        self.molecule_id = molecule_id
        self.local_iso_id = local_iso_id
        self.nu = nu
        self.sw = sw
        self.a = a
        self.gamma_air = gamma_air
        self.gamma_self = gamma_self
        self.elower = elower
        self.n_air = n_air
        self.delta_air = delta_air
        self.global_upper_quanta = global_upper_quanta
        self.global_lower_quanta = global_lower_quanta
        self.local_upper_quanta = local_upper_quanta
        self.local_lower_quanta = local_lower_quanta
        self.ierr = ierr
        self.iref = iref
        self.line_mixing_flag = line_mixing_flag
        self.gp = gp
        self.gpp = gpp

    def __str__(self):
        """
        Overload of the str function. It is used to print instance variables of the class to be sure everything was stored as it should be.
        You cannot overload functions in pur C as it is not object-oriented (you can in C++ but it can be hard). However in C, reading delimited data like
        this is way easier, so verifying the output is not mandatory
        :return: string object representing the class
        """
        #  initialisation of the result string
        res = ''
        for i in self.__dict__.keys():  # looping over every instance variable stored in __dict__.keys()
            res += i + ' = ' + str(self.__dict__[i]) + '\n'
        return(res)

class HITRANData(list):
    #TODO: en vrai sw dépende de la tempéraature et de la pression. il serait interessant de le prendre en compte
    """
    High level class for handling HITRAN data. It inherits form list object. It actualy is a list of HITRANStructure objects
    """
    def __init__(self, specie, isIsotpologueData=False):
        """
        Constructor of the class HITRANData
        :param specie (str): Name of the specie represented by HITRANData object. It is mainly used for labeling graphs
        :param isIsotpologueData (bool): Flag to know if HITRANData contains multiple isotopologues or not. It should be understood that it IS NOT
            the result of a data analysis but it is just declared at class initialisation
        """
        list.__init__(self)
        self.isIsotpologueData = isIsotpologueData
        self.specie = specie

    def getIsotopologueByID(self, id):
        """
        Extracts data from an isotopologue given its local HITRAN id
        :param id (int): Local HITRAN id
        :return (HITRANData): Isotopologue data
        """
        res = HITRANData(self.specie, isIsotpologueData=True)
        for line in self:
            if line.local_iso_id == id:
                res.append(line)
        return(res)

    def getData(self, data):
        """
        Extract a given column from the data
        :param data (str): Name of the data to extract, as given in HITRAN headers
        :return (list): List of the values ordered as in the file containing data
        """
        dat = []
        for line in self:
            dat.append(getattr(line, data))
        return(dat)


    def getTau(self, nu, T, P, specieFraction):
        """
        Computes optical depth per cm given a list of frequencies at which to compute values
        It can be computationaly optimized by reducing the number of large array operations

        This function neglects variation of line intensity with T.
        It only considers Lorentzian profile (pressure-broadening). For low pressure, a dopler-broadened (Gaussian) profile is more accurate

        :param nu (np.array, one dimensional): List of wavenumber at which to compute values (cm-1)
        :param T (float): Temperature (K)
        :param P (float): Pressure (atm)
        :param specieFraction (float): Specie mol fraction
        :return (float): Absorption coefficient
        """

        if not self.isIsotpologueData:
            raise ValueError('Impossible to compute absorption coefficient for all isotopologues at once in this version')


        T0 = 296  # reference temperature
        P0 = 1.0  # reference pressure

        gamma_air = np.array(self.getData('gamma_air'))
        gamma_self = np.array(self.getData('gamma_self'))
        n_air = np.array(self.getData('n_air'))
        delta = np.array(self.getData('delta_air'))
        nui = np.array(self.getData('nu'))
        Sw = np.array(self.getData('sw'))

        # limitaiton : constant line intensity with T
        n = n_air # or 0.5
        gamma = (T0 / T) ** n * (gamma_air * P * (1 - specieFraction) + gamma_self * P * specieFraction)

        fL = 1 / np.pi * gamma / (gamma ** 2 + (nu[:, np.newaxis] - (nui + delta * P)) ** 2)

        #alphaD = nui / c * np.sqrt(2 * Na * kb * T * np.log(2) / molarMass)

        #fG = np.sqrt(np.log(2) / np.pi / alphaD ** 2) * np.exp(- (nu[:, np.newaxis] - nui) ** 2 * np.log(2) / (alphaD ** 2))


        t = np.sum(Sw * fL, axis=1) * (101325 * P / kb / T) * (1e-2) ** 3 * specieFraction
        #t = np.sum(Sw * f, axis=1)

        return(t)


    def plotData(self, y='linear'):
        """
        Creates a graph object using matplotlib.pyplot.plot function (does not show it. you still need to use matplotlib.pyplot.show() afterwards)
        :param y (str): Argument directly passed to matplotlib.pyplot.yscale. See matplotlib documentation for more insight. Two common values are 'linear' and 'log'
        :return: None
        """
        if not self.isIsotpologueData:
            raise ValueError('Impossible to plot data because it contains multiple isotopologues. It will be corrected in next versions of the code')
        nu = []
        Sw = []
        for line in self:
            nu.append(line.nu)
            Sw.append(line.sw)
        plt.plot(nu, Sw)
        plt.yscale(y)
        plt.title(f'Absorption coefficient of {self.specie}')
        plt.xlabel('Frequency [cm-1]')
        plt.ylabel('Absorption coefficient, isotopologue abundance acounted for [cm-1 / (molecule * cm-2)]')



class DataHandler:
    """
    Handles all operations related to stored data: downloading and reading
    """

    API_FORMAT = 'https://hitran.org/lbl/api?iso_ids_list={}&numin={}&numax={}'

    def download(self, specie, nu_min, nu_max):
        """
        Function used to download data from HITRAN via API. data stored in data folder and using the naming pattern 'bandes_{specie}.par'
        :param specie (str): Name of the specie to download (all known isotopologs will be downloaded)
            The specie name must be in convertion table ISO
        :param nu_min (float): lower limit of frequency to download
        :param nu_max (float): upper limit of frequency to download
        :return: None
        """
        ids = []
        for val in corr.ISO.values():
            if val[-1] == specie:
                ids.append(val[0])
        if len(ids) == 0:
            raise ValueError(f'No data found for specie {specie}')
        else:
            iso_ids_list = ''
            for id in ids:
                iso_ids_list += str(id) + ','
            iso_ids_list = iso_ids_list[:-1]  # remove last comma
            data = requests.get(DataHandler.API_FORMAT.format(iso_ids_list, nu_min, nu_max)).text
            with open(f'data/bandes_{specie}.par', 'w') as file:
                file.write(data)

    def readData(self, specie):
        """
        Function used to read data sotored in the data folder
        :param specie (str): Name of the specie to read (all known isotopologues will be read)
        :return (HITRANData): HITRANData object representing specie data
        """
        res = HITRANData(specie)

        with open(f'data/bandes_{specie}.par', 'r') as file:
            data = file.readlines()

        for line in data:
            res.append(HITRANStructure(
                molecule_id=int(line[0:2]),
                local_iso_id=int(line[2:3]),
                nu=float(line[3:15]),
                sw=float(line[15:25]),
                a=float(line[25:35]),
                gamma_air=float(line[35:40]),
                gamma_self=float(line[40:45]),
                elower=float(line[45:55]),
                n_air=float(line[55:59]),
                delta_air=float(line[59:67]),
                global_upper_quanta=line[67:82],
                global_lower_quanta=line[82:97],
                local_upper_quanta=line[97:112],
                local_lower_quanta=line[112:127],
                ierr=line[127:128],
                iref=line[128:130],
                line_mixing_flag=line[130:131],
                gp=line[131:138],
                gpp=line[138:145],
            ))
        return(res)


if __name__ == '__main__':
    a = DataHandler()
    CO2 = a.readData('CO2')
    CO2_iso_1 = CO2.getIsotopologueByID(1)
    for l in CO2_iso_1:
        print(l.local_iso_id)
    CO2_iso_1.plotData('log')
    plt.show()
    #a.download('H2O', 200, 1600)
    H2O = a.readData('H2O')
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