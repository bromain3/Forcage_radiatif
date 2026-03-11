import math


class AtmosphericLayer:
    """
    Data structure defining all variable describing an atmospheric layer.
    Can be replicated by a struct type in C
    """

    def __init__(self, h_min, h_max, T_gradient, T_layer_base, P_layer_base):
        """
        Initialization of the struct. No need for that in C
        :param h_min (float): Minimal geopotential height of the layer in meter
        :param h_max (float): Maximal geopotential height of the layer in meter
        :param T_gradient (float): Temperature gradient of the layer in Kelvin/kilometer
        :param T_layer_base (float): Temperature of the layer base in Kelvin
        :param P_layer_base (float): Pressure of the layer base in Pa
        """
        self.h_min = h_min
        self.h_max = h_max
        self.T_gradient = T_gradient
        self.T_layer_base = T_layer_base
        self.P_layer_base = P_layer_base


class StandardAtm:
    """
    Class containing all functions related to computing Standard atmosphere properties as defined in US Standard Atmosphere, 1976.
    In C, as you don't have classes you have to dissociate every functions by passing variables contained in self as inputs
    """
    T_STD = 288.15  # mean see level temperature (K)
    P_STD = 101325  # mean see level pressure (Pa)
    RHO_STD = 1.225  # mean see level density (kg/m3)
    g = 9.80665  # gravitational constant (m/s2)
    r = 287.05  # specific gaz constant for air (J/kg.K)
    Re = 6_356_766  # Mean Earth radius (m)
    LAYERS = [
        AtmosphericLayer(0, 11000, -0.0065, 288.15, 101325),  # Troposphère
        AtmosphericLayer(11000, 20000, 0.0000, 216.65, 22632),  # Stratosphère 1
        AtmosphericLayer(20000, 32000, 0.0010, 216.65, 5474.9),  # Stratosphère 2
        AtmosphericLayer(32000, 47000, 0.0028, 228.65, 868.02),  # Stratosphère 3
        AtmosphericLayer(47000, 51000, 0.0000, 270.65, 110.91),  # Stratosphère 4
        AtmosphericLayer(51000, 71000, -0.0028, 270.65, 66.94),  # Mésosphère 1
        AtmosphericLayer(71000, 86000, -0.0020, 214.65, 3.96)  # Mésosphère 2
    ]

    def h2z(self, h):
        """
        Converts geopotential height to geometric altitude
        :param h (float): Geopotential height in meter
        :return (float): Geometric altitude in meter
        """
        return(StandardAtm.Re * h / (StandardAtm.Re - h))

    def z2h(self, z):
        """
        Converts geometric altitude to geopotential height
        :param z (float): Geometric altitude in meter
        :return (float): Geopotential altitude in meter
        """
        return(StandardAtm.Re * z / (StandardAtm.Re + z))

    def CCO2(self, h):
        """
        Model for CO2 mol fraction. Based on constant concentration profile
        :param h: Geopotential altitude in meter
        :return: CO2 mol fraction
        """
        """
        C0 = 420e-6  # ground concentration of CO2
        H = 8.5e3
        return(C0 * math.exp(- h / H))
        """
        return(420e-6)


    def CCH4(self, h):
        """
        Model for CH4 mol fraction. Based on constant concentration profile
        :param h: Geopotential altitude in meter
        :return: CH4 mol fraction
        """
        """
        C0 = 1.9e-6  # ground concentration of CO2
        H = 25e3
        return(C0 * math.exp(- h / H))
        """
        return(1.9e-6)

    def CN2O(self, h):
        """
        Model for N2O mol fraction. Based on constant concentration profile
        :param h: Geopotential altitude in meter
        :return: N2O mol fraction
        """
        """
        C0 = 0.33e-6  # ground concentration of CO2
        H = 20e3
        return(C0 * math.exp(- h / H))
        """
        return(0.33e-6)

    def CH2O(self, h):
        """
        Model for N2O mol fraction. Based on constant concentration profile
        :param h: Geopotential altitude in meter
        :return: N2O mol fraction
        """
        """
        C0 = 0.01  # ground concentration of CO2
        H = 2.5e3
        return(C0 * math.exp(- h / H))
        """
        return(0.01)


    def getAtmProperties(self, z):
        """
        Computes atmosphere properties given geometric altitude
        :param z (float): Geometric altitude in meter
        :return tuple[float, float, float]: Temperature, pressure and density of the atmosphere at geometric altitude computed with US Standard Atmosphere model
        """
        h = self.z2h(z)

        if h < 0:
            raise ValueError(f'Impossible to compute standard atmosphere values for negative geopotential altitudes (h = {h})')
        elif h > 86_000:
            raise ValueError(f'Impossible to compute standard atmosphere values above a geopotential altitude of 86000 m (h = {h}) with the US standard atmosphere model')


        for n, layer in enumerate(StandardAtm.LAYERS):
            if layer.h_min <= h < layer.h_max:
                break
        # This way not be clean but the variable layer is still accessible outside the for loop

        dh = h - layer.h_min

        T = layer.T_layer_base + layer.T_gradient * dh

        if layer.T_gradient == 0:  # If the layer is isothermal
            P = layer.P_layer_base * math.exp(- StandardAtm.g * dh / (StandardAtm.r * layer.T_layer_base))
        else:
            P = layer.P_layer_base * (T / layer.T_layer_base) ** (- StandardAtm.g / (layer.T_gradient * StandardAtm.r))

        rho = P / (StandardAtm.r * T)

        speciesConcentration = {'CO2': self.CCO2(h), 'CH4': self.CCH4(h), 'N2O': self.CN2O(h), 'H2O': self.CH2O(h)}
        return(T, P, rho, speciesConcentration)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    a = StandardAtm()
    z = np.linspace(0, 86000, 1000)
    res = {'T': [], 'P': [], 'rho': []}
    for alt in z:
        tmp = a.getAtmProperties(alt)
        print(f"Altitude: {alt:6.0f} m | Température: {tmp[0]:6.2f} K | Pression: {tmp[1]:10.2f} Pa | Densité: {tmp[2]:.5f} kg/m³")
        res['T'].append(tmp[0])
        res['P'].append(tmp[1])
        res['rho'].append(tmp[2])

    z = z / 1000

    plt.subplot(1,3, 1)
    plt.plot(res['T'], z)
    plt.title('T')
    plt.subplot(1, 3, 2)
    plt.plot(res['P'], z)
    plt.title('P')
    plt.subplot(1, 3, 3)
    plt.plot(res['rho'], z)
    plt.title('rho')
    plt.show()