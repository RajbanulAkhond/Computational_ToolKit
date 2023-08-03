''' In this code, we use NumPy for array manipulations, and SciPy's integrate.quad function to calculate the definite integral. 
    Make sure to replace 'dos_data.csv' with the correct file path to your CSV file containing the Density of States (DOS) data in two columns, 
    the first column representing energies and the second column representing DOS values. 
    The rest of the code should work as expected, plotting the Quantum Capacitance against the External Potential.
'''

import numpy as np
import pandas as pd
import scipy.constants as const
from scipy import integrate
import matplotlib.pyplot as plt

# Constants
e = const.elementary_charge
k_B = const.Boltzmann
T = 300  # temperature in Kelvin
mass_mos2 = 2.12488e-21
mass_nc = 3.18161e-21

# Read Density of States (DOS) from CSV file
dos_data = pd.read_csv('dos_data.csv')
nc_o_dos = dos_data.to_numpy()

# Define the thermal broadening function
def thermal_broadening(E, Vex, T):
    return np.power(1 / np.cosh((e * (E - Vex)) / (2 * k_B * T)), 2)

# Define the range of external potentials
Vext_range = np.arange(-0.5, 0.51, 0.01)

# Initialize an array to store the capacitance values
Cq = np.zeros(len(Vext_range))

# Sweep the external potential
for i, Vex in enumerate(Vext_range):
    # Calculate the product of DOS and thermal broadening function
    dosF_T = nc_o_dos[:, 1] * thermal_broadening(nc_o_dos[:, 0], Vex, T)
    # Define the integral function
    integrand = lambda E: np.interp(E, nc_o_dos[:, 0], dosF_T, left=0, right=0)
    # Calculate the definite integral using quad function from SciPy
    Cq[i], _ = integrate.quad(integrand, -np.inf, np.inf)

Cq *= e**2 * (4 * k_B * T)**-1

# Plot the results
plt.plot(Vext_range, Cq / mass_nc)
plt.xlabel('External Potential (V)')
plt.ylabel('Quantum Capacitance (F/g)')
plt.show()
