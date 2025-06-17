import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

# Constants and parameters
mu0 = 4 * np.pi * 1e-7  # Vacuum permeability [H/m]
I = 57.0                # Current of coil A [A]
beta = 1.0              # Ratio IB / IA
IA = I
IB = beta * I
r = 0.002               # Radius of coil [m]
d = 0.011               # Air gap at equilibrium [m]
L = 4.0                 # Width of the guideway [m]

# Dipole interaction function
def f(y):
    return 3 * mu0 * np.pi * r**4 * y / (2 * (r**2 + y**2)**(5/2))

# Forces from each side
def FA(y):
    return IA**2 * f(2*d - 2*y) + IA*IB * f(L + 2*d) + IA*IB * f(L - 2*y)

def FB(y):
    return - (IB**2 * f(2*d + 2*y) + IA*IB * f(L + 2*d) + IA*IB * f(L + 2*y))

def Ft(y):
    return FA(y) + FB(y)

# Displacement range
y_vals = np.linspace(-d + d/10, d - d/10, 1000)
F_vals = Ft(y_vals)

# Numerical integration for potential energy
U_vals = np.concatenate(([0], cumtrapz(-F_vals, y_vals)))

# Plotting
plt.figure(figsize=(10, 4))

plt.subplot(2, 1, 1)
plt.plot(y_vals, F_vals, label='Total force')
plt.axhline(0, color='gray', ls='--', lw=0.5)
plt.xlabel('Lateral displacement y (m)')
plt.ylabel('Force (N)')
plt.title('Lateral Force vs Displacement (Method of Images – 1 Coil)')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(y_vals, U_vals, color='orange')
plt.xlabel('Lateral displacement y (m)')
plt.ylabel('Potential Energy U (J)')
plt.title('Potential Energy Profile')
plt.grid(True)

plt.tight_layout()
plt.show()
