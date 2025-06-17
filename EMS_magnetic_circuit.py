import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

# Physical constants / parameters (example values)
mu0 = 4 * np.pi * 1e-7  # Vacuum permeability [H/m]
I = 57.0                # Current of coil A [A]
alph = 1.0              # Current of coil B [A] is alph * I

R = 0.06                # Radius [m]
d = 0.011               # Distance parameter [m]
N = 210                 # Number of coils
A = np.pi * R**2        # Cross-sectional area [m^2]

def F(dist):
    return (mu0 * A * N**2) / 4 * (I / dist)**2

# Composite forces
def FA(x):
    return F(d - x)

def FB(x):
    return -alph**2 * F(d + x)

def Ft(x):
    return FA(x) + FB(x)

# Set up x grid
x_vals = np.linspace(-d + d / 10, d - d / 10, 1000)

# Compute on that grid
F_vals = Ft(x_vals)
FA_vals = FA(x_vals)
FB_vals = FB(x_vals)

# Numerically integrate -F to get U, with U(-d) = 0
U_vals = np.concatenate(([0], cumtrapz(-F_vals, x_vals)))

# Plotting
plt.figure(figsize=(10, 4))

plt.subplot(2, 1, 1)
plt.plot(x_vals, F_vals, label='F(x)')
plt.plot(x_vals, FA_vals, '--', label='F_A(x)')
plt.plot(x_vals, FB_vals, '--', label='F_B(x)')
plt.axhline(0, color='gray', ls='--', lw=0.5)
plt.xlabel('Lateral displacement y (m)')
plt.ylabel('Force (N)')
plt.title('Forces vs. lateral displacement')
plt.legend()
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(x_vals, U_vals)
plt.xlabel('Lateral displacement y (m)')
plt.ylabel('Potential energy U (J)')
plt.title('Potential energy vs. lateral displacement')
plt.grid(True)

plt.tight_layout()
plt.show()
