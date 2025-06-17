import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

# Parameters
mu0 = 4 * np.pi * 1e-7  # Permeability of free space
I = 57.0                # Current in coil A
N = 210                 # Number of superposed loops
r = 0.06                # Radius of coils (m)
d = 0.011               # Gap between train and railway
L = 4.0                 # Width between guideway sides
l = 0.1                 # Distance between the two coils of the electromagnet

# Laplace force expression (simplified 2D model from appendix)
def Fy(rho, y):
    return -mu0 * I**2 * N**2 * r * (rho / ((rho**2 + y**2)**(3/2)))

# Total force function
def F_total(y):
    dy_left = 2*d - 2*y
    dy_right = 2*d + 2*y
    return 2 * (Fy(0, dy_left) - Fy(l, dy_left)) + 2 * (Fy(l, dy_right) - Fy(0, dy_right))

# Evaluate over range of y
y_vals = np.linspace(-0.01, 0.01, 1000)
F_vals = np.vectorize(F_total)(y_vals)
U_vals = np.concatenate(([0], cumtrapz(-F_vals, y_vals)))

# Plotting
plt.figure(figsize=(10, 4))

plt.subplot(2, 1, 1)
plt.plot(y_vals, F_vals, label='Total force')
plt.axhline(0, color='gray', linestyle='--')
plt.xlabel('Lateral displacement y (m)')
plt.ylabel('Force (N)')
plt.title('Force vs Displacement – Method of Images (Two Coils per Electromagnet)')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(y_vals, U_vals, color='orange')
plt.xlabel('Lateral displacement y (m)')
plt.ylabel('Potential Energy (J)')
plt.title('Potential Energy Profile')
plt.grid(True)

plt.tight_layout()
plt.show()
