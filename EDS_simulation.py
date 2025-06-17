import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.integrate import cumtrapz

# Constants
mu0 = 4 * np.pi * 1e-7

# Parameters
R_train = 1.5     # radius of train coil (m)
R_guide = 1.0     # radius of guideway coil (m)
L = 1.2           # lateral spacing between guideway coils (m)
z0 = 0.0          # guideway z-position
z1 = 0.7          # train z-position
I_train = 7e5     # current in train coil (A)
RE = 10.0         # resistance of guideway coil (Ohm)
v = 160.0         # train speed (m/s)

# Time & space setup
T = 0.5
dt = 1e-3
t = np.arange(0, T, dt)
x_train = v * t
y = 0.0  # lateral offset
Mh = lambda sign: sign * I_train * np.pi * R_train**2  # Dipole moment

# Compute distance vector from train coil to guideway coil (top and bottom loops)
def r_vector(xj, xh, yj, zloop):
    return np.array([xj - xh, yj - y, zloop - z1])

# B field from dipole at point r
def B_field(M, r):
    r_mag = np.linalg.norm(r)
    if r_mag == 0:
        return np.zeros(3)
    term1 = 3 * np.dot(M, r) * r / r_mag**5
    term2 = M / r_mag**3
    return (mu0 / (4 * np.pi)) * (term1 - term2)

# Guideway coil positions (2 per side, 2 loops per coil)
xj_vals = [R_guide, 3*R_guide]
sides = [-L/2, L/2]
loops = [+R_guide, -R_guide]

# Train has 4 coils with alternating dipoles
def train_pos(h):
    return v * t + (2*h - 1)*R_train

# Compute total flux over time for one 8-shape coil
def compute_flux():
    flux_t = np.zeros_like(t)
    for h in range(1, 5):
        M_vec = np.array([0, (-1)**h * Mh(1), 0])  # Dipole in y-direction
        xh = train_pos(h)
        for side in sides:
            for xj in xj_vals:
                for zloop in loops:
                    for i, xt in enumerate(xh):
                        r = r_vector(xj, xt, side, z0 + zloop)
                        B = B_field(M_vec, r)
                        flux_t[i] += B[1] * np.pi * R_guide**2 * (1 if zloop > 0 else -1)
    return flux_t

# Compute induced current using Faraday + Ohm's law
flux = compute_flux()
emf = -np.gradient(flux, dt)
current = emf / RE

# FFT for frequency analysis
freqs = fftfreq(len(t), dt)
fft_vals = fft(current)

# Plot current
plt.figure(figsize=(10, 4))
plt.plot(t, current)
plt.title("Induced Current in Guideway Coil")
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")
plt.grid(True)
plt.tight_layout()
plt.show()

# Compute potential energy as time-averaged dipole interaction
def compute_energy(y_disp=0.0):
    U_t = np.zeros_like(t)
    for h in range(1, 5):
        Mp = np.array([0, (-1)**h * Mh(1), 0])
        xh = train_pos(h)
        for side in sides:
            for xj in xj_vals:
                for zloop in loops:
                    sign = 1 if zloop > 0 else -1
                    for i, xt in enumerate(xh):
                        r = r_vector(xj, xt, side, z0 + zloop)
                        B = B_field(np.array([0, sign * current[i] * np.pi * R_guide**2, 0]), r)
                        U_t[i] += -np.dot(Mp, B)
    return U_t

U_t = compu_
