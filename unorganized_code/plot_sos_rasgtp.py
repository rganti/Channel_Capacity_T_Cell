import matplotlib.pyplot as plt
import numpy as np

sos = np.loadtxt("Sos")
ras = np.loadtxt("RasGTP")
plt.plot(sos, ras, linestyle="None", marker='o', label="Equal Rates")
plt.xlabel("Sos")
plt.ylabel("Ras-GTP")
plt.autoscale()

plt.savefig("sos_ras", format="pdf")
