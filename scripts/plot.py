import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

df = pd.read_csv('../logs/log.csv', delimiter=';')

cinematic_energy = df["cinetic_energy"]
potential_energy = df["ljs_potential"]
temperature = df["temp"]
sum_forces = df["sum_forces"]
x = [ i for i in range( len(temperature))]

# Create a figure with subplots
plt.figure(figsize=(10, 6))

# Plot the first subplot
plt.subplot(2, 2, 1)
plt.plot(x, potential_energy, label="potential")
plt.plot(x, cinematic_energy, label="cinematic")
plt.plot(x, cinematic_energy+potential_energy, label="total")
plt.title("energy")
plt.legend()

# Plot the second subplot
plt.subplot(2, 2, 2)
plt.plot(x, temperature)
plt.title("temperature")

# Plot the third subplot
plt.subplot(2, 2, 3)
plt.plot(x, sum_forces)
plt.title("sum of forces")

# Show the plots
plt.tight_layout()
plt.show()
plt.savefig("../img/result.png")



