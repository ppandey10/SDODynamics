import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("~/custom.mplstyle")

# Define the path to the text file
file_path = 'mpc-data-scattering-centaurs.txt'

# Read the text file, skipping the first line and using fixed-width format
column_names = [
    "Designation", "Prov", "Des", "q", "Q", "H", "Epoch", "M", "Peri", 
    "Node", "Incl", "e", "a", "Opps", "Ref", "Designation_name", "Discovery_date_site_discoverer"
]

# Define the fixed widths for each column based on the structure of the file
column_widths = [
    27, 5, 7, 7, 10, 7, 10, 7, 6, 6, 6, 6, 9, 9, 19, 33, 48
]

# Read the file using pandas.read_fwf
df = pd.read_fwf(file_path, widths=column_widths, names=column_names, skiprows=1)

fig, ax = plt.subplots()
ax.scatter(df["a"], df["e"], s=3)
ax.set_xscale("log")
plt.show()

# TODO: remove centaurs, calculate perihelion for 3:2, 2:1 resonances
