import pandas as pd
import numpy as np
import csv
import re
import xarray as xr

A = [0, 0.25, 0.8]
mu0 = np.linspace(0.02, 1, 50)
mu = mu0
phi = np.arange(0, 181, 30)
Tau = [1, 2, 4, 8, 16, 32, 64, 100, 128, 256, 512, 1024]
Stokes = ["I", "Q", "U"]
for s in Stokes:
    for t in Tau:
        # Open and Read
        file = open(f"Sasktran2 Addition/Stokes.tar/{s}_UP_TAU_{t}")
        reader = csv.reader(file, delimiter=f" ")
        reader = [[[re.sub(r'[^\x00-\x7F]+', '-', k) for k in num] for num in line if num != ''] for line in reader if
                  len(line) >= 11 if line[3] != '0.00']
        Benchmark_lis = [[i for i in line] for line in reader]
        Benchmark_lis = [["".join(i) for i in lis] for lis in Benchmark_lis]
        file.close()

        # Sort data wanted into list and then into df
        templis = []
        for lis in Benchmark_lis:
            if lis[0] != 'TAU' and lis[0] != 'mu0':
                lis = [float(st) for st in lis]
                templis.append(lis[2::])

        cols = [f"phi={j}" for j in np.arange(0, 181, 30)]
        df1 = pd.DataFrame(templis, columns=cols)

        if t == 1:
            df = df1
        else:
            df = df._append(df1, ignore_index=True)

        # turn df to numpy array then to data array
        if t == 1024 and s == "I":
            na = df.to_numpy()
            na = na.reshape(12, 3, 50, 50, 7)
            naI = na
            daI = xr.DataArray(naI,
                               dims=["Tau", "A", "Mu0", "Mu", "Phi"],
                               coords={"Tau": Tau, "A": A, "Mu0": mu0, "Mu": mu, "Phi": phi}
                               )

        if t == 1024 and s == "Q":
            na = df.to_numpy()
            na = na.reshape(12, 3, 50, 50, 7)
            naQ = na * -1  # fix signs
            daQ = xr.DataArray(naQ,
                               dims=["Tau", "A", "Mu0", "Mu", "Phi"],
                               coords={"Tau": Tau, "A": A, "Mu0": mu0, "Mu": mu, "Phi": phi}
                               )

        if t == 1024 and s == "U":
            na = df.to_numpy()
            na = na.reshape(12, 3, 50, 50, 7)
            naU = na * -1  # fix signs
            daU = xr.DataArray(naU,
                               dims=["Tau", "A", "Mu0", "Mu", "Phi"],
                               coords={"Tau": Tau, "A": A, "Mu0": mu0, "Mu": mu, "Phi": phi}
                               )


# Add all data arrays to a data set
Benchmark_ds = xr.Dataset({"I": daI, "Q": daQ, "U": daU})
# print(Benchmark_ds, '\n\n')

Benchmark_ds.to_netcdf("Stokes_Benchmark_Data_Set.nc")
