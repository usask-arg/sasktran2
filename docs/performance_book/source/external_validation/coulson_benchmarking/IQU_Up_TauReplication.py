import xarray as xr
import numpy as np

# Takes about 45 minutes to run

# Open Benchmark data file
# if working in python console, need "SaskTran2 Testing/2009JournalDataReplication/Sasktran2 Addition/Stokes_Benchmark_Data_Set.nc"
Benchmark_ds = xr.open_dataset("Stokes_Benchmark_Data_Set.nc")
Calc_ds = xr.open_dataset("Sk2_Stokes_Calc_Data_Set.nc")

# Create Percent Difference Dataset (doesn't account for division by zero)
P_Diff_ds = xr.Dataset({"I": ((Benchmark_ds.I - Calc_ds.I)/Benchmark_ds.I) * 100,
                        "Q": ((Benchmark_ds.Q - Calc_ds.Q)/Benchmark_ds.Q) * 100,
                        "U": ((Benchmark_ds.U - Calc_ds.U)/Benchmark_ds.U) * 100})

# Create ds to replace division by zero values
Mask_Diff_ds = xr.Dataset({"I": Calc_ds.I,
                           "Q": (Benchmark_ds.Q - Calc_ds.Q)/Calc_ds.I,
                           "U": (Benchmark_ds.U - Calc_ds.U)/Calc_ds.I})
# print(Mask_Diff_ds.I, '\n\n')
# print(Mask_Diff_ds.Q, '\n\n')
# print(Mask_Diff_ds.U, '\n\n')
# Go through all P_Diff to clean up nan, inf, -inf
P_Diff_ds = P_Diff_ds.where(P_Diff_ds.map(np.isfinite)).fillna(Mask_Diff_ds)

# print("PI", P_Diff_ds.I, '\n\n')
# print("PQ", P_Diff_ds.Q, '\n\n')
# print("PU", P_Diff_ds.U, '\n\n')

