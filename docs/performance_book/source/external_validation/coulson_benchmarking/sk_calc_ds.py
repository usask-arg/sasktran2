import xarray as xr
import sasktran2 as sk
import numpy as np
import time
start = time.time()

Albedo = [0, 0.25, 0.8]
Mu0 = np.linspace(0.02, 1, 50)
Mu = Mu0
Phi = np.arange(0, 181, 30)
Tau = [1, 2, 4, 8, 16, 32, 64, 100, 128, 256, 512, 1024]
Stokes = ["I", "Q", "U"]

print(f"Time to run up to sasktran2 (min): {(time.time() - start)/60}")

# Initialize Storage of Stokes numpy arrays calculated by sasktran2
Temp_na_I = np.zeros(shape=(12, 3, 50, 50, 7))
Temp_na_Q = np.zeros(shape=(12, 3, 50, 50, 7))
Temp_na_U = np.zeros(shape=(12, 3, 50, 50, 7))

# Set up sasktran2
config = sk.Config()
config.multiple_scatter_source = sk.MultipleScatterSource.DiscreteOrdinates
config.single_scatter_source = sk.SingleScatterSource.DiscreteOrdinates
config.input_validation_mode = sk.InputValidationMode.Disabled
config.num_streams = 40
config.num_singlescatter_moments = 40
config.num_stokes = 3

count_tau = 0
for tau in Tau:
    count_a = 0
    for a in Albedo:
        print(f"New A! =>{a} for Tau={tau}\n", (time.time() - start)/60, '\n')  # Checkpoint times
        count_cos_sza = 0
        for cos_sza in Mu0:
            model_geometry = sk.Geometry1D(cos_sza=cos_sza,
                                           solar_azimuth=0,
                                           earth_radius_m=6372000,
                                           altitude_grid_m=np.array([0, 1.0]),
                                           interpolation_method=sk.InterpolationMethod.LinearInterpolation,
                                           geometry_type=sk.GeometryType.PlaneParallel)

            # Viewing Geometry Config
            viewing_geo = sk.ViewingGeometry()

            # Look at all values in the tables
            for mu in Mu:
                for phi in Phi:
                    viewing_geo.add_ray(sk.GroundViewingSolar(cos_sza, np.deg2rad(phi), mu, 2.0))

            # Setting up Atmosphere
            atmosphere = sk.Atmosphere(model_geometry, config, calculate_derivatives=False, numwavel=1)

            atmosphere.storage.total_extinction[:] = tau  # Tau Value
            atmosphere.storage.ssa[:] = 1

            # Legendre coefficient expansion for Rayleigh scattering
            atmosphere.leg_coeff.a1[0, :, 0] = 1  # Const term
            atmosphere.leg_coeff.a1[2, :, 0] = 0.5  # cos^2 term

            # Polarization of light
            atmosphere.leg_coeff.a2[2] = 3
            atmosphere.leg_coeff.b1[2] = -np.sqrt(6.0) / 2

            atmosphere.surface.albedo[:] = a  # Surface Reflectance

            # Run Calculation
            engine = sk.Engine(config, model_geometry, viewing_geo)
            radiance = engine.calculate_radiance(atmosphere)

            radiance["radiance"] *= np.pi  # these all upwelling at TOA

            # Sort Radiance into numpy arrays of the 3 Stokes Parameters
            temp_lis_I = [float(radiance["radiance"][0][i][0].values) for i in range(len(radiance["radiance"][0]))]
            # for i in range(len(temp_lis_I)):  # Rounding to zero
            #     if 1e-17 > i > -1e-17:
            #         temp_lis_I[i] = 0
            Temp_na_I[count_tau][count_a][count_cos_sza] = np.array(temp_lis_I).reshape(50, 7)

            temp_lis_Q = [float(radiance["radiance"][0][i][1].values) for i in range(len(radiance["radiance"][0]))]
            # for i in range(len(temp_lis_Q)):  # Rounding to zero
            #     if 1e-17 > i > -1e-17:
            #         temp_lis_Q[i] = 0
            Temp_na_Q[count_tau][count_a][count_cos_sza] = np.array(temp_lis_Q).reshape(50, 7)

            temp_lis_U = [float(radiance["radiance"][0][i][2].values) for i in range(len(radiance["radiance"][0]))]
            # for i in range(len(temp_lis_U)):  # Rounding to zero
            #     if 1e-17 > i > -1e-17:
            #         temp_lis_U[i] = 0
            Temp_na_U[count_tau][count_a][count_cos_sza] = np.array(temp_lis_U).reshape(50, 7)

            count_cos_sza += 1
        count_a += 1
    count_tau += 1


print(f"Time to run sasktran2 software (min): {(time.time() - start)/60} \n\n")

# temp list np arrays to data arrays (da) for each Stokes
Calc_da_I = xr.DataArray(Temp_na_I,
                         dims=["Tau", "A", "Mu0", "Mu", "Phi"],
                         coords={"Tau": Tau, "A": Albedo, "Mu0": Mu0, "Mu": Mu, "Phi": Phi}
                         )
Calc_da_Q = xr.DataArray(Temp_na_Q,
                         dims=["Tau", "A", "Mu0", "Mu", "Phi"],
                         coords={"Tau": Tau, "A": Albedo, "Mu0": Mu0, "Mu": Mu, "Phi": Phi}
                         )
Calc_da_U = xr.DataArray(Temp_na_U,
                         dims=["Tau", "A", "Mu0", "Mu", "Phi"],
                         coords={"Tau": Tau, "A": Albedo, "Mu0": Mu0, "Mu": Mu, "Phi": Phi}
                         )

# Create Dataset for da's and save to file
Calc_ds = xr.Dataset({"I": Calc_da_I, "Q": Calc_da_Q, "U": Calc_da_U})
Calc_ds.to_netcdf("Sk2_Stokes_Calc_Data_Set.nc")

print(Calc_ds, '\n\n')
print(Calc_ds.I, '\n\n')
print(Calc_ds.Q, '\n\n')
print(Calc_ds.U, '\n\n')
print("End Time: ", (time.time() - start)/60)
