from __future__ import annotations

import contextlib
from pathlib import Path

import numpy as np

import sasktran2 as sk
from sasktran2.constituent import VMRAltitudeAbsorber
from sasktran2.database.web import StandardDatabase
from sasktran2.optical.base import OpticalProperty


def _atm_file_path(folder_name: str, file_name: str) -> Path:
    """
    Determine the location of a specified '.atm' file.

    Parameters
    ----------
    folder_name : str
        Name of folder with the file, e.g. 'fascode', 'mipas_1998', or 'mipas_2001'. Do not include slashes
    file_name : str
        Name of the file, e.g. 'std.atm', 'tro.atm', ...

    Returns
    -------
    Path
        Absolute path to the file.
    """
    return StandardDatabase().path(
        (Path("climatology") / folder_name / file_name).as_posix()
    )


def _atm_reader(atm_file: str) -> dict:
    """
    Reads in a file with the extension '.atm' which contains information about atmospheric profiles and returns a
    dictionary with the atmospheric parameters as keys and the profiles as values. Refer to the input atm file for
    units of each profile but typically heights are in km, temperatures in K, pressures in mb, and species populations
    in ppmv.

    Parameters
    ----------
    atm_file : str

    Returns
    -------
    dict of [species: str, profile: np.array]
    """
    profiles = (
        {}
    )  # dictionary where keys are the species tags and values are the profiles
    with Path.open(atm_file) as f:
        num_levels_set = False
        cur_profile = ""
        for line in f:
            if line[0] == "!":  # line is a comment
                pass
            elif line[0] == "*":  # line is a profile header
                cur_profile = ""
                for letter in line[1::]:
                    if letter == " ":
                        break  # end of species name
                    cur_profile += letter
                if cur_profile == "END\n":
                    break  # end of file
                profiles[cur_profile.upper()] = np.array([])
            elif (
                not num_levels_set
            ):  # first uncommented line gives number of levels in each profile
                num_levels_set = True
            else:  # line contains profile values
                no_space = list(filter(None, line.split(" ")))
                no_space_or_comma = []
                for item in no_space:
                    no_space_or_comma.append(next(filter(None, item.split(","))))
                for item in no_space_or_comma:
                    with contextlib.suppress(ValueError):
                        if item != "\n":
                            profiles[cur_profile.upper()] = np.append(
                                profiles[cur_profile.upper()], float(item)
                            )
    return profiles


def constituent(
    species: str,
    optical_property: OpticalProperty,
    dataset: str = "fascode",
    climatology: str = "std",
) -> VMRAltitudeAbsorber:
    """
    Creates a VMRAltitudeAbsorber constituent for the given species and profile from several reference atmosphere models available at
    http://eodg.atm.ox.ac.uk/RFM/atm/

    Note: Species in the "Minor" category do not have a climatology associated VMR.  Instead they have single VMR profile

    The following table lists the species supported by each climatology.
    Each dataset contains several different climatologies. Each climatology contains unique profiles of temperature,
    pressure, and VMRs of the major species. The VMRs for the minor species are the same for all climatologies in each
    dataset.

    +------------+---------------------------------+---------------------------------+---------------------------------+
    | Dataset    | Climatologies                   | Major Species                   | Minor Species                   |
    +============+=================================+=================================+=================================+
    | fascode    | tro (Tropical),                 | H2O, CO2, O3, N2O, CO, CH4, O2  | NO, SO2, NO2, NH3, HNO3, OH,    |
    |            | mls (Mid-Latitude Summer),      |                                 | HF, HCl, HBr, HI, ClO, OCS,     |
    |            | mlw (Mid-Latitude Winter),      |                                 | H2CO, HOCl, N2, HCN, CH3Cl,     |
    |            | sas (Sub-Arctic Summer),        |                                 | H2O2, C2H2, C2H6, PH3, COF2,    |
    |            | saw (Sub-Arctic Winter),        |                                 | SF6, H2S, CFCl3, CF2Cl2, CClF3, |
    |            | std (US Standard Atmosphere),   |                                 | CF4, CHCl2F, CHClF2, C2Cl3F3,   |
    |            |                                 |                                 | C2Cl2F4, C2ClF5, CCl4, ClONO2,  |
    |            |                                 |                                 | N2O5, HNO4                      |
    +------------+---------------------------------+---------------------------------+---------------------------------+
    | mipas_1998 | day_imk (Mid-Latitude Day)      | N2, O2, O3P, CO2, O3, H2O, CH4, | HCN, H2O2, F12, F14, F22, COF2, |
    |            | ngt_imk (Mid-Latitude Night)    | N2O, HNO3, CO, NO2, N2O5, ClO,  | OCS, NH3, SO2, CFCl3, C2H2,     |
    |            | win_imk (Polar Winter)          | HOCl, ClONO2, NO                | C2H6, CCl4, SF6, HNO4, CH3Cl,   |
    |            | sum_imk (Polar Summer)          |                                 | CClF3, CHCl2F, C2Cl3F3, C2Cl2F4 |
    +------------+---------------------------------+---------------------------------+---------------------------------+
    | mipas_2001 | day (Mid-Latitude Day)          | N2, O2, CO2, O3, H2O, CH4, N2O, | CClF3, CHCl2F, C2Cl3F3,         |
    |            | ngt (Mid-Latitude Night)        | HNO3, CO, NO2, N2O5, ClO, HOCl, | C2Cl2F4, C2ClF5, CH3Cl, H2S     |
    |            | win (Polar Winter)              | ClONO2, NO, HNO4, HCN, NH3,     |                                 |
    |            | sum (Polar Summer)              | F11, F12, F14, F22, CCl4, COF2, |                                 |
    |            | equ (Equatorial Day)            | H2O2, C2H2, C2H6, OCS, SO2, SF6 |                                 |
    +------------+---------------------------------+---------------------------------+---------------------------------+

    Parameters
    ----------
    species : str
        The species to include, see the table below for supported species
    optical_property : OpticalProperty
        The optical property to use for the given species
    dataset : str, optional
        Dataset to use, see table below for options, by default "fascode"
    climatology : str, optional
        Climatology to use, see table below for options, by default "std"


    Returns
    -------
    sk.constituent.VMRAltitudeAbsorber
        The resulting constituent for the given species and profile
    """
    atm_file = climatology if climatology.endswith(".atm") else climatology + ".atm"
    file_path = _atm_file_path(dataset, atm_file)
    data = _atm_reader(file_path)

    if species.upper() in data:
        species_vmr = data[species.upper()] / 1.0e6  # Convert from ppm to vmr
        species_heights_m = data["HGT"] * 1000.0  # Convert from km to m
    else:
        # This is a minor species, need to load the minor species data
        if dataset == "mipas_2001":
            file_path = _atm_file_path(dataset, "extra.atm")
        elif dataset == "mipas_1998":
            file_path = _atm_file_path(dataset, "extra_imk.atm")
        elif dataset == "fascode":
            file_path = _atm_file_path(dataset, "minor.atm")
        else:
            msg = f"{dataset} is not a valid dataset"
            raise ValueError(msg)
        data_minor = _atm_reader(file_path)
        species_vmr = data_minor[species.upper()] / 1.0e6  # convert from ppm to vmr
        species_heights_m = data_minor["HGT"] * 1000.0  # convert from km to m

    return VMRAltitudeAbsorber(optical_property, species_heights_m, species_vmr)


def add_to_atmosphere(
    atmosphere: sk.Atmosphere,
    species: dict,
    dataset: str = "fascode",
    climatology: str = "std",
    set_pressure_temperature: bool = True,
):
    """
    Adds multiple species objects to an atmosphere object using the reference atmospheric profiles from http://eodg.atm.ox.ac.uk/RFM/atm/.

    See :py:func:`constituent` for more information on the species, dataset, and climatology options.

    Parameters
    ----------
    atmosphere : sk.Atmosphere
        The atmosphere to add the species objects to
    species : dict
        A dictionary where keys are Species (e.g. "O3") and values are the optical properties associated with them.  For example
        Species = {"O3": sk.optical.O3DBM(), "NO2": sk.optical.NO2Vandaele()}
    dataset : str, optional
        Reference dataset to use see :py:func:`constituent`, by default "fascode"
    climatology : str, optional
        Climatology to use, see :py:func:`constituent`, by default "std"
    set_pressure_temperature: bool, optional
        If true the atmospheric pressure/temperature will also be set by the climatology, by default True
    """
    for specie, optical_property in species.items():
        atmosphere[specie] = constituent(specie, optical_property, dataset, climatology)

    if set_pressure_temperature:
        atm_file = climatology if climatology.endswith(".atm") else climatology + ".atm"
        file_path = _atm_file_path(dataset, atm_file)
        data = _atm_reader(file_path)

        altitudes_m = data["HGT"] * 1000.0  # Convert from km to m
        temperature_k = data["TEM"]
        pressure_pa = data["PRE"] * 100.0

        atmosphere.temperature_k = np.interp(
            atmosphere.model_geometry.altitudes(), altitudes_m, temperature_k
        )
        # Interpolate pressure in log space
        atmosphere.pressure_pa = np.exp(
            np.interp(
                atmosphere.model_geometry.altitudes(), altitudes_m, np.log(pressure_pa)
            )
        )
