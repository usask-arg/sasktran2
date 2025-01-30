
(_changelog)=
# Changelog

## Unreleased

## 2025.01.2
- Adds manylinux aarch64 wheel builds
- More graceful error handling when homogenous solution fails
- Enable OMP support in Mac wheels

## 2025.01.1
- Fixes a bug in ecef_to_sasktran2_ray

## 2025.01.0
- Added capability for thermal emissions, see [Atmospheric and Surface Emissions](users_guide/emissions.md)
- Several accuracy and speed improvements
- Dropped support for Python 3.10

## 2024.11.0
- Added {py:class}`sasktran2.viewinggeo.LimbVertical`
- Added {py:class}`sasktran2.solar.SolarGeometryHandlerAstropy`

## 2024.10.2
- Internal re-release with wheel upload fixed

## 2024.10.1
- Add internal backend for Mie calculations
- Add capability to include refractive raytracing
- Add wheel builds for Python 3.13

## 2024.10.0
- Add python 3.13 conda builds
- Add linux arm64 builds

## 2024.07.0
- Add `config.do_backprop` option for discrete ordinates source
- Add `freeze` method to particle size distributions

## 2024.03.0
- Fix binary wheel uploads on PyPI

## 2024.02.1
- Fix the license classifier on PyPI

## 2024.02.0
- Templated number of streams has been enabled in DO, resulting in large speedups for 2-stream calculations
- Added support for Mie calculations (Using SASKTRAN Legacy)
- Added support for HITRAN cross sections (Using SASKTRAN Legacy)

## 2024.01.1
- Model has been relicensed under the MIT licenses

## 2024.01.0
- Add support for Python 3.12 wheels
- Remove support for Python versions below 3.10
- SASKTRAN2 is now available in conda-forge
- Added the {py:class}`sasktran2.LinearizedMie` object to perform Mie calculations
- Added support for `sk.GeometryType.PlaneParallel` and `sk.GeometryType.PseudoSpherical`

## 2023.12.0
- First official release
