# CHANGELOG

This file explain main changes between different GRASP versions.
To know how to update your settings files from a version to run in different one,
please, check the file [UPGRADE.md](UPGRADE.md)


## v0.7.4

- Adding a setting to dump initial guess information: input.imagedat.dump
- Important bug was fixed that mainly affect to retrieval with pixels with different number of wavelengths.
- Introduced module report stop: It allows to don't stop entire processing if single segment is wrong.
- Added two extra wildcards to streams: segment_first_date, segment_last_date and iinversion
- Added a new default transformer: segment_imagedat to set up different initial guess for each pixel in a single segment retrieval
- Minor improvements, bugs fixed, code more cleaned, documentation improved and code stability improved

## v0.7.3

- Models approach is integrated and can be loaded
- Refactor of SSA output
- Old bug fixed: Sometimes appear a pixels with all values set to 0 in the output
- Minor improvements, bugs fixed, code more cleaned, documentation improved and code stability improved

## v0.7.2

- Settings module allow to set default values for arrays with different values. This change affects to the settings interface of extension modules.
- Moved loops in inversion routine to their own module for better abstraction and easier parallelization.
- Minor improvements, bugs fixed, code more cleaned, documentation improved and code stability improved

## v0.7.1

- Refactoring for vertical profile part; discret vertical profile for given number of altitudes is passes as input to radiative transfer routine. 
- Possibility to run 2 aerosol mode retrieval with precalculated lognormal bin size distribution model
- Particle component volume fractions linear mixture model of refractive index contains: 
  carbon, dust, iron and water components for 1 aerosol mode retrieval;  
  black carbon, brown carbon, quartz and water (fine mode) and quartz, iron and water (coarse mode) for 2 aerosol mode retrieval.
- Other minor improvements, bugs fixed, code more cleaned, documentation improved and code stability improved

## v0.7.0

- Added license agreement
- Added official examples
- Minor improvements and bugs fixed


## v0.6.8

- Code has been cleaned removing many external code (extensions, modules, kernels...)
- grasp-manager can manage modules, kernels and constants sets.
- Minor improvements, bugs fixed, code more cleaned, improved documentation and improved code stability


## v0.6.7

- An important bug introduced in version v0.6.6 has been fixed. It affects to retrieval over ocean surfaces.
- Minor improvements, bugs fixed, code more cleaned, improved documentation and improved code stability


## v0.6.6

- Minor improvements, bugs fixed, code more cleaned, improved documentation and improved code stability


## v0.6.5

- Minor improvements and bugs fixed


## v0.6.4

- Added new wildcards to grasp output streams
- Many new actions added to grasp-manager script
- Code has been speeded up in not-mixed surface (pure land or pure ocean)
- Minor improvements, bugs fixed, code more cleaned, improved documentation and improved code stability

## v0.6.3

- Parameter **retrieval.constraints.characteristic[1].mode[1].single_pixel.a_priori_estimates.difference_order** has been removed and
  **retrieval.constraints.characteristic[].mode[].single_pixel.a_priori_estimates.lagrange_multiplier** accept an array instead of single value.
  It allows to define a different constraint to each parameter in the initial guess.  
- New option was added allowing to use in current segment retrieval retrieved parameters of pixels-neighbors from preceding retrievals.
  In order to use this option additional settings for edges have to be provided in settings file providing the number 
  of pixels-edges for longitude (X), latitude (Y) and time (T) coordinates to be used in current segment inversion (by default 0 is set). 
- Particulate matter (air quality) product added
- Aerosol typing added (aerosol type indicies: 
    0 – Complex mixture, 
    1 – Background Aerosol
    2 – Water/Maritime
    3 — Urban Polluted
    4 – Mixed aerosol
    5 – Urban Clean
    6 – Smoke Smoldering
    7 – Smoke flaming
    8 – Mineral dust)
- Option to select aerosol vertical profile type between normal and exponential added

```yml
retrieval:
    edges_size:
        x: 2
        y: 2
        t: 2
```

- Minor improvements, bugs fixed, code more cleaned, improved documentation and improved code stability

## v0.6.2

- Added a parameter in settings files: **retrieval.product_configuration.wavelength_indices_for_ndvi**: Indices of wavelengths which will be used to calculate NDVI if it is calculated
- Added "valgrind" constants set
- Minor improvements, bugs fixed, code more cleaned, improved documentation and improved code stability


## v0.6.1

- Added more actions to grasp-manager script
- Minor improvements, bugs fixed, code more cleaned, and improved code stability

## v0.6.0

