# UPGRADE

This file explain how to update your settings files from a version to run in different one.
If you want to know the changes in code (bug-fixed, new features, ...) please, check 
the file [CHANGELOG.md](CHANGELOG.md)


## v0.7.4

nothing to do to update to this version from v0.7.3

## v0.7.3

nothing to do to update to this version from v0.7.2

## v0.7.2

nothing to do to update to this version from v0.7.1 for users. Internally, the interface
of settings module has changed so new extensions have to have it into account.


## v0.7.1

nothing to do to update to this version from v0.7.0


## v0.7.0

nothing to do to update to this version from v0.6.8


## v0.6.8

- A characteristic type has been renamed: **size_distribution_bimodal** now is **size_distribution_lognormal**
- Kernels now are not available by default and the user has to install them via grasp-manager. The name of kernels has changed and it will affect to **kernels_folder** parameter
- The retrieval.kernel section has been reorganized

Before:
```yml
retrieval:
    kernel:
        phase_matrix_package:
        size_binning_method_for_triangle_bins: logarithm
            elements_of_phase_matrix: 4
            path_to_kernels: "AERONET_n22k16_83/"     
        radius:
            mode[1]:
                min: 0.05
                max: 15.0
```

After:
```yml
retrieval:
    phase_matrix:
        size_binning_method_for_triangle_bins: logarithm
        number_of_elements: 4
        kernels_folder: "aeronet_n22k16_83"     
        radius:
            mode[1]:
                min: 0.05
                max: 15.0
```

- Some input settings related with polder inversion has been moved to polder driver. If you don't use them you can remove them. If you work with polder driver you'll need to update following section:

Before:
```yml
input:
    gas_correction:
        enable: true
        path_no2: <valid path>
        path_o3: <valid path>
    land_ocean_filter: none
```

After:
```yml
input: 
    driver_settings:
        polder:
            gas_correction:
                enable: true
                path_no2: <valid path>
                path_o3: <valid path>
            land_ocean_filter: none
```


## v0.6.7

nothing to do to update to this version from v0.6.6


## v0.6.6

- **retrieval.debug.silent** is renamed by **retrieval.debug.verbose** and the value is reverse comparing to before value. True is print. Before it was false.
- **retrieval.debug.retrieval_iteration_information** is renamed by **retrieval.debug.iteration_information** 


## v0.6.5

nothing to do to update to this version from v0.6.4

## v0.6.4

### Some characteristic names has been renamed

Some of the possible values of retrieval.constraints.characteristic[].type has been change:
- **complex_refractive_index_chemistry** is renamed by **particle_component_fractions_chemical_mixture**
- **sphericity_fraction** is renamed by **sphere_fraction**
- **sphericity_distribution** is renamed by **aspect_ratio_distribution**

## v0.6.3

### single_pixel.a_priori_estimates have changed to an array.

The parameter **retrieval.constraints.characteristic[].mode[].single_pixel.a_priori_estimates.difference_order** has been removed
and **retrieval.constraints.characteristic[].mode[].single_pixel.a_priori_estimates.lagrange_multiplier** have change to an array of values.
Before, same value was applied to the characteristic but now it is possible to modify it for each parameter. By defualt these values are 0 so 
possibly you need just to remove these options:

```yml
a_priori_estimates:
    difference_order:   0
    lagrange_multiplier: 0.0
```

after each characteristic type. If you had not defined the value like 0 you have to add an array with as many elements as the characteristic
type with same value that you had before.

### Default value of retrieval.kernel.phase_matrix_package.size_binning_method_for_triangle_bins has changed

Before, the default value of size_binning_method_for_triangle_bins was absolute. If you have a 
settings file which you are using implicity this value, now you'll need to write it explicitly
because by default 'logarithm' will be used.

## v0.6.2

### New setting wavelength_indices_for_ndvi

A new setting has been added:

retrieval.product_configuration.wavelength_indices_for_ndvi: Indices of wavelengths which will be used to calculate NDVI if it is calculated

Before, this value was hardcoded to work with POLDER data. The settings files to
process polder data have been updated with:

```yml
retrieval:
    product_configuration:
        wavelength_indices_for_ndvi: [4, 5]
```


### The measurements 11 (TOD) and 12 (AOD) have been renamed

To clarify the code some constants have been renamed. TAU and EXT now are AOD and TOD. 
This change affects settings file in the description of noises. Specifically the parameter 
retrieval.noises.noise[].measurement_type[].type can be set to "aod" and "tod" instead 
of "tau" and "ext". So now retrieval.noises.noise[].measurement_type[].type can
be: tod, aod, p11, p12, p22, p33, p34, p44, ls, dp, rl, i, q or u .


## v0.6.1

### General block in settings has been refactored

The old general block has been reorganized completely. Here the list of changes that you have to do
in order to adapt a settings file from v0.6.0:

- **retrieval.general.path_to_internal_files** moved to **retrieval.general.path_to_internal_files**
- **retrieval.general.stop_before_performing_retrieval** moved to **retrieval.convergence.stop_before_performing_retrieval** 
- **retrieval.general.minimization_convention** moved to **retrieval.convergence.minimization_convention**
- **retrieval.general.maximum_iterations_of_Levenberg-Marquardt** moved to **retrieval.convergence.maximum_iterations_of_Levenberg-Marquardt** 
- **retrieval.general.threshold_for_stopping** moved to **retrieval.convergence.threshold_for_stopping**  
- **retrieval.general.threshold_for_stopping_Q_iterations** moved to **retrieval.convergence.threshold_for_stopping_Q_iterations** 
- **retrieval.general.scale_for_finite_difference** moved to **retrieval.convergence.scale_for_finite_difference** 
- **retrieval.general.normal_system_solver** moved to **retrieval.convergence.normal_system_solver** 
- **retrieval.general.shift_for_applying_logarithm_to_negative_values** moved to **retrieval.convergence.shift_for_applying_logarithm_to_negative_values**  
- **retrieval.general.maximum_iterations_for_stopping** moved to **retrieval.convergence.maximum_iterations_for_stopping** 
- **retrieval.general.regime_of_measurement_fitting** moved to **retrieval.regime_of_measurement_fitting.polarization** 
- **retrieval.general.regime_of_radiance_measurement_fitting** moved to **retrieval.regime_of_measurement_fitting.radiance** 
- **retrieval.general.regime_of_multipixel_constraints.inversion_regime** moved to **retrieval.regime_of_multipixel_constraints.inversion_regime** 
- **retrieval.general.regime_of_multipixel_constraints.time_space_groups** moved to **retrieval.regime_of_multipixel_constraints.time_space_groups**
- **retrieval.general.regime_of_multipixel_constraints.time-scale** moved to **retrieval.regime_of_multipixel_constraints.time-scale** 
- **retrieval.general.regime_of_multipixel_constraints.x-scale** moved to **retrieval.regime_of_multipixel_constraints.x-scale**
- **retrieval.general.regime_of_multipixel_constraints.y-scale** moved to **retrieval.regime_of_multipixel_constraints.y-scale**
- **retrieval.general.wavelength_indices_for_angstrom** moved to **retrieval.product_configuration.wavelength_indices_for_angstrom** 
- **retrieval.general.binning_method** moved to **retrieval.kernel.phase_matrix_package.size_binning_method_for_triangle_bins** 


Following example shows an update procedure:

Original code (v0.6.0)

```yml
retrieval:    
    general:
        minimization_convention: logarithm
        stop_before_performing_retrieval: false
        binning_method: logarithm
        maximum_iterations_of_Levenberg-Marquardt: 15
        threshold_for_stopping: 1.0e-5
        path_to_internal_files: "../../src/retrieval/internal_files/"
        regime_of_measurement_fitting: linear_polarization
        normal_system_solver: sparse_matrix_solver
        threshold_for_stopping_Q_iterations: 1e-12
        scale_for_finite_difference: 1.0e-5  
        wavelength_indices_for_angstrom: [4, 5]
        regime_of_multipixel_constraints: 
            inversion_regime: single_pixel 

```
Code after update:

```yml
retrieval:    
    general:
        path_to_internal_files: "../../src/retrieval/internal_files/"
        
    convergence:
        stop_before_performing_retrieval: false
        minimization_convention: logarithm
        maximum_iterations_of_Levenberg-Marquardt: 15
        threshold_for_stopping: 1.0e-5
        threshold_for_stopping_Q_iterations: 1e-12
        scale_for_finite_difference: 1.0e-5  
        normal_system_solver: sparse_matrix_solver
                
    regime_of_measurement_fitting:
        polarization: linear_polarization                
        
    product_configuration:
        wavelength_indices_for_angstrom: [4, 5]
    
    regime_of_multipixel_constraints: 
        inversion_regime: single_pixel 

    kernel:
        phase_matrix_package:
            size_binning_method_for_triangle_bins: logarithm
```


## v0.6.0
