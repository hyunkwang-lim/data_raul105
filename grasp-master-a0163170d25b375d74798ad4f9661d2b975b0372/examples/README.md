EXAMPLES
========

This directory includes subdirectories with examples for several GRASP applications. Subdirectories contain examples of GRASP settings and measurement data (sdata) input files together with sample of the code output corresponding to these files. 

- **sunphotometer:**
    * **settings\_example\_sunphotometer\_inversion.yml** is an example of settings that use sdata driver and is set to invert sunphotomter total optical depth and radinaces measurements and to retrieve size distribution, spectral dependent real and imaginary parts of refractive index and sphericity fraction for sigle mode of aerosol.  
    * **example\_sunphotometer.sdat** is an example of synthetic sunphotomter total optical depth and radinaces measurements.
    * **example\_sunphotometer\_inversion\_output.txt** is an example of GRASP output provided after performing the inversion using two files described above.
- **lidar\_and\_sunphotometer:** 
    * **settings\_example\_lidar\_and\_sunphotometer.yml** is an example of settings that use sdata driver and is set to invert simultaneously sunphotometer and lidar measurements and to retrieve size distribution, spectral dependent real and imaginary parts of refractive index, aerosol vertical distribution profile distinguished between fine and coarse aerosol modes and same sphericity fraction for both modes. 
    * **example\_lidar\_and\_sunphotometer.sdat** is an example of combined synthetic sunphotomter and lidar measurements.
    * **example\_sunphotometer\_and\_lidar\_inversion\_output.txt** is an example of GRASP output provided after performing the inversion using two files described above.
- **polder:** 
    * **settings\_example\_polder\_inversion.yml** is an example of settings that use sdata driver and is set to invert simultaneously multiple pixel polder measurements reading them from data input file and to retrieve size distribution (**16 triangular bins**), spectral dependent real and imaginary parts of refractive index, sphericity fraction and aerosol vertical distribution scale altitude for sigle mode of aerosol together with the Ross-Li and Maignan-Breon properties for the underlying land surface and Cox-Munk properties for underlying ocean surface for each of the measurement pixels. 
    * **settings\_example\_polder\_inversion_LB.yml** is an example of settings that use sdata driver and is set to invert simultaneously multiple pixel polder measurements reading them from data input file and to retrieve size distribution (**5 precalculated lognormal bins**), spectral dependent real and imaginary parts of refractive index, sphericity fraction and aerosol vertical distribution scale altitude for sigle mode of aerosol together with the Ross-Li and Maignan-Breon properties for the underlying land surface and Cox-Munk properties for underlying ocean surface for each of the measurement pixels.
    * **example\_polder.sdat** is an example of multi-pixel (2x2) multi-temporal (30 days) synthetic POLDER measurements of stockes vector (I,Q and U) over a region with both land, ocean and mixed regions.
    * **example\_polder\_inversion\_output.yml** is an example of GRASP output provided after performing the inversion using example\_polder\_inversion.yml and example\_polder.sdat.
    * **example\_polder\_inversion\_output\_LB.yml** is an example of GRASP output provided after performing the inversion using example\_polder\_inversion_LB.yml and example\_polder.sdat.





