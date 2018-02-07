/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/*      MODULE mod_par_DLS
! ********************************************************** c
! **  Parameters for matrix_...f                          ** c
! ********************************************************** c
! **  KN1par - number of grid radii in original kernels   ** c
! **           and fixed kernels                          ** c
! **  KM1par - number of scattering angles in org.kernels ** c
! **  KR1par - number of grid aspect ratios for original  ** c
! **                                            kernels   ** c
! **  KREpar - number of refractive index real parts      ** c
! **  KIMpar - number of refractive index imaginary parts ** c
! **                                                      ** c
! **  KNpar  - number of grid radii for optical charact.  ** c
! **  KRpar  - number of grid aspect ratios for axial     ** c
! **           ratio distribution                         ** c
! **  KMpar  - number of scattering angles in fixed       ** c
! **           kernels and for opt.characteristics        ** c
! **  KMD    - number of modes in size distribution       ** c
! **           (up to 2)                                  ** c
! **  rootdir - kernel directory name                     ** c
! ********************************************************** c
*/                      

/* parameters for Original or Fixed Kernels */

#define  _KN1par 41
#define  _KR1par 2  
#define  _KM1par 181 
#define  _KREpar 22
#define  _KIMpar 16 


/* parameters for Optical Characteristics */

#define _KNpar 22
#define _KRpar 2
#define _KMpar 181 


/* parameters for Size Distribution */

#define _KMD 2 
#define _rootdir "./"

