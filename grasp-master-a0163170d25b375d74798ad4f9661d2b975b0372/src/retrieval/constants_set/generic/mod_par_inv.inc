/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */


/* RELATED TO MEASUREMENTS SPECIFICATIONS */

/* KWM: max number of wave lengths for measurement */
#define _KWM 12 

/* KW: number of wavelengths */
#define _KW 12 

/* KSHAPE: number of aerosol components  with different shape distributions */
#define _KSHAPE 2 

/* KBF: number of parameters in BRF model */
#define _KBF 4 

/* RELATED TO INVERSION ORGANIZATION */

/* KIX: number of pixels corresponding with different X - coordinates */
#define _KIX 2

/* KIY: number of pixels corresponding with different Y - coordinates */
#define _KIY 2

/* KITIME: number of pixels corresponding to different times */
/* ! 92x2x2; 25x5x5 - Parasol data */
#define _KITIME 30

/* KIMAGE: number of pixels */
#define _KIMAGE (_KITIME*_KIX*_KIY)

/* KIEDGE: max number of edge pixels for different X,Y,T - coordinates of inverted segment edges */
#define _KIEDGE 2

/* KDIM1: number of retrieved characteristics (SD,REAL,IMAG,SHAPE, height ) */
#define _KIDIM1 10 

/* KDIM2: number of sub-component in KIDIM1  */
#define _KIDIM2 4

/* KDIM3: number of sub/sub/sub-component in KIDIM1 */
/* KIDIM3=22  KIDIM3=61 */
#define _KIDIM3 61

/* KPARS: number of retrieved parameters for each pixel */
#define _KPARS 209

/* KKNOISE: number of measured optical characteristics with different noise distribution  */
#define _KKNOISE 6 

/* KVERTM: maximum of number of measured profile vert. heights and retrieved profile vert. heights + 2 */
#define _KVERTM 63

/* KNBVM: parameter for size of measurement arrays (max of number of angles or heights) */
/* max(_KVERTM,_NBVM) */
#define _KNBVM 300

/* KIP: max number of meas types for single pixel */
#define _KIP 4 

/*??? ask Oleg _KMESS (_KNBVM*_KWM*_KIP) */
#define _KMESS (_KNBVM*_KW*_KIP) 

/* */
#define _KPAR  (_KPARS*_KIMAGE)

/* nnz_par: parameter for sparse matrix */
#define _nnz_par (_KPAR*(_KPARS+_KIMAGE-1))

/* KMPSM: parameter for multi pixel smoothness */
#define _KMPSM 3*(_KITIME*_KIX+_KITIME*_KIY+_KIX*_KIY)

/* KDIF: parameter for differences in smoothness KDIF=MAX of (KPARS,KTIME,KIX,KIY,_KIEDGE) */
#define _KDIF  (max(_KPARS,_KITIME,_KIX,_KIY,_KIEDGE))

/* */
#define _KSURF 1  

/* #define _KGAS  1 */
 

