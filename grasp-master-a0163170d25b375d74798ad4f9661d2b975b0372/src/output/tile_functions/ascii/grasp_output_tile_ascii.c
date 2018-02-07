/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include "grasp_output_tile_ascii.h"
#include "../../../global/grasp_retrieval_characteristic_type.h"
#include <grasp/utils.h>
#include <inttypes.h>
#include <grasp/utils.h>

#include "yamlsettings/yamlsettings.h"
#include "yamlsettings/yamlsettings_dictionary.h"
#include "yamlsettings/yamlsettings_assign_data.h"
#include "yamlsettings/yamlsettings_validators.h"
#include "grasp_output_tile_ascii_settings.h"

grasp_output_tile_function_t grasp_output_tile_function_ascii(){
    grasp_output_tile_function_t x;
    
    x.init=grasp_output_tile_function_ascii_init;
    x.function=grasp_output_tile_function_ascii_process;
    x.close=grasp_output_tile_function_ascii_close;
    
    return x;
}

int grasp_output_tile_function_ascii_init(grasp_settings *settings, grasp_tile_description_t *input_information){
    return 0;
}

int grasp_output_tile_function_ascii_close(void){
    return 0;
}

int grasp_output_tile_function_ascii_process(grasp_output_stream *stream, grasp_settings *settings, grasp_tile_description_t *tile_description, grasp_results_t *results){
    int itime, ix, iy, iwln,inoise, iparam, isd, ipixel, ihlv, sd_index;
    float aerosol_concentration;
    char str_date_time[20 + 1];
    char param_title[64];
    FILE *f;

    f=grasp_output_stream_open(stream, settings, NULL, NULL, &tile_description->dimensions, -1, -1, -1);
    if(grasp_output_stream_writable(stream)==false){ // If it is not writable, close stream and finish function
        grasp_output_stream_close(stream);
        grasp_output_stream_destroy(stream); 
        return 0;
    }
    
    if(grasp_output_tile_information_tile_npixels(results)<=0){
        fprintf(stderr, "The tile is empty. Results with ascii tile output function will not be printed in %s.\n", stream->filename );
    }else{    
        fprintf(f,"row col datetime unixtimestamp lon lat niterations");


        for (inoise = 0; inoise < settings->retrieval.NOISE.INOISE; inoise++) {
            fprintf(f," residual_relative_noise%d", inoise);
        }    

        fprintf(f, " land_percent");        
        fprintf(f, " cloud_mask");

        if(grasp_output_tile_products_retrieval_par(results)){
            for (iparam = 0; iparam < settings->retrieval.KNSING; iparam++) {
                grasp_parameters_get_characteric_type_pretty_name_by_parameter_number(&(settings->retrieval.NDIM), iparam, false, settings->retrieval.WAVE, settings->retrieval.IWW_SINGL, param_title, 64);
                fprintf(f, " %s", param_title);
            }               
        }

        if(grasp_output_tile_products_aerosol_opt(results)){
            fprintf(f," AExp");
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," tau%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," ssa%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," aaod%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }            
            
            if(settings->retrieval.NSD>1){
                for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                    for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                        fprintf(f," tau%d_%d", (int)(settings->retrieval.WAVE[iwln]*1000),isd);
                    }
                }
                for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                    for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                        fprintf(f," ssa%d_%d", (int)(settings->retrieval.WAVE[iwln]*1000),isd);
                    }
                }
            }
        }
                
        if(grasp_output_tile_products_aerosol_sd2m_ext(results)){
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," tauF%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }

            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," tauC%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }  
        }
        
        if(grasp_output_tile_products_aerosol_lidar(results)){
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," lidar_ratio%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }
        }
        
        if(grasp_output_tile_products_forcing_bbflux(results)){
           ipixel=0;
           while (results->tile_result_map[ipixel] == NULL) ipixel++;
           for (ihlv = 0; ihlv < grasp_output_tile_forcing_bbflux_nhlv(results,0,0,ipixel); ihlv++) {
              fprintf(f," bbufx0%0.0fkm", grasp_output_tile_forcing_bbflux_hlv(results,0,0,ipixel, ihlv));
           }
           for (ihlv = 0; ihlv < grasp_output_tile_forcing_bbflux_nhlv(results,0,0,ipixel); ihlv++) {
              fprintf(f," bbdfx0%0.0fkm", grasp_output_tile_forcing_bbflux_hlv(results,0,0,ipixel, ihlv));
           }
           for (ihlv = 0; ihlv < grasp_output_tile_forcing_bbflux_nhlv(results,0,0,ipixel); ihlv++) {
              fprintf(f," bbufxa%0.0fkm", grasp_output_tile_forcing_bbflux_hlv(results,0,0,ipixel, ihlv));
           }
           for (ihlv = 0; ihlv < grasp_output_tile_forcing_bbflux_nhlv(results,0,0,ipixel); ihlv++) {
              fprintf(f," bbdfxa%0.0fkm", grasp_output_tile_forcing_bbflux_hlv(results,0,0,ipixel, ihlv));
           }
        }
        
        if(grasp_output_tile_products_forcing_forcing(results)){
           ipixel=0;
           while (results->tile_result_map[ipixel] == NULL) ipixel++;
           for (ihlv = 0; ihlv < grasp_output_tile_forcing_forcing_nhlv(results,0,0,ipixel); ihlv++) {
              fprintf(f," netforc%0.0fkm", grasp_output_tile_forcing_forcing_hlv(results,0,0,ipixel, ihlv));
           }
           for (ihlv = 0; ihlv < grasp_output_tile_forcing_forcing_nhlv(results,0,0,ipixel); ihlv++) {
              fprintf(f," forceff%0.0fkm", grasp_output_tile_forcing_forcing_hlv(results,0,0,ipixel, ihlv));
           }
        }

        if(grasp_output_tile_products_surface_surf(results)){
            fprintf(f," ndvi");
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," salbedo%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }
        }        
        
        if(grasp_output_tile_products_aerosol_pm(results)){
            for (iparam = 0; iparam < settings->retrieval.nPM_diam; iparam++) {
                fprintf(f," PM(%f)", settings->retrieval.PM_diam[iparam]);
            }
        }
        
        if(grasp_output_tile_products_aerosol_types(results)){
            fprintf(f," Aerosol_type_index");
        }
        
        if(grasp_output_tile_products_aerosol_rind(results)){
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                    fprintf(f," reff_index_real%d_%d", (int)(settings->retrieval.WAVE[iwln]*1000),isd);
                }
            }
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                    fprintf(f," reff_index_imag%d_%d", (int)(settings->retrieval.WAVE[iwln]*1000),isd);
                }
            }
        }
        
        if(grasp_output_tile_products_aerosol_chem(results)){
            for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                fprintf(f," chem_aer_relative_humidity_%d", isd);
            }
            for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                fprintf(f," chem_aer_water_fraction_%d", isd);
            }
            for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                fprintf(f," chem_aer_soluble_fraction_%d", isd);
            }
            for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                fprintf(f," chem_aer_insoluble_fraction_%d", isd);
            }
            for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                fprintf(f," chem_aer_soot_fraction_%d", isd);
            }
            for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                fprintf(f," chem_aer_iron_fraction_%d", isd);
            }

/* modify by lei on 15/11/2016 */
            for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                fprintf(f," chem_aer_brc_fraction_%d", isd);
            }
/* modify by lei on 15/11/2016 */



            if(settings->output.tile_function.ascii.chemical_concentration==true){
                for(isd = 0; isd<grasp_parameters_number_of_modes_of_parameter(&(settings->retrieval.NDIM), par_type_CXRI_nmix);isd++){ 
                    for (iparam = 0; iparam < grasp_parameters_number_of_elements_of_parameter(&(settings->retrieval.NDIM), par_type_CXRI_nmix, isd+1); iparam++) {
                        fprintf(f," chemical_elements_concentration%d_%d", iparam+1,isd+1); 
                    }
                }
            }

/* modify by lei on 15/11/2016 */
            if(settings->output.tile_function.ascii.chemical_concentration==true){

                      for(isd = 0; isd<grasp_parameters_number_of_modes_of_parameter(&(settings->retrieval.NDIM), par_type_CXRI_chem);isd++){

                       fprintf(f," chem_aer_water_volume_concentration_%d", isd);
                       fprintf(f," chem_aer_soluble_volume_concentration_%d", isd);
                       fprintf(f," chem_aer_insoluble_volume_concentration_%d", isd);
                       fprintf(f," chem_aer_soot_volume_concentration_%d", isd);
                       fprintf(f," chem_aer_iron_volume_concentration_%d", isd);
                       fprintf(f," chem_aer_brc_volume_concentration_%d", isd);
                }
            }
/* modify by lei on 15/11/2016 */

            
            
        }
        
        fprintf(f,"\n");

        for (itime = 0; itime < tile_description->dimensions.tile_nt; itime++) {
            for (ix = 0; ix < tile_description->dimensions.tile_nx; ix++) {
                for (iy = 0; iy < tile_description->dimensions.tile_ny; iy++) {
                    if( results->tile_result_map[index3D(itime,ix,iy,tile_description->dimensions.tile_nx,tile_description->dimensions.tile_ny)] != NULL ){
                        time_to_string_r((time_t) grasp_output_tile_pixel_information_time(results,itime,ix,iy), "%FT%H:%M:%SZ", sizeof (str_date_time), str_date_time);
                        fprintf(f, "%d %d %s ", grasp_output_tile_pixel_information_out_x(results,itime,ix,iy), grasp_output_tile_pixel_information_out_y(results,itime,ix,iy), str_date_time);
                        fprintf(f, "%" PRId64 " %f %f", grasp_output_tile_pixel_information_time(results,itime,ix,iy), grasp_output_tile_pixel_information_longitude(results,itime,ix,iy),grasp_output_tile_pixel_information_latitude(results,itime,ix,iy));
                        fprintf(f, " %d", grasp_output_tile_retrieval_res_niter(results,itime,ix,iy));

                        for (inoise = 0; inoise < settings->retrieval.NOISE.INOISE; inoise++) {
                            fprintf(f," %f", grasp_output_tile_retrieval_res_resr(results,itime,ix,iy,inoise));
                        }        

                        fprintf(f, " %f", grasp_output_tile_pixel_information_land_percent(results,itime,ix,iy));
                        fprintf(f, " %d", grasp_output_tile_pixel_information_cloud_flag(results,itime,ix,iy));

                        if(grasp_output_tile_products_retrieval_par(results)){
                            for (iparam = 0; iparam < settings->retrieval.KNSING; iparam++) {
                                fprintf(f, " %f", grasp_output_tile_retrieval_par_parameters(results,itime,ix,iy,iparam));
                            }  
                        }
                        
                        if(grasp_output_tile_products_aerosol_opt(results)){
                            fprintf(f, " %f", grasp_output_tile_aerosol_opt_aexp(results, itime,ix,iy));
                            for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_opt_extt(results, itime,ix,iy,iwln));
                            }
                            for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_opt_ssat(results, itime,ix,iy,iwln));
                            }
                            for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_opt_aextt(results, itime,ix,iy,iwln));
                            }
                            if(settings->retrieval.NSD>1){
                                for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                    for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                        fprintf(f, " %f", grasp_output_tile_aerosol_opt_ext(results, itime,ix,iy,iwln,isd));   
                                    }   
                                }
                                for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                    for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                        fprintf(f, " %f", grasp_output_tile_aerosol_opt_ssa(results, itime,ix,iy,iwln,isd));   
                                    }   
                                }
                            }
                        }
                        if(grasp_output_tile_products_aerosol_sd2m_ext(results)){
                            for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_sd2m_opt_ext(results, itime,ix,iy,iwln,0)); 
                            }        
                            for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_sd2m_opt_ext(results, itime,ix,iy,iwln,1));
                            } 
                        }
                        if(grasp_output_tile_products_aerosol_lidar(results)){ 
                            for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                fprintf(f," %f", grasp_output_tile_aerosol_lidar_lrt(results, itime,ix,iy,iwln));
                            }
                        }
                        if(grasp_output_tile_products_forcing_bbflux(results)){
                            for (ihlv = 0; ihlv < grasp_output_tile_forcing_bbflux_nhlv(results, itime,ix,iy); ihlv++) {
                                fprintf(f, " %f", grasp_output_tile_forcing_bbflux_bbufx0(results, itime,ix,iy,ihlv));
                            }
                        }
                        if(grasp_output_tile_products_forcing_bbflux(results)){
                            for (ihlv = 0; ihlv < grasp_output_tile_forcing_bbflux_nhlv(results, itime,ix,iy); ihlv++) {
                                fprintf(f, " %f", grasp_output_tile_forcing_bbflux_bbdfx0(results, itime,ix,iy,ihlv));
                            }
                        }
                        if(grasp_output_tile_products_forcing_bbflux(results)){
                            for (ihlv = 0; ihlv < grasp_output_tile_forcing_bbflux_nhlv(results, itime,ix,iy); ihlv++) {
                                fprintf(f, " %f", grasp_output_tile_forcing_bbflux_bbufxa(results, itime,ix,iy,ihlv));
                            }
                        }
                        if(grasp_output_tile_products_forcing_bbflux(results)){
                            for (ihlv = 0; ihlv < grasp_output_tile_forcing_bbflux_nhlv(results, itime,ix,iy); ihlv++) {
                                fprintf(f, " %f", grasp_output_tile_forcing_bbflux_bbdfxa(results, itime,ix,iy,ihlv));
                            }
                        }
                        if(grasp_output_tile_products_forcing_forcing(results)){
                            for (ihlv = 0; ihlv < grasp_output_tile_forcing_forcing_nhlv(results, itime,ix,iy); ihlv++) {
                                fprintf(f, " %f", grasp_output_tile_forcing_forcing_netforc(results, itime,ix,iy,ihlv));
                            }
                        }
                        if(grasp_output_tile_products_forcing_forcing(results)){
                            for (ihlv = 0; ihlv < grasp_output_tile_forcing_forcing_nhlv(results, itime,ix,iy); ihlv++) {
                                fprintf(f, " %f", grasp_output_tile_forcing_forcing_forceff(results, itime,ix,iy,ihlv));
                            }
                        }
                        if(grasp_output_tile_products_surface_surf(results)){
                            fprintf(f," %f", grasp_output_tile_surface_ndvi(results, itime,ix,iy));
                            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                                fprintf(f," %f", grasp_output_tile_surface_salbedo(results, itime,ix,iy,iwln));
                            }
                        }
                        
                        if(grasp_output_tile_products_aerosol_pm(results)){
                            for (iparam = 0; iparam < grasp_output_tile_information_npm_diam(results); iparam++) {
                                fprintf(f," %f", grasp_output_tile_aerosol_pm_pm(results, itime,ix,iy,iparam));
                            }
                        }
                        
                        if(grasp_output_tile_products_aerosol_types(results)){
                            fprintf(f," %d", grasp_output_tile_aerosol_types_index(results, itime,ix,iy));
                        }
                        
                        if(grasp_output_tile_products_aerosol_rind(results)){
                            for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                    fprintf(f, " %f", grasp_output_tile_aerosol_rind_mreal(results, itime,ix,iy,iwln,isd));
                                }
                            }
                            for (iwln = 0; iwln < grasp_output_tile_pixel_information_nwl(results, itime,ix,iy); iwln++) {
                                for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                    fprintf(f, " %f", grasp_output_tile_aerosol_rind_mimag(results, itime,ix,iy,iwln,isd));
                                }
                            }
                        }
                        
                        if(grasp_output_tile_products_aerosol_chem(results)){
                            for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_chem_rh(results, itime,ix,iy,isd));
                            }
                            for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_chem_fwtr(results, itime,ix,iy,isd));
                            }
                            for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_chem_fslbl(results, itime,ix,iy,isd));
                            }
                            for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_chem_finslbl(results, itime,ix,iy,isd));
                            }
                            for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_chem_fsoot(results, itime,ix,iy,isd));
                            }
                            for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_chem_firon(results, itime,ix,iy,isd));
                            }


/* modify by lei on 15/11/2016 */
                            for (isd = 0; isd < grasp_output_tile_information_nsd(results); isd++) {
                                fprintf(f, " %f", grasp_output_tile_aerosol_chem_fbrc(results, itime,ix,iy,isd));
                            }
/* modify by lei on 15/11/2016 */


                            
                            
                            if(settings->output.tile_function.ascii.chemical_concentration==true){
                                // For each mode
                                for(isd = 0; isd<grasp_parameters_number_of_modes_of_parameter(&(settings->retrieval.NDIM), par_type_CXRI_nmix);isd++){   
                                    // Calculation of aerosol concentratio of this mode
                                    sd_index=grasp_parameters_characteristic_code_present_of_kind_of_parameter(&(settings->retrieval.NDIM) , par_type_SD_beg, par_type_SD_end);
                                    aerosol_concentration=0;
                                    for (iparam = 0; iparam < grasp_parameters_number_of_elements_of_parameter(&(settings->retrieval.NDIM), sd_index, isd+1); iparam++) {
                                        aerosol_concentration+=grasp_output_tile_retrieval_par_parameters (results, itime,ix,iy, grasp_parameters_get_position(&(settings->retrieval.NDIM), sd_index, isd+1 , iparam));
                                    }
                                    // And for each element of par_type_CXRI_nmix
                                    for (iparam = 0; iparam < grasp_parameters_number_of_elements_of_parameter(&(settings->retrieval.NDIM), par_type_CXRI_nmix, isd+1); iparam++) {
                                        // We print it
                                        fprintf(f," %f", aerosol_concentration*grasp_output_tile_retrieval_par_parameters (results, itime,ix,iy, grasp_parameters_get_position(&(settings->retrieval.NDIM), par_type_CXRI_nmix, isd+1 , iparam))); 
                                    }
                                }
                            }
                            

/* modify by lei on 15/11/2016 */
                            if(settings->output.tile_function.ascii.chemical_concentration==true){
                                // For each mode
                                for(isd = 0; isd<grasp_parameters_number_of_modes_of_parameter(&(settings->retrieval.NDIM), par_type_CXRI_chem);isd++){
                                    // Calculation of aerosol concentratio of this mode
                                    sd_index=grasp_parameters_characteristic_code_present_of_kind_of_parameter(&(settings->retrieval.NDIM) , par_type_SD_beg, par_type_SD_end);
                                    aerosol_concentration=0;

                                    for (iparam = 0; iparam < grasp_parameters_number_of_elements_of_parameter(&(settings->retrieval.NDIM), sd_index, isd+1); iparam++) {
                                        aerosol_concentration+=grasp_output_tile_retrieval_par_parameters (results, itime,ix,iy, grasp_parameters_get_position(&(settings->retrieval.NDIM), sd_index, isd+1 , iparam));
                                    }


                                        fprintf(f, " %f", aerosol_concentration*grasp_output_tile_aerosol_chem_fwtr(results, itime,ix,iy, isd));
                                        fprintf(f, " %f", aerosol_concentration*grasp_output_tile_aerosol_chem_fslbl(results, itime,ix,iy,isd));
                                        fprintf(f, " %f", aerosol_concentration*grasp_output_tile_aerosol_chem_finslbl(results, itime,ix,iy,isd));
                                        fprintf(f, " %f", aerosol_concentration*grasp_output_tile_aerosol_chem_fsoot(results, itime,ix,iy,isd));
                                        fprintf(f, " %f", aerosol_concentration*grasp_output_tile_aerosol_chem_firon(results, itime,ix,iy,isd));
                                        fprintf(f, " %f", aerosol_concentration*grasp_output_tile_aerosol_chem_fbrc(results, itime,ix,iy,isd));

                                }
                            }
                            
/* modify by lei on 15/11/2016 */
                            
                            
                            
                        }
                        
                        fprintf(f, "\n");
                    }
                }
            }
        }

        grasp_output_stream_close(stream);
    }    
    
    return 0;
}

grasp_settings_parameter_array *grasp_output_tile_function_settings_ascii(grasp_settings *settings){
    int i;
    
    // Static definition of a dictionary
    yamlsettings_parameter parameters[]= {
         //                 name                                                                                                                                                                                   , memory direction                                , counter memory direction              , func to set variable type (length of value)             , initial value                               ,  number of elements          , allow_array             , parameter description                                                                                                                                                                                                                                                                                                                ,                 validator 1                 ,                   validator 2                         ,                    validator 3                            , {input function, output function, assigned}   
        {"chemical_concentration"      , &settings->output.tile_function.ascii.chemical_concentration             , NULL       , YS_DATA_TYPE_BOOLEAN              , {ys_d_all,{"false"}}          , {0,{}}     , YS_PARAM_SCALAR         , "Calculate and print chemical concentration"                                        , {  },YAMLSETTINGS_END_VAR  },                  
    };
    grasp_settings_parameter_array *result;
    
    result = (grasp_settings_parameter_array *) malloc(sizeof (grasp_settings_parameter_array));
    
    result->nparameters=sizeof(parameters)/sizeof(yamlsettings_parameter);
    
    result->parameters = (yamlsettings_parameter *) malloc(sizeof (yamlsettings_parameter)*result->nparameters);
    
    for(i=0;i<sizeof(parameters)/sizeof(yamlsettings_parameter);i++){
       yamlsettings_copy_parameter(&(parameters[i]),&(result->parameters[i]));
    }
    
    return result;   
}

