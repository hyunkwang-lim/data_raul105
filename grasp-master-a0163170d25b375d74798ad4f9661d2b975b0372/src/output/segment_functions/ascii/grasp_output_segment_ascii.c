/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include "grasp_output_segment_ascii.h"
#include "../../grasp_output_stream.h"
#include <grasp/utils.h>
#include <inttypes.h>


grasp_output_segment_function_t grasp_output_segment_function_ascii(){
    grasp_output_segment_function_t x;
    
    x.init=grasp_output_segment_function_ascii_init;
    x.function=grasp_output_segment_function_ascii_process;
    x.close=grasp_output_segment_function_ascii_close;
    
    return x;
}

int grasp_output_segment_function_ascii_init(grasp_settings *settings, grasp_tile_description_t *input_information){
    return 0;
}

int grasp_output_segment_function_ascii_close(void){
    return 0;
}

int grasp_output_segment_function_ascii_process(grasp_output_stream *stream, grasp_settings *settings, grasp_segment_t *segment, output_segment_general *output, grasp_tile_description_t *tile_description,int icol,int irow,int itime){
    int ipixel, iwln,inoise, iparam, ihlv, isd;
    char str_date_time[20 + 1];
    char param_title[64];
    FILE *f;

    f=grasp_output_stream_open(stream, settings, segment, output, &tile_description->dimensions, icol, irow, itime);
    if(grasp_output_stream_writable(stream)==false){ // If it is not writable, close stream and finish function
        grasp_output_stream_close(stream);
        grasp_output_stream_destroy(stream); 
        return 0;
    }
        
    if(segment->sdata.npixels<=0){
        fprintf(stderr, "The segment in position time=%d, x=%d and y=%d has no pixels. Results with ascii segment output function will not be printed in %s.\n", itime, icol, irow, stream->filename );
    }else{    

        fprintf(f,"row col datetime unixtimestamp lon lat niteraciones");// res(single_pixel)");

        for (inoise = 0; inoise < settings->retrieval.NOISE.INOISE; inoise++) {
            fprintf(f," residual_relative_noise%d", inoise);
        }

        fprintf(f, " land_percent");
        fprintf(f, " cloud_mask");

        if(settings->retrieval.products.retrieval.par==true && grasp_output_segment_products_retrieval_par(output)==true){
            for (iparam = 0; iparam < settings->retrieval.KNSING; iparam++) {
                grasp_parameters_get_characteric_type_pretty_name_by_parameter_number(&(settings->retrieval.NDIM), iparam, false, settings->retrieval.WAVE, settings->retrieval.IWW_SINGL, param_title, 64);
                fprintf(f, " %s", param_title);
            }             
        }
        
        if(settings->retrieval.products.aerosol.opt==true && grasp_output_segment_products_aerosol_opt(output)==true){
            fprintf(f," AExp");
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," tau%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }   
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," ssa%d",(int)(settings->retrieval.WAVE[iwln]*1000));
            }  
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," aaod%d",(int)(settings->retrieval.WAVE[iwln]*1000));
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
        
        if(settings->retrieval.products.aerosol.sd2m_ext==true && grasp_output_segment_products_aerosol_sd2m_ext(output)==true){
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," tauF%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," tauC%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }    
        }
        
        if(settings->retrieval.products.aerosol.lidar && grasp_output_segment_products_aerosol_lidar(output)){ 
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," lidar_ratio%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }
        }
        
        if(settings->retrieval.products.forcing.bbflux && grasp_output_segment_products_forcing_bbflux(output)){
           for (ihlv = 0; ihlv < grasp_output_segment_forcing_bbflux_nhlv(output, 0); ihlv++) {
              fprintf(f," bbufx0%0.0fkm", grasp_output_segment_forcing_bbflux_hlv(output, 0, ihlv));
           }
           for (ihlv = 0; ihlv < grasp_output_segment_forcing_bbflux_nhlv(output, 0); ihlv++) {
              fprintf(f," bbdfx0%0.0fkm", grasp_output_segment_forcing_bbflux_hlv(output, 0, ihlv));
           }
           for (ihlv = 0; ihlv < grasp_output_segment_forcing_bbflux_nhlv(output, 0); ihlv++) {
              fprintf(f," bbufxa%0.0fkm", grasp_output_segment_forcing_bbflux_hlv(output, 0, ihlv));
           }
           for (ihlv = 0; ihlv < grasp_output_segment_forcing_bbflux_nhlv(output, 0); ihlv++) {
              fprintf(f," bbdfxa%0.0fkm", grasp_output_segment_forcing_bbflux_hlv(output, 0, ihlv));
           }
        }
        
        if(settings->retrieval.products.forcing.forcing && grasp_output_segment_products_forcing_forcing(output)){
           for (ihlv = 0; ihlv < grasp_output_segment_forcing_forcing_nhlv(output, 0); ihlv++) {
              fprintf(f," netforc%0.0fkm", grasp_output_segment_forcing_forcing_hlv(output, 0, ihlv));
           }
           for (ihlv = 0; ihlv < grasp_output_segment_forcing_forcing_nhlv(output, 0); ihlv++) {
              fprintf(f," forceff%0.0fkm", grasp_output_segment_forcing_forcing_hlv(output, 0, ihlv));
           }
        }    
                
        if(settings->retrieval.products.surface.surf && grasp_output_segment_products_surface_surf(output)){
            fprintf(f," ndvi");
            for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                fprintf(f," salbedo%d", (int)(settings->retrieval.WAVE[iwln]*1000));
            }
        }
        
        if(settings->retrieval.products.aerosol.pm && grasp_output_segment_products_aerosol_pm(output)){
            for (iparam = 0; iparam < settings->retrieval.nPM_diam; iparam++) {
                fprintf(f," PM(%f)", settings->retrieval.PM_diam[iparam]);
            }
        }
        
        if(settings->retrieval.products.aerosol.types && grasp_output_segment_products_aerosol_types(output)){
            fprintf(f," Aerosol_type_index");
        }
        
        if(grasp_output_segment_products_aerosol_rind(output)){
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
        
        if(grasp_output_segment_products_aerosol_chem(output)){
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
        }
        
        fprintf(f,"\n");

        for (ipixel = 0; ipixel < segment->sdata.npixels; ipixel++) {
            time_to_string_r((time_t)segment->sdata.pixel[ipixel].t, "%FT%H:%M:%SZ" , sizeof(str_date_time), str_date_time);        
            fprintf(f,"%d %d %s ", segment->sdata.pixel[ipixel].irow, segment->sdata.pixel[ipixel].icol, str_date_time);
            fprintf(f,"%" PRId64 " %f %f %d", segment->sdata.pixel[ipixel].t, segment->sdata.pixel[ipixel].x, segment->sdata.pixel[ipixel].y, grasp_output_segment_retrieval_res_niter(output));        

            for (inoise = 0; inoise < settings->retrieval.NOISE.INOISE; inoise++) {
                fprintf(f," %f", grasp_output_segment_retrieval_res_pixel_resr(output,ipixel,inoise));
            }        
            
            fprintf(f, " %f", segment->sdata.pixel[ipixel].land_percent);
            if(segment->sdata.pixel[ipixel].cloudy==0){
                fprintf(f, " 1");
            }else if(segment->sdata.pixel[ipixel].cloudy==1){
                fprintf(f, " 0");
            }else{
                fprintf(f, " %d", segment->sdata.pixel[ipixel].cloudy);
            }
            
            if(settings->retrieval.products.retrieval.par==true && grasp_output_segment_products_retrieval_par(output)==true){
                for (iparam = 0; iparam < settings->retrieval.KNSING; iparam++) {
                    fprintf(f, " %f", grasp_output_segment_retrieval_par_parameters(output,ipixel,iparam));
                }  
            }            
            
            if(settings->retrieval.products.aerosol.opt==true && grasp_output_segment_products_aerosol_opt(output)==true){
                fprintf(f," %f", grasp_output_segment_aerosol_opt_aexp(output,ipixel));
                for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                    fprintf(f," %f", grasp_output_segment_aerosol_opt_extt(output,ipixel,iwln));
                }
                for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                    fprintf(f," %f", grasp_output_segment_aerosol_opt_ssat(output,ipixel,iwln));
                }   
                for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                    fprintf(f," %f", grasp_output_segment_aerosol_opt_aextt(output,ipixel,iwln));
                }     
                if(settings->retrieval.NSD>1){
                    for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                        for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                            fprintf(f, " %f", grasp_output_segment_aerosol_opt_ext(output,ipixel,iwln,isd));   
                        }   
                    }
                    for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                        for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                            fprintf(f, " %f", grasp_output_segment_aerosol_opt_ssa(output,ipixel,iwln,isd));   
                        }   
                    }
                }
            }
            
            if(settings->retrieval.products.aerosol.sd2m_ext==true && grasp_output_segment_products_aerosol_sd2m_ext(output)){
                for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                    fprintf(f," %f", grasp_output_segment_aerosol_sd2m_opt_ext(output,ipixel, iwln, 0) );
                }        
                for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                    fprintf(f," %f", grasp_output_segment_aerosol_sd2m_opt_ext(output,ipixel, iwln, 1));
                }        
            }
            
            if(settings->retrieval.products.aerosol.lidar && grasp_output_segment_products_aerosol_lidar(output)){ 
                for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                    fprintf(f," %f", grasp_output_segment_aerosol_lidar_lrt(output,ipixel,iwln));
                }
            }
            
            if(settings->retrieval.products.forcing.bbflux==true && grasp_output_segment_products_forcing_bbflux(output)==true){
               for (ihlv = 0; ihlv < grasp_output_segment_forcing_forcing_nhlv(output, 0); ihlv++) {
                    fprintf(f, " %f",  grasp_output_segment_forcing_bbflux_bbufx0(output,ipixel,ihlv));
               }
               for (ihlv = 0; ihlv < grasp_output_segment_forcing_forcing_nhlv(output, 0); ihlv++) {
                    fprintf(f, " %f",  grasp_output_segment_forcing_bbflux_bbdfx0(output,ipixel,ihlv));
               }
               for (ihlv = 0; ihlv < grasp_output_segment_forcing_forcing_nhlv(output, 0); ihlv++) {
                    fprintf(f, " %f",  grasp_output_segment_forcing_bbflux_bbufxa(output,ipixel,ihlv));
               }
               for (ihlv = 0; ihlv < grasp_output_segment_forcing_forcing_nhlv(output, 0); ihlv++) {
                    fprintf(f, " %f",  grasp_output_segment_forcing_bbflux_bbdfxa(output,ipixel,ihlv));
               }
            }
            
            if(settings->retrieval.products.forcing.forcing==true && grasp_output_segment_products_forcing_forcing(output)==true){
               for (ihlv = 0; ihlv < grasp_output_segment_forcing_forcing_nhlv(output, 0); ihlv++) {
                    fprintf(f, " %f",  grasp_output_segment_forcing_forcing_netforc(output,ipixel,ihlv));
               }
               for (ihlv = 0; ihlv < grasp_output_segment_forcing_forcing_nhlv(output, 0); ihlv++) {
                    fprintf(f, " %f",  grasp_output_segment_forcing_forcing_forceff(output,ipixel,ihlv));
               }
            }   
            
            if(settings->retrieval.products.surface.surf && grasp_output_segment_products_surface_surf(output)){
                fprintf(f," %f", grasp_output_segment_surface_ndvi(output, ipixel));
                for (iwln = 0; iwln < settings->retrieval.NW; iwln++) {
                    fprintf(f," %f", grasp_output_segment_surface_salbedo(output, ipixel, iwln));
                }
            }

            if(settings->retrieval.products.aerosol.pm && grasp_output_segment_products_aerosol_pm(output)){
                for (iparam = 0; iparam < settings->retrieval.nPM_diam; iparam++) {
                    fprintf(f," %f", grasp_output_segment_aerosol_pm_pm(output, ipixel, iparam));
                }
            }
        
            if(settings->retrieval.products.aerosol.types && grasp_output_segment_products_aerosol_types(output)){
                fprintf(f," %d", grasp_output_segment_aerosol_types_index(output, ipixel));
            }

            if(grasp_output_segment_products_aerosol_rind(output)){
                for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                    for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                        fprintf(f, " %f", grasp_output_segment_aerosol_rind_mreal(output, ipixel,iwln,isd));
                    }
                }
                for (iwln = 0; iwln < segment->sdata.pixel[ipixel].nwl; iwln++) {
                    for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                        fprintf(f, " %f", grasp_output_segment_aerosol_rind_mimag(output, ipixel,iwln,isd));
                    }
                }
            }

            if(grasp_output_segment_products_aerosol_chem(output)){
                for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                    fprintf(f, " %f", grasp_output_segment_aerosol_chem_rh(output, ipixel,isd));
                }
                for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                    fprintf(f, " %f", grasp_output_segment_aerosol_chem_fwtr(output, ipixel,isd));
                }
                for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                    fprintf(f, " %f", grasp_output_segment_aerosol_chem_fslbl(output, ipixel,isd));
                }
                for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                    fprintf(f, " %f", grasp_output_segment_aerosol_chem_finslbl(output, ipixel,isd));
                }
                for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                    fprintf(f, " %f", grasp_output_segment_aerosol_chem_fsoot(output, ipixel,isd));
                }
                for (isd = 0; isd < settings->retrieval.NSD; isd++) {
                    fprintf(f, " %f", grasp_output_segment_aerosol_chem_firon(output, ipixel,isd));
                }
            }
                        
            
            fprintf(f,"\n");
        }


        grasp_output_stream_close(stream);

    }
        
    return 0;
}


grasp_settings_parameter_array *grasp_output_segment_function_settings_ascii(grasp_settings *settings){
    return NULL;
}
