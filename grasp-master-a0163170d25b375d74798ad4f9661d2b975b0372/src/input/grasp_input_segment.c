/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include "grasp_input_segment.h"


void grasp_segment_remove_pixel(grasp_segment_t *segment, int ipixel){
    int i;
    int pixel_time;
    int pixels_in_time=0;
    bool last_time=false;
    
    // Check if nt is still valid or we have emptied a specific time and we have to reduce it in 1.
    pixel_time=segment->sdata.pixel[ipixel].it;
    for(i=0;i<segment->sdata.npixels;i++){
        if(segment->sdata.pixel[i].it==pixel_time){
            pixels_in_time++;
        }
    }
    assert(pixels_in_time>=1);
    if(pixels_in_time==1){ // The pixel we are going to remove is the only pixel in the current time so we have to remove completely this time
        last_time=true;
        if(segment->edges.N_I_EDGE>0){ 
            // WARNING: This function does not support edges yet.
            assert(1!=1);
        }
        segment->sdata.nt--;
    }
    
    // Remove the pixel
    segment->sdata.npixels--;
    
    for(i=ipixel;i<segment->sdata.npixels;i++){
        memcpy(&(segment->sdata.pixel[i]), &(segment->sdata.pixel[i+1]), sizeof(pixel_t));
        if(last_time==true){
            segment->sdata.pixel[i].it--;
        }
        memcpy(segment->iguess[i], segment->iguess[i+1], sizeof(float)*_KPARS);
    }
}

void grasp_segment_remove_pixels(grasp_segment_t *segment, int npixels, int *pixels){
    int i;
    
    for(i=0;i<npixels-1;i++){ // To ensure monitonic incrasing order
        assert(pixels[i]<pixels[i+1]);
    }
    
    for(i=0;i<npixels;i++){
        grasp_segment_remove_pixel(segment,pixels[i]-i);
    }
}

void grasp_segment_set_nt(sensor_data_t *sdata, int nt){
    assert(nt>=0);
    assert(nt<=_KITIME);
    
    sdata->nt=nt;
}

void grasp_segment_set_nx(sensor_data_t *sdata, int nx){
    assert(nx>=0);
    assert(nx<=_KIX);
    
    sdata->nx=nx;    
}

void grasp_segment_set_ny(sensor_data_t *sdata, int ny){
    assert(ny>=0);
    assert(ny<=_KIY);
    
    sdata->ny=ny;    
}

void grasp_segment_set_npixels(sensor_data_t *sdata, int npixels){
    assert(npixels>=0);
    assert(npixels<=_KIMAGE);
    
    sdata->npixels=npixels;
}

void grasp_segment_set_pixel_it(sensor_data_t *sdata,int ipixel,int it){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE);
    assert(it>0);
    assert(it<=_KITIME);
    
    sdata->pixel[ipixel].it=it;
}

void grasp_segment_set_pixel_ix(sensor_data_t *sdata,int ipixel,int ix){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE);
    assert(ix>0);
    assert(ix<=_KIX);
    
    sdata->pixel[ipixel].ix=ix;    
}

void grasp_segment_set_pixel_iy(sensor_data_t *sdata,int ipixel,int iy){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE);
    assert(iy>0);
    assert(iy<=_KIY);
    

    sdata->pixel[ipixel].iy=iy;    
}

void grasp_segment_set_pixel_cloudy(sensor_data_t *sdata,int ipixel,int cloudy){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE);    
    assert(cloudy==0 || cloudy==1);
    
    sdata->pixel[ipixel].cloudy=cloudy;
}

void grasp_segment_set_pixel_irow(sensor_data_t *sdata,int ipixel,int irow){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE);  
    
    sdata->pixel[ipixel].irow=irow;
}

void grasp_segment_set_pixel_icol(sensor_data_t *sdata,int ipixel,int icol){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE);    
    
    sdata->pixel[ipixel].icol=icol;
}

void grasp_segment_set_pixel_file_index(sensor_data_t *sdata,int ipixel,int file_index){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(file_index>=-1);
    
    sdata->pixel[ipixel].file_index=file_index;
}

void grasp_segment_set_pixel_x(sensor_data_t *sdata,int ipixel,float x){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(x>=-180.0);
    assert(x<=360.0);
    
    sdata->pixel[ipixel].x=x;  
}

void grasp_segment_set_pixel_y(sensor_data_t *sdata,int ipixel,float y){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(y>=-90.0);
    assert(y<=90.0);
    
    sdata->pixel[ipixel].y=y;    
}

void grasp_segment_set_pixel_t(sensor_data_t *sdata,int ipixel,int64_t t){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    
    
    sdata->pixel[ipixel].t=t;    
}
void grasp_segment_set_pixel_out_t(sensor_data_t *sdata,int ipixel,int out_t){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    
    sdata->pixel[ipixel].out_t=out_t;    
}

void grasp_segment_set_pixel_out_x(sensor_data_t *sdata,int ipixel,int out_x){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    
    sdata->pixel[ipixel].out_x=out_x;    
}

void grasp_segment_set_pixel_out_y(sensor_data_t *sdata,int ipixel,int out_y){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    
    sdata->pixel[ipixel].out_y=out_y;    
}

void grasp_segment_set_pixel_masl(sensor_data_t *sdata,int ipixel,float masl){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    
    sdata->pixel[ipixel].masl=masl;
}

void grasp_segment_set_pixel_hobs(sensor_data_t *sdata,int ipixel,float hobs){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    
    sdata->pixel[ipixel].hobs=hobs;
}

void grasp_segment_set_pixel_land_percent(sensor_data_t *sdata,int ipixel,float land_percent){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(land_percent>=0.0);
    assert(land_percent<=100.0);
    
    sdata->pixel[ipixel].land_percent=land_percent;
}

void grasp_segment_set_pixel_nwl(sensor_data_t *sdata,int ipixel,int nwl){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(nwl>0);  
    assert(nwl<=_KWM);
    
    sdata->pixel[ipixel].nwl=nwl;
}

void grasp_segment_set_pixel_ifgas(sensor_data_t *sdata,int ipixel,int ifgas){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(ifgas==0 || ifgas==1);
    
    sdata->pixel[ipixel].ifgas=ifgas;
}

void grasp_segment_set_pixel_hvp(sensor_data_t *sdata,int ipixel, int ivm, float hvp){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(ivm>=0);
    assert(ivm<=_KVERTM);
    
    sdata->pixel[ipixel].hvp[ivm]=hvp;
}

void grasp_segment_set_pixel_meas_wl(sensor_data_t *sdata,int ipixel, int iwl, float wl){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM); 
    
    sdata->pixel[ipixel].meas[iwl].wl=wl;     
}

void grasp_segment_set_pixel_meas_ind_wl(sensor_data_t *sdata,int ipixel, int iwl, float ind_wl){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM); 
    
    sdata->pixel[ipixel].meas[iwl].ind_wl=ind_wl;     
}

void grasp_segment_set_pixel_meas_nsurf(sensor_data_t *sdata,int ipixel, int iwl, int nsurf){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);     
    assert(nsurf>=0);
    assert(nsurf<=_KSURF);
    
    sdata->pixel[ipixel].meas[iwl].nsurf=nsurf;
}

void grasp_segment_set_pixel_meas_gaspar(sensor_data_t *sdata,int ipixel, int iwl, float gaspar){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);     
    
    sdata->pixel[ipixel].meas[iwl].gaspar=gaspar;
}

void grasp_segment_set_pixel_meas_sza(sensor_data_t *sdata,int ipixel, int iwl, float sza){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);     
    
    sdata->pixel[ipixel].meas[iwl].sza=sza;
}

void grasp_segment_set_pixel_meas_groundpar(sensor_data_t *sdata,int ipixel, int iwl, int isurf, float groundpar){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(isurf>=0);
    assert(isurf<_KSURF);
    
    sdata->pixel[ipixel].meas[iwl].groundpar[isurf]=groundpar;
}

void grasp_segment_set_pixel_meas_nip(sensor_data_t *sdata,int ipixel, int iwl, int nip){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(nip>0);    
    assert(nip<=_KIP);
    
    sdata->pixel[ipixel].meas[iwl].nip=nip;
}

void grasp_segment_set_pixel_meas_meas_type(sensor_data_t *sdata,int ipixel, int iwl, int ip, int meas_type){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ip>=0);    
    assert(ip<_KIP);

    sdata->pixel[ipixel].meas[iwl].meas_type[ip]=meas_type;
}

void grasp_segment_set_pixel_meas_nbvm(sensor_data_t *sdata,int ipixel, int iwl, int ip, int nbvm){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ip>=0);    
    assert(ip<_KIP);
    assert(nbvm>=0);
    assert(nbvm<=_KNBVM);
    
    sdata->pixel[ipixel].meas[iwl].nbvm[ip]=nbvm;
}

void grasp_segment_set_pixel_meas_ifcov(sensor_data_t *sdata,int ipixel, int iwl, int ip, int ifcov){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ip>=0);    
    assert(ip<_KIP);   
    assert(ifcov==0 || ifcov==1);
    
    sdata->pixel[ipixel].meas[iwl].ifcov[ip]=ifcov;
}

void grasp_segment_set_pixel_meas_ifmp(sensor_data_t *sdata,int ipixel, int iwl, int ip, int ifmp){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ip>=0);    
    assert(ip<_KIP);
    assert(ifmp==0 || ifmp==1);
    
    sdata->pixel[ipixel].meas[iwl].ifmp[ip]=ifmp;
}

void grasp_segment_set_pixel_meas_thetav(sensor_data_t *sdata,int ipixel, int iwl, int ip, int ivalidmeas, float thetav){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ip>=0);    
    assert(ip<_KIP);
    assert(ivalidmeas>=0);
    assert(ivalidmeas<sdata->pixel[ipixel].meas[iwl].nbvm[ip]);    
    
    sdata->pixel[ipixel].meas[iwl].thetav[ip][ivalidmeas]=thetav;
}

void grasp_segment_set_pixel_meas_phi(sensor_data_t *sdata,int ipixel, int iwl, int ip, int ivalidmeas, float phi){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ip>=0);    
    assert(ip<_KIP);
    assert(ivalidmeas>=0);
    assert(ivalidmeas<sdata->pixel[ipixel].meas[iwl].nbvm[ip]);
    
    sdata->pixel[ipixel].meas[iwl].phi[ip][ivalidmeas]=phi;
}

void grasp_segment_set_pixel_meas_cmtrx(sensor_data_t *sdata,int ipixel, int iwl, int ip, int ivalidmeas, float cmtrx){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ip>=0);    
    assert(ip<_KIP);
    assert(ivalidmeas>=0);
    assert(ivalidmeas<sdata->pixel[ipixel].meas[iwl].nbvm[ip]);
    
    sdata->pixel[ipixel].meas[iwl].cmtrx[ip][ivalidmeas]=cmtrx;
}

void grasp_segment_set_pixel_meas_mprof(sensor_data_t *sdata,int ipixel, int iwl, int ip, int ivalidmeas, float mprof){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ip>=0);    
    assert(ip<_KIP);
    assert(ivalidmeas>=0);
    assert(ivalidmeas<sdata->pixel[ipixel].meas[iwl].nbvm[ip]);  
    assert(ivalidmeas<_KVERTM);

    sdata->pixel[ipixel].meas[iwl].mprof[ip][ivalidmeas]=mprof;
}

void grasp_segment_set_pixel_meas_tau(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float tau){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].tau[ivalidmeas]=tau;
}

void grasp_segment_set_pixel_meas_p11(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float p11){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].p11[ivalidmeas]=p11;
}

void grasp_segment_set_pixel_meas_p12(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float p12){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
        
    sdata->pixel[ipixel].meas[iwl].p12[ivalidmeas]=p12;
}

void grasp_segment_set_pixel_meas_p22(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float p22){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].p22[ivalidmeas]=p22;    
}

void grasp_segment_set_pixel_meas_p33(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float p33){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    assert(p33>=0.0);    
    
    sdata->pixel[ipixel].meas[iwl].p33[ivalidmeas]=p33;
}

void grasp_segment_set_pixel_meas_p34(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float p34){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].p34[ivalidmeas]=p34;
}

void grasp_segment_set_pixel_meas_p44(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float p44){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].p44[ivalidmeas]=p44;    
}

void grasp_segment_set_pixel_meas_i(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float i){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].i[ivalidmeas]=i;    
}

void grasp_segment_set_pixel_meas_q(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float q){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].q[ivalidmeas]=q;    
}


void grasp_segment_set_pixel_meas_u(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float u){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].u[ivalidmeas]=u;    
}

void grasp_segment_set_pixel_meas_p(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float p){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_NBVM);  
    
    sdata->pixel[ipixel].meas[iwl].p[ivalidmeas]=p;    
}

void grasp_segment_set_pixel_meas_ls(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float ls){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_KVERTM);  
    
    sdata->pixel[ipixel].meas[iwl].ls[ivalidmeas]=ls;      
}

void grasp_segment_set_pixel_meas_dp(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float dp){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_KVERTM);  
    
    sdata->pixel[ipixel].meas[iwl].dp[ivalidmeas]=dp;  
}

void grasp_segment_set_pixel_meas_rl(sensor_data_t *sdata,int ipixel, int iwl, int ivalidmeas, float rl){
    assert(ipixel>=0);
    assert(ipixel<_KIMAGE); 
    assert(iwl>=0);
    assert(iwl<_KWM);  
    assert(ivalidmeas>=0);
    assert(ivalidmeas<_KVERTM);  
    
    sdata->pixel[ipixel].meas[iwl].rl[ivalidmeas]=rl;  
}
