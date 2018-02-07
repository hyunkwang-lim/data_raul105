/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "grasp_output_tile_result.h"
#include <grasp/utils.h>

#define ipixel_ index3D(it,ix,iy,output->information.tile_npixels_x,output->information.tile_npixels_y)
#define NSD_ output->information.nsd
#define NMPAR_ output->information.nmpar
#define NWL_ output->tile_result_map[ipixel_]->information.nwl

bool grasp_output_tile_is_pixel(const grasp_results_t *output, int it, int ix, int iy){
    if(output->tile_result_map[ipixel_]!=NULL){
        return true;
    }else{
        return false;
    }
}

bool grasp_output_tile_products_retrieval_res (const grasp_results_t *output) {
    return output->products.retrieval.res;
}

bool grasp_output_tile_products_retrieval_par (const grasp_results_t *output) {
    return output->products.retrieval.par;
}

bool grasp_output_tile_products_retrieval_fit (const grasp_results_t *output) {
    return output->products.retrieval.fit;
}

bool grasp_output_tile_products_aerosol_opt (const grasp_results_t *output) {
    return output->products.aerosol.opt;
}

bool grasp_output_tile_products_aerosol_rind (const grasp_results_t *output) {
    return output->products.aerosol.rind;
}

bool grasp_output_tile_products_aerosol_chem (const grasp_results_t *output) {
    return output->products.aerosol.chem;
}

bool grasp_output_tile_products_aerosol_phmx (const grasp_results_t *output) {
    return output->products.aerosol.phmx;
}

bool grasp_output_tile_products_aerosol_lidar (const grasp_results_t *output) {
    return output->products.aerosol.lidar;
}

bool grasp_output_tile_products_aerosol_sd2m_mph (const grasp_results_t *output) {
    return output->products.aerosol.sd2m_mph;
}

bool grasp_output_tile_products_aerosol_sd2m_ext (const grasp_results_t *output) {
    return output->products.aerosol.sd2m_ext;
}

bool grasp_output_tile_products_aerosol_pm (const grasp_results_t *output) {
    return output->products.aerosol.pm;
}

bool grasp_output_tile_products_aerosol_types (const grasp_results_t *output) {
    return output->products.aerosol.types;
}

bool grasp_output_tile_products_clouds_opt (const grasp_results_t *output) {
    return output->products.clouds.opt;
}

bool grasp_output_tile_products_clouds_rind (const grasp_results_t *output) {
    return output->products.clouds.rind;
}

bool grasp_output_tile_products_clouds_chem (const grasp_results_t *output) {
    return output->products.clouds.chem;
}

bool grasp_output_tile_products_clouds_phmx (const grasp_results_t *output) {
    return output->products.clouds.phmx;
}

bool grasp_output_tile_products_clouds_lidar (const grasp_results_t *output) {
    return output->products.clouds.lidar;
}

bool grasp_output_tile_products_clouds_sd2m_mph (const grasp_results_t *output) {
    return output->products.clouds.sd2m_mph;
}

bool grasp_output_tile_products_clouds_sd2m_ext (const grasp_results_t *output) {
    return output->products.clouds.sd2m_ext;
}

bool grasp_output_tile_products_clouds_pm (const grasp_results_t *output) {
    return output->products.clouds.pm;
}

bool grasp_output_tile_products_clouds_types (const grasp_results_t *output) {
    return output->products.clouds.types;
}

bool grasp_output_tile_products_surface_surf (const grasp_results_t *output) {
    return output->products.surface.surf;
}

bool grasp_output_tile_products_errest_par (const grasp_results_t *output) {
    return output->products.errest.par;
}

bool grasp_output_tile_products_errest_aerosol_opt (const grasp_results_t *output) {
    return output->products.errest.aerosol.opt;
}

bool grasp_output_tile_products_errest_aerosol_lidar (const grasp_results_t *output) {
    return output->products.errest.aerosol.lidar;
}

bool grasp_output_tile_products_errest_clouds_opt (const grasp_results_t *output) {
    return output->products.errest.clouds.opt;
}

bool grasp_output_tile_products_errest_clouds_lidar (const grasp_results_t *output) {
    return output->products.errest.clouds.lidar;
}

bool grasp_output_tile_products_forcing_bbflux (const grasp_results_t *output) {
    return output->products.forcing.bbflux;
}

bool grasp_output_tile_products_forcing_forcing (const grasp_results_t *output) {
    return output->products.forcing.forcing;
}

int grasp_output_tile_information_tile_npixels(const grasp_results_t *output){
    return output->information.tile_npixels;
}

int grasp_output_tile_information_tile_npixels_t(const grasp_results_t *output){
    return output->information.tile_npixels_t;
}

int grasp_output_tile_information_tile_npixels_x(const grasp_results_t *output){
    return output->information.tile_npixels_x;
}

int grasp_output_tile_information_tile_npixels_y(const grasp_results_t *output){
    return output->information.tile_npixels_y;
}

int grasp_output_tile_information_npars(const grasp_results_t *output){
    return output->information.npars;
}

int grasp_output_tile_information_nrr(const grasp_results_t *output){
    return output->information.nrr;
}

int grasp_output_tile_information_nrc(const grasp_results_t *output){
    return output->information.nrc;
}

int grasp_output_tile_information_nmpar(const grasp_results_t *output){
    return output->information.nmpar;
}

int grasp_output_tile_information_nsd(const grasp_results_t *output){
    return output->information.nsd;
}

float grasp_output_tile_information_retrieval_par_radius(const grasp_results_t *output, int irr){
    return output->information.retrieval_par_radius[irr];
}

float grasp_output_tile_information_retrieval_par_SDL(const grasp_results_t *output, int irr, int irc){
    return output->information.retrieval_par_SDL[irc][irr];
}

int grasp_output_tile_information_npm_diam(const grasp_results_t *output){
    return output->information.nPM_diam;
}

int grasp_output_tile_information_phmx_angle(const grasp_results_t *output, int iangle){
    return output->information.phmx_angle[iangle];
}

int grasp_output_tile_information_nnoises(const grasp_results_t *output){
    return output->information.nnoises;
}

int grasp_output_tile_pixel_information_segment_time(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.segment_time;
}

int grasp_output_tile_pixel_information_segment_col(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.segment_col;
}

int grasp_output_tile_pixel_information_segment_row(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.segment_row;
}

int grasp_output_tile_pixel_information_it(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.it;
}

int grasp_output_tile_pixel_information_ix(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.ix;
}

int grasp_output_tile_pixel_information_iy(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.iy;
}

int grasp_output_tile_pixel_information_out_x(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.out_x;
}

int grasp_output_tile_pixel_information_out_y(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.out_y;
}

int grasp_output_tile_pixel_information_out_t(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.out_t;
}

float grasp_output_tile_pixel_information_latitude(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.latitude;
}

float grasp_output_tile_pixel_information_longitude(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.longitude;
}

int grasp_output_tile_pixel_information_grid_col(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.grid_col;
}

int grasp_output_tile_pixel_information_grid_row(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.grid_row;
}

int64_t grasp_output_tile_pixel_information_time(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.time;
}

int grasp_output_tile_pixel_information_nwl(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.nwl;
}

int grasp_output_tile_pixel_information_cloud_flag(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.cloud_flag;
}

float grasp_output_tile_pixel_information_land_percent(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.land_percent;
}

int grasp_output_tile_pixel_information_file_index(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.file_index;
}

float grasp_output_tile_pixel_information_masl(const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->information.masl;
}

float grasp_output_tile_pixel_information_sza(const grasp_results_t *output,int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->information.sza[iwl];
}

int grasp_output_tile_pixel_segment_npixels(const grasp_results_t *output, int it, int ix, int iy, int nrCol, int nrRow){
    int index = (output->tile_result_map[ipixel_]->information.segment_time * nrCol*nrRow +
                 output->tile_result_map[ipixel_]->information.segment_col * nrRow +
                 output->tile_result_map[ipixel_]->information.segment_row);
    return output->segment_result[index].npixel;
}

int grasp_output_tile_retrieval_res_niter (const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->retrieval_res->niter;
}

float grasp_output_tile_retrieval_res_rest (const grasp_results_t *output,int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->retrieval_res->rest;
}

float grasp_output_tile_retrieval_res_resa (const grasp_results_t *output,int it, int ix, int iy, int inoise){
    return output->tile_result_map[ipixel_]->retrieval_res->resa[inoise];
}

float grasp_output_tile_retrieval_res_resr (const grasp_results_t *output,int it, int ix, int iy, int inoise){
    return output->tile_result_map[ipixel_]->retrieval_res->resr[inoise];
}

float grasp_output_tile_retrieval_par_parameters (const grasp_results_t *output, int it, int ix, int iy, int ipar){
    return output->tile_result_map[ipixel_]->retrieval_par->parameters[ipar];
}

const pixel_t *grasp_output_tile_retrieval_fit_pixel_original (const grasp_results_t *output, int it, int ix, int iy){
    return &output->tile_result_map[ipixel_]->retrieval_fit->original_pixel;
}

const pixel_t *grasp_output_tile_retrieval_fit_pixel_fit (const grasp_results_t *output, int it, int ix, int iy){
    return &output->tile_result_map[ipixel_]->retrieval_fit->fitting_pixel;
}

float grasp_output_tile_aerosol_opt_aexp (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_opt->Aexp;
}

float grasp_output_tile_aerosol_opt_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->aerosol_opt->extt[iwl];
}

float grasp_output_tile_aerosol_opt_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->aerosol_opt->ssat[iwl];
}

float grasp_output_tile_aerosol_opt_aextt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->aerosol_opt->aext[iwl];
}

float grasp_output_tile_aerosol_opt_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->aerosol_opt->ext[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_aerosol_opt_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->aerosol_opt->ssa[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_aerosol_opt_aext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->aerosol_opt->aext[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_aerosol_rind_mreal (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->aerosol_rind->mreal[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_aerosol_rind_mimag (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->aerosol_rind->mimag[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_aerosol_phmx_ph11 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->aerosol_phmx->ph11[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_ph12 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->aerosol_phmx->ph12[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_ph22 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->aerosol_phmx->ph22[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_ph33 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->aerosol_phmx->ph33[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_ph34 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->aerosol_phmx->ph34[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_ph44 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->aerosol_phmx->ph44[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_pht11 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->aerosol_phmx->pht11[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_pht12 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->aerosol_phmx->pht12[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_pht22 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->aerosol_phmx->pht22[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_pht33 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->aerosol_phmx->pht33[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_pht34 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->aerosol_phmx->pht34[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_aerosol_phmx_pht44 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->aerosol_phmx->pht44[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_aerosol_lidar_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->aerosol_lidar->lr[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_aerosol_lidar_ldpr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->aerosol_lidar->ldpr[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_aerosol_lidar_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->aerosol_lidar->lrt[iwl];
}

float grasp_output_tile_aerosol_lidar_ldprt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->aerosol_lidar->ldprt[iwl];
}

float grasp_output_tile_aerosol_sd2m_mph_cv (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->cv[i];
}

float grasp_output_tile_aerosol_sd2m_mph_std (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->std[i];
}

float grasp_output_tile_aerosol_sd2m_mph_rm (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->rm[i];
}

float grasp_output_tile_aerosol_sd2m_mph_reff (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->reff[i];
}

float grasp_output_tile_aerosol_sd2m_opt_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int i){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_ext->ext[index2D(iwl,i,2)];
}

float grasp_output_tile_aerosol_sd2m_mph_cv_fine_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->cv[0];
}

float grasp_output_tile_aerosol_sd2m_mph_std_fine_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->std[0];
}

float grasp_output_tile_aerosol_sd2m_mph_rm_fine_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->rm[0];
}

float grasp_output_tile_aerosol_sd2m_mph_reff_fine_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->reff[0];
}

float grasp_output_tile_aerosol_sd2m_opt_ext_fine_mode (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_ext->ext[index2D(iwl,0,2)];
}

float grasp_output_tile_aerosol_sd2m_mph_cv_coarse_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->cv[1];
}

float grasp_output_tile_aerosol_sd2m_mph_std_coarse_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->std[1];
}

float grasp_output_tile_aerosol_sd2m_mph_rm_coarse_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->rm[1];
}

float grasp_output_tile_aerosol_sd2m_mph_reff_coarse_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_mph->reff[1];
}

float grasp_output_tile_aerosol_sd2m_opt_ext_coarse_mode (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->aerosol_sd2m_ext->ext[index2D(iwl,1,2)];
}

float grasp_output_tile_aerosol_chem_rh (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->aerosol_chem->rh[isd];
}

float grasp_output_tile_aerosol_chem_fwtr (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->aerosol_chem->fwrt[isd];
}

float grasp_output_tile_aerosol_chem_fslbl (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->aerosol_chem->fslbl[isd];
}

float grasp_output_tile_aerosol_chem_finslbl (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->aerosol_chem->finslbl[isd];
}

float grasp_output_tile_aerosol_chem_fsoot (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->aerosol_chem->fsoot[isd];
}

float grasp_output_tile_aerosol_chem_firon (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->aerosol_chem->firon[isd];
}

/* modify by lei on 15/11/2016 */
float grasp_output_tile_aerosol_chem_fbrc (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->aerosol_chem->fbrc[isd];
}
/* modify by lei on 15/11/2016 */


float grasp_output_tile_aerosol_pm_pm (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->aerosol_pm->pm[i];
}

int grasp_output_tile_aerosol_types_index (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->aerosol_types->index;
}


/////////

float grasp_output_tile_clouds_opt_aexp (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_opt->Aexp;
}

float grasp_output_tile_clouds_opt_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->clouds_opt->extt[iwl];
}

float grasp_output_tile_clouds_opt_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->clouds_opt->ssat[iwl];
}

float grasp_output_tile_clouds_opt_aextt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->clouds_opt->aext[iwl];
}

float grasp_output_tile_clouds_opt_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->clouds_opt->ext[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_clouds_opt_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->clouds_opt->ssa[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_clouds_opt_aext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->clouds_opt->aext[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_clouds_rind_mreal (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->clouds_rind->mreal[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_clouds_rind_mimag (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->clouds_rind->mimag[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_clouds_phmx_ph11 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->clouds_phmx->ph11[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_ph12 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->clouds_phmx->ph12[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_ph22 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->clouds_phmx->ph22[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_ph33 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->clouds_phmx->ph33[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_ph34 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->clouds_phmx->ph34[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_ph44 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd){
    return output->tile_result_map[ipixel_]->clouds_phmx->ph44[index3D(iwl,isd,impar,NSD_,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_pht11 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->clouds_phmx->pht11[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_pht12 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->clouds_phmx->pht12[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_pht22 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->clouds_phmx->pht22[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_pht33 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->clouds_phmx->pht33[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_pht34 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->clouds_phmx->pht34[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_clouds_phmx_pht44 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar){
    return output->tile_result_map[ipixel_]->clouds_phmx->pht44[index2D(iwl,impar,NMPAR_)];
}

float grasp_output_tile_clouds_lidar_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->clouds_lidar->lr[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_clouds_lidar_ldpr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->clouds_lidar->ldpr[index2D(iwl,isd,NSD_)];
}

float grasp_output_tile_clouds_lidar_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->clouds_lidar->lrt[iwl];
}

float grasp_output_tile_clouds_lidar_ldprt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->clouds_lidar->ldprt[iwl];
}

float grasp_output_tile_clouds_sd2m_mph_cv (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->cv[i];
}

float grasp_output_tile_clouds_sd2m_mph_std (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->std[i];
}

float grasp_output_tile_clouds_sd2m_mph_rm (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->rm[i];
}

float grasp_output_tile_clouds_sd2m_mph_reff (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->reff[i];
}

float grasp_output_tile_clouds_sd2m_opt_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int i){
    return output->tile_result_map[ipixel_]->clouds_sd2m_ext->ext[index2D(iwl,i,2)];
}

float grasp_output_tile_clouds_sd2m_mph_cv_fine_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->cv[0];
}

float grasp_output_tile_clouds_sd2m_mph_std_fine_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->std[0];
}

float grasp_output_tile_clouds_sd2m_mph_rm_fine_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->rm[0];
}

float grasp_output_tile_clouds_sd2m_mph_reff_fine_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->reff[0];
}

float grasp_output_tile_clouds_sd2m_opt_ext_fine_mode (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->clouds_sd2m_ext->ext[index2D(iwl,0,2)];
}

float grasp_output_tile_clouds_sd2m_mph_cv_coarse_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->cv[1];
}

float grasp_output_tile_clouds_sd2m_mph_std_coarse_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->std[1];
}

float grasp_output_tile_clouds_sd2m_mph_rm_coarse_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->rm[1];
}

float grasp_output_tile_clouds_sd2m_mph_reff_coarse_mode (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_sd2m_mph->reff[1];
}

float grasp_output_tile_clouds_sd2m_opt_ext_coarse_mode (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->clouds_sd2m_ext->ext[index2D(iwl,1,2)];
}

float grasp_output_tile_clouds_chem_rh (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->clouds_chem->rh[isd];
}

float grasp_output_tile_clouds_chem_fwtr (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->clouds_chem->fwrt[isd];
}

float grasp_output_tile_clouds_chem_fslbl (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->clouds_chem->fslbl[isd];
}

float grasp_output_tile_clouds_chem_finslbl (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->clouds_chem->finslbl[isd];
}

float grasp_output_tile_clouds_chem_fsoot (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->clouds_chem->fsoot[isd];
}

float grasp_output_tile_clouds_chem_firon (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->clouds_chem->firon[isd];
}

float grasp_output_tile_clouds_chem_fbrc (const grasp_results_t *output, int it, int ix, int iy, int isd){
    return output->tile_result_map[ipixel_]->clouds_chem->fbrc[isd];
}


float grasp_output_tile_clouds_pm_pm (const grasp_results_t *output, int it, int ix, int iy, int i){
    return output->tile_result_map[ipixel_]->clouds_pm->pm[i];
}

int grasp_output_tile_clouds_types_index (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->clouds_types->index;
}

/////////

float grasp_output_tile_surface_ndvi (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->surface_surf->ndvi;
}

float grasp_output_tile_surface_salbedo (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->surface_surf->salbedo[iwl];
}


float grasp_output_tile_errest_par_errp (const grasp_results_t *output, int it, int ix, int iy, int ipar){
    return output->tile_result_map[ipixel_]->errest_par->ERRP[ipar];
}

float grasp_output_tile_errest_par_biasp (const grasp_results_t *output, int it, int ix, int iy, int ipar){
    return output->tile_result_map[ipixel_]->errest_par->BIASP[ipar];
}

float grasp_output_tile_errest_aerosol_opt_err_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_aerosol_opt->ERR_ext[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_aerosol_opt_bias_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_aerosol_opt->BIAS_ext[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_aerosol_opt_err_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_aerosol_opt->ERR_extt[iwl];
}

float grasp_output_tile_errest_aerosol_opt_bias_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_aerosol_opt->BIAS_extt[iwl];
}

float grasp_output_tile_errest_aerosol_opt_err_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_aerosol_opt->ERR_ssa[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_aerosol_opt_bias_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_aerosol_opt->BIAS_ssa[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_aerosol_opt_err_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_aerosol_opt->ERR_ssat[iwl];
}

float grasp_output_tile_errest_aerosol_opt_bias_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_aerosol_opt->BIAS_ssat[iwl];
}

float grasp_output_tile_errest_aerosol_lidar_err_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_aerosol_lidar->ERR_lr[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_aerosol_lidar_bias_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_aerosol_lidar->BIAS_lr[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_aerosol_lidar_err_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_aerosol_lidar->ERR_lrt[iwl];
}

float grasp_output_tile_errest_aerosol_lidar_bias_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_aerosol_lidar->BIAS_lrt[iwl];
}
/////////

float grasp_output_tile_errest_clouds_opt_err_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_clouds_opt->ERR_ext[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_clouds_opt_bias_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_clouds_opt->BIAS_ext[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_clouds_opt_err_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_clouds_opt->ERR_extt[iwl];
}

float grasp_output_tile_errest_clouds_opt_bias_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_clouds_opt->BIAS_extt[iwl];
}

float grasp_output_tile_errest_clouds_opt_err_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_clouds_opt->ERR_ssa[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_clouds_opt_bias_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_clouds_opt->BIAS_ssa[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_clouds_opt_err_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_clouds_opt->ERR_ssat[iwl];
}

float grasp_output_tile_errest_clouds_opt_bias_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_clouds_opt->BIAS_ssat[iwl];
}

float grasp_output_tile_errest_clouds_lidar_err_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_clouds_lidar->ERR_lr[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_clouds_lidar_bias_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd){
    return output->tile_result_map[ipixel_]->errest_clouds_lidar->BIAS_lr[index2D(isd,iwl,NWL_)];
}

float grasp_output_tile_errest_clouds_lidar_err_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_clouds_lidar->ERR_lrt[iwl];
}

float grasp_output_tile_errest_clouds_lidar_bias_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl){
    return output->tile_result_map[ipixel_]->errest_clouds_lidar->BIAS_lrt[iwl];
}
///////////////

int grasp_output_tile_forcing_bbflux_nhlv (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->forcing_bbflux->nhlv;
}

float grasp_output_tile_forcing_bbflux_bbufx0 (const grasp_results_t *output, int it, int ix, int iy, int iknt){
    return output->tile_result_map[ipixel_]->forcing_bbflux->bbufx0[iknt];
}

float grasp_output_tile_forcing_bbflux_bbdfx0 (const grasp_results_t *output, int it, int ix, int iy, int iknt){
    return output->tile_result_map[ipixel_]->forcing_bbflux->bbdfx0[iknt];
}

float grasp_output_tile_forcing_bbflux_bbufxa (const grasp_results_t *output, int it, int ix, int iy, int iknt){
    return output->tile_result_map[ipixel_]->forcing_bbflux->bbufxa[iknt];
}

float grasp_output_tile_forcing_bbflux_bbdfxa (const grasp_results_t *output, int it, int ix, int iy, int iknt){
    return output->tile_result_map[ipixel_]->forcing_bbflux->bbdfxa[iknt];
}

float grasp_output_tile_forcing_bbflux_hlv (const grasp_results_t *output, int it, int ix, int iy, int iknt){
    return output->tile_result_map[ipixel_]->forcing_bbflux->hlv[iknt];
}

int grasp_output_tile_forcing_forcing_nhlv (const grasp_results_t *output, int it, int ix, int iy){
    return output->tile_result_map[ipixel_]->forcing_forcing->nhlv;
}

float grasp_output_tile_forcing_forcing_netforc (const grasp_results_t *output, int it, int ix, int iy, int iknt){
    return output->tile_result_map[ipixel_]->forcing_forcing->netforc[iknt];
}

float grasp_output_tile_forcing_forcing_forceff (const grasp_results_t *output, int it, int ix, int iy, int iknt){
    return output->tile_result_map[ipixel_]->forcing_forcing->forceff[iknt];
}

float grasp_output_tile_forcing_forcing_hlv (const grasp_results_t *output, int it, int ix, int iy, int iknt){
    return output->tile_result_map[ipixel_]->forcing_forcing->hlv[iknt];
}
