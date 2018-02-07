/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <grasp/utils.h>
#include "grasp_output_segment_result.h"


/// NOTE: A developer can take advantage of the this private access to fields
///       to implement assertion verifications (values inside valid ranges of arrays)

const float *grasp_output_segment_parameters(const output_segment_general *output, int ipix){
    return output->retrieval.par.pixel[ipix].par;
}

bool grasp_output_segment_products_retrieval_res (const output_segment_general *output) {
    return output->products.retrieval.res;
}

bool grasp_output_segment_products_retrieval_par (const output_segment_general *output) {
    return output->products.retrieval.par;
}

bool grasp_output_segment_products_retrieval_fit (const output_segment_general *output) {
    return output->products.retrieval.fit;
}

bool grasp_output_segment_products_aerosol_opt (const output_segment_general *output) {
    return output->products.aerosol.opt;
}

bool grasp_output_segment_products_aerosol_rind (const output_segment_general *output) {
    return output->products.aerosol.rind;
}

bool grasp_output_segment_products_aerosol_chem (const output_segment_general *output) {
    return output->products.aerosol.chem;
}

bool grasp_output_segment_products_aerosol_phmx (const output_segment_general *output) {
    return output->products.aerosol.phmx;
}

bool grasp_output_segment_products_aerosol_lidar (const output_segment_general *output) {
    return output->products.aerosol.lidar;
}

bool grasp_output_segment_products_aerosol_sd2m_mph (const output_segment_general *output) {
    return output->products.aerosol.sd2m_mph;
}

bool grasp_output_segment_products_aerosol_sd2m_ext (const output_segment_general *output) {
    return output->products.aerosol.sd2m_ext;
}

bool grasp_output_segment_products_aerosol_pm (const output_segment_general *output) {
    return output->products.aerosol.pm;
}

bool grasp_output_segment_products_aerosol_types (const output_segment_general *output) {
    return output->products.aerosol.types;
}

bool grasp_output_segment_products_clouds_opt (const output_segment_general *output) {
    return output->products.clouds.opt;
}

bool grasp_output_segment_products_clouds_rind (const output_segment_general *output) {
    return output->products.clouds.rind;
}

bool grasp_output_segment_products_clouds_chem (const output_segment_general *output) {
    return output->products.clouds.chem;
}

bool grasp_output_segment_products_clouds_phmx (const output_segment_general *output) {
    return output->products.clouds.phmx;
}

bool grasp_output_segment_products_clouds_lidar (const output_segment_general *output) {
    return output->products.clouds.lidar;
}

bool grasp_output_segment_products_clouds_sd2m_mph (const output_segment_general *output) {
    return output->products.clouds.sd2m_mph;
}

bool grasp_output_segment_products_clouds_sd2m_ext (const output_segment_general *output) {
    return output->products.clouds.sd2m_ext;
}

bool grasp_output_segment_products_clouds_pm (const output_segment_general *output) {
    return output->products.clouds.pm;
}

bool grasp_output_segment_products_clouds_types (const output_segment_general *output) {
    return output->products.clouds.types;
}

bool grasp_output_segment_products_surface_surf (const output_segment_general *output) {
    return output->products.surface.surf;
}

bool grasp_output_segment_products_errest_par (const output_segment_general *output) {
    return output->products.errest.par;
}

bool grasp_output_segment_products_errest_aerosol_opt (const output_segment_general *output) {
    return output->products.errest.aerosol.opt;
}

bool grasp_output_segment_products_errest_aerosol_lidar (const output_segment_general *output) {
    return output->products.errest.aerosol.lidar;
}

bool grasp_output_segment_products_errest_clouds_opt (const output_segment_general *output) {
    return output->products.errest.clouds.opt;
}

bool grasp_output_segment_products_errest_clouds_lidar (const output_segment_general *output) {
    return output->products.errest.clouds.lidar;
}

bool grasp_output_segment_products_forcing_bbflux (const output_segment_general *output) {
    return output->products.forcing.bbflux;
}

bool grasp_output_segment_products_forcing_forcing (const output_segment_general *output) {
    return output->products.forcing.forcing;
}

///////////////////

int grasp_output_segment_retrieval_res_niter (const output_segment_general *output) {
    return output->retrieval.res.niter;
}

float grasp_output_segment_retrieval_res_rest (const output_segment_general *output) {
    return output->retrieval.res.rest;
}

float grasp_output_segment_retrieval_res_resat (const output_segment_general *output, int inoise) {
    return output->retrieval.res.resat[inoise];
}

float grasp_output_segment_retrieval_res_resrt (const output_segment_general *output, int inoise) {
    return output->retrieval.res.resrt[inoise];
}

int grasp_output_segment_retrieval_res_pixel_niter (const output_segment_general *output, int ipix) {
    return output->retrieval.res.pixel[ipix].niter;
}

float grasp_output_segment_retrieval_res_pixel_res (const output_segment_general *output, int ipix) {
    return output->retrieval.res.pixel[ipix].res;
}

float grasp_output_segment_retrieval_res_pixel_resa (const output_segment_general *output, int ipix, int inoise) {
    return output->retrieval.res.pixel[ipix].resa[inoise];
}

float grasp_output_segment_retrieval_res_pixel_resr (const output_segment_general *output, int ipix, int inoise) {
    return output->retrieval.res.pixel[ipix].resr[inoise];
}

int grasp_output_segment_retrieval_par_ngrid (const output_segment_general *output) {
    return output->retrieval.par.ngrid;
}

float grasp_output_segment_retrieval_par_radius (const output_segment_general *output, int irr) {
    return output->retrieval.par.radius[irr];
}

float grasp_output_segment_retrieval_par_sdl (const output_segment_general *output, int irr, int irc) {
    return output->retrieval.par.SDL[irc][irr];
}

float grasp_output_segment_retrieval_par_parameters (const output_segment_general *output, int ipix, int ipar) {
    return output->retrieval.par.pixel[ipix].par[ipar];
}

const sensor_data_t *grasp_output_segment_retrieval_fit_segment_fit (const output_segment_general *output) {
    return &output->retrieval.fit.segment_fit;
}

///////////////

float grasp_output_segment_aerosol_opt_aexp (const output_segment_general *output, int ipix) {
    return output->aerosol.opt.pixel[ipix].Aexp;
}

float grasp_output_segment_aerosol_opt_extt (const output_segment_general *output, int ipix, int iwl) {
    return output->aerosol.opt.pixel[ipix].wl[iwl].extt;
}

float grasp_output_segment_aerosol_opt_ssat (const output_segment_general *output, int ipix, int iwl) {
    return output->aerosol.opt.pixel[ipix].wl[iwl].ssat;
}

float grasp_output_segment_aerosol_opt_aextt (const output_segment_general *output, int ipix, int iwl) {
    return output->aerosol.opt.pixel[ipix].wl[iwl].aextt;
}

float grasp_output_segment_aerosol_opt_ext (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->aerosol.opt.pixel[ipix].wl[iwl].ext[isd];
}

float grasp_output_segment_aerosol_opt_ssa (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->aerosol.opt.pixel[ipix].wl[iwl].ssa[isd];
}

float grasp_output_segment_aerosol_opt_aext (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->aerosol.opt.pixel[ipix].wl[iwl].aext[isd];
}

float grasp_output_segment_aerosol_rind_mreal (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->aerosol.rind.pixel[ipix].wl[iwl].mreal[isd];
}

float grasp_output_segment_aerosol_rind_mimag (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->aerosol.rind.pixel[ipix].wl[iwl].mimag[isd];
}

int grasp_output_segment_aerosol_phmx_nangle (const output_segment_general *output){
    return output->aerosol.phmx.nangle;
}

float grasp_output_segment_aerosol_phmx_angle (const output_segment_general *output, int impar){
    return output->aerosol.phmx.angle[impar];
}

float grasp_output_segment_aerosol_phmx_ph11 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].ph11[isd][impar];
}

float grasp_output_segment_aerosol_phmx_ph12 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].ph12[isd][impar];
}

float grasp_output_segment_aerosol_phmx_ph22 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].ph22[isd][impar];
}

float grasp_output_segment_aerosol_phmx_ph33 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].ph33[isd][impar];
}

float grasp_output_segment_aerosol_phmx_ph34 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].ph34[isd][impar];
}

float grasp_output_segment_aerosol_phmx_ph44 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].ph44[isd][impar];
}

float grasp_output_segment_aerosol_phmx_pht11 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].pht11[impar];
}

float grasp_output_segment_aerosol_phmx_pht12 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].pht12[impar];
}

float grasp_output_segment_aerosol_phmx_pht22 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].pht22[impar];
}

float grasp_output_segment_aerosol_phmx_pht33 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].pht33[impar];
}

float grasp_output_segment_aerosol_phmx_pht34 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].pht34[impar];
}

float grasp_output_segment_aerosol_phmx_pht44 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->aerosol.phmx.pixel[ipix].wl[iwl].pht44[impar];
}

float grasp_output_segment_aerosol_lidar_lr (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->aerosol.lidar.pixel[ipix].wl[iwl].lr[isd];
}

float grasp_output_segment_aerosol_lidar_ldpr (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->aerosol.lidar.pixel[ipix].wl[iwl].ldpr[isd];
}

float grasp_output_segment_aerosol_lidar_lrt (const output_segment_general *output, int ipix, int iwl){
    return output->aerosol.lidar.pixel[ipix].wl[iwl].lrt;
}

float grasp_output_segment_aerosol_lidar_ldprt (const output_segment_general *output, int ipix, int iwl){
    return output->aerosol.lidar.pixel[ipix].wl[iwl].ldprt;
}

float grasp_output_segment_aerosol_sd2m_mph_cv (const output_segment_general *output, int ipix, int i){
    return output->aerosol.sd2m.mph.pixel[ipix].cv[i];
}

float grasp_output_segment_aerosol_sd2m_mph_std (const output_segment_general *output, int ipix, int i){
    return output->aerosol.sd2m.mph.pixel[ipix].std[i];
}

float grasp_output_segment_aerosol_sd2m_mph_rm (const output_segment_general *output, int ipix, int i){
    return output->aerosol.sd2m.mph.pixel[ipix].rm[i];
}

float grasp_output_segment_aerosol_sd2m_mph_reff (const output_segment_general *output, int ipix, int i){
    return output->aerosol.sd2m.mph.pixel[ipix].reff[i];
}

float grasp_output_segment_aerosol_sd2m_opt_ext (const output_segment_general *output, int ipix, int iwl, int i){
    return output->aerosol.sd2m.opt.pixel[ipix].ext[i][iwl];
}

float grasp_output_segment_aerosol_sd2m_mph_cv_fine_mode (const output_segment_general *output, int ipix){
    return output->aerosol.sd2m.mph.pixel[ipix].cv[0];
}

float grasp_output_segment_aerosol_sd2m_mph_std_fine_mode (const output_segment_general *output, int ipix){
    return output->aerosol.sd2m.mph.pixel[ipix].std[0];
}

float grasp_output_segment_aerosol_sd2m_mph_rm_fine_mode (const output_segment_general *output, int ipix){
    return output->aerosol.sd2m.mph.pixel[ipix].rm[0];
}

float grasp_output_segment_aerosol_sd2m_mph_reff_fine_mode (const output_segment_general *output, int ipix){
    return output->aerosol.sd2m.mph.pixel[ipix].reff[0];
}

float grasp_output_segment_aerosol_sd2m_opt_ext_fine_mode (const output_segment_general *output, int ipix, int iwl){
    return output->aerosol.sd2m.opt.pixel[ipix].ext[0][iwl];
}

float grasp_output_segment_aerosol_sd2m_mph_cv_coarse_mode (const output_segment_general *output, int ipix){
    return output->aerosol.sd2m.mph.pixel[ipix].cv[1];
}

float grasp_output_segment_aerosol_sd2m_mph_std_coarse_mode (const output_segment_general *output, int ipix){
    return output->aerosol.sd2m.mph.pixel[ipix].std[1];
}

float grasp_output_segment_aerosol_sd2m_mph_rm_coarse_mode (const output_segment_general *output, int ipix){
    return output->aerosol.sd2m.mph.pixel[ipix].rm[1];
}

float grasp_output_segment_aerosol_sd2m_mph_reff_coarse_mode (const output_segment_general *output, int ipix){
    return output->aerosol.sd2m.mph.pixel[ipix].reff[1];
}

float grasp_output_segment_aerosol_sd2m_opt_ext_coarse_mode (const output_segment_general *output, int ipix, int iwl){
    return output->aerosol.sd2m.opt.pixel[ipix].ext[1][iwl];
}

float grasp_output_segment_aerosol_chem_rh (const output_segment_general *output, int ipix, int isd){
    return output->aerosol.chem.pixel[ipix].rh[isd];
}

float grasp_output_segment_aerosol_chem_fwtr (const output_segment_general *output, int ipix, int isd){
    return output->aerosol.chem.pixel[ipix].fwtr[isd];
}

float grasp_output_segment_aerosol_chem_fslbl (const output_segment_general *output, int ipix, int isd){
    return output->aerosol.chem.pixel[ipix].fslbl[isd];
}

float grasp_output_segment_aerosol_chem_finslbl (const output_segment_general *output, int ipix, int isd){
    return output->aerosol.chem.pixel[ipix].finslbl[isd];
}

float grasp_output_segment_aerosol_chem_fsoot (const output_segment_general *output, int ipix, int isd){
    return output->aerosol.chem.pixel[ipix].fsoot[isd];
}

float grasp_output_segment_aerosol_chem_firon (const output_segment_general *output, int ipix, int isd){
    return output->aerosol.chem.pixel[ipix].firon[isd];
}

/* add by lei on 15/11/2016 */
float grasp_output_segment_aerosol_chem_fbrc (const output_segment_general *output, int ipix, int isd){
    return output->aerosol.chem.pixel[ipix].fbrc[isd];
}


float grasp_output_segment_aerosol_pm_pm (const output_segment_general *output, int ipix, int i){
    return output->aerosol.pm.pixel[ipix].pm[i];
}

int grasp_output_segment_aerosol_types_index (const output_segment_general *output, int ipix){
    return output->aerosol.types.pixel[ipix].index;
}

///////////////

float grasp_output_segment_clouds_opt_aexp (const output_segment_general *output, int ipix) {
    return output->clouds.opt.pixel[ipix].Aexp;
}

float grasp_output_segment_clouds_opt_extt (const output_segment_general *output, int ipix, int iwl) {
    return output->clouds.opt.pixel[ipix].wl[iwl].extt;
}

float grasp_output_segment_clouds_opt_ssat (const output_segment_general *output, int ipix, int iwl) {
    return output->clouds.opt.pixel[ipix].wl[iwl].ssat;
}

float grasp_output_segment_clouds_opt_aextt (const output_segment_general *output, int ipix, int iwl) {
    return output->clouds.opt.pixel[ipix].wl[iwl].aextt;
}

float grasp_output_segment_clouds_opt_ext (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->clouds.opt.pixel[ipix].wl[iwl].ext[isd];
}

float grasp_output_segment_clouds_opt_ssa (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->clouds.opt.pixel[ipix].wl[iwl].ssa[isd];
}

float grasp_output_segment_clouds_opt_aext (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->clouds.opt.pixel[ipix].wl[iwl].aext[isd];
}

float grasp_output_segment_clouds_rind_mreal (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->clouds.rind.pixel[ipix].wl[iwl].mreal[isd];
}

float grasp_output_segment_clouds_rind_mimag (const output_segment_general *output, int ipix, int iwl, int isd) {
    return output->clouds.rind.pixel[ipix].wl[iwl].mimag[isd];
}

int grasp_output_segment_clouds_phmx_nangle (const output_segment_general *output){
    return output->clouds.phmx.nangle;
}

float grasp_output_segment_clouds_phmx_angle (const output_segment_general *output, int impar){
    return output->clouds.phmx.angle[impar];
}

float grasp_output_segment_clouds_phmx_ph11 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->clouds.phmx.pixel[ipix].wl[iwl].ph11[isd][impar];
}

float grasp_output_segment_clouds_phmx_ph12 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->clouds.phmx.pixel[ipix].wl[iwl].ph12[isd][impar];
}

float grasp_output_segment_clouds_phmx_ph22 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->clouds.phmx.pixel[ipix].wl[iwl].ph22[isd][impar];
}

float grasp_output_segment_clouds_phmx_ph33 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->clouds.phmx.pixel[ipix].wl[iwl].ph33[isd][impar];
}

float grasp_output_segment_clouds_phmx_ph34 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->clouds.phmx.pixel[ipix].wl[iwl].ph34[isd][impar];
}

float grasp_output_segment_clouds_phmx_ph44 (const output_segment_general *output, int ipix, int iwl, int impar, int isd){
    return output->clouds.phmx.pixel[ipix].wl[iwl].ph44[isd][impar];
}

float grasp_output_segment_clouds_phmx_pht11 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->clouds.phmx.pixel[ipix].wl[iwl].pht11[impar];
}

float grasp_output_segment_clouds_phmx_pht12 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->clouds.phmx.pixel[ipix].wl[iwl].pht12[impar];
}

float grasp_output_segment_clouds_phmx_pht22 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->clouds.phmx.pixel[ipix].wl[iwl].pht22[impar];
}

float grasp_output_segment_clouds_phmx_pht33 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->clouds.phmx.pixel[ipix].wl[iwl].pht33[impar];
}

float grasp_output_segment_clouds_phmx_pht34 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->clouds.phmx.pixel[ipix].wl[iwl].pht34[impar];
}

float grasp_output_segment_clouds_phmx_pht44 (const output_segment_general *output, int ipix, int iwl, int impar){
    return output->clouds.phmx.pixel[ipix].wl[iwl].pht44[impar];
}

float grasp_output_segment_clouds_lidar_lr (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->clouds.lidar.pixel[ipix].wl[iwl].lr[isd];
}

float grasp_output_segment_clouds_lidar_ldpr (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->clouds.lidar.pixel[ipix].wl[iwl].ldpr[isd];
}

float grasp_output_segment_clouds_lidar_lrt (const output_segment_general *output, int ipix, int iwl){
    return output->clouds.lidar.pixel[ipix].wl[iwl].lrt;
}

float grasp_output_segment_clouds_lidar_ldprt (const output_segment_general *output, int ipix, int iwl){
    return output->clouds.lidar.pixel[ipix].wl[iwl].ldprt;
}

float grasp_output_segment_clouds_sd2m_mph_cv (const output_segment_general *output, int ipix, int i){
    return output->clouds.sd2m.mph.pixel[ipix].cv[i];
}

float grasp_output_segment_clouds_sd2m_mph_std (const output_segment_general *output, int ipix, int i){
    return output->clouds.sd2m.mph.pixel[ipix].std[i];
}

float grasp_output_segment_clouds_sd2m_mph_rm (const output_segment_general *output, int ipix, int i){
    return output->clouds.sd2m.mph.pixel[ipix].rm[i];
}

float grasp_output_segment_clouds_sd2m_mph_reff (const output_segment_general *output, int ipix, int i){
    return output->clouds.sd2m.mph.pixel[ipix].reff[i];
}

float grasp_output_segment_clouds_sd2m_opt_ext (const output_segment_general *output, int ipix, int iwl, int i){
    return output->clouds.sd2m.opt.pixel[ipix].ext[i][iwl];
}

float grasp_output_segment_clouds_sd2m_mph_cv_fine_mode (const output_segment_general *output, int ipix){
    return output->clouds.sd2m.mph.pixel[ipix].cv[0];
}

float grasp_output_segment_clouds_sd2m_mph_std_fine_mode (const output_segment_general *output, int ipix){
    return output->clouds.sd2m.mph.pixel[ipix].std[0];
}

float grasp_output_segment_clouds_sd2m_mph_rm_fine_mode (const output_segment_general *output, int ipix){
    return output->clouds.sd2m.mph.pixel[ipix].rm[0];
}

float grasp_output_segment_clouds_sd2m_mph_reff_fine_mode (const output_segment_general *output, int ipix){
    return output->clouds.sd2m.mph.pixel[ipix].reff[0];
}

float grasp_output_segment_clouds_sd2m_opt_ext_fine_mode (const output_segment_general *output, int ipix, int iwl){
    return output->clouds.sd2m.opt.pixel[ipix].ext[0][iwl];
}

float grasp_output_segment_clouds_sd2m_mph_cv_coarse_mode (const output_segment_general *output, int ipix){
    return output->clouds.sd2m.mph.pixel[ipix].cv[1];
}

float grasp_output_segment_clouds_sd2m_mph_std_coarse_mode (const output_segment_general *output, int ipix){
    return output->clouds.sd2m.mph.pixel[ipix].std[1];
}

float grasp_output_segment_clouds_sd2m_mph_rm_coarse_mode (const output_segment_general *output, int ipix){
    return output->clouds.sd2m.mph.pixel[ipix].rm[1];
}

float grasp_output_segment_clouds_sd2m_mph_reff_coarse_mode (const output_segment_general *output, int ipix){
    return output->clouds.sd2m.mph.pixel[ipix].reff[1];
}

float grasp_output_segment_clouds_sd2m_opt_ext_coarse_mode (const output_segment_general *output, int ipix, int iwl){
    return output->clouds.sd2m.opt.pixel[ipix].ext[1][iwl];
}

float grasp_output_segment_clouds_chem_rh (const output_segment_general *output, int ipix, int isd){
    return output->clouds.chem.pixel[ipix].rh[isd];
}

float grasp_output_segment_clouds_chem_fwtr (const output_segment_general *output, int ipix, int isd){
    return output->clouds.chem.pixel[ipix].fwtr[isd];
}

float grasp_output_segment_clouds_chem_fslbl (const output_segment_general *output, int ipix, int isd){
    return output->clouds.chem.pixel[ipix].fslbl[isd];
}

float grasp_output_segment_clouds_chem_finslbl (const output_segment_general *output, int ipix, int isd){
    return output->clouds.chem.pixel[ipix].finslbl[isd];
}

float grasp_output_segment_clouds_chem_fsoot (const output_segment_general *output, int ipix, int isd){
    return output->clouds.chem.pixel[ipix].fsoot[isd];
}

float grasp_output_segment_clouds_chem_firon (const output_segment_general *output, int ipix, int isd){
    return output->clouds.chem.pixel[ipix].firon[isd];
}


/* add by lei on 15/11/2016 */
float grasp_output_segment_clouds_chem_fbrc (const output_segment_general *output, int ipix, int isd){
    return output->clouds.chem.pixel[ipix].fbrc[isd];
}


float grasp_output_segment_clouds_pm_pm (const output_segment_general *output, int ipix, int i){
    return output->clouds.pm.pixel[ipix].pm[i];
}

int grasp_output_segment_clouds_types_index (const output_segment_general *output, int ipix){
    return output->clouds.types.pixel[ipix].index;
}


///////////////

float grasp_output_segment_surface_ndvi (const output_segment_general *output, int ipix){
    return output->surface.pixel[ipix].ndvi;
}

float grasp_output_segment_surface_salbedo (const output_segment_general *output, int ipix, int iwl){
    return output->surface.pixel[ipix].wl[iwl].salbedo;
}

/////////////////

float grasp_output_segment_errest_par_errp (const output_segment_general *output, int ipix, int ipar){
    return output->errest.par.pixel[ipix].ERRP[ipar];
}

float grasp_output_segment_errest_par_biasp (const output_segment_general *output, int ipix, int ipar){
    return output->errest.par.pixel[ipix].BIASP[ipar];
}

float grasp_output_segment_errest_aerosol_opt_err_ext (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.aerosol.opt.pixel[ipix].ERR_ext[isd][iwl];
}

float grasp_output_segment_errest_aerosol_opt_bias_ext (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.aerosol.opt.pixel[ipix].BIAS_ext[isd][iwl];
}

float grasp_output_segment_errest_aerosol_opt_err_extt (const output_segment_general *output, int ipix, int iwl){
    return output->errest.aerosol.opt.pixel[ipix].ERR_extt[iwl];
}

float grasp_output_segment_errest_aerosol_opt_bias_extt (const output_segment_general *output, int ipix, int iwl){
    return output->errest.aerosol.opt.pixel[ipix].BIAS_extt[iwl];
}

float grasp_output_segment_errest_aerosol_opt_err_ssa (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.aerosol.opt.pixel[ipix].ERR_ssa[isd][iwl];
}

float grasp_output_segment_errest_aerosol_opt_bias_ssa (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.aerosol.opt.pixel[ipix].BIAS_ssa[isd][iwl];
}

float grasp_output_segment_errest_aerosol_opt_err_ssat (const output_segment_general *output, int ipix, int iwl){
    return output->errest.aerosol.opt.pixel[ipix].ERR_ssat[iwl];
}

float grasp_output_segment_errest_aerosol_opt_bias_ssat (const output_segment_general *output, int ipix, int iwl){
    return output->errest.aerosol.opt.pixel[ipix].BIAS_ssat[iwl];
}

float grasp_output_segment_errest_aerosol_lidar_err_lr (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.aerosol.lidar.pixel[ipix].ERR_lr[isd][iwl];
}

float grasp_output_segment_errest_aerosol_lidar_bias_lr (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.aerosol.lidar.pixel[ipix].BIAS_lr[isd][iwl];
}

float grasp_output_segment_errest_aerosol_lidar_err_lrt (const output_segment_general *output, int ipix, int iwl){
    return output->errest.aerosol.lidar.pixel[ipix].ERR_lrt[iwl];
}

float grasp_output_segment_errest_aerosol_lidar_bias_lrt (const output_segment_general *output, int ipix, int iwl){
    return output->errest.aerosol.lidar.pixel[ipix].BIAS_lrt[iwl];
}

float grasp_output_segment_errest_clouds_opt_err_ext (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.clouds.opt.pixel[ipix].ERR_ext[isd][iwl];
}

float grasp_output_segment_errest_clouds_opt_bias_ext (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.clouds.opt.pixel[ipix].BIAS_ext[isd][iwl];
}

float grasp_output_segment_errest_clouds_opt_err_extt (const output_segment_general *output, int ipix, int iwl){
    return output->errest.clouds.opt.pixel[ipix].ERR_extt[iwl];
}

float grasp_output_segment_errest_clouds_opt_bias_extt (const output_segment_general *output, int ipix, int iwl){
    return output->errest.clouds.opt.pixel[ipix].BIAS_extt[iwl];
}

float grasp_output_segment_errest_clouds_opt_err_ssa (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.clouds.opt.pixel[ipix].ERR_ssa[isd][iwl];
}

float grasp_output_segment_errest_clouds_opt_bias_ssa (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.clouds.opt.pixel[ipix].BIAS_ssa[isd][iwl];
}

float grasp_output_segment_errest_clouds_opt_err_ssat (const output_segment_general *output, int ipix, int iwl){
    return output->errest.clouds.opt.pixel[ipix].ERR_ssat[iwl];
}

float grasp_output_segment_errest_clouds_opt_bias_ssat (const output_segment_general *output, int ipix, int iwl){
    return output->errest.clouds.opt.pixel[ipix].BIAS_ssat[iwl];
}

//
float grasp_output_segment_errest_clouds_lidar_err_lr (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.clouds.lidar.pixel[ipix].ERR_lr[isd][iwl];
}

float grasp_output_segment_errest_clouds_lidar_bias_lr (const output_segment_general *output, int ipix, int iwl, int isd){
    return output->errest.clouds.lidar.pixel[ipix].BIAS_lr[isd][iwl];
}

float grasp_output_segment_errest_clouds_lidar_err_lrt (const output_segment_general *output, int ipix, int iwl){
    return output->errest.clouds.lidar.pixel[ipix].ERR_lrt[iwl];
}

float grasp_output_segment_errest_clouds_lidar_bias_lrt (const output_segment_general *output, int ipix, int iwl){
    return output->errest.clouds.lidar.pixel[ipix].BIAS_lrt[iwl];
}

///////////////

int grasp_output_segment_forcing_bbflux_nhlv (const output_segment_general *output, int ipix){
    return output->forcing.bbflux.pixel[ipix].nhlv;
}

float grasp_output_segment_forcing_bbflux_bbufx0 (const output_segment_general *output, int ipix, int iknt){
    return output->forcing.bbflux.pixel[ipix].bbufx0[iknt];
}

float grasp_output_segment_forcing_bbflux_bbdfx0 (const output_segment_general *output, int ipix, int iknt){
    return output->forcing.bbflux.pixel[ipix].bbdfx0[iknt];
}

float grasp_output_segment_forcing_bbflux_bbufxa (const output_segment_general *output, int ipix, int iknt){
    return output->forcing.bbflux.pixel[ipix].bbufxa[iknt];
}

float grasp_output_segment_forcing_bbflux_bbdfxa (const output_segment_general *output, int ipix, int iknt){
    return output->forcing.bbflux.pixel[ipix].bbdfxa[iknt];
}

float grasp_output_segment_forcing_bbflux_hlv (const output_segment_general *output, int ipix, int iknt){
    return output->forcing.bbflux.pixel[ipix].hlv[iknt];
}

int grasp_output_segment_forcing_forcing_nhlv (const output_segment_general *output, int ipix){
    return output->forcing.forcing.pixel[ipix].nhlv;
}

float grasp_output_segment_forcing_forcing_netforc (const output_segment_general *output, int ipix, int iknt){
    return output->forcing.forcing.pixel[ipix].netforc[iknt];
}

float grasp_output_segment_forcing_forcing_forceff (const output_segment_general *output, int ipix, int iknt){
    return output->forcing.forcing.pixel[ipix].forceff[iknt];
}

float grasp_output_segment_forcing_forcing_hlv (const output_segment_general *output, int ipix, int iknt){
    return output->forcing.forcing.pixel[ipix].hlv[iknt];
}






