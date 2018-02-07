/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

 /**
 * @file grasp_output_stream.h
 * @author David Fuertes
 * @date 11 Jun 2014
 * @brief This library defines a general way to work with the output
 *
 * Technical documentation available: (http://www.grasp-open.com/tech-doc/chap02.php#chap02030101)
 * User documentation: (http://www.grasp-open.com/doc/ch04.php#settings-file)
 */

#ifndef GRASP_STREAM_H
#define	GRASP_STREAM_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "../settings/grasp_settings.h"
#include "../input/grasp_input.h"
#include "../output/grasp_output.h"    



    /**
     * @brief Return a allocated string created from original_name replacing the special pattern wildcards. 
     * 
     * Return a allocated string created from original_name replacing the special pattern wildcards. 
     * In addition, gs structure will be change to contains the information about if the stream is writable. 
     * The following list show all available wildcards that can be used for creating dynamic output filenames:
     * - auto(N): #grasp_output_stream_filename_generate_by_auto
     * - icol(N): #grasp_output_stream_filename_generate_by_icol
     * - irow(N): #grasp_output_stream_filename_generate_by_irow
     * - itime(N): #grasp_output_stream_filename_generate_by_itime
     * - segment_nx(N): #grasp_output_stream_filename_generate_by_segment_nx
     * - segment_ny(N): #grasp_output_stream_filename_generate_by_segment_ny
     * - segment_nt(N): #grasp_output_stream_filename_generate_by_segment_nt
     * - tile_from(FORMAT): #grasp_output_stream_filename_generate_by_tile_from
     * - tile_to(FORMAT): #grasp_output_stream_filename_generate_by_tile_to
     * - tile_corner_column(N): #grasp_output_stream_filename_generate_by_tile_corner_column
     * - tile_corner_row(N): #grasp_output_stream_filename_generate_by_tile_corner_row
     * - tile_center_longitude(FORMAT): #grasp_output_stream_filename_generate_by_tile_center_longitude
     * - tile_center_latitude(FORMAT): #grasp_output_stream_filename_generate_by_tile_center_latitude
     * - tile_coordinate_x(I): #grasp_output_stream_filename_generate_by_tile_coordinate_x
     * - tile_coordinate_y(I): #grasp_output_stream_filename_generate_by_tile_coordinate_y
     * - tile_width(N): #grasp_output_stream_filename_generate_by_tile_width
     * - tile_height(N): #grasp_output_stream_filename_generate_by_tile_height
     * - segment_corner_column(N): #grasp_output_stream_filename_generate_by_segment_corner_column
     * - segment_corner_row(N): #grasp_output_stream_filename_generate_by_segment_corner_row
     * - settings_filename: #grasp_output_stream_filename_generate_by_settings_filename
     * - version: #grasp_output_stream_filename_generate_by_version
     * - branch: #grasp_output_stream_filename_generate_by_branch
     * - commit: #grasp_output_stream_filename_generate_by_commit
     * - constants_set: #grasp_output_stream_filename_generate_by_constants_set
     * - pwd: #grasp_output_stream_filename_generate_by_pwd
     * - yml: #grasp_output_stream_filename_generate_by_yml
     * @param gs definition of grasp output stream 
     * @param settings current settings used
     * @param segment current processed segment
     * @param output output of processed segment
     * @param tile_description general description of tile
     * @param icol column position of the segment in the tile
     * @param irow row position of the segment in the tile
     * @param itime time position of the segment in the tile
     * @return Return a allocated string created from original_name replacing the special pattern wildcards
     */
    char *grasp_output_stream_filename(grasp_output_stream *gs, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief Validate if filename is well-formed
     * @param original_name Pattern to check if it is a valid stream pattern
     * @return true if the patern is a valid one
     */
    bool grasp_output_stream_filename_validation(const char *original_name);


    /**
     * @brief Initilize gs stream argument with filename string.
     * @param[in] filename a valid grasp output stream pattern
     * @param[out] gs grasp output stream structure initialized
     */
    void grasp_output_stream_initialize(const char *filename, grasp_output_stream *gs);

    /**
     * @brief Open a stream returning FILE pointer.
     * Open a stream returning FILE pointer. You should use grasp_stream_writable for knowing if you can write in it. If you don't want to add extra information you can give null and -1 values to extra information arguments
     * @param gs definition of grasp output stream 
     * @param settings current settings used
     * @param segment current processed segment
     * @param output output of processed segment
     * @param tile_description general description of tile
     * @param icol column position of the segment in the tile
     * @param irow row position of the segment in the tile
     * @param itime time position of the segment in the tile
     * @return FILE pointer to stream.
     */
    FILE *grasp_output_stream_open(grasp_output_stream *gs, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief Return true if the stream is writable
     * @param gs definition of grasp output stream 
     * @return true if the stream is writable (there was not error)
     */
    bool grasp_output_stream_writable(const grasp_output_stream *gs);

    /**
     * @brief Return true if the stream is writable file (only if is a file)
     * @param gs definition of grasp output stream 
     * @return rue if the stream is writable FILE (there was not error)
     */
    bool grasp_output_stream_writable_file(const grasp_output_stream *gs);

    /**
     * @brief  Close current open stream
     * @param gs definition of grasp output stream which will be closed
     */
    void grasp_output_stream_close(grasp_output_stream *gs);

    /**
     * @brief It is an implementation of fprinf which uses directly grasp_output_streams
     */
#define gos_fprintf(f_, s_, ...) \
            do { \
               if((f_)->writable==true) fprintf(((f_)->file), (s_), ## __VA_ARGS__); \
            } while (0)

    /**
     * @brief It is an implementation of fprinf which uses directly grasp_output_streams and after printing it force a flush statament
     */
#define gos_fprintf_flushed(f_, s_, ...) \
            do { \
               if((f_)->writable==true) fprintf(((f_)->file), (s_), ## __VA_ARGS__); \
               fflush(((f_)->file)); \
            } while (0)                

    /**
     * @brief It is an implementation of fprinf which uses directly grasp_output_streams and benchmark structures. It is thought to print benchmarks quickly
     */
#define gos_benchmark_print_stream(f_, label_, benchmark_) \
            do { \
               if((f_)->writable==true)  fprintf(((f_)->file), ("%s:%d: %s: clock_time: %f cpu_time: %f\n"),  __FILE__,__LINE__,label_,(benchmark_)->delta_ct, (benchmark_)->delta_ut); \
               fflush(((f_)->file)); \
            } while (0)  


    /**
     * @brief Deallocate memory taken by grasp output stream
     * @param gs definition of grasp output stream
     */
    void grasp_output_stream_destroy(grasp_output_stream *gs);

    /**
     * @brief Dump debug information of current grasp output stream
     * @param gs definition of grasp output stream
     * @param f pointer to file where the information will be dumped
     */
    void grasp_output_stream_debug(grasp_output_stream *gs, FILE *f);

    /**
     * Returns the name of a file applying substitutions to a pattern with wildcards
     * @param original_name complete pattern of the filename
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return A string with all wildcards replaced
     */
    char *grasp_output_stream_filename_generator(const char *original_name, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * typedef definition of a function which implements a specific translation from a token and return the substitution string
     */
    typedef char *(*grasp_output_stream_filename_generate_by)(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    typedef struct grasp_output_stream_filename_generate_by_template_ {
        const char *token;
        grasp_output_stream_filename_generate_by function;
    } grasp_output_stream_filename_generate_by_template;

    /**
     * @brief auto(N): itimexicolxirow with N zeros at the left
     * @param token Token to be translated: auto
     * @param format Format of the token read from pattern: "%dx%dx%d" by default but if N is a number it will represent number of zeros at left of each number. Example for N=3: 000111222
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated
     */
    char *grasp_output_stream_filename_generate_by_auto(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief icol(N): current column number with N zeros at the left
     * @param token Token to be translated: icol
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated
     */
    char *grasp_output_stream_filename_generate_by_icol(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief irow(N): current column number with N zeros at the left
     * @param token Token to be translated: irow
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated
     */
    char *grasp_output_stream_filename_generate_by_irow(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief itime(N): current time number with N zeros at the left
     * @param token Token to be translated: itime
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated
     */
    char *grasp_output_stream_filename_generate_by_itime(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief iinversion(N): current inversion id with N zeros at the left
     * @param token Token to be translated: itime
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated
     */
    char *grasp_output_stream_filename_generate_by_iinversion(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief segment_nx(N): number of X elements per segment with N zeros at the left
     * @param token Token to be translated segment_nx
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_segment_nx(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief segment_ny(N): number of Y elements per segment with N zeros at the left
     * @param token Token to be translated: segment_ny
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_segment_ny(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief segment_nt(N): number of T elements per segment with N zeros at the left
     * @param token Token to be translated: segment_nt
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_segment_nt(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief tile_from(FORMAT): start tile date in FORMAT. By default FORMAT is %FT%H:%M:%SZ
     * @param token Token to be translated tile_from
     * @param format Format of the token read from pattern. By default FORMAT is %FT%H:%M:%SZ
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_from(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief tile_to(FORMAT): final tile date in FORMAT. By default FORMAT is %FT%H:%M:%SZ
     * @param token Token to be translated. tile_to
     * @param format Format of the token read from pattern. By default FORMAT is %FT%H:%M:%SZ
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_to(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief tile_corner_column(N): number of the corner (column) of the tile defined in settings file. Requirement: Input data have to be defined using input.corner instead of input.center
     * @param token Token to be translated. tile_corner_column
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_corner_column(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief tile_corner_row(N): number of the corner of (row) the tile defined in settings file. Requirement: Input data have to be defined using input.corner instead of input.center
     * @param token Token to be translated tile_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_corner_row(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief tile_center_longitude(FORMAT): longitude of the center of the tile defined in settings file. Requirement: Input data have to be defined using input.center instead of input.corner
     * @param token Token to be translated. tile_center_longitude
     * @param format Format of the token read from pattern: F. It has the format of %f C wildcard. %FORMATf 
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_center_longitude(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief tile_center_latitude(FORMAT): latitude of the center of the tile defined in settings file. Requirement: Input data have to be defined using input.center instead of input.corner
     * @param token Token to be translated tile_center_latitude
     * @param format Format of the token read from pattern: F. It has the format of %f C wildcard. %FORMATf
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_center_latitude(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);
    
    /**
     * @brief tile_coordinate_x(I): x input reference of center of the tile defined in settings file. It can be defined by corner or latitude. I is N in case it was defined by corner or 0.I in case it was defined like center
     * @param token Token to be translated. tile_center_longitude
     * @param format Format of the token read from pattern: It can be defined by corner or latitude. I is N in case it was defined by corner or 0.I in case it was defined like center
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_coordinate_x(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
    * @brief tile_coordinate_y(I): y input reference of center of the tile defined in settings file. It can be defined by corner or latitude. I is N in case it was defined by corner or 0.I in case it was defined like center
     * @param token Token to be translated tile_center_latitude
     * @param format Format of the token read from pattern: It can be defined by corner or latitude. I is N in case it was defined by corner or 0.I in case it was defined like center
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_coordinate_y(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);        
    
    /**
     * @brief tile_width(N): Number of X elements in tile with N zeros at the left
     * @param token Token to be translated: tile_width
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_tile_width(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief tile_height(N): Number of Y elements in tile with N zeros at the left
     * @param token Token to be translated tile_height
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated
     */
    char *grasp_output_stream_filename_generate_by_tile_height(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief segment_corner_column(N): number of column of the segment corner with N zeros at the left. Requirement: Input data have to be defined using input.corner instead of input.center
     * @param token Token to be translated segment_corner_column
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_segment_corner_column(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief segment_corner_row(N): number of row of the segment corner with N zeros at the left Requirement: Input data have to be defined using input.corner instead of input.center
     * @param token Token to be translated: segment_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_segment_corner_row(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief settings_filename: the name of settings file used to run the retrieval
     * @param token Token to be translated: segment_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_settings_filename(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief branch: git branch of grasp if it is compiled with saving this information
     * @param token Token to be translated: segment_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */    
    char *grasp_output_stream_filename_generate_by_branch(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief commit: Return commit reference of grasp if it is compiled with saving this information
     * @param token Token to be translated: segment_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */    
    char *grasp_output_stream_filename_generate_by_commit(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);
        
    /**
     * @brief constants_set: Return constants set used in compilation time
     * @param token Token to be translated: segment_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */    
    char *grasp_output_stream_filename_generate_by_constants_set(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);
        
    /**
     * @brief version: Return version of grasp if it is compiled with saving this information
     * @param token Token to be translated: segment_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */    
    char *grasp_output_stream_filename_generate_by_version(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);
        
    /**
     * @brief version: Return first date present in specific segment
     * @param token Token to be translated: segment_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */    
    char *grasp_output_stream_filename_generate_by_segment_first_date(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);
    
    
    /**
     * @brief version: Return last date present in specific segment
     * @param token Token to be translated: segment_corner_row
     * @param format Format of the token read from pattern: N. It has to be a integer number which represents number of 0s at left
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */    
    char *grasp_output_stream_filename_generate_by_segment_last_date(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);
    
    /**
     * @brief pwd: this is replaced by current folder and is only valid at the beginning of the stream definition
     * @param token Token to be translated: pwd
     * @param format Format of the token read from pattern. It is ignored
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_pwd(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief yml: this is replaced by current folder of main configuration file and it is only valid at the beginning of the stream definition    
     * @param token Token to be translated: yml
     * @param format Format of the token read from pattern. It is ignored
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_yml(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);

    /**
     * @brief When the token is known this function is called and will return a copy of input token
     * @param token Token to be translated. It is the value read from pattern
     * @param format Format of the token read from pattern. It is ignored
     * @param settings Current settings
     * @param segment Input segment
     * @param output Output current segment
     * @param tile_description Tile description (dimensions)
     * @param icol Number of column. If it is known you have to use -1
     * @param irow Number of current row of the segment
     * @param itime Number of current time of the segment
     * @return the wildcard translated 
     */
    char *grasp_output_stream_filename_generate_by_token(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime);


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_STREAM_H */

