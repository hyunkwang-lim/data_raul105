/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include "grasp_output_stream.h"
#include <grasp/utils.h>
#include "../global/grasp_runtime_information.h"
#include "../global/grasp_compilation_information.h"
#include "../retrieval/constants_set/mod_globals.inc"


grasp_output_stream_filename_generate_by_template grasp_output_stream_wildcard_dictionary[] = {
    {"auto", grasp_output_stream_filename_generate_by_auto},
    {"icol", grasp_output_stream_filename_generate_by_icol},
    {"irow", grasp_output_stream_filename_generate_by_irow},
    {"itime", grasp_output_stream_filename_generate_by_itime},
    {"iinversion", grasp_output_stream_filename_generate_by_iinversion},
    {"segment_nx", grasp_output_stream_filename_generate_by_segment_nx},
    {"segment_ny", grasp_output_stream_filename_generate_by_segment_ny},
    {"segment_nt", grasp_output_stream_filename_generate_by_segment_nt},
    {"tile_from", grasp_output_stream_filename_generate_by_tile_from},
    {"tile_to", grasp_output_stream_filename_generate_by_tile_to},
    {"tile_corner_column", grasp_output_stream_filename_generate_by_tile_corner_column},
    {"tile_corner_row", grasp_output_stream_filename_generate_by_tile_corner_row},
    {"tile_center_longitude", grasp_output_stream_filename_generate_by_tile_center_longitude},
    {"tile_center_latitude", grasp_output_stream_filename_generate_by_tile_center_latitude},
    {"tile_coordinate_x", grasp_output_stream_filename_generate_by_tile_coordinate_x},
    {"tile_coordinate_y", grasp_output_stream_filename_generate_by_tile_coordinate_y},    
    {"tile_width", grasp_output_stream_filename_generate_by_tile_width},
    {"tile_height", grasp_output_stream_filename_generate_by_tile_height},
    {"segment_corner_column", grasp_output_stream_filename_generate_by_segment_corner_column},
    {"segment_corner_row", grasp_output_stream_filename_generate_by_segment_corner_row},
    {"segment_first_date", grasp_output_stream_filename_generate_by_segment_first_date},
    {"segment_last_date", grasp_output_stream_filename_generate_by_segment_last_date},
    {"settings_filename", grasp_output_stream_filename_generate_by_settings_filename},
    {"version", grasp_output_stream_filename_generate_by_version}, 
    {"branch", grasp_output_stream_filename_generate_by_branch}, 
    {"commit", grasp_output_stream_filename_generate_by_commit}, 
    {"constants_set", grasp_output_stream_filename_generate_by_constants_set}, 
    {"pwd", grasp_output_stream_filename_generate_by_pwd},
    {"yml", grasp_output_stream_filename_generate_by_yml},
};

char *grasp_output_stream_filename_generator(const char *original_name, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *copy, *tokens[200];
    int total_length = 0;
    char *result = NULL;
    int i, ntokens;
    // For strtok_t
    char *rest, *token, *ptr; // params for main string tokenizer
    char *wildcard, *format, *rest2, *ptr2; // params for secondary string tokenizer
    char *token_result;
    grasp_output_stream_filename_generate_by function;
    int iknownwildcard, nknownwildcards;
    
    // If pattern is not valid we will return an empty stream
    if(grasp_output_stream_filename_validation(original_name)==false){
        result = (char *) trackmem_malloc(sizeof (char)*2);
        assert(result!=NULL);
        strcpy(result, "");
        return result;
    }
    
    copy = (char *) trackmem_malloc(sizeof (char)*(strlen(original_name) + 1));  
    assert(copy!=NULL);
    strcpy(copy, original_name);

    nknownwildcards = sizeof (grasp_output_stream_wildcard_dictionary) / sizeof (grasp_output_stream_filename_generate_by_template);

    ptr = copy; // initialize main string tokenizer
    i = 0;
    // loop over main tokens in segment_filename
    while ((token = strtok_r(ptr, "{}", &rest))) {
        ptr2 = token; // Initilizating secondary string tokenizer   
        wildcard = strtok_r(ptr2, "()", &rest2);
        ptr2 = rest2;
        format = strtok_r(ptr2, "()", &rest2);
        // Transform wildcard to its corresponding value
        function = NULL;
        for (iknownwildcard = 0; iknownwildcard < nknownwildcards; iknownwildcard++) {
            if (strcmp(wildcard, grasp_output_stream_wildcard_dictionary[iknownwildcard].token) == 0) {
                function = grasp_output_stream_wildcard_dictionary[iknownwildcard].function;
                break;
            }
        }
        if (function == NULL) { // If there was not match we'll assign default one
            function = grasp_output_stream_filename_generate_by_token;
        }

        token_result = function(wildcard, format, settings, segment, output, tile_description, icol, irow, itime);
        // Save token
        tokens[i] = token_result;
        // Analize length
        total_length += strlen(tokens[i]);
        // Prepare to next loop
        ptr = rest;
        i++;
    }
    // Save total of tokens
    ntokens = i;

    // result length is equal to total length plus end of string '\0' character
    total_length += 1;

    // Allocate result and initialize it
    result = (char *) trackmem_malloc(sizeof (char)*total_length);
    assert(result != NULL);
    strcpy(result, "");

    // Set result and clean memory
    for (i = 0; i < ntokens; i++) {
        strcat(result, tokens[i]);
        trackmem_free(tokens[i]);
    }
    trackmem_free(copy);

    return result;
}

char *grasp_output_stream_filename(grasp_output_stream *gs, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *filename;

    filename = grasp_output_stream_filename_generator(gs->filename, settings, segment, output, tile_description, icol, irow, itime);

    // Check if file exists. 
    if (access(filename, F_OK) != -1) {
        fprintf(stderr, "WARNING: File %s already exists and it will be overwritten\n", filename);
    }

    gs->writable = true;

    return filename;
}

bool grasp_output_stream_filename_validation(const char *original_name) {
    char *copy;
    //char *tmp;
    int init_wildcard, end_wildcard, total_token_length,total_length, i, tmp, c;
    // For strtok_t
    char *rest, *token, *ptr; // params for main string tokenizer
    char *wildcard, *rest2, *ptr2; //, *format; // params for secondary string tokenizer
    int iknownwildcard, nknownwildcards;
    bool found;
    assert(original_name != NULL);
    total_length = strlen(original_name);
    nknownwildcards = sizeof (grasp_output_stream_wildcard_dictionary) / sizeof (grasp_output_stream_filename_generate_by_template);
    copy = (char *) trackmem_malloc(sizeof (char)*(strlen(original_name) + 1));
    assert( copy!= NULL);
    strcpy(copy, original_name);

    // Validate that the brackets open and close properly
    tmp = 0;
    c = 0;
    for (i = 0; i < total_length; i++) {
        if (original_name[i] == '{') {
            tmp++;
            c = 0;
        }
        if (original_name[i] == '}') {
            tmp--;
            if (c < 2) {
                return false;
            }
            c = 0;
        }
        if (tmp < 0 || tmp > 1) {
            return false;
        }
        c++;
    }
    if (tmp != 0) {
        return false;
    }

    ptr = copy; // initialize main string tokenizer
    // loop over main tokens in segment_filename
    while ((token = strtok_r(ptr, "{}", &rest))) {
        tmp = 0;
        c = 0;
        total_token_length = strlen(token);
        for (i = 0; i < total_token_length; i++) {
            if (token[i] == '(') {
                tmp++;
                c = 0;
                ;
            }
            if (token[i] == ')') {
                tmp--;
                if (c < 2) {
                    return false;
                }
                c = 0;
            }
            if (tmp < 0 || tmp > 1) {
                return false;
            }
            c++;
        }
        if (tmp != 0) {
            return false;
        }
        //   validate if wildcard exists. We try to find  '{' and '}' symbols that should be in indexes init_wildcard and end_wildcard  
        init_wildcard = 0;
        while (&copy[init_wildcard] != &token[0]) {
            init_wildcard++; // Looking for the position of the char
        }
        end_wildcard = strlen(token);
        init_wildcard--;
        end_wildcard++;
        end_wildcard = end_wildcard + init_wildcard;

        ptr2 = token; // Initilizating secondary string tokenizer   
        wildcard = strtok_r(ptr2, "()", &rest2);
        ptr2 = rest2;
        //format =   strtok_r(ptr2, "()", &rest2);

        //check than format is a number
        /* I remove this because also could be date string format
         * if(format!=NULL){
            for(i=0;i<strlen(format);i++){
                if(format[i]!='0' &&
                   format[i]!='1' &&
                   format[i]!='2' &&
                   format[i]!='3' &&
                   format[i]!='4' &&
                   format[i]!='5' &&
                   format[i]!='6' &&
                   format[i]!='7' &&
                   format[i]!='8' &&
                   format[i]!='9' 
                        ){
                    return false; 
                }
            }
        }*/
        // printf("evaluating %s with %d=%c and %d=%c\n", wildcard,init_wildcard, original_name[init_wildcard], end_wildcard,original_name[end_wildcard] );
        if (init_wildcard >= 0 && end_wildcard < total_length && original_name[init_wildcard] == '{' && original_name[end_wildcard] == '}') {
            found=false;
            for (iknownwildcard = 0; iknownwildcard < nknownwildcards; iknownwildcard++) {
                if (strcmp(token, grasp_output_stream_wildcard_dictionary[iknownwildcard].token) == 0) {
                    found=true;
                }
            }
            if(found==false){
                return false;
            }
            // the wilcards pwd and yml have to be at the beginning of the pattern. We will check it
            if (strcmp(wildcard, "pwd") == 0) {
                if (strlen(original_name) > 5 && original_name[0] == '{' && original_name[1] == 'p' && original_name[2] == 'w' && original_name[3] == 'd' && original_name[4] == '}') {
                    // ok
                } else {
                    return false;
                }
            } 
            if (strcmp(wildcard, "yml") == 0) {
                if (strlen(original_name) > 5 && original_name[0] == '{' && original_name[1] == 'y' && original_name[2] == 'm' && original_name[3] == 'l' && original_name[4] == '}') {
                    // ok
                } else {
                    return false;
                }
            }
        }
        ptr = rest;
    }

    trackmem_free(copy);

    return true;
}

void grasp_output_stream_initialize(const char *filename, grasp_output_stream *gs) {
    char *c;

    assert(filename != NULL && strcmp(filename, "") != 0);

    gs->filename = (char *) trackmem_malloc(sizeof (char)*(strlen(filename) + 1));
    assert( gs->filename!= NULL);
    strcpy(gs->filename, filename);
    gs->file = NULL;
    c = strtolower(filename);

    if (strcmp(c, "screen") == 0 || strcmp(c, "stdout") == 0 || strcmp(c, "true") == 0 || strcmp(c, "t") == 0 || strcmp(c, "1") == 0) {
        gs->none = false;
        gs->screen = true;
        gs->writable = false;
        gs->open = false;
    } else if (strcmp(c, "none") == 0 || strcmp(c, "null") == 0 || strcmp(c, "false") == 0 || strcmp(c, "f") == 0 || strcmp(c, "0") == 0) {
        gs->none = true;
        gs->screen = false;
        gs->writable = false;
        gs->open = false;
    } else {
        gs->none = false;
        gs->screen = false;
        gs->writable = false;
        gs->open = false;
    }

    trackmem_free(c);
}

FILE *grasp_output_stream_open(grasp_output_stream *gs, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *filename;
    assert(gs->open == false);
    if (gs->none == false) {
        gs->writable = true;
        gs->open = true;
        if (gs->screen == true) {
            gs->file = stdout;
        } else {
            filename = grasp_output_stream_filename(gs, settings, segment, output, tile_description, icol, irow, itime);
            gs->file = fopen(filename, "w");
            if (gs->file == NULL) {
                gs->writable = false;
                fprintf(stderr, "ERROR: File %s can not be opened\n", filename);
            }
            trackmem_free(filename);
        }
    } else {
        gs->writable = false;
        gs->open = true;
        gs->file = NULL;
    }

    return gs->file;
}

bool grasp_output_stream_writable(const grasp_output_stream *gs) {
    return gs->writable;
}

bool grasp_output_stream_writable_file(const grasp_output_stream *gs) {
    if (gs->writable == true && gs->screen == false) {
        return true;
    } else {
        return false;
    }
}

void grasp_output_stream_close(grasp_output_stream *gs) {
    assert(gs->open == true);

    if (gs->writable == true) {
        gs->writable = false;
        gs->open = false;
        if (gs->screen == false) {
            fclose(gs->file);
        }
    }
}

void grasp_output_stream_destroy(grasp_output_stream *gs) {
    trackmem_free(gs->filename);
}

void grasp_output_stream_debug(grasp_output_stream *gs, FILE *f) {
    fprintf(f, "filename=%s; none=%d; screen=%d; file=%p; writable=%d; open=%d\n"
            , gs->filename
            , gs->none
            , gs->screen
            , gs->file
            , gs->writable
            , gs->open);
}

char *grasp_output_stream_filename_generate_by_auto(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*20);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (icol < 0 || irow < 0 || itime < 0 || strcmp(token, "auto") != 0) {
        return token_result;
    }

    if (format == NULL) {
        strcpy(used_format, "%dx%dx%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
        strcat(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
        strcat(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, itime, icol, irow);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_icol(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    strcpy(token_result, "");
    if (icol < 0 || strcmp(token, "icol") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue        
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, icol);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_irow(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (irow < 0 || strcmp(token, "irow") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue          
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, irow);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_itime(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (itime < 0 || strcmp(token, "itime") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue          
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, itime);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_iinversion(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (segment==NULL || segment->sdata.id_inversion < 0 || strcmp(token, "iinversion") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue          
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, segment->sdata.id_inversion);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_segment_nx(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || strcmp(token, "segment_nx") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue           
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, settings->input.segment_nx);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_segment_ny(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || strcmp(token, "segment_ny") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue 
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, settings->input.segment_ny);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_segment_nt(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || strcmp(token, "segment_nt") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue 
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, settings->input.segment_nt);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_tile_from(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    time_t time;
    int err;

    token_result = (char *) trackmem_malloc(sizeof (char)*35);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || strcmp(token, "tile_from") != 0) {
        return token_result;
    }
    err = convert_string_to_time(settings->input.time_from, &time, TIMEFMT_ISO8601);
    if (err != 0) {
        printf("%s:%d: Wrong format input date\n", __FILE__, __LINE__);
        abort();
    }
    if (format == NULL) {
        strcpy(used_format, "%FT%H:%M:%SZ");
    } else {
        strcpy(used_format, format);
    }
    time_to_string_r(time, used_format, sizeof (char)*35, token_result);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_tile_to(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    time_t time;
    int err;

    token_result = (char *) trackmem_malloc(sizeof (char)*35);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || strcmp(token, "tile_to") != 0) {
        return token_result;
    }
    err = convert_string_to_time(settings->input.time_to, &time, TIMEFMT_ISO8601);
    if (err != 0) {
        printf("%s:%d: Wrong format input date\n", __FILE__, __LINE__);
        abort();
    }
    if (format == NULL) {
        strcpy(used_format, "%FT%H:%M:%SZ");
    } else {
        strcpy(used_format, format);
    }
    time_to_string_r(time, used_format, sizeof (char)*35, token_result);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_tile_corner_column(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || settings->input.coordinates_grid.col < 0 || strcmp(token, "tile_corner_column") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue 
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, settings->input.coordinates_grid.col - settings->input.grid_offset.col);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_tile_corner_row(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || settings->input.coordinates_grid.row < 0 || strcmp(token, "tile_corner_row") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue 
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, settings->input.coordinates_grid.row - settings->input.grid_offset.row);

    return token_result;
}


char *grasp_output_stream_filename_generate_by_tile_center_longitude(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    float satof, satof_r;
    
    token_result = (char *) trackmem_malloc(sizeof (char)*20);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || settings->input.coordinates_grid.col > 0 || strcmp(token, "tile_center_longitude") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%f");
    } else {
        // Validate format
        satof_r = safe_atofloat(format, &satof);
        if (satof_r != 0) {
            return token_result;
        }
        // If it valid we continue 
        strcpy(used_format, "%");
        strcat(used_format, format);
        strcat(used_format, "f");
    }
    sprintf(token_result, used_format, (settings->input.coordinates.lon));
    
    return token_result;
}

char *grasp_output_stream_filename_generate_by_tile_center_latitude(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    float satof, satof_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*20);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || settings->input.coordinates_grid.row > 0 || strcmp(token, "tile_center_latitude") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%f");
    } else {
        // Validate format
        satof_r = safe_atofloat(format, &satof);
        if (satof_r != 0) {
            return token_result;
        }
        // If it valid we continue 
        strcpy(used_format, "%");
        strcat(used_format, format);
        strcat(used_format, "f");
    }
    sprintf(token_result, used_format, (settings->input.coordinates.lat));

    return token_result;
}


char *grasp_output_stream_filename_generate_by_tile_coordinate_x(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result, *new_format, *token_corner="tile_corner_column", *token_center="tile_center_longitude";
    
    if (settings == NULL || strcmp(token, "tile_coordinate_x") != 0) {
        token_result = (char *) trackmem_malloc(sizeof (char)*5);
        assert( token_result!= NULL);
        strcpy(token_result, "");
        return token_result;
    }
    
    if(settings->input.coordinates_grid.col >= 0 ){
        return grasp_output_stream_filename_generate_by_tile_corner_column(token_corner,format,settings,segment,output,tile_description,icol,irow,itime);
    }else{
        if(format!=NULL){
            new_format = (char *) trackmem_malloc(sizeof (char)*(strlen(format)+3));
            assert( new_format!= NULL);
            strcpy(new_format,".");
            strcat(new_format,format);
        }else{
            new_format=NULL;
        }
        token_result= grasp_output_stream_filename_generate_by_tile_center_longitude(token_center,new_format,settings,segment,output,tile_description,icol,irow,itime);    
        trackmem_free(new_format);
        return token_result;
    }
}

char *grasp_output_stream_filename_generate_by_tile_coordinate_y(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result, *new_format, *token_corner="tile_corner_row", *token_center="tile_center_latitude";
    
    if (settings == NULL || strcmp(token, "tile_coordinate_y") != 0) {
        token_result = (char *) trackmem_malloc(sizeof (char)*5);
        assert( token_result!= NULL);
        strcpy(token_result, "");
        return token_result;
    }
    
    if(settings->input.coordinates_grid.row >= 0 ){
        return grasp_output_stream_filename_generate_by_tile_corner_row(token_corner,format,settings,segment,output,tile_description,icol,irow,itime);
    }else{
        if(format!=NULL){
            new_format = (char *) trackmem_malloc(sizeof (char)*(strlen(format)+3));
            assert( new_format!= NULL);
            strcpy(new_format,".");
            strcat(new_format,format);
        }else{
            new_format=NULL;
        }
        token_result= grasp_output_stream_filename_generate_by_tile_center_latitude(token_center,new_format,settings,segment,output,tile_description,icol,irow,itime);    
        trackmem_free(new_format);
        return token_result;
    }
}

char *grasp_output_stream_filename_generate_by_tile_width(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || settings->input.area_width < 0 || strcmp(token, "tile_width") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue 
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, (settings->input.area_width));

    return token_result;
}

char *grasp_output_stream_filename_generate_by_tile_height(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || settings->input.area_height < 0 || strcmp(token, "tile_height") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue 
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, (settings->input.area_height));

    return token_result;
}

char *grasp_output_stream_filename_generate_by_segment_corner_column(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || settings->input.coordinates_grid.col < 0 || icol < 0 || strcmp(token, "segment_corner_column") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, (settings->input.coordinates_grid.col + (icol * settings->input.segment_nx) - settings->input.grid_offset.col));

    return token_result;
}

char *grasp_output_stream_filename_generate_by_segment_corner_row(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    int satoi, satoi_r;

    token_result = (char *) trackmem_malloc(sizeof (char)*10);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || settings->input.coordinates_grid.row < 0 || icol < 0 || strcmp(token, "segment_corner_row") != 0) {
        return token_result;
    }
    if (format == NULL) {
        strcpy(used_format, "%d");
    } else {
        // Validate format
        satoi_r = safe_atoi(format, &satoi);
        if (satoi_r != 0) {
            return token_result;
        }
        // If it valid we continue
        strcpy(used_format, "%0");
        strcat(used_format, format);
        strcat(used_format, "d");
    }
    sprintf(token_result, used_format, (settings->input.coordinates_grid.row + (irow * settings->input.segment_ny) - settings->input.grid_offset.row));

    return token_result;
}

char *grasp_output_stream_filename_generate_by_segment_first_date(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    
    token_result = (char *) trackmem_malloc(sizeof (char)*35);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || strcmp(token, "segment_first_date") != 0) {
        return token_result;
    }

    if (format == NULL) {
        strcpy(used_format, "%FT%H:%M:%SZ");
    } else {
        strcpy(used_format, format);
    }
    time_to_string_r(segment->sdata.pixel[0].t, used_format, sizeof (char)*35, token_result);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_segment_last_date(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;
    char used_format[50];
    
    token_result = (char *) trackmem_malloc(sizeof (char)*35);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (settings == NULL || strcmp(token, "segment_last_date") != 0) {
        return token_result;
    }

    if (format == NULL) {
        strcpy(used_format, "%FT%H:%M:%SZ");
    } else {
        strcpy(used_format, format);
    }
    time_to_string_r(segment->sdata.pixel[segment->sdata.npixels-1].t, used_format, sizeof (char)*35, token_result);

    return token_result;
}


char *grasp_output_stream_filename_generate_by_settings_filename(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime){
    char *token_result;

    if (strcmp(token, "settings_filename") != 0) {
        token_result = (char *) trackmem_malloc(sizeof (char)* _GBL_FILE_PATH_LEN);
        assert( token_result!= NULL);
        strcpy(token_result, "");        
        return token_result;
    }
    token_result=name_of_file_without_extension(grasp_main_settings_file);

    return token_result; // WARNING: We are returning memory allocated without trackmem library
}

char *grasp_output_stream_filename_generate_by_version(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime){
    char *token_result;
     
    if (strcmp(token, "version") != 0) {
        token_result = (char *) trackmem_malloc(sizeof (char)* 3);
        assert( token_result!= NULL);
        strcpy(token_result, "");           
        return token_result;
    }

    return grasp_compilation_information_version(); 
}

char *grasp_output_stream_filename_generate_by_branch(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime){
    char *token_result;
  
    if (strcmp(token, "branch") != 0) {
        token_result = (char *) trackmem_malloc(sizeof (char)* 3);
        assert( token_result!= NULL);
        strcpy(token_result, "");           
        return token_result;
    }
    
    return grasp_compilation_information_branch_name(); 
}

char *grasp_output_stream_filename_generate_by_commit(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime){
    char *token_result;
  
    if (strcmp(token, "commit") != 0) {
        token_result = (char *) trackmem_malloc(sizeof (char)* 3);
        assert( token_result!= NULL);
        strcpy(token_result, "");           
        return token_result;
    }

    return grasp_compilation_information_commit_ref(); 
}
        
char *grasp_output_stream_filename_generate_by_constants_set(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime){
    char *token_result;

    if (strcmp(token, "constants_set") != 0) {
        token_result = (char *) trackmem_malloc(sizeof (char)* 3);
        assert( token_result!= NULL);
        strcpy(token_result, "");           
        return token_result;
    }

    return grasp_compilation_information_constants_set(); 
}
        
char *grasp_output_stream_filename_generate_by_pwd(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;

    token_result = (char *) trackmem_malloc(sizeof (char)* _GBL_FILE_PATH_LEN);
    assert( token_result!= NULL);
    strcpy(token_result, "");
    if (strcmp(token, "pwd") != 0) {
        return token_result;
    }
    strcpy(token_result, grasp_current_path);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_yml(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;

    if (strcmp(token, "yml") != 0) {
        token_result = (char *) trackmem_malloc(sizeof (char)*5);
        assert( token_result!= NULL);
        strcpy(token_result, "");
        return token_result;
    }
    token_result = pathoffile(grasp_main_settings_file);

    return token_result;
}

char *grasp_output_stream_filename_generate_by_token(const char *token, const char *format, const grasp_settings *settings, const grasp_segment_t *segment, const output_segment_general *output, const grasp_tile_dimensions_t *tile_description, int icol, int irow, int itime) {
    char *token_result;

    if (format != NULL) {
        token_result = (char *) trackmem_malloc(sizeof (char)*(strlen(token) + strlen(format) + 1));
        assert( token_result!= NULL);
        sprintf(token_result, "%s(%s)", token, format);
    } else {
        token_result = (char *) trackmem_malloc(sizeof (char)*(strlen(token) + 1));
        assert( token_result!= NULL);
        strcpy(token_result, token);
    }

    return token_result;
}
