/* 
 * File:   yamlsettings_parse_settings_file_mode.h
 * Author: fuertes
 *
 * Created on March 7, 2015, 4:14 PM
 */

#ifndef YAMLSETTINGS_PARSE_SETTINGS_FILE_MODE_H
#define	YAMLSETTINGS_PARSE_SETTINGS_FILE_MODE_H

#ifdef	__cplusplus
extern "C" {
#endif


// This enumeration defines the way to parse settings. If YS_PARSE_SETTINGS_FILE_MANDATORY
// this file have to be defined like first argument. If YS_PARSE_SETTINGS_FILE_FORBIDDEN
// first argument never will be read like input settings file. If YS_PARSE_SETTINGS_FILE_OPTIONAL
// first argument will be readed like settings file but if the file can not be opened 
// it will try to read like a parameter. Warning: last option could be a problem if there is a
// file with same name than a parameter...
typedef enum  { YS_PARSE_SETTINGS_FILE_MANDATORY , YS_PARSE_SETTINGS_FILE_FORBIDDEN, YS_PARSE_SETTINGS_FILE_OPTIONAL } yamlsettings_parser_settings_file_mode;
 

#ifdef	__cplusplus
}
#endif

#endif	/* YAMLSETTINGS_PARSE_SETTINGS_FILE_MODE_H */

