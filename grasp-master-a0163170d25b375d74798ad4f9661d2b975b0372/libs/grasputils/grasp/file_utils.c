/**
 * @file:  file_utils.c
 * @Author: Fabrice Ducos <fabrice.ducos@univ-lille1.fr>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strncpy */
#include <glob.h>
#include <assert.h>

#include "utils/file_utils.h"

int find_files_from_pattern(const char *pattern, filepath_t **matches) {
  glob_t glob_obj;
  int err;
  size_t num_matches;
  filepath_t *matches_;
  size_t i;
  
  
  err = glob(pattern, 0, NULL, &glob_obj);

  if (err != 0) {
    switch (err) {
    case GLOB_NOMATCH:
      //fprintf(stderr, "%s: no match\n", pattern);
      return -1;
    case GLOB_ABORTED:
      //fprintf(stderr, "%s: glob aborted", pattern);
      //perror(pattern);
      return -2;
    case GLOB_NOSPACE:
      //fprintf(stderr, "fatal error: not enough memory for glob");
      //abort();
      return -3;
    default:
      //fprintf(stderr, "fatal error: unexpected return value (%d) from glob", err);
      //abort();
      return -4;
    }
    return 0;
  }
  
  num_matches = glob_obj.gl_pathc;
  matches_ = malloc(num_matches * sizeof(filepath_t));
  for (i = 0 ; i < num_matches ; i++) {
    strncpy(matches_[i], glob_obj.gl_pathv[i], FILEPATH_LEN_);
  }
  
  globfree(&glob_obj);

  *matches = matches_;
  return num_matches;
}

/* checks if a file exists and is readable */
bool is_readable_file(const char *fname) {
  FILE *fp;

  fp = fopen(fname, "r");
  if (fp != NULL) {
    fclose(fp);
    return true;
  }
  return false;
}

char *pathoffile(const char *file) {
    char *result;
    char *tmp;
    int resultlength;
    
    if (file!=NULL){
        tmp = strrchr(file, '/');
        if (tmp == NULL) {
            result = (char *) malloc(sizeof (char)*1);
            result[0]='\0';
        } else {
            // Allocate space for path string
            resultlength=((tmp-file)/sizeof(char)) +2;
            result = (char *) malloc(sizeof (char)* resultlength);
            
            // Copy length-1 characters
            strncpy(result,file,resultlength-1);
            
            // Stablish end of line
            result[resultlength-1]='\0';
        }
    }else{
        result = (char *) malloc(sizeof (char)*1);
        result[0]='\0';
    }

    return result;
}

char *name_of_file(const char *file){
    char *result=NULL;
    const char *ptr;
    assert(file!=NULL);
    
    result=malloc(strlen(file)+1);
    
    ptr=strrchr(file, '/');
    if(ptr==NULL){
        ptr=file;
    }else{
        ptr++;
    }
    strcpy(result,ptr);
    
    return result;
}

char *name_of_file_without_extension(const char *file){
    char *result=NULL;
    char *ptr;
    assert(file!=NULL);
        
    result=name_of_file(file);
  
    ptr=strrchr(result, '.');
    if(ptr==NULL){
        // Nothing to do
    }else{
        *ptr='\0';
    }
    
    return result;
}

bool isabsolute(const char *file) {
	if(file!=NULL){
		if (file[0] == '/') {
			return true;
		} else {
			return false;
		}
	}else{
		 return false;
	}
}

/*
 Resources for creating recursively directories
 
 http://nion.modprobe.de/blog/archives/357-Recursive-directory-creation.html
 
 snippet:
 
For ii I searched for a method to create directories recursive with all parents. 
 There was no standard c function to do this. mkdir(char *path, mode_t mode) only 
 creates one directory, not the whole path tree. This is what we implemented:
 
static void _mkdir(const char *dir) {
        char tmp[256];
        char *p = NULL;
        size_t len;
 
        snprintf(tmp, sizeof(tmp),"%s",dir);
        len = strlen(tmp);
        if(tmp[len - 1] == '/')
                tmp[len - 1] = 0;
        for(p = tmp + 1; *p; p++)
                if(*p == '/') {
                        *p = 0;
                        mkdir(tmp, S_IRWXU);
                        *p = '/';
                }
        mkdir(tmp, S_IRWXU);
}

Now _mkdir ("/home/nion/irc/server/channel"); works. If someone knows a better and maybe shorter variant please let me know.
 */
