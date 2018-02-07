/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */


/* File Path Length */

#define _GBL_FILE_PATH_LEN 255

#define G_ERROR(msg) \
  write(error_unit,'(/,"A new error occured at: ",A,":",I0,":",A,A)') __FILE__, __LINE__, NEW_LINE('A'), trim(msg) ; \
  if (stop_report%status) then ; \
    write(error_unit,'("Fatal error: it masks an old error at: ",A,":",I0,":",A,A)') \
	  trim(stop_report%file), stop_report%line, NEW_LINE('A'), trim(stop_report%message) ; \
	stop ; \
  end if ; \
  stop_report%status = .true. ; \
  stop_report%file = __FILE__ ; \
  stop_report%line = __LINE__ ; \
  stop_report%message = msg ; \
  return
