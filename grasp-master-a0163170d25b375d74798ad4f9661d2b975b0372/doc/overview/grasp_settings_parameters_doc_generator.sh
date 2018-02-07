#!/bin/bash

# This script update settings_parameters.dbk file which is used in chapter 5 to show
# all grasp settings parameters. It should be run regularly

# Define path to grasp executable
grasp_exec="../../build/nb/bin/grasp"

# Define temporal file path (then it is complete step by step)
mkdir doc 2> /dev/null
settings_parameters="doc/help_output.txt"

# Obtaining settings parameters from grasp
# Cleaning first and last line
  # Removing the two white-spaces at the beginning of the line
  # Adding </row><row> at the beginning of the line
  # Replacing first ":" by </row><row>
  # Remove first </entry></row>
$grasp_exec help | \
    grep "^  " | \
    sed -e 's/^  //'  -e 's/^/<\/entry><\/row><row><entry>/' -e 's/:/<\/entry><entry>/' -e '1s/^<\/entry><\/row>//' -e 's/\./\. /g'  > $settings_parameters
     

cat src/settings_parameters_template.dbk | sed -e "s/{file_content}/$(<doc/help_output.txt sed -e 's/[\&/]/\\&/g'  | tr -d '\n')/" > src/settings_parameters.dbk



