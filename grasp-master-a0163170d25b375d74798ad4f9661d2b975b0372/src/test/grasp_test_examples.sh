# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

#This script will run all tests and check that they exit returning 0 code.

grasp_app="./../../build/bin/grasp"
grasp_args=""
#grasp_args="controller.debug.perform_retrieval=false"

if [ ! -e $grasp_app ]
then
    echo "You have to have GRASP compiled to run this script"
    exit -1;
fi

files=$(find ./../../examples -name settings_example_*.yml)
files+=" "$(find ./../../src -name settings_example_*.yml)
files+=" "$(find ./../../internal_examples -name settings_example_*.yml)

nfiles=$(echo $files | wc -w)

files_with_problems=""
problems=0

i=0

for file in $files
do
    percent=$(echo "scale=2; $i*100/$nfiles" | bc)
    printf "%6s" "$percent"
    echo -n "% Testing "$file" "
    $grasp_app $file $grasp_args > /dev/null 2>&1
    if [ $? -ne 0 ]
    then
       files_with_problems+=" "$file 
       let problems++
       echo "FAIL"
    else 
       echo "OK"
    fi
    let i++
done

if [ $problems -ne 0 ]
then
    #echo "List of problematic files"
    #for file in $files_with_problems
    #do
    #    echo " "$file
    #done
    echo "There were "$problems" unsuccessful tests"
    exit -1
else
    echo "All tests were run successfully"
    exit 0
fi

