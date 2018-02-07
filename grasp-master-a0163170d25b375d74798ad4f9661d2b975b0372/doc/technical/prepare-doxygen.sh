#!/bin/bash

if [ "$#" != "1" ]
then
    echo "Usage: prepare-doxygen.sh [format=html|php]"
fi

format=$1


# PREPARING DOXYFILE FROM TEMPLATE

cp Doxyfile_template Doxyfile
sed -i -e 's/{{format}}/'${format}'/g' Doxyfile



# PREPARING MARKDOWN FILES FOR DOXYGEN

input_dir='markdown'
output_dir='markdown-doxygen'

if [ ! -d "$input_dir" ] ; then
  echo "$0: $input_dir: no such directory" 1>&2
  exit 1
fi

if [ -d "$output_dir" ] ; then
  rm -rf "$output_dir"
fi

cp -r $input_dir $output_dir

find $output_dir -name '*.md' | while read input_file
do
  sed -i -e '
  s/```c/~~~\{.c\}/g;s/```/~~~/
  s/```c/~~~\{.c\}/g;s/```/~~~/' $input_file
done # input_file

exit 0