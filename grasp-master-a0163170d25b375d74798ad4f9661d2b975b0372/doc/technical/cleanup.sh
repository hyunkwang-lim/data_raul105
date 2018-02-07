#!/bin/bash
grep -v "8md\.xml" docbook/index.xml > temp.xml
mv temp.xml docbook/index.xml
