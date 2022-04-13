#!/bin/bash
rm -v *iwfmt
for file in *d1trfl.*
do
    if [ -f $file ]; then
       echo $file
       echo -e "$file\n 1\n" |$COLUMBUS/iwfmt.x > $file.iwfmt 2> err.ignore
    fi
done