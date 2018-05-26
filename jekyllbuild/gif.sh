#!/bin/bash

#We go in the img folder, and take the list of gif file to convert
my_path=./img/
cd $my_path;
all_files=$( find . -type f -name '*.gif' )

#We go back at the root, and process to the modificaiton, file per file
cd ./../;

while read line
do
    echo "About to convert $line ..."
    b=$(basename $line)
    # DEBUG # echo "file name : ${b%.*} " 
    
    # We can resize pictures
    #gifsicle -b $my_path/$line -O3 --resize 450x450 -o ./media/compressed/crops/450x450/${b%.*}-450x450.gif
    #gifsicle -b $my_path/$line -O3 --resize 450x253 -o ./media/compressed/crops/450x253/${b%.*}-450x253.gif
    
    # We can crop picture instead of resizing (deformation free :) )
    gifsicle -b $my_path/$line -O3 --crop 0,0-450,450 -o ./media/compressed/crops/450x450/${b%.*}-450x450.gif
    gifsicle -b $my_path/$line -O3 --crop 0,0-450,253 -o ./media/compressed/crops/450x253/${b%.*}-450x253.gif

#More information : https://www.lcdf.org/gifsicle/man.html

    echo "Done! writen to : ./media/compressed/crops/450x253/${b%.*}-450xXXX.gif "
done <<< "$all_files"