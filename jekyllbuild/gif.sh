#!/bin/bash

#We go in the img folder, and take the list of gif file to convert
my_path=./img/
cd $my_path;
all_files=$( find . -type f -name '*.gif' )

#We go back at the root, and process to the modificaiton, file per file
cd ./../;

while read line
do
    echo -e "\nAbout to convert $line ..."
    b=$(basename $line)
    # DEBUG # echo "file name : ${b%.*} " 
    
    # We can resize pictures
    #gifsicle -b $my_path/$line -O3 --resize 450x450 -o ./media/compressed/crops/450x450/${b%.*}-450x450.gif
    #gifsicle -b $my_path/$line -O3 --resize 450x253 -o ./media/compressed/crops/450x253/${b%.*}-450x253.gif
    
    # We can crop picture instead of resizing (deformation free :) )
    if [ ! -f ./media/compressed/crops/450x450/${b%.*}-450x450.gif ]; then
        gifsicle -b $my_path/$line -O3 --crop 0,0-450,450 -o ./media/compressed/crops/450x450/${b%.*}-450x450.gif
        echo "1. Done! writen to : ./media/compressed/crops/450x450/${b%.*}-450x450.gif "
    else
        echo "1. File ./media/compressed/crops/450x450/${b%.*}-450x450.gif already exist."
    fi
    if [ ! -f ./media/compressed/crops/450x253/${b%.*}-450x253.gif ]; then
        gifsicle -b $my_path/$line -O3 --crop 0,0-450,253 -o ./media/compressed/crops/450x253/${b%.*}-450x253.gif
        echo "2. Done! writen to : ./media/compressed/crops/450x253/${b%.*}-450x253.gif "
    else
        echo "2. File ./media/compressed/crops/450x253/${b%.*}-450x253.gif already exist."
    fi
    
#More information : https://www.lcdf.org/gifsicle/man.html
done <<< "$all_files"