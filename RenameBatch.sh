for f in *.png; 
    do 
    b=$(basename $f)
    # DEBUG # echo "file name : ${b%.*} " 
    cp "$f" "${f#*.}"; 
done
