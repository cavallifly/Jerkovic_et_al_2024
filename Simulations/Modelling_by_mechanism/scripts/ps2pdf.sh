#!/usr/bin/env bash

files=$(ls *.ps 2> /dev/null); 
if [[ $1 != "" ]];
then
    files=$(ls $1 2> /dev/null)
fi

for file in ${files};
do 
    echo "Converting "$file

    ps2pdf $file
    
    ~/pdfcrop.pl --margin 15 ${file%.ps}.pdf _${file%.ps}.pdf #&> /dev/null
    mv _${file%.ps}.pdf ${file%.ps}.pdf
    convert -background white -alpha remove -density 600 ${file%ps}pdf ${file%ps}png
    
    #rm $file ${file%ps}pdf
    rm $file
done
