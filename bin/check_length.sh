for i in $@;
do 
ffmpeg -i $i 2>&1 | grep Duration | cut -d ' ' -f 4 | sed s/,//
done
