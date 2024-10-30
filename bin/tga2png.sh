for i in *.tga; do name=$(echo "$i" | cut -f 1 -d '.'); convert $name.tga $name.png; done
