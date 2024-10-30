awk -v a="$2" '{filename = "split." int((NR-1)/a) ".xyz"; print >> filename}' $1
