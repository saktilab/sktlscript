nAt=$(head -n 1 $1)
grep -A $nAt -B 1 "energy" $1 | grep -v -- "^--$"
