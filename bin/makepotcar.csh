#!/bin/csh

set MET=`grep -i "method" $1 | cut -d "=" -f 2`
set FUNC=`grep -i "functional" $1 | cut -d "=" -f 2`
set ELEMENTS=`grep "element" $1 | cut -d "=" -f 2 | tr "," " "`
set POTCAR=POTCAR

set VASP=/users/ishi/VASP
set PPDIR=pot${MET}_${FUNC}

touch $POTCAR

foreach ELEM ($ELEMENTS)
	set POTDIR=${VASP}/${PPDIR}/${ELEM}
	if ( -e $POTDIR/POTCAR ) then
		cat $POTDIR/POTCAR >> $POTCAR
	else if ( -e $POTDIR/POTCAR.Z ) then
		zcat $POTDIR/POTCAR.Z >> $POTCAR
	else
		echo "POTCAR does not exist for $POTDIR"
	endif
end
