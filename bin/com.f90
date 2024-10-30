subroutine CenterOfMass( Length, X, Y, Z, TYPEg, Ngrs, MASS, xcom, ycom, zcom )

implicit none

! This routine calculates the center of mass of a chain of molecules.

! Length is the number of beads in the molecule.
! X, Y, and Z are the coordinates of the molecules.
! TYPEg contains the group identity of each molecule.

integer, intent(in)								:: Length
real, dimension(Length), intent(in)				:: X, Y, Z
integer, dimension(Length), intent(in)			:: TYPEg

! Ngrs is the number of groups within the molecule.
! MASS is the mass of each group.

integer, intent(in)								:: Ngrs
real, dimension(Ngrs), intent(in)				:: MASS

! xcom, ycom, and zcom are the coordinates of the center of mass of the molecule.

real, intent(out)								:: xcom, ycom, zcom

! Local Variables

integer											:: i
real											:: MassSum
real											:: XMSum, YMSum, ZMSum

MassSum = 0.0
XMSum = 0.0
YMSum = 0.0
ZMSum = 0.0

do i = 1, Length
	
	MassSum = MassSum + MASS( TYPEg(i) )

	XMSum = XMSum + X(i) * MASS( TYPEg(i) )
	YMSum = YMSum + Y(i) * MASS( TYPEg(i) )
	ZMSum = ZMSum + Z(i) * MASS( TYPEg(i) )

end do

xcom = XMSum / MassSum
ycom = YMSum / MassSum
zcom = ZMSum / MassSum

return

end	subroutine CenterOfMass
