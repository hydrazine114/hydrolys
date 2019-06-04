Subroutine DihedralAngle(RA,RB,RC,RD,iRad,Theta)

Implicit Real(8) (A-H,O-Z)

Real(8) RA(3),RB(3),RC(3),RD(3),A(3),B(3),C(3),D(3),X(3),Y(3),Z(3),U(3)

!
! Subroutine calculates dihedral angle for 4 atoms given by their ccordinates RA,RB,RC,RD (A-B-C-D)
! with normal conventions of QC: A - atom X, BCD - reference atoms
!

A=RA-RB
B=RC-RB
C=-B
D=RD-RC

Call CrossProduct(A,B,X,XNorm)
Call CrossProduct(C,D,Y,YNorm)

X=X/XNorm
Y=Y/YNorm
ct=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)

Call CrossProduct(B,X,Z,ZNorm)
Z=Z/ZNorm
st=Z(1)*Y(1)+Z(2)*Y(2)+Z(3)*Y(3)

If (DABS(DABS(ct)-1.d0)<1.d-6) ct=DSIGN(1.d0,ct)
If (iRad>0) Then
	Theta=DACOS(ct)
Else
	Theta=DACOSD(ct)
Endif
Theta1=Theta
If (st<0.d0) Theta=-Theta

End
