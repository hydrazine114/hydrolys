Subroutine RotAroundRvec(R,Angle,iRad,R0,RX)
Implicit Real(8) (A-H,O-Z)

Real(8) R0(3),RX(3),R(3),U(3,3),RR(3)

!
! RotAroundRvec rotates vector R0 around arbitrary vector R by Angle (radians if iRad==1, degrees if iRad==0)
! Rx = U * R0  (U - rotation matrix) 
! Algorithm is based on the quaternion algebra 
!

Pi=3.141592653589793d0
ang=Angle
If (iRad==0) ang=Angle*Pi/180.d0

Rnorm=DSQRT(R(1)**2+R(2)**2+R(3)**2)
RR=R/Rnorm

!
! Construct quaternion Q=(x,y,z,w)
!
sa = DSIN(0.5d0*ang)
ca = DCOS(0.5d0*ang)

x  = RR(1) * sa
y  = RR(2) * sa
z  = RR(3) * sa
w  = ca

!
!quaternion_normalise( q );
!
qnorm=DSQRT(X*X+Y*Y+Z*Z+W*W)
x=x/qnorm
y=y/qnorm
z=z/qnorm
w=w/qnorm


! Quaternion transformation to the 3x3 rotation matrix
!
!        |       2     2                                |
!        | 1 - 2Y  - 2Z    2XY - 2ZW      2XZ + 2YW     |
!        |                                              |
!        |                       2     2                |
!    U = | 2XY + 2ZW       1 - 2X  - 2Z   2YZ - 2XW     |
!        |                                              |
!        |                                      2     2 |
!        | 2XZ - 2YW       2YZ + 2XW      1 - 2X  - 2Y  |
!        |                                              |

!
! Rotation matrix
!
U(1,1)=0.5d0-Y*Y-Z*Z
U(1,2)=X*Y-Z*W
U(1,3)=X*Z+Y*W
U(2,1)=X*Y+Z*W
U(2,2)=0.5d0-X*X-Z*Z
U(2,3)=Y*Z-X*W
U(3,1)=X*Z-Y*W
U(3,2)=Y*Z+X*W
U(3,3)=0.5d0-X*X-Y*Y

U=U*2.d0

!
! Rotation
!
!Rx=matmul(U,R0)
!
Do k=1,3
	RX(k)=U(k,1)*R0(1)+U(k,2)*R0(2)+U(k,3)*R0(3)
Enddo

End
