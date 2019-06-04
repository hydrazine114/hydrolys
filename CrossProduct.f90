Subroutine CrossProduct(X,Y,Z,ZNorm)

Implicit Real(8) (A-H,O-Z)

Real(8) X(3),Y(3),Z(3)

Z(1)=X(2)*Y(3)-Y(2)*X(3)
Z(2)=X(3)*Y(1)-X(1)*Y(3)
Z(3)=X(1)*Y(2)-Y(1)*X(2)

ZNorm=DSQRT(Z(1)**2+Z(2)**2+Z(3)**2)

End
