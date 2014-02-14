C
C     ************************************************************
C
      SUBROUTINE CHOSOL(A,B,N,X)
      REAL L(3,3),A(3,3),BL(3),B(3),X(3)
      L(1,1)=SQRT(A(1,1))
      DO 35 I=2,N
      L(I,1)=A(I,1)/L(1,1)
      DO 35 J=2,I
      Z=0.
      DO 36 K=1,J-1
36    Z=Z+L(I,K)*L(J,K)
      Z=A(I,J)-Z
C       TRAP FOR NEGATIVE SQUARE ROOT
      IF(Z.LE.0.)Z=1.E-10
      IF(J.EQ.I)L(I,J)=SQRT(Z)
      IF(J.NE.I)L(I,J)=Z/L(J,J)
35    CONTINUE
      BL(1)=B(1)/L(1,1)
      DO 37 I=2,N
      Z=0.
      DO 38 K=1,I-1
38    Z=Z+L(I,K)*BL(K)
37    BL(I)=(B(I)-Z)/L(I,I)
      X(N)=BL(N)/L(N,N)
      DO 39 I1=1,N-1
      I=N-I1
      Z=0.
      DO 40 K=I+1,N
40    Z=Z+L(K,I)*X(K)
39    X(I)=(BL(I)-Z)/L(I,I)
      RETURN
      END
