C
C     ************************************************************
C
      SUBROUTINE TROPUS(NGROUPS,NPTS,N,P,OX,X)
      INTEGER P
      REAL X(4096),OX(NPTS,NGROUPS)
      REAL Y(68000)
      common/savetime/ngo,i2pwr
      common/pulseshape/convolR(8192),CONVOLI(8192)
      COMMON/MISSCHANNELS/MM
      COMMON/DETECT/A,B,E,c,d
      REAL A(96),B(96),c(96),d(96),E(8192)
8     DO 1 I=1,68000
1     Y(I)=0.0
      DO 2 K=1,npts
      SR=0.
      SI=0.
      DO 3 J=1,NGROUPS
      SR=SR+A(J)*OX(K,J)
3     SI=SI+B(J)*OX(K,J)
      Y(K+K-1)=SR*E(K)
      Y(K+K)=SI*E(K)
2     CONTINUE
      CALL FFT(i2pwr+1,Y,-1.)
      DO 4 I=1,N
4     X(I)=(Y(I+I-1)*convolR(I)+Y(I+I)*CONVOLI(I))
      RETURN
      END
