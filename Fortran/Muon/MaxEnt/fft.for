C
C*****************************************************************
C
      SUBROUTINE FFT(KK,A,SN)
      COMPLEX A(1),B,EQ,EL,EN
      DIMENSION INU(16)
      DATA INU/1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,
     + 16384,32768/
      N=INU(KK+1)
      N2=N/2
      E=SN*3.1415926536
      J=1+N2
      DO 1 I=2,N
      K=KK-1
      M=N2
      IF(I-J)2,4,4
2     B=A(J)
      A(J)=A(I)
      A(I)=B
4     L=J-M
      IF(L)5,5,3
3     J=L
      M=INU(K)
      K=K-1
      IF(K)5,5,4
5     J=J+M
1     CONTINUE
      K=1
6     IF(K-N)9,8,8
9     L=K+K
      C2=COS(E)
      EQ=CMPLX(C2,SIN(E))
      EL=(1.,0.)
      C2=C2+C2
      DO 10 M=1,K
      DO 7 I=M,N,L
      J=I+K
      B=EL*A(J)
      A(J)=A(I)-B
      A(I)=A(I)+B
7     CONTINUE
      EN=EQ*C2-EL
      X=REAL(EN)
      Y=AIMAG(EN)
      EN=EN*.5*(3.-X*X-Y*Y)
      EL=EQ
      EQ=EN
10    CONTINUE
      K=L
      E=0.5*E
      GO TO 6
8     CONTINUE
      RETURN
      END
