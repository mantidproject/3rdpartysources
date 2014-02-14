C
C     ************************************************************
C
      FUNCTION CHINOW(AX,M)
      REAL A(3,3),B(3),C1(3),C2(3,3),S1(3),S2(3,3),BETA(3)
      COMMON/SPACE/CHISQ,CHTARG,CHIZER,XSUM,BETA,C1,C2,S1,S2,BLANK
      BX=1.-AX
      DO 28 K=1,M
      DO 29 L=1,M
29    A(K,L)=BX*C2(K,L)-AX*S2(K,L)
28    B(K)=-(BX*C1(K)-AX*S1(K))
      CALL CHOSOL(A,B,M,BETA)
      W=0.
      DO 31 K=1,M
      Z=0.
      DO 32 L=1,M
32    Z=Z+C2(K,L)*BETA(L)
31    W=W+BETA(K)*(C1(K)+0.5*Z)
      CHINOW=1.+W
      RETURN
      END
