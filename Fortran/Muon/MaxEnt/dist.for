C
C     ************************************************************
C
      FUNCTION DIST(M)
      REAL BETA(3),S1(3),S2(3,3),C1(3),C2(3,3)
      COMMON/SPACE/CHISQ,CHTARG,CHIZER,XSUM,BETA,C1,C2,S1,S2,BLANK
      W=0.
      DO 26 K=1,M
      Z=0.
      DO 27 L=1,M
27    Z=Z-S2(K,L)*BETA(L)
26    W=W+BETA(K)*Z
      DIST=W
      RETURN
      END
