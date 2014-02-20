C
C     ************************************************************
C
      SUBROUTINE MAXENT(NGROUPS,NPTS,P,DATUM,
     +SIGMA,FLAT,BASE,MAX,SUMFIX)
      INTEGER P
      REAL F(68000),BASE(68000),OX(68000),CGRAD(68000),
     +SGRAD(68000),XI(68000,3),ETA(68000,3),S1(3),S2(3,3),
     +C1(3),C2(3,3),BETA(3),SIGMA(P),DATUM(P)
	 character*2048 str
      COMMON/SPACE/CHISQ,CHTARG,CHIZER,XSUM,BETA,C1,C2,S1,S2,BLANK
      common/savetime/ngo,i2pwr
      common/fac/factor,facdef,facfake,ratio
      common/heritage/iter
	common/MaxPage/n,f
      LOGICAL SUMFIX
c     
      BLANK=FLAT
      IF(BLANK.EQ.0.)GOTO 7
      DO 3 I=1,N
3     BASE(I)=BLANK
      GOTO 5
7     SUM=0.
      DO 4 I=1,N
4     SUM=SUM+BASE(I)
      BLANK=SUM/N
5     CHIZER=FLOAT(P)
      CHTARG=CHIZER
      ITER=0
      M=3
      if( ngo.gt.0 ) then
        iter=1
        goto 6
      endif
      DO 8 I=1,N
8     F(I)=BASE(I)
6     CALL OPUS(NGROUPS,NPTS,N,P,F,OX)
      CHISQ=0.
      DO 10 J=1,P
      A=OX(J)-DATUM(J)
      CHISQ=CHISQ+A*A/(SIGMA(J)*SIGMA(J))
10    OX(J)=2.*A/(SIGMA(J)*SIGMA(J))
      CALL TROPUS(NGROUPS,NPTS,N,P,OX,CGRAD)
      XSUM=0.
      TEST=0.
      SNORM=0.
      CNORM=0.
      TNORM=0.
      DO 12 I=1,N
      XSUM=XSUM+F(I)
      SGRAD(I)=-ALOG(F(I)/BASE(I))/BLANK
      SNORM=SNORM+SGRAD(I)*SGRAD(I)*F(I)
      CNORM=CNORM+CGRAD(I)*CGRAD(I)*F(I)
12    TNORM=TNORM+SGRAD(I)*CGRAD(I)*F(I)
      SNORM=SQRT(SNORM)
      CNORM=SQRT(CNORM)
      A=1.
      B=1./CNORM
      C=1./CNORM
      IF(ITER.NE.0)TEST=SQRT(0.5*abs(1.-TNORM/(SNORM*CNORM)))
      if(test.lt..0000001) test=.0000001
      IF(ITER.NE.0)A=1./(SNORM*2.*TEST)
      IF(ITER.NE.0)B=1./(CNORM*2.*TEST)
      DO 13 I=1,N
      XI(I,1)=F(I)*C*CGRAD(I)
13    XI(I,2)=F(I)*(A*SGRAD(I)-B*CGRAD(I))
      IF(SUMFIX)CALL PROJECT(1,N,XI)
      IF(SUMFIX)CALL PROJECT(2,N,XI)
      CALL OPUS(NGROUPS,NPTS,N,P,XI(1,1),ETA(1,1))
      CALL OPUS(NGROUPS,NPTS,N,P,XI(1,2),ETA(1,2))
      DO 14 J=1,P
14    OX(J)=ETA(J,2)/(SIGMA(J)*SIGMA(J))
      CALL TROPUS(NGROUPS,NPTS,N,P,OX,XI(1,3))
      A=0.
      DO 15 I=1,N
      B=F(I)*XI(I,3)
      A=A+B*XI(I,3)
15    XI(I,3)=B
      A=1./SQRT(A)
      DO 16 I=1,N
16    XI(I,3)=A*XI(I,3)
      IF(SUMFIX)CALL PROJECT(3,N,XI)
      CALL OPUS(NGROUPS,NPTS,N,P,XI(1,3),ETA(1,3))
      DO 17 K=1,M
      S1(K)=0.
      C1(K)=0.
      DO 18 I=1,N
      S1(K)=S1(K)+XI(I,K)*SGRAD(I)
18    C1(K)=C1(K)+XI(I,K)*CGRAD(I)
17    C1(K)=C1(K)/CHISQ
      DO 19 K=1,M
      DO 19 L=1,K
      S2(K,L)=0.
      C2(K,L)=0.
      DO 20 I=1,N
20    S2(K,L)=S2(K,L)-XI(I,K)*XI(I,L)/F(I)
      DO 21 J=1,P
21    C2(K,L)=C2(K,L)+ETA(J,K)*ETA(J,L)/(SIGMA(J)*SIGMA(J))
      S2(K,L)=S2(K,L)/BLANK
19    C2(K,L)=2.*C2(K,L)/CHISQ
      C2(1,2)=C2(2,1)
      C2(1,3)=C2(3,1)
      C2(2,3)=C2(3,2)
      S2(1,2)=S2(2,1)
      S2(1,3)=S2(3,1)
      S2(2,3)=S2(3,2)
      S=0.
      DO 22 I=1,N
22    S=S-F(I)*ALOG(F(I)/(BASE(I)*2.7182818285))/(BLANK*2.7182818285)
      A=S*BLANK*2.7182818285/XSUM
      
      write(str,103)ITER,TEST,S,CHTARG,CHISQ,XSUM
      call print_log_msg("notice", TRIM(str))

103   FORMAT(I3,4X,5(E10.4,2X))
      BETA(1)=-0.5*C1(1)/C2(1,1)
      BETA(2)=0.
      BETA(3)=0.
      IF(ITER.NE.0)CALL MOVE(3,sigma,p)
      A=0.
      DO 23 I=1,N
      DO 24 K=1,M
24    F(I)=F(I)+BETA(K)*XI(I,K)
      IF(F(I).LT.0.)F(I)=1.E-3*BLANK
      A=A+F(I)
23    CONTINUE
      IF(.NOT.SUMFIX)GOTO 50
      DO 51 I=1,N
51    F(I)=F(I)/A
50    ITER=ITER+1
      if(iter.le.1) goto 6
      IF(TEST.LT.0.02.AND.ABS(CHISQ/CHIZER-1.).LT.0.01.OR.ITER.GT.
     *MAX) THEN
c        write(99,98)
c        write(99,99)(F(I),I=1,N,4)
        RETURN
      ENDIF
99    FORMAT(1X,5E13.4)
98    FORMAT(1X,//,1X,'FREQUENCY SPECTRUM')
      GOTO 6
      END
