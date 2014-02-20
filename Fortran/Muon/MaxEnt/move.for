C
C     ************************************************************
C
      SUBROUTINE MOVE(M,sigma,p)
      REAL BETA(3),S1(3),C1(3),S2(3,3),C2(3,3),sigma(1)
	character*255 str
      COMMON/SPACE/CHISQ,CHTARG,CHIZER,XSUM,BETA,C1,C2,S1,S2,BLANK
      common/fac/factor,facdef,facfake,ratio
      integer p
      CMIN=CHINOW(0.,M)
      IF(CMIN*CHISQ.GT.CHIZER)CTARG=0.5*(1.+CMIN)
      IF(CMIN*CHISQ.LE.CHIZER)CTARG=CHIZER/CHISQ
      A1=0.
      A2=1.
      jtest=0
      F1=CMIN-CTARG
      F2=CHINOW(1.,M)-CTARG
41    ANEW=0.5*(A1+A2)
      FX=CHINOW(ANEW,M)-CTARG
      IF(F1*FX.GT.0.)A1=ANEW
      IF(F1*FX.GT.0.)F1=FX
      IF(F2*FX.GT.0.)A2=ANEW
      IF(F2*FX.GT.0.)F2=FX
      IF(ABS(FX).GE.1.E-3) then
        jtest=jtest+1
        if(jtest.gt.10000) then

          write(str,*) ' stuck in MOVE : chi**2 not tight enough'
          call print_log_msg("notice", TRIM(str))

          do 222 i=1,p
222       sigma(i)=sigma(i)*.99
          factor=factor*.99
          facdef=factor
          facfake=facfake*.99
          
          write(str,333) factor
333       format(' tightening looseness factor by 1 % to:',f5.3)
          call print_log_msg("notice", TRIM(str))

          goto 444
        endif
        GOTO 41
      endif
444      W=DIST(M)
      IF(W.LE.0.1*XSUM/BLANK)GOTO 42
      DO 44 K=1,M
44    BETA(K)=BETA(K)*SQRT(0.1*XSUM/(BLANK*W))
42    CHTARG=CTARG*CHISQ
      RETURN
      END
