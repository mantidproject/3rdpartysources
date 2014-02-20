
C
C ************************************************************
C
      SUBROUTINE OUTSPEC(NGROUPS,NPTS,N,P,DATUM,F,SIGMA,datt)
      INTEGER P
      REAL DATUM(NPTS,NGROUPS),SIGMA(NPTS,NGROUPS),F(68000)
      REAL SF(9000),GUESS(9000,64),TEST(9000,64),chi(64)
      REAL A(64),B(64),c(64),d(64),E(8192)
      real datt(npts,ngroups)
      CHARACTER*1 ANS,fitdead,fixphase,fitamp
      CHARACTER*1 AA
	character*2046 str
      REAL ZR(8192),ZI(8192)  
      COMMON/FLAGS/FITDEAD,FIXPHASE,FITAMP
      COMMON/RUNDATA/RES,FRAMES,FNORM,IRUNNO,HISTS
      COMMON/MISSCHANNELS/MM
      common/channels/itzero,ipulse1,ipulse2,i1stgood,itotal
      common/pulseshape/convolR(8192),CONVOLI(8192)
      common/fac/factor,facdef,facfake,ratio
      COMMON/DETECT/A,B,E,c,d
      CALL ZFT(P,N,F,ZR,ZI)
      do 998 j=1,NGROUPS
      chi(j)=0.0
      ibads=0
      ilast=itotal-itzero
      if(npts.lt.ilast) ilast=npts
      do 999 i=1,ilast
      guess(i,j)=(A(J)*ZR(I)+B(J)*ZI(I))
      test(i,j)=(datum(i,j)-guess(i,j))*ratio/sigma(i,j)
      chi(j)=chi(j)+(test(i,j)/ratio)**2
      if (abs(test(i,j)).lt.5.0) goto 999
      ibads=ibads+1
      if (ibads.gt.10) goto 999

      write(str,555)i,j
      call print_log_msg("notice", TRIM(str))
555   format('Devn > 5 std devs on point:',i4,' of group:',i2)

      if (ibads.eq.10) then

        write(str,*) ' ... lots of baddies in group ',j,' ...'
        call print_log_msg("notice", TRIM(str))

	    endif
999   continue
      chi(j)=chi(j)/float(npts)
998   continue

      write(str,*) 'Contributions to chi**2 from each group:'
      call print_log_msg("notice", TRIM(str))

      write(str,*) (chi(J),J=1,NGROUPS)
      call print_log_msg("notice", TRIM(str))

      DO 99 I=npts/2,NPTS
       DO 99 J=1,NGROUPS
        IF (SIGMA(I,J).GT.1.E3) GOTO 99
        SUM=SUM+(GUESS(I,J)/SIGMA(I,J))**2
        NSUM=NSUM+1
99    CONTINUE
      if(nsum.eq.0) then
        
        write(str,*)'No data in last half, so last half sum = 0'
        call print_log_msg("notice", TRIM(str))

        goto 87
      endif
      BUM=SUM/NSUM

      write(str,*) 'Points:',npts
      call print_log_msg("notice", TRIM(str))

      write(str,*) 'Last half sum:',BUM
      call print_log_msg("notice", TRIM(str))

87    write(str,*) 'Run: ',irunno
      call print_log_msg("notice", TRIM(str))

!      write(99,'('' Do you want to see some time diagnostics? [n]'',$)')
!      read(5,84) aa
!84    format(a1)
!      if(aa.ne.'y'.and.aa.ne.'Y') goto 85
!      write(99,'('' What group to be outputted as a spectrum? [1]'',$)')
!      read(5,11) J
!11    FORMAT(I)
!C12    rewind(5)
!      IF(J.EQ.0)J=1
!      DO 4 K=1,ilast
!      SF(K)=guess(k,j)
!4     continue
!      OPEN(8,FILE='OUTSPEC',STATUS='UNKNOWN')    
!      WRITE(8,2)(SF(I),I=1,ilast)
!      close(8)
!2     FORMAT (E13.4)   
!      DO 10 I=1,ilast
!      SF(I)=datum(I,j)
!10    CONTINUE
!      OPEN(8,FILE='OUTSPEC',STATUS='UNKNOWN')    
!      WRITE(8,2)(SF(I),I=1,ilast)
!      CLOSE(8)
!      OPEN(8,FILE='OUTSPEC',STATUS='UNKNOWN')    
!      WRITE(8,2)(test(I,j),I=1,ilast)
!      CLOSE(8)
c      write(99,*) ' dodgey point guess, test, sigma and datum (1005,1) '
c      write(99,*) guess(1005,1),test(1005,1),sigma(1005,1),datum(1005,1)
c      write(99,*) ' dodgey point datt, exp, fnorm (1005,1) '
c      write(99,*) datt(1005,1),d(1)*e(1005),fnorm
c      write(99,*) ' good point guess, test, sigma and datum (1003,1) '
c      write(99,*) guess(1003,1),test(1003,1),sigma(1003,1),datum(1003 ,1)
c      write(99,*) ' good point datt, exp, fnorm (1003,1) '
c      write(99,*) datt(1003,1),d(1)*e(1003),fnorm
      NSUM=0
      SUM=0.0
85    RETURN
      END
