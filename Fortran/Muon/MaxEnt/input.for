C
C ***************************************************************
C
      SUBROUTINE INPUT(NHISTS,NGROUPS,NPTS,P,DATUM,SIGMA,CORR,DATT)
      INTEGER P,NPTS,NGROUPS
      common/channels/itzero,ipulse1,ipulse2,i1stgood,itotal
      COMMON/FLAGS/FITDEAD,FIXPHASE,FITAMP
      COMMON/MISSCHANNELS/MM
      COMMON/RUNDATA/RES,FRAMES,FNORM,IRUNNO,HISTS
      COMMON/GROUPING/GROUP
      COMMON/MACH/MACHINE
      common/sense/phi(96),TAUD(96)
      common/datall/rdata(96,4096)
      common/datin/runin(96,4096),histco(96),hblock(15)
      integer rdata,runin
      character*80 hblock,histco
      COMMON/FILE/NAME
      common/fac/factor,facdef,facfake,ratio
      common/pulses/npulse,DEF
      CHARACTER*40 TITLE
      CHARACTER*20 NAME
      CHARACTER*1 ANS,fitdead,fixphase,fitamp
      CHARACTER*3 MACHINE
	character*2046 str
      REAL CORR(NPTS,NGROUPS),DATT(NPTS,NGROUPS),HISTS(96)
      REAL DATUM(NPTS,NGROUPS),SIGMA(NPTS,NGROUPS),TAUD,NFRAM
      INTEGER TEMP(9000),ITOTAL,NHISTOS,HISTO,IDEST
      INTEGER COLDSTORE(9000,96),GROUP(96)
c      OPEN(5,FILE=NAME,STATUS='OLD')
      do 110 i=1,15
c      read(5,11)title
c      read(hblock(i),11)title
110   continue
      TOTFCOUNTS=0.0
      TOTHCOUNTS=0.0
c	write(99,*) ngroups,nhists
      DO 101 J=1,NGROUPS
      DO 101 I=1,9000
      COLDSTORE(I,J)=0
101   CONTINUE
C       LOOP TO LOAD AND REGROUP DATA
      DO 1 J=1,NHISTS
c      READ(5,11) TITLE
11    FORMAT(A40)
      IDEST=GROUP(J)
c      READ(5,*) (TEMP(K),K=1,ITOTAL)
      do 345 k=1,itotal
345   temp(k)=rdata(j,k)
      IF (IDEST.EQ.0)GOTO 121
      DO 2 K=1,ITOTAL
      COLDSTORE(K,IDEST)=COLDSTORE(K,IDEST)+TEMP(K)
2     CONTINUE
      HISTS(IDEST)=HISTS(IDEST)+1.
c      write(99,*) ' histo no, idest, temp(100), coldstore(100,idest)'
c      write(99,*) j,idest,temp(100),coldstore(100,idest)
c121   READ(5,'(5X,F9.0)') HCOUNTS
C     write(69,*) histco(j)
121   READ(histco(j),'(5X,F9.0)') HCOUNTS
C     write(69,*) hcounts
c	write(69,*) i1stgood
      TOTHCOUNTS=TOTHCOUNTS+HCOUNTS
c      READ(5,11) TITLE
1     CONTINUE
      MM=0
      DO 944 J=1,NGROUPS
      IF (HISTS(J).EQ.0.)MM=MM+1
944   CONTINUE
C       LOOP TO CREATE 'DATUM(K,I)'

      call print_log_msg("notice", "Group, no. of hists...")

      DO 3 I=1,NGROUPS
      IF(HISTS(I).EQ.0.)GOTO 55
      DO 6 J=1,I1STGOOD-1
6     COLDSTORE(J,I)=0
      DO 5 K=1,NPTS
      DATUM(K,I)=COLDSTORE(K-1+ITZERO,I)
c      IF (DATUM(K,I).GT.1.E6)write(99,*)K,I,DATUM(K,I)
      TOTFCOUNTS=TOTFCOUNTS+DATUM(K,I)
5     CONTINUE
55    write(str,*) i,hists(i)

      call print_log_msg("notice", TRIM(str))
  
3     CONTINUE

      write(str,*) 'Total counts in all histos=',TOTHCOUNTS/1.e3
      call print_log_msg("notice", TRIM(str))

   	write(str,*) 'Good muons (used for normn. to 10 Mev)=',TOTFCOUNTS
      call print_log_msg("notice", TRIM(str))

      icount=0
C     NORMALISE TO 10 MEGAEVENTS
      FNORM=10000000./TOTFCOUNTS
      DO 4 J=1,NGROUPS
      DO 4 I=1,npts
      IF(HISTS(J).EQ.0.) THEN
        SIGMA(I,J)=1.E15
        GOTO 4
      ENDIF
      ifirst=i1stgood-itzero+1
      ilast=itotal-itzero+1
      IF (I.ge.ifirst.and.I.le.ilast) THEN
         SIGMA(I,J)=SQRT(DATUM(I,J)+2)*FNORM
         DATUM(I,J)=DATUM(I,J)*FNORM
         DATT(I,J)=DATUM(I,J)
C            DATT STORES ORIGINAL DATA (NORMALISED)
         icount=icount+1
         ELSE
         SIGMA(I,J)=1.e15
         datum(i,j)=0.0
         datt(i,j)=0.0
      ENDIF
      TAU=TAUD(J)/(RES*HISTS(J)*FRAMES*fnorm)
c       fnorm allows for scaling of data to 10 Mev
      if(datum(i,j)*tau.ge.1.0) then

      write(str,*) i,j,datum(i,j),taud(j)
      call print_log_msg("notice", TRIM(str))

	endif
      CORR(I,J)=TAU*DATUM(I,J)*DATUM(I,J)
      IF (DATUM(I,J).GT.0.5) THEN
         SIGMA(I,J)=SIGMA(I,J)*(1.+.5*CORR(I,J)/DATUM(I,J))
      ENDIF
      DATUM(I,J)=DATUM(I,J)+CORR(I,J)
4     CONTINUE
c      write(99,*) ' DATT,DATUM AND CORR AT 100th POINT for 4 GROUPS:'
c      write(99,*) (DATT(100,l),l=1,4)
c      write(99,*) (DATUM(100,l),l=1,4)
c      write(99,*) (corr(100,l),l=1,4)
      ratio=SQRT(float(icount)/FLOAT(npts*(NGROUPS-MM)))
c adjusts sigma for missing counts or groups
      facfake=factor*ratio
      do  i=1,npts
       do  j=1,NGROUPS
        sigma(i,j)=sigma(i,j)*facfake
        if(sigma(i,j).eq.0.0) then

         write(str,'('' ZERO SIG!!!'',i,i)') i,j
         call print_log_msg("notice", TRIM(str))

        endif
	 enddo
	enddo
c      CLOSE(5)
      RETURN
      END
