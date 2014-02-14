c
C***********************************************************************
C           SUBROUTINE TO READ RUN DETAILS FROM OPENGENIE
C************************************************************************
C

      subroutine opengenie_maxent(pars_in,pars_out)
c     implicit none
        
      external pars_in,pars_out

	character*255 datafilename,readdata
      logical onepulse, varydt, fixphasein, retlog
      integer nptsin, zeroch, istgoodch, ierr, debug
      real deflevel, Sigmaloss

      INTEGER P,FIELD,ig
      CHARACTER*1 ANS,fitdead,fixphase,fitamp,firstgo
	common/MaxPage/n,f
      COMMON/ANS/ANS,firstgo
      COMMON/FASE/PHASE,phshift
      common/fac/factor,facfake,ratio
      common/channels/itzero,ipulse1,ipulse2,i1stgood,itotal
      common/points/npts,ngroups,nhists
      common/pulses/npulse,DEF
      COMMON/RUNDATA/RES,FRAMES,FNORM,IRUNNO,HISTS
      common/datall/rdata(64,4096)
      common/framin/iframnew
      common/datin/runin(64,4096),histco(64),hblock(15)
      integer rdata,runin
      character*80 hblock,histco
      COMMON/FLAGS/FITDEAD,FIXPHASE,FITAMP
      COMMON/FILE/NAME
      COMMON/GROUPING/GROUP(64)
      COMMON/MACH/MACHINE
	common/sense/phi(64),TAUD(64)
      common/savetime/ngo,i2pwr
      common/moments/bmindef
      CHARACTER*40 TITLE,NAME 
      character*3 machine
      character*1 mulrun
      CHARACTER*40 GRPTITLE
      CHARACTER*80 COMMENTLINE
	character*2046 str
      REAL HISTS(64)
      REAL NFRAM, F(68000),fchan(68000),fperchan
      INTEGER ITOTAL,NHISTOS,HISTO,DEST
      INTEGER GROUP, tempdata(256000)
      REAL RES


	integer histlen, nhisto

      call module_get_int(pars_in,'RunNo',Irunno)
      call module_get_int(pars_in,'frames',Ifram)
      call module_get_int(pars_in,'res',Ires)
	res=real(Ires)/1e6
c	write(99,*) res
c      call module_print("Hello")

      call module_get_int(pars_in,'Tzeroch',zeroch)
	call module_get_int(pars_in,'firstgoodch',istgood)
	call module_get_int(pars_in,'fitphase',fph)
	if (fph.eq.1) then
	  fixphasein=.true.
	else
	  fixphasein=.false.
	endif
	call module_get_int(pars_in,'fitdt',fdt)
	if (fph.eq.1) then
	  varydt=.true.
	else
	  varydt=.false.
	endif

	call module_get_real(pars_in,'deflevel',deflevel)
	onepulse=.true.

	call module_get_real(pars_in,'sigloose',Sigmaloss)

	call module_get_int(pars_in,'ptstofit',nptsin)

	call module_get_int(pars_in,'histolen',histlen)
	itotal=histlen
	call module_get_int(pars_in,'nhisto',nhisto)
	nhists=nhisto
	write(10,*) varydt, fixphasein, onepulse, deflevel,sigmaloss
     +,nptsin,histlen,nhisto

      call module_get_int_array
     +              (pars_in,'counts',tempdata,histlen*nhisto)
c	 write(69,*)tempdata
	do i=1,nhisto
	 do j=1,histlen
	   rdata(i,j)=tempdata((i-1)*histlen+j)
	 enddo
	enddo










      
c  default values of various consts collected for easy changing
      if (firstgo.ne.'T') goto 998
c for second time round, default value is prev. chosen
      irundef=-1
      npulsdef=1
      fitdead='Y'
      fixphase='Y'
      fitamp='Y'
      mulrun='N'
      defdef=0.1
      nptsdef=4096
      bmindef=25.
c itzero and i1stgood for 1/2 pulses, res= 8/16 ns
      itzdef18=24
      itzdef28=44
      i1stdef8=80
      itzdef116=39
      itzdef216=17
      i1stdef16=50
      facdef=1.04
998   mulrun='N'
997   write(99,*) ' Please enter -1 for currant bun ...'
!      pause
      write(99,'('' Run number ['',i5,'']:'',$)')irundef
!	pause
!      read(5,5,end=1) irunno
!      read(5,5) irunno
!5     format(i5)
      datafilename=""
      write(99,*) datafilename
!1     rewind(5)
      if(irunno.eq.0)irunno=irundef
c re-read the same data
      if(irunno.lt.0)irunno=0
c -ve means read the currant bun, which goes to multibunrd as 0
      ihaveaproblem=0
c      call multibunrd(datafilename,ihaveaproblem)
      if (ihaveaproblem.gt.0) then
        write(99,*) ' error in reading file....'
        goto 997
      endif
      irundef=irunno  
      firstgo='F'
      write(99,'( '' Resolution: '',f8.4)')res
      FRAMES=IFRAM
      write(69,*) frames
994   write(99,993) mulrun
993   FORMAT(1X,'Do you want to add in more runs? [',a1,']:',$)
      ans=' '
!      read(5,992) ans
      ans='N'
992   FORMAT(A)  
      if (ans.eq.'y'.or.ans.eq.'Y') then
        mulrun='Y'
        ans=' '
      endif
      if (ans.eq.' '.and.mulrun.eq.'Y') then
        mulrun='Y'
c        irundef=irundef+1
        goto 997 
      else
        mulrun='N'
      endif
      mulrun='N'
      write(99,'('' Single (1) or double (2) pulse? ['',
     + i1,''] :'',$)')npulsdef
C      read(5,82,end=80) npulse
!       read(5,82) npulse
      IF (onepulse) npulse=1
	IF (.not.onepulse) npulse=2

C80    rewind(5)
82    format(i)
      if(npulse.lt.1.or.npulse.gt.2) npulse=npulsdef
      npulsdef=npulse
	write(99,*) npulse
      write(99,4)fitdead 
4     FORMAT(1X,'Fit dead time correction ? [',a1,']:',$)
      IF (varydt) ans='Y'
	write(99,*) ans
!      read(5,3)ans
      IF (ans.eq.'Y'.or.ans.eq.'y') then
        fitdead='Y'
      ELSEIF (ans.ne.' ') then 
        fitdead='N'    
      ENDIF
101   write(99,2)fixphase      
2     FORMAT(1X,'Do you want to keep the phases fixed ['
     + ,a1,']:',$)
      IF (fixphasein) then 
	 ans='Y'
	else
	 ans='N'
	endif
	write(99,*) ans
!      read(5,3) ans
3     FORMAT(A)  
      if (ans.eq.'y'.or.ans.eq.'Y') then
          fixphase='Y' 
      elseif (ans.ne.' ') then
          fixphase='N'
      endif
      write(99,'('' What default level? ['',f5.2,''] :'',$)')defdef
      def=deflevel
	write(99,*) deflevel
!      read(5,122) def
C120   rewind(5)
122    format(g6.3)
      if(def.le.0.0) def=defdef
      defdef=def
125    write(99,'('' Points to fit? ['',i4,''] :'',$)')nptsdef
      npts=nptsin
	write(99,*) npts
!      read(5,123) npts
C124    rewind(5)
      if(npts.gt.8192) then
        write(99,*) ' Max points at present is 8192'
        goto 125
      endif
      if(npts.eq.0) npts=nptsdef
123    format(I4)
      tn=log(float(npts))/log(2.)
	write(99,*) tn
      if(abs(tn-float(int(tn+.5))).gt.1.e-5) goto 125
      i2pwr=int(tn+.5)
      nptsdef=npts
	write(99,*) res
      n=1024*npts/(2048*int(.016/res+.5))
      write(99,*)' No of frequency spectrum points chosen as:',n
c this'll give approx 1000 gauss max, for any res or npts
      if(n.gt.npts) n=npts
      if(n.lt.256) n=256
c will do for now till we make it like the psi program
      IF(NPULSE.EQ.2) THEN
c assuming res= .008 < .010 or .016 > .010
        if (res.lt..010) then
          write(99,'('' Enter zero-time channel ['',
     + i2,'']:'',$)')itzdef28
          itzero=zeroch
	    write(99,*) itzero
!          read(5,20) itzero
c231       rewind(5)
          IF(ITZERO.EQ.0) ITZERO=itzdef28
          itzdef28=itzero
        else
          write(99,'('' Enter zero-time channel ['',
     + i2,'']:'',$)')itzdef216
          itzero=zeroch
	    write(99,*) itzero
!          read(5,20) itzero
c232       rewind(5)
          IF(ITZERO.EQ.0) ITZERO=itzdef216
          itzdef216=itzero
        endif
      ELSEIF(NPULSE.EQ.1) THEN
        if (res.lt..010) then
          write(99,'('' Enter zero-time channel ['',
     + i2,'']:'',$)')itzdef18
          itzero=zeroch
	    write(99,*) itzero
!          read(5,20) itzero
C233       rewind(5)
          IF(ITZERO.EQ.0) ITZERO=itzdef18
          itzdef18=itzero
        else
          write(99,'('' Enter zero-time channel ['',
     + i2,'']:'',$)')itzdef116
          itzero=zeroch
	    write(99,*) itzero
!          read(5,20) itzero
C234       rewind(5)
          IF(ITZERO.EQ.0) ITZERO=itzdef116
          itzdef116=itzero
        endif
      ENDIF 
      if (res.lt..010) then
         write(99,'('' Enter first good channel ['',
     + i3,'']:'',$)')i1stdef8
      i1stgood=istgood
	write(99,*) "hello"
	write(99,*) i1stgood
!         read(5,20) i1stgood
C241      rewind(5)
         IF(i1stgood.EQ.0) i1stgood=i1stdef8
         i1stdef8=i1stgood
      else
         write(99,'('' Enter first good channel ['',
     + i3,'']:'',$)')i1stdef16
      i1stgood=istgood
	write(99,*) "goodbye"
	write(99,*) i1stgood
!         read(5,20) i1stgood
C242      rewind(5)
         IF(i1stgood.EQ.0) i1stgood=i1stdef16
         i1stdef16=i1stgood
      endif
      write(99,'('' Enter sigma looseness factor  ['',
     + f4.2,'']'',$)')facdef
      factor=sigmaloss
	write(99,*) factor
!      read(5,21) factor
C26    rewind(5)
      IF(FACTOR.LE.0.0) FACTOR=facdef
      facdef=factor
20    FORMAT(I2)
21    FORMAT(F4.2)
22    FORMAT(F5.1)
!      CLOSE(5)
C       LOADS IN GROUPING 
C       NGROUPS MAY CHANGE AT THIS POINT
      OPEN(20,FILE='STDGRP.PAR',STATUS='OLD')
      READ(20,300) GRPTITLE
      write(99,*) ' Heading of grouping table read in: ',GRPTITLE
300   FORMAT(1X,A40)
      READ(20,301) COMMENTLINE
301   FORMAT(1X,A80)
      READ(20,302) NHISTOS,NGROUPS
302   FORMAT(1X,I2,1X,I2)
      READ(20,301) COMMENTLINE
      READ(20,301) COMMENTLINE
      DO 900 J=1,nhistos
      READ(20,302) HISTO,GROUP(J)
900   CONTINUE
      close(20)
      DO 945 J=1,NGROUPS
      HISTS(J)=0
945   CONTINUE
C      This allows for the removal of a whole group, but...
c      scaling will be buggered if all groups are not the same size.
      DO 55 I=1,64
        TAUD(I)=0
55    CONTINUE
C       LOADS IN DEADTIMES
      OPEN(16,FILE='TAUD.DAT',STATUS='OLD')  
      DO 50 I=1,NGROUPS
      READ(16,*) TAUD(I)
50    CONTINUE
      CLOSE(16)
      write(str,*) ' Initial values of dead times'
c	call module_print(TRIM(str))
      write(str,*) (taud(i),i=1,ngroups)
c	call module_print(TRIM(str))
 	tcounts=0
	write(69,*) nhists, histlen
      do ig=1,NHISTS
	  sum=0 !integer sum of all data in this group
	  do ic=1,histlen
	    sum=sum+int(rdata(ig,ic))
		if(tcounts.lt.2000000000)tcounts=tcounts+rdata(ig,ic)

	  enddo
	  write(histco(ig),1978) sum

1978  FORMAT('End: ',i10)
      enddo
	write(69,*) histco
	write(69,*) "tcounts",tcounts,sum
      call multimax
	fperchan=1./(RES*float(npts)*2.)
      do  i=1,n
       fchan(i)=float(i-1)*fperchan/135.5e-4
	enddo


      call module_put_real_array(
     +    pars_out,"f",f,n)
	 call module_put_real_array(
     +    pars_out,"fchan",fchan,n)
	  call module_put_real_array(
     +    pars_out,"taud",taud,ngroups)
	 call module_put_real_array(
     +    pars_out,"phi",phi,ngroups)




      RETURN
      END
