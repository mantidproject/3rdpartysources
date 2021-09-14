C**<main program>*******************************************************
C      PROGRAM QUASI_STEXP_2D
C     -----------------------
C
C-----------------------------------------------------------------------
C  This program tries to estimate the parameters of an elstic peak and 
C  a single "stretched exponential" pertaining to some relevant data.
C  It was modified from QUASI_LINES_2D, a 2-D Bayesian Quasi-elastic 
C  line-fitting program. The resolution function is given on a constant 
C  X-binning grid, and is assumed to be invariant. The background is 
C  assumed to be linear and the data should be of an intermediate Genie
C  binary format, such as that output by batch-mode ICON.
C-----------------------------------------------------------------------
C  Written by: D.S. Sivia, Rutherford Appleton Lab., Oxfordshire, England.
C-----------------------------------------------------------------------
C Initial release for trials by CJC, MAA, WSH.        (DSS: 27-FEB-1992)
C Changed file access from sequential to direct.      (DSS:  5-MAR-1992)
C Changed STEPSZ reduction from 0.85 to 0.6 .         (DSS:  7-MAR-1992)
C Can read in a detector normalisation file.          (DSS: 12-MAR-1992)
C Readonly for resol. file; silent-running modified.  (DSS: 14-MAY-1992)
C WBF array-size increased in subroutine INBUFF.      (DSS: 28-MAY-1992)
C Fixed bug in CLRBUF; problem for lots of detectors. (DSS: 29-MAY-1992)
C Hard-wired meV as input and ueV as output.          (DSS:  2-OCT-1992)
C Put in STEPSZ safety-valve in case Chisq stuck!     (DSS: 27-OCT-1992)
C Took out small and negative number safety-valve.    (DSS: 17-DEC-1992)
C Fixed bug in DATIN1 (no more large error-bars!).    (DSS: 27-SEP-1993)
C Modification of QUASI_LINES_2D.FOR .                (DSS: 26-APR-1994)
C Do "stretched-exp" quasi_lines properly, I hope!    (DSS:  9-JUN-1999)
C Corrected minor crashable bug!                      (DSS:  7-JUL-1999)
C-----------------------------------------------------------------------
C
      SUBROUTINE QLstexp(numb,x_in,y_in,e_in,reals,opft,
     1 XD_in,XB_in,YB_in,Wy_in,We_in,dtn,xsc,sfile,rfile,l_fn,
     2 nd_out,xout,yout,eout,yfit,yprob)
      INCLUDE 'mod_data.f90'
      INCLUDE 'mod_files.f90'
      INCLUDE 'options.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(10),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /WRKCOM/ WORK(m_d2,2)
      COMMON /STEXP/  BETEXP,XBETA(101),BETCHI(101),NBTOT
      COMMON/ModPars/NBIN,IMIN,IMAX,RSCL,BNORM
      real x_in(m_d), y_in(m_d), e_in(m_d)
cf2py intent(in) :: x_in, y_in, e_in                    !sample data
      real reals(4)
cf2py intent(in) :: reals                               !real parameters
      real XD_in(m_d), XB_in(m_d), YB_in(m_d)
cf2py intent(in) :: XD_in, XB_in, YB_in                 !sample xrange, res data (blur)
      real Wy_in(m_sp), We_in(m_sp), dtn(m_sp), xsc(m_sp)
cf2py intent(in) :: Wy_in, We_in, dtn, xsc              !fixed width data
      integer numb(9)
cf2py intent(in) :: numb                                !integer parameters
      integer opft(4)
cf2py intent(in) :: opft                                 !options parameters
      integer l_fn
cf2py intent(in) :: l_fn                                 !length of filenames
      character*140 sfile, rfile
cf2py intent(in) :: sfile, rfile                         !sample & res filenames
      integer nd_out
cf2py intent(out) :: nd_out                              !number of output points
      real xout(m_d), yout(m_d), eout(m_d)
cf2py intent(out) :: xout, yout, eout                    !data values
      real yfit(4*m_d)
cf2py intent(out) :: yfit                                !fit values
      real yprob(2)
cf2py intent(out) :: yprob                               !probability values
      integer l_title,l_user
      character*80 title,user
      real XBLR(m_d),YBLR(m_d)
      REAL GRAD(m_p),HESS(m_p,m_p),DPAR(m_p),COVAR(m_p,m_p)
      REAL SIGPAR(m_p),FITPSV(m_p),PARS(m_p),DPARS(m_p)
      REAL DTNORM(m_sp),XSCALE(m_sp)
      REAL PRBSV(4,m_sp),POUT(4,m_sp),
     1 PRMSV(7,4,m_sp),SIGSV(7,4,m_sp)
      INTEGER INDX(m_p)
      prog='s'
      CALL init_paras
      do n=1,m_p
       GRAD(n)=0.0
       DPAR(n)=0.0
       SIGPAR(n)=0.0
       FITPSV(n)=0.0
       do m=1,m_p
        HESS(n,m)=0.0
        COVAR(n,m)=0.0
       end do
      end do
      do n=1,m_sp
       do n1=1,4
        PRBSV(n1,n)=0.0
        POUT(n1,n)=0.0
        do n2=1,7
         PRMSV(n2,n1,n)=0.0
         SIGSV(n2,n1,n)=0.0
        end do
       end do
      end do
c     numb = [ngrp, nsp, ntc, Ndat, nbin, Imin, Imax, Nb, nrbin]
      NSPEC=numb(1)                                      !no. of groups
      ISP=numb(2)                                   !group number
      ISPEC=ISP
      ntc=numb(3)                                     !no. of points
      NDAT=numb(4)
      NBIN=numb(5)
      IMIN=numb(6)
      IMAX=numb(7)
      NB=numb(8)
      nrbin=numb(9)
c     reals = [efix, theta[isp], rscl, bnorm]
      efix=reals(1)
      theta(ISP)=reals(2)
      RSCL=reals(3)                                      
      BNORM=reals(4)
      do n=1,m_d
       xin(n)=x_in(n)
       yin(n)=y_in(n)
       ein(n)=e_in(n)
       XDAT(n)=XD_in(n)
       XBLR(n)=XB_in(n)
       YBLR(n)=YB_in(n)
      end do
      o_el=opft(1)
      o_bgd=opft(2)
      o_w1=opft(3)
      do n=1,m_sp
       DTNORM(n)=dtn(n)                      !DTNORM
       XSCALE(n)=xsc(n)                      !XSCALE
       yprob=0.0
      end do
      lptfile = ''
      fileout1 = ''
      fileout2 = ''
      fileout3 = ''
      l_lpt=l_fn+8
      lptfile(1:l_lpt)=sfile(1:l_fn)//'_Qse.lpt'
      l_file=l_fn+8
      fileout1(1:l_lpt)=sfile(1:l_fn)//'_Qse.qse'
      l_title=9
      title='<unknown>'
      l_user=9
      user='<unknown>'
      if(ISP.eq.1)then                          !print info	
        call open_f(53,lptfile)
        WRITE(53,1107)sfile
1107    format(' Sample file: ',a140)
        WRITE(53,1108)rfile
1108    format(' Resolution file: ',a140)
        write(53,1110)xin(imin),xin(imax)
1110    format(' Energy range: ',f8.3,' to ',f8.3,' meV')
        if (o_el.eq.0) write(53,1111)
        if (o_el.eq.1) write(53,1112)
        if (o_bgd.eq.2) write(53,1113)
        if (o_bgd.eq.1) write(53,1114)
        if (o_bgd.eq.0) write(53,1115)
1111    format(' Elastic option : NO peak')
1112    format(' Elastic option : WITH peak')
1113    format(' Background option : sloping')
1114    format(' Background option : flat')
1115    format(' Background option : zero')
       if(o_w1.eq.1)then
        write(53,1116)
       else
        write(53,1117)
       endif
1116   format(' Width option : fixed from file ')
1117   format(' Width option : free')
       close(unit=53)
       call open_f(1,fileout1)
       n=1
        WRITE(n,1107)sfile
        write(n,1121)title(1:l_title)
1121    format(' Title : ',a)
        write(n,1122)user(1:l_user)
1122    format(' User : ',a)
        write(n,120)NSPEC,NDAT,xin(imin),xin(imax)
 120    FORMAT(2X,2I10,2x,2f10.3,' ueV')
        WRITE(n,*)' -------------------------------------------------'
        WRITE(n,1108)rfile
        WRITE(n,*) '  1'
        WRITE(n,*)' -------------------------------------------------'
        close(unit=n)
      endif                                            !end print
c
      CALL BLRINT(XBLR,YBLR,NB,1.0,.FALSE.)
      NDAT=ntc-1
      CALL DPINIT
      CALL GDINIT
      NDAT=ntc-1
      CALL DATIN(ISP,DTNORM,IDUF)
      CALL BLRINT(XBLR,YBLR,NB,XSCALE(ISP),.TRUE.)
      CALL DPINIT
      CALL PRINIT(FITP,3,NFEW,1)
      CALL FileInit(1,ISP)
      CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
      GOTO 2
    1 CALL SEARCH(GRAD,HESS,DPAR,NFEW,INDX,COVAR,FITP)
    2 NPARMS=4+2*NFEW
      BETEXP=1.0
      CHIOLD=CBCHI(FITP)
      CALL VCOPY(FITP,FITPSV,NPARMS)
      STEPSZ=0.3
      IF (NFEW.GE.1) STEPSZ=STEPSZ/10.0
      IAGAIN=0
      CDIFMN=0.003
      do I=1,200
        CALL REFINE(GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,STEPSZ)
        CALL NEWEST(COVAR,GRAD,NPARMS,NFEW,DPAR,FITP)
        CNORM=CBCHI(FITP)
        IF (CNORM.LE.CHIOLD) THEN
          CHIDIF=(CHIOLD-CNORM)/CNORM
          IF (ABS(CHIDIF).LE.CDIFMN) THEN
            IF (IAGAIN.EQ.0) THEN
              CDIFMN=0.00005
              STEPSZ=0.15
              IAGAIN=1
            ELSE
	          GOTO 3
	        ENDIF
          ENDIF
          CHIOLD=CNORM
          CALL VCOPY(FITP,FITPSV,NPARMS)
        ELSE
          CALL VCOPY(FITPSV,FITP,NPARMS)
          STEPSZ=STEPSZ*0.6
          IF (STEPSZ.LT.1.0E-10) GOTO 3
        ENDIF
      end do
   3  CALL REFINE(GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,0.25)
      CNORM=CBCHI(FITP)
      CALL REFINE(GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,0.7)
      CALL ERRBAR(COVAR,NPARMS,SIGPAR)
      IF (NFEW.GT.0) THEN
        CALL VCOPY(FITP,PARS,NPARMS)
        PARS(NPARMS+1)=BETEXP
        CALL VCOPY(SIGPAR,DPARS,NPARMS)
        DPARS(NPARMS+1)=0.1
        CALL SIMOPT(PARS,DPARS,COVAR,NPARMS+1)
        CALL VCOPY(PARS,FITP,NPARMS+1)
        CALL ERRBAR(COVAR,NPARMS+1,SIGPAR)
        CNORM=CCHI(FITP)/FLOAT(NDAT)
        do I=1,NPARMS+1
          DPARS(I)=0.4*SIGPAR(I)
        end do
        CALL SIMOPT(PARS,DPARS,COVAR,NPARMS+1)
        CALL VCOPY(PARS,FITP,NPARMS+1)
        CALL ERRBAR(COVAR,NPARMS+1,SIGPAR)
        CNORM=CCHI(FITP)/FLOAT(NDAT)
      ENDIF
      CALL SEEFIT(SIGPAR,CNORM,PRMSV(1,NFEW+1,ISP),
     *              SIGSV(1,NFEW+1,ISP),NPARMS+1)
      CALL OUTPRM(FITP,COVAR,NPARMS+1,NFEW,CNORM)
      call PROBN(CNORM,numb(4),DETLOG,NFEW,PROBLG)
      noff=NDAT*NFEW
      do n=1,NDAT
       yfit(noff+n)=FIT(n)
      end do 
      NFEW=NFEW+1
      if (o_el.eq.0) then               !no peak
       FITP(3)=0.0
       FITPSV(3)=0.0
      endif
      IF (NFEW.LE.1) GOTO 1
      nd_out=NDAT
      do n=1,nd_out
       xout(n)=XDAT(n)
       yout(n)=DAT(n)
       if(SIG(n).lt.1.0e-10)then
        eout(n)=0.0
       else
        eout(n)=SQRT(2.0/SIG(n))
       endif
      end do
      CALL PRBOUT(PRBSV,2,ISP,POUT)
      do l=1,2
       yprob(l)=POUT(l,isp)
      end do
      RETURN
      END
