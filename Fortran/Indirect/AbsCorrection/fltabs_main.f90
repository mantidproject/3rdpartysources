C CALCULATE FLAT PLATE ABSORPTION FACTORS
C
c  Input parameters :
c  sampars - list of sample scattering & absorption cross-sections and density
c  canpars - list of can scattering & absorption cross-sections and density
c  ncan_in - =0 no can, >1 with can
c  thick - list of thicknesses ts,t1,t2
c  angles - list of angles
c  waves - list of wavelengths
c  Output parameters :
c  kill - =0 if succesful, else error code
c  A1 - Ass
c  A2 - Assc
c  A3 - Acsc
c  A4 - Acc
c
      SUBROUTINE FLTABS(ncan_in, thick, density, sigs, siga,
     1 angles, waves, n_sp, fname, l_fn,
     2 kill, A1, A2, A3, A4)
      parameter (mlam=10)
      real density(3),sigs(3),siga(3),thick(3)
cf2py intent(in) :: density, sigs, siga, thick
      real angles(2)
cf2py intent(in) :: angles
      real waves(mlam)
cf2py intent(in) :: waves
      integer ncan_in, n_sp, l_fn
cf2py intent(in) :: ncan_in, n_sp, l_fn
      character*120 fname
cf2py intent(in) :: fname
      integer kill
cf2py intent(out) :: kill
      real A1(mlam),A2(mlam),A3(mlam),A4(mlam)
cf2py intent(out) :: A1, A2, A3, A4
      real AMUS1(mlam),AMUC1(mlam)
      integer l_lpt
      character*120 lptfile
C
      PI=3.141592653
      PICONV=PI/180.
      do n=1,mlam
       AMUS1(n)=0.0
       AMUC1(n)=0.0
       A1(n)=1.0
       A2(n)=1.0
       A3(n)=1.0
       A4(n)=1.0
      end do
      ssigs=sigs(1)                            !sam scatt x-sect
      ssiga=siga(1)                             !sam abs x-sect
      rhos=density(1)                           !sam density
      ncan=ncan_in                              !no.of sam + cans
      ts=thick(1)                               !sam thicknes
      t1=thick(2)                               !can thickness 1
      t2=thick(3)                               !can thickness 2
      csigs=sigs(2)                              !can scatt x-sect
      csiga=siga(2)                         !can abs x-sect
      rhoc=density(2)                           !can density
      l_fn=0
      if(n_sp.eq.1)then                     !if first spectrum, lpt output
       if(l_fn.gt.0)then                    !if len>0 output
c      l_lpt=l_fn+8
C      lptfile(1:l_lpt)=fname(1:l_fn)//'_flt.lpt'
        lptfile='G:/Science/Mantid/IRIS/Abs_flt.lpt'
        call open_f(6,'G:/Science/Mantid/IRIS/Abs_flt.lpt')
        write(6,1001)ssiga,ssigs
1001    format(' Sample cross-sections : abs = ',f10.5,' ; scatt = ',
     1  f10.5)
        write(6,1002)ts,rhos
1002    format(' Sample thickness = ',f10.5,' ; density = ',f10.5)
        write(6,1003)angles(1)
1003    format(' Can angle = ',f10.5)
        if(ncan.gt.1)then
         write(6,1004)csiga,csigs
1004     format(' Can cross-sections : abs = ',f10.5,' ; scatt = ',
     1   f10.5)
         write(6,1005)t1,rhoc
1005     format(' Can thickness = ',f10.5,' ; density = ',f10.5)
        endif
        close(unit=6)
       endif                                 !end output
      endif                                  !end 1st spectrum
c
      if(l_fn.gt.0)then
       OPEN(UNIT=6,FILE='G:/Science/Mantid/IRIS/Abs_flt.lpt',
     1  STATUS='old',FORM='formatted',access='append')
       write(6,1006)n_sp,angles(2)
1006   format(' Spectrum = ',i4,' ; detector angle = ',f10.5)
       close(unit=6)
      endif
      TCAN1=angles(1)                     !angle can to beam
      TCAN=TCAN1*PICONV
      THETA1=angles(2)                          !THETAB value - detector angle
      THETA=PICONV*THETA1
c read sample cross sections
      do n=1,mlam
       AMUS1(n)=ssigs + ssiga*waves(n)/1.8
       AMUS1(n)=AMUS1(n)*rhos
      end do
      if(ncan.gt.1) then
       do n=1,mlam
        AMUC1(n)=csigs + csiga*waves(n)/1.8
        AMUC1(n)=AMUC1(n)*rhoc
       end do

      else
       rhoc=0.
      endif

C CALCULATE CORRECTIONS AT EACH WAVELENGTH
C
      SEC1=1./COS(TCAN)
      TSEC=THETA1-TCAN1
C
C TSEC IS THE ANGLE THE SCATTERED BEAM MAKES WITH THE NORMAL TO THE SAMPLE
C SURFACE.  IF ABS(TSEC) IS CLOSE TO 90 DEG. CALCULATION OF ABSORPTION
C COEFFICIENTS IS UNRELIABLE
C
      if(ABS(ABS(TSEC)-90.0).LT.1.0)then
c case where TSEC is close to 90
       ass=1.0
       do K=1,mlam                            !start loop over wavelengths
        A1(K)=ASS
        A2(K)=ASS
        A3(K)=ASS
        A4(K)=ASS
       end do
      else
       TSEC=TSEC*PICONV
       SEC2=1./COS(TSEC)
       do K=1,mlam                            !start loop over wavelengths
        AMUS=AMUS1(K)
        FS=F(AMUS,TS,SEC1,SEC2)
        ES1=AMUS*TS*SEC1
        ES2=AMUS*TS*SEC2
        if(ncan.gt.1)then
         AMUC=AMUC1(K)
         F1=F(AMUC,T1,SEC1,SEC2)
         F2=F(AMUC,T2,SEC1,SEC2)
         E11=AMUC*T1*SEC1
         E12=AMUC*T1*SEC2
         E21=AMUC*T2*SEC1
         E22=AMUC*T2*SEC2
        endif
        if(SEC2.LT.0.)then
         ASS=FS/TS
         if(ncan.gt.1)then
          ASSC=ASS*EXP(-(E11-E12))
          ACC1=F1
          ACC2=F2*EXP(-(E11-E12))
          ACSC1=ACC1
          ACSC2=ACC2*EXP(-(ES1-ES2))
         else
          ASSC=1.0
          ACSC=1.0
          ACC=1.0
         endif
        else
         ASS=EXP(-ES2)*FS/TS
         if(ncan.gt.1)then
          ASSC=EXP(-(E11+E22))*ASS
          ACC1=EXP(-(E12+E22))*F1
          ACC2=EXP(-(E11+E22))*F2
          ACSC1=ACC1*EXP(-ES2)
          ACSC2=ACC2*EXP(-ES1)
         else
          ASSC=1.0
          ACSC=1.0
          ACC=1.0
         endif
        endif
        tsum=t1+t2
        if(tsum.gt.0.) then
         ACC=(ACC1+ACC2)/tsum
         ACSC=(ACSC1+ACSC2)/tsum
        else
         ACC=1.0
         ACSC=1.0
        endif
        A1(K)=ASS
        A2(K)=ASSC
        A3(K)=ACSC
        A4(K)=ACC
        if(l_fn.gt.0)then
         OPEN(UNIT=6,FILE='G:/Science/Mantid/IRIS/Abs_flt.lpt',
     1   STATUS='old',FORM='formatted',access='append')
         WRITE(6,1010) waves(K),ASS,ASSC,ACSC,ACC
1010     FORMAT(5F10.5)
         close(unit=6)
        endif
       end do                                       !end loop over wavelengths
      endif
      kill=0
      RETURN
      END
