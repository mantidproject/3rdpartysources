C CALCULATE CYLINDER ABSORPTION FACTORS
C
c  Input parameters :
c  nprof_in - number of profile values
c  profile - list of profile values
c  astep_in - step size
c  beam - list of sample & beam parameters
c  nan_in -number of annuli
c  radii - list of radii (for each annulus)
c  density - list of densities (for each annulus)
c  sigs - list of scattering cross-sections (for each annulus)
c  siga - list of absorption cross-sections (for each annulus)
c  angle - list of angles
c  wavelas - elastic wavelength
c  waves - list of wavelengths
c  Output parameters :
c  kill - =0 if succesful, else error code
c  Abs1 - Ass
c  Abs2 - Assc
c  Abs3 - Acsc
c  Abs4 - Acc
c
      SUBROUTINE CYLABS(astep_in,beam,nan_in,radii,
     1 density,sigs,siga,angle,wavelas,waves,n_sp, fname, l_fn,
     2 kill, Abs1, Abs2, Abs3, Abs4)
      INCLUDE 'mod_files.inc'
      parameter (mlam=10,mann=3,mbank=20)
      real beam(9)
cf2py intent(in) :: beam
      real radii(mann+1)
cf2py intent(in) :: radii
      real density(mann),sigs(mann),siga(mann)
cf2py intent(in) :: density, sigs, siga
      real wavelas, angle, astep_in
cf2py intent(in) :: wavelas, angle, astep_in
      real waves(mlam)
cf2py intent(in) :: waves
      integer nan_in, n_sp, l_fn
cf2py intent(in) :: nan_in, n_sp, l_fn
      character*120 fname
cf2py intent(in) :: fname
      integer kill
cf2py intent(out) :: kill
      real Abs1(mlam),Abs2(mlam),Abs3(mlam),Abs4(mlam)
cf2py intent(out) :: Abs1, Abs2, Abs3, Abs4
C
      REAL abscs(mann),assc(mann)
      COMMON/PROFBL/PROFIL(50),NPROF,PRSTEP
      COMMON NAN,RAD(mann+1),THETA,PI,A,B,AMUS(mann),AMUT(mann),
     1 DEN(mann),amuse(mann),amute(mann),scatxs(mann)
      REAL MUS(mann,mlam),MUT(mann,mlam),SIGTL(mann,mlam),
     1 SIGSL(mann,mlam),sigse(mann),sigte(mann)
      INTEGER lpt

      PI=3.141592653
      PICONV=PI/180.
      l_lpt=33                              !this works
      lptfile(1:l_lpt)='G:/Science/Mantid/IRIS/cylabs.lpt'
c      l_lpt=35
c      lptfile(1:l_lpt)='G:/Science/Mantid/IRIS/irs26176.lpt'
c      l_lpt=l_fn+10
c      lptfile(1:l_lpt)=fname(1:l_fn)//'cylabs.lpt'
c      l_fn=l_lpt
      lpt=0

      NLAMB=mlam
      do J=1,NLAMB                                !loop over wavelengths
       Abs1(J)=1.0
       Abs2(J)=1.0
       Abs3(J)=1.0
       Abs4(J)=1.0
      end do
      nprof=2                     ! NO. OF PROFILE VALUES AND VALUES
      do n=1,nprof
       PROFIL(n)=1.0
      end do
      AM=0.                                  !NORMALIZE PROFILE VALUES
      do I=1,NPROF
       AM=AMAX1(AM,PROFIL(I))
      end do
      do I=1,NPROF
       PROFIL(I)=PROFIL(I)/AM
      end do
      ASTEP=astep_in                          !INTEGRATION PARAMETERS

C  READ HEIGHT OF SAMPLE AND BEAM PARAMETERS
      HEIGHT=beam(1)
      A=beam(2)
      B=beam(3)
      A1=beam(4)
      B1=beam(5)
      HDOWN=beam(6)
      HUP=beam(7)
      HSDOWN=beam(8)
      HSUP=beam(9)
C
C  SET UP EQUIVALENT BEAM - INTEGRATE THE PROFILE FUNCTION AND
C  GENERATE AN EQUIVALENT WIDTH TO FIND THE CENTRE OF THE
C  EQUIVALENT BEAM RELATIVE TO THE CENTRE OF THE SAMPLE
C
      PSUM=0.
      WSUM=0.
      PRSTEP=(A-B)/(NPROF-1)
      PSUM=PROFIL(1)/2.
      WSUM=A*PSUM
      do I=2,NPROF
       X=A-(I-1)*PRSTEP
       PSUM=PSUM+PROFIL(I)
       WSUM=WSUM+X*PROFIL(I)
      end do
      PSUM=PSUM-PROFIL(NPROF)*0.5
      WSUM=WSUM/PSUM
      WIDTH=PSUM*PRSTEP
C
C DEFINE VALUES A2,B2 TO GIVE WIDTH AND POSITION OF EQUIVALENT
C BEAM 0F EQUAL INTENSITY BUT RECTANGULAR PROFILE -WSUM
C GIVES THE POSITION OF CENTRE OF NEW BEAM AND WIDTH GIVES
C WIDTH OF THE BEAM
C
      A2=WSUM+WIDTH/2
      B2=WSUM-WIDTH/2
C
      NAN=nan_in                             !NO. OF ANNULI
      IF (NAN.GT.mann) NAN=mann              !PRESENT MAXIMUM IS mann
      NAN1=NAN+1                             !NAN+1 RAD VALUES
      do n=1,NAN1
       RAD(n)=radii(n)
      end do
      do n=1,NAN
       DEN(n)=density(n)                     !density
       ABSCS(n)=siga(n)                      !absorption Xsec
       scatxs(n)=sigs(n)                      !scattering Xsec
       sigte(n)=scatxs(n)+siga(n)*wavelas/1.7979
      end do
      if(n_sp.eq.1)then
       if(lpt.gt.0)then
        call open_f(6,lptfile)
        write(6,1111)fname(1:l_fn)
1111    format(' filename : ',a)
        write(6,1001)siga(1),sigs(1)
1001    format(' Sample cross-sections : abs = ',f10.5,' ; scatt = ',
     1 f10.5)
        write(6,1002)radii(1),radii(2)
1002    format(' Sample radius1 = ',f10.5,' ; radius2 = ',f10.5)
        write(6,1003)density(1)
1003    format(' Sample density = ',f10.5)
        if(nan.eq.2)then
         write(6,1004)siga(2),sigs(2)
1004     format(' Can cross-sections : abs = ',f10.5,' ; scatt = ',
     1   f10.5)
         write(6,1005)radii(3)
1005     format(' Sample radius3 = ',f10.5)
         write(6,1006)density(2)
1006     format(' Can density = ',f10.5)
        endif
        close(unit=6)
       endif
      endif
c
      if(lpt.gt.0)then
       OPEN(UNIT=6,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
       write(6,1009)n_sp,angle
1009   format(' Spectrum = ',i4,' ; Detector angle = ',f10.5)
      endif
      do I=1,NAN                              !loop over annuli
       do IL=1,NLAMB                          !CALC SCATTERING C/S AT EACH WAVELENGTH 
        SIGSL(I,IL)=scatxs(I)
        SIGTL(I,IL)=scatxs(I)+ABSCS(I)*waves(IL)/1.7979
        MUS(I,IL)=DEN(I)*SIGSL(I,IL)          !Mu for scattering
        MUT(I,IL)=DEN(I)*SIGTL(I,IL)          !Mu for total
       end do                                 !end loop IL
       sigse(i)=scatxs(i)
       amuse(i)=den(i)*sigse(i)
       amute(i)=den(i)*sigte(i)
      end do                                   !end loop I over annuli

      MS=(RAD(2)-RAD(1)+0.0001)/ASTEP
      IF (MS.LT.1)MS=1
C
C CALCULATE CORRECTIONS AT EACH WAVELENGTH
C
      THETA=angle*PICONV
      do J=1,NLAMB                                !loop over wavelengths
       do IR=1,NAN                                !loop over annuli
        AMUS(IR)=MUS(IR,J)
        AMUT(IR)=MUT(IR,J)
       end do                                     !end loop IR
       nanw=nan-1
       CALL ACYL(MS,ASS,assc,ACSC,ACC,AREAS,AREAC)
       Abs1(J)=ASS
       Abs2(J)=ASSC(1)
       Abs3(J)=ACSC
       Abs4(J)=ACC
       if(lpt.gt.0)WRITE(6,1010) waves(J),ASS,ASSC(1),ACSC,ACC
1010   FORMAT(5F10.5)
      end do                                       !end loop J
      if(lpt.gt.0)close(unit=6)
      kill=0
      RETURN
      END
