c   Version of DISCUS  for use with Quasi-elastic scattering
c   Modified by WSH         Version March 2012
c
      SUBROUTINE MUSCAT_data(IDET,sfile,l_sf,rfile,l_rf,rinstr,nran,
     1 ijeom,rgeom,sam,ims,dqw,Q_in,S_in,
     2 kill,totals,iw,energy,scat1,scatm,RR,S_out)
      include 'muscat_com.f90'
      parameter (pi=3.14159265,d2r=pi/180.,E0=25.2429,AH=0.4818025)
      real*4 QSSUM(10),TOT(mw2,6)
      real*4 RE(mw2)
      character*140 filelpt
c
      integer IDET
cf2py intent(in) :: IDET                     !detecter/angle number
      integer l_sf, l_rf
cf2py intent(in) :: l_sf, l_rf                !length of filenames
      character*140 sfile, rfile
cf2py intent(in) :: sfile, rfile              !sample & spe filenames
      real rinstr(3)
cf2py intent(in) :: rinstr                     !instrument real parameters
      integer nran(4)
cf2py intent(in) :: nran                      !random no. parameters
      integer ijeom(2)
cf2py intent(in) :: ijeom                      !geometry integer parameters
      real rgeom(3)
cf2py intent(in) :: rgeom                     !geometry real parameters
      real sam(4)
cf2py intent(in) :: sam                       !sample parameters
      integer ims(5)
cf2py intent(in) :: ims                      !mscat integer parameters
      real dqw(2)
cf2py intent(in) :: dqw                       !dQ,dw parameters
      real Q_in(mq)
cf2py intent(in) :: Q_in                       !Q values
      real S_in(mq,mw1)
cf2py intent(in) :: S_in                     !detecter/angle number
      integer kill, iw
cf2py intent(out) :: kill, iw                 !error code, no.energies
      real TOTALS(5)
cf2py intent(out) :: totals                   !total mscat values
      real energy(mw2),scat1(mw2),scatm(mw2),RR(mw2)
cf2py intent(out) :: energy,scat1,scatm,RR    !Scat values
      real S_out(mw2)
cf2py intent(out) :: S_out                   !S(Q,w) values

      PKHT=1.0
      PKSH=0.0
      call INIT_SQW
      if(l_sf.gt.0)then
c       filelpt(1:l_sf)=sfile(1:l_sf)           !LPT o/p
       filelpt(1:l_sf)='G:/Science/Mantid/IRIS/Muscat_d.lpt'
      endif
      l_sf=0
c     rinstr = [efixed, theta, alfa]
      E2=rinstr(1)
      theta=rinstr(2)
      alfa=rinstr(3)
      if(E2.lt.0.00001)then
       kill=4  
       return
      else
       wave=1.8*SQRT(E0/E2)
      endif
c     nran = [NRUN1, NRUN2, JRAND, MRAND]
      NRUN1=nran(1)                        !no. of ntns 1
      NRUN2=nran(2)                        !no. of ntns 2
      JRAND=nran(3)                        !seed 1
      MRAND=nran(4)                        !seed 2
c     ijeom = [jeom, Jann]
      JEOM=ijeom(1)
      Jann=ijeom(2)
      if(JEOM.eq.2)then
c     rgeom = [thick, width, height]
       THICK=rgeom(1)
       WIDTH=rgeom(2)
      endif
      if(JEOM.eq.3)then
c     rgeom = [rad1, rad2, height]
       WIDTH=rgeom(1)
       WIDTH2=rgeom(2)
      endif
      HEIGHT=rgeom(3)
c     sam = [temp, dens, siga, sigb]
      TEMP=sam(1)
      DENS=sam(2)
      SIGA=sam(3)
      SIGB=sam(4)
c     ims = [NMST, NQ, NW, KR1]
      NMST=ims(1)
      NQ=ims(2)
      NW=ims(3)
      Nel=ims(4)
      KR1=ims(5)
      if(nq.gt.mq)then
       kill=5  
       return
      endif   
      if(nw.gt.mw1)then
       kill=6  
       return
      endif   
c     dqw = [DQ, DW]
      DQ=dqw(1)
      DW=dqw(2)
      N2=1
      alpha=alfa -90.             !change definition of alpha
      Vkx=COS(d2r*alpha)          !calc direction cosines
      Vky=SIN(d2r*alpha)
      Vkz=0.                      !default vertical
      VKINC=2*pi/wave             !inc wave.vector
      WINC = (VKINC**2)/AH

      if(KR1.lt.0)then                         !full LPT o/p
       np_q=1                                  !all values of array
       np_w=1
       np_r=1
      else
       np_q=2                                  !restricted values
       np_w=5
       np_r=2
      endif
      np_r=1
      KR1=IABS(KR1)
C
      do numq=1,NQ
       Qq(numq)=Q_in(numq)
       do numw=1,NW
        S(numq,numw)=S_in(numq,numw)
       end do
      end do
      QMIN=0.25*DQ                !must be non-zero
      QMAX = Qq(NQ)             !max-Q input
      Qran=QMAX-QMIN             !Q-range for random generator
      do numw=1,NW
       S_out(numw)=S(IDET,numw)
      end do
c
      if(IDET.eq.1)then
       if(l_sf.gt.0)then
        call open_f(2,filelpt)
        write(2,104)
104     FORMAT(' MuScat Data - DISCUS VERSION FOR microEV ENERGIES :',
     1  '  MODIFIED March 2012  : WSH')
        write(2,108)sfile(1:l_sf)
108     FORMAT(/'  Runnumber : ',a)
        write(2,117) TEMP,DENS,SIGA,SIGB
117     FORMAT(/'   Temperature (K)           = ',F10.3
     2 '   Atom Density (10-24)cm-3 = ',F10.5
     3 /'   Absorption Xsect at inc.K = ',F10.3
     4 '   Scattering Cross Section = ',F10.3)
        write(2,159)alfa,alpha
159     format(/'   Sample angle : input (deg) =',f8.2,
     1 ' ; calc (deg) =',f8.2)
        call set_geom(JEOM,Jann,WIDTH,WIDTH2,HEIGHT,THICK,.true.)
        write(2,160)wave,VKINC
160     format('   Incident wavelength (A) =',f8.3,
     1 ' ; Momentum (A-1) =',F8.3,' ; Energy (meV) =',F8.3)
        write(2,113) VKX,VKY,VKZ
113     FORMAT('   Incident Beam Direction Cosines :',3F10.5)
        write(2,1029)rfile(1:l_rf)
1029    format(/'  S(Q,w) from file : 'a)
        write(2,118) NQ,DQ,Qq(1),NW,DW
118     FORMAT(/'  Q MESH :',I5,' points at interval (DQ) =',f8.3,
     1 ' starting at ',F8.3,
     2 /'  W MESH :',I5,' points at interval (DW) =',f8.3,' meV')
        write(2,1105)
1105    FORMAT(/'    S(Q,W)  Array  in (meV)-1 '/)
        do I=1,NQ,np_q              !print S(Q,w)
         write(2,1106)I,Qq(I)
1106     FORMAT(' Q (',i3,') = ',f8.3)
         write(2,1107)(S(I,J),J=1,NW,np_w)
1107     FORMAT(1X,10e11.3)
        end do
        write(2,115) NRUN1,NRUN2
115     FORMAT(/'   Number of Incident Neutrons :',2I8)
        close(unit=2)
       else
        call set_geom(JEOM,Jann,WIDTH,WIDTH2,HEIGHT,THICK,.false.)
       endif
      endif
C
C     VKINC IS THE INCIDENT NEUTRON MOMENTUM IN A-1
C     THICK IS THE THICKNESS OF THE SAMPLE
C     DENS = NO. OF ATOMS A-3
C     SIGA = ABSORPTION XSECT AT INCIDENT K
C     NW = NO. OF W USED IN THE Q,W, MESH
C     NQ = NO. OF Q  DITTO
C     DQ,DW = STEP SIZE IN THE Q,W MESH
C     W = ENERGY GAINED BY THE NEUTRON
C     HOVERT = PLANK'S H OVER KB*TEMP IN MEV-1 UNITS
C     AH = 2*M/HCROSS IN A-2.MEV-1
      HOVERT = 11.6087/TEMP
      call LOG_SIG
      DO I=1,10
       QSSUM(I) = 1.0
      end do
      VKM=VKINC
      CALL NEWV
      if(JEOM.LT.3)then           !if plate sample
       X=0.
       Y=0.
       Z=0.
       Vx=Vkx
       Vy=Vky
       Vz=Vkz
       isurf=0
       call DTOEX
       DL=exdist
       TRANS=EXP(-Vmu*DL)         !calc transmission
       if(IDET.eq.1)then
        if(l_sf.gt.0)then
         OPEN(UNIT=2,FILE=filelpt,STATUS='old',FORM='formatted',
     1   access='append')
         write(2,139)TRANS
139      FORMAT(/'   Transmission      =',F10.5)
         close(unit=2)
        end if
       end if
      end if
      SIGAC = SIGA
      call unit_vector(VKX,VKY,VKZ)       !mod of inc wave.vector
      IEL = mw1                           !centre point of range mw2=2*mw1

c      do IDET = 1,NDET            !start loop over dets
       SIGA = SIGAC              !keep value 'cos for NE=1 it's set to 0
       beta=alpha-theta          !angle 'tween sample & det
       DKX=COS(d2r*beta)       !calc detector direction cosines
       DKY=SIN(d2r*beta)
       DKZ=0.                  !default vertical
       Qel=2.0*VKinc*SIN(0.5*d2r*theta)
       if(l_sf.gt.0)then
        OPEN(UNIT=2,FILE=filelpt,STATUS='old',FORM='formatted',
     1  access='append')
        write(2,135) IDET
135     FORMAT(//5X,'DETECTOR',I5/)
        write(2,161)theta,Qel
161     format('   Scattering angle (deg) = ',f8.2,' ; Q elastic = ',
     1  f8.3)
        write(2,136) beta
136     FORMAT('   Angle between Incident Beam and Detector',
     1 ' Direction (deg) =',F8.2)
        write(2,114) DKX,DKY,DKZ
114     FORMAT('   Detector Direction Cosines :',3F10.5)
       endif
       CALL FA01DS(JRAND,MRAND)    !INITIALISE THE RANDOM NO. GENERATOR
       call unit_vector(DKX,DKY,DKZ)   !mod of scattered wave.vector
       do I=1,6                !initialise array
        do J=1,mw2
         TOT(J,I) = 0.0
        end do
       end do

       NMSTP1 = NMST + 1           !no. of scattering events
       do NE=1,NMSTP1         !start loop over NE
        SIGA = SIGAC
        IF(NE.EQ.1) SIGA = 0.0
        QSSUM(NE) = 0.0
        VKM = VKINC
        NOINC = 0
1000    VKM = VKINC                !NEW NEUTRONS START HERE
        UKX = VKX
        UKY = VKY
        UKZ = VKZ
        CALL RONE(NOINC,JEOM)

        if(NE.GT.2) then          !scattering events only
         NEM1 = NE-1
         do IMST= 2,NEM1      !loop over scattering events (IMST)
          B9 = B9*SIGS1
          call Q_dir(Ukx,Uky,Ukz,QSS)
          QSSUM(NE) = QSSUM(NE) + QSS
          CALL NEWV           !NOW HAVE NEW NEUTRON DIRECTION
          Vx=Ukx
          Vy=Uky
          Vz=Ukz
          isurf=0
          call DTOEX
          DL=exdist
          B4 = 1.0 -EXP(-Vmu*DL)
          VL = -(VMFP*ALOG(1.0 - FA01AS(1)*B4))
          B9 = B9*B4/SIGT
          call inc_xyz(VL)
         end do         !end loop over scattering events (IMST)
        end if          !end scattering events (NE>2)

        VKMC = VKM
        Vx=Dkx
        Vy=Dky
        Vz=Dkz
        isurf=0
        call DTOEX
        DL=exdist
        IF(NE.EQ.1) DL = 0.0
        do IIW=1,NW                   ! LOOP OVER ENERGY TRANSFERS
         W =(IIW-Nel)*DW               !energy range of S(Q,w)
         TKM2 = VKMC*VKMC + AH*W
         if(TKM2.GT.0.0)then
          VKM = SQRT(TKM2)
          CALL NEWV
          QX = VKM*DKX - VKMC*UKX
          QY = VKM*DKY - VKMC*UKY
          QZ = VKM*DKZ - VKMC*UKZ
          Qv = SQRT(QX*QX + QY*QY + QZ*QZ)
          IF(Qv.LE.Qq(NQ)) then
           CALL SINT
           WEIGHT = B9*EXP(-Vmu*DL)*SQW*VKM/VKMC
           WEIGHT = WEIGHT*EXP(-HOVERT*W)*SIGB/(4*PI)
           WF = (VKM**2)/AH
           IWF = NINT((WF-WINC)/DW)  +  IEL
           IF(IWF.GE.1.AND.IWF.LE.mw2)then
            TOT(IWF,NE) = TOT(IWF,NE) + WEIGHT
           endif
          endif
         endif
        end do
        if(NE.LE.2) then
         if(NOINC.le.NRUN1) goto 1000
        else
         if(NOINC.le.NRUN2) goto 1000
        endif
       end do                             !end loop over NE

       CALL FA01CS(JRAND,MRAND)
       if(l_sf.gt.0)then
        write(2,130) JRAND,MRAND
130     FORMAT('   Random Number Starting Values :',2I8)
        write(2,106)
106     FORMAT(/12X,'Energy ',7X,'J1*',7X,'J1',8X,'J2',8X,
     1 'J3',8X,'J4',8X,'J5',6X,'Total',6X,'Mult',8X,'R1',8X,'R1*'/
     2 12X,'(meV)',66X,'Scat'/)
       endif

       DO IA=1,4               !initialise array
        TOTALS(IA)=0.
       end do
       nw0=(NW-1)/2            !#energies each side of zero
       iwmin=IEL-nw0           !energy range of input Q,w
       iwmax=IEL+nw0
       iwtot=iwmax-iwmin+1
       iw=0
       do I=iwmin,iwmax             !start loop over nwt2
        W = (I-IEL)*DW             !energy transfer
        iw=iw+1
        energy(iw)=W
        D0 = TOT(I,1)/NRUN1
        D1 = TOT(I,2)/NRUN1
        D2 = TOT(I,3)/QSSUM(3)
        D3 = TOT(I,4)*4*NRUN2/(QSSUM(4)*QSSUM(4))
        D4 = TOT(I,5)*27*NRUN2*NRUN2/(QSSUM(5)**3)
        D5 = TOT(I,6)*16*16*NRUN2*NRUN2*NRUN2/(QSSUM(6)**4)
        DMUL=D2+D3+D4+D5           !mult scatt
        DSUM=D1+DMUL               !total scatt
        if(DSUM.LT.0.00001)Then       !check dsum>0
         R1STAR = 0.0
         R1 = 0.0
        else                  !for division
         R1STAR = D0/DSUM
         R1 = D1/DSUM
        endif
          numb1=I/20
          numb2=20*numb1
          if(I.eq.numb2)then
           if(l_sf.gt.0)then
            write(2,105)I,W,D0,D1,D2,D3,D4,D5,DSUM,DMUL,R1,R1STAR
105         FORMAT(1X,I5,F12.4,3X,1P8E10.3,2F10.4)
           endif
          endif
         TOTALS(1)=TOTALS(1)+D0
         TOTALS(2)=TOTALS(2)+D1
         TOTALS(3)=TOTALS(3)+D2
         TOTALS(4)=TOTALS(4)+D3
         TOTALS(5)=TOTALS(5)+DSUM
        IF(KR1.EQ.1)RR(iw)=R1
        IF(KR1.EQ.2)RR(iw)=R1STAR
        WT=1000/NRUN2
        RE(N)=(0.01*RR(iw))**2*WT
        scat1(iw)=D1
        scatm(iw)=DMUL
       end do                  !end loop over nwt2

       if(l_sf.gt.0)then
        write(2,155)
155     FORMAT(16X,30('-'),30X,10('-'))
        write(2,156)TOTALS
156     FORMAT(16X,1P3E10.3,30X,E10.3)
        close(unit=2)
       endif
c      end do              !end loop over detectors

      RETURN
      END
