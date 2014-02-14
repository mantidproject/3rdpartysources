C
C**************************************************************
C
      SUBROUTINE START(NGROUPS,NPTS,P)
      REAL AMP(64),PH(64),PHASE(64)
      common/pulses/npulse,DEF
      COMMON/FASE/PHASE,phshift
      COMMON/DETECT/A,B,E,c,d
      REAL A(64),B(64),c(64),d(64),E(8192),AA(8192),GW(8192),WW(8192)
      COMMON/MISSCHANNELS/MM
      COMMON/RUNDATA/RES,FRAMES,FNORM,IRUNNO,HISTS
      COMMON/MACH/MACHINE
      common/pulseshape/convolR(8192),CONVOLI(8192)
      common/savetime/ngo,i2pwr
      common/channels/itzero,ipulse1,ipulse2,i1stgood,itotal
	common/MaxPage/n,f
      real hpulsew,pulse2t,pi
	real f(68000)
      TAUlife=2.19704/RES
      DO 2 I=1,npts
2     E(I)=EXP(-float(I-1)/TAUlife)
c
c     create convolution factor for finite pulse width and 2 pulse
c
      hpulsew=.05
      pulse2t=.324
      IF(NPULSE.EQ.1) PULSE2T=0.0
      fperchan=1./(RES*float(npts)*2.)
      pi=3.1415926535
      TPION=0.026
      TMUON=2.19704
      do 4 i=1,N
C     
      ww(i)=2*pi*fperchan*float(i-1)
C            GW(I) is parabolic proton pulse ft GW(0)=1
      IF(I.NE.1)THEN
      GW(I)=(3/HPULSEW)*(SIN(WW(I)*HPULSEW)/(HPULSEW**2*WW(I)**3)-
     +COS(WW(I)*HPULSEW)/(HPULSEW*WW(I)**2))
      ELSE
      GW(I)=1.0
      ENDIF
      AA(I)=GW(I)/(1+(WW(I)*TPION)**2)
C
c      IF(NPULSE.EQ.1) AA(I)=-AA(I)
      convolR(i)=AA(I)*(cos(ww(i)*pulse2t/2)-
     +TANH(PULSE2T/(2*TMUON))*SIN(WW(I)*PULSE2T/2)*WW(I)*TPION)
C
      CONVOLI(I)=-AA(I)*(TANH(PULSE2T/(2*TMUON))*SIN(WW(I)*PULSE2T/2)+
     +COS(pulse2t*WW(I)/2)*WW(I)*TPION)
4     continue
      RETURN
3     END
