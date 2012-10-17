c  minus parameters
      parameter (mq=500,mw1=1000,mw2=2000)
      parameter (msurf=15,mreg=10)
      LOGICAL IN
      COMMON/TES/IN
      COMMON/DTO/ISURF,EXDIST
      COMMON/DIR/VX,VY,VZ
      COMMON/DEF/X,Y,Z,IREG
      COMMON/GEOM1/NSURF,NREG,IGEOM(msurf,mreg)
      COMMON/GEOM2/A(mreg),B(mreg),C(mreg),D(mreg),E(mreg),F(mreg),
     1 G(mreg),P(mreg),Q(mreg),R(mreg)
      COMMON/XSEC/DENS,SIGA,SIGB,SIGT,SIGT1,SIGS1
      COMMON/SofQW/NQ,DQ,NW,Nel,DW,en(mw1),Qq(mq),Ww(mw1),S(mq,mw1),
     1 HOVERT
      COMMON/Vector/VKINC,Vkm,Vkm2,VMFP,Vmu,VKX,VKY,VKZ,B1,B9
      COMMON/KEW/Qv,QMIN,Qran,QMAX,IIW,SQW
