      parameter (mpt=1000,mpar=13)
      integer l_lpt
      character*140 lptfile
      COMMON/Files/l_lpt,lptfile
      COMMON/DATCOM/ X(mpt),Y(mpt),E(mpt),NDAT,MORDER,FIT(mpt)
