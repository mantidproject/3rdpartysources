C**********************************************************************
C* Subroutine to accept data/parameters from F2Py (Mantid), and put   *
C* it all into the common blocks. We're using three different         *
C* callbacks to pass in the parameters; this will save us from having *
C* to have a huge amount of arguments listed in the subroutine        *
C* definition.  This method is preferable to just passing in arrays of*
C* parameters since it's easy to get indices mixed up.  This way is   *
C* almost like dictionary access.                                     *
C**********************************************************************

      subroutine mantid_maxent(data_in, groups_in, taud_in, phase_in,
     +f_out, fchan_out, taud_out, phi_out)

C     Correct way to use STDOUT pipe.  For debug statements.
      use iso_fortran_env, only stdout => output_unit

C     Various array inputs.
      integer data_in(262144)
Cf2py intent(in) data_in
      integer groups_in(64)
Cf2py intent(in) groups_in
      real taud_in(64)
Cf2py intent(in) taud_in
      real phase_in(64)
Cf2py intent(in) phase_in

C     Declare python callbacks, which we use to look up input parameters of
C     various types.
Cf2py character(100) par_name
Cf2py intent(callback) int_par
      external int_par
      integer int_par
Cf2py integer int_value
Cf2py int_value = int_par(par_name)
Cf2py intent(callback) real_par
      external real_par
      real real_par
Cf2py real real_value
Cf2py real_value = real_par(par_name)
Cf2py intent(callback) bool_par
      external bool_par
      logical bool_par
Cf2py logical bool_value
Cf2py bool_value = bool_par(par_name)

C     Declare python callback for logging.  Arbitrarily large size for
C     log_message parameter.
Cf2py character(20)   log_priority
Cf2py character(1000) log_message
Cf2py intent(callback) logger
      external logger
      logical logger
      logical void
Cf2py void = logger(log_priority, log_message)

C     Results.
      real f_out(68000)
Cf2py intent(out) f_out
      real fchan_out(68000)
Cf2py intent(out) fchan_out
      real taud_out(64)
Cf2py intent(out) taud_out
      real phi_out(64)
Cf2py intent(out) phi_out

C     Same variables as declared in the original opengenie_maxent subroutine.
      character(255) readdata
      logical onepulse, varydt, fitdead, fixphase, retlog
      integer nptsin, zeroch, istgoodch, ierr, debug
      real deflevel, sigmaloss
      integer p,field,ig
      character(1) ans,fitamp,firstgo
      common/MaxPage/n,f
      common/ans/ans,firstgo
      common/fase/phase,phshift
      common/fac/factor,facfake,ratio
      common/channels/itzero,ipulse1,ipulse2,i1stgood,itotal
      common/points/npts,ngroups,nhists
      common/pulses/npulse,def
      common/rundata/res,frames,fnorm,irunno,hists
      common/datall/rdata(64,4096)
      common/framin/iframnew
      common/datin/runin(64,4096),histco(64),hblock(15)
      integer rdata,runin
      character(80) hblock,histco
      common/flags/fitdead,fixphase,fitamp
      common/file/name
      common/grouping/group(64)
      common/mach/machine
      common/sense/phi(64),taud(64),phases(64)
      common/savetime/ngo,i2pwr
      common/moments/bmindef
      character(40) title,name
      character(3) machine
      character(1) mulrun
      character(40) grptitle
      character(80) commentline
      character(2046) str
      real hists(64)
      real nfram,f(68000),fperchan
      integer itotal,nhistos,histo,dest,histlen
      integer group
      real res, num, demon

      call print_log_msg("debug", "Starting underlying Fortran.")

C     Assign input parameters to common blocks in roughly the same fashion as
C     the original opengenie_maxent subroutine.  I've removed all of the code
C     to do with reading from STDIN and tried to simplify things, but there
C     may be some redundancy still.

      irunno = int_par("RunNo")
      ifram = int_par("frames")
      zeroch = int_par("Tzeroch")
      istgood = int_par("firstgoodch")
      nptsin = int_par("ptstofit")
      histlen = int_par("histolen")
      ires = int_par("res")
      nhistos = int_par("nhisto")
      ngroups = int_par("n_groups")

      res = real(ires)/1e6

      deflevel = real_par("deflevel")
      Sigmaloss = real_par("sigloose")

C     In the OpenGENIE version, fixphase was (incorrectly?) used to set varydt.
C     I'm assuming this way is correct.
      fitdead = bool_par("fitdt")
      fixphase = bool_par("fixphase")

      onepulse=.true.


      itotal = histlen
      nhists = nhistos

      do i=1,nhistos
      do j=1,histlen
          rdata(i,j)=data_in((i-1)*histlen+j)
      enddo
      enddo
      
      fitamp='Y'
      mulrun='N'
      nptsdef=4096
      bmindef=25.
c     itzero and i1stgood for 1/2 pulses, res= 8/16 ns
      itzdef18=24
      itzdef28=44
      i1stdef8=80
      itzdef116=39
      itzdef216=17
      i1stdef16=50
      facdef=1.04

997   firstgo='f'
      frames=ifram

      ans='N'

C     npulsdef = single (1) or double (2) pulse
      if (onepulse) npulse=1
      if (.not.onepulse) npulse=2

      def=deflevel

      if(def.le.0.0) def=0.1

      npts=nptsin
      if(npts.eq.0) npts=nptsdef
      tn=log(float(npts))/log(2.)
      i2pwr=int(tn+.5)
      nptsdef=npts

C     Number of frequency spectrum points
      num = 1024*npts
      demon = 2048*int(.016/res+.5)

C      n=1024*npts/(2048*int(.016/res+.5))
      n = int(num / demon)

c     this'll give approx 1000 gauss max, for any res or npts
      if(n.gt.npts) n=npts
      if(n.lt.256) n=256
c     will do for now till we make it like the psi program

      if(npulse.eq.2) then         
c         c assuming res= .008 < .010 or .016 > .010
          if (res.lt..010) then
              itzero=zeroch
              if(itzero.eq.0) itzero=itzdef28
              itzdef28=itzero
          else
              itzero=zeroch
              if(itzero.eq.0) itzero=itzdef216
              itzdef216=itzero
          endif
      elseif(npulse.eq.1) then
          if (res.lt..010) then
              itzero=zeroch
              if(itzero.eq.0) itzero=itzdef18
              itzdef18=itzero
          else
              itzero=zeroch
              if(itzero.eq.0) itzero=itzdef116
              itzdef116=itzero
          endif
      endif 

      if (res.lt..010) then
          i1stgood=istgood
          if(i1stgood.eq.0) i1stgood=i1stdef8
          i1stdef8=i1stgood
      else
          i1stgood=istgood
          if(i1stgood.eq.0) i1stgood=i1stdef16
          i1stdef16=i1stgood
      endif

      factor=sigmaloss
      if(factor.le.0.0) factor=facdef

      group = groups_in

      do j = 1, ngroups
          hists(j)=0
      enddo

c     deadtimes
      taud = taud_in
      phases = phase_in

      tcounts=0
      do ig=1,nhists
          sum=0 !integer sum of all data in this group
          do ic=1,histlen
              sum=sum+int(rdata(ig,ic))
              if(tcounts.lt.2000000000)tcounts=tcounts+rdata(ig,ic)

          enddo
          write(histco(ig),1978) sum
1978      format('end: ',i10)
      enddo

      call multimax
      fperchan=1./(res*float(npts)*2.)
      do  i=1,n
          fchan_out(i)=float(i-1)*fperchan/135.5e-4
      enddo

      f_out = f
      if (fixphase) then
        phi_out = phase_in
      else
        phi_out = phi
      endif
      taud_out = taud

      call print_log_msg("debug", "Leaving Fortran.")

      end subroutine

      subroutine print_log_msg(log_priority, log_message)
C     A wrapper for logger function pointer, which Fortran will only ever allow
C     us to call if we store it's return value.  Using this subroutine we can
C     now just do 'call print_log_msg("debug", "message")'
      logical void
      void = logger(log_priority, log_message)
      end