      program topo
c
c-----------------------------------------------------------------------
c
c     prepare topography & geometry for MOM by linearly interpolating
c     Scripps one degree topography data onto the MOM "t" grid to
c     produce "kmt"
c
c     note: this program is meant for internal GFDL ocean group use 
c           only and should not be used as a black box. Each user
c           must verify that the results are correct for his case.
c
c-----------------------------------------------------------------------
c
c
c     iw   = # of longitues for Scripps 1 deg topography (including
c            cyclic condition)
c     jw   = # of latitudes for Scripps 1 deg topography 
c     d    = array of Scripps 1 deg topography (meters)
c     cx   = longitudes of Scripps data
c     cy   = latitudes of Scripps data
c     dt   = "d" translated in longitude so model grid lies entirely
c             within data grid 
c     cxt  = "cx" translated in longitude so model grid lies entirely
c            within data grid 
c     a    = Scripts topography on MOM "t" grid
c     nfil = number of filter applications
c     kmt  = Scripts topography discretized to model levels on "t" grid
c     dmin = miniumum allowable depth for ocean model points (meters)
c
c     am   = lateral mixing coeff in cm**2/sec (for Killworth stability
c            analysis)
c     dtuv = time step in sec (for Killworth stability analysis)
c     alonis = longitude(s) for island seed(s). see index.h
c     alatis = latitude(s) for island seed(s). see index.h
c     nudata = unit to read Scripps 1 deg topography
c     nutopo = unit to write model depths at "t" points
c     nukmt  = unit to write model levels at "t" points
c
c     slat   = starting latitude of window for viewing plots
c     slon   = starting longitude of window for viewing plots
c     elat   = ending latitude of window for viewing plots
c     elon   = ending longitude of window for viewing plots
c
      external blkdta
      character*60 name
#include "param.h"
      parameter (iw=362, jw=180)
#include "coord.h"
#include "grdvar.h"
#include "index.h"
#include "levind.h"
#ifdef plots
#include "pltvar.h"
#endif
#include "scalar.h"
      dimension d(iw,jw), cx(iw), cy(jw)
      dimension dt(iw,jw), cxt(iw)
      dimension scr(iw,jw)
      dimension a(imt,jmt), b(imt,jmt), c(imt,jmt)
c
      data  dmin /100.0/, nutopo/83/, nukmt/84/, nrdata/22/
      data slon /0.0/, elon/360.0/, slat/-90.0/, elat/90.0/
      data nfil /0/
      data  name/' program topo '/
c
      namelist /input/ dmin, nutopo, nukmt, name, slon, elon
     $,                slat, elat, nfil, am, dtuv, alonis, alatis
c
#ifdef cyclic
      write (stdout,*) '  ==> cyclic option set'
#endif
#ifdef removebays
      write (stdout,*) '  ==> remove isolated bays option set'
#endif
#if defined islands && defined rigidlid
      write (stdout,*) '  ==> rigidlid & islands option set'
#endif
#ifdef plots
      write (stdout,*) '  ==> plots option set'
      call start
#endif
      rewind stdin
      read  (stdin,input)
      write (stdout,input)
c
c-----------------------------------------------------------------------
c     read  Scripps 1 deg topography data
c-----------------------------------------------------------------------
c
      open (unit=nrdata,access='sequential', form='unformatted')
      read (nrdata) lons, lats, cx, cy, d
      close (unit=nrdata, status='keep')
      write (stdout,*) ' '
      write (stdout,*) ' Scripps 1 deg Topography read from unit '
     $,                 nrdata 
c
c-----------------------------------------------------------------------
c     smooth Scripps topography here (retain coastlines)
c-----------------------------------------------------------------------
c
      if (nfil .ne. 0) then
        do 20 j=1,jw
          do 10 i=1,iw
            dt(i,j) = c1
            if (d(i,j) .ge. c0) then
              d(i,j)  = c0
              dt(i,j) = c0
            endif
10        continue
20      continue
c
c       choose an area to smooth
c
        is = indp (0.5,   cx, lons)
        ie = indp (359.5, cx, lons)
        js = indp (-89.5, cy, lats)
        je = indp ( 89.5, cy, lats)
        call sfltr (d, iw, iw, jw, scr, dt, is, ie, js, je, nfil)
      endif
c
c-----------------------------------------------------------------------
c     set up the model grid coordinates
c-----------------------------------------------------------------------
c
      write (stdout,*) ' '
      write (stdout,*) ' Model grid point definition:'
c
c     set constants needed by "grids"
c
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
c
      call grids
c
c-----------------------------------------------------------------------
c     translate data grid so model grid lies entirely within data grid
c     interpolate topography data to MOM "t" grid
c-----------------------------------------------------------------------
c
#ifdef cyclic
      call tranlon (d, iw, iw-2, jw, dt, cx, xt(2), imt-2, cxt)
      call intrp (dt, iw, iw, jw, a(2,1), imt, imt-2, jmt
     $,           cxt, cy, xt(2), yt)
c
      do 30 j=1,jmt
        a(1,j)   = a(imtm1,j)
        a(imt,j) = a(2,j)
30    continue
#else
      call tranlon (d, iw, iw-2, jw, dt, cx, xt, imt, cxt)
      call intrp (dt, iw, iw, jw, a, imt, imt, jmt
     $,           cxt, cy, xt, yt)
#endif
c
      write(stdout,'(///,50x, a30,/)') 'Scripps 1 deg topography grid'
      write (stdout,'(//,1x,a30,/)') '            latitudes:'
      write (stdout,'(1x,10f10.5)') cy
      write (stdout,'(//,1x,a30,/)') '  original longitudes:'
      write (stdout,'(1x,10f10.5)') cx
      write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
      write (stdout,'(1x,10f10.5)') cxt
c
c-----------------------------------------------------------------------
c     convert depths from meters to cm and change sign. set mask
c     for filtering
c-----------------------------------------------------------------------
c
      dptmin = dmin*100.0
      do 60 j=1,jmt
        do 58 i=1,imt
          if (a(i,j) .lt. c0) then
            a(i,j) = max(dptmin,-a(i,j)*100.0)
            b(i,j) = c1
          else
            a(i,j) = c0
            b(i,j) = c0
          endif
58      continue
60    continue
c
c-----------------------------------------------------------------------
c     set topography modifications on model grid here
c-----------------------------------------------------------------------
c
c     set artificial walls at northern & southern boundaries
c
      do 70 i=1,imt
        a(i,1 )  = c0
        a(i,jmt) = c0
70    continue
c
#ifdef plots
c
c-----------------------------------------------------------------------
c     allow view of sub region of model domain
c-----------------------------------------------------------------------
c
      n = nchar(name)
      doc(1:) = ' topo: '//name(1:n)
      if (nfil .ne. 0) then
        n = nchar(doc)
        write (doc(n+1:),'(a2,i2,a11)') ' (', nfil, ' filter(s)))'
      endif
      ndoc = nchar(doc)
      doc = doc(1:ndoc)//' (Scripps data in model domain)'
      is = indp (slon, xt, imt)
      ie = indp (elon, xt, imt)
      js = indp (slat, yt, jmt)
      je = indp (elat, yt, jmt)
c
      call foverh (a, b, is, ie, js, je)
#endif
c
c-----------------------------------------------------------------------
c     discretize topography to  "kmt" (min kmt = 2) & print "kmt"
c-----------------------------------------------------------------------
c
      do 110 j=1,jmt
        do 100 i=1,imt
          kmt(i,j) = 0
          if (a(i,j) .ne. c0) then
            kmt(i,j) = indp (a(i,j), zw, km)
            a(i,j)   = zw(kmt(i,j))
            if (kmt(i,j) .lt. 2) write(6,9008) i, j, kmt(i,j), a(i,j)
          endif
100     continue
110   continue
c
c-----------------------------------------------------------------------
c     compute number of vertical levels on the "u" grid
c-----------------------------------------------------------------------
c
      do 310 j=1,jmt
        kmu(imt,j) = 0
310   continue
c
      do 320 i=1,imt
        kmu(i,jmt) = 0
320   continue
c
      do 340 j=1,jmtm1
        do 330 i=1,imtm1
         kmu(i,j) = min (kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
330     continue
340   continue
#ifdef cyclic
      do 350 j=1,jmt
        kmu(imt,j) = kmu(2,j)
350   continue
#endif
#ifdef symmetry
      do 360 i=1,imt
        kmu(i,jmt) = kmu(i,jmtm2)
360   continue
#endif
c
c-----------------------------------------------------------------------
c     search for isolated bays... "t" grid boxes at the surface which
c     can't be influenced by advection
c-----------------------------------------------------------------------
c     
      do 400 j=2,jmtm1
        do 390 i=2,imtm1
          if (kmt(i,j) .ne. 0) then
            if (kmu(i,j) .eq. 0 .and. kmu(i-1,j) .eq. 0 .and.
     $        kmu(i,j-1) .eq. 0 .and. kmu(i-1,j-1) .eq. 0) then
#ifdef removebays
              write (stdout,'(10x,a42,i4,a1,i4,a9,i3,a20)')
     $        '==>   Warning: isolated "kmt" at (i,j) = (',i,',',j
     $,       '), kmt = ', kmt(i,j),' is being reset to 0'
              kmt(i,j) = 0
#else
              write (stdout,'(10x,a42,i4,a1,i4,a9,i3)')
     $        '==>   Warning: isolated "kmt" at (i,j) = (',i,',',j
     $,       '), kmt = ', kmt(i,j)
#endif
            endif
          endif
390     continue
400   continue
c
c-----------------------------------------------------------------------
c     search for topographic instabilities (based  on the  work of
c     Peter Killworth  ...  eqn 11 from ocean modeling nov 1987)
c-----------------------------------------------------------------------
c
      write (stdout,'(/1x,a59,a26)/')
     $ '==> Warning:  Make sure "am" & "dtuv" are set to the values'
     $, ' intended for use with MOM'
c
      num   = 50
      do 500 j=2,jmtm1
        dysq = dyt(j)**2
        do 490 i=2,imtm1
          if (kmt(i+1,j-1) .ne. 0 .and. kmt(i+1,j) .ne. 0) then
            dxsq = (dxt(i)*cst(j))**2
            clix = am*dtuv/dxsq
            h1   = zw(kmt(i+1,j-1))
            h2   = zw(kmt(i+1,j))
            hx   = (8.0*h1*h2/(h1+h2)**2 + dxsq/dysq)/(4.0 + dxsq/dysq)
            if (clix .ge. hx .and. num .ge. 0) then
              num = num - 1
              write (stdout,'(1x,a56,a31,i4,a1,i4,a1)')
     $        '==>   Warning: Killworth topographic roughness condition'
     $,       ' exceeded at location (i,j) = (',i+1,',',j,')'
              if (num .eq. 0) then
                write (stdout,9000)
     $         '==> Warning: msgs terminated after 50 cases were found '
              endif
            endif
          endif
490     continue
500   continue
c
c-----------------------------------------------------------------------
c     print "kmt" field 
c-----------------------------------------------------------------------
c
      write (stdout,'(50x,a28)') '"kmt" field for use with MOM'
      call pr (kmt, imt, jmt)
c
c-----------------------------------------------------------------------
c     write out h & kmt fields (note: h is not used by MOM)
c-----------------------------------------------------------------------
c
      open  (unit=nutopo, access='sequential', form='unformatted')
      write (nutopo) a
      close (unit=nutopo, status='keep')
c
      open  (unit=nukmt ,access='sequential', form='unformatted')
      write (nukmt) kmt
      close (unit=nukmt, status='keep')
c
#ifdef plots
c
c-----------------------------------------------------------------------
c     allow view of sub region after discretization to model levels
c-----------------------------------------------------------------------
c
      doc = doc(1:ndoc)
     $      //'(Scripps data discretized to model levels)'
      is = indp (slon, xt, imt)
      ie = indp (elon, xt, imt)
      js = indp (slat, yt, jmt)
      je = indp (elat, yt, jmt)
c
      call foverh (a, b, is, ie, js, je)
      call last  
#endif
#ifdef islands
# ifdef removebays
c
c-----------------------------------------------------------------------
c     compute number of vertical levels on the "u" grid
c-----------------------------------------------------------------------
c
      do 510 j=1,jmt
        kmu(imt,j) = 0
510   continue
c
      do 520 i=1,imt
        kmu(i,jmt) = 0
520   continue
c
      do 600 j=1,jmtm1
        do 590 i=1,imtm1
         kmu(i,j) = min (kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
590     continue
600   continue
#  ifdef cyclic
      do 710 j=1,jmt
        kmu(imt,j) = kmu(2,j)
710   continue
#  endif
#  ifdef symmetry
      do 720 i=1,imt
        kmu(i,jmt) = kmu(i,jmtm2)
720   continue
#  endif
# endif
c
c-----------------------------------------------------------------------
c     find island perimeter points
c-----------------------------------------------------------------------
c
      call iperim (b)
#endif
      stop
9000  format (/,(1x,a80))
 9008 format (/,10x,' !!!! too few levels at i=',i4, ', j=',i4,
     1        ', levels=',i3,', depth=',e14.7,'  !!!!',/)
      end
      subroutine pr (kmt, im, jm)
c
c-----------------------------------------------------------------------
c        print map of levels
c-----------------------------------------------------------------------
c
      dimension kmt(im,jm)
      inc=40
      do 20 l=0,im,inc
        incr = min(inc,im-l)
        write (6,9300) (l+i,i=1,incr,2)
        do 10 jj=1,jm
          j = jm + 1 - jj
          write (6,9100) j, (kmt(l+i,j),i=1,incr)
10      continue
20    continue
      return
9100  format (1x,i3,1x,40i3)
9300  format (///,2x,20i6/)
      end
      subroutine foverh (d, foh, is, ie, js, je)
#ifdef plots
c
c-----------------------------------------------------------------------
c     construct f/h field & plot it
c-----------------------------------------------------------------------
c
#include "param.h"
#include "pltvar.h"
      dimension d(imt,jmt), f(jmt), foh(imt,jmt)
#include "coord.h"
      data omega/7.2921e-5/, degrad/0.0174533/
c
      do 10 j=1,jmt
        f(j) = c2*omega*sin(yt(j)*degrad)
10    continue
c
c-----------------------------------------------------------------------
c     limit contours to open ocean values
c-----------------------------------------------------------------------
c
      zmin=1.5e5
      do 20 j=1,jmt
        do 20 i=1,imt
          foh(i,j) = f(j)/max(d(i,j),zmin)
20      continue
30    continue
c
c-----------------------------------------------------------------------
c     plot h & f/h fields
c-----------------------------------------------------------------------
c
      jl = je - js + 1
      il = ie - is + 1
      sq = (yt(je)-yt(js))/(xt(ie)-xt(is))
      nc=10
      shadng=.true.
      extshd=.true.
      shadev=0.0
      numpaf=1
      proj='rectangular'
      mppfx = 0
      if (imt .gt. 100) mppfx = imt/2
      fnflag=1.e-50
      xint = 30.0
      do 50 j=1,jmt
        do 40 i=1,imt
         if (d(i,j) .eq. 0.0) then
           d(i,j)   = fnflag
           foh(i,j) = fnflag
         endif
40      continue
50    continue
      call splt (d(is,js), imt, il, jl, xt(is), yt(js), sq, nc
     $, 'lon', 'lat', ' Topography')
      call splt (foh(is,js), imt, il, jl, xt(is), yt(js), sq, 2*nc
     $, 'lon', 'lat', ' f/h')
      do 70 j=1,jmt
        do 60 i=1,imt
         if (d(i,j) .eq. fnflag) then
           d(i,j)   = 0.0
           foh(i,j) = 0.0
         endif
60      continue
70    continue
#endif
      return
      end
      subroutine sfltr (a, im, il, jl, s, f, i1, i2, j1, j2, nfil)
c
c-----------------------------------------------------------------------
c     this 2 dimensional symmetric filter is based on a weighting of
c     1/4  1/2  1/4. 
c     f is the mask field = 0 over land & 1 in ocean
c     a is the data to be smooth. It must be zero over land !
c     s is a scratch array
c     nfil is the number of smoothing passes
c     i1 .. i2 are the points along the 1st dimension to smooth
c     j1 .. j2 are the points along the 2nd dimension to smooth
c     this filter behaves as if there is a no flux boundary condition
c     across coastlines. It therefore preserves coasts.
c     note: this filter does not take into account variable grids
c-----------------------------------------------------------------------
c           
      dimension a(im,jl), s(im,jl), f(im,jl)
c
      if (nfil .eq. 0) return
c
      do 20 j=1,jl
        do 10 i=1,il
          if (f(i,j) .eq. 0.0 .and. a(i,j) .ne. 0.0) then
            print *,' "a" is not zero but "f" is zero in "sfltr"'
            print *,' i=',i,' j=',j,' a(i,j)=',a(i,j)
            stop
          endif
10      continue
20    continue
c
      is = max (i1,2)
      ie = min (i2,il-1)
      js = max (j1,2)
      je = min (j2,jl-1)
      print *,' '
      print *,' ==> Applying smoothing filter ', nfil,' times over'
      print *,
     $'     domain defined by is=',is,', ie=',ie,', js=',js,', je=',je
      print *,' '
c
      do 100 n=1,nfil
        do 50 j=js,je
          do 40 i=2,il-1
            if (i .ge. is .and. i .le. ie) then
              x1p = 0.25*f(i,j)
              x2p = 1.0-0.25*(f(i-1,j+1)+f(i+1,j+1))
              x3p = 0.25*f(i,j)
c            
              x1c = 0.25*f(i,j)
              x2c = 1.0-0.25*(f(i-1,j)+f(i+1,j))
              x3c= 0.25*f(i,j)
c
              x1m = 0.25*f(i,j)
              x2m = 1.0-0.25*(f(i-1,j-1)+f(i+1,j-1))
              x3m = 0.25*f(i,j)
c
              y1m = 0.25*f(i,j)
              y2m = 1.0-0.25*(f(i-1,j+1)+f(i-1,j-1))
              y3m = 0.25*f(i,j)
c
              y1c = 0.25*f(i,j)
              y2c = 1.0-0.25*(f(i,j+1)+f(i,j-1))
              y3c = 0.25*f(i,j)
c
              y1p = 0.25*f(i,j)
              y2p = 1.0-0.25*(f(i+1,j+1)+f(i+1,j-1))
              y3p = 0.25*f(i,j)
c
              s(i,j) = 
     $        x1p*y1m*a(i-1,j+1) + x2p*y1c*a(i,j+1) + x3p*y1p*a(i+1,j+1)
     $      + x1c*y2m*a(i-1,j)   + x2c*y2c*a(i,j)   + x3c*y2p*a(i+1,j)
     $      + x1m*y3m*a(i-1,j-1) + x2m*y3c*a(i,j-1) + x3m*y3p*a(i+1,j-1)
            else
              s(i,j) = a(i,j)
            endif
40        continue  
50      continue
#ifdef cyclic
        do 60 j=js,je
          s(1,j)  = s(il-1,j)
          s(il,j) = s(2,j)
60      continue
#else
        do 60 j=js,je
          s(1,j)  = 0.0
          s(il,j) = 0.0
60      continue
#endif
        do 80 j=js,je
          do 70 i=1,il
            a(i,j) = s(i,j)
70        continue  
80      continue
100   continue
      return
      end
