      program sbc
c
c-----------------------------------------------------------------------
c
c     sbc prepares surface boundary conditions using ...
c     Hellerman wind stress monthly climatology and 
c     Oort air temperature monthly climatology
c     for MOM
c
c     note: this program is meant for internal GFDL ocean group use 
c           only and should not be used as a black box. Each user
c           must verify that the results are correct for his case.
c
c-----------------------------------------------------------------------
c
c
c     unit 39 is for calculated monthly curl tau
c     unit 41 is for taux
c     unit 42 is for tauy
c     unit 43 is for oort's air temp
c
c     slat   = starting latitude of window for viewing plots
c     slon   = starting longitude of window for viewing plots
c     elat   = ending latitude of window for viewing plots
c     elon   = ending longitude of window for viewing plots
c
      external blkdta
      character*60 name
#include "param.h"
#include "coord.h"
#include "grdvar.h"
#include "scalar.h"
      data  name/' MOM boundary conditions '/
      data slon, elon, slat, elat/ 0.0, 360.0, -90.0, 90.0/
      namelist /input/ name, slon, elon, slat, elat
#ifdef cyclic
      write (stdout,*) '  ==> cyclic option set'
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
c     set up the MOM grid coordinates
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
c     interpolate sol's monthly wind stress to MOM "u,v" grid
c-----------------------------------------------------------------------
c
      call stress (slon, elon, slat, elat)
c
c-----------------------------------------------------------------------
c     interpolate oorts monthly air temp to MOM "t" grid
c-----------------------------------------------------------------------
c
      call airtem (slon, elon, slat, elat)
c
#ifdef plots
      call last
#endif
      stop
      end 
      subroutine stress (slon, elon, slat, elat)
c
c-----------------------------------------------------------------------
c     linearly interpolate Sol Hellerman's monthly wind stress
c      climatology to the MOM grid
c-----------------------------------------------------------------------
c
      character *70 title
#include "param.h"
      parameter (ib=182, jb=90)
      parameter (imjmh=imt*jmt, lentau=3+ib+jb+12*2*ib*jb)
      common/tds/ h(imt,jmt), d(imt,jmt), a(ib,jb), b(ib,jb)
#include "grdvar.h"
#include "coord.h"
#include "levind.h"
#ifdef plots
#include "pltvar.h"
#endif
#include "scalar.h"
      common /tau/ months, lons, lats, cx(ib), cy(jb)
     $,            clmtau(ib,jb,12,2)
      dimension at(ib,jb), cxt(ib)
      dimension tx(imt,jmt), ty(imt,jmt), curl(imt,jmt)
      dimension dxr(imt), dyr(jmt), cs(jmt), csm(jmt), fcr(jmt)
#ifdef cyclic
      write (stdout,*) ' '
      write (stdout,*) ' ==> cyclic option turned on in stress'
      write (stdout,*) ' '
#endif
c
c-----------------------------------------------------------------------
c     read in Sol Hellerman's monthly wind stress climatology
c-----------------------------------------------------------------------
c
      open (unit=15,access='sequential', form='unformatted')
      read (15) months, lons, lats, cx, cy, clmtau 
      close (unit=15, status='keep')
c
c-----------------------------------------------------------------------
c     open unit 39 for MOM monthly wind stress curl
c     open unit 41 for MOM monthly taux component of wind stress
c     open unit 42 for MOM monthly tauy component of wind stress
c-----------------------------------------------------------------------
c
      open (unit=39,access='direct', form='unformatted', recl=imjmh*8)
      open (unit=41,access='direct', form='unformatted', recl=imjmh*8)
      open (unit=42,access='direct', form='unformatted', recl=imjmh*8)
c
      do 100 n=1,months
        do 90 m=1,2
c
c-----------------------------------------------------------------------
c         translate data grid so it entirely encloses the MOM "u,v"
c         grid. interpolate stress data to model "u,v" grid.
c         set last row (=jmt) to zero stress
c-----------------------------------------------------------------------
c
#ifdef cyclic
          call tranlon (clmtau(1,1,n,m), ib, ib-2, jb, at
     $,                 cx, xu(2), imt-2, cxt)
          call intrp (at, ib, ib, jb, d(2,1), imt, imt-2, jmt-1
     $,               cxt, cy, xu(2), yu)
c
          do 70 j=1,jmt
            d(1,j)   = d(imtm1,j)
            d(imt,j) = d(2,j)
70        continue
#else
          call tranlon (clmtau(1,1,n,m), ib, ib-2, jb, at
     $,                cx, xu, imt, cxt)
          call intrp (at, ib, ib, jb, d, imt, imt, jmt-1
     $,               cxt, cy, xu, yu)
#endif
          do 80 i=1,imt
            d(i,jmt) = c0
80        continue
c
          if (n .eq. 1 .and. m .eq. 1) then
            write(stdout,'(///,1x, a38,/)') 
     $       'Hellerman wind stress grid definition:'
            write (stdout,'(//,1x,a30,/)') '            latitudes:'
            write (stdout,'(1x,10f10.5)') cy
            write (stdout,'(//,1x,a30,/)') '  original longitudes:'
            write (stdout,'(1x,10f10.5)') cx
            write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
            write (stdout,'(1x,10f10.5)') cxt
          endif
c
          if (m .eq. 1) then
            nv = 41
            write (nv, rec=n) d
            write (6,*) ' ==> mean taux for month ',n,' = '
     $,     avgit(d,imt,imt,jmt,jmt,1),',  taux written to unit ', nv
          else
            nv = 42
            write (nv, rec=n) d
            write (6,*) ' ==> mean tauy for month ',n,' = '
     $,     avgit(d,imt,imt,jmt,jmt,1),',  tauy written to unit ', nv
          endif
90      continue
100   continue
c
c-----------------------------------------------------------------------
c     calculate wind stress curl
c-----------------------------------------------------------------------
c
      fact = 180.0/(c2*pi*radius)
      do 200 j=2,jmt
        cs(j)  = cos(yu(j)/radian)
        csm(j) = cos(yu(j-1)/radian)
        fcr(j) = c1/cos(yt(j)/radian)
        dyr(j) = fact/(yu(j)-yu(j-1))
200   continue
      do 300 i=2,imt
        dxr(i) = fact/(xu(i)-xu(i-1))
300   continue
c
      do 400 n=1,months
        read (41, rec=n) tx
        read (42, rec=n) ty
        do 350 j=1,jmt
          do 340 i=1,imt
            curl(i,j) = c0
340       continue
350     continue
        do 370 j=2,jmt
          do 360 i=2,imt
            curl(i,j)=fcr(j)*dxr(i)*(
     $              ty(i,j)+ty(i,j-1)-ty(i-1,j)-ty(i-1,j-1))
360       continue
370     continue
        do 390 j=2,jmt
          do 380 i=2,imt
            curl(i,j)=curl(i,j)-fcr(j)*dyr(j)*(cs(j)*
     $              (tx(i,j)+tx(i-1,j))-csm(j)*(tx(i,j-1)+tx(i-1,j-1)))
380       continue
390     continue
c
        nv = 39
        write (nv, rec=n) curl
        write (6,*) ' ==> mean curl tau for month ',n,' = '
     $, avgit(curl,imt,imt,jmt,jmt,1),', curl tau written to unit ', nv
400   continue
c
c
c-----------------------------------------------------------------------
c     annual taux & tauy & curl tau on record 13
c-----------------------------------------------------------------------
c
      do 480 nn=1,3
        if (nn .eq. 1) nv = 41
        if (nn .eq. 2) nv = 42
        if (nn .eq. 3) nv = 39
        do 420 j=1,jmt
          do 410 i=1,imt
            h(i,j) = c0
410       continue
420     continue
        do 450 n=1,12
          read(unit=nv, rec=n) d
          do 440 j=1,jmt
            do 430 i=1,imt
              h(i,j) = h(i,j) + d(i,j)
430         continue
440       continue
450     continue
        r12 = c1/12.0
        do 470 j=1,jmt
          do 460 i=1,imt
            h(i,j) = r12*h(i,j)
460       continue
470     continue
        write (unit=nv, rec=13) h
480   continue
#ifdef plots
c
c     bring in MOM kmt field
c
      nukmt = 18
      open  (unit=nukmt ,access='sequential', form='unformatted')
      rewind nukmt
      read (nukmt) kmt
c
c-----------------------------------------------------------------------
c     allow view of sub region
c-----------------------------------------------------------------------
c
      is = indp (slon, xt, imt)
      ie = indp (elon, xt, imt)
      js = indp (slat, yt, jmt)
      je = indp (elat, yt, jmt)
c
      jl = je - js + 1
      il = ie - is + 1
      sq = (yt(je)-yt(js))/(xt(ie)-xt(is))
      nc = 15
      shadng=.true.
      extshd=.false.
      shadev=0.0
      numpaf=1
      proj='rectangular'
      mppfx = 0
      if (imt .gt. 100) mppfx = imt/2
      fnflag=1.e-50
      xint = 30.0
c
      do 600 n=1,13
c
        read (unit=41, rec=n) d
        do 520 j=js,je
          do 510 i=is,ie
            if (kmt(i,j) .eq. 0) d(i,j) = fnflag
510       continue
520     continue
        if (n .le. 12) then
          write (title,'(a45,i2)') 
     $ 'Hellerman taux (dynes/cm**2) for MOM:  month=',n
        else
          write (title,'(a45)') 
     $ 'Hellerman taux (dynes/cm**2) for MOM:  Annual'
        endif
        call splt (d(is,js), imt, il, jl, xt(is), yt(js), sq, nc
     $, 'lon', 'lat', title)
c
        read (unit=42, rec=n) d
        do 540 j=js,je
          do 530 i=is,ie
            if (kmt(i,j) .eq. 0) d(i,j) = fnflag
530       continue
540     continue
        if (n .le. 12) then
          write (title,'(a45,i2)') 
     $ 'Hellerman tauy (dynes/cm**2) for MOM:  month=',n
        else
          write (title,'(a45)') 
     $ 'Hellerman tauy (dynes/cm**2) for MOM:  Annual'
        endif
        call splt (d(is,js), imt, il, jl, xt(is), yt(js), sq, nc
     $, 'lon', 'lat', title)
600   continue
#endif
      close (unit=39, status='keep')
      close (unit=41, status='keep')
      close (unit=42, status='keep')
c
      return
      end
      subroutine airtem (slon, elon, slat, elat)
c
c-----------------------------------------------------------------------
c     interp Oort's climatological air temperature to MOM "t" grid.
c-----------------------------------------------------------------------
c
      character *70 title
#include "param.h"
#include "coord.h"
#include "levind.h"
#ifdef plots
#include "pltvar.h"
#endif
      parameter (iw=73, jw=73)
      parameter (imjmh=imt*jmt)
      common /airs/ months, lons, lats, cx(iw), cy(jw)
     $,             clmair(iw,jw,12)
      dimension  d(imt,jmt), airt(iw,jw), adata(jw,iw)
      dimension airtp(iw,jw), cxt(iw), h(imt,jmt)
#ifdef cyclic
      write (stdout,*) ' '
      write (stdout,*) ' ==> cyclic option turned on in airtem'
      write (stdout,*) ' '
#endif
c
c-----------------------------------------------------------------------
c     read in Oort's monthly air temperature climatology
c-----------------------------------------------------------------------
c
      open (unit=50,access='sequential', form='unformatted')
      read (50) months, lons, lats, cx, cy, clmair 
      close (unit=50, status='keep')
c
c-----------------------------------------------------------------------
c     open unit 43 for MOM monthly air temperature
c-----------------------------------------------------------------------
c
      nw=43
      open (unit=nw,access='direct', form='unformatted', recl=imjmh*8)
c
      do 50 n=1,months
c
c-----------------------------------------------------------------------
c       translate data then interp to MOM "t" grid
c-----------------------------------------------------------------------
c
#ifdef cyclic
        call tranlon (clmair(1,1,n), iw, iw-1, jw, airtp
     $,               cx, xt(2), imt-2, cxt)
        call intrp (airtp, iw ,iw, jw ,d(2,1), imt, imt-2, jmt
     $,           cxt, cy, xt(2), yt)
c
        do 40 j=1,jmt
          d(1,j)   = d(imtm1,j)
          d(imt,j) = d(2,j)
40      continue
#else
        call tranlon (clmair(1,1,n), iw, iw-1, jw, airtp
     $,               cx, xt, imt, cxt)
        call intrp (airtp, iw ,iw, jw ,d, imt, imt, jmt
     $,           cxt, cy, xt, yt)
#endif
c
        if (n .eq. 1) then
          write(stdout,'(///,1x, a30,/)') 
     $     'Oort air temp grid definition:'
          write (stdout,'(//,1x,a30,/)') '            latitudes:'
          write (stdout,'(1x,10f10.5)') cy
          write (stdout,'(//,1x,a30,/)') '  original longitudes:'
          write (stdout,'(1x,10f10.5)') cx
          write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
          write (stdout,'(1x,10f10.5)') cxt
        endif
        write (nw, rec=n) d
        write (6,*) ' ==> mean air temp for month ',n,' = '
     $, avgit(d,imt,imt,jmt,jmt,1),',  air temp written to unit ', nw
50    continue
c
c-----------------------------------------------------------------------
c     annual air temp on record 13
c-----------------------------------------------------------------------
c
      do 420 j=1,jmt
        do 410 i=1,imt
          h(i,j) = 0.0
410     continue
420   continue
      do 450 n=1,12
        read(unit=nw, rec=n) d
        do 440 j=1,jmt
          do 430 i=1,imt
            h(i,j) = h(i,j) + d(i,j)
430       continue
440     continue
450   continue
      r12 = c1/12.0
      do 470 j=1,jmt
        do 460 i=1,imt
          h(i,j) = r12*h(i,j)
460     continue
470   continue
      write (unit=nw, rec=13) h
#ifdef plots
c
c     bring in MOM kmt field
c
      nukmt = 18
      open  (unit=nukmt ,access='sequential', form='unformatted')
      rewind nukmt
      read (nukmt) kmt
c
c-----------------------------------------------------------------------
c     allow view of sub region
c-----------------------------------------------------------------------
c
      is = indp (slon, xt, imt)
      ie = indp (elon, xt, imt)
      js = indp (slat, yt, jmt)
      je = indp (elat, yt, jmt)
c
      jl = je - js + 1
      il = ie - is + 1
      sq = (yt(je)-yt(js))/(xt(ie)-xt(is))
      nc=15
      shadng=.true.
      extshd=.false.
      shadev=0.0
      numpaf=1
      proj='rectangular'
      mppfx = 0
      if (imt .gt. 100) mppfx = imt/2
      fnflag=1.e-50
      xint = 30.0
c
      do 600 n=1,13
c
        read (unit=nw, rec=n) d
        do 520 j=js,je
          do 510 i=is,ie
            if (kmt(i,j) .eq. 0) d(i,j) = fnflag
510       continue
520     continue
        if (n .le. 12) then
          write (title,'(a45,i2)')
     $ 'Oort Air Temperature (deg C) for MOM:  month=',n
        else
          write (title,'(a45)')
     $ 'Oort Air Temperature (deg C) for MOM:  Annual'
        endif
        call splt (d(is,js), imt, il, jl, xt(is), yt(js), sq, nc
     $, 'lon', 'lat', title)
c
600   continue
#endif
      close (unit=43, status='keep')
      return
      end

