      program ic
c
c-----------------------------------------------------------------------
c     prepare initial temperature & salt for MOM based on
c     monthly Levitus climatology
c
c     note: this program is meant for internal GFDL ocean group use 
c           only and should not be used as a black box. Each user
c           must verify that the results are correct for his case.
c
c-----------------------------------------------------------------------
c
#include "param.h"
      parameter (imd=362, jmd=180, kd=33)
      parameter (lend=imd*jmd+imd+jmd+kd+4, long=8*lend)
      parameter (ni=15, nj=15, nk=15)
      parameter (flag=1.e24)
      parameter (imkmh=imt*km)
c
      character*4 tors
      character*3 months(13)
      character*3 curmon
      character*36 title
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
      common /bigone/ month, lons, lats, depths, cx(imd), cy(jmd)
     $,               dd(kd), data(imd,jmd)
c
#ifdef plots
#include "pltvar.h"
#endif
      dimension record(lend)
      equivalence (record, month)
      common /large/ block(imt,jmt,kd)
c
      common /bufrs/ s(km), djk(jmt,km), dkj(km,jmt), tm(imt,km),
     1               tij(imt,jmt), tdv(imt,kd)
c
      dimension jv(nj), iv(ni), kv(nk), td(imd,jmd)
      dimension tdp(imd,jmd), cxt(imd)
      dimension logi(kd), rlat(nj), rlon(ni), rdepth(nk)
c
      data months /'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul'
     $,            'aug', 'sep', 'oct', 'nov', 'dec', 'ann'/
c
c     tors = 'temp' or 'salt' to choose temperature or salinity data
c
      data tors/'temp'/
c
c     rlat   = latitude(s) at which plots are desired
c     rlon   = longitude(s) at which plots are desired
c     rdepth = depth(s) at which plots are desired
c
      data rlat /nj*flag/, rlon /ni*flag/, rdepth /nk*flag/
      data  nuik, nujk, nuij, nukmt, nuann, nuseas, nusurf
     1     /23, 24, 25, 18, 20, 19, 44/
c
c     slat   = starting latitude of window for viewing plots
c     slon   = starting longitude of window for viewing plots
c     elat   = ending latitude of window for viewing plots
c     elon   = ending longitude of window for viewing plots
c
      data slon, elon, slat, elat/ 0.0, 360.0, -90.0, 90.0/
c
      namelist /input/  tors, nukmt, nuik, nujk, nuij, nusurf
     $,                 nuann, nuseas, rlat, rlon, rdepth,  curmon
     $,                 slon, elon, slat, elat
#ifdef cyclic
      write (stdout,*) '  ==> cyclic option set'
#endif
#ifdef plots
      write (stdout,*) '  ==> plots option set'
      call start
#endif
      rewind stdin
      read  (stdin, input)
      write (stdout,input)
c
      open (unit=nuann,  access='direct', form='unformatted', recl=long)
      open (unit=nuseas, access='direct', form='unformatted', recl=long)
c
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
c     read  MOM "kmt" field 
c-----------------------------------------------------------------------
c
      open (unit=nukmt, access='sequential', form='unformatted')
      read (nukmt) kmt
      write (stdout,*) ' '
      write (stdout,*) ' MOM levels read from unit ', nukmt
      close (unit=nukmt, status='keep')
      open (unit=nuik
     $,           access='direct', form='unformatted', recl=imt*km*8)
c
c-----------------------------------------------------------------------
c     bring in the grid definition for levitus climatology
c-----------------------------------------------------------------------
c
      read (nuann, rec=1) record
#ifdef cyclic
      call tranlon (data, imd, imd-2, jmd, tdp, cx, xt(2), imt-2, cxt)
#else
      call tranlon (data, imd, imd-2, jmd, tdp, cx, xt, imt, cxt)
#endif
c
      write(stdout,'(///,1x, a38,/)')
     $   'Levitus grid definition:'
      write (stdout,'(//,1x,a30,/)') '            latitudes:'
      write (stdout,'(1x,10f10.5)') cy
      write (stdout,'(//,1x,a30,/)') '  original longitudes:'
      write (stdout,'(1x,10f10.5)') cx
      write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
      write (stdout,'(1x,10f10.5)') cxt
      write (stdout,'(//,1x,a30,/)') '               depths:'
      write (stdout,'(1x,10g10.2)') dd
c
c-----------------------------------------------------------------------
c     identify which data levels are above & below MOM levels to
c     eliminate unnecessary horizontal interpolations
c-----------------------------------------------------------------------
c
      do 40 k=1,kd
        logi(k) = 0
40    continue
      do 60 k=1,km
        do 50 l=1,kd
          if (dd(l) .ge. zt(k)) then
            logi(l)          = 1
            logi(max(1,l-1)) = 1
            if (dd(l) .eq. zt(k)) logi(min(kd,l+1))=1
            go to 60
          endif
50      continue
60    continue
c
      write (6,9000)
      do 70 k=1,kd
        if (logi(k) .eq. 1) write (6,9001) k, dd(k)
70    continue
c
c-----------------------------------------------------------------------
c     horizontally interpolate levitus data to MOM "t" grid
c-----------------------------------------------------------------------
c
c     temperature varies monthly down to level 19
c     salinity varies monthly down to level 24
c
      lev = 19
      if (tors .ne. 'temp') lev = 24
      print *,' => lev = ',lev
c
      do 160 k=1,kd
c
        if (logi(k) .eq. 0) then
          do 120 j=1,jmt
            do 110 i=1,imt
              block(i,j,k) = -10000.
110         continue
120       continue
          go to 160
        endif
c
c       read in Levitus data at level k
c
        if (k .le. lev) then
          read (nuseas, rec=k) record
        else
          read (nuann, rec=k) record
        endif
c
c       translate levitus grid & data on level k
c       horizontally interpolate to the MOM "t" grid 
c
#ifdef cyclic
        call tranlon (data, imd, imd-2, jmd, tdp, cx, xt(2), imt-2, cxt)
        call intrp (tdp, imd, imd, jmd, block(2,1,k), imt, imt-2, jmt,
     1              cxt, cy, xt(2), yt)
c
          do 130 j=1,jmt
            block(1,j,k)   = block(imtm1,j,k)
            block(imt,j,k) = block(2,j,k)
130       continue
#else
        call tranlon (data, imd, imd-2, jmd, tdp, cx, xt, imt, cxt)
        call intrp (tdp, imd, imd, jmd, block(1,1,k), imt, imt, jmt,
     1              cxt, cy, xt, yt)
#endif
c
        if (tors .ne. 'temp') then
          do 150 j=1,jmt
            do 140 i=1,imt
              block(i,j,k) = 1.e-3*block(i,j,k)-.035
140         continue
150       continue
        endif
160   continue
c
c-----------------------------------------------------------------------
c    interpolate vertically to MOM levels for i-k slabs
c-----------------------------------------------------------------------
c
      do 300 j=1,jmt
        do 210 k=1,kd
          do 200 i=1,imt
            tdv(i,k) = block(i,j,k)
200       continue
210     continue
c
        call intrp (tdv, imt, imt, kd, tm, imt, imt, km, xt, dd, xt, zt)
c
c       set data to zero in land & save surface data
c
        do 230 i=1,imt
          kz = kmt(i,j)
          do 220 k=1,km
            if (k .gt. kz) then
              tm(i,k) = c0
            endif
220       continue
          tij(i,j) = tm(i,1)
230     continue
c
        write (nuik, rec=j) tm
c
	if (j .eq. jmt/2) then
	  write (6,'(//"j=",i2/10(1pe10.3,1x))')
     .          j,(tm(i,1),i=1,imt)
	  write (6,'(//)')
	endif
300   continue
c
c-----------------------------------------------------------------------
c     store surface data on MOM "t" grid for each month (annual is 13)
c-----------------------------------------------------------------------
c
      open (unit=nusurf, access='direct', form='unformatted'
     $,     recl=imt*jmt*8)
      n = 0
      do 310 nn=1,13
        if (curmon .eq. months(nn)) then
          n = nn
          go to 320
        endif
310   continue
320   continue
      if (n .eq. 0) then
        print *,' stop in ic.F because curmon=',curmon
        stop 
      endif
      write (unit=nusurf, rec=n) tij
      print *,' ==> writing ',tors,' surface data for month = '
     $ ,n,' on unit =',nusurf
#ifdef plots
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
      nc=20
      shadng=.true.
      extshd=.false.
      shadev=0.0
      numpaf=1
      proj='rectangular'
      mppfx = 0
      if (imt .gt. 100) mppfx = imt/2
      fnflag=1.e-50
      xint = 30.0
      yint = 15.0
      do 420 j=js,je
        do 410 i=is,ie
          if (kmt(i,j) .eq. 0) tij(i,j) = fnflag
410     continue
420   continue
      if (tors .eq. 'temp') then
        write (title,'(a28,i2,a2,a3,a1)')
     $     'Levitus SST for MOM: month =',n, ' (',curmon,')'
      else
        write (title,'(a28,i2,a2,a3,a1)')
     $     'Levitus SSS for MOM: month =',n, ' (',curmon,')'
      endif
      call splt (tij(is,js), imt, il, jl, xt(is), yt(js), sq, nc
     $, 'lon', 'lat', title)
c     
c
c     open units for j-k & i-j storage
c
      open (unit=nujk, status='scratch'
     $,       access='direct', form='unformatted', recl=jmt*km*8)
      open (unit=nuij, status='scratch'
     $,       access='direct', form='unformatted', recl=imt*jmt*8)
c
c-----------------------------------------------------------------------
c     get the latitudes, longitudes & depths for sections
c-----------------------------------------------------------------------
c
      do 400 j=1,nj
        jv(j) = -1
        if (rlat(j) .ne. flag) then
          jv(j)   = indp(rlat(j), yt, jmt)
          rlat(j) = yt(jv(j))
        endif
400   continue
      do 500 i=1,ni
        iv(i) = -1
        if (rlon(i) .ne. flag) then
          iv(i) = indp(rlon(i), xt, imt)
          rlon(i) = xt(iv(i))
        endif
500   continue
      do 600 k=1,nk
        kv(k) = -1
        if (rdepth(k) .ne. flag) then
          kv(k) = indp(rdepth(k), zt, km)
          rdepth(k) = zt(kv(k))
        endif
600   continue
c
c-----------------------------------------------------------------------
c     fill the 3d memory block with data
c-----------------------------------------------------------------------
c
      do 700 j=1,jmt
        read (nuik, rec=j) tm
        do 690 i=1,imt
          do 680 k=1,km
            block(i,j,k) = tm(i,k)
680       continue
690     continue
700   continue
c
c-----------------------------------------------------------------------
c     create i-j storage
c-----------------------------------------------------------------------
c
      do 800 k=1,km
        do 790 j=1,jmt
          do 780 i=1,imt
            tij(i,j) = block(i,j,k)
780       continue
790     continue
        write (nuij, rec= k) tij
800   continue
c
c-----------------------------------------------------------------------
c     create j-k storage
c-----------------------------------------------------------------------
c
      do 900 i=1,imt
        do 890 k=1,km
          do 880 j=1,jmt
            djk(j,k) = block(i,j,k)
880       continue
890     continue
        write (nujk, rec=i) djk
900   continue
c
c-----------------------------------------------------------------------
c     inspect sections
c-----------------------------------------------------------------------
c
      call sample (nuik,  tm, imt, km, jv, nj, zt, s, nuik, nujk
     $,           nuij, tors, is, ie, js, je, curmon)
      call sample (nujk, djk, jmt, km, iv, ni, zt, s, nuik, nujk
     $,           nuij, tors, is, ie, js, je, curmon)
      call sample (nuij,tij, imt, jmt, kv, nk, zt, tm, nuik, nujk
     $,           nuij, tors, is, ie, js, je, curmon)
      call last
#endif
      close (unit=nusurf, status='keep')
      close (unit=nuik, status='keep')
      stop
9000  format (1x,'data levels for horizontal interpolations',/)
9001  format (5(3x,i3,e14.7))
      end
      subroutine sample (nu, a, lx, ly, ind, ni, y, s, nuik, nujk, nuij
     $,                  tors, is, ie, js, je, curmon)
#ifdef plots
#include "param.h"
      character*36 text
      character*(*) tors
      character*3 curmon
#include "coord.h"
#include "grdvar.h"
#include "pltvar.h"
      dimension a(lx,ly), ind(ni), y(ly), s(ly)
      shadng=.true.
      extshd=.true.
      shadev=0.0
      numpaf=1
      proj='rectangular'
      mppfx=0
      fnflag=1.e-50
      yint = 15.0
      xint = 30.0
      zint = 250.0e2
      nc = 20
      do 100 n=1,ni
        idd=ind(n)
        if (idd.le.0) return
c
        read (nu, rec=idd) a
c
        do 60 jj=1,ly
          do 50 ii=1,lx
            if (a(ii,jj) .eq. c0) a(ii,jj) = fnflag
50        continue
60      continue
c
        if (nu .eq. nuik) then
          write (text,9001)
     $      curmon,' levitus ', tors, yt(idd)
        endif
        if (nu .eq. nujk) then
          write (text,9000) 
     $      curmon,' levitus ', tors, xt(idd)
        endif
        if (nu .eq. nuij) then
          write (text,9002)
     $      curmon,' levitus ', tors, y(idd)
        endif
        if (nu .ne. nuij) then
          do 1 l=1,lx
            do 2 jj=1,ly
              s(ly+1-jj) = a(l,jj)
2           continue
            do 3 j=1,ly
              a(l,j) = s(j)
3           continue
1         continue
          do 4 jj=1,ly
            s(jj) = y(ly+1-jj)
4         continue
          if (nu .eq. nuik) then
            mppfx = 0
            if (imt .gt. 100) mppfx = imt/2
            il = ie - is + 1
            call splt (a(is,1), lx, il, ly, xt(is), s,
     $                              1.0, nc, 'lon', 'depth', text)
          endif
          if (nu .eq. nujk) then
            mppfx = 0
            jl = je - js + 1
            call splt (a, lx, jl, ly, yt(js), s,
     $                              1.0, nc, 'lat', 'depth', text)
          endif
        else
          mppfx = 0
          if (imt .gt. 100) mppfx = imt/2
          il = ie - is + 1
          sq = (yt(je)-yt(js))/(xt(ie)-xt(is))
          call splt (a(is,js), lx, il, ly, xt(is), yt(js), sq, nc
     $,      'lon', 'lat', text)
        endif
        do 90 jj=1,ly
          do 80 ii=1,lx
            if (a(ii,jj) .eq. fnflag) a(ii,jj) = c0
80        continue
90      continue
100   continue
#endif
      return
 9000 format (1x,a3,a9,a4,' lon = ',f8.2)
 9001 format (1x,a3,a9,a4,' lat = ',f8.2)
 9002 format (1x,a3,a9,a4,' depth =',f8.0,'cm')
      end
