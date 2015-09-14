      subroutine tranlon (c, ic, il, jl, t, cx, fx, ifl, tx)
c
c-----------------------------------------------------------------------
c     translate longitudes  "cx" to "tx" so that tx(i) i=1..ic
c     completely encloses longitudes  fx(i) i=1..ifl
c     note that "tx" may extend beyond 360 degrees to contain "fx".
c     use same mapping to translate data in "c" to "t"
c-----------------------------------------------------------------------
c
c     c  = original data array
c     t  = translated (in longitudinal) data array
c     cx = original data longitudes
c     tx = translated data longitudes
c     fx = model longitudes
c
      dimension c(ic,jl), t(ic,jl), tx(ic), cx(ic), fx(ifl)
c
c-----------------------------------------------------------------------
c     find starting index for translation "iw"
c-----------------------------------------------------------------------
c
      if (cx(1) .le. fx(1) .and. cx(ic) .ge. fx(ifl)) then
        do 10 i=1,ic
          tx(i) = cx(i)
10      continue
        do 30 j=1,jl
          do 20 i=1,ic
            t(i,j) = c(i,j)
20        continue
30      continue
        return
      else if (cx(1) .le. fx(1)) then
        iw = 1
        do 40 i=2,ic
          if (cx(i) .le. fx(1)) then
            iw = i
          else
            go to 50
          endif
40      continue
      else
        write (6,996) ic, ifl
        write (6,998) 'cx= ',(cx(i),i=1,ic)
        write (6,998) 'fx= ',(fx(i),i=1,ifl)
        stop
      endif
50    if (iw .eq. 1 .and. cx(ic) .ge. fx(ifl)) then
        return
      endif
c
c-----------------------------------------------------------------------
c     translate longitudes such that tx(1) = cx(iw)
c-----------------------------------------------------------------------
c
      do 60 i=1,ic
        tx(i) = cx(mod(i+iw-2,il) + 1)
        im1   = max(1,i-1)
        if (tx(i) .lt. tx(im1)) tx(i) = tx(i) + 360.0
60    continue
      if (fx(ifl) .gt. tx(ic)) then
        write (6,997) iw, ic, ifl
        write (6,998) 'tx= ',(tx(i),i=1,ic)
        write (6,998) 'fx= ',(fx(i),i=1,ifl)
        stop
      endif
c
c-----------------------------------------------------------------------
c     translate data in a similar way
c-----------------------------------------------------------------------
c
      do 80 j=1,jl
        do 70 i=1,ic
          t(i,j) = c(mod(i+iw-2,il) + 1,j)
70      continue
80    continue
c
      return
996   format (1x, ' ===>  cx(1) > fx(1) in tranlon. ic=',i6,' ifl=',i6)
997   format (1x, ' ===>  tx(ic) < fx(ifl) in tranlon. iw=',i6,
     1         ' il=',i6,' ifl=',i6)
998   format (1x,a4,(5x,10e11.4))
      end
      subroutine intrp (cn,imc,ic,jc,fn,imf,if,jf,cx,cy,fx,fy)
      dimension cn(imc,jc),cx(ic),cy(jc),fn(imf,jf),fx(if),fy(jf)
c
c-----------------------------------------------------------------------
c     intrp interpolates data from grid cn to grid fn where the
c     gridspace cn completely overlies the gridspace fn.
c     imc= 1st dimension of cn in calling routine
c     ic= ending i index on cn
c     jc= ending j index on cn
c     imf= 1st dimension of fn in calling routine
c     if= ending i index on fn
c     jf= ending j index on fn
c     cx= distribution of grid points in i direction on cn
c     cy= distribution of grid points in j direction on cn
c     fx= distribution of grid points in i direction on fn
c     fy= distribution of grid points in j direction on fn
c     the distributions cx and fx (cy and fy) are arbitrary except that
c     they be monotonic increasing and consistent with each other.
c     i. e. if cx(i)=i say, then fx(i) will take on values consistent
c     with this distribution and determined by the placement of the
c     points of fn with respect to that of the points of cn.
c     note that cx(1) should be less than or equal to fx(1)
c               cy(1) should be less than or equal to fy(1)
c               cx(ic) should be greater than or equal to fx(ic)
c               cy(ic) should be greater than or equal to fy(ic)
c-----------------------------------------------------------------------
c
      if (.not.(cx(1).le.fx(1).and.cx(ic).ge.fx(if).and.
     .          cy(1).le.fy(1).and.cy(jc).ge.fy(jf))) go to 21
      do 12 i=2,ic
        if (cx(i).lt.cx(i-1)) go to 20
 12   continue
      do 14 j=2,jc
        if (cy(j).lt.cy(j-1)) go to 20
 14   continue
      do 16 i=2,if
        if (fx(i).lt.fx(i-1)) go to 20
 16   continue
      do 18 j=2,jf
        if (fy(j).lt.fy(j-1)) go to 20
 18   continue
c
      do 100 j=1,jf
        do 110 l=1,jc
          if (cy(l).ge.fy(j)) then
            jn=l-1
            js=l
            if (jn.eq.0) js=2
            if (jn.eq.0) jn=1
            go to 111
          endif
 110    continue
 111    continue
        if (jf*jc.eq.1) then
          c=0.5
          d=0.5
          ty=1.0
          js=1
          jn=1
        else
          c=fy(j)-cy(jn)
          d=cy(js)-fy(j)
          ty=cy(js)-cy(jn)
        endif
        do 101 i=1,if
          do 120 l=1,ic
            if (cx(l).ge.fx(i)) then
              ie=l
              iw=l-1
              if (iw.eq.0) ie=2
              if (iw.eq.0) iw=1
              go to 121
            endif
 120      continue
 121      continue
          if (if*ic.eq.1) then
            a=0.5
            b=0.5
            tx=1.0
            ie=1
            iw=1
          else
            a=fx(i)-cx(iw)
            b=cx(ie)-fx(i)
            tx=cx(ie)-cx(iw)
          endif
          fn(i,j)=(c*(a*cn(ie,js)+b*cn(iw,js))+
     .             d*(a*cn(ie,jn)+b*cn(iw,jn)))/(tx*ty)
  101   continue
  100 continue
c
      return
c
   20 continue
      write (6,9998)
      go to 30
   21 continue
      write (6,9999)
   30 continue
c
      write (6,9994) 'cx='
      write (6,9993) (cx(i),i=1,ic)
      write (6,9994) 'cy='
      write (6,9993) (cy(j),j=1,jc)
      write (6,9994) 'fx='
      write (6,9993) (fx(i),i=1,if)
      write (6,9994) 'fy='
      write (6,9993) (fy(j),j=1,jf)
c
      stop
 9999 format (1x,'grid interpolated to exceeds grid interploated from')
 9998 format (1x,'cx, cy, fx, or fy not monotonic')
 9994 format (1x,a4)
 9993 format (1x,10e13.6)
c
      end
      function avgit (a, id, il, jd, jl, kl)
      dimension a(id,jd,kl)
c
c     a function to return a simple mean of the array a(id,jd,kl)
c     over the region defined by 1:il, 1:jl, 1:kl.
c
      sum = 0.0
      do 100 k=1,kl
        do 90 j=1,jl
          do 80 i=1,il
            sum = sum + a(i,j,k)
80        continue
90      continue
100   continue
c
c       calculate and return the mean
c
      avgit = sum / real(il * jl * kl)
c
      return
      end
