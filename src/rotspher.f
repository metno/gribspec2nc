      SUBROUTINE ancorr(delta,rlon,rlat,xcen,ycen)
c
c     Direction angle correction, from rotated Wam grid to geographic.
c
c     CODED BY: J.E.HAUGEN, A.FOSS   - DNMI/R&D  APRIL 1998
c
      real delta,rlon,rlat,xcen,ycen
c
      zpir18 = 2.0*asin(1.0)/180.
c
      xca = xcen*zpir18
      yca = ycen*zpir18
      x3  = rlon*zpir18
      y3  = rlat*zpir18
c
      call sphrot(x2,y2,x3,y3,1,xca,yca,-1)
c
      zsyca = sin(yca)
      zcyca = cos(yca)
c
      zsxsph = sin(x2)
      zcxsph = cos(x2)
      zsysph = sin(y2)
      zcysph = cos(y2)
      zxmxc  = x2 - xca
      zsxmxc = sin(zxmxc)
      zcxmxc = cos(zxmxc)
      zsxrot = sin(x3)
      zcxrot = cos(x3)
      zsyrot = sin(y3)
      zcyrot = cos(y3)
      za1 = zcxmxc*zcxrot + zcyca*zsxmxc*zsxrot
      za2 = zcyca*zsxmxc*zcxrot*zsyrot + zsyca*zsxmxc*zcyrot -
     +      zcxmxc*zsxrot*zsyrot
      za3 =-zsyca*zsxrot/zcysph
      za4 = (zcyca*zcyrot - zsyca*zcxrot*zsyrot)/zcysph
c
      dda = 45.
      ua  = -1.
      va  = -1.
c
      u = za1*ua + za2*va
      v = za3*ua + za4*va
c
      dd=270.-atan2(v,u)/zpir18
      if(dd.gt.+180.) dd=dd-360.
      if(dd.lt.-180.) dd=dd+360.
c
      delta=dd-dda
      if(delta.gt.+180.) delta=delta-360.
      if(delta.lt.-180.) delta=delta+360.
c
      return
      end
c
c
c
      subroutine sphrot(xsph,ysph,xrot,yrot,n,
     +                  xcen,ycen,icall)
c
c  conversion between spherical (xsph,ysph) and spherical rotated
c  (xrot,yrot) coordinates. (xcen,ycen) is the position of the
c  rotated equator/greenwich in terms of (longitude,latitude).
c  all values are given in radians.
c
      real xsph(n), ysph(n), xrot(n), yrot(n)
c
      zsycen = sin(ycen)
      zcycen = cos(ycen)
c
      if (icall.eq.1) then
c
c  compute spherical rotated coordinates as function of
c  spherical coordinates
c
      do j = 1,n
         zxmxc  = xsph(j) - xcen
         zsxmxc = sin(zxmxc)
         zcxmxc = cos(zxmxc)
         zsysph = sin(ysph(j))
         zcysph = cos(ysph(j))
         zsyrot = zcycen*zsysph - zsycen*zcysph*zcxmxc
         zsyrot = max(zsyrot,-1.0)
         zsyrot = min(zsyrot,+1.0)
         yrot(j) = asin(zsyrot)
         zcyrot = cos(yrot(j))
         zcxrot = (zcycen*zcysph*zcxmxc +
     +             zsycen*zsysph)/zcyrot
         zcxrot = max(zcxrot,-1.0)
         zcxrot = min(zcxrot,+1.0)
         zsxrot = zcysph*zsxmxc/zcyrot
         xrot(j) = acos(zcxrot)
         if (zsxrot.lt.0.0) xrot(j) = -xrot(j)
      enddo
c
      elseif (icall.eq.-1) then
c
c  compute spherical coordinates as function of
c  spherical rotated coordinates
c
      do j = 1,n
         zsxrot = sin(xrot(j))
         zcxrot = cos(xrot(j))
         zsyrot = sin(yrot(j))
         zcyrot = cos(yrot(j))
         zsysph = zcycen*zsyrot + zsycen*zcyrot*zcxrot
         zsysph = max(zsysph,-1.0)
         zsysph = min(zsysph,+1.0)
         ysph(j) = asin(zsysph)
         zcysph = cos(ysph(j))
         zcxmxc = (zcycen*zcyrot*zcxrot -
     +             zsycen*zsyrot)/zcysph
         zcxmxc = max(zcxmxc,-1.0)
         zcxmxc = min(zcxmxc,+1.0)
         zsxmxc = zcyrot*zsxrot/zcysph
         zxmxc  = acos(zcxmxc)
         if (zsxmxc.lt.0.0) zxmxc = -zxmxc
         xsph(j) = zxmxc + xcen
      enddo
c
      else
c
      write(*,'(1x,''invalid icall in sphrot'')')
      stop
c
      endif
c
      return
      end
