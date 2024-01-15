      program photosphere

c v3.0 CCE (May 2017): added ability to switch between KH96 and PM13 
c  PM13: A0 to M5 stars, U through K band 
c  KH96: A0 to M6 stars, U through M band
c  also, changed to find standard from Teff, within 25K

c v2.0, Nuria (Dec 2016): Converts Mbol to M using V
c       and BC (in table A5) for Tef > 5770 and eq.
c       (A1) for cooler stars.
c
c
c	reads Teff and R and calculates photospheric fluxes using table
c	A5 from KH95 at the photometric bands.
c	Then interpolates in the wavelength scale where disk and wall
c	fluxes are calculated. Assumes Rayleigh-Jeans for the long
c wavelength extension


      parameter(nf1=12,nf2=12,ndate=3000,nkh=55,npm=80,nbess=18)
      parameter(nwl0=1000)
      real*4 linter
      character*200 linea,filewave,tableA5kh95,table5pm13
      character*200 fileout, table
      character*2 sistema(ndate),banda(nf1),bandapm(nf1)
      character*2 tipokh(nkh),nada,tipopm(npm)
      character*2 tipoesp
      common/kh95/ tefkh(nkh),bckh(nkh),umvkh(nkh),bmvkh(nkh),
     .  vmrckh(nkh),
     .  vmrjkh(nkh),vmickh(nkh),vmijkh(nkh),vmjkh(nkh),
     .  vmhkh(nkh),vmkkh(nkh),vmlkh(nkh),vmmkh(nkh)
      common/pm13/ tefpm(npm),bcpm(npm),umbpm(npm),bmvpm(npm),
     .  vmrcpm(npm),vmicpm(npm),vmjpm(npm),
     .  vmhpm(npm),vmkpm(npm),
     .  xkmw1pm(npm),xkmw2pm(npm),xkmw3pm(npm),xkmw4pm(npm)
      dimension wef(nf1),zp(nf1),fluxb(nf1)
      dimension wefpm(nf1),zppm(nf1)
      dimension wl(nwl0),flux(nwl0)
      dimension xmo(nf1),pro(nf1)


      data banda/'U','B','V','Rc','Rj','Ic','Ij','J','H',
     .	'K','L','M'/
c	1:U,2:B,3:V,4:Rc,5:Rj,6:Ic,7:Ij,8:J,9:H,10:K,11:L,12:M
	    data wef/0.36,0.44,0.55,0.64,0.7,0.79,0.9,1.22,1.63,2.19,3.45,
     .  4.75/
c       zero points - erg/cm2/s/Hz
c       from Bessell 1979,Johson66(Jsystem),BB88 nir
      data zp/1.81e-20,4.26e-20,3.64e-20,3.08e-20,3.01e-20,2.55e-20,
     .  2.43e-20,1.57e-20,1.02e-20,6.36e-21,2.81e-21,1.54e-21/
     
        data bandapm/'U','B','V','Rc','Ic','J','H',
     .	'K','W1','W2','W3','W4'/
c	1:U,2:B,3:V,4:Rc,5:Ic,6:J,7:H,8:K,9:W1,10:W2,11:W4,12:W4
	    data wefpm/0.36,0.44,0.55,0.64,0.79,1.22,1.63,2.19,3.4,4.6,
     .  12,22/
c       zero points - erg/cm2/s/Hz
c       from Bessell 1979,BB88 nir, Jarrett et al. 2011 WISE 
c		for Jarrett convert Flam to Fnu i.e., (*lam^2 in micron/c in micron)*1e7 (W->erg/s)
      data zppm/1.81e-20,4.26e-20,3.64e-20,3.08e-20,2.55e-20,
     .  1.57e-20,1.02e-20,6.36e-21,3.10e-21,1.72e-21,3.17e-22,8.36e-23/
        
      data xmsol/1.99e33/,xlsol/3.83e33/,rsol/6.96e10/
      data spa/3.15e7/,c/2.99793e10/,gg/6.67e-8/
      data sigmar/5.6e-5/
      data h/6.6262E-27/,bk/1.38062E-16/
      data pc/3.09e18/,au/1.496e13/
      pi=acos(-1.0)

c--------------------------------------------
c	input
c	reads Teff, R, distance in parsecs, and wavelength file name
      read(5,*)teff
      read(5,*)radius
      read(5,*)distance
      read(5,'(a)')table
      read(5,'(a)')filewave
      read(5,'(a)')tableA5kh95
      read(5,'(a)')table5pm13
      read(5,'(a)')fileout
      open(31,file=fileout,status='unknown')
c----------------------------------------------
c	reads wavelengths
      open(50,file=filewave,status='old')
      read(50,*)nwl
      do i=1,nwl
        read(50,*)wl(i)
      enddo
      close(50)
      
c---------------------------------------------
c	luminosity in solar units
      xlum=4.*pi*((radius*rsol)**2/xlsol)*sigmar*teff**4
c	bolometric magnitude
      xmbol=4.75-2.5*alog10(xlum)
      dmodulus=5.*alog10(distance)-5.

c---------------------------------------------
c
	if(table == "kh95") then
	
c	values of corresponding standard
c	interpolate in table A5 of KH95

c	read table A5 de KH95
c       open(unit=52,file='/Users/ncalvet/datos/tabla_a5_kh95',
      open(unit=52,file=tableA5kh95,
     .	status='old')
c	solo a partir de A0, para tener al menos L, y no lee N,12
      nskip=11
      do i=1,nskip
        read(52,401)linea
      enddo
401   format(a)
      do i=1,nkh-nskip+1
        read(52,400,end=1000)tipokh(i),nada,itefkh,bckh(i),umvkh(i),
     .  bmvkh(i),vmrckh(i),
     .  vmrjkh(i),vmickh(i),vmijkh(i),vmjkh(i),vmhkh(i),vmkkh(i),
     .	vmlkh(i)
     . 	,vmmkh(i)

400   format(a2,a2,i6,12f6.2)
      tefkh(i)=itefkh
      enddo
1000	continue
      nkh1=nkh-nskip+1
      close(52)

      iesp=0
      diff=1.0e10
c	find closest standard from Teff
      do ik=1,nkh1
        if (abs(teff - tefkh(ik)).lt.diff) then
          iesp=ik
          diff = abs(teff - tefkh(ik))
        endif
      enddo


c	spectral type and intrinsic colors U - L KH95
      tipoesp=tipopm(iesp)
      umv0=umvkh(iesp)
      bmv0=bmvkh(iesp)
      umb0=umv0-bmv0
      vmr0j=vmrjkh(iesp)
      vmr0c=vmrckh(iesp)
      vmi0j=vmijkh(iesp)
      vmi0c=vmickh(iesp)
      vmj0=vmjkh(iesp)
      vmh0=vmhkh(iesp)
      vmk0=vmkkh(iesp)
      vml0=vmlkh(iesp)
      vmm0=vmmkh(iesp)
c      write(*,*)umv0,bmv0,vmr0j,vmr0c,vmi0j
      
c v2 using Tef of G5 in Table 5 of KH95 to separate
      if(teff.lt.5770.) then
c	Mj, using eq A1 in KH95
        xmj=xmbol-0.10*vmk0-1.17
c	J magnitude
        xj=xmj+dmodulus
c       magnitudes scaled at J
        xmo(8)=xj
        xmo(3)=vmj0+xj
      else
c       Mv, using BC in Table A5 of KH95
        xmv=xmbol-bckh(iesp)
        v=xmv+dmodulus
c       magnitudes scaled at V
        xmo(3) = v
        xmo(8) = v - vmj0
      endif

      xmo(2)=bmv0+xmo(3)
      xmo(1)=umb0+xmo(2)
      xmo(4)=xmo(3)-vmr0c
      xmo(5)=xmo(3)-vmr0j
      xmo(6)=xmo(3)-vmi0c
      xmo(7)=xmo(3)-vmi0j
      xmo(9)=xmo(3)-vmh0
      xmo(10)=xmo(3)-vmk0
      xmo(11)=xmo(3)-vml0
      xmo(12)=xmo(3)-vmm0
c     write(*,*)bmv0
      


c	Average magnitudes
      do if=1,nf1
        pro(if)=xmo(if)
c        write(*,*)pro(if)
      enddo

c	write(*,*)(pro(i),i=1,nf1)

      nf=nf1
      do i=1,nf
        fnu=zp(i)*10.**(-pro(i)/2.5)
        xnu=c*1e4/wef(i)
        wangs=wef(i)*1.e4
c	flambda per A
        fl=xnu*fnu/wangs
        if (i.eq.1)flu=fl
        fluxb(i)=alog10(xnu*fnu)
c	write(*,*)wef(i)*1.e4,xnu*fnu
      enddo

c	photosphere in adopted wl scale
      i0=1
      do i=1,nwl
        if(wl(i).lt.wef(nf1)) then
c	interpolate in the log
          flux(i)=linter(wl(i),wef,fluxb,nf1,i0)
          flux(i)=10.**flux(i)
        else
          flux(i)=10.**fluxb(nf1)*(wef(nf1)/wl(i))**3
        endif
c	flux(i)=alog10(flux(i))
        write(31,*)wl(i),flux(i)
      enddo
      
      endif
c
c---------------------------------------------

c---------------------------------------------
c
	if(table == "pm13") then
	
c	values of corresponding standard
c	interpolate in table 5 of PM13

c	read table 5 of PM13
      open(unit=52,file=table5pm13,
     .	status='old')

c	reading F0 to M5 
      nskip=51
      do i=1,nskip
        read(52,402)linea
      enddo

402   format(a)
      do i=1,npm-nskip+1
        read(52,403,end=1001)tipopm(i),nada,itefpm,bcpm(i),umbpm(i),
     .  bmvpm(i),vmrcpm(i),vmicpm(i),
     .  vmjpm(i),vmhpm(i),vmkpm(i),
     .	xkmw1pm(i),xkmw2pm(i),xkmw3pm(i),xkmw4pm(i)

c403   format(a2,a1,i8,1f6.2,7f7.3)
403   format(a2,a1,i8,1f6.2,7f7.3,2f6.3,1f7.3,1f6.3)
        tefpm(i)=itefpm    
        
c        write(*,*)tipopm(i),nada,itefpm,bcpm(i),umbpm(i),
c     .  bmvpm(i),vmrcpm(i),vmicpm(i),
c     .  vmjpm(i),vmhpm(i),vmkpm(i)
          
      enddo
1001	continue
      npm1=npm-nskip+1
      close(52)
      
      iesp=0
      diff=1.0e10
c	find closest standard from Teff
      do ik=1,npm1
        if (abs(teff - tefpm(ik)).lt.diff) then
          iesp=ik
          diff = abs(teff - tefpm(ik))
        endif
      enddo

c	spectral type and intrinsic colors U - K PM13
      tipoesp=tipopm(iesp)
      umb0=umbpm(iesp)
      bmv0=bmvpm(iesp)
      vmr0c=vmrcpm(iesp)
      vmi0c=vmicpm(iesp)
      vmj0=vmjpm(iesp)
      vmh0=vmhpm(iesp)
      vmk0=vmkpm(iesp)
      xkmw10=xkmw1pm(iesp)
      xkmw20=xkmw2pm(iesp)
      xkmw30=xkmw3pm(iesp)
      xkmw40=xkmw4pm(iesp)
      
c v2 using Tef of G5 in Table 5 of PM13 to separate
      if(teff.lt.5660.) then
c	Mj, using eq A1 in KH95
        xmj=xmbol-0.10*vmk0-1.17
c	J magnitude
        xj=xmj+dmodulus
c       magnitudes scaled at J
        xmo(6)=xj
        xmo(3)=vmj0+xj
      else
c       Mv, using BC in Table 5 of PM13
        xmv=xmbol-bcpm(iesp)
        v=xmv+dmodulus
c       magnitudes scaled at V
        xmo(3) = v
        xmo(6) = v - vmj0
      endif

      xmo(2)=bmv0+xmo(3)
      xmo(1)=umb0+xmo(2)
      xmo(4)=xmo(3)-vmr0c
      xmo(5)=xmo(3)-vmi0c
      xmo(7)=xmo(3)-vmh0
      xmo(8)=xmo(3)-vmk0
      xmo(9)=xmo(8)-xkmw10
      xmo(10)=xmo(8)-xkmw20
      xmo(11)=xmo(8)-xkmw30
      xmo(12)=xmo(8)-xkmw40


c	Average magnitudes
      do if=1,nf2
        pro(if)=xmo(if)
      enddo

c	write(*,*)(pro(i),i=1,nf2)

      nf=nf2
      do i=1,nf
        fnu=zppm(i)*10.**(-pro(i)/2.5)
        xnu=c*1e4/wefpm(i)
        wangs=wefpm(i)*1.e4
c	flambda per A
        fl=xnu*fnu/wangs
        if (i.eq.1)flu=fl
        fluxb(i)=alog10(xnu*fnu)
c	write(*,*)wefpm(i)*1.e4,xnu*fnu
      enddo

c--------------------------------
c	photosphere in adopted wl scale
      i0=1
      do i=1,nwl
        if(wl(i).lt.wefpm(nf2)) then
c	interpolate in the log
          flux(i)=linter(wl(i),wefpm,fluxb,nf2,i0)
          flux(i)=10.**flux(i)
        else
          flux(i)=10.**fluxb(nf2)*(wefpm(nf2)/wl(i))**3
        endif
c	flux(i)=alog10(flux(i))
        write(31,*)wl(i),flux(i)
      enddo

      endif
c      
c--------------------------------    
  
      end

c---------------------------
      real*4 function linter(x0, x, y, n,i0)
c	linear interpolation at x0 on arrays x and y
      dimension x(n), y(n)
      if (x(2) .lt. x(1)) then
c	x en orden decreciente
        do while (x0 .lt. x(i0))
          i0 = i0 + 1
        end do
      else
c	x en orden creciente
        do while (x0 .gt. x(i0))
          i0 = i0 + 1
        end do
      end if
      if (i0 .gt. 1) i0 = i0 - 1
      linter = ((y(i0) * (x(i0 + 1) - x0)) + (y(i0 + 1) * (x0 - x(i0))))
     & / (x(i0 + 1) - x(i0))
      return
      end
