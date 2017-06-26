 	Program solution	
	  
	IMPLICIT NONE
	INTEGER(8), PARAMETER :: nstep=500000
	INTEGER, PARAMETER ::lx=72, ly=72, lz=72, nnd=24, nnd1=18, nkind=3, d=2
	INTEGER :: imove, nmove, ti, good
	INTEGER :: ntot, ntotm, ntotp, ntotc, nntotc, nndx, nnd_2, lxm, lym, lzm, nx, nx0, nx2, nn1, renx,&
                      modnnd, chnx, nn0, ntote, ntrnt, nacc,&
		     	   ncc,  nnx, nndxx, mii,  nchax, nny, nx1, nxx, ichg, nstar, nend, nnum, nni, nii, nnj, nnk, njj1, njj, mm1,&
			   iia, iib, ichanx, ichanx2, ichaa, ichab,ncc1, ncc2, ichnx1, jcc, jcc2,mmp, jcc1, njcc, njcc1, nccp,icnccp,&
			   ichncc, ichaap, ichabp, iip, nccpp,nccp2, ichnccp, nnb, iistend, nb0, nb,nnd2,nxb,lmove,il
	INTEGER :: i, j, k, i1, i2, i3, ii, jj, kk, kkk, ic, iz, m, n, x1, y1, z1, x2, y2, z2, cbond 
	INTEGER  :: ichapp(nnd+1), mb0(3), mb(3)
	DOUBLE PRECISION :: eab(6,6)
	INTEGER, ALLOCATABLE :: atom(:,:), nna(:,:), nn(:), icha(:), ichaP(:), nchain(:,:), ichb(:),ichbp(:)
	DOUBLE PRECISION, ALLOCATABLE :: bond(:,:)
	DOUBLE PRECISION :: lx_2, ly_2, lz_2, rrx, rry, rrz, rr, xb, yb, zb
	DOUBLE PRECISION :: cons, T0, E, EP, EO, EN, ET,dE, percent, pacc, e1, e2
	INTEGER :: idum, ISEED
 	COMMON /CSEED/ ISEED
	REAL :: ranf
	
	OPEN (UNIT=2,FILE='atom.txt') 

	OPEN (UNIT=3,FILE='chain.txt')  
	OPEN (UNIT=4,FILE='anneal.txt')  

	idum=-29
	!ISEED = 137
	ntot = lx*ly*lz     !!!!!!!!!!!!总格点的个数
	ntotc = int(0.06*ntot/nnd)
	nndx = nnd-nnd1
    !nnd2 = int((nnd+1)/2.D0)
    nnd_2 = int(nnd/2.D0) 
	lx_2 = 0.5D0 * lx
	ly_2 = 0.5D0 * ly
	lz_2 = 0.5D0 * lz
	lxm = 2.d0*lx
	lym = 2.d0*ly
	lzm = 2.d0*lz
	ntotm = lxm*lym*lzm    !!!!!!!!!!!!包括中点在内的格点的总数目
	
	ALLOCATE (atom(ntot,3), nna(ntot,18), nn(ntot), icha(ntot), ichaP(ntot), nchain(ntotc,nnd)) 
    ALLOCATE (bond(ntotm,3), ichb(ntotm), ichbp(ntotm))
	WRITE(*,*)'initializing start'
	
	do i = 1,nkind  
		do j = 1, nkind
			eab(i,j) = 0.0
		end do
	end do			  
	eab(1,2)=1.0
	eab(2,1)=eab(1,2)
	eab(1,3)=1.0
	eab(3,1)=eab(1,3)
	eab(2,3)=-1.0
	eab(3,2)=eab(2,3)

	atom(1,1)=0.0
	atom(1,2)=0.0
	atom(1,3)=0.0
	kkk=0
	  do i1=1,lx
		do i2=1,ly	 
			do i3=1,lz
				kkk=kkk+1
				atom(kkk,1)=atom(1,1)+(i1-1)			
				atom(kkk,2)=atom(1,2)+(i2-1)
				atom(kkk,3)=atom(1,3)+(i3-1)
			end do
		end do
	end do
	bond(1,1)=0.0
	bond(1,2)=0.0
	bond(1,3)=0.0
	kkk=0
	do i1=1,lxm
		do i2=1,lym	 
			do i3=1,lzm
				kkk=kkk+1
				bond(kkk,1)=bond(1,1)+(i1-1)*0.5d0			
				bond(kkk,2)=bond(1,2)+(i2-1)*0.5d0    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!计算键的个数
				bond(kkk,3)=bond(1,3)+(i3-1)*0.5d0
			end do
		end do
	end do
	
	do i=1,ntot
		icha(i)=3
	end do	
	do i=1,ntotm      
		ichb(i)=0
	end do	

	
	do i=1,ntot
  		 nn(i)=0
	end do
   	 
	nna(ntot,18) = 0
	do 11 i=1,ntot
		j=i
		do  while (j.lt.ntot)
			j=j+1
			rrx=atom(i,1)-atom(j,1)
			rry=atom(i,2)-atom(j,2)
			rrz=atom(i,3)-atom(j,3)
		
			if(rrx.gt.lx_2)then
				rrx=-lx+rrx
			else if(rrx.lt.-lx_2)then
				rrx=lx+rrx
			endif
			
			if(rry.gt.ly_2)then
				rry=-ly+rry
			else if(rry.lt.-ly_2)then
				rry=ly+rry
			endif
			if(rrz.gt.lz_2)then
				rrz=-lz+rrz
			else if(rrz.lt.-lz_2)then
				rrz=lz+rrz
			endif
			
			rr=rrx*rrx+rry*rry+rrz*rrz
			
			if(rr.LT. 2.5)then

				nn(i)=nn(i)+1
				nn(j)=nn(j)+1
				nna(i,nn(i))=j
				nna(j,nn(j))=i
			endif
		enddo 
11	continue
	do i = 1, ntotc
		iz = int(lz/nnd_2)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!判断在一个平面内可以排几层链
		m = int((i-1)/(ly/(d+2)*iz))!!!!!!!!!!!!!!!!!!!链子可以排几个平面
		n = mod(i,ly/4)    !!!!!!!!! 在一层内可以排个链子
		
		if (n.eq.0) then
			n = ly/2
		end if
		nx0 = int((i-m*ly/(d+2)*iz-1)/(ly/(d+2)))*nnd_2+m*ly*lz      !!!!!!!!!!!!!!!!!链子底部的第一个点
		do j=1, nnd_2
			if (n .eq. 1) then
				nx= 1+nx0+j-1
			else	
				nx= 1+nx0+(n-1)*(2+d)*lz+j-1
			end if
			nchain(i,j)=nx
		
		end do    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!左侧的链子上的点由下往上排
		do j=1, nnd_2
			if (n .eq. 1) then
				nx= 1+nx0+lz+j-1
			else
				nx= 1+nx0+(n-1)*(2+d)*lz+lz+j-1     !!!!!!!!!表示左侧和右侧相差一个lz的间距
			end if
			nchain(i,nnd+1-j)=nx
		
		end do    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!右侧链子上的点由下往上排
!	write(*,*)i, (nchain(i,j),j=1,nnd)
!!!!!!!!!!!!!以上是排链子的过程
	end do
	do i=1,ntotc
		do j=1,nnd1
			icha(nchain(i,j))=1
		end do
		do j=nnd1+1,nnd
			icha(nchain(i,j))=2
		end do
	end do   	
    percent=1.d0*ntotc*nnd/ntot     !!!!!!!!!!表示链的浓度
    write(*,*)ntotc, percent 

	do i=1,ntotc
		do j=1,nnd
			x1 = atom(nchain(i,j),1)
			y1 = atom(nchain(i,j),2)
			z1 = atom(nchain(i,j),3)
			if(j.eq.nnd)then
				nx0=nchain(i,1)
			else
				nx0=nchain(i,j+1)
			end if
			x2 = atom(nx0,1)
			y2 = atom(nx0,2)
			z2 = atom(nx0,3)
			if(x1 .eq.0 .and. x2 .eq. lx-1) then   !如果键连接的两个单体坐标（x,y,z），位于盒子的两个边界上,把位于0边界的单体加周期边界条件，再计算键的坐标位置。
			   x1 = lx               
	
			else if(x2 .eq.0 .and. x1 .eq. lx-1)then
				x2 = lx
	
			end if
			if(y1 .eq.0 .and. y2 .eq. ly-1) then  
			   y1 = ly               
	
			else if(y2 .eq.0 .and. y1 .eq. ly-1)then
				y2 = ly
	
			end if
			if(z1 .eq.0 .and. z2 .eq. lz-1) then  
			   z1 = lz               
	
			else if(z2 .eq.0 .and. z1 .eq. lz-1)then
				z2 = lz
	
			end if					
			xb = x1+x2
			yb = y1+y2
			zb = z1+z2
			nxb = 1+zb + lzm*yb + lzm*lym*xb    
			ichb(nxb) = 1
		end do
	end do
	
	E=0.
	do i=1,ntot
		nnb = nn(i)
		do j=1,nnb
			ic=icha(nna(i,j))
			if(ICHA(i).lt.ic)then
				E=E+eab(icha(i),ic)
			endif
		end do
	end do
    write(*,*)'E0=', E
	
	do i1=1,ntot
		  ichap(i1)=icha(i1)
	end do
	
	WRITE(*,*)'initializing done'	 
	nb0=0
	do i=1,ntotm
	   if (ichb(i)==1)then
		nb0=nb0+1
		end if
	end do
!goto 9999
!  MC MOVE	
	EO=1D7
	EN=E
	lmove=ntotc*nnd
	T0=70.
	
	do 999 ti =1,75
	!	if(EO-EN.gt.500.)then
			T0=T0*0.94d0
	!	else
		!	T0=T0*0.92d0
	!	endif
		EO=EN
		EP=0.
	!	imove=0
		nacc=0
   		cons=1./T0
		do 777 i=1,nstep
		do 888 il=1,lmove
			ncc=1
! hoping start			
			nnx=INT(1.0+ntotc*ranf(idum))
			nxx=INT(1.0+nnd*ranf(idum))
			nx=nchain(nnx,nxx)
			nnb = nn(nx)
			nx1=int(1.0+nnb*ranf(idum))
			nx2=nna(nx,nx1)
			if(icha(nx2) .eq. 3)then
				ichg=1
				if (nxx.eq.1)then
					nstar=nnd
					nend=nxx+1
				else if(nxx.eq.nnd)then
					nstar=nxx-1
					nend=1
				else
					nstar=nxx-1
					nend=nxx+1
				endif
 !*****************记录键的原始位置************************
				k = 0
				do ii=nstar,nend,nend-nstar 
					k = k+1
					nni=nchain(nnx,ii)
					x1 = atom(nx,1)
					y1 = atom(nx,2)
					z1 = atom(nx,3)
					x2 = atom(nni,1)
					y2 = atom(nni,2)
					z2 = atom(nni,3)
					if(x1 .eq.0 .and. x2 .eq. lx-1) then   !如果键连接的两个单体坐标（x,y,z），位于盒子的两个边界上,把位于0边界的单体加周期边界条件，再计算键的坐标位置。
						 x1 = lx               
					else if(x2 .eq.0 .and. x1 .eq. lx-1)then
						 x2 = lx
					end if
					if(y1 .eq.0 .and. y2 .eq. ly-1) then  
						y1 = ly               
					else if(y2 .eq.0 .and. y1 .eq. ly-1)then
						y2 = ly
					end if
					if(z1 .eq.0 .and. z2 .eq. lz-1) then  
						z1 = lz               
					else if(z2 .eq.0 .and. z1 .eq. lz-1)then
						z2 = lz
					end if					
					xb = x1+x2
					yb = y1+y2
					zb = z1+z2
					mb0(k) = 1 + zb + lzm*yb + lzm*lym*xb 
				end do	
 !***********************************************************	    				
 !***************判断链是否断掉******************************				
					
				do 66 ii=nstar,nend,nend-nstar
					iistend=0
					nni=nchain(nnx,ii)
					kk=1
					do while(iistend.eq.0 .and. kk.le.18)
						if(nna(nni,kk) .eq. nx2) then
							iistend=1
						else
							kk=kk+1		
						end if
					end do
					if(iistend.eq.0)then	
						ichg=ichg+1				
						nii=ii					
						njj1=nni
					end if		
66				continue
			
!************************************************************			
!*************如果一侧链断开***********************************
				njj=nx
				if(ichg.eq.2)then
					if(nii.eq.nstar) then
						mm1=-1
						mb0(1) = mb0(2) 
					else if(nii.eq.nend) then
						mm1=1
					endif
					cbond = 0
					x1 = atom(nx2,1)
					y1 = atom(nx2,2)
					z1 = atom(nx2,3)
					x2 = atom(nx,1)
					y2 = atom(nx,2)
					z2 = atom(nx,3)
					if(x1 .eq.0 .and. x2 .eq. lx-1) then   !如果键连接的两个单体坐标（x,y,z），位于盒子的两个边界上,把位于0边界的单体加周期边界条件，再计算键的坐标位置。
						x1 = lx               
					else if(x2 .eq.0 .and. x1 .eq. lx-1)then
						x2 = lx
					end if
					if(y1 .eq.0 .and. y2 .eq. ly-1) then  
						y1 = ly               
					else if(y2 .eq.0 .and. y1 .eq. ly-1)then
						y2 = ly
					end if
					if(z1 .eq.0 .and. z2 .eq. lz-1) then  
						z1 = lz               
					else if(z2 .eq.0 .and. z1 .eq. lz-1)then
						z2 = lz
					end if				
					
					xb = x1+x2
					yb = y1+y2
					zb = z1+z2
					mb(1) = 1 + zb + lzm*yb + lzm*lym*xb   
					if(ichb(mb(1)) == 1) then        !判断是否键交叉
						cbond = 1         
					end if	
					if(cbond .eq. 1)then	
						goto 888                      !禁止键交叉
					endif	
					
					if(mm1.eq.-1)then
					   nni = nchain(nnx,nend)		
					else
					   nni = nchain(nnx,nstar)	
					endif
					x1 = atom(nx2,1)
					y1 = atom(nx2,2)
					z1 = atom(nx2,3)
					x2 = atom(nni,1)
					y2 = atom(nni,2)
					z2 = atom(nni,3)
					if(x1 .eq.0 .and. x2 .eq. lx-1) then   !如果键连接的两个单体坐标（x,y,z），位于盒子的两个边界上,把位于0边界的单体加周期边界条件，再计算键的坐标位置。
						x1 = lx               
					else if(x2 .eq.0 .and. x1 .eq. lx-1)then
						x2 = lx
					end if
					if(y1 .eq.0 .and. y2 .eq. ly-1) then  
						y1 = ly               
					else if(y2 .eq.0 .and. y1 .eq. ly-1)then
						y2 = ly
					end if
					if(z1 .eq.0 .and. z2 .eq. lz-1) then  
						z1 = lz               
					else if(z2 .eq.0 .and. z1 .eq. lz-1)then
						z2 = lz
					end if			
					xb = x1+x2
					yb = y1+y2
					zb = z1+z2
					mb(2) = 1 + zb + lzm*yb + lzm*lym*xb  
					if(ichb(mb(2)) == 1) then        !判断是否键交叉
						cbond = 1         
					end if	
					if(cbond .eq. 1)then	
						goto 888                      !禁止键交叉
					endif	
					ncc=1
					do while(ichg.eq.2 .and. ncc.ne.nnd)
						ncc=ncc+1
						nii=nii+mm1
						do while (nii.eq.0 .or. nii.eq.nnd+1)
							nii=nii-mm1*nnd
						end do
						if(nii.eq.nxx)then
							nni=nx2
						else
 							nni=nchain(nnx,nii)
						end if
						kk=1
						do while(ichg.ne.1.and.kk.le.18)
							if(nna(nni,kk).eq.njj)then		
								ichg=1
							else
								kk=kk+1
							endif
						enddo
						njj=njj1
						njj1=nni
					enddo
					if(ncc.eq.nnd)then		
						if (ichg.eq.2)then	
							goto 888
						endif
					endif	
					nni = nxx+mm1*ncc-mm1*2
					nnj = nni+mm1
					nnk = nnj+mm1
					if(nni.le.0 .or. nni.ge.(nnd+1))then
							nni=nni-mm1*nnd 
					end if
					if(nnj.le.0 .or. nnj.ge.(nnd+1))then
							nnj=nnj-mm1*nnd 
					end if
					if(nnk.le.0 .or. nnk.ge.(nnd+1))then
							nnk=nnk-mm1*nnd 
					end if
				

!***************记录键的原始位置******************************						
					x1 = atom(nchain(nnx,nni),1)
					y1 = atom(nchain(nnx,nni),2)
					z1 = atom(nchain(nnx,nni),3)
					x2 = atom(nchain(nnx,nnj),1)
					y2 = atom(nchain(nnx,nnj),2)
					z2 = atom(nchain(nnx,nnj),3)
					if(x1 .eq.0 .and. x2 .eq. lx-1) then   !如果键连接的两个单体坐标（x,y,z），位于盒子的两个边界上,把位于0边界的单体加周期边界条件，再计算键的坐标位置。
						x1 = lx               
					else if(x2 .eq.0 .and. x1 .eq. lx-1)then
						x2 = lx
					end if
					if(y1 .eq.0 .and. y2 .eq. ly-1) then  
						y1 = ly               
					else if(y2 .eq.0 .and. y1 .eq. ly-1)then
						y2 = ly
					end if
					if(z1 .eq.0 .and. z2 .eq. lz-1) then  
						z1 = lz               
					else if(z2 .eq.0 .and. z1 .eq. lz-1)then
						z2 = lz
					end if		
					xb = x1+x2
					yb = y1+y2
					zb = z1+z2
					mb0(2) = 1 + zb + lzm*yb + lzm*lym*xb 
					x1 = atom(nchain(nnx,nnk),1)
					y1 = atom(nchain(nnx,nnk),2)
					z1 = atom(nchain(nnx,nnk),3)
					x2 = atom(nchain(nnx,nnj),1)
					y2 = atom(nchain(nnx,nnj),2)
					z2 = atom(nchain(nnx,nnj),3)
					if(x1 .eq.0 .and. x2 .eq. lx-1) then   !如果键连接的两个单体坐标（x,y,z），位于盒子的两个边界上,把位于0边界的单体加周期边界条件，再计算键的坐标位置。
						x1 = lx               
					else if(x2 .eq.0 .and. x1 .eq. lx-1)then
						x2 = lx
					end if
					if(y1 .eq.0 .and. y2 .eq. ly-1) then  
						y1 = ly               
					else if(y2 .eq.0 .and. y1 .eq. ly-1)then
						y2 = ly
					end if
					if(z1 .eq.0 .and. z2 .eq. lz-1) then  
						z1 = lz               
					else if(z2 .eq.0 .and. z1 .eq. lz-1)then
						z2 = lz
					end if				
					xb = x1+x2
					yb = y1+y2
					zb = z1+z2
					mb0(3) = 1 + zb + lzm*yb + lzm*lym*xb 		
						
!***************************************************************	

!***********************键的新位置******************************

				
				
					x1 = atom(nchain(nnx,nni),1)
					y1 = atom(nchain(nnx,nni),2)
					z1 = atom(nchain(nnx,nni),3)
					if(ncc .eq. nnd)then
						x2 = atom(nx2,1)
						y2 = atom(nx2,2)
						z2 = atom(nx2,3)
					else
						x2 = atom(nchain(nnx,nnk),1)
						y2 = atom(nchain(nnx,nnk),2)
						z2 = atom(nchain(nnx,nnk),3)
					end if
					if(x1 .eq.0 .and. x2 .eq. lx-1) then   !如果键连接的两个单体坐标（x,y,z），位于盒子的两个边界上,把位于0边界的单体加周期边界条件，再计算键的坐标位置。
						x1 = lx               
					else if(x2 .eq.0 .and. x1 .eq. lx-1)then
						x2 = lx
					end if
					if(y1 .eq.0 .and. y2 .eq. ly-1) then  
						y1 = ly               
					else if(y2 .eq.0 .and. y1 .eq. ly-1)then
						y2 = ly
					end if
					if(z1 .eq.0 .and. z2 .eq. lz-1) then  
						z1 = lz               
					else if(z2 .eq.0 .and. z1 .eq. lz-1)then
						z2 = lz
					end if	
					xb = x1+x2
					yb = y1+y2
					zb = z1+z2
					mb(3) = 1 + zb + lzm*yb + lzm*lym*xb     
					if(ichb(mb(3)) == 1) then        !判断是否键交叉
						cbond = 1        
					end if	
					if(cbond .eq. 1 .or. mb(3) .eq. mb(1) .or. mb(3) .eq. mb(2))then	
						goto 888                      !禁止键交叉
					endif	
				    if(ncc .eq. nnd)then
						mb(2) = mb(3)
				    end if		
				else
!***********************键的新位置******************************
					!判断是否键交叉
					cbond = 0 
					k = 0
					do ii=nstar,nend,nend-nstar
						k = k + 1
						nni=nchain(nnx,ii)
						x1 = atom(nx2,1)
						y1 = atom(nx2,2)
						z1 = atom(nx2,3)
						x2 = atom(nni,1)
						y2 = atom(nni,2)
						z2 = atom(nni,3)
						if(x1 .eq.0 .and. x2 .eq. lx-1) then  
							x1 = lx               
						else if(x2 .eq.0 .and. x1 .eq. lx-1)then
							x2 = lx
						end if
						if(y1 .eq.0 .and. y2 .eq. ly-1) then  
							y1 = ly               
						else if(y2 .eq.0 .and. y1 .eq. ly-1)then
							y2 = ly
						end if
						if(z1 .eq.0 .and. z2 .eq. lz-1) then  
							z1 = lz               
						else if(z2 .eq.0 .and. z1 .eq. lz-1)then
							z2 = lz
						end if		
						xb = x1+x2
						yb = y1+y2
						zb = z1+z2
						mb(k) = 1 + zb + lzm*yb + lzm*lym*xb      
						if(ichb(mb(k)) == 1) then
							cbond = 1         !cbond: crossing bond
						end if	
					end do
					if (cbond.eq.1)then	              !禁止键交叉
						goto 888
					endif		
				end if	
!****************************************************************		

!*************如果符合键的限制（包括直接交换与蛇形运动后两种情况）**********************************			
						
				if(ichg.eq.1)then
    				dE=0.
					if(ncc.eq.1)then         !情况1：单体与空位直接交换，无蛇形运动
						nnb = nn(nx)
						do 77 ii=1,nnb
							iia=nna(nx,ii)
							iib=nna(nx2,ii)
							ichanx=icha(nx)
							ichanx2=icha(nx2)
							if(iia.ne.nx2)then   
								ichaa=icha(iia)
								ichaap=ichap(iia)
								de=de+eab(ichanx2,ichaap)-eab(ichanx,ichaa)
							endif
							if(iib.ne.nx)then   
								ichab=icha(iib)
								ichabp=ichap(iib)
								de=de+eab(ichanx,ichabp)-eab(ichanx2,ichab)
							endif
77						continue
												
						good = 1
						if (dE > 0)  then
							pacc=exp(-de*cons)
							if (Pacc < ranf(idum))  good = 0
						end if
						if (good == 1) then
							E=E+dE	
							icha(nx2)=ichanx
							icha(nx)=ICHANX2
							ichap(nx2)=ichanx
							ichap(nx)=ICHANX2
							nacc=nacc+ncc
							nchain(nnx,nxx)=nx2
							do k =1, 2           !更新键中点的属性
								ichb(mb0(k))=0
								ichb(mb(k))=1
							end do
						endif
					else 
						ncc1=ncc-1              !情况2：蛇形运动
						ncc2=ncc+1
						nx1=nx
						ichnx1=ichap(nx1)
						ichapp(1)=nx2
						mmp=nxx
						DO 28 JCC=2,NCC2
							ichapp(jcc)=nx1
							if(jcc.ne.ncc2)then
								mmp=mmp+mm1
								if(mmp.eq.0 .or. mmp.eq.(nnd+1))then
									mmp=mmp-mm1*nnd 
								end if
								nx1=NCHAIN(NNX,mmP)
							endif
28						continue
						DO 38 JCC=1,NCC
							jcc1=jcc+1
							njcc=ichapp(jcc)
							njcc1=ichapp(jcc1)
							ichap(njcc)=icha(njcc1)
38						continue
						ichap(ichapp(ncc2))=ICHA(NX2)
						e1=0                        !新的能量算法
						e2=0
						do 48 jcc=1,ncc2
							nccp=ichapp(jcc)
							do  ii=1,18
								iia=nna(nccp,ii)
								e1=e1+eab(icha(nccp),icha(iia))
								e2=e2+eab(ichap(nccp),ichap(iia))
							end do
48						continue
						do 58 jcc=1,ncc
							nccp=ichapp(jcc)
							jcc2=jcc+1
							do while(jcc2.le.ncc2)
								nccp2=ichapp(jcc2)
								do ii=1,18
									if(nna(nccp2,ii).eq.nccp)then
										e1=e1-eab(icha(nccp2),icha(nccp))
										e2=e2-eab(ichap(nccp2),ichap(nccp))
									endif
								end do
								jcc2=jcc2+1
							enddo
58						continue
						
						de=e2-e1

						mmp=nxx
						good = 1
						if (dE > 0)  then
							pacc=exp(-de*cons)
							if (Pacc < ranf(idum))  good = 0
						end if
						if (good == 1) then
							e=e+de	
!	WRITE(*,*)E, DE
							DO 78 JCC=1,NCC

								nccp=ichapp(jcc)
								icha(nccp)=ichap(nccp)
								nchain(nnx,mmp)=nccp
								mmp=mmp+mm1
								if(mmp.eq.0 .or. mmp.eq.(nnd+1))then
									mmp=mmp-mm1*nnd
								end if
78							continue
							nccp=ichapp(ncc2)
							icha(nccp)=ichap(nccp)
							nacc=nacc+ncc
							
							if(ncc .eq. nnd)then
								do k =1, 2
									ichb(mb0(k))=0
									ichb(mb(k))=1
								end do
							else			    	
								do k =1, 3           !更新键中点的属性
									ichb(mb0(k))=0
									ichb(mb(k))=1
						 		end do
							end if
						else
							DO 99 JCC=1,NCC2
								njcc=ichapp(jcc)
								ichap(njcc)=icha(njcc)
99							continue
						endif
 					endif
				endif
			endif
			EP=EP+E
		!	imove=imove+1
888	    end do
777		enddo	 
		EP=EP/lmove
		EP=EP/nstep
		write(*,*)'ti=', ti, 'nacc=', nacc, 'E=', E ,'lmove=', lmove,'nstep=',nstep
		write(4,*)ti, nacc, cons , EP
	!	EN=EP
999 continue	
 
!**************** check ****************** 
8888 ET=0.
	do i=1,ntot
		nnb = nn(i)
		do j=1,nnb
			ic=icha(nna(i,j))
			if(ICHA(i).le.ic)then
				ET=ET+eab(icha(i),ic)
			endif
		end do
	end do
    write(*,*)'E0=',E,'ET=',ET
	if(E-et.gt.1.0.or.et-E.gt.1.0)then
	 endif
 
 	 do 9091 i=1,ntotc
		do 9092 j=1,nnd
		    if(j.eq.nnd)then
			rrx=atom(nchain(i,nnd),1)-atom(nchain(i,1),1)
			rry=atom(nchain(i,nnd),2)-atom(nchain(i,1),2)
			rrz=atom(nchain(i,nnd),3)-atom(nchain(i,1),3)
			else
			rrx=atom(nchain(i,j),1)-atom(nchain(i,j+1),1)
			rry=atom(nchain(i,j),2)-atom(nchain(i,j+1),2)
			rrz=atom(nchain(i,j),3)-atom(nchain(i,j+1),3)
			endif

			if(rrx.gt.lx_2) then
	          rrx=-lx+rrx
	        else if(rrx.lt.-lx_2) then
	          rrx=lx+rrx
	        end if
	        if(rry.gt.ly_2) then
	          rry=-ly+rry
	        else if(rry.lt.-ly_2) then
	          rry=ly+rry
	        end if
			if(rrz.gt.lz_2) then
	          rrz=-lz+rrz
	        else if(rrz.lt.-lz_2) then
	          rrz=lz+rrz
	        end if
			rr=rrx*rrx+rry*rry+rrz*rrz

			if(rr.gt.2.5)then
				write(*,*) 'error',i,j,j+1
			endif
9092	continue
9091 continue
	nb=0
	do i=1,ntotm
	   if (ichb(i)==1)then
		nb=nb+1
		end if
	end do
	if(nb .ne. nb0)then
		write(*,*)'nb0=',nb0, 'nb=', nb
	end if



write(*,*)'移动无错误,输出移动后结果'


9999	do kkk=1,ntot
   		write(2,*)atom(kkk,1),atom(kkk,2),atom(kkk,3),icha(kkk),ichap(kkk)
	end do
	do kkk=1,ntotc
		write(3,*)(nchain(kkk,j),j=1,nnd)
	end do	 
	write(*,*)'e=',e
	e=0.
	do i=1,ntot
		nnb = nn(i)
		do j=1,nnb
			ic=icha(nna(i,j))
			if(ICHA(i).lt.ic)then
 				e=e+eab(icha(i),ic)
			endif
		end do	
	end do
    write(*,*)'ep=',e 

END

 
 !  DOUBLE PRECISION FUNCTION RANF()
 !  DATA IA/16807/,IC/2147483647/, IQ/127773/,IR/2836/
 !  COMMON /CSEED/ ISEED
 !  IH=ISEED/IQ
 !  IL=MOD(ISEED,IQ)
 !  IT=IA*IL-IR*IH
 !  IF(IT.GT.0) THEN
 !  ISEED=IT
 !  ELSE
 !  ISEED=IC+IT
 !  END IF
 !  RANF=ISEED/FLOAT(IC)
 !  RETURN 
 !  END
   


!********************************************************************************************************
! "Minimal" random number generator of Park and Miller combined with a Marsaglia shift sequence.
! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the end point values,(0,1)). This 
! fully portable, scalar generator has the "traditional"(not Fortran 90) calling sequence with a
! random deviate as the returned function value: call with idum a negative integer to initialize;
! thereafter, do not alter idum except to reinitialize. The period of this generotor is about 3.1*10^^18.
!*********************************************************************************************************
FUNCTION ranf(idum)
IMPLICIT NONE
INTEGER, PARAMETER :: K4B=selected_int_kind(9)
INTEGER(K4B), INTENT(INOUT) :: idum
REAL :: ranf
INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
REAL, SAVE :: am
INTEGER(K4B), SAVE :: ix=-1, iy=-1, k
if (idum <= 0 .or. iy < 0) then					! Initialize 
	am=nearest(1.0,-1.0)/IM
	iy=ior(ieor(888889999, abs(idum)),1)
	ix=ieor(777755555, abs(idum))
	idum=abs(idum)+1							! Set idum positive
end if
ix=ieor(ix, ishft(ix, 13))						! Marsaglia shift sequence with period 2^^32-1
ix=ieor(ix, ishft(ix, -17))
ix=ieor(ix, ishft(ix, 5))
k=iy/IQ											! Park-Miller sequence by Schrage's method, period 2^^31-1
iy=IA*(iy-k*IQ)-IR*k
if (iy < 0 ) iy=iy+IM
ranf=am*ior(iand(IM, ieor(ix,iy)),1)			! Combine the two generators with masking to ensure nonzero value
END FUNCTION ranf

	 
