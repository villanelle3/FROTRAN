module avanco_temporal
use constantes
use termos_rhs
! real(pr),parameter,dimension(6)::alfa=(/ 0.0_pr, -0.737101392796_pr, -1.634740794341_pr, -0.744739003780_pr,  &
! &  -1.469897351522_pr, -2.813971388035_pr /), beta=(/ 0.032918605146_pr, 0.823256998200_pr, 0.381530948900_pr, &
! &   0.200092213184_pr, 1.718581042715_pr, 0.27_pr /), c=(/ 0.0_pr, 0.032918605146_pr, 0.249351723343_pr, & 
! &   0.466911705055_pr, 0.582030414044_pr, 0.847252983783_pr /) !Berland

real(pr),parameter,dimension(6)::alfa=(/ 0.0_pr, -0.691750960670_pr, -1.727127405211_pr, -0.694890150986_pr,  &
&  -1.039942756197_pr, -1.531977447611_pr /), beta=(/ 0.122_pr, 0.477263056358_pr, 0.381941220320_pr, &
&   0.447757195744_pr, 0.498614246822_pr, 0.186648570846_pr /), c=(/ 0.0_pr, 0.122_pr, 0.269115878630_pr, & 
&   0.447717183551_pr, 0.749979795490_pr, 0.898555413085_pr /) !Allampalli 6 estagios

! real(pr),parameter,dimension(7)::alfa=(/ 0.0_pr, -0.647900745934_pr, -2.704760863204_pr, -0.460080550118_pr,  &
! &  -0.500581787785_pr, -1.906532255913_pr, -1.45_pr /), beta=(/ 0.117322146869_pr, 0.503270262127_pr, 0.233663281658_pr, &
! &   0.283419634625_pr, 0.540367414023_pr, 0.371499414620_pr, 0.136670099385_pr /), c=(/ 0.0_pr, 0.117322146869_pr, 0.294523230758_pr, & 
! &   0.305658622131_pr, 0.582864148403_pr, 0.858664273599_pr, 0.868664273599_pr /) !Allampalli 7 estagios

integer(pint)::ii

contains

subroutine RK46DF !Berland et al., 2006
implicit none
complex(pr),dimension(nx,ny)::viscx,viscy,tnlxe,tnlye,ZBx,ZBy,k1x,k1y,RHSx,RHSy
k1x=0.0_pr
k1y=0.0_pr
do ii=1,6 !Runge-Kutta de quarta ordem com 6 passos otimizado (Berland et al., 2006)
	viscx=ni*k2*u !Calculo do termo viscoso na direcao x
	viscy=ni*k2*v !Calculo do termo viscoso na direcao y
	ut=u
	vt=v
	utfis=ufis
	vtfis=vfis
	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
	!t3=((nstep)+c(ii))*dt
	!call fonte

	RHSx=-tnlxe-viscx !+fontex   !-ZBx calculo do lado direito da equacao em x
	RHSy=-tnlye-viscy !+fontey   !-ZBy calculo do lado direito da equacao em y
	call projecao(RHSx,RHSy)

	!coeficientes do RK4
	k1x=alfa(ii)*k1x+dt*RHSx
	k1y=alfa(ii)*k1y+dt*RHSy

	!atualizacao da predicao das velocidades
	u=u+k1x*beta(ii)
	v=v+k1y*beta(ii)
	
	ufis=u
	vfis=v
	call ZFFT2D(ufis,NX,NY,1)
	call ZFFT2D(vfis,NX,NY,1)
	ufis=real(ufis)
	vfis=real(vfis)
end do !fim do laco do Runge-Kutta de quarta ordem com 6 passos otimizado
end subroutine RK46DF

subroutine RK4DF !Berland et al., 2006
implicit none
complex(pr),dimension(nx,ny)::viscx,viscy,tnlxe,tnlye,ZBx,ZBy,k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,knx,kny,RHSx,RHSy
	
	knx=u
	kny=v
	t3=t

	viscx=ni*k2*knx !Calculo do termo viscoso na direcao x
	viscy=ni*k2*kny !Calculo do termo viscoso na direcao y
	ut=knx
	vt=kny
	utfis=ufis
	vtfis=vfis
	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
	call fonte
	!Zona de amortecimentao (Buffer Zone)
	ZBx=phi*(ufis-Qtx)
	ZBy=phi*(vfis-Qty)
	call ZFFT2D(ZBx,NX,NY,-1)
	call ZFFT2D(ZBy,NX,NY,-1)

	RHSx=-tnlxe-viscx-ZBx+fontex   !calculo do lado direito da equacao em x
	RHSy=-tnlye-viscy-ZBy+fontey   !calculo do lado direito da equacao em y
	call projecao(RHSx,RHSy)
	!coeficientes do RK4
	k1x=RHSx
	k1y=RHSy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	knx=u+0.5_pr*dt*k1x
	kny=v+0.5_pr*dt*k1y
	t3=t+0.5_pr*dt
	viscx=ni*k2*knx !Calculo do termo viscoso na direcao x
	viscy=ni*k2*kny !Calculo do termo viscoso na direcao y
	ut=knx
	vt=kny
	utfis=ufis
	vtfis=vfis
	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
	call fonte
	!Zona de amortecimentao (Buffer Zone)
	ZBx=phi*(ufis-Qtx)
	ZBy=phi*(vfis-Qty)
	call ZFFT2D(ZBx,NX,NY,-1)
	call ZFFT2D(ZBy,NX,NY,-1)

	RHSx=-tnlxe-viscx-ZBx+fontex   !calculo do lado direito da equacao em x
	RHSy=-tnlye-viscy-ZBy+fontey   !calculo do lado direito da equacao em y
	call projecao(RHSx,RHSy)
	!coeficientes do RK4
	k2x=RHSx
	k2y=RHSy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	knx=u+0.5_pr*dt*k2x
	kny=v+0.5_pr*dt*k2y
	t3=t+0.5_pr*dt
	viscx=ni*k2*knx !Calculo do termo viscoso na direcao x
	viscy=ni*k2*kny !Calculo do termo viscoso na direcao y
	ut=knx
	vt=kny
	utfis=ufis
	vtfis=vfis
	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
	call fonte
	!Zona de amortecimentao (Buffer Zone)
	ZBx=phi*(ufis-Qtx)
	ZBy=phi*(vfis-Qty)
	call ZFFT2D(ZBx,NX,NY,-1)
	call ZFFT2D(ZBy,NX,NY,-1)

	RHSx=-tnlxe-viscx-ZBx+fontex   !calculo do lado direito da equacao em x
	RHSy=-tnlye-viscy-ZBy+fontey   !calculo do lado direito da equacao em y
	call projecao(RHSx,RHSy)
	!coeficientes do RK4
	k3x=RHSx
	k3y=RHSy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	knx=u+dt*k3x
	kny=v+dt*k3y
	t3=t+dt
	viscx=ni*k2*knx !Calculo do termo viscoso na direcao x
	viscy=ni*k2*kny !Calculo do termo viscoso na direcao y
	ut=knx
	vt=kny
	utfis=ufis
	vtfis=vfis
	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
	call fonte
	!Zona de amortecimentao (Buffer Zone)
	ZBx=phi*(ufis-Qtx)
	ZBy=phi*(vfis-Qty)
	call ZFFT2D(ZBx,NX,NY,-1)
	call ZFFT2D(ZBy,NX,NY,-1)

	RHSx=-tnlxe-viscx-ZBx+fontex   !calculo do lado direito da equacao em x
	RHSy=-tnlye-viscy-ZBy+fontey   !calculo do lado direito da equacao em y
	call projecao(RHSx,RHSy)
	!coeficientes do RK4
	k4x=RHSx
	k4y=RHSy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!atualizacao da predicao das velocidades
	u=u+1.0_pr/6.0_pr*dt*(k1x+2.0_pr*k2x+2.0_pr*k3x+k4x)
	v=v+1.0_pr/6.0_pr*dt*(k1y+2.0_pr*k2y+2.0_pr*k3y+k4y)
	
	ufis=u
	vfis=v
	call ZFFT2D(ufis,NX,NY,1)
	call ZFFT2D(vfis,NX,NY,1)
	ufis=real(ufis)
	vfis=real(vfis)
end subroutine RK4DF

subroutine eulerDF !Berland et al., 2006
implicit none
complex(pr),dimension(nx,ny)::viscx,viscy,tnlxe,tnlye,ZBx,ZBy,RHSx,RHSy

	viscx=ni*k2*u !Calculo do termo viscoso na direcao x
	viscy=ni*k2*v !Calculo do termo viscoso na direcao y
	ut=u
	vt=v
	utfis=ufis
	vtfis=vfis
	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
	t3=t
	call fonte
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Zona de amortecimentao (Buffer Zone)
	ZBx=phi*(ufis-Qtx)
	ZBy=phi*(vfis-Qty)
	call ZFFT2D(ZBx,NX,NY,-1)
	call ZFFT2D(ZBy,NX,NY,-1)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	RHSx=-tnlxe-viscx-ZBx+fontex   !calculo do lado direito da equacao em x
	RHSy=-tnlye-viscy-ZBy+fontey   !calculo do lado direito da equacao em y
	call projecao(RHSx,RHSy)
	
	!atualizacao da predicao das velocidades
	u=u+dt*RHSx
	v=v+dt*RHSy
		
	ufis=u
	vfis=v
	call ZFFT2D(ufis,NX,NY,1)
	call ZFFT2D(vfis,NX,NY,1)
	ufis=real(ufis)
	vfis=real(vfis)
end subroutine eulerDF

subroutine euler
implicit none
complex(pr),dimension(nx,ny)::viscx,viscy,tnlxe,tnlye,ZBx,ZBy,RHSx,RHSy
	viscx=ni*k2*u !Calculo do termo viscoso na direcao x
	viscy=ni*k2*v !Calculo do termo viscoso na direcao y
	ut=u
	vt=v
	utfis=ufis
	vtfis=vfis
	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
	call forca !Calcula a forca da fronteira imersa via MFV (Mariano, 2007)
	RHSx=(-tnlxe+fx)  ! calculo do lado direitao da equacao em x
	RHSy=(-tnlye+fy)  ! calculo do lado direitao da equacao em y
	call projecao(RHSx,RHSy)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Zona de amortecimentao (Buffer Zone)
	ZBx=phi*(ufis-Qtx)
	ZBy=phi*(vfis-Qty)
	
	call ZFFT2D(ZBx,NX,NY,-1)
	call ZFFT2D(ZBy,NX,NY,-1)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	RHSx=(RHSx-viscx-ZBx)  !calculo do lado direitao da equacao em x
	RHSy=(RHSy-viscy-ZBy)  !calculo do lado direitao da equacao em y

	u=u+dt*RHSx
	v=v+dt*RHSy

	call projecao(u,v)
		
	ufis=u
	vfis=v
	call ZFFT2D(ufis,NX,NY,1)
	call ZFFT2D(vfis,NX,NY,1)
	ufis=real(ufis)
	vfis=real(vfis)
end subroutine euler

subroutine RK46 !Berland et al., 2007
implicit none
complex(pr),dimension(nx,ny)::viscx,viscy,tnlxe,tnlye,ZBx,ZBy,k1x,k1y,RHSx,RHSy
k1x=0.0_pr
k1y=0.0_pr
do ii=1,6 !Runge-Kutta de quarta ordem com 6 passos otimizado
	viscx=ni*k2*u !Calculo do termo viscoso na direcao x
	viscy=ni*k2*v !Calculo do termo viscoso na direcao y
	ut=u
	vt=v
	utfis=ufis
	vtfis=vfis
	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
	call forca !Calcula a forca da fronteira imersa via MFV (Mariano, 2007)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Zona de amortecimentao (Buffer Zone)
	ZBx=phi*(ufis-Qtx)
	ZBy=phi*(vfis-Qty)
	call ZFFT2D(ZBx,NX,NY,-1)
	call ZFFT2D(ZBy,NX,NY,-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	RHSx=(-tnlxe-viscx-ZBx+fx)  ! calculo do lado direito da equacao em x
	RHSy=(-tnlye-viscy-ZBy+fy)  ! calculo do lado direito da equacao em y
	call projecao(RHSx,RHSy)

	!coeficientes do RK4
	k1x=alfa(ii)*k1x+dt*RHSx
	k1y=alfa(ii)*k1y+dt*RHSy

	!atualizacao da predicao das velocidades
	u=u+k1x*beta(ii)
	v=v+k1y*beta(ii)
	
	ufis=u
	vfis=v
	call ZFFT2D(ufis,NX,NY,1)
	call ZFFT2D(vfis,NX,NY,1)
	ufis=real(ufis)
	vfis=real(vfis)
end do !Runge-Kutta de quarta ordem com 6 passos otimizado
	call projecao(u,v)
end subroutine RK46

subroutine imp_visc

implicit none
complex(pr),dimension(nx,ny)::tnlxe,tnlye,ZBx,ZBy,RHSx,RHSy
	ut=u
	vt=v
	utfis=ufis
	vtfis=vfis
 	call tnlinear(tnlxe,tnlye) !Calculo do termo nao-linear
 	call forca !Calcula a forca da fronteira imersa via MFV (Mariano, 2007)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	!Zona de amortecimentao (Buffer Zone)
! 	ZBx=phi*(ufis-Qtx)
! 	ZBy=phi*(vfis-Qty)
! 	call ZFFT2D(ZBx,NX,NY,-1)
! 	call ZFFT2D(ZBy,NX,NY,-1)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	RHSx=(-tnlxe+fx)  !-ZBx calculo do lado direitao da equacao em x
 	RHSy=(-tnlye+fy)  !-ZBy calculo do lado direitao da equacao em y
 	call projecao(RHSx,RHSy)

	u=(u+dt*RHSx)*aux
	v=(v+dt*RHSy)*aux
end subroutine imp_visc

subroutine RK46DFBIF !Berland et al., 2006
implicit none
complex(pr),dimension(nx,ny)::tnlxe,tnlye,ZBx,ZBy,k1x,k1y,RHSx,RHSy
k1x=0.0_pr
k1y=0.0_pr
do ii=1,6 !Runge-Kutta de quarta ordem com 6 passos otimizado (Berland et al., 2006)

	ut=u
	vt=v
	utfis=ufis
	vtfis=vfis

	dudx=0.0_pr
	dudy=0.0_pr
	dvdx=0.0_pr
	dvdy=0.0_pr

	call derivadas(ut,vt)
	call visc
	call tnlinear(tnlx,tnly) !Calculo do termo nao-linear
	!t3=((nstep)+c(ii))*dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Zona de amortecimentao (Buffer Zone)
! 	ZBx=phi*(ufis-Qtx)
! 	ZBy=phi*(vfis-Qty)
! 	call ZFFT2D(ZBx,NX,NY,-1)
! 	call ZFFT2D(ZBy,NX,NY,-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	RHSx=-viscx-tnlx+ffx+empx   !calculo do lado direito da equacao em x
	RHSy=-viscy-tnly+ffy+empy   !calculo do lado direito da equacao em y

	call projecao(RHSx,RHSy)

	!coeficientes do RK4
	k1x=alfa(ii)*k1x+dt*RHSx
	k1y=alfa(ii)*k1y+dt*RHSy

	!atualizacao da predicao das velocidades
	u=u+k1x*beta(ii)
	v=v+k1y*beta(ii)

	ufis=u
	vfis=v
	call ZFFT2D(ufis,nx,ny,1)
	call ZFFT2D(vfis,nx,ny,1)
	ufis=real(ufis)
	vfis=real(vfis)
end do !fim do laco do Runge-Kutta de quarta ordem com 6 passos otimizado
end subroutine RK46DFBIF

end module avanco_temporal