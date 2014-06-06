! This file contains functions (to be wrapped with f2py) which are
! used with pyMission


subroutine getFinDiffCoef(M, N, x0, nCoord, coord, d)
  ! subroutine used to compute the finite element coefficients
  ! that are multiplied with discrete values to provide estimates
  ! of the derivatives
  ! these coefficients can be generated for uniform and non-uniform
  ! spacings

  implicit none

  !f2py intent(in) M, N, x0, nCoord, coord
  !f2py intent(out) d
  !f2py depend(M,N) d
  !f2py depend(nCoord) coord

  ! input/output
  integer, intent(in) :: M, N
  double precision, intent(in) :: x0
  integer, intent(in) :: nCoord
  double precision, dimension(0:nCoord), intent(in) :: coord
  double precision, dimension(0:M,0:N-1,0:N-1), intent(out) :: d

  ! variable initialization
  integer :: i = 0, j = 0, k = 0
  double precision :: c1 = 0.0, c2 = 0.0, c3 = 0.0

  do i = 0, M
     do j = 0, N-1
        do k = 0, N-1
           d(i,j,k) = 0.0
        end do
     end do
  end do

  d(0,0,0) = 1.0
  c1 = 1.0

  do j = 1, N-1
     c2 = 1.0

     do k = 0, j-1
        c3 = coord(j)-coord(k)
        c2 = c2*c3
        
        do i = 0, min(j,M)
           if (i == 0) then
              d(i,j,k) = ((coord(j)-x0)*d(i,j-1,k))/c3
           else
              d(i,j,k) = ((coord(j)-x0)*d(i,j-1,k)-i*d(i-1,j-1,k))/c3
           end if
        end do
     end do

     do i = 0, min(j,M)
        if (i == 0) then
           d(i,j,j) = -(c1/c2)*(coord(j-1)-x0)*d(i,j-1,j-1)
        else
           d(i,j,j) = (c1/c2)*(i*d(i-1,j-1,j-1)-(coord(j-1)-x0)*d(i,j-1,j-1))
        end if
     end do
     
     c1 = c2
  end do

end subroutine getFinDiffCoef

subroutine get_alpha(numElem, CL, eta, alpha)
  !f2py intent(in) numElem
  !f2py intent(in) CL, eta
  !f2py intent(out) alpha
  !f2py depend(numElem) CL, eta, alpha

  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(in) :: CL, eta
  double precision, dimension(0:numElem), intent(out) :: alpha

  alpha = (1/4.24)*(CL-0.26-0.27*eta)

end subroutine get_alpha

subroutine get_alpha_d(numElem, eta, CL, deta, dCL)
  ! compute the derivatives of alpha wrt alpha, eta, and CL

  !f2py intent(in) numElem
  !f2py intent(in) eta, CL
  !f2py intent(out) deta, dCL
  !f2py depend(numElem) eta, CL, deta, dCL

  ! Input/Output
  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(in) :: eta, CL
  double precision, dimension(0:numElem), intent(out) :: deta, dCL

  integer :: i = 0

  do i = 0,numElem
     deta(i) = -0.27/4.24
     dCL(i) = 1/4.24
  enddo

end subroutine get_alpha_d

subroutine get_CD(numElem, AR, e, CL, CD)
  ! compute CD using simple drag polar

  !f2py intent(in) numElem
  !f2py intent(in) AR, e
  !f2py intent(in) CL
  !f2py intent(out) CD
  !f2py depend(numElem) CL, CD

  ! Input/Output
  integer, intent(in) :: numElem
  double precision, intent(in) :: AR, e
  double precision, dimension(0:numElem), intent(in) :: CL
  double precision, dimension(0:numElem), intent(out) :: CD
  
  double precision, parameter :: PI = 3.14159265359

  CD = (0.018 + CL**2/(PI*AR*e))

end subroutine get_CD

subroutine get_CD_d(numElem, AR, e, CL, dAR, de, dCL)
  !f2py intent(in) numElem
  !f2py intent(in) AR, e
  !f2py intent(in) CL
  !f2py intent(out) dAR, de, dCL
  !f2py depend(numElem) CL, dAR, de, dCL

  integer, intent(in) :: numElem
  double precision, intent(in) :: AR, e
  double precision, dimension(0:numElem), intent(in) :: CL
  double precision, dimension(0:numElem), intent(out) :: dAR, de, dCL
  
  integer :: i = 0
  double precision, parameter :: PI = 3.14159254359

  do i = 0,numElem
     dAR(i) = -CL(i)**2/(PI*e*AR**2)
     de(i) = -CL(i)**2/(PI*AR*e**2)
     dCL(i) = 2*CL(i)/(PI*AR*e)
  enddo

end subroutine get_CD_d

subroutine get_eta(numElem, CM, alpha, eta)
  ! compute eta (tail rotation angle)

  !f2py intent(in) numSeg
  !f2py intent(in) CM, alpha
  !f2py intent(out) eta
  !f2py depend(numElem) CM, alpha, eta

  !Input/Output
  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(in) :: CM, alpha
  double precision, dimension(0:numElem), intent(out) :: eta

  eta = (1/1.06)*(0.63*alpha - CM)

end subroutine get_eta

subroutine get_eta_d(numElem, alpha, CM, dalpha, dCM)
  ! compute the derivatives of eta wrt alpha, eta, CM

  !f2py intent(in) numElem
  !f2py intent(in) alpha, CM
  !f2py intent(out) dalpha, dCM
  !f2py depend(numElem) alpha, CM, dalpha, dCM

  !Input/Output
  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(in) :: alpha, CM
  double precision, dimension(0:numElem), intent(out) :: dalpha, dCM

  integer :: i = 0

  do i = 0,numElem
     dalpha(i) = 0.63/1.06
     dCM(i) = -1/1.06
  enddo

end subroutine get_eta_d

subroutine get_sfc(numElem, SFCSL, h, SFC)
  !f2py intent(in) numElem
  !f2py intent(in) SFCSL
  !f2py intent(in) h
  !f2py intent(out) SFC
  !f2py depend(numElem) h, SFC

  integer, intent(in) :: numElem
  double precision, intent(in) :: SFCSL
  double precision, dimension(0:numElem), intent(in) :: h
  double precision, dimension(0:numElem), intent(out) :: SFC

  SFC = SFCSL + (6.39e-13)*h

end subroutine get_sfc

subroutine get_sfc_d(numElem, dSFC)
  !f2py intent(in) numElem
  !f2py intent(out) dSFC
  !f2py depend(numElem) dSFC

  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(out) :: dSFC

  integer :: i
  
  do i=0,numElem
     dSFC(i) = 6.39e-13
  enddo

end subroutine get_sfc_d

subroutine get_gamma(numElem, h, x, gamma)
  !f2py intent(in) numElem
  !f2py intent(in) h, x
  !f2py intent(out) gamma
  !f2py depend(numElem) h, x, gamma

  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(in) :: h, x
  double precision, dimension(0:numElem), intent(out) :: gamma

  integer :: i
  double precision :: dx

  dx = x(1)-x(0)
  do i=0,1
     gamma(i) = -(3.0/2.0)*h(i) + 2.0*h(i+1) - (1.0/2.0)*h(i+2)
     gamma(i) = gamma(i) / dx
  enddo

  do i=2,numElem-2
     gamma(i) = (1.0/12.0)*h(i-2) - (2.0/3.0)*h(i-1) + (2.0/3.0)*h(i+1) - (1.0/12.0)*h(i+2)
     gamma(i) = gamma(i) / dx
  enddo

  do i=numElem-2,numElem
     gamma(i) = (1.0/2.0)*h(i-2) - 2.0*h(i-1) + (3.0/2.0)*h(i)
     gamma(i) = gamma(i) / dx
  enddo
end subroutine get_gamma

subroutine get_gamma_d(numElem, h, x, dhgamma1, dhgamma2, dhgamma3, dhgamma4, &
     dhgamma5, dxgamma0, dxgamma1)
  !f2py intent(in) numElem
  !f2py intent(in) h, x
  !f2py intent(out) dhgamma1, dhgamma2, dhgamma3, dhgamma4, dhgamma5, dxgamma
  !f2py depend(numElem) h, x, dhgamma1, dhgamma2, dhgamma3, dhgamma4, dhgamma5, dxgamma0, dxgamma1

  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(in) :: h, x
  double precision, dimension(0:numElem), intent(out) :: dhgamma1, dhgamma2, &
       dhgamma3, dhgamma4, dhgamma5, dxgamma0, dxgamma1

  integer :: i
  double precision :: dx, dx0, dx1

  dx = x(1)-x(0)
  dx0 = -1/numElem
  dx1 = 1/numElem

  do i=0,1
     dhgamma1(i) = 0.0
     dhgamma2(i) = 0.0
     dhgamma3(i) = -3.0/2.0/dx
     dhgamma4(i) = 2.0/dx
     dhgamma5(i) = -1.0/2.0/dx
     dxgamma0(i) = (3.0/2.0)*h(i) - 2.0*h(i+1) + (1.0/2.0)*h(i+2)
     dxgamma0(i) = dxgamma0(i) / dx**2
     dxgamma1(i) = dxgamma0(i) * dx1
     dxgamma0(i) = dxgamma0(i) * dx0
  enddo

  do i=2,numElem-2
     dhgamma1(i) = 1.0/12.0/dx
     dhgamma2(i) = -2.0/3.0/dx
     dhgamma3(i) = 0.0
     dhgamma4(i) = 2.0/3.0/dx
     dhgamma5(i) = -1.0/12.0/dx
     dxgamma0(i) = -(1.0/12.0)*h(i-2) + (2.0/3.0)*h(i-1) - (2.0/3.0)*h(i+1) + (1.0/12.0)*h(i+2)
     dxgamma0(i) = dxgamma0(i) / dx**2
     dxgamma1(i) = dxgamma0(i) * dx1
     dxgamma0(i) = dxgamma0(i) * dx0
  enddo

  do i=numElem-2,numElem
     dhgamma1(i) = 1.0/2.0/dx
     dhgamma2(i) = -2.0/dx
     dhgamma3(i) = 3.0/2.0/dx
     dhgamma4(i) = 0.0
     dhgamma5(i) = 0.0
     dxgamma0(i) = -(1.0/2.0)*h(i-2) + 2.0*h(i-1) - (3.0/2.0)*h(i)
     dxgamma0(i) = dxgamma0(i) / dx**2
     dxgamma1(i) = dxgamma0(i) * dx1
     dxgamma0(i) = dxgamma0(i) * dx0
  enddo

end subroutine get_gamma_d

subroutine get_Temp(numElem, h, Temp)
  !f2py intent(in) numElem
  !f2py intent(in) h
  !f2py intent(out) Temp
  !f2py depend(numElem) h, Temp

  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(in) :: h
  double precision, dimension(0:numElem), intent(out) :: Temp

  integer :: i

  Temp = 288.16 - (6.5e-3)*h

  do i=0,numElem
     if (Temp(i) < 1) then
        Temp(i) = 1.0
     endif
  enddo

end subroutine get_Temp

subroutine get_Temp_d(numElem, h, dTemp)
  !f2py intent(in) numElem
  !f2py intent(in) h
  !f2py intent(out) dTemp
  !f2py depend(numElem) h, dTemp

  integer, intent(in) :: numElem
  double precision, dimension(0:numElem), intent(in) :: h
  double precision, dimension(0:numElem), intent(out) :: dTemp

  integer :: i
  double precision, dimension(0:numElem) :: Temp

  Temp = 288.16 - (6.5e-3)*h

  do i=0,numElem
     dTemp(i) = -6.5e-3
     if (Temp(i) < 1) then
        dTemp(i) = 0.0
     endif
  enddo
end subroutine get_Temp_d

subroutine get_rho(numElem, g, Temp, rho)
  !f2py intent(in) numElem
  !f2py intent(in) g
  !f2py intent(in) Temp
  !f2py intent(out) rho
  !f2py depend(numElem) Temp, rho

  integer, intent(in) :: numElem
  double precision, intent(in) :: g
  double precision, dimension(0:numElem), intent(in) :: Temp
  double precision, dimension(0:numElem), intent(out) :: rho

  integer :: i

  rho = 1.225*(Temp/288.16)**(-((9.81/((-6.5e-3)*287))+1))

  do i=0,numElem
     if (rho(i) < 0.01) then
        rho(i) = 0.01
     endif
  enddo

end subroutine get_rho

subroutine get_rho_d(numElem, g, Temp, dRho)
  !f2py intent(in) numElem
  !f2py intent(in) g
  !f2py intent(in) Temp
  !f2py intent(out) dRho
  !f2py depend(numElem) Temp, dRho

  integer, intent(in) :: numElem
  double precision, intent(in) :: g
  double precision, dimension(0:numElem), intent(in) :: Temp
  double precision, dimension(0:numElem), intent(out) :: dRho
  
  integer :: i
  double precision, dimension(0:numElem) :: rho

  rho = 1.225*(Temp/288.16)**(-((9.81/((-6.5e-3)*287))+1))

  do i=0,numElem
     dRho(i) = 1.225*(Temp(i)/288.16)**(-((9.81/((-6.5e-3)*287))+2))
     dRho(i) = dRho(i) * (-9.81/((-6.5e-3)*287)-1)*(1/288.16)
     if (rho(i) < 0.01) then
        dRho(i) = 0.0
     endif
  enddo
end subroutine get_rho_d

subroutine get_v(numElem, v_ends, v)
  !f2py intent(in) numElem
  !f2py intent(in) v_ends
  !f2py intent(out) v
  !f2py depend(numElem) v

  integer, intent(in) :: numElem
  double precision, dimension(0:1), intent(in) :: v_ends
  double precision, dimension(0:numElem), intent(out) :: v

  integer :: i
  double precision :: dv

  dv = (v_ends(1)-v_ends(0))/numElem
  do i=0,numElem
     v(i) = dv*i+v_ends(0)
  enddo
end subroutine get_v

subroutine get_v_d(numElem, v_ends, dv0, dv1)
  !f2py intent(in) numElem
  !f2py intent(in) v_ends
  !f2py intent(out) dv
  !f2py depend(numElem) dv

  integer, intent(in) :: numElem
  double precision, dimension(0:1), intent(in) :: v_ends
  double precision, dimension(0:numElem), intent(out) :: dv0, dv1

  integer :: i

  do i=0,numElem
     dv0(i) = -real(i)/numElem + 1
     dv1(i) = real(i)/numElem
  enddo
end subroutine get_v_d

subroutine get_CT_Init(numElem, tau, S, cThrustSL, rho, v, h, CT)
  !f2py intent(in) numElem
  !f2py intent(in) tau, S, cThrustSL
  !f2py intent(in) rho, v, h
  !f2py intent(out) CT
  !f2py depend(numElem) rho, v, h, CT

  integer, intent(in) :: numElem
  double precision, intent(in) :: tau, S, cThrustSL
  double precision, dimension(0:numElem), intent(in) :: rho, v, h
  double precision, dimension(0:numElem), intent(out) :: CT

  integer :: i
  double precision, dimension(0:numElem) :: cThrust

  do i=0,numElem
     cThrust(i) = cThrustSL - 0.072*h(i)
     CT(i) = (2/(rho(i)*v(i)**2*S))*cThrust(i)*tau
  enddo
end subroutine get_CT_Init

subroutine get_CT_Init_d(numElem, tau, S, cThrustSL, rho, v, h, &
     dcThrustSL, dh, drho, dv, dS, dtau)
  !f2py intent(in) numElem
  !f2py intent(in) tau, S, cThrustSL
  !f2py intent(in) rho, v, h
  !f2py intent(out) dcThrustSL, dh, drho, dv, dS, dtau
  !f2py depend(numElem) rho, v, h, dcThrustSL, dh, drho, dv, dS, dtau

  integer, intent(in) :: numElem
  double precision, intent(in) :: tau, S, cThrustSL
  double precision, dimension(0:numElem), intent(in) :: rho, v, h
  double precision, dimension(0:numElem), intent(out) :: dcThrustSL, &
       dh, drho, dv, dS, dtau

  integer :: i
  double precision :: cThrust

  do i=0,numElem
     cThrust = cThrustSL - 0.072*h(i)
     dcThrustSL(i) = 2*tau/(rho(i)*v(i)**2*S)
     dh(i) = -0.144*tau/(rho(i)*v(i)**2*S)
     drho(i) = -2*cThrust*tau/(rho(i)**2*v(i)**2*S)
     dv(i) = -4*cThrust*tau/(rho(i)*v(i)**3*S)
     dS(i) = -2*cThrust*tau/(rho(i)*v(i)**2*S**2)
     dtau(i) = 2*cThrust/(rho(i)*v(i)**2*S)
  enddo
end subroutine get_CT_Init_d

subroutine get_CT(numElem, numInt, S, Wac, x, alpha, rho, v, CD, Wf, gamma, CT, CTRes)
  !f2py intent(in) numElem, numInt
  !f2py intent(in) S, Wac
  !f2py intent(in) x, alpha, rho, v, CD, Wf, gamma, CT
  !f2py intent(out) CTRes
  !f2py depend(numElem) x, alpha, rho, v, CD, Wf, gamma, CT, CTRes

  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: S, Wac
  double precision, dimension(0:numElem), intent(in) :: x, alpha, rho, v, CD, Wf, &
       gamma, CT
  double precision, dimension(0:numElem), intent(out) :: CTRes

  integer :: i, j
  double precision :: dx, temp
  double precision, parameter :: param_zero = 0.0, param_one = 1.0
  double precision, dimension(0:numInt-1) :: R1, R2
  double precision, dimension(0:numInt-1) :: xTemp, vTemp, rhoTemp, QTemp
  double precision, dimension(0:numInt-1) :: gammaTemp, sinGamma, cosAlpha
  double precision, dimension(0:numInt-1) :: CTTemp, WTemp, aTemp, CDTemp

  do i = 0, numElem
     CTRes(i) = 0.0
  enddo

  call linspace(numInt, param_one, param_zero, R1)
  call linspace(numInt, param_zero, param_one, R2)

  do i = 0, numElem-1
     call linspace(numInt, rho(i), rho(i+1), rhoTemp)
     call linspace(numInt, v(i), v(i+1), vTemp)
     dx = (x(i+1)-x(i))/numInt
     
     do j = 0, numInt-1
        QTemp(j) = 0.5*rhoTemp(j)*vTemp(j)**2*S
        xTemp(j) = dx
     enddo

     xTemp(0) = 0.5*dx
     xTemp(numInt-1) = 0.5*dx
     
     call linspace(numInt, gamma(i), gamma(i+1), gammaTemp)
     sinGamma = sin(gammaTemp)
     
     call linspace(numInt, CD(i), CD(i+1), CDTemp)
     call linspace(numInt, Wf(i), Wf(i+1), WTemp)
     call linspace(numInt, alpha(i), alpha(i+1), aTemp)
     call linspace(numInt, CT(i), CT(i+1), CTTemp)

     cosAlpha = cos(aTemp)

     do j = 0, numInt-1
        WTemp(j) = WTemp(j) + Wac
     enddo

     do j = 0, numInt-1
        temp = -QTemp(j)*CTTemp(j)*cosAlpha(j) + &
             QTemp(j)*CDTemp(j) + WTemp(j)*sinGamma(j)
        temp = temp*R1(j)*xTemp(j)
        CTRes(i) = CTRes(i) + temp
        temp = -QTemp(j)*CTTemp(j)*cosAlpha(j) + &
             QTemp(j)*CDTemp(j) + WTemp(j)*sinGamma(j)
        temp = temp*R2(j)*xTemp(j)
        CTRes(i+1) = CTRes(i+1) + temp
     enddo
  enddo
end subroutine get_CT

subroutine get_CT_d(numElem, numInt, S, Wac, x, alpha, rho, v, CD, Wf, gamma, CT, &
     dCTdS, dCTdWac, dCTdAlpha1, dCTdAlpha2, dCTdAlpha3, dCTdRho1, dCTdRho2, &
     dCTdRho3, dCTdV1, dCTdV2, dCTdV3, dCTdCD1, dCTdCD2, dCTdCD3, dCTdWf1, &
     dCTdWf2, dCTdWf3, dCTdGamma1, dCTdGamma2, dCTdGamma3, dCTdCT1, dCTdCT2, &
     dCTdCT3)
  !f2py intent(in) numElem, numInt
  !f2py intent(in) S, Wac
  !f2py intent(in) x, alpha, rho, v, CD, Wf, gamma, CT
  !f2py intent(out) dCTdS, dCTdWac, dCTdAlpha1, dCTdAlpha2, dCTdAlpha3, dCTdRho1, dCTdRho2, dCTdRho3, dCTdV1, dCTdV2, dCTdV3, dCTdCD1, dCTdCD2, dCTdCD3, dCTdWf1, dCTdWf2, dCTdWf3, dCTdGamma1, dCTdGamma2, dCTdGamma3, dCTdCT1, dCTdCT2, dCTdCT3
  !f2py depend(numElem) x, alpha, rho, v, CD, Wf, gamma, CT, dCTdS, dCTdWac, dCTdAlpha1, dCTdAlpha2, dCTdAlpha3, dCTdRho1, dCTdRho2, dCTdRho3, dCTdV1, dCTdV2, dCTdV3, dCTdCD1, dCTdCD2, dCTdCD3, dCTdWf1, dCTdWf2, dCTdWf3, dCTdGamma1, dCTdGamma2, dCTdGamma3, dCTdCT1, dCTdCT2, dCTdCT3

  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: S, Wac
  double precision, dimension(0:numElem), intent(in) :: x, alpha, rho, v, CD, Wf, &
       gamma, CT
  double precision, dimension(0:numElem), intent(out) :: dCTdS, dCTdWac, dCTdAlpha1, & 
       dCTdAlpha2, dCTdAlpha3, dCTdRho1, dCTdRho2, dCTdRho3, dCTdV1, dCTdV2, dCTdV3, &
       dCTdCD1, dCTdCD2, dCTdCD3, dCTdWf1, dCTdWf2, dCTdWf3, dCTdGamma1, dCTdGamma2, &
       dCTdGamma3, dCTdCT1, dCTdCT2, dCTdCT3

  integer :: i, j
  double precision :: temp, deltax
  double precision, parameter :: param_zero = 0.0, param_one = 1.0
  double precision, dimension(0:numInt-1) :: R1, R2
  double precision, dimension(0:numInt-1) :: xTemp, aTemp, rhoTemp, vTemp, CDTemp, &
       WTemp, gammaTemp, CTTemp, sinGamma, cosAlpha, QTemp, dQTempdRho1, &
       dQTempdRho2, dQTempdV1, dQTempdV2, dQTempdS, dRhoTemp1, dRhoTemp2, dVTemp1, &
       dVTemp2, dGammaTemp1, dGammaTemp2, dSinGamma1, dSinGamma2, dCDTemp1, &
       dCDTemp2, dWTempdWf1, dWTempdWf2, dWTempdWac, dATemp1, dATemp2, dCTTemp1, &
       dCTTemp2, dCosAlpha1, dCosAlpha2

  call linspace(numInt, param_one, param_zero, R1)
  call linspace(numInt, param_zero, param_one, R2)

  do i = 0, numElem
     dCTdWac(i) = 0.0
     dCTdS(i) = 0.0
     dCTdV1(i) = 0.0
     dCTdV2(i) = 0.0
     dCTdV3(i) = 0.0
     dCTdRho1(i) = 0.0
     dCTdRho2(i) = 0.0
     dCTdRho3(i) = 0.0
     dCTdCD1(i) = 0.0
     dCTdCD2(i) = 0.0
     dCTdCD3(i) = 0.0
     dCTdWf1(i) = 0.0
     dCTdWf2(i) = 0.0
     dCTdWf3(i) = 0.0
     dCTdGamma1(i) = 0.0
     dCTdGamma2(i) = 0.0
     dCTdGamma3(i) = 0.0
     dCTdCT1(i) = 0.0
     dCTdCT2(i) = 0.0
     dCTdCT3(i) = 0.0
     dCTdAlpha1(i) = 0.0
     dCTdAlpha2(i) = 0.0
     dCTdAlpha3(i) = 0.0
  enddo

  do i = 0, numElem-1
     call linspace(numInt, rho(i), rho(i+1), rhoTemp)
     call linspace(numInt, v(i), v(i+1), vTemp)
     call dlinspace(numInt, rho(i), rho(i+1), dRhoTemp1, dRhoTemp2)
     call dlinspace(numInt, v(i), v(i+1), dVTemp1, dVTemp2)
     deltax = (x(i+1)-x(i))/numInt
     
     do j = 0, numInt-1
        QTemp(j) = 0.5*rhoTemp(j)*vTemp(j)**2*S
        dQTempdRho1(j) = 0.5*vTemp(j)**2*S*dRhoTemp1(j)
        dQTempdRho2(j) = 0.5*vTemp(j)**2*S*dRhoTemp2(j)
        dQTempdV1(j) = rhoTemp(j)*vTemp(j)*S*dVTemp1(j)
        dQTempdV2(j) = rhoTemp(j)*vTemp(j)*S*dVTemp2(j)
        dQTempdS(j) = 0.5*rhoTemp(j)*vTemp(j)**2
        xTemp(j) = deltax
     enddo

     xTemp(0) = 0.5*deltax
     xTemp(numInt-1) = 0.5*deltax
     
     call linspace(numInt, gamma(i), gamma(i+1), gammaTemp)
     call dlinspace(numInt, gamma(i), gamma(i+1), dGammaTemp1, dGammaTemp2)
     sinGamma = sin(gammaTemp)
     dSinGamma1 = cos(gammaTemp)*dGammaTemp1
     dSinGamma2 = cos(gammaTemp)*dGammaTemp2
     
     call linspace(numInt, CD(i), CD(i+1), CDTemp)
     call dlinspace(numInt, CD(i), CD(i+1), dCDTemp1, dCDTemp2)
     call linspace(numInt, alpha(i), alpha(i+1), aTemp)
     call dlinspace(numInt, alpha(i), alpha(i+1), dATemp1, dATemp2)
     call linspace(numInt, CT(i), CT(i+1), CTTemp)
     call dlinspace(numInt, CT(i), CT(i+1), dCTTemp1, dCTTemp2)
     call linspace(numInt, Wf(i), Wf(i+1), WTemp)
     call dlinspace(numInt, Wf(i), Wf(i+1), dWTempdWf1, dWTempdWf2)

     cosAlpha = cos(aTemp)
     dCosAlpha1 = -sin(aTemp)*dATemp1
     dCosAlpha2 = -sin(aTemp)*dATemp2

     do j = 0, numInt-1
        WTemp(j) = WTemp(j) + Wac
        dWTempdWac(j) = 1.0
     enddo

     do j = 0, numInt-1
        dCTdWac(i) = dCTdWac(i) + R1(j)*xTemp(j)*(sinGamma(j))
        dCTdS(i) = dCTdS(i) + R1(j)*xTemp(j)*(-dQTempdS(j)*CTTemp(j)* &
             cosAlpha(j)+dQTempdS(j)*CDTemp(j))
        dCTdV2(i) = dCTdV2(i) + R1(j)*xTemp(j)*(-dQTempdV1(j)*CTTemp(j)* &
             cosAlpha(j)+dQTempdV1(j)*CDTemp(j))
        dCTdV3(i) = dCTdV3(i) + R1(j)*xTemp(j)*(-dQTempdV2(j)*CTTemp(j)* &
             cosAlpha(j)+dQTempdV2(j)*CDTemp(j))
        dCTdRho2(i) = dCTdRho2(i) + R1(j)*xTemp(j)*(-dQTempdRho1(j)* &
             CTTemp(j)*cosAlpha(j)+dQTempdRho1(j)*CDTemp(j))
        dCTdRho3(i) = dCTdRho3(i) + R1(j)*xTemp(j)*(-dQTempdRho2(j)* &
             CTTemp(j)*cosAlpha(j)+dQTempdRho2(j)*CDTemp(j))
        dCTdCD2(i) = dCTdCD2(i) + R1(j)*xTemp(j)*(QTemp(j)*dCDTemp1(j))
        dCTdCD3(i) = dCTdCD3(i) + R1(j)*xTemp(j)*(QTemp(j)*dCDTemp2(j))
        dCTdWf2(i) = dCTdWf2(i) + R1(j)*xTemp(j)*(dWTempdWf1(j)* &
             sinGamma(j))
        dCTdWf3(i) = dCTdWf3(i) + R1(j)*xTemp(j)*(dWTempdWf2(j)* &
             sinGamma(j))
        dCTdGamma2(i) = dCTdGamma2(i) + R1(j)*xTemp(j)*(WTemp(j)* &
             dSinGamma1(j))
        dCTdGamma3(i) = dCTdGamma3(i) + R1(j)*xTemp(j)*(WTemp(j)* &
             dSinGamma2(j))
        dCTdCT2(i) = dCTdCT2(i) + R1(j)*xTemp(j)*(-QTemp(j)*cosAlpha(j)* &
             dCTTemp1(j))
        dCTdCT3(i) = dCTdCT3(i) + R1(j)*xTemp(j)*(-QTemp(j)*cosAlpha(j)* &
             dCTTemp2(j))
        dCTdAlpha2(i) = dCTdAlpha2(i) + R1(j)*xTemp(j)*(-QTemp(j)* &
             CTTemp(j)*dCosAlpha1(j))
        dCTdAlpha3(i) = dCTdAlpha3(i) + R1(j)*xTemp(j)*(-QTemp(j)* &
             CTTemp(j)*dCosAlpha2(j))

        dCTdWac(i+1) = dCTdWac(i+1) + R2(j)*xTemp(j)*(sinGamma(j))
        dCTdS(i+1) = dCTdS(i+1) + R2(j)*xTemp(j)*(-dQTempdS(j)*CTTemp(j)* &
             cosAlpha(j)+dQTempdS(j)*CDTemp(j))
        dCTdV1(i+1) = dCTdV1(i+1) + R2(j)*xTemp(j)*(-dQTempdV1(j)* &
             CTTemp(j)*cosAlpha(j)+dQTempdV1(j)*CDTemp(j))
        dCTdV2(i+1) = dCTdV2(i+1) + R2(j)*xTemp(j)*(-dQTempdV2(j)* &
             CTTemp(j)*cosAlpha(j)+dQTempdV2(j)*CDTemp(j))
        dCTdRho1(i+1) = dCTdRho1(i+1) + R2(j)*xTemp(j)*(-dQTempdRho1(j)* &
             CTTemp(j)*cosAlpha(j)+dQTempdRho1(j)*CDTemp(j))
        dCTdRho2(i+1) = dCTdRho2(i+1) + R2(j)*xTemp(j)*(-dQTempdRho2(j)* &
             CTTemp(j)*cosAlpha(j)+dQTempdRho2(j)*CDTemp(j))
        dCTdCD1(i+1) = dCTdCD1(i+1) + R2(j)*xTemp(j)*(QTemp(j)*dCDTemp1(j))
        dCTdCD2(i+1) = dCTdCD2(i+1) + R2(j)*xTemp(j)*(QTemp(j)*dCDTemp2(j))
        dCTdWf1(i+1) = dCTdWf1(i+1) + R2(j)*xTemp(j)*(dWTempdWf1(j)* &
             sinGamma(j))
        dCTdWf2(i+1) = dCTdWf2(i+1) + R2(j)*xTemp(j)*(dWTempdWf2(j)* &
             sinGamma(j))
        dCTdGamma1(i+1) = dCTdGamma1(i+1) + R2(j)*xTemp(j)*(WTemp(j)* &
             dSinGamma1(j))
        dCTdGamma2(i+1) = dCTdGamma2(i+1) + R2(j)*xTemp(j)*(WTemp(j)* &
             dSinGamma2(j))
        dCTdCT1(i+1) = dCTdCT1(i+1) + R2(j)*xTemp(j)*(-QTemp(j)* &
             cosAlpha(j)*dCTTemp1(j))
        dCTdCT2(i+1) = dCTdCT2(i+1) + R2(j)*xTemp(j)*(-QTemp(j)* &
             cosAlpha(j)* dCTTemp2(j))
        dCTdAlpha1(i+1) = dCTdAlpha1(i+1) + R2(j)*xTemp(j)*(-QTemp(j)* &
             CTTemp(j)*dCosAlpha1(j))
        dCTdAlpha2(i+1) = dCTdAlpha2(i+1) + R2(j)*xTemp(j)*(-QTemp(j)* &
             CTTemp(j)*dCosAlpha2(j))
     enddo
  enddo
end subroutine get_CT_d

subroutine get_h(numElem, numInt, S, Wac, x_ends, h_ends, Wf, CT, alpha, CD, rho, v, h)
  !f2py intent(in) numElem, numInt
  !f2py intent(in) S, Wac
  !f2py intent(in) x_ends, h_ends
  !f2py intent(in) Wf, CT, alpha, CD, rho, v
  !f2py intent(out) h
  !f2py depend(numElem) Wf, CT, alpha, CD, rho, v, h

  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: S, Wac
  double precision, dimension(0:1), intent(in) :: x_ends, h_ends
  double precision, dimension(0:numElem), intent(in) :: Wf, &
       CT, alpha, CD, rho, v
  double precision, dimension(0:numElem), intent(out) :: h

  integer :: i
  double precision :: dx, dh
  double precision, dimension(0:numElem) :: x
  double precision, dimension(0:numInt) :: xTemp, WfTemp, CTTemp, aTemp, &
       CDTemp, rhoTemp, vTemp

  dx = (x_ends(1)-x_ends(0))/numElem

  do i=0,numElem
     x(i) = i*dx+x_ends(0)
     h(i) = 0
  enddo

  h(0) = h_ends(0)
  do i=0,numElem-1
     h(i+1) = h(i)
     call linspace(numInt+1,x(i),x(i+1),xTemp)
     call linspace(numInt+1,Wf(i),Wf(i+1),WfTemp)
     call linspace(numInt+1,CT(i),CT(i+1),CTTemp)
     call linspace(numInt+1,alpha(i),alpha(i+1),aTemp)
     call linspace(numInt+1,CD(i),CD(i+1),CDTemp)
     call linspace(numInt+1,rho(i),rho(i+1),rhoTemp)
     call linspace(numInt+1,v(i),v(i+1),vTemp)

     do j=0,numInt-1
        dh = ((xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac))*(0.5*rhoTemp(j)*vTemp(j)**2*S)
        dh = dh*(CTTemp(j)*COS(aTemp(j))-CDTemp(j))
        h(i+1) = h(i+1)+dh
     enddo
  enddo

end subroutine get_h

subroutine get_h_d(numElem, numInt, S, Wac, x_ends, h_ends, Wf, CT, alpha, CD, rho, &
     v, dhdS, dhdWac, dhdx0, dhdx1, dhdWf1, dhdWf2, dhdRho1, dhdRho2, dhdV1, dhdV2, dhdCT1, &
     dhdCT2, dhda1, dhda2, dhdCD1, dhdCD2)
  !f2py intent(in) numElem, numInt
  !f2py intent(in) S, Wac
  !f2py intent(in) x_ends, h_ends
  !f2py intent(in) Wf, CT, alpha, CD, rho, v
  !f2py intent(out) dhdS, dhdWac, dhdx0, dhdx1, dhdWf1, dhdWf2, dhdRho1, dhdRho2, dhdV1, dhdV2, dhdCT1, dhdCT2, dhda1, dhda2, dhdCD1, dhdCD2
  !f2py depend(numElem) Wf, CT, alpha, CD, rho, v, dhdS, dhdWac, dhdx0, dhdx1, dhdWf1, dhdWf2, dhdRho1, dhdRho2, dhdV1, dhdV2, dhdCT1, dhdCT2, dhda1, dhda1, dhdCD1, dhdCD2

  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: S, Wac
  double precision, dimension(0:1), intent(in) :: x_ends, h_ends
  double precision, dimension(0:numElem), intent(in) :: Wf, &
       CT, alpha, CD, rho, v
  double precision, dimension(0:numElem), intent(out) :: dhdS, dhdWac, dhdx0, dhdx1, &
       dhdWf1, dhdWf2, dhdRho1, dhdRho2, dhdV1, dhdV2, dhdCT1, dhdCT2, &
       dhda1, dhda2, dhdCD1, dhdCD2

  double precision :: cosa, delta_x, temp
  double precision, dimension(0:numElem) :: delxTemp0, delxTemp1, x
  double precision, dimension(0:numInt) :: dxTemp0, dxTemp1, dWfTemp1, &
       dWfTemp2, dCTTemp1, dCTTemp2, daTemp1, daTemp2, dCDTemp1, dCDTemp2, &
       drhoTemp1, drhoTemp2, dvTemp1, dvTemp2, xTemp, WfTemp, CTTemp, aTemp, &
       CDTemp, rhoTemp, vTemp

  delta_x = (x_ends(1)-x_ends(0))/numElem

  do i=0,numElem
     x(i) = i*delta_x+x_ends(0)
     delxTemp0(i) = -i/numElem+1
     delxTemp1(i) = i/numElem
     dhdS(i) = 0.0
     dhdWac(i) = 0.0
     dhdx0(i) = 0.0
     dhdx1(i) = 0.0
     dhdWf1(i) = 0.0
     dhdWf2(i) = 0.0
     dhdRho1(i) = 0.0
     dhdRho2(i) = 0.0
     dhdV1(i) = 0.0
     dhdV2(i) = 0.0
     dhdCT1(i) = 0.0
     dhdCT2(i) = 0.0
     dhda1(i) = 0.0
     dhda2(i) = 0.0
     dhdCD1(i) = 0.0
     dhdCD2(i) = 0.0
  enddo

  do i=0,numElem-1

     call dlinspace(numInt+1,x(i),x(i+1),dxTemp0,dxTemp1)
     call dlinspace(numInt+1,Wf(i),Wf(i+1),dWfTemp1,dWfTemp2)
     call dlinspace(numInt+1,CT(i),CT(i+1),dCTTemp1,dCTTemp2)
     call dlinspace(numInt+1,alpha(i),alpha(i+1),daTemp1,daTemp2)
     call dlinspace(numInt+1,CD(i),CD(i+1),dCDTemp1,dCDTemp2)
     call dlinspace(numInt+1,rho(i),rho(i+1),drhoTemp1,drhoTemp2)
     call dlinspace(numInt+1,v(i),v(i+1),dvTemp1,dvTemp2)

     call linspace(numInt+1,x(i),x(i+1),xTemp)
     call linspace(numInt+1,Wf(i),Wf(i+1),WfTemp)
     call linspace(numInt+1,CT(i),CT(i+1),CTTemp)
     call linspace(numInt+1,alpha(i),alpha(i+1),aTemp)
     call linspace(numInt+1,CD(i),CD(i+1),CDTemp)
     call linspace(numInt+1,rho(i),rho(i+1),rhoTemp)
     call linspace(numInt+1,v(i),v(i+1),vTemp)

     do j=0,numInt-1
        cosa = COS(aTemp(j))
        
        temp = ((xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac))
        temp = temp * (0.5*rhoTemp(j)*vTemp(j)**2)
        temp = temp * (CTTemp(j)*cosa-CDTemp(j))
        dhdS(i+1) = dhdS(i+1) + temp

        temp = dxTemp1(j+1)*delxTemp0(i+1)+dxTemp0(j+1)*delxTemp0(i)
        temp = temp - dxTemp1(j)*delxTemp0(i+1)-dxTemp0(j)*delxTemp0(i)
        temp = temp / (WfTemp(j)+Wac) * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdx0(i+1) = dhdx0(i+1) + temp

        temp = dxTemp1(j+1)*delxTemp1(i+1)+dxTemp0(j+1)*delxTemp1(i)
        temp = temp - dxTemp1(j)*delxTemp1(i+1)-dxTemp0(j)*delxTemp1(i)
        temp = temp / (WfTemp(j)+Wac) * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdx1(i+1) = dhdx1(i+1) + temp

        temp = -(xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)**2
        temp = temp * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdWac(i+1) = dhdWac(i+1) + temp

        temp = -(xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)**2
        temp = temp * dWfTemp1(j)*(0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdWf1(i+1) = dhdWf1(i+1) + temp
        
        temp = -(xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)**2
        temp = temp * dWfTemp2(j)*(0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdWf2(i+1) = dhdWf2(i+1) + temp
        
        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (0.5*vTemp(j)**2*S)*drhoTemp1(j)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdRho1(i+1) = dhdRho1(i+1) + temp
        
        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (0.5*vTemp(j)**2*S)*drhoTemp2(j)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdRho2(i+1) = dhdRho2(i+1) + temp

        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (rhoTemp(j)*vTemp(j)*dvTemp1(j)*S)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdV1(i+1) = dhdV1(i+1) + temp
        
        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (rhoTemp(j)*vTemp(j)*dvTemp2(j)*S)
        temp = temp * (CTTemp(j)*cosa - CDTemp(j))
        dhdV2(i+1) = dhdV2(i+1) + temp

        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (dCTTemp1(j)*cosa)
        dhdCT1(i+1) = dhdCT1(i+1) + temp
        
        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (dCTTemp2(j)*cosa)
        dhdCT2(i+1) = dhdCT2(i+1) + temp

        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (-CTTemp(j)*SIN(aTemp(j))*daTemp1(j))
        dhda1(i+1) = dhda1(i+1) + temp

        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (-CTTemp(j)*SIN(aTemp(j))*daTemp2(j))
        dhda2(i+1) = dhda2(i+1) + temp

        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (-dCDTemp1(j))
        dhdCD1(i+1) = dhdCD1(i+1) + temp
        
        temp = (xTemp(j+1)-xTemp(j))/(WfTemp(j)+Wac)
        temp = temp * (0.5*rhoTemp(j)*vTemp(j)**2*S)
        temp = temp * (-dCDTemp2(j))
        dhdCD2(i+1) = dhdCD2(i+1) + temp
     enddo
  enddo
end subroutine get_h_d

subroutine linspace(n, x0, x1, y)
  ! fortran implementation of linspace function

  !Input/Output
  integer, intent(in) :: n
  double precision, intent(in) :: x0, x1
  double precision, dimension(0:n-1), intent(out) :: y

  double precision :: dx
  integer :: i

  dx = (x1-x0)/(n-1)

  do i = 0,n-1
     y(i) = x0 + dx*i
  enddo

end subroutine linspace

subroutine dlinspace(n, x0, x1, dy1, dy2)
  ! computes the derivatives of the output of linspace function
  ! wrt to the two inputs x0, and x1

  !Input/Output
  integer, intent(in) :: n
  double precision, intent(in) :: x0, x1
  double precision, dimension(0:n-1), intent(out) :: dy1, dy2

  integer :: i

  do i = 0,n-1
     dy1(i) = (1-real(i)/real(n-1))
     dy2(i) = real(i)/real(n-1)
  enddo

end subroutine dlinspace

subroutine get_tau(numElem, cThrustSL, S, CT, rho, v, h, tau)
  !f2py intent(in) numElem
  !f2py intent(in) cThrustSL, S
  !f2py intent(in) CT, rho, v, h
  !f2py intent(out) tau
  !f2py depend(numElem) CT, rho, v, h, tau

  integer, intent(in) :: numElem
  double precision, intent(in) :: cThrustSL, S
  double precision, dimension(0:numElem), intent(in) :: CT, rho, v, h
  double precision, dimension(0:numElem), intent(out) :: tau

  integer :: i = 0
  double precision :: a
  double precision, dimension(0:numElem) :: cThrust, Thrust

  a = 2.0

  do i = 0,numElem
     cThrust(i) = cThrustSL - 0.072*h(i)
     Thrust(i) = 0.5*rho(i)*v(i)**2*S*CT(i)
     !tau(i) = Thrust(i)/cThrust(i)
     tau(i) = 1 - (1/a)*LOG((1-EXP(a))*Thrust(i)/cThrust(i)+EXP(a))
  enddo

end subroutine get_tau

subroutine get_tau_d(numElem, cThrustSL, S, h, CT, rho, v, dCThrustSL, dh, dCT, &
     dRho, dV, dS)
  ! computes the derivatives of throttle setting wrt
  ! propulsion parameters, h, tau, and thrust

  !f2py intent(in) numElem
  !f2py intent(in) cThrustSL, S
  !f2py intent(in) h, CT, rho, v
  !f2py intent(out) dCThrustSL, dh, dCT, dRho, dV, dS
  !f2py depend(numElem) h, CT, rho, v, dCThrustSL, dh, dCT, dRho, dV, dS

  integer, intent(in) :: numElem
  double precision, intent(in) :: cThrustSL, S
  double precision, dimension(0:numElem), intent(in) :: h, CT, rho, v
  double precision, dimension(0:numElem), intent(out) :: dCThrustSL, dh, dCT, dRho, dV, dS

  integer :: i = 0
  double precision :: a
  double precision, dimension(0:numElem) :: cThrust

  a = 2.0

  do i = 0,numElem
     !dRho(i) = (0.5*v(i)**2*S*CT(i))/(cThrustSL-0.072*h(i))
     !dV(i) = (rho(i)*v(i)*S*CT(i))/(cThrustSL-0.072*h(i))
     !dS(i) = (0.5*rho(i)*v(i)**2*CT(i))/(cThrustSL-0.072*h(i))
     !dCT(i) = (0.5*rho(i)*v(i)**2*S)/(cThrustSL-0.072*h(i))
     !dCThrustSL(i) = -(0.5*rho(i)*v(i)**2*S*CT(i))/(cThrustSL-0.072*h(i))**2
     !dh(i) = 0.072*(0.5*rho(i)*v(i)**2*S*CT(i))/(cThrustSL-0.072*h(i))**2
     cThrust(i) = cThrustSL - 0.072*h(i)
     dRho(i) = -(1/a)*(1/(0.5*(1-EXP(a))*rho(i)*v(i)**2*S*CT(i)/cThrust(i)+EXP(a)))* &
          (0.5*(1-EXP(a))*v(i)**2*S*CT(i)/cThrust(i))
     dV(i) = -(1/a)*(1/(0.5*(1-EXP(a))*rho(i)*v(i)**2*S*CT(i)/cThrust(i)+EXP(a)))* &
          ((1-EXP(a))*rho(i)*v(i)*S*CT(i)/cThrust(i))
     dS(i) = -(1/a)*(1/(0.5*(1-EXP(a))*rho(i)*v(i)**2*S*CT(i)/cThrust(i)+EXP(a)))* &
          (0.5*(1-EXP(a))*rho(i)*v(i)**2*CT(i)/cThrust(i))
     dCT(i) = -(1/a)*(1/(0.5*(1-EXP(a))*rho(i)*v(i)**2*S*CT(i)/cThrust(i)+EXP(a)))* &
          (0.5*(1-EXP(a))*rho(i)*v(i)**2*S/cThrust(i))
     dCThrustSL(i) = -(1/a)*(1/(0.5*(1-EXP(a))*rho(i)*v(i)**2*S*CT(i)/cThrust(i)+EXP(a)))* &
          (-0.5*(1-EXP(a))*rho(i)*v(i)**2*S*CT(i)/cThrust(i)**2)
     dh(i) = -(1/a)*(1/(0.5*(1-EXP(a))*rho(i)*v(i)**2*S*CT(i)/cThrust(i)+EXP(a)))* &
          (-0.5*(1-EXP(a))*rho(i)*v(i)**2*S*CT(i)*(-0.072)/cThrust(i)**2)
  enddo

end subroutine get_tau_d

subroutine get_Wf(numInt, numElem, x, v, gamma, CT, SFC, rho, WfIn, g, WfSeg, S, Wf)
  ! computes fuel weight for each segment
  ! using a forward Euler scheme implemented backwards from the
  ! end of the mission
  ! the current initial condition is set to be 0 fuel at the end
  ! of the mission

  !f2py intent(in) numInt, numElem
  !f2py intent(in) x, v, gamma, CT, SFC, rho, WfIn
  !f2py intent(in) g, WfSeg, S
  !f2py intent(out) Wf
  !f2py depend(numElem) x, v, gamma, CT, SFC, rho, WfIn, Wf

  !Input/Output
  integer, intent(in) :: numInt, numElem
  double precision, dimension(0:numElem), intent(in) :: x, v, gamma
  double precision, dimension(0:numElem), intent(in) :: CT, SFC, rho, WfIn
  double precision, intent(in) :: g, WfSeg, S
!  double precision, intent(in) :: Wff
  double precision, dimension(0:numElem), intent(out) :: Wf

  double precision, dimension(0:numInt-1) :: vTemp, xTemp, SFCTemp, gammaTemp, cosGamma, rhoTemp, QTemp, CTTemp
  double precision :: dx, WfTemp
  integer :: i = 0, j = 0, k = 0, l = 0

  do i = 0,numElem
     Wf(i) = 0.0
  enddo

  do i = 0,numElem-1

     j = numElem-1-i
     call linspace(numInt, v(j), v(j+1), vTemp)
     dx = (x(j+1)-x(j))/real(numInt)

     do k = 0,numInt-1
        xTemp(k) = dx
     enddo

     xTemp(0) = 0.5*dx
     xTemp(numInt-1) = 0.5*dx
     call linspace(numInt, SFC(j), SFC(j+1), SFCTemp)
     call linspace(numInt, gamma(j), gamma(j+1), gammaTemp)
     call linspace(numInt, CT(j), CT(j+1), CTTemp)
     call linspace(numInt, rho(j), rho(j+1), rhoTemp)
     QTemp = 0.5*rhoTemp*vTemp**2*S
     cosGamma = cos(gammaTemp)
     
     do k = 0,numInt-1
        WfTemp = SFCTemp(k)*CTTemp(k)*QTemp(k)
        WfTemp = WfTemp*xTemp(k)/(vTemp(k)*cosGamma(k))
        !WfTemp = vTemp(k)*1e4
        Wf(j) = Wf(j) + WfTemp*9.81
     enddo
  enddo

  do i = 0,numElem-1
     j = numElem-1-i
     Wf(j) = Wf(j) + Wf(j+1)
  enddo

  do i = 0,numElem
     Wf(i) = Wf(i) + WfSeg
  enddo

end subroutine get_Wf

subroutine get_Wf_d(numElem, numInt, g, WfSeg, S, x, v, gamma, CT, SFC, rho, &
     dCT1, dCT2, dSFC1, dSFC2, dV1, dV2, dGamma1, dGamma2, dRho1, dRho2, dS)
  ! computes the derivatives of fuel weights wrt
  ! propulsion parameters, thrust, v, gamma, and fuel weight
  ! the resultant jacobian is bi-diagonal, and is stored in 2
  ! vectors for each derivative, ie, d___1 and d___2

  !f2py intent(in) numElem, numInt
  !f2py intent(in) g, WfSeg, S
  !f2py intent(in) x, v, gamma, CT, SFC, rho
  !f2py intent(out) dCT1, dCT2, dSFC1, dSFC2, dV1, dV2, dGamma1, dGamma2, dRho1, dRho2, dS
  !f2py depend(numElem) x, v, gamma, CT, SFC, rho, dCT1, dCT2, dSFC1, dSFC2, dV1, dV2, dGamma1, dGamma2, dRho1, dRho2, dS

  !Input/Output
  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: g, WfSeg, S
  double precision, dimension(0:numElem), intent(in) :: x, v, gamma
  double precision, dimension(0:numElem), intent(in) :: CT, SFC, rho
  double precision, dimension(0:numElem), intent(out) :: dCT1, dCT2, dSFC1, dSFC2
  double precision, dimension(0:numElem), intent(out) :: dV1, dV2, dGamma1, dGamma2, dRho1, dRho2, dS

  integer :: i = 0, j = 0, k = 0
  double precision :: deltax
  double precision :: dWfTempdSFC1, dWfTempdSFC2, dWfTempdCT1, dWfTempdCT2
  double precision :: dWfTempdV1, dWfTempdV2
  double precision :: dWfTempdGamma1, dWfTempdGamma2, dWfTempdRho1, dWfTempdRho2, dWfTempdS
  double precision, dimension(0:numInt-1) :: vTemp, xTemp, SFCTemp, gammaTemp, cosGamma, CTTemp, rhoTemp
  double precision, dimension(0:numInt-1) :: dVTemp1, dVTemp2, dSFCTemp1, dSFCTemp2, dGammaTemp1, dGammaTemp2
  double precision, dimension(0:numInt-1) :: dCosGamma1, dCosGamma2, dCTTemp1, dCTTemp2, dRhoTemp1, dRhoTemp2
  double precision, dimension(0:numInt-1) :: QTemp, dQTempdRho1, dQTempdRho2, dQTempdV1, dQTempdV2, dQTempdS

  do i = 0,numElem
     dSFC1(i) = 0.0
     dSFC2(i) = 0.0
     dCT1(i) = 0.0
     dCT2(i) = 0.0
     dV1(i) = 0.0
     dV2(i) = 0.0
     dGamma1(i) = 0.0
     dGamma2(i) = 0.0
     dRho1(i) = 0.0
     dRho2(i) = 0.0
     dS(i) = 0.0
  enddo

  do i = 0,numElem-1

     j = numElem-1-i
     call linspace(numInt, v(j), v(j+1), vTemp)
     call dlinspace(numInt, v(j), v(j+1), dVTemp1, dVTemp2)
     deltax = (x(j+1)-x(j))/real(numInt)

     do k = 0,numInt-1
        xTemp(k) = deltax
     enddo

     xTemp(0) = 0.5*deltax
     xTemp(numInt-1) = 0.5*deltax

     call linspace(numInt, SFC(j), SFC(j+1), SFCTemp)
     call dlinspace(numInt, SFC(j), SFC(j+1), dSFCTemp1, dSFCTemp2)
     call linspace(numInt, gamma(j), gamma(j+1), gammaTemp)
     call dlinspace(numInt, gamma(j), gamma(j+1), dGammaTemp1, dGammaTemp2)
     call linspace(numInt, CT(j), CT(j+1), CTTemp)
     call dlinspace(numInt, CT(j), CT(j+1), dCTTemp1, dCTTemp2)
     call linspace(numInt, rho(j), rho(j+1), rhoTemp)
     call dlinspace(numInt, rho(j), rho(j+1), dRhoTemp1, dRhoTemp2)
     QTemp = 0.5*rhoTemp*vTemp**2*S
     dQTempdRho1 = 0.5*vTemp**2*S*dRhoTemp1
     dQTempdRho2 = 0.5*vTemp**2*S*dRhoTemp2
     dQTempdV1 = rhoTemp*vTemp*S*dVTemp1
     dQTempdV2 = rhoTemp*vTemp*S*dVTemp2
     dQTempdS = 0.5*rhoTemp*vTemp**2
     cosGamma = cos(gammaTemp)
     dCosGamma1 = -sin(gammaTemp)*dGammaTemp1
     dCosGamma2 = -sin(gammaTemp)*dGammaTemp2

     do k = 0,numInt-1

        dWfTempdSFC1 = (CTTemp(k)*QTemp(k)*xTemp(k))/(vTemp(k)*cosGamma(k))
        dWfTempdSFC2 = dWfTempdSFC1
        dWfTempdSFC1 = dWfTempdSFC1*dSFCTemp1(k)
        dWfTempdSFC2 = dWfTempdSFC2*dSFCTemp2(k)
        dSFC1(j) = dSFC1(j)+dWfTempdSFC1*g
        dSFC2(j) = dSFC2(j)+dWfTempdSFC2*g

        dWfTempdCT1 = (SFCTemp(k)*QTemp(k)*xTemp(k))/(vTemp(k)*cosGamma(k))
        dWfTempdCT2 = dWfTempdCT1
        dWfTempdCT1 = dWfTempdCT1*dCTTemp1(k)
        dWfTempdCT2 = dWfTempdCT2*dCTTemp2(k)
        dCT1(j) = dCT1(j)+dWfTempdCT1*g
        dCT2(j) = dCT2(j)+dWfTempdCT2*g

        dWfTempdRho1 = (SFCTemp(k)*CTTemp(k)*xTemp(k))/(vTemp(k)*cosGamma(k))
        dWfTempdRho2 = dWfTempdRho1
        dWfTempdRho1 = dWfTempdRho1*dQTempdRho1(k)
        dWfTempdRho2 = dWfTempdRho2*dQTempdRho2(k)
        dRho1(j) = dRho1(j)+dWfTempdRho1*g
        dRho2(j) = dRho2(j)+dWfTempdRho2*g
        
        dWfTempdV1 = -(SFCTemp(k)*CTTemp(k)*QTemp(k)*xTemp(k))/(cosGamma(k)*vTemp(k)**2)
        dWfTempdV2 = dWfTempdV1
        dWfTempdV1 = dWfTempdV1*dVTemp1(k)
        dWfTempdV2 = dWfTempdV2*dVTemp2(k)
        dWfTempdV1 = dWfTempdV1 + (SFCTemp(k)*CTTemp(k)*dQTempdV1(k)*xTemp(k))/(vTemp(k)*cosGamma(k))
        dWfTempdV2 = dWfTempdV2 + (SFCTemp(k)*CTTemp(k)*dQTempdV2(k)*xTemp(k))/(vTemp(k)*cosGamma(k))
        dV1(j) = dV1(j)+dWfTempdV1*g
        dV2(j) = dV2(j)+dWfTempdV2*g

        dWfTempdGamma1 = -(SFCTemp(k)*CTTemp(k)*QTemp(k)*xTemp(k))/(vTemp(k)*cosGamma(k)**2)
        dWfTempdGamma2 = dWfTempdGamma1
        dWfTempdGamma1 = dWfTempdGamma1*dCosGamma1(k)
        dWfTempdGamma2 = dWfTempdGamma2*dCosGamma2(k)
        dGamma1(j) = dGamma1(j)+dWfTempdGamma1*g
        dGamma2(j) = dGamma2(j)+dWfTempdGamma2*g
        
        dWfTempdS = (SFCTemp(k)*CTTemp(k)*dQTempdS(k)*xTemp(k))/(vTemp(k)*cosGamma(k))
        dS(j) = dS(j)+dWfTempdS*g
     enddo
  enddo
       
end subroutine get_Wf_d

  ! compute the residuals for CL using the governing flight equation
  ! the flight equation is integrated, and the resultant value is 
  ! determined to be the residual for alpha
  
subroutine get_CL(numElem, numInt, Wac, S, g, x, v, rho, CL, Wf, gamma, &
     CT, alpha, CLRes)

  !f2py intent(in) numElem, numInt
  !f2py intent(in) Wac, S, g
  !f2py intent(in) x, v, rho, CL, Wf, gamma, CT, alpha
  !f2py intent(out) CLRes
  !f2py depend(numElem) x, v, rho, CL, Wf, gamma, CT, alpha

  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: Wac, S, g
  double precision, dimension(0:numElem), intent(in) :: x, v, rho, CL, Wf, gamma, CT, alpha
  double precision, dimension(0:numElem), intent(out) :: CLRes

  integer :: i = 0, j = 0, k = 0
  double precision :: temp, dx
  double precision, parameter :: param_zero = 0.0, param_one = 1.0
  double precision, dimension(0:numInt-1) :: R1, R2
  double precision, dimension(0:numInt-1) :: xTemp, vTemp, rhoTemp, QTemp
  double precision, dimension(0:numInt-1) :: gammaTemp, cosGamma
  double precision, dimension(0:numInt-1) :: CLTemp, WTemp, aTemp, tTemp
  double precision, dimension(0:numElem,0:numInt-1) :: alphaTemp

  do i = 0,numElem
     do j = 0,numInt-1
        alphaTemp(i,j) = 0.0
     enddo
     CLRes(i) = 0.0
  enddo

  call linspace(numInt, param_one, param_zero, R1)
  call linspace(numInt, param_zero, param_one, R2)

  do i = 0,numElem-1
     call linspace(numInt, rho(i), rho(i+1), rhoTemp)
     call linspace(numInt, v(i), v(i+1), vTemp)
     dx = (x(i+1)-x(i))/numInt

     do j = 0,numInt-1
        QTemp(j) = 0.5*rhoTemp(j)*vTemp(j)**2*S
        xTemp(j) = dx
     enddo

     xTemp(0) = 0.5*dx
     xTemp(numInt-1) = 0.5*dx

     call linspace(numInt, gamma(i), gamma(i+1), gammaTemp)
     cosGamma = cos(gammaTemp)

     call linspace(numInt, CL(i), CL(i+1), CLTemp)
     call linspace(numInt, Wf(i), Wf(i+1), WTemp)
     call linspace(numInt, alpha(i), alpha(i+1), aTemp)
     call linspace(numInt, CT(i), CT(i+1), tTemp)

     do j = 0,numInt-1
        WTemp(j) = WTemp(j) + Wac
     enddo

     do j = 0,numInt-1
        temp = -QTemp(j)*CLTemp(j)+WTemp(j)*cosGamma(j) &
             -sin(aTemp(j))*tTemp(j)*QTemp(j)
        temp = temp*R1(j)*xTemp(j)
        CLRes(i) = CLRes(i) + temp
        temp = -QTemp(j)*CLTemp(j)+WTemp(j)*cosGamma(j) &
             -sin(aTemp(j))*tTemp(j)*QTemp(j)
        temp = temp*R2(j)*xTemp(j)
        CLRes(i+1) = CLRes(i+1) + temp
     enddo
  enddo
end subroutine get_CL

subroutine get_CL_d(numElem, numInt, Wac, S, g, x, v, rho, CL, Wf, gamma, &
     CT, alpha, dCLdWac, dCLdS, dCLdV1, dCLdV2, dCLdV3, dCLdRho1, dCLdRho2, dCLdRho3, &
     dCLdCL1, dCLdCL2, dCLdCL3, dCLdWf1, dCLdWf2, dCLdWf3, dCLdGamma1, dCLdGamma2, dCLdGamma3, dCLdThrust1, &
     dCLdThrust2, dCLdThrust3, dCLdAlpha1, dCLdAlpha2, dCLdAlpha3)
  ! compute the derivatives of the residuals of CL wrt gamma,
  ! empty weight, S, v, rho, CL, fuel weight, gamma, thrust,
  ! alpha, dgamma/dx
  ! the jacobian structure is tri-diagonal, and is stored in 3
  ! vectors for each derivative

  !f2py intent(in) numElem, numInt
  !f2py intent(in) Wac, S, g
  !f2py intent(in) x, v, rho, CL, Wf, gamma, CT, alpha
  !f2py intent(out) dCLdWac, dCLdS, dCLdV1, dCLdV2, dCLdV3, dCLdRho1, dCLdRho2, dCLdRho3, dCLdCL1, dCLdCL2, dCLdCL3, dCLdGamma1, dCLdGamma2, dCLdGamma3, dCLdThrust1, dCLdThrust2, dCLdThrust3, dCLdAlpha1, dCLdAlpha2, dCLdAlpha3
  !f2py depend(numElem) x, v, rho, CL, Wf, gamma, CT, alpha, dCLdWac, dCLdS, dCLdV1, dCLdV2, dCLdV3, dCLdRho1, dCLdRho2, dCLdRho3, dCLdCL1, dCLdCL2, dCLdCL3, dCLdWf1, dCLdWf2, dCLdWf3, dCLdGamma1, dCLdGamma2, dCLdGamma3, dCLdThrust1, dCLdThrust2, dCLdThrust3, dCLdAlpha1, dCLdAlpha2, dCLdAlpha3

  ! Input/Output
  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: Wac, S, g
  double precision, dimension(0:numElem), intent(in) :: x, v, rho, CL, Wf, gamma, CT, alpha
  double precision, dimension(0:numElem), intent(out) :: dCLdWac, dCLdS, dCLdV1, dCLdV2, dCLdV3, dCLdRho1, dCLdRho2, dCLdRho3
  double precision, dimension(0:numElem), intent(out) :: dCLdCL1, dCLdCL2, dCLdCL3, dCLdWf1, dCLdWf2, dCLdWf3, dCLdGamma1
  double precision, dimension(0:numElem), intent(out) :: dCLdGamma2, dCLdGamma3
  double precision, dimension(0:numElem), intent(out) :: dCLdThrust1, dCLdThrust2, dCLdThrust3, dCLdAlpha1, dCLdAlpha2, dCLdAlpha3

  integer :: i = 0, j = 0, k = 0
  double precision :: temp, deltax
  double precision, parameter :: param_zero = 0.0, param_one = 1.0
  double precision, dimension(0:numInt-1) :: R1, R2
  double precision, dimension(0:numInt-1) :: xTemp, vTemp, rhoTemp, QTemp
  double precision, dimension(0:numInt-1) :: gammaTemp, cosGamma
  double precision, dimension(0:numInt-1) :: CLTemp, WTemp, aTemp, CTTemp
  double precision, dimension(0:numInt-1) :: dQTempdRho1, dQTempdRho2
  double precision, dimension(0:numInt-1) :: dQTempdV1, dQTempdV2
  double precision, dimension(0:numInt-1) :: dQTempdS
  double precision, dimension(0:numInt-1) :: dRhoTemp1, dRhoTemp2
  double precision, dimension(0:numInt-1) :: dVTemp1, dVTemp2
  double precision, dimension(0:numInt-1) :: dGammaTemp1, dGammaTemp2
  double precision, dimension(0:numInt-1) :: dCosGamma1, dCosGamma2
  double precision, dimension(0:numInt-1) :: dCLTemp1, dCLTemp2
  double precision, dimension(0:numInt-1) :: dWTempdWf1, dWTempdWf2
  double precision, dimension(0:numInt-1) :: dWTempdWac, dATemp1, dATemp2
  double precision, dimension(0:numInt-1) :: dCTTemp1, dCTTemp2

  call linspace(numInt, param_one, param_zero, R1)
  call linspace(numInt, param_zero, param_one, R2)

  do i = 0,numElem
     dCLdWac(i) = 0.0
     dCLdS(i) = 0.0
     dCLdV1(i) = 0.0
     dCLdV2(i) = 0.0
     dCLdV3(i) = 0.0
     dCLdRho1(i) = 0.0
     dCLdRho2(i) = 0.0
     dCLdRho3(i) = 0.0
     dCLdCL1(i) = 0.0
     dCLdCL2(i) = 0.0
     dCLdCL3(i) = 0.0
     dCLdWf1(i) = 0.0
     dCLdWf2(i) = 0.0
     dCLdWf3(i) = 0.0
     dCLdGamma1(i) = 0.0
     dCLdGamma2(i) = 0.0
     dCLdGamma3(i) = 0.0
     dCLdThrust1(i) = 0.0
     dCLdThrust2(i) = 0.0
     dCLdThrust3(i) = 0.0
     dCLdAlpha1(i) = 0.0
     dCLdAlpha2(i) = 0.0
     dCLdAlpha3(i) = 0.0
  enddo

  do i = 0,numElem-1
     call linspace(numInt, rho(i), rho(i+1), rhoTemp)
     call linspace(numInt, v(i), v(i+1), vTemp)
     call dlinspace(numInt, rho(i), rho(i+1), dRhoTemp1, dRhoTemp2)
     call dlinspace(numInt, v(i), v(i+1), dVTemp1, dVTemp2)
     deltax = (x(i+1)-x(i))/(numInt)
     
     do j = 0,numInt-1
        QTemp(j) = 0.5*rhoTemp(j)*vTemp(j)*vTemp(j)*S
        dQTempdRho1(j) = 0.5*vTemp(j)**2*S*dRhoTemp1(j)
        dQTempdRho2(j) = 0.5*vTemp(j)**2*S*dRhoTemp2(j)
        dQTempdV1(j) = rhoTemp(j)*vTemp(j)*dVTemp1(j)*S
        dQTempdV2(j) = rhoTemp(j)*vTemp(j)*dVTemp2(j)*S
        dQTempdS(j) = 0.5*rhoTemp(j)*vTemp(j)**2
        xTemp(j) = deltax
     enddo

     xTemp(0) = 0.5*deltax
     xTemp(numInt-1) = 0.5*deltax

     call linspace(numInt, gamma(i), gamma(i+1), gammaTemp)
     call dlinspace(numInt, gamma(i), gamma(i+1), dGammaTemp1, dGammaTemp2)
     cosGamma = cos(gammaTemp)
     dCosGamma1 = -sin(gammaTemp)*dGammaTemp1
     dCosGamma2 = -sin(gammaTemp)*dGammaTemp2

     call linspace(numInt, CL(i), CL(i+1), CLTemp)
     call dlinspace(numInt, CL(i), CL(i+1), dCLTemp1, dCLTemp2)
     call linspace(numInt, alpha(i), alpha(i+1), aTemp)
     call dlinspace(numInt, alpha(i), alpha(i+1), dATemp1, dATemp2)
     call linspace(numInt, CT(i), CT(i+1), CTTemp)
     call dlinspace(numInt, CT(i), CT(i+1), dCTTemp1, dCTTemp2)
     call linspace(numInt, Wf(i), Wf(i+1), WTemp)
     call dlinspace(numInt, Wf(i), Wf(i+1), dWTempdWf1, dWTempdWf2)
     
     do j = 0,numInt-1
        WTemp(j) = WTemp(j) + Wac
        dWTempdWac(j) = 1.0
     enddo

     do j = 0,numInt-1
        dCLdWac(i) = dCLdWac(i)+R1(j)*xTemp(j)*(cosGamma(j))
        dCLdS(i) = dCLdS(i)-R1(j)*xTemp(j)*(dQTempdS(j)*CLTemp(j) &
             +SIN(aTemp(j))*CTTemp(j)*dQTempdS(j))
        dCLdV2(i) = dCLdV2(i)+R1(j)*xTemp(j)*(-dQTempdV1(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdV1(j))
        dCLdV3(i) = dCLdV3(i)+R1(j)*xTemp(j)*(-dQTempdV2(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdV2(j))
        dCLdRho2(i) = dCLdRho2(i)+R1(j)*xTemp(j)*(-dQTempdRho1(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdRho1(j))
        dCLdRho3(i) = dCLdRho3(i)+R1(j)*xTemp(j)*(-dQTempdRho2(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdRho2(j))
        dCLdCL2(i) = dCLdCL2(i)+R1(j)*xTemp(j)*(-QTemp(j)*dCLTemp1(j))
        dCLdCL3(i) = dCLdCL3(i)+R1(j)*xTemp(j)*(-QTemp(j)*dCLTemp2(j))
        dCLdWf2(i) = dCLdWf2(i)+R1(j)*xTemp(j)*(dWTempdWf1(j)*cosGamma(j))
        dCLdWf3(i) = dCLdWf3(i)+R1(j)*xTemp(j)*(dWTempdWf2(j)*cosGamma(j))
        dCLdGamma2(i) = dCLdGamma2(i)+R1(j)*xTemp(j)*(WTemp(j)*dCosGamma1(j))
        dCLdGamma3(i) = dCLdGamma3(i)+R1(j)*xTemp(j)*(WTemp(j)*dCosGamma2(j))
        dCLdThrust2(i) = dCLdThrust2(i)+R1(j)*xTemp(j)*(-sin(aTemp(j)) &
             *dCTTemp1(j)*QTemp(j))
        dCLdThrust3(i) = dCLdThrust3(i)+R1(j)*xTemp(j)*(-sin(aTemp(j)) &
             *dCTTemp2(j)*QTemp(j))
        dCLdAlpha2(i) = dCLdAlpha2(i)+R1(j)*xTemp(j)*(-cos(aTemp(j)) &
             *CTTemp(j)*dATemp1(j)*QTemp(j))
        dCLdAlpha3(i) = dCLdAlpha3(i)+R1(j)*xTemp(j)*(-cos(aTemp(j)) &
             *CTTemp(j)*dATemp2(j)*QTemp(j))

        dCLdWac(i+1) = dCLdWac(i+1)+R2(j)*xTemp(j)*(cosGamma(j))
        dCLdS(i+1) = dCLdS(i+1)+R2(j)*xTemp(j)*(-dQTempdS(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdS(j))
        dCLdV1(i+1) = dCLdV1(i+1)+R2(j)*xTemp(j)*(-dQTempdV1(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdV1(j))
        dCLdV2(i+1) = dCLdV2(i+1)+R2(j)*xTemp(j)*(-dQTempdV2(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdV2(j))
        dCLdRho1(i+1) = dCLdRho1(i+1)+R2(j)*xTemp(j)*(-dQTempdRho1(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdRho1(j))
        dCLdRho2(i+1) = dCLdRho2(i+1)+R2(j)*xTemp(j)*(-dQTempdRho2(j)*CLTemp(j) &
             -SIN(aTemp(j))*CTTemp(j)*dQTempdRho2(j))
        dCLdCL1(i+1) = dCLdCL1(i+1)+R2(j)*xTemp(j)*(-QTemp(j)*dCLTemp1(j))
        dCLdCL2(i+1) = dCLdCL2(i+1)+R2(j)*xTemp(j)*(-QTemp(j)*dCLTemp2(j))
        dCLdWf1(i+1) = dCLdWf1(i+1)+R2(j)*xTemp(j)*(dWTempdWf1(j)*cosGamma(j))
        dCLdWf2(i+1) = dCLdWf2(i+1)+R2(j)*xTemp(j)*(dWTempdWf2(j)*cosGamma(j))
        dCLdGamma1(i+1) = dCLdGamma1(i+1)+R2(j)*xTemp(j)*(WTemp(j)*dCosGamma1(j))
        dCLdGamma2(i+1) = dCLdGamma2(i+1)+R2(j)*xTemp(j)*(WTemp(j)*dCosGamma2(j))
        dCLdThrust1(i+1) = dCLdThrust1(i+1)+R2(j)*xTemp(j)*(-sin(aTemp(j)) &
             *dCTTemp1(j)*QTemp(j))
        dCLdThrust2(i+1) = dCLdThrust2(i+1)+R2(j)*xTemp(j)*(-sin(aTemp(j)) &
             *dCTTemp2(j)*QTemp(j))
        dCLdAlpha1(i+1) = dCLdAlpha1(i+1)+R2(j)*xTemp(j)*(-cos(aTemp(j)) &
             *CTTemp(j)*dATemp1(j)*QTemp(j))
        dCLdAlpha2(i+1) = dCLdAlpha2(i+1)+R2(j)*xTemp(j)*(-cos(aTemp(j)) &
             *CTTemp(j)*dATemp2(j)*QTemp(j))
     enddo
  enddo
end subroutine get_CL_d

subroutine get_CM(numElem, numInt, S, chord, x, v, rho, CM, CMRes)
  ! compute the residuals of CM by integrating the moment flight equation

  !f2py intent(in) numElem, numInt
  !f2py intent(in) S, chord
  !f2py intent(in) x, v, rho, CM
  !f2py intent(out) CMRes
  !f2py depend(numElem) x, v, rho, CM, CMRes

  !Input/Output
  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: S, chord
  double precision, dimension(0:numElem), intent(in) :: x, v, rho
  double precision, dimension(0:numElem), intent(in) :: CM
  double precision, dimension(0:numElem), intent(out) :: CMRes

  integer :: i = 0, j = 0, k = 0
  double precision :: temp, dx
  double precision, parameter :: param_zero = 0.0, param_one = 1.0
  double precision, dimension(0:numInt-1) :: R1, R2
  double precision, dimension(0:numInt-1) :: xTemp, vTemp, rhoTemp, QTemp
  double precision, dimension(0:numInt-1) :: CMTemp
  double precision, dimension(0:numElem,0:numInt-1) :: etaTemp

  do i = 0,numElem
     do j = 0,numInt-1
        etaTemp(i,j) = 0.0
     enddo
     CMRes(i) = 0.0
  enddo

  call linspace(numInt, param_one, param_zero, R1)
  call linspace(numInt, param_zero, param_one, R2)
  
  do i = 0,numElem-1
     call linspace(numInt, rho(i), rho(i+1), rhoTemp)
     call linspace(numInt, v(i), v(i+1), vTemp)
     dx = (x(i+1)-x(i))/numInt

     do j = 0,numInt-1
        QTemp(j) = 0.5*rhoTemp(j)*vTemp(j)*vTemp(j)*S
        xTemp(j) = dx
     enddo

     xTemp(0) = 0.5*dx
     xTemp(numInt-1) = 0.5*dx

     do j = 0,numInt-1
        CMTemp(j) = CM(i)*R1(j) + CM(i+1)*R2(j)
     enddo

     do j = 0,numInt-1
        temp = QTemp(j)*chord*CMTemp(j)
        temp = temp*R1(j)*xTemp(j)
        etaTemp(i,j) = etaTemp(i,j) + temp
     enddo
  enddo

  do i = 0,numElem-1
     call linspace(numInt, rho(i), rho(i+1), rhoTemp)
     call linspace(numInt, v(i), v(i+1), vTemp)
     dx = (x(i+1)-x(i))/numInt
     
     do j = 0,numInt-1
        QTemp(j) = 0.5*rhoTemp(j)*vTemp(j)*vTemp(j)*S
        xTemp(j) = dx
     enddo

     xTemp(0) = 0.5*dx
     xTemp(numInt-1) = 0.5*dx
     
     do j = 0,numInt-1
        CMTemp(j) = CM(i)*R1(j) + CM(i+1)*R2(j)
     enddo

     do j = 0,numInt-1
        temp = QTemp(j)*chord*CMTemp(j)
        temp = temp*R2(j)*xTemp(j)
        etaTemp(i+1,j) = etaTemp(i+1,j) + temp
     enddo
  enddo

  do i = 0,numElem
     do j = 0,numInt-1
        CMRes(i) = CMRes(i) + etaTemp(i,j)
     enddo
  enddo

end subroutine get_CM

subroutine get_CM_d(numElem, numInt, S, chord, x, v, rho, &
     CM, dCMdS, dCMdC, dCMdV1, dCMdV2, dCMdV3, dCMdRho1, dCMdRho2, dCMdRho3, &
     dCMdCM1, dCMdCM2, dCMdCM3)
  ! compute the derivatives of residuals of CM wrt S, chord,
  ! inertia, v, rho, gamma, dgamma/dx, d2gamma/dx2, dv/dx, CM
  ! the jacobians are tri-diagonal, and are stored in 3 vectors
  
  !f2py intent(in) numElem, numInt
  !f2py intent(in) S, chord
  !f2py intent(in) x, v, rho, CM
  !f2py intent(out) dCMdS, dCMdC, dCMdV1, dCMdV2, dCMdV3, dCMdRho1, dCMdRho2, dCMdRho3, dCMdCM1, dCMdCM2, dCMdCM3
  !f2py depend(numElem) x, v, rho, CM, dCMdS, dCMdC, dCMdV1, dCMdV2, dCMdV3, dCMdRho1, dCMdRho2, dCMdRho3, dCMdCM1, dCMdCM2, dCMdCM3

  !Input/Output
  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: S, chord
  double precision, dimension(0:numElem), intent(in) :: x, v, rho, CM
  double precision, dimension(0:numElem), intent(out) :: dCMdS, dCMdC
  double precision, dimension(0:numElem), intent(out) :: dCMdV1, dCMdV2, dCMdV3, dCMdRho1
  double precision, dimension(0:numElem), intent(out) :: dCMdRho2, dCMdRho3
  double precision, dimension(0:numElem), intent(out) :: dCMdCM1, dCMdCM2, dCMdCM3

  integer :: i = 0, j = 0, k = 0
  double precision :: deltax
  double precision, parameter :: param_zero = 0.0, param_one = 1.0
  double precision, dimension(0:numInt-1) :: R1, R2
  double precision, dimension(0:numInt-1) :: xTemp, vTemp, rhoTemp, QTemp
  double precision, dimension(0:numInt-1) :: CMTemp
  double precision, dimension(0:numInt-1) :: dRhoTemp1, dRhoTemp2, dVTemp1
  double precision, dimension(0:numInt-1) :: dVTemp2, dQTempdRho1, dQTempdRho2
  double precision, dimension(0:numInt-1) :: dQTempdV1, dQTempdV2, dQTempdS
  double precision, dimension(0:numInt-1) :: dCMTemp1, dCMTemp2
  
  do i = 0,numElem
     dCMdS(i) = 0.0
     dCMdC(i) = 0.0
     dCMdV1(i) = 0.0
     dCMdV2(i) = 0.0
     dCMdV3(i) = 0.0
     dCMdRho1(i) = 0.0
     dCMdRho2(i) = 0.0
     dCMdRho3(i) = 0.0
     dCMdCM1(i) = 0.0
     dCMdCM2(i) = 0.0
     dCMdCM3(i) = 0.0
  enddo

  call linspace(numInt, param_one, param_zero, R1)
  call linspace(numInt, param_zero, param_one, R2)

  do i = 0,numElem-1
     call linspace(numInt, rho(i), rho(i+1), rhoTemp)
     call linspace(numInt, v(i), v(i+1), vTemp)
     call dlinspace(numInt, rho(i), rho(i+1), dRhoTemp1, dRhoTemp2)
     call dlinspace(numInt, v(i), v(i+1), dVTemp1, dVTemp2)
     deltax = (x(i+1)-x(i))/numInt
     
     do j = 0,numInt-1
        QTemp(j) = 0.5*rhoTemp(j)*vTemp(j)**2*S
        dQTempdRho1(j) = 0.5*dRhoTemp1(j)*vTemp(j)**2*S
        dQTempdRho2(j) = 0.5*dRhoTemp2(j)*vTemp(j)**2*S
        dQTempdV1(j) = rhoTemp(j)*vTemp(j)*dVTemp1(j)*S
        dQTempdV2(j) = rhoTemp(j)*vTemp(j)*dVTemp2(j)*S
        dQTempdS(j) = 0.5*rhoTemp(j)*vTemp(j)**2
        xTemp(j) = deltax
     enddo

     xTemp(0) = 0.5*deltax
     xTemp(numInt-1) = 0.5*deltax

     do j = 0,numInt-1
        CMTemp(j) = CM(i)*R1(j) + CM(i+1)*R2(j)
        dCMTemp1(j) = R1(j)
        dCMTemp2(j) = R2(j)
     enddo

     do j = 0,numInt-1
        dCMdS(i) = dCMdS(i)+R1(j)*xTemp(j)*(dQTempdS(j)*chord*CMTemp(j))
        dCMdC(i) = dCMdC(i)+R1(j)*xTemp(j)*(QTemp(j)*CMTemp(j))
        dCMdV2(i) = dCMdV2(i)+R1(j)*xTemp(j)*(dQTempdV1(j)*chord*CMTemp(j))
        dCMdV3(i) = dCMdV3(i)+R1(j)*xTemp(j)*(dQTempdV2(j)*chord*CMTemp(j))
        dCMdRho2(i) = dCMdRho2(i)+R1(j)*xTemp(j)*(dQTempdRho1(j)*chord*CMTemp(j))
        dCMdRho3(i) = dCMdRho3(i)+R1(j)*xTemp(j)*(dQTempdRho2(j)*chord*CMTemp(j))
        dCMdCM2(i) = dCMdCM2(i)+R1(j)*xTemp(j)*QTemp(j)*chord*dCMTemp1(j)
        dCMdCM3(i) = dCMdCM3(i)+R1(j)*xTemp(j)*QTemp(j)*chord*dCMTemp2(j)

        dCMdS(i+1) = dCMdS(i+1)+R2(j)*xTemp(j)*(dQTempdS(j)*chord*CMTemp(j))
        dCMdC(i+1) = dCMdC(i+1)+R2(j)*xTemp(j)*(QTemp(j)*CMTemp(j))
        dCMdV1(i+1) = dCMdV1(i+1)+R2(j)*xTemp(j)*(dQTempdV1(j)*chord*CMTemp(j))
        dCMdV2(i+1) = dCMdV2(i+1)+R2(j)*xTemp(j)*(dQTempdV2(j)*chord*CMTemp(j))
        dCMdRho1(i+1) = dCMdRho1(i+1)+R2(j)*xTemp(j)*(dQTempdRho1(j)*chord*CMTemp(j))
        dCMdRho2(i+1) = dCMdRho2(i+1)+R2(j)*xTemp(j)*(dQTempdRho2(j)*chord*CMTemp(j))
        dCMdCM1(i+1) = dCMdCM1(i+1)+R2(j)*xTemp(j)*QTemp(j)*chord*dCMTemp1(j)
        dCMdCM2(i+1) = dCMdCM2(i+1)+R2(j)*xTemp(j)*QTemp(j)*chord*dCMTemp2(j)
     enddo
  enddo

end subroutine get_CM_d

subroutine get_h_des(numElem,numInt,h_dot,x_ends,v_ends,h_ends,h)
  
  !f2py intent(in) numElem, numInt
  !f2py intent(in) h_dot
  !f2py intent(in) x_ends, v_ends, h_ends
  !f2py intent(out) h
  !f2py depend(numElem) h
  
  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: h_dot
  double precision, dimension(0:1), intent(in) :: x_ends, v_ends, h_ends
  double precision, dimension(0:numElem), intent(out) :: h

  integer :: i, j
  double precision :: dx, dv, dh
  double precision, dimension(0:numElem) :: x, v
  double precision, dimension(0:numInt) :: xTemp

  dx = (x_ends(1)-x_ends(0))/numElem
  dv = (v_ends(1)-v_ends(0))/numElem

  do i = 0,numElem
     x(i) = i*dx+x_ends(0)
     v(i) = i*dv+v_ends(0)
     h(i) = 0
  enddo

  h(0) = h_ends(0)
  do i = 0,numElem-1
     h(i+1) = h(i)
     call linspace(numInt+1,x(i),x(i+1),xTemp)
     do j = 0,numInt-1
        dh = h_dot/SQRT((v(i)+(v(i+1)-v(i))*xTemp(j)/(x(i+1)-x(i)))**2-h_dot**2)
        h(i+1) = h(i+1)+dh
     enddo
     h(i+1) = h(i+1)*(x(i+1)-x(i))/numInt
  enddo
endsubroutine get_h_des

subroutine get_h_des_d(numElem,numInt,h_dot,x_ends,v_ends,h_ends,dh_dhEnds0,dh_dhEnds1, &
     dh_dhDot, dh_dxEnds0, dh_dxEnds1, dh_dvEnds0, dh_dvEnds1)

  !f2py intent(in) numElem, numInt
  !f2py intent(in) h_dot
  !f2py intent(in) x_ends, v_ends, h_ends
  !f2py intent(out) dh_dhEnds0, dh_dhEnds1, dh_dhDot, dh_dxEnds0, dh_dxEnds1, dh_dvEnds0, dh_dvEnds1
  !f2py depend(numElem) dh_dhEnds0, dh_dhEnds1, dh_dhDot, dh_dxEnds0, dh_dxEnds1, dh_dvEnds0, dh_dvEnds1

  integer, intent(in) :: numElem, numInt
  double precision, intent(in) :: h_dot
  double precision, dimension(0:1), intent(in) :: x_ends, v_ends, h_ends
  double precision, dimension(0:numElem), intent(out) :: dh_dhEnds0, &
       dh_dhEnds1, dh_dhDot, dh_dxEnds0, dh_dxEnds1, dh_dvEnds0, dh_dvEnds1

  integer :: i, j
  double precision :: dx, dv, sqrtTerm, delta
  double precision :: ddh_dhDot, ddh_dh0, ddh_dh1, ddh_dx0, ddh_dx1
  double precision :: ddh_dv0, ddh_dv1
  double precision, dimension(0:numElem) :: x, v, h
  double precision, dimension(0:numElem) :: dx_dx0, dx_dx1, dv_dv0, dv_dv1
  double precision, dimension(0:numInt) :: xTemp, dxTemp1, dxTemp2

  dx = (x_ends(1)-x_ends(0))/numElem
  dv = (v_ends(1)-v_ends(0))/numElem

  do i = 0,numElem
     dx_dx0(i) = 1 - i/numElem
     dx_dx1(i) = i/numElem
     dv_dv0(i) = 1 - i/numElem
     dv_dv1(i) = i/numElem
     x(i) = i*dx+x_ends(0)
     v(i) = i*dv+v_ends(0)
     dh_dhEnds0(i) = 1
     dh_dhEnds1(i) = 0
     dh_dhDot(i) = 0
     dh_dxEnds0(i) = 0
     dh_dxEnds1(i) = 0
     dh_dvEnds0(i) = 0
     dh_dvEnds1(i) = 0
  enddo

  do i = 0,numElem-1

     call dlinspace(numInt+1,x(i),x(i+1),dxTemp1,dxTemp2)
     do j = 0,numInt-1
        sqrtTerm = SQRT((v(i)+(v(i+1)-v(i))*xTemp(j)/(x(i+1)-x(i)))**2-h_dot**2)
        ddh_dhDot = 1/sqrtTerm
        ddh_dhDot = ddh_dhDot + h_dot**2/(sqrtTerm**3)
        dh_dhDot(i+1) = dh_dhDot(i+1) + ddh_dhDot

        ddh_dx0 = (v(i+1)-v(i))*(dxTemp1(j)*dx_dx0(i)+dxTemp2(j)*dx_dx0(i+1))/(x(i+1)-x(i))
        ddh_dx0 = ddh_dx0 - ((v(i+1)-v(i))*xTemp(j)/(x(i+1)-x(i))**2)*(dx_dx0(i+1)-dx_dx0(i))
        ddh_dx0 = ddh_dx0 * (v(i)+(v(i+1)-v(i))*xTemp(j)/(x(i+1)-x(i)))
        ddh_dx0 = ddh_dx0 * (-h_dot/(sqrtTerm**3))
        dh_dxEnds0(i+1) = dh_dxEnds0(i+1) + ddh_dx0
        
        ddh_dx1 = (v(i+1)-v(i))*(dxTemp1(j)*dx_dx1(i)+dxTemp2(j)*dx_dx1(i+1))/(x(i+1)-x(i))
        ddh_dx1 = ddh_dx1 - ((v(i+1)-v(i))*xTemp(j)/(x(i+1)-x(i))**2)*(dx_dx1(i+1)-dx_dx1(i))
        ddh_dx1 = ddh_dx1 * (v(i)+(v(i+1)-v(i))*xTemp(j)/(x(i+1)-x(i)))
        ddh_dx1 = ddh_dx1 * (-h_dot/(sqrtTerm**3))
        dh_dxEnds1(i+1) = dh_dxEnds1(i+1) + ddh_dx1

        ddh_dv0 = (xTemp(j)/(x(i+1)-x(i)))*(dv_dv0(i+1)-dv_dv0(i))
        ddh_dv0 = ddh_dv0 + dv_dv0(i)
        ddh_dv0 = ddh_dv0 * (v(i)+(v(i+1)-v(i))*xTemp(j)/(x(i+1)-x(i)))
        ddh_dv0 = ddh_dv0 * (-h_dot/(sqrtTerm**3))
        dh_dvEnds0(i+1) = dh_dvEnds0(i+1) + ddh_dv0
        
        ddh_dv1 = (xTemp(j)/(x(i+1)-x(i)))*(dv_dv1(i+1)-dv_dv1(i))
        ddh_dv1 = ddh_dv1 + dv_dv1(i)
        ddh_dv1 = ddh_dv1 * (v(i)+(v(i+1)-v(i))*xTemp(j)/(x(i+1)-x(i)))
        ddh_dv1 = ddh_dv1 * (-h_dot/(sqrtTerm**3))
        dh_dvEnds1(i+1) = dh_dvEnds1(i+1) + ddh_dv1

        dh = h_dot/sqrtTerm
        h(i+1) = h(i+1)+dh
     enddo

     delta = (x(i+1)-x(i))/numInt
     dh_dhDot(i+1) = dh_dhDot(i+1) * delta
     dh_dvEnds0(i+1) = dh_dvEnds0(i+1) * delta
     dh_dvEnds1(i+1) = dh_dvEnds1(i+1) * delta
     dh_dxEnds0(i+1) = h(i+1)*(dx_dx0(i+1)-dx_dx0(i))/numInt + dh_dxEnds0(i+1)*delta
     dh_dxEnds1(i+1) = h(i+1)*(dx_dx1(i+1)-dx_dx1(i))/numInt + dh_dxEnds1(i+1)*delta
  enddo
endsubroutine get_h_des_d

