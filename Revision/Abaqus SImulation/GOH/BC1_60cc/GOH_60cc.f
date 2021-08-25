!****************************************************************************
C Column Numbers:
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
c...  ------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      dimension statev(nstatv)

      statev(1)=1.0d0    ! lam1g
      statev(2)=1.0d0    ! lam1e
      statev(3)=1.0d0    ! lam1
      statev(4)=1.0d0    ! lam2g
      statev(5)=1.0d0    ! lam2e
      statev(6)=1.0d0    ! lam2

      statev(7)=1.0d0    ! xn(1)
      statev(8)=1.0d0    ! xn(2)
      statev(9)=1.0d0    ! xn(3)
      statev(10)=1.0d0    ! f1(1)
      statev(11)=1.0d0    ! f1(2)
      statev(12)=1.0d0    ! f1(3)

      statev(13)=1.0d0    ! f0(1)
      statev(14)=1.0d0    ! f0(2)
      statev(15)=1.0d0    ! f0(3)
      statev(16)=1.0d0    ! s0(1)
      statev(17)=1.0d0    ! s0(2)
      statev(18)=1.0d0    ! s0(3)
      statev(19)=1.0d0    ! s1(1)
      statev(20)=1.0d0    ! s1(2)
      statev(21)=1.0d0    ! s1(3)

      return
      end

c...  ------------------------------------------------------------------
c     user amplitude subroutine
      Subroutine UAMP(
     *     ampName, time, ampValueOld, dt, nProps, props, nSvars, svars,
     *     lFlagsInfo, nSensor, sensorValues, sensorNames,
     *     jSensorLookUpTable,
C          to be defined
     *     ampValueNew,
     *     lFlagsDefine,
     *     AmpDerivative, AmpSecDerivative, AmpIncIntegral,
     *     AmpIncDoubleIntegral)

      include 'aba_param.inc'

C     svars - additional state variables, similar to (V)UEL

C     time indices
      parameter (iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
C     flags passed in for information
      parameter (iInitialization   = 1,
     *           iRegularInc       = 2,
     *           nFlagsInfo        = 2)
C     optional flags to be defined
      parameter (iComputeDeriv       = 1,
     *           iComputeSecDeriv    = 2,
     *           iComputeInteg       = 3,
     *           iComputeDoubleInteg = 4,
     *           iStopAnalysis       = 5,
     *           iConcludeStep       = 6,
     *           nFlagsDefine        = 6)

      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine)
      dimension jSensorLookUpTable(*)
      common /pressure/ porepres
      common /porepres_pre/ porepres_pre
      common /volume/ porevol
      common /preKiterm/ preKiterm
      common /preerr/ preerr
      common /tStart_before/ tStart_before

      double precision porepres,porevol,preKiterm,
     #                 preerr,tStart_before,porepres_pre
      dimension sensorValues(nSensor), svars(nSvars), props(nProps)

      character*80 sensorNames(nSensor)
      character*80 ampName

c     get sensor values first
      porevol = GetSensorValue('VOL_EXP',jSensorLookUpTable,
     *                           sensorValues)
      porepres_pre = GetSensorValue('PRES_EXP',jSensorLookUpTable,
     *                           sensorValues)

      ampValueOld = 0.00003

      if (ampName(1:7) .eq. 'AMP_PID' ) then

C        User code to compute  ampValue = F(sensors)
         if (lFlagsInfo(iInitialization).eq.1) then
c	    tim = time(iStepTime)
            tStart_before = time(iStepTime)
            ampValueNew  = ampValueOld
            porepres = ampValueOld

         else
c           Example: f(t) = t

            tim = time(iStepTime)
            tStart = tim - dt
            tEnd   = tim
            call pid_control(dt,tim,tStart)
            ampValueNew = porepres


C           stop the analysis, if desired
            if (rfAve.gt.1000000) lFlagsDefine(iStopAnalysis)=1
         end if

      end if

      return
      end

c...  ------------------------------------------------------------------
      subroutine pid_control(dt,tim,tStart)
      implicit none
      common /pressure/ porepres
      common /porepres_pre/ porepres_pre
      common /volume/ porevol
      common /preKiterm/ preKiterm
      common /preerr/ preerr
      common /tStart_before/ tStart_before


      real*8  tim,tStart, tStart_before
      real*8  prevol, curvol, err, SP, Kp, Ki, Kd, dt
      real*8  output,porepres,porevol,preKiterm,preerr,porepres_pre


      curvol = porevol


      SP = 55000.d0

c      Kp = 2.5d0
c      Ki = 2.0d0
c      Kd = 0.001d0

      Kp = 5.0d0 ! 3 works
      Ki = 10.0d0
      Kd = 0.001d0 ! 0.001 works

      err = SP - curvol

      if (tStart_before.eq.tStart) then
        porevol = porevol
        tStart_before = tStart
      else


        output = Kp*err + preKiterm + Ki*dt*err+ Kd*(err-preerr)/dt

        if (tim.gt.0.01d0) then
          porepres= porepres_pre + output/1000000000.d0 ! 200000000.d0 works
        else
          porepres= porepres_pre + output/1000000000.d0
        end if

        print *,'output,preerr,err,err_preerr,preKiterm,Kp*err,Ki*dt*err,Kd*err-preerr/dt',
     #           output,preerr,err,err-preerr,preKiterm,Kp*err,Ki*dt*err,Kd*(err-preerr)/dt

        porevol = curvol + output ! new volume
        preKiterm = preKiterm + Ki*dt*err ! new Ki term
        preerr = err ! new error term

        tStart_before = tStart

      end if


      RETURN
      END
c...  ------------------------------------------------------------------


c...  ------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     #rpl,ddsddt,drplde,drpldt,
     #stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     #ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     #celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     #ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     #stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     #props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      call umat_area_stretch(stress,statev,ddsdde,sse,
     #                       time,dtime,coords,props,dfgrd1,
     #                       ntens,ndi,nshr,nstatv,nprops,
     #                       noel,npt,kstep,kinc)

      return
      end

c...  ------------------------------------------------------------------
      subroutine umat_area_stretch(stress,statev,ddsdde,sse,
     #                             time,dtime,coords,props,dfgrd1,
     #                             ntens,ndi,nshr,nstatv,nprops,
     #                             noel,npt,kstep,kinc)
c...  ------------------------------------------------------------------
c * " This UMAT is written for a GOH constitutive material.
c * " It implements area-stretch-driven area growth in the plane normal to xn0.
c * " It is not suitable for use with local material directions?

c * " Written by Maria Holland in Aug 2012
c * " Updated by Maria Holland in Sep 2017
c * " Edited by Taeksang Lee in Mar 2019
c...  ------------------------------------------------------------------

      implicit none

c...  variables to be defined
      real*8  stress(ntens), ddsdde(ntens,ntens), statev(nstatv), sse

c...  variables passed in for information
      real*8  time(2), dtime, coords(3), props(nprops), dfgrd1(3,3)
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt, kstep, kinc

c...  material properties
      real*8 blk, mu, kappa, k1, k2, xn0(3), f0(3), s0(3)
      real*8 alpha, tcr, tmax, gam
      real*8 kg1, kg2, mg1, mg2, ng1, ng2
      real*8 lamg1, lamg2, kkg1, kkg2, tcr1, tcr2
      real*8 tcrg1, tcrg2, phig, dot_lamg1, dot_lamg2
      real*8 dot_tlamg1_dlamg1, dot_tlamg2_dlamg2
      real*8 dot_lamg1_dlamg1, dot_lamg2_dlamg2


c...  local variables
ccccc indices, counters, integer stuff
      integer ii, jj, kk, ll, nitl, iii, jjj
ccccc deformation tensors (F,C,Fe,Fg,b,Cinv,Ceinv...)
      real*8 finv(3,3), fe(3,3), be(6), b(6), ce(6)
      real*8  cinv(6), cinvmat(3,3), ceinv(6), feinv(3,3)
      real*8 ceinvmat(3,3), cemat(3,3)
      real*8 fg(6), fginv(6), fginvmat(3,3), fgmat(3,3), c(6)
ccccc scalars derived from deformation, and growth variables
      real*8 detfe, norm, the, theg, theg_n, detf, lnJe, detce
      real*8 Ebar, Ebar_el, E, E_el
      real*8  I1, I1bar, I1_el, I1bar_el, I4, I4bar, I4_el, I4bar_el
      real*8 lam1, lam1g, lam1g_n, lam1e, lam2, lam2g, lam2g_n, lam2e
ccccc stress and tangent tensors, including growth tangent tensors
      real*8 stressmat(3,3),  ddsedde(6,6), ddsgdde(6,6)
      real*8  pk2_isc_iso(3,3), pk2_isc_aniso(3,3), pk2_vol(3,3)
      real*8 pk2(3,3)
      real*8 cg_ij
      real*8 cg_kl(6), cg_ij_1(6), cg_kl_1(6), cg_ij_2(6), cg_kl_2(6)
      real*8 cg_ijmat(3,3), cg_ijmat_1(3,3), cg_ijmat_2(3,3)
      real*8 dFgdtheg(3,3)
ccccc fiber direction stuff
      real*8 xn(3), ftn0(3), normf, f1(3), ff0(6), ff1(6)
      real*8 xn0xn0(3,3), fbar(3), fbar_e(3), ffbar(6), ff0mat(3,3)
      real*8 ss0(6), ss0mat(3,3), s1(3), ss1(6)
      real*8 H(6)
ccccc local newton iteration stuff
      real*8  res1, dres1, res2, dres2, fac1, fac2, xtol
ccccc temporal variables, no clear name
      real*8  coef1, coef2, coef3, coef4, coef5, coef6
      real*8  tmp1(3,3), tmp2(3,3), tmp3(3,3)
      real*8 const1, const2(6), const2mat(3,3)
      real*8  tmp1_1(3,3), tmp2_1(3,3), tmp3_1(3,3)
      real*8  tmp1_2(3,3), tmp2_2(3,3), tmp3_2(3,3)
      real*8 const1_1, const2_1(6), const2mat_1(3,3)
      real*8 const1_2, const2_2(6), const2mat_2(3,3)
ccccc constants
      real*8  xi(6), ximat(3,3), xixi(6,6)


ccccc ###############################################
ccccc ABT: what is this data? -> TL: 2nd order identity tensor
ccccc in Voigt form. I think the name of data does not have
ccccc meaning of the actual data, but just initial definition.
ccccc ###############################################
c      print *, 'dfgrd1', dfgrd1(1,1), dfgrd1(1,2), dfgrd1(1,3)
c      print *, 'dfgrd1', dfgrd1(2,1), dfgrd1(2,2), dfgrd1(2,3)
c      print *, 'dfgrd1', dfgrd1(3,1), dfgrd1(3,2), dfgrd1(3,3)

      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
      xtol = 1.d-12

c...  initialize material parameters
      blk    = props(1)   ! Bulk modulus
      mu     = props(2)   ! shear modulus
      kappa = props(3)    ! kappa, fiber dispersion: 0 to 1/3
      k1     = props(4)   ! nonlinear constant in goh model
      k2     = props(5)   ! nonlinear constant in goh model
      f0(1) = props(6)    ! initial fiber orientation[1]
      f0(2) = props(7)    ! initial fiber orientation[2]
      f0(3) = props(8)    ! initial fiber orientation[3]
      xn0(1) = props(9)   ! n0[1]
      xn0(2) = props(10)  ! n0[2]
      xn0(3) = props(11)  ! n0[3]
      kkg1 = props(12)       ! hill function parameter for direction 1
      kkg2 = props(13)       ! hill function parameter for direction 2
      mg1 = props(14)       ! hill function parameter for direction 1
      mg2 = props(15)       ! hill function parameter for direction 2
      ng1=props(16)         ! hill function parameter for direction 1
      ng2=props(17)         ! hill function parameter for direction 2
      tcr1    = props(18)  ! critical stretch 1
      tcr2    = props(19)  ! critical stretch 2

ccccc ###############################################
ccccc calculate s0; cross product(n0,f0)
ccccc s0 is orthogonal to fiber and normal, CHECK
      s0(1) = xn0(2)*f0(3)-xn0(3)*f0(2)
      s0(2) = xn0(3)*f0(1)-xn0(1)*f0(3)
      s0(3) = xn0(1)*f0(2)-xn0(2)*f0(1)
      statev(16) = s0(1)
      statev(17) = s0(2)
      statev(18) = s0(3)

c      s0(1) = 0.d0
c      s0(2) = 1.d0
c      s0(3) = 0.d0
ccccc ###############################################

c...  make xi to matrix
      call move6to33(xi,ximat)

c...  calculate 4th order identity
      xixi(1,1) = 0.5d0*(xi(1)*xi(1) + xi(1)*xi(1))
      xixi(2,2) = 0.5d0*(xi(2)*xi(2) + xi(2)*xi(2))
      xixi(3,3) = 0.5d0*(xi(3)*xi(3) + xi(3)*xi(3))
      xixi(1,2) = 0.5d0*(xi(4)*xi(4) + xi(4)*xi(4))
      xixi(1,3) = 0.5d0*(xi(5)*xi(5) + xi(5)*xi(5))
      xixi(1,4) = 0.5d0*(xi(1)*xi(4) + xi(4)*xi(1))
      xixi(1,5) = 0.5d0*(xi(1)*xi(5) + xi(5)*xi(1))
      xixi(1,6) = 0.5d0*(xi(4)*xi(5) + xi(5)*xi(4))
      xixi(2,4) = 0.5d0*(xi(4)*xi(2) + xi(2)*xi(4))
      xixi(2,5) = 0.5d0*(xi(4)*xi(6) + xi(6)*xi(4))
      xixi(2,6) = 0.5d0*(xi(2)*xi(6) + xi(6)*xi(2))
      xixi(3,4) = 0.5d0*(xi(5)*xi(6) + xi(6)*xi(5))
      xixi(3,5) = 0.5d0*(xi(5)*xi(3) + xi(3)*xi(5))
      xixi(3,6) = 0.5d0*(xi(6)*xi(3) + xi(3)*xi(6))
      xixi(4,4) = 0.5d0*(xi(1)*xi(2) + xi(4)*xi(4))
      xixi(4,5) = 0.5d0*(xi(1)*xi(6) + xi(5)*xi(4))
      xixi(4,6) = 0.5d0*(xi(4)*xi(6) + xi(5)*xi(2))
      xixi(5,5) = 0.5d0*(xi(1)*xi(3) + xi(5)*xi(5))
      xixi(5,6) = 0.5d0*(xi(4)*xi(3) + xi(5)*xi(6))
      xixi(6,6) = 0.5d0*(xi(2)*xi(3) + xi(6)*xi(6))

c...  calculate deformed elastic normal
      xn(1) = (dfgrd1(1,1)*xn0(1)+dfgrd1(1,2)*xn0(2)+dfgrd1(1,3)*xn0(3))
      xn(2) = (dfgrd1(2,1)*xn0(1)+dfgrd1(2,2)*xn0(2)+dfgrd1(2,3)*xn0(3))
      xn(3) = (dfgrd1(3,1)*xn0(1)+dfgrd1(3,2)*xn0(2)+dfgrd1(3,3)*xn0(3))
c      norm = sqrt(xn(1)*xn(1) + xn(2)*xn(2) + xn(3)*xn(3))

c...  update state variable history
      statev(7) = xn(1)
      statev(8) = xn(2)
      statev(9) = xn(3)

      statev(13) = f0(1)
      statev(14) = f0(2)
      statev(15) = f0(3)

c... calculate xn0 otimes xn0
      do ii = 1,3
        do jj = 1,3
          xn0xn0(ii,jj) = xn0(ii)*xn0(jj)
        end do
      end do

c...  ff0 = f0 \otimes f0
      ff0(1) = f0(1)*f0(1)
      ff0(2) = f0(2)*f0(2)
      ff0(3) = f0(3)*f0(3)
      ff0(4) = f0(1)*f0(2)
      ff0(5) = f0(1)*f0(3)
      ff0(6) = f0(2)*f0(3)

c     ...

c...  ss0 = s0 \otimes s0
      ss0(1) = s0(1)*s0(1)
      ss0(2) = s0(2)*s0(2)
      ss0(3) = s0(3)*s0(3)
      ss0(4) = s0(1)*s0(2)
      ss0(5) = s0(1)*s0(3)
      ss0(6) = s0(2)*s0(3)

c...  ff0mat
      ff0mat(1,1) = ff0(1)
      ff0mat(2,2) = ff0(2)
      ff0mat(3,3) = ff0(3)
      ff0mat(1,2) = ff0(4)
      ff0mat(1,3) = ff0(5)
      ff0mat(2,3) = ff0(6)
      ff0mat(2,1) = ff0(4)
      ff0mat(3,1) = ff0(5)
      ff0mat(3,2) = ff0(6)

c...  ss0mat
      ss0mat(1,1) = ss0(1)
      ss0mat(2,2) = ss0(2)
      ss0mat(3,3) = ss0(3)
      ss0mat(1,2) = ss0(4)
      ss0mat(1,3) = ss0(5)
      ss0mat(2,3) = ss0(6)
      ss0mat(2,1) = ss0(4)
      ss0mat(3,1) = ss0(5)
      ss0mat(3,2) = ss0(6)

c...  structure tensor, symmetric
      do ii = 1,6
        H(ii) = kappa*xi(ii) + (1.d0 - 3.d0*kappa)*ff0(ii)
      end do

c...  calculate determinant of deformation gradient
      detf = +dfgrd1(1,1)*(dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #       -dfgrd1(1,2)*(dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #       +dfgrd1(1,3)*(dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))

c...  calculate area stretch
      finv(1,1) = (+dfgrd1(2,2)*dfgrd1(3,3) - dfgrd1(2,3)*dfgrd1(3,2))/detf
      finv(1,2) = (-dfgrd1(1,2)*dfgrd1(3,3) + dfgrd1(1,3)*dfgrd1(3,2))/detf
      finv(1,3) = (+dfgrd1(1,2)*dfgrd1(2,3) - dfgrd1(1,3)*dfgrd1(2,2))/detf
      finv(2,1) = (-dfgrd1(2,1)*dfgrd1(3,3) + dfgrd1(2,3)*dfgrd1(3,1))/detf
      finv(2,2) = (+dfgrd1(1,1)*dfgrd1(3,3) - dfgrd1(1,3)*dfgrd1(3,1))/detf
      finv(2,3) = (-dfgrd1(1,1)*dfgrd1(2,3) + dfgrd1(1,3)*dfgrd1(2,1))/detf
      finv(3,1) = (+dfgrd1(2,1)*dfgrd1(3,2) - dfgrd1(2,2)*dfgrd1(3,1))/detf
      finv(3,2) = (-dfgrd1(1,1)*dfgrd1(3,2) + dfgrd1(1,2)*dfgrd1(3,1))/detf
      finv(3,3) = (+dfgrd1(1,1)*dfgrd1(2,2) - dfgrd1(1,2)*dfgrd1(2,1))/detf

      ftn0(1) = finv(1,1)*xn0(1)+finv(2,1)*xn0(2)+finv(3,1)*xn0(3)
      ftn0(2) = finv(1,2)*xn0(1)+finv(2,2)*xn0(2)+finv(3,2)*xn0(3)
      ftn0(3) = finv(1,3)*xn0(1)+finv(2,3)*xn0(2)+finv(3,3)*xn0(3)

      the = detf*sqrt(ftn0(1)*ftn0(1) + ftn0(2)*ftn0(2) + ftn0(3)*ftn0(3))

c...  calculate deformed fiber orientation: f1 = F*f0
      f1(1) = dfgrd1(1,1)*f0(1)
     #      + dfgrd1(1,2)*f0(2)
     #      + dfgrd1(1,3)*f0(3)
      f1(2) = dfgrd1(2,1)*f0(1)
     #      + dfgrd1(2,2)*f0(2)
     #      + dfgrd1(2,3)*f0(3)
      f1(3) = dfgrd1(3,1)*f0(1)
     #      + dfgrd1(3,2)*f0(2)
     #      + dfgrd1(3,3)*f0(3)

c...  calculate deformed s0 orientation: s1 = F*s0
      s1(1) = dfgrd1(1,1)*s0(1)
     #      + dfgrd1(1,2)*s0(2)
     #      + dfgrd1(1,3)*s0(3)
      s1(2) = dfgrd1(2,1)*s0(1)
     #      + dfgrd1(2,2)*s0(2)
     #      + dfgrd1(2,3)*s0(3)
      s1(3) = dfgrd1(3,1)*s0(1)
     #      + dfgrd1(3,2)*s0(2)
     #      + dfgrd1(3,3)*s0(3)

      statev(19) = s1(1)
      statev(20) = s1(2)
      statev(21) = s1(3)

c...  ff1 = f1 \otimes f1
      ff1(1) = f1(1)*f1(1)
      ff1(2) = f1(2)*f1(2)
      ff1(3) = f1(3)*f1(3)
      ff1(4) = f1(1)*f1(2)
      ff1(5) = f1(1)*f1(3)
      ff1(6) = f1(2)*f1(3)

c...  ss1 = s1 \otimes s1
      ss1(1) = s1(1)*s1(1)
      ss1(2) = s1(2)*s1(2)
      ss1(3) = s1(3)*s1(3)
      ss1(4) = s1(1)*s1(2)
      ss1(5) = s1(1)*s1(3)
      ss1(6) = s1(2)*s1(3)

c...  update state variable history
      statev(10) = f1(1)
      statev(11) = f1(2)
      statev(12) = f1(3)

ccccc ###############################################
ccccc calculate lam1 and lam2, total stretch, CHECK
      lam1 = sqrt(ff1(1)+ff1(2)+ff1(3))
      lam2 = sqrt(ss1(1)+ss1(2)+ss1(3))
ccccc ###############################################
c      print *, 'lam2',lam2
c...  obtain state variable history
      lam1g_n = statev(1)
      lam1g  = lam1g_n
      lam2g_n = statev(4)
      lam2g  = lam2g_n


c...  ------------------------------------------------------------------
c...  local newton iteration lambda1
c...  ------------------------------------------------------------------
      nitl = 0

      phig = (lam1/lam1g) - tcr1
c      print *, 'lam1',lam1
c      print *, 'lam1g',lam1g
c      print *, 'tcr1',tcr1

      if (phig.gt.0) then
  200    continue
         if (time(2).gt.0) then
           continue
           nitl = nitl + 1
           dot_lamg1 = kkg1*(lam1/lam1g-tcr1)**ng1/(mg1**ng1+(lam1/lam1g-tcr1)**ng1)
           dot_lamg1_dlamg1 = (kkg1*ng1*(lam1/lam1g-tcr1)**(ng1-1.d0)*(mg1**ng1+(lam1/lam1g-tcr1)**ng1)
     #     -kkg1*(lam1/lam1g-tcr1)**ng1*ng1*(lam1/lam1g-tcr1)**(ng1-1.d0))
     #      /(mg1**ng1 + (lam1/lam1g-tcr1)**ng1)**2.d0
     #     *(-lam1/lam1g**2.d0)

           res1  = lam1g - lam1g_n - dot_lamg1*dtime
           dres1 = 1.d0 - dot_lamg1_dlamg1*dtime ! =K_theta
           lam1g = lam1g - res1 / dres1

           if ((nitl.lt.20).and.(dabs(res1).gt.xtol)) go to 200
           if (nitl.eq.20) print *, 'no local convergence! |r|=',dabs(res1)

           fac1 = -dtime/dres1*dot_lamg1_dlamg1/lam1/lam1g
         else
           lam1g = lam1g_n
         end if
      else
        lam1g = lam1g_n
        fac1 =0.d0
      end if
c...  ------------------------------------------------------------------
c...  local newton iteration lambda2
c...  ------------------------------------------------------------------
      nitl = 0

      phig = (lam2/lam2g) - tcr2

      if (phig.gt.0) then
  201    continue
         if (time(2).gt.0) then
           continue
           nitl = nitl + 1
           dot_lamg2 = kkg2*(lam2/lam2g-tcr2)**ng2/(mg2**ng2+(lam2/lam2g-tcr2)**ng2)
           dot_lamg2_dlamg2 = (kkg2*ng2*(lam2/lam2g-tcr2)**(ng2-1.d0)*(mg2**ng2+(lam2/lam2g-tcr2)**ng2)
     #     -kkg2*(lam2/lam2g-tcr2)**ng2*ng2*(lam2/lam2g-tcr2)**(ng2-1.d0))
     #      /(mg2**ng2 + (lam2/lam2g-tcr2)**ng2)**2.d0
     #     *(-lam2/lam2g**2.d0)

           res2  = lam2g - lam2g_n - dot_lamg2*dtime
           dres2 = 1.d0 - dot_lamg2_dlamg2*dtime ! =K_theta
           lam2g = lam2g - res2 / dres2

           if ((nitl.lt.20).and.(dabs(res2).gt.xtol)) go to 201
           if (nitl.eq.20) print *, 'no local convergence! |r|=',dabs(res2)

           fac2 = -dtime/dres2*dot_lamg2_dlamg2/lam2/lam2g
         else
           lam2g=lam2g_n
         end if
      else
        lam2g = lam2g_n
        fac2 =0.d0
      end if
c...  ------------------------------------------------------------------
c...  end local newton iteration
c      print *, 'res2,dres2,dot_lamg2,dot_lamg2_dlamg2', res2,dres2,dot_lamg2,dot_lamg2_dlamg2

c...  update state variables
      statev(1) = lam1g
      statev(2) = lam1/lam1g
      statev(3) = lam1
      statev(4) = lam2g
      statev(5) = lam2/lam2g
      statev(6) = lam2


ccccc ###################################################
ccccc compute new growth tensor Fg, its inverse, and Fe
ccccc CHECK
      fg(1) = lam1g*ff0(1) + lam2g*ss0(1) + xn0xn0(1,1)
      fg(2) = lam1g*ff0(2) + lam2g*ss0(2) + xn0xn0(2,2)
      fg(3) = lam1g*ff0(3) + lam2g*ss0(3) + xn0xn0(3,3)
      fg(4) = lam1g*ff0(4) + lam2g*ss0(4) + xn0xn0(1,2)
      fg(5) = lam1g*ff0(5) + lam2g*ss0(5) + xn0xn0(1,3)
      fg(6) = lam1g*ff0(6) + lam2g*ss0(6) + xn0xn0(2,3)
      call move6to33(fg,fgmat)
      fginv(1) = 1.d0/lam1g*ff0(1) + 1.d0/lam2g*ss0(1) + xn0xn0(1,1)
      fginv(2) = 1.d0/lam1g*ff0(2) + 1.d0/lam2g*ss0(2) + xn0xn0(2,2)
      fginv(3) = 1.d0/lam1g*ff0(3) + 1.d0/lam2g*ss0(3) + xn0xn0(3,3)
      fginv(4) = 1.d0/lam1g*ff0(4) + 1.d0/lam2g*ss0(4) + xn0xn0(1,2)
      fginv(5) = 1.d0/lam1g*ff0(5) + 1.d0/lam2g*ss0(5) + xn0xn0(1,3)
      fginv(6) = 1.d0/lam1g*ff0(6) + 1.d0/lam2g*ss0(6) + xn0xn0(2,3)
      call move6to33(fginv,fginvmat)
      fe(1,1)= dfgrd1(1,1)*fginv(1)+ dfgrd1(1,2)*fginv(4)+ dfgrd1(1,3)*fginv(5)
      fe(1,2)= dfgrd1(1,1)*fginv(4)+ dfgrd1(1,2)*fginv(2)+ dfgrd1(1,3)*fginv(6)
      fe(1,3)= dfgrd1(1,1)*fginv(5)+ dfgrd1(1,2)*fginv(6)+ dfgrd1(1,3)*fginv(3)
      fe(2,1)= dfgrd1(2,1)*fginv(1)+ dfgrd1(2,2)*fginv(4)+ dfgrd1(2,3)*fginv(5)
      fe(2,2)= dfgrd1(2,1)*fginv(4)+ dfgrd1(2,2)*fginv(2)+ dfgrd1(2,3)*fginv(6)
      fe(2,3)= dfgrd1(2,1)*fginv(5)+ dfgrd1(2,2)*fginv(6)+ dfgrd1(2,3)*fginv(3)
      fe(3,1)= dfgrd1(3,1)*fginv(1)+ dfgrd1(3,2)*fginv(4)+ dfgrd1(3,3)*fginv(5)
      fe(3,2)= dfgrd1(3,1)*fginv(4)+ dfgrd1(3,2)*fginv(2)+ dfgrd1(3,3)*fginv(6)
      fe(3,3)= dfgrd1(3,1)*fginv(5)+ dfgrd1(3,2)*fginv(6)+ dfgrd1(3,3)*fginv(3)
c      fe(1,1)= 1.d0/lam1g*f1(1)*f0(1) +1.d0/lam2g*s1(1)*s0(1) +xn(1)*xn0(1)
c      fe(1,2)= 1.d0/lam1g*f1(1)*f0(2) +1.d0/lam2g*s1(1)*s0(2) +xn(1)*xn0(2)
c      fe(1,3)= 1.d0/lam1g*f1(1)*f0(3) +1.d0/lam2g*s1(1)*s0(3) +xn(1)*xn0(3)
c      fe(2,1)= 1.d0/lam1g*f1(2)*f0(1) +1.d0/lam2g*s1(2)*s0(1) +xn(2)*xn0(1)
c      fe(2,2)= 1.d0/lam1g*f1(2)*f0(2) +1.d0/lam2g*s1(2)*s0(2) +xn(2)*xn0(2)
c      fe(2,3)= 1.d0/lam1g*f1(2)*f0(3) +1.d0/lam2g*s1(2)*s0(3) +xn(2)*xn0(3)
c      fe(3,1)= 1.d0/lam1g*f1(3)*f0(1) +1.d0/lam2g*s1(3)*s0(1) +xn(3)*xn0(1)
c      fe(3,2)= 1.d0/lam1g*f1(3)*f0(2) +1.d0/lam2g*s1(3)*s0(2) +xn(3)*xn0(2)
c      fe(3,3)= 1.d0/lam1g*f1(3)*f0(3) +1.d0/lam2g*s1(3)*s0(3) +xn(3)*xn0(3)
ccccc ###################################################

c     Je
c...  calculate determinant of Fe
      detfe = +fe(1,1) * (fe(2,2)*fe(3,3)-fe(2,3)*fe(3,2))
     #        -fe(1,2) * (fe(2,1)*fe(3,3)-fe(2,3)*fe(3,1))
     #        +fe(1,3) * (fe(2,1)*fe(3,2)-fe(2,2)*fe(3,1))
      lnJe = dlog(detfe)

c     fiber in the intermediate configuration
c     this fiber should be normalized because cauchy stress
c     is calculated from the intermediate configuration.
c     It is kind of new initial fiber normalized in the intermediate
c.    Continuous change of fbar in the intermediate configuration makes big issue
      fbar(1) = f0(1)
      fbar(2) = f0(2)
      fbar(3) = f0(3)

c      fbar(1) = fg(1)*f0(1) + fg(4)*f0(2) + fg(5)*f0(3)
c      fbar(2) = fg(4)*f0(1) + fg(2)*f0(2) + fg(6)*f0(3)
c      fbar(3) = fg(5)*f0(1) + fg(6)*f0(2) + fg(3)*f0(3)
 
c      fbar(1) = fbar(1)/sqrt(fbar(1)**2.d0 + fbar(2)**2.d0 + fbar(3)**2.d0)
c      fbar(2) = fbar(2)/sqrt(fbar(1)**2.d0 + fbar(2)**2.d0 + fbar(3)**2.d0)
c      fbar(3) = fbar(3)/sqrt(fbar(1)**2.d0 + fbar(2)**2.d0 + fbar(3)**2.d0)

c     fiber in the current configuration
c     It should not be normalized to capture the elastic deformation
c     along the fiber and itself.

      fbar_e(1) = fe(1,1)*fbar(1) + fe(1,2)*fbar(2) + fe(1,3)*fbar(3)
      fbar_e(2) = fe(2,1)*fbar(1) + fe(2,2)*fbar(2) + fe(2,3)*fbar(3)
      fbar_e(3) = fe(3,1)*fbar(1) + fe(3,2)*fbar(2) + fe(3,3)*fbar(3)

c...  ffbar = fbar \otimes fbar
      ffbar(1) = fbar_e(1)*fbar_e(1)
      ffbar(2) = fbar_e(2)*fbar_e(2)
      ffbar(3) = fbar_e(3)*fbar_e(3)
      ffbar(4) = fbar_e(1)*fbar_e(2)
      ffbar(5) = fbar_e(1)*fbar_e(3)
      ffbar(6) = fbar_e(2)*fbar_e(3)

c     C^-1
c...  calculate inverse of total right cauchy-green deformation tensor
      cinv(1) = finv(1,1)*finv(1,1)+finv(1,2)*finv(1,2)+finv(1,3)*finv(1,3)
      cinv(2) = finv(2,1)*finv(2,1)+finv(2,2)*finv(2,2)+finv(2,3)*finv(2,3)
      cinv(3) = finv(3,1)*finv(3,1)+finv(3,2)*finv(3,2)+finv(3,3)*finv(3,3)
      cinv(4) = finv(1,1)*finv(2,1)+finv(1,2)*finv(2,2)+finv(1,3)*finv(2,3)
      cinv(5) = finv(1,1)*finv(3,1)+finv(1,2)*finv(3,2)+finv(1,3)*finv(3,3)
      cinv(6) = finv(2,1)*finv(3,1)+finv(2,2)*finv(3,2)+finv(2,3)*finv(3,3)

      call move6to33(cinv,cinvmat)

c...  Ce
c...  calculate elastic right cauchy-green deformation tensor ce = fe^t * fe
      ce(1) = fe(1,1)*fe(1,1) + fe(2,1)*fe(2,1) + fe(3,1)*fe(3,1)
      ce(2) = fe(1,2)*fe(1,2) + fe(2,2)*fe(2,2) + fe(3,2)*fe(3,2)
      ce(3) = fe(1,3)*fe(1,3) + fe(2,3)*fe(2,3) + fe(3,3)*fe(3,3)
      ce(4) = fe(1,1)*fe(2,1) + fe(2,1)*fe(2,2) + fe(3,1)*fe(2,3)
      ce(5) = fe(1,1)*fe(3,1) + fe(2,1)*fe(3,2) + fe(3,1)*fe(3,3)
      ce(6) = fe(1,2)*fe(3,1) + fe(2,2)*fe(3,2) + fe(3,2)*fe(3,3)

c.... cemat
      call move6to33(ce,cemat)

      detce = +ce(1) * (ce(2)*ce(3)-ce(6)*ce(6))
     #        -ce(4) * (ce(4)*ce(3)-ce(6)*ce(5))
     #        +ce(5) * (ce(4)*ce(6)-ce(2)*ce(5))

c...  C
c...  calculate right cauchy-green deformation tensor c = f^t * f
      c(1) = dfgrd1(1,1)*dfgrd1(1,1)+dfgrd1(2,1)*dfgrd1(2,1)
     #       +dfgrd1(3,1)*dfgrd1(3,1)
      c(2) = dfgrd1(1,2)*dfgrd1(1,2)+dfgrd1(2,2)*dfgrd1(2,2)
     #       +dfgrd1(3,2)*dfgrd1(3,2)
      c(3) = dfgrd1(1,3)*dfgrd1(1,3)+dfgrd1(2,3)*dfgrd1(2,3)
     #        +dfgrd1(3,3)*dfgrd1(3,3)
      c(4) = +dfgrd1(1,1)*dfgrd1(1,2)
     #       +dfgrd1(2,1)*dfgrd1(2,2)
     #       +dfgrd1(3,1)*dfgrd1(3,2)
      c(5) = dfgrd1(1,1)*dfgrd1(1,3)
     #      +dfgrd1(2,1)*dfgrd1(2,3)
     #      +dfgrd1(3,1)*dfgrd1(3,3)
      c(6) = dfgrd1(1,2)*dfgrd1(1,3)
     #      +dfgrd1(2,2)*dfgrd1(2,3)
     #      +dfgrd1(3,2)*dfgrd1(3,3)

c...  calculate elastic left cauchy-green deformation tensor be = fe * fe^t
      be(1) = fe(1,1)*fe(1,1) + fe(1,2)*fe(1,2) + fe(1,3)*fe(1,3)
      be(2) = fe(2,1)*fe(2,1) + fe(2,2)*fe(2,2) + fe(2,3)*fe(2,3)
      be(3) = fe(3,1)*fe(3,1) + fe(3,2)*fe(3,2) + fe(3,3)*fe(3,3)
      be(4) = fe(1,1)*fe(2,1) + fe(1,2)*fe(2,2) + fe(1,3)*fe(2,3)
      be(5) = fe(1,1)*fe(3,1) + fe(1,2)*fe(3,2) + fe(1,3)*fe(3,3)
      be(6) = fe(2,1)*fe(3,1) + fe(2,2)*fe(3,2) + fe(2,3)*fe(3,3)

c...  calculate left cauchy-green deformation tensor b = f * f^t
      b(1) = dfgrd1(1,1)*dfgrd1(1,1)
     #     + dfgrd1(1,2)*dfgrd1(1,2)
     #     + dfgrd1(1,3)*dfgrd1(1,3)
      b(2) = dfgrd1(2,1)*dfgrd1(2,1)
     #     + dfgrd1(2,2)*dfgrd1(2,2)
     #     + dfgrd1(2,3)*dfgrd1(2,3)
      b(3) = dfgrd1(3,1)*dfgrd1(3,1)
     #     + dfgrd1(3,2)*dfgrd1(3,2)
     #     + dfgrd1(3,3)*dfgrd1(3,3)
      b(4) = dfgrd1(1,1)*dfgrd1(2,1)
     #     + dfgrd1(1,2)*dfgrd1(2,2)
     #     + dfgrd1(1,3)*dfgrd1(2,3)
      b(5) = dfgrd1(1,1)*dfgrd1(3,1)
     #     + dfgrd1(1,2)*dfgrd1(3,2)
     #     + dfgrd1(1,3)*dfgrd1(3,3)
      b(6) = dfgrd1(2,1)*dfgrd1(3,1)
     #     + dfgrd1(2,2)*dfgrd1(3,2)
     #     + dfgrd1(2,3)*dfgrd1(3,3)

c... Fe^-1
c    Fe^-1 = Fg*F^-1
      feinv = matmul(fgmat,finv)

c... Ce^-1 = fe^-1 * fe^-t
c... calculate inverse of elastic right cauchy-green deformation tensor

      ceinv(1) = feinv(1,1)*feinv(1,1)
     #         + feinv(1,2)*feinv(1,2)
     #         + feinv(1,3)*feinv(1,3)
      ceinv(2) = feinv(2,1)*feinv(2,1)
     #         + feinv(2,2)*feinv(2,2)
     #         + feinv(2,3)*feinv(2,3)
      ceinv(3) = feinv(3,1)*feinv(3,1)
     #         + feinv(3,2)*feinv(3,2)
     #         + feinv(3,3)*feinv(3,3)
      ceinv(4) = feinv(1,1)*feinv(2,1)
     #         + feinv(1,2)*feinv(2,2)
     #         + feinv(1,3)*feinv(2,3)
      ceinv(5) = feinv(1,1)*feinv(3,1)
     #         + feinv(1,2)*feinv(3,2)
     #         + feinv(1,3)*feinv(3,3)
      ceinv(6) = feinv(2,1)*feinv(3,1)
     #         + feinv(2,2)*feinv(3,2)
     #         + feinv(2,3)*feinv(3,3)

c...  convert viogt to full matrix of cinv
      call move6to33(ceinv,ceinvmat)

c     stretch square in the intermediate configuration: |Fe*f0|^2; strecth^2 or I4
c     I1 and I4
      I1 = c(1) + c(2) + c(3)
      I1bar = detf**(-2.d0/3.d0)*I1

      I1_el = ce(1) + ce(2) + ce(3)
      I1bar_el = detfe**(-2.d0/3.d0)*I1_el

      I4 = ff1(1) + ff1(2) + ff1(3)
      I4bar = detf**(-2.d0/3.d0)*I4

      I4_el = fbar_e(1)*fbar_e(1)+fbar_e(2)*fbar_e(2)+fbar_e(3)*fbar_e(3)
      I4bar_el = detfe**(-2.d0/3.d0)*I4_el

c...  Ebar = kappa*(I1-3.d0) + (1.d0 - 3.d0*kappa)*(I4-1.d0)
      E = kappa*(I1-3.d0) + (1.d0-3.d0*kappa)*(I4-1.d0)
      Ebar = kappa*(I1bar-3.d0) + (1.d0-3.d0*kappa)*(I4bar-1.d0)

c...  E_el = kappa*I1_el + (1-3*kappa)*I4_el
      E_el = kappa*(I1_el-3.d0) + (1.d0-3.d0*kappa)*(I4_el-1.d0)
      Ebar_el = kappa*(I1bar_el-3.d0) + (1.d0-3.d0*kappa)*(I4bar_el-1.d0)

c...  only tension can generate contribution of fiber
      if (E.lt.0) then
        E = 0.d0
      end if

      if (Ebar.lt.0) then
        Ebar = 0.d0
      end if

      if (E_el.lt.0) then
        E_el = 0.d0
      end if

      if (Ebar_el.lt.0) then
        Ebar_el = 0.d0
      end if

c...  calculate PK2 stress
c     pk2 = pk2_vol + pk2_isc
c         = pk2_vol + pk2_isc_iso + pk2_isc_aniso

      pk2_vol = blk/2.d0*(detfe**2.d0-1.d0)*ceinvmat

      pk2_isc_iso = mu/detfe**(2.d0/3.d0)
     #            *(ximat - I1/3.d0*ceinvmat)

      pk2_isc_aniso = 2.d0*kappa*matmul(fginvmat,transpose(fginvmat))
     #              +2.d0*(1.d0-3.d0*kappa)
     #              *matmul(matmul(fginvmat,ff0mat),transpose(fginvmat))
     #              -2.d0/3.d0*kappa*I1_el*cinvmat
     #              -2.d0/3.d0*(1.d0-3.d0*kappa)*I4_el*cinvmat

      pk2_isc_aniso = k1*Ebar_el*exp(k2*Ebar_el**2.d0)*detfe**(-2.d0/3.d0)
     #              *pk2_isc_aniso

      pk2 = detf/detfe*(pk2_vol + pk2_isc_iso + pk2_isc_aniso)


      coef1 = blk*detfe**2.d0 + 2.d0/9.d0*mu*detfe**(-2.d0/3.d0)*I1_el
     #        +4.d0/9.d0*k1*exp(k2*Ebar_el**2.d0)*detfe**(-4.d0/3.d0)
     #          *(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el)**2.d0
     #        +8.d0/9.d0*k1*k2*Ebar_el**2.d0*exp(k2*Ebar_el**2.d0)
     #        *detfe**(-4.d0/3.d0)
     #          *(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el)**2.d0
     #        +4.d0/9.d0*k1*Ebar_el*exp(k2*Ebar_el**2.d0)*detfe**(-2.d0/3.d0)
     #          *(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el)
      coef1 = coef1/detfe

      coef2 = -2.d0/3.d0*mu*detfe**(-2.d0/3.d0)
      coef2 = coef2/detfe

      coef3 = -blk/2.d0*(detfe**2.d0-1.d0) + mu/3.d0*detfe**(-2.d0/3.d0)*I1_el
     #        +2.d0/3.d0*k1*Ebar_el*exp(k2*Ebar_el**2.d0)*detfe**(-2.d0/3.d0)
     #          *(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el)
      coef3 = coef3/detfe

      coef4 = -4.d0/3.d0*k1*exp(k2*Ebar_el**2.d0)
     #        *detfe**(-4.d0/3.d0)*(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el)
     #        -8.d0/3.d0*k1*k2*Ebar_el**2.d0*exp(k2*Ebar_el**2.d0)
     #          *detfe**(-4.d0/3.d0)*(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el)
     #        -4.d0/3.d0*k1*Ebar_el*exp(k2*Ebar_el**2.d0)*detfe**(-2.d0/3.d0)
      coef4 = coef4/detfe

      coef5 = 4.d0*k1*exp(k2*Ebar_el**2.d0)*detfe**(-4.d0/3.d0)
     #        +8.d0*k1*k2*Ebar_el**2.d0*exp(k2*Ebar_el**2.d0)
     #        *detfe**(-4.d0/3.d0)
      coef5 = coef5/detfe

      coef6 = -4.d0/3.d0*k1*Ebar_el*exp(k2*Ebar_el**2.d0)*detfe**(-2.d0/3.d0)
     #        -4.d0/3.d0*k1*exp(k2*Ebar_el**2.d0)*detfe**(-4.d0/3.d0)
     #          *(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el)
     #        -8.d0/3.d0*k1*k2*Ebar_el**2.d0*exp(k2*Ebar_el**2.d0)
     #        *detfe**(-4.d0/3.d0)
     #          *(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el)
      coef6 = coef6/detfe


c.... when fiber does not have tension.
      if (Ebar_el.eq.0.d0) then
        coef1 = blk*detfe**2.d0 + 2.d0/9.d0*mu*detfe**(-2.d0/3.d0)*I1_el
        coef1 = coef1/detfe

        coef2 = -2.d0/3.d0*mu*detfe**(-2.d0/3.d0)
        coef2 = coef2/detfe

        coef3 = -blk/2.d0*(detfe**2.d0-1.d0) + mu/3.d0*detfe**(-2.d0/3.d0)*I1_el
        coef3 = coef3/detfe

        coef4 = 0.d0
        coef4 = coef4/detfe

        coef5 = 0.d0
        coef5 = coef5/detfe

        coef6 = 0.d0
        coef6 = coef6/detfe
      end if
c      print *, 'coef1, coef6', coef1, coef6


c...  calculate elastic tangent, ddsedde
      ddsedde(1,1) = coef1*xi(1)*xi(1)
     # + coef2*be(1)*xi(1) + coef2*xi(1)*be(1)
     #              +coef3*2.d0*xixi(1,1)
     #     +coef4*xi(1)*(kappa*be(1) + (1.d0-3.d0*kappa)*fbar_e(1)*fbar_e(1))
     #     +coef5*(kappa*be(1) + (1.d0-3.d0*kappa)*fbar_e(1)*fbar_e(1))
     #    *(kappa*be(1) + (1.d0-3.d0*kappa)*fbar_e(1)*fbar_e(1))
     #    +coef6*(kappa*be(1) + (1.d0-3.d0*kappa)*fbar_e(1)*fbar_e(1))*xi(1)

      ddsedde(2,2) = coef1*xi(2)*xi(2)
     # + coef2*be(2)*xi(2) + coef2*xi(2)*be(2)
     #              +coef3*2.d0*xixi(2,2)
     #     +coef4*xi(2)*(kappa*be(2) + (1.d0-3.d0*kappa)*fbar_e(2)*fbar_e(2))
     #   +coef5*(kappa*be(2) + (1.d0-3.d0*kappa)*fbar_e(2)*fbar_e(2))
     #   *(kappa*be(2) + (1.d0-3.d0*kappa)*fbar_e(2)*fbar_e(2))
     #   +coef6*(kappa*be(2) + (1.d0-3.d0*kappa)*fbar_e(2)*fbar_e(2))*xi(2)

      ddsedde(3,3) = coef1*xi(3)*xi(3)
     # + coef2*be(3)*xi(3) + coef2*xi(3)*be(3)
     #              +coef3*2.d0*xixi(3,3)
     #    +coef4*xi(3)*(kappa*be(3) + (1.d0-3.d0*kappa)*fbar_e(3)*fbar_e(3))
     #   +coef5*(kappa*be(3) + (1.d0-3.d0*kappa)*fbar_e(3)*fbar_e(3))
     #   *(kappa*be(3) + (1.d0-3.d0*kappa)*fbar_e(3)*fbar_e(3))
     #  +coef6*(kappa*be(3) + (1.d0-3.d0*kappa)*fbar_e(3)*fbar_e(3))*xi(3)

c     ii=1,jj=1,kk=2,ll=2; iii=1,jjj=2
      ii=1
      jj=1
      kk=2
      ll=2
      iii=1
      jjj=2
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #   +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #   +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #   *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=1,jjj=3; ii=1,jj=1,kk=3,ll=3;
      iii=1
      jjj=3
      ii=1
      jj=1
      kk=3
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #  +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #  *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     # +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=2,jjj=3; ii=2,jj=2,kk=3,ll=3;
      iii=2
      jjj=3
      ii=2
      jj=2
      kk=3
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #  +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #   +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #  *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     # +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=1,jjj=4; ii=1,jj=1,kk=1,ll=2;
      iii=1
      jjj=4
      ii=1
      jj=1
      kk=1
      ll=2
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #   +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #   *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=2,jjj=4; ii=2,jj=2,kk=1,ll=2;
      iii=2
      jjj=4
      ii=2
      jj=2
      kk=1
      ll=2
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #  +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #  *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=3,jjj=4; ii=3,jj=3,kk=1,ll=2;
      iii=3
      jjj=4
      ii=3
      jj=3
      kk=1
      ll=2
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #  +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #  *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     # +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=1,jjj=5; ii=1,jj=1,kk=1,ll=3;
      iii=1
      jjj=5
      ii=1
      jj=1
      kk=1
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #  +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #  *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=2,jjj=5; ii=2,jj=2,kk=1,ll=3;
      iii=2
      jjj=5
      ii=2
      jj=2
      kk=1
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     # +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #   *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     # +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=3,jjj=5; ii=3,jj=3,kk=1,ll=3;
      iii=3
      jjj=5
      ii=3
      jj=3
      kk=1
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #   +coef3*2.d0*xixi(iii,jjj)
     #   +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #    +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #   *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #   +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=4,jjj=5; ii=1,jj=2,kk=1,ll=3;
      iii=4
      jjj=5
      ii=1
      jj=2
      kk=1
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #  +coef3*2.d0*xixi(iii,jjj)
     #   +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #   +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #   *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=4,jjj=4; ii=1,jj=2,kk=1,ll=2;
      iii=4
      jjj=4
      ii=1
      jj=2
      kk=1
      ll=2
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #  +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #   *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=5,jjj=5; ii=1,jj=3,kk=1,ll=3;
      iii=5
      jjj=5
      ii=1
      jj=3
      kk=1
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #   +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #    +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #   *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     # +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=6,jjj=6; ii=2,jj=3,kk=2,ll=3;
      iii=6
      jjj=6
      ii=2
      jj=3
      kk=2
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     #   + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #  +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #   +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #  *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=1,jjj=6; ii=1,jj=1,kk=2,ll=3;
      iii=1
      jjj=6
      ii=1
      jj=1
      kk=2
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     # + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #   +coef3*2.d0*xixi(iii,jjj)
     #   +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #   +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #    *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #   +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=2,jjj=6; ii=2,jj=2,kk=2,ll=3;
      iii=2
      jjj=6
      ii=2
      jj=2
      kk=2
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     #     + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #    +coef3*2.d0*xixi(iii,jjj)
     #   +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #   *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=3,jjj=6; ii=3,jj=3,kk=2,ll=3;
      iii=3
      jjj=6
      ii=3
      jj=3
      kk=2
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     #     + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #  +coef3*2.d0*xixi(iii,jjj)
     #  +coef4*xi(iii)*(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #   +coef5*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #    *(kappa*be(jjj) + (1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #  +coef6*(kappa*be(iii) + (1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=4,jjj=6; ii=1,jj=2,kk=2,ll=3;
      iii=4
      jjj=6
      ii=1
      jj=2
      kk=2
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     #        + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #        +coef3*2.d0*xixi(iii,jjj)
     #     +coef4*xi(iii)*(kappa*be(jjj)+(1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #     +coef5*(kappa*be(iii)+(1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #      *(kappa*be(jjj)+(1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #    +coef6*(kappa*be(iii)+(1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)

c     iii=5,jjj=6; ii=1,jj=3,kk=2,ll=3;
      iii=5
      jjj=6
      ii=1
      jj=3
      kk=2
      ll=3
      ddsedde(iii,jjj) = coef1*xi(iii)*xi(jjj)
     #              + coef2*be(iii)*xi(jjj) + coef2*xi(iii)*be(jjj)
     #              +coef3*2.d0*xixi(iii,jjj)
     #              +coef4*xi(iii)*(kappa*be(jjj)
     #                +(1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #              +coef5*(kappa*be(iii)+(1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))
     #                *(kappa*be(jjj)+(1.d0-3.d0*kappa)*fbar_e(kk)*fbar_e(ll))
     #     +coef6*(kappa*be(iii)+(1.d0-3.d0*kappa)*fbar_e(ii)*fbar_e(jj))*xi(jjj)


c...  use symmetry to fill in the rest
      do ii = 2,ntens
        do jj = 1,ii-1
          ddsedde(ii,jj) = ddsedde(jj,ii)
        end do
      end do

c...  calculate Cauchy stress
c        stress(i) = ((lam*lnJe-mu)*xi(i) + mu*be(i))/detfe
      do ii = 1,ntens
        stress(ii) = (blk/2.d0*(detfe**2.d0-1.d0)
     #                -mu/3.d0*detfe**(-2.d0/3.d0)*I1_el
     #                -2.d0/3.d0*k1*Ebar_el*exp(k2*Ebar_el**2.d0)
     #                *detfe**(-2.d0/3.d0)
     #                *(kappa*I1_el+(1.d0-3.d0*kappa)*I4_el))*xi(ii)
     #               +(mu*detfe**(-2.d0/3.d0)
     #                +2.d0*k1*Ebar_el*exp(k2*Ebar_el**2.d0)*detfe**(-2.d0/3.d0)
     #                *kappa)*be(ii)
     #               +2.d0*k1*Ebar_el*exp(k2*Ebar_el**2.d0)*detfe**(-2.d0/3.d0)
     #                 *(1.d0-3.d0*kappa)*ffbar(ii)
        stress(ii) = stress(ii)/detfe
      end do

c... make stress to matrix foam
      call move6to33(stress,stressmat)

c     Fg^-1:dFg/dtheg
      do ii = 1,3
        do jj = 1,3
          dFgdtheg(ii,jj) = (ximat(ii,jj)-xn0xn0(ii,jj))/2.d0/theg**0.5d0
        end do
      end do

      tmp1 = matmul(matmul(fe,dFgdtheg),matmul(pk2,transpose(dfgrd1)))
      tmp2 = matmul(matmul(dfgrd1,pk2),matmul(dFgdtheg,transpose(fe)))
      tmp3 = matmul(fe,matmul(dFgdtheg,finv))
      call doubledot(fginvmat,dFgdtheg,const1)
      call doubledot42(ddsedde,tmp3,const2)
      call move6to33(const2,const2mat)

ccccc ###################################################
ccccc Similar to the isotropic case, but for anisotropic case we
ccccc need dFgdlam1g and dFgdlam2g, however, not that due to the shape
ccccc of the growth tensor we simply get
ccccc dFgdlam1g = ff0mat
ccccc dFgdlam2g = ss0mat

      tmp1_1 = matmul(matmul(fe,ff0mat),matmul(pk2,transpose(dfgrd1)))
      tmp2_1 = matmul(matmul(dfgrd1,pk2),matmul(ff0mat,transpose(fe)))
      tmp3_1 = matmul(fe,matmul(ff0mat,finv))
      call doubledot(fginvmat,ff0mat,const1_1)
      call doubledot42(ddsedde,tmp3_1,const2_1)
      call move6to33(const2_1,const2mat_1)

      tmp1_2 = matmul(matmul(fe,ss0mat),matmul(pk2,transpose(dfgrd1)))
      tmp2_2 = matmul(matmul(dfgrd1,pk2),matmul(ss0mat,transpose(fe)))
      tmp3_2 = matmul(fe,matmul(ss0mat,finv))
      call doubledot(fginvmat,ff0mat,const1_2)
      call doubledot42(ddsedde,tmp3_2,const2_2)
      call move6to33(const2_2,const2mat_2)
ccccc ###################################################


c...  calculate growth tangent
      do ii =1,3
        do jj = 1,3
          cg_ijmat(ii,jj) = +stressmat(ii,jj)*const1
     #                   -(tmp1(ii,jj) + tmp2(ii,jj))/detf
     #                   -const2mat(ii,jj)
        end do
      end do
      call move33to6(cg_ijmat,cg_ij)

ccccc ###################################################
ccccc Growth tangent for the two directions, like before
ccccc for separate for each of the directions
      do ii =1,3
        do jj = 1,3
          cg_ijmat_1(ii,jj) = +stressmat(ii,jj)*const1_1
     #                   -(tmp1_1(ii,jj) + tmp2_1(ii,jj))/detf
     #                   -const2mat_1(ii,jj)
        end do
      end do
      call move33to6(cg_ijmat_1,cg_ij_1)

      do ii =1,3
        do jj = 1,3
          cg_ijmat_2(ii,jj) = +stressmat(ii,jj)*const1_2
     #                   -(tmp1_2(ii,jj) + tmp2_2(ii,jj))/detf
     #                   -const2mat_2(ii,jj)
        end do
      end do
      call move33to6(cg_ijmat_2,cg_ij_2)
ccccc ###################################################

c...  Second part of the growth tangent, this comes from dthedC
      cg_kl(1) = the*xi(1) - detf*detf/the*ftn0(1)*ftn0(1)
      cg_kl(2) = the*xi(2) - detf*detf/the*ftn0(2)*ftn0(2)
      cg_kl(3) = the*xi(3) - detf*detf/the*ftn0(3)*ftn0(3)
      cg_kl(4) = the*xi(4) - detf*detf/the*ftn0(1)*ftn0(2)
      cg_kl(5) = the*xi(5) - detf*detf/the*ftn0(1)*ftn0(3)
      cg_kl(6) = the*xi(6) - detf*detf/the*ftn0(2)*ftn0(3)

ccccc ###################################################
ccccc We need different cg_kl for lam1 and lam 2
ccccc CHECK cg_kl wrt lam1 -> dlam1/dC
ccccc CHECK cg_kl wrt lam2 -> dlam2/dC
ccccc TL: I think it is correct
      cg_kl_1(1) = 1.d0/2.d0/lam1*ff0(1)
      cg_kl_1(2) = 1.d0/2.d0/lam1*ff0(2)
      cg_kl_1(3) = 1.d0/2.d0/lam1*ff0(3)
      cg_kl_1(4) = 1.d0/2.d0/lam1*ff0(4)
      cg_kl_1(5) = 1.d0/2.d0/lam1*ff0(5)
      cg_kl_1(6) = 1.d0/2.d0/lam1*ff0(6)

      cg_kl_2(1) = 1.d0/2.d0/lam2*ss0(1)
      cg_kl_2(2) = 1.d0/2.d0/lam2*ss0(2)
      cg_kl_2(3) = 1.d0/2.d0/lam2*ss0(3)
      cg_kl_2(4) = 1.d0/2.d0/lam2*ss0(4)
      cg_kl_2(5) = 1.d0/2.d0/lam2*ss0(5)
      cg_kl_2(6) = 1.d0/2.d0/lam2*ss0(6)
ccccc ###################################################

c...  geometric tangent
      do ii=1,3
        ddsgdde(ii,ii) = 2.d0*stress(ii)
      end do
      ddsgdde(1,2) = 0.d0
      ddsgdde(1,3) = 0.d0
      ddsgdde(2,3) = 0.d0

      ddsgdde(1,4) = stress(4)
      ddsgdde(1,5) = stress(5)
      ddsgdde(1,6) = 0.d0
      ddsgdde(2,4) = stress(4)
      ddsgdde(2,5) = 0.d0
      ddsgdde(2,6) = stress(6)
      ddsgdde(3,4) = 0.d0
      ddsgdde(3,5) = stress(5)
      ddsgdde(3,6) = stress(6)
      ddsgdde(4,5) = stress(6)/2.d0
      ddsgdde(4,6) = stress(5)/2.d0
      ddsgdde(5,6) = stress(4)/2.d0
      ddsgdde(4,4) = (stress(2) + stress(1))/2.d0
      ddsgdde(5,5) = (stress(3) + stress(1))/2.d0
      ddsgdde(6,6) = (stress(3) + stress(2))/2.d0

c... use symmetric condition to fill other parts of geometric tangent
      do ii = 2,ntens
        do jj = 1,ii-1
          ddsgdde(ii,jj) = ddsgdde(jj,ii)
        end do
      end do

ccccc ###################################################
ccccc CHANGING here for the anisotropic growth, adding two
ccccc contributions for the tangent instead of one
c...  compile tangent
      do ii = 1,ntens
        do jj = 1,ntens
          ddsdde(ii,jj)=ddsedde(ii,jj)+ddsgdde(ii,jj)
     #                  +fac1*cg_ij_1(ii)*cg_kl_1(jj)
     #                  +fac2*cg_ij_2(ii)*cg_kl_2(jj)
        end do
      end do
ccccc ###################################################
c      print *, 'ddsdde(1,1),ddsdde(2,2),ddsdde(3,3)',ddsdde(1,1),ddsdde(2,2),ddsdde(3,3)
c...  calculate strain energy
      sse = +mu/2.d0*(I1bar_el-3.d0)
     #      +blk/2.d0*((detfe**2.d0-1.d0)/2.d0-lnJe)
     #      +0.5d0*k1/k2*(exp(k2*Ebar_el**2.d0)-1.d0)
      return
      end

c...  ------------------------------------------------------------------
c...  Subroutine for doubledot product of two 2nd order tensor
c...  ------------------------------------------------------------------
      subroutine doubledot(matA,matB,prd1)
      implicit none
      real*8   matA(3,3), matB(3,3)
      real*8   prd1

      prd1=+matA(1,1)*matB(1,1)+matA(2,2)*matB(2,2)+matA(3,3)*matB(3,3)
     #       +matA(1,2)*matB(1,2)+matA(1,3)*matB(1,3)+matA(2,3)*matB(2,3)
     #       +matA(2,1)*matB(2,1)+matA(3,1)*matB(3,1)+matA(3,2)*matB(3,2)

      return
      end

c...  ------------------------------------------------------------------
c...  Subroutine for doubledot product of symmetric 4th order and 2nd order tensor
c...  ------------------------------------------------------------------
      subroutine doubledot42(matA,matB,voigtC)
      implicit none
      real*8   matA(6,6), matB(3,3), voigtC(6)
      integer  ii

      do ii = 1,6
        voigtC(ii)=+matA(ii,1)*matB(1,1)+matA(ii,2)
     #             *matB(2,2)+matA(ii,3)*matB(3,3)
     #             +matA(ii,4)*matB(1,2)+matA(ii,5)
     #             *matB(1,3)+matA(ii,6)*matB(2,3)
     #             +matA(ii,4)*matB(2,1)+matA(ii,5)
     #             *matB(3,1)+matA(ii,6)*matB(3,2)
      end do

      return
      end

c...  ------------------------------------------------------------------
c...  Subroutine for moving 6x1 tensor to 3x3 tensor
c...  ------------------------------------------------------------------
      subroutine move6to33(a6,a33)
      implicit none
      real*8   a6(6),a33(3,3)

      a33(1,1) = a6(1)
      a33(2,2) = a6(2)
      a33(3,3) = a6(3)
      a33(1,2) = a6(4)
      a33(1,3) = a6(5)
      a33(2,3) = a6(6)
      a33(2,1) = a33(1,2)
      a33(3,1) = a33(1,3)
      a33(3,2) = a33(2,3)

      return
      end

c...  ------------------------------------------------------------------
c...  Subroutine for moving 3x3 tensor to 6x1 tensor
c...  ------------------------------------------------------------------
      subroutine move33to6(a33,a6)
      implicit none
      real*8   a6(6),a33(3,3)

      a6(1) = a33(1,1)
      a6(2) = a33(2,2)
      a6(3) = a33(3,3)
      a6(4) = a33(1,2)
      a6(5) = a33(1,3)
      a6(6) = a33(2,3)

      return
      end
c...  ------------------------------------------------------------------
