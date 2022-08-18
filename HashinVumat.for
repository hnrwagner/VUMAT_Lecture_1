      subroutine vumat( 
c Read only - 
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, 
     2  stepTime, totalTime, dt, cmname, coordMp, charLength, 
     3  props, density, strainInc, relSpinInc, 
     4  tempOld, stretchOld, defgradOld, fieldOld, 
     5  stressOld, stateOld, enerInternOld, enerInelasOld, 
     6  tempNew, stretchNew, defgradNew, fieldNew, 
c Write only - 
     7  stressNew, stateNew, enerInternNew, enerInelasNew ) 
c 
      include 'vaba_param.inc' 
c 
c 3D Orthotropic Elasticity with Hashin 3d Failure criterion 
c 
c The state variables are stored as: 
c    state(*,1)   = material point status 
c    state(*,2:7) = damping stresses 
c 
c User defined material properties are stored as 
c  * First line: 
c     props(1) --> Young's modulus in 1-direction, E1 
c     props(2) --> Young's modulus in 2-direction, E2 
c     props(3) --> Young's modulus in 3-direction, E3 
c     props(4) --> Poisson's ratio, nu12 
c     props(5) --> Poisson's ratio, nu13 
c     props(6) --> Poisson's ratio, nu23 
c     props(7) --> Shear modulus, G12 
c     props(8) --> Shear modulus, G13 
c 
c  * Second line: 
c     props(9)  --> Shear modulus, G23 
c     props(10) --> beta damping parameter 
c     props(11) --> "not used" 
c     props(12) --> "not used" 
c     props(13) --> "not used" 
c     props(14) --> "not used" 
c     props(15) --> "not used" 
c     props(16) --> "not used" 
c 
c  * Third line: 
c     props(17) --> Ultimate tens stress in 1-direction, sigu1t 
c     props(18) --> Ultimate comp stress in 1-direction, sigu1c 
c     props(19) --> Ultimate tens stress in 2-direction, sigu2t 
c     props(20) --> Ultimate comp stress in 2-direction, sigu2c 
c     props(21) --> Ultimate tens stress in 2-direction, sigu3t 
c     props(22) --> Ultimate comp stress in 2-direction, sigu3c 
c     props(23) --> "not used" 
c     props(24) --> "not used" 
c 
c  * Fourth line: 
c     props(25) --> Ultimate shear stress, sigu12 
c     props(26) --> Ultimate shear stress, sigu13 
c     props(27) --> Ultimate shear stress, sigu23 
c     props(28) --> "not used" 
c     props(29) --> "not used" 
c     props(30) --> "not used" 
c     props(31) --> "not used" 
c     props(32) --> "not used" 
c 

      dimension props(nprops), density(nblock), 
     1  coordMp(nblock,*), 
     2  charLength(*), strainInc(nblock,ndir+nshr), 
     3  relSpinInc(nblock,nshr), tempOld(nblock), 
     4  stretchOld(nblock,ndir+nshr), defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr), 
     6  stateOld(nblock,nstatev), enerInternOld(nblock), 
     7  enerInelasOld(nblock), tempNew(*), 
     8  stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv), stressNew(nblock,ndir+nshr), 
     1  stateNew(nblock,nstatev), 
     2  enerInternNew(nblock), enerInelasNew(nblock) 
* 
      character*80 cmname 
* 
      parameter( zero = 0.d0, one = 1.d0, two = 2.d0, half = .5d0 ) 
* 
      parameter( 
     *     i_svd_DmgFiberT   = 1, 
     *     i_svd_DmgFiberC   = 2, 
     *     i_svd_DmgMatrixT  = 3, 
     *     i_svd_DmgMatrixC  = 4, 
     *     i_svd_statusMp   = 5, 
     *     i_svd_dampStress = 6, 
c     *    i_svd_dampStressXx = 6, 
c     *    i_svd_dampStressYy = 7, 
c     *    i_svd_dampStressZz = 8, 
c     *    i_svd_dampStressXy = 9, 
c     *    i_svd_dampStressYz = 10, 
c     *    i_svd_dampStressZx = 11, 
     *     i_svd_Strain   = 12, 
c     *    i_svd_StrainXx = 12, 
c     *    i_svd_StrainYy = 13, 
c     *    i_svd_StrainZz = 14, 
c     *    i_svd_StrainXy = 15, 
c     *    i_svd_StrainYz = 16, 
c     *    i_svd_StrainZx = 17, 
     *     n_svd_required = 17 ) 
* 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6 ) 
* 
* Structure of property array 
      parameter ( 
     *     i_pro_E1    = 1, 
     *     i_pro_E2    = 2, 
     *     i_pro_E3    = 3, 
     *     i_pro_nu12  = 4, 
     *     i_pro_nu13  = 5, 
     *     i_pro_nu23  = 6, 
     *     i_pro_G12   = 7, 
     *     i_pro_G13   = 8, 
     *     i_pro_G23   = 9, 
* 
     *     i_pro_beta  = 19, 
* 
     *     i_pro_sigu1t = 10, 
     *     i_pro_sigu1c = 11, 
     *     i_pro_sigu2t = 12, 
     *     i_pro_sigu2c = 13, 
     *     i_pro_sigu3t = 14, 
     *     i_pro_sigu3c = 15, 
     *     i_pro_sigu12 = 16, 
     *     i_pro_sigu13 = 17, 
     *     i_pro_sigu23 = 18 ) 
* Temporary arrays 
      dimension eigen(maxblk*3) 
* 
* Read material properties 
* 
      E1 = props(i_pro_E1) 
      E2 = props(i_pro_E2) 
      E3 = props(i_pro_E3) 
      xnu12 = props(i_pro_nu12) 
      xnu13 = props(i_pro_nu13) 
      xnu23 = props(i_pro_nu23) 
      G12 = props(i_pro_G12) 
      G13 = props(i_pro_G13) 
      G23 = props(i_pro_G23) 
* 
      xnu21 = xnu12 * E2 / E1 
      xnu31 = xnu13 * E3 / E1 
      xnu32 = xnu23 * E3 / E2 
* 
* 
* Compute terms of stiffness matrix 
      gg = one / ( one - xnu12*xnu21 - xnu23*xnu32 - xnu31*xnu13 
     *     - two*xnu21*xnu32*xnu13 ) 
      C11  = E1 * ( one - xnu23*xnu32 ) * gg 
      C22  = E2 * ( one - xnu13*xnu31 ) * gg 
      C33  = E3 * ( one - xnu12*xnu21 ) * gg 
      C12  = E1 * ( xnu21 + xnu31*xnu23 ) * gg 
      C13  = E1 * ( xnu31 + xnu21*xnu32 ) * gg 
      C23  = E2 * ( xnu32 + xnu12*xnu31 ) * gg 
* 
      f1t = props(i_pro_sigu1t) 
      f1c = props(i_pro_sigu1c) 
      f2t = props(i_pro_sigu2t) 
      f2c = props(i_pro_sigu2c) 
      f3t = props(i_pro_sigu3t) 
      f3c = props(i_pro_sigu3c) 
      f12 = props(i_pro_sigu12) 
      f13 = props(i_pro_sigu13) 
      f23 = props(i_pro_sigu23) 
* 
      beta = props(i_pro_beta) 
* 
* Assume purely elastic material at the beginning of the analysis 
*       
      if ( totalTime .eq. zero ) then 
         if (nstatev .lt. n_svd_Required) then 
            call xplb_abqerr(-2,'Subroutine VUMAT requires the '// 
     *           'specification of %I state variables. Check the '// 
     *           'definition of *DEPVAR in the input file.', 
     *           n_svd_Required,zero,' ') 
            call xplb_exit 
         end if 
         call OrthoEla3dExp ( nblock, 
     *        stateOld(1,i_svd_DmgFiberT), 
     *        stateOld(1,i_svd_DmgFiberC), 
     *        stateOld(1,i_svd_DmgMatrixT), 
     *        stateOld(1,i_svd_DmgMatrixC), 
     *        C11, C22, C33, C12, C23, C13, G12, G23, G13, 
     *        strainInc, 
     *        stressNew ) 
         return 
      end if 
* 
*  Update total elastic strain 
      call strainUpdate ( nblock, strainInc, 
     *     stateOld(1,i_svd_strain), stateNew(1,i_svd_strain) ) 
* 
* Stress update 
      call OrthoEla3dExp ( nblock, 
     *     stateOld(1,i_svd_DmgFiberT), 
     *     stateOld(1,i_svd_DmgFiberC), 
     *     stateOld(1,i_svd_DmgMatrixT), 
     *     stateOld(1,i_svd_DmgMatrixC), 
     *     C11, C22, C33, C12, C23, C13, G12, G23, G13, 
     *     stateNew(1,i_svd_strain), 
     *     stressNew ) 
* 
* Failure evaluation 
* 
      call copyr ( nblock, 
     *     stateOld(1,i_svd_DmgFiberT), stateNew(1,i_svd_DmgFiberT) ) 
      call copyr ( nblock, 
     *     stateOld(1,i_svd_DmgFiberC), stateNew(1,i_svd_DmgFiberC) ) 
      call copyr ( nblock, 
     *     stateOld(1,i_svd_DmgMatrixT), stateNew(1,i_svd_DmgMatrixT) ) 
      call copyr ( nblock, 
     *     stateOld(1,i_svd_DmgMatrixC), stateNew(1,i_svd_DmgMatrixC) ) 
      nDmg = 0 
      call eig33Anal ( nblock, stretchNew, eigen ) 
      call Hashin3d  ( nblock, nDmg, 
     *     f1t, f2t, f3t, f1c, f2c, f3c, f12, f23, f13, 
     *     stateNew(1,i_svd_DmgFiberT), 
     *     stateNew(1,i_svd_DmgFiberC), 
     *     stateNew(1,i_svd_DmgMatrixT), 
     *     stateNew(1,i_svd_DmgMatrixC), 
     *     stateNew(1,i_svd_statusMp), 
     *     stressNew, eigen ) 
*     -- Recompute stresses if new Damage is occurring 
      if ( nDmg .gt. 0 ) then 
         call OrthoEla3dExp ( nblock, 
     *        stateNew(1,i_svd_DmgFiberT), 
     *        stateNew(1,i_svd_DmgFiberC), 
     *        stateNew(1,i_svd_DmgMatrixT), 
     *        stateNew(1,i_svd_DmgMatrixC), 
     *        C11, C22, C33, C12, C23, C13, G12, G23, G13, 
     *        stateNew(1,i_svd_strain), 
     *        stressNew ) 
      end if 
* 
* Beta damping 
      if ( beta .gt. zero ) then 
         call betaDamping3d ( nblock, 
     *        beta, dt, strainInc, 
     *        stressOld, stressNew, 
     *        stateNew(1,i_svd_statusMp), 
     *        stateOld(1,i_svd_dampStress), 
     *        stateNew(1,i_svd_dampStress) ) 
      end if 
* 
* Integrate the internal specific energy (per unit mass) 
* 
      call EnergyInternal3d ( nblock, stressOld, stressNew, 
     *   strainInc, density, enerInternOld, enerInternNew ) 
* 
      return 
      end 


************************************************************ 
*   OrthoEla3dExp: Orthotropic elasticity - 3d             * 
************************************************************ 
      subroutine OrthoEla3dExp ( nblock, 
     *     dmgFiberT, dmgFiberC, dmgMatrixT, dmgMatrixC, 
     *     C11, C22, C33, C12, C23, C13, G12, G23, G13, 
     *     strain, stress ) 
* 
      include 'vaba_param.inc' 

*  Orthotropic elasticity, 3D case - 
* 
      parameter( zero = 0.d0, one = 1.d0, two = 2.d0) 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      dimension  strain(nblock,n_s33_Car), 
     *     dmgFiberT(nblock), dmgFiberC(nblock), 
     *     dmgMatrixT(nblock), dmgMatrixC(nblock), 
     *     stress(nblock,n_s33_Car) 
*     -- shear fraction in matrix tension and compression mode 
      parameter ( smt = 0.9d0, smc = 0.5d0 ) 
* 
      do k = 1, nblock 
*     -- Compute damaged stiffness 
         dft = dmgFiberT(k) 
         dfc = dmgFiberC(k) 
         dmt = dmgMatrixT(k) 
         dmc = dmgMatrixC(k) 
         df = one - ( one - dft ) * ( one - dfc ) 
* 
         dC11 = ( one - df ) * C11 
         dC22 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C22 
         dC33 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C33 
         dC12 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C12 
         dC23 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C23 
         dC13 = ( one - df ) * ( one - dmt ) * ( one - dmc ) * C13 
         dG12 = ( one - df ) 
     *        * ( one - smt*dmt ) * ( one - smc*dmc ) * G12 
         dG23 = ( one - df ) 
     *        * ( one - smt*dmt ) * ( one - smc*dmc ) * G23 
         dG13 = ( one - df ) 
     *        * ( one - smt*dmt ) * ( one - smc*dmc ) * G13 
*     -- Stress update 
         stress(k,i_s33_Xx) = dC11 * strain(k,i_s33_Xx) 
     *        + dC12 * strain(k,i_s33_Yy) 
     *        + dC13 * strain(k,i_s33_Zz) 
         stress(k,i_s33_Yy) = dC12 * strain(k,i_s33_Xx) 
     *        + dC22 * strain(k,i_s33_Yy) 
     *        + dC23 * strain(k,i_s33_Zz) 
         stress(k,i_s33_Zz) = dC13 * strain(k,i_s33_Xx) 
     *        + dC23 * strain(k,i_s33_Yy) 
     *        + dC33 * strain(k,i_s33_Zz) 
         stress(k,i_s33_Xy) = two * dG12 * strain(k,i_s33_Xy) 
         stress(k,i_s33_Yz) = two * dG23 * strain(k,i_s33_Yz) 
         stress(k,i_s33_Zx) = two * dG13 * strain(k,i_s33_Zx) 
      end do 
*     
      return 
      end 

************************************************************ 
*   strainUpdate: Update total strain                      * 
************************************************************ 
      subroutine strainUpdate ( nblock, 
     *     strainInc, strainOld, strainNew ) 
* 
      include 'vaba_param.inc' 
* 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      dimension strainInc(nblock,n_s33_Car), 
     *     strainOld(nblock,n_s33_Car), 
     *     strainNew(nblock,n_s33_Car) 
* 
      do k = 1, nblock 
         strainNew(k,i_s33_Xx)= strainOld(k,i_s33_Xx) 
     *                        + strainInc(k,i_s33_Xx) 
         strainNew(k,i_s33_Yy)= strainOld(k,i_s33_Yy) 
     *                        + strainInc(k,i_s33_Yy) 
         strainNew(k,i_s33_Zz)= strainOld(k,i_s33_Zz) 
     *                        + strainInc(k,i_s33_Zz) 
         strainNew(k,i_s33_Xy)= strainOld(k,i_s33_Xy) 
     *                        + strainInc(k,i_s33_Xy) 
         strainNew(k,i_s33_Yz)= strainOld(k,i_s33_Yz) 
     *                        + strainInc(k,i_s33_Yz) 
         strainNew(k,i_s33_Zx)= strainOld(k,i_s33_Zx) 
     *                        + strainInc(k,i_s33_Zx) 
      end do 
* 
      return 
      end 


************************************************************ 
*   Hashin3d w/ Modified Puck: Evaluate Hashin 3d failure  * 
*   criterion for fiber, Puck for matrix                   * 
************************************************************ 
      subroutine Hashin3d ( nblock, nDmg, 
     *     f1t, f2t, f3t, f1c, f2c, f3c, f12, f23, f13, 
     *     dmgFiberT, dmgFiberC, dmgMatrixT, dmgMatrixC, 
     *     statusMp, stress, eigen ) 
* 
      include 'vaba_param.inc' 

      parameter( zero = 0.d0, one = 1.d0, half = 0.5d0, three =3.d0 ) 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      parameter(i_v3d_X=1,i_v3d_Y=2,i_v3d_Z=3 ) 
      parameter(n_v3d_Car=3 ) 
* 
      parameter ( eMax = 1.00d0, eMin = -0.8d0 ) 
* 
      dimension  dmgFiberT(nblock), dmgFiberC(nblock), 
     *     dmgMatrixT(nblock), dmgMatrixC(nblock), 
     *     stress(nblock,n_s33_Car), 
     *     eigen(nblock,n_v3d_Car), 
     *     statusMp(nblock) 
* 
      f1tInv = zero 
      f2tInv = zero 
      f3tInv = zero 
      f1cInv = zero 
      f2cInv = zero 
      f3cInv = zero 
      f12Inv = zero 
      f23Inv = zero 
      f13Inv = zero 
* 
      if ( f1t .gt. zero ) f1tInv = one / f1t 
      if ( f2t .gt. zero ) f2tInv = one / f2t 
      if ( f3t .gt. zero ) f3tInv = one / f3t 
      if ( f1c .gt. zero ) f1cInv = one / f1c 
      if ( f2c .gt. zero ) f2cInv = one / f2c 
      if ( f3c .gt. zero ) f3cInv = one / f3c 
      if ( f12 .gt. zero ) f12Inv = one / f12 
      if ( f23 .gt. zero ) f23Inv = one / f23 
      if ( f13 .gt. zero ) f13Inv = one / f13 
* 
      do k = 1, nblock 
         if ( statusMp(k) .eq. one ) then 
*     
         lFail = 0 
* 
         s11 = stress(k,i_s33_Xx) 
         s22 = stress(k,i_s33_Yy) 
         s33 = stress(k,i_s33_Zz) 
         s12 = stress(k,i_s33_Xy) 
         s23 = stress(k,i_s33_Yz) 
         s13 = stress(k,i_s33_Zx) 
* 
*     Evaluate Fiber modes 
         if ( s11 .gt. zero ) then 
*     -- Tensile Fiber Mode 
            rft = (s11*f1tInv )**2 + (s12*f12Inv )**2 + (s13*f13Inv )**2
            if ( rft .ge. one ) then 
               lDmg = 1 
               dmgFiberT(k) = one 
            end if 
         else if ( s11 .lt. zero ) then 
*     -- Compressive Fiber Mode 
            rfc = abs(s11) * f1cInv 
            if ( rfc .ge. one ) then 
               lDmg = 1 
               dmgFiberC(k) = one 
            end if 
         end if 
* 
*     Evaluate Matrix Modes 
         if ( ( s22 + s33 ) .gt. zero ) then 
*     -- Tensile Matrix mode 
            rmt = ( s11 * half * f1tInv )**2 
     *           + ( s22**2 * abs(f2tInv * f2cInv) ) 
     *           + ( s12 * f12Inv )**2 
     *           + ( s22 * (f2tInv + f2cInv) ) 
            if ( rmt .ge. one ) then 
               lDmg = 1 
               dmgMatrixT(k) = one 
            end if 
         else if ( ( s22 + s33 ) .lt. zero ) then 
*     -- Compressive Matrix Mode 
            rmc = ( s11 * half * f1tInv )**2 
     *           + ( s22**2 * abs(f2tInv * f2cInv) ) 
     *           + ( s12 * f12Inv )**2 
     *           + ( s22 * (f2tInv + f2cInv) ) 
            if ( rmc .ge. one ) then 
               lDmg = 1 
               dmgMatrixC(k) = one 
            end if 
         end if 
* 
         eigMax=max(eigen(k,i_v3d_X),eigen(k,i_v3d_Y),eigen(k,i_v3d_Z)) 
         eigMin=min(eigen(k,i_v3d_X),eigen(k,i_v3d_Y),eigen(k,i_v3d_Z)) 
         enomMax = eigMax - one 
         enomMin = eigMin - one 
* 
         if ( enomMax .gt. eMax .or. 
     *        enomMin .lt. eMin .or. 
     *        dmgFiberT(k) .eq. one ) then 
            statusMp(k) = zero 
         end if 
* 
         nDmg = nDmk + lDmg 
* 
         end if 
* 
      end do 
* 
      return 
      end 


************************************************************ 
*   betaDamping: Add beta damping                          * 
************************************************************ 
      subroutine betaDamping3d ( nblock, 
     *     beta, dt, strainInc, sigOld, sigNew, 
     *     statusMp, sigDampOld, sigDampNew ) 
* 
      include 'vaba_param.inc' 
* 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      dimension  sigOld(nblock,n_s33_Car), 
     *     sigNew(nblock,n_s33_Car), 
     *     strainInc(nblock,n_s33_Car), 
     *     statusMp(nblock), 
     *     sigDampOld(nblock,n_s33_Car), 
     *     sigDampNew(nblock,n_s33_Car)       
* 
      parameter ( zero = 0.d0, one = 1.d0, two=2.0d0, 
     *     half = 0.5d0, third = 1.d0/3.d0 ) 
      parameter ( asmall = 1.d-16 ) 
*     
      betaddt =  beta / dt 
* 
      do k =1 , nblock 
         sigDampNew(k,i_s33_Xx) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Xx) 
     *        - ( sigOld(k,i_s33_Xx) - sigDampOld(k,i_s33_Xx) ) ) 
         sigDampNew(k,i_s33_Yy) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Yy) 
     *        - ( sigOld(k,i_s33_Yy) - sigDampOld(k,i_s33_Yy) ) ) 
         sigDampNew(k,i_s33_Zz) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Zz) 
     *        - ( sigOld(k,i_s33_Zz) - sigDampOld(k,i_s33_Zz) ) ) 
         sigDampNew(k,i_s33_Xy) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Xy) 
     *        - ( sigOld(k,i_s33_Xy) - sigDampOld(k,i_s33_Xy) ) ) 
         sigDampNew(k,i_s33_Yz) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Yz) 
     *        - ( sigOld(k,i_s33_Yz) - sigDampOld(k,i_s33_Yz) ) ) 
         sigDampNew(k,i_s33_Zx) = betaddt * statusMp(k) * 
     *        ( sigNew(k,i_s33_Zx) 
     *        - ( sigOld(k,i_s33_Zx) - sigDampOld(k,i_s33_Zx) ) ) 
* 
         sigNew(k,i_s33_Xx) = sigNew(k,i_s33_Xx)+sigDampNew(k,i_s33_Xx) 
         sigNew(k,i_s33_Yy) = sigNew(k,i_s33_Yy)+sigDampNew(k,i_s33_Yy) 
         sigNew(k,i_s33_Zz) = sigNew(k,i_s33_Zz)+sigDampNew(k,i_s33_Zz) 
         sigNew(k,i_s33_Xy) = sigNew(k,i_s33_Xy)+sigDampNew(k,i_s33_Xy) 
         sigNew(k,i_s33_Yz) = sigNew(k,i_s33_Yz)+sigDampNew(k,i_s33_Yz) 
         sigNew(k,i_s33_Zx) = sigNew(k,i_s33_Zx)+sigDampNew(k,i_s33_Zx) 
* 
      end do 
*     
      return 
      end 


************************************************************ 
*   EnergyInternal3d: Compute internal energy for 3d case  * 
************************************************************ 
      subroutine EnergyInternal3d(nblock, sigOld, sigNew , 
     *   strainInc, curDensity, enerInternOld, enerInternNew) 
* 
      include 'vaba_param.inc' 
* 
      parameter( 
     *     i_s33_Xx = 1, 
     *     i_s33_Yy = 2, 
     *     i_s33_Zz = 3, 
     *     i_s33_Xy = 4, 
     *     i_s33_Yz = 5, 
     *     i_s33_Zx = 6, 
     *     n_s33_Car = 6 ) 
* 
      parameter( two = 2.d0, half = .5d0 ) 
* 
      dimension sigOld (nblock,n_s33_Car), sigNew (nblock,n_s33_Car), 
     *     strainInc (nblock,n_s33_Car), curDensity (nblock), 
     *     enerInternOld(nblock), enerInternNew(nblock) 
* 
      do k = 1, nblock 
         stressPower  = half * ( 
     *        ( sigOld(k,i_s33_Xx) + sigNew(k,i_s33_Xx) ) 
     *        * ( strainInc(k,i_s33_Xx) ) 
     *        +       ( sigOld(k,i_s33_Yy) + sigNew(k,i_s33_Yy) ) 
     *        * ( strainInc(k,i_s33_Yy)) 
     *        +       ( sigOld(k,i_s33_Zz) + sigNew(k,i_s33_Zz) ) 
     *        * ( strainInc(k,i_s33_Zz)) 
     *        + two * ( sigOld(k,i_s33_Xy) + sigNew(k,i_s33_Xy) ) 
     *        * strainInc(k,i_s33_Xy) 
     *        + two * ( sigOld(k,i_s33_Yz) + sigNew(k,i_s33_Yz) ) 
     *        * strainInc(k,i_s33_Yz) 
     *        + two * ( sigOld(k,i_s33_Zx) + sigNew(k,i_s33_Zx) ) 
     *        * strainInc(k,i_s33_Zx) ) 
*     
         enerInternNew(k) = enerInternOld(k) + stressPower/curDensity(k)
      end do 
*     
      return   
      end   

************************************************************ 
*   CopyR: Copy from one array to another                  * 
************************************************************ 
      subroutine CopyR(nCopy, from, to ) 
* 
      include 'vaba_param.inc' 
* 
      dimension from(nCopy), to(nCopy) 
* 
      do k = 1, nCopy 
         to(k) = from(k) 
      end do 
* 
      return 
      end 

********************************************************************* 
***** 
* eig33Anal: Compute eigen values of a 3x3 symmetric matrix analytically * 
********************************************************************* 
***** 
      subroutine eig33Anal( nblock, sMat, eigVal ) 
* 
      include 'vaba_param.inc' 
* 
      parameter(i_s33_Xx=1,i_s33_Yy=2,i_s33_Zz=3 ) 
      parameter(i_s33_Xy=4,i_s33_Yz=5,i_s33_Zx=6 ) 
      parameter(i_s33_Yx=i_s33_Xy ) 
      parameter(i_s33_Zy=i_s33_Yz ) 
      parameter(i_s33_Xz=i_s33_Zx,n_s33_Car=6 ) 
* 
      parameter(i_v3d_X=1,i_v3d_Y=2,i_v3d_Z=3 ) 
      parameter(n_v3d_Car=3 ) 
* 
      parameter ( zero = 0.d0, one = 1.d0, two = 2.d0, 
     *     three = 3.d0, half = 0.5d0, third = one / three, 
     *     pi23 = 2.094395102393195d0, 
     *     fuzz = 1.d-8, 
     *     preciz = fuzz * 1.d4 ) 
* 
      dimension eigVal(nblock,n_v3d_Car), sMat(nblock,n_s33_Car) 
* 
      do k = 1, nblock 
        sh  = third*(sMat(k,i_s33_Xx)+sMat(k,i_s33_Yy)+sMat(k,i_s33_Zz))
        s11 = sMat(k,i_s33_Xx) - sh 
        s22 = sMat(k,i_s33_Yy) - sh 
        s33 = sMat(k,i_s33_Zz) - sh 
        s12 = sMat(k,i_s33_Xy) 
        s13 = sMat(k,i_s33_Xz) 
        s23 = sMat(k,i_s33_Yz) 
* 
        fac  = max(abs(s11), abs(s22), abs(s33)) 
        facs = max(abs(s12), abs(s13), abs(s23)) 
        if( facs .lt. (preciz*fac) ) then 
          eigVal(k,i_v3d_X) = sMat(k,i_s33_Xx) 
          eigVal(k,i_v3d_Y) = sMat(k,i_s33_Yy) 
          eigVal(k,i_v3d_Z) = sMat(k,i_s33_Zz) 
        else 
          q = third*((s12**2+s13**2+s23**2)+half*(s11**2+s22**2+s33**2))
          fac = two * sqrt(q) 
          if( fac .gt. fuzz ) then 
            ofac = two/fac 
          else 
            ofac = zero 
          end if 
          s11 = ofac*s11 
          s22 = ofac*s22 
          s33 = ofac*s33 
          s12 = ofac*s12 
          s13 = ofac*s13 
          s23 = ofac*s23 
          r = s12*s13*s23 
     *         + half*(s11*s22*s33-s11*s23**2-s22*s13**2-s33*s12**2) 
          if( r .ge. one-fuzz ) then 
            cos1 = -half 
            cos2 = -half 
            cos3 = one 
          else if( r .le. fuzz-one ) then 
            cos1 = -one 
            cos2 = half 
            cos3 = half 
          else 
            ang = third * acos(r) 
            cos1 = cos(ang) 
            cos2 = cos(ang+pi23) 
            cos3 =-cos1-cos2 
          end if 
          eigVal(k,i_v3d_X) = sh + fac*cos1 
          eigVal(k,i_v3d_Y) = sh + fac*cos2 
          eigVal(k,i_v3d_Z) = sh + fac*cos3 
        end if 
      end do 
* 
      return 
      end