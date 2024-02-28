  	  subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C      
      integer readflag
      real*4  edgedim(350000,3)
      common /integerbuf/readflag 
      common / realbuf/edgedim

        
        IF ( cmname(1:8).eq.'MATFIBRE' ) THEN
C		
            CALL  vumatfibre( 
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
	 
		ELSEIF ( cmname(1:6).eq.'MATCOH' ) THEN
 
            CALL vumatcoh(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
     
		ENDIF
C         
      RETURN
      END       

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C  subroutine: vumatfibre                                            C
C  function: Fibre failure activation using maximum stress criteria  C
C            and evolution                                           C 
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C SDV1  - Strain 11
C SDV2  - Strain 22
C SDV3  - Strain 33
C SDV4  - Strain 12, engineering shear strain gamma
C SDV5  - Temporary strain 12
C SDV6  - Element deformation gradient, F
C SDV7  - Failure index under tension, FIt
C SDV8  - Failure index under compression, FIc
C SDV9  - Damage variable, D
C SDV10 - Element Delete flag; 0-delete
C SDV11 - Stress 11 at fibre damage initiation under tension, sig_f0_T     
C SDV12 - Strain 11 at fibre damage initiation under tension, eps_f0_T
C SDV13 - Strain 11 at fibre damage completion under tension, eps_ff_T
C SDV14 - Characteristic length, L
C SDV15 - Fibre failure falg: 0-elastic, 1-softening
C SDV16 - Old d1Plus
C SDV17 - Old D1Minus
C SDV18 - Stress 11 at fibre damage initiation under compression, sig_f0_C
C SDV19 - Strain 11 at fibre damage initiation under compression, eps_f0_C
C SDV20 - Strain 11 at fibre damage completion under compression, eps_ff_C
C SDV21 - Fibre failure flag under tension: 0-elastic, 1-softening
C SDV22 - Fibre failure flag under compression: 0-elastic, 1-softening


#include 'matrixUtil.for'
#include 'stress.for'

		 subroutine vumatfibre(
! Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
! Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
	  
	  Use matrixAlgUtil_Mod
	  Use stress_Mod
      
	  include 'vaba_param.inc'
	  
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
	 
       parameter(zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0, 
     * third = one/three, half = 0.5d0, twothird = two/three,
     * threehalf = 1.5d0, safety = 1.d-20, tolSgn = 1.d-20, S_res=50.d0)              
     
       character*80 cmname,cpname
       character*256 outdir,fullpath
       integer intnum,locnum,jrcd,lenoutdir

	!  Material properties
       Double Precision :: E11,E22,E33,v12,v21,v13,v31, 
     *				       v23,v32,G12,G13,G23
	   Double Precision :: Xt,Xc,Yc,Yt,Sl
	   Double Precision :: Gft,Gfc

	!  Stiffness matrix
	   Double Precision :: n,m,delta
	   Double Precision :: C(3+nshr, 3+nshr)

	!  Nonlinear shear parameters 
       Double Precision :: sgn4,A,B
	   
	!  DGD related parameters
	   Double Precision :: F(3,3)
	   Double Precision :: F_old(3,3)
	   Double Precision :: U(3,3)
	
	!  DGD Step-1
	   Double Precision :: eye(ndir, ndir) ! Identity
	   Double Precision :: GLStrain(3, 3)
	
	!  DGD Step-2
	   Double Precision :: stress(3,3)
	   Double Precision :: stressABQ(3,3)

	!  DGD Step-3
	   Double Precision :: Cauchy(3,3)
	
	!  DGD Step-4
	   Double Precision :: Rot(3,3) 
	   Double Precision :: CauchyABQ(3,3)
	   
	   ! Double Precision :: eps(ndir, ndir), eps_old(ndir, ndir)
       ! Double Precision :: FIt, FIc, FIt_old, FIc_old, FIt_max, FIc_max
       ! Double Precision :: sig_f0, eps_f0, eps_f
	   ! Double Precision :: sig_f0_T, eps_f0_T, eps_ff_T, eps_f_T
	   ! Double Precision :: sig_f0_C, eps_f0_C, eps_ff_C, eps_f_C
	   ! Double Precision :: D, d1Plus, d1Minus
	   
	   integer i,j
C
C INITIALISATION: MATERIAL CARD PARAMETERS
C =========================================
C
      E11 = props(1) !E11 = 61646
      E22 = props(2) !E33 = 13368
      E33 = props(3) !E22 = 61646
      v12 = props(4) !V13 = 0.3070
      v13 = props(5) !V12 = 0.3187
      v23 = props(6) !V32 = 0.0667 
      G12 = props(7) !G13 = 4575
      G13 = props(8) !G12 = 23373
      G23 = props(9) !G32 = 4575  
C	  
	  NLS = props(10)
	  A   = props(11) 
	  B   = props(12) 
C	 
	  Xt  = props(13)  
	  Xc  = props(14)  
C	  
      Gft = props(15)  
	  Gfc = props(16)	  
C
	  DelCriterion = props(17)

C
C AXISYMMETRIC MODEL STIFFNESS MATRIX 
C ====================================
C
      v21 = v12*E22/E11
      v31 = v13*E33/E11
      v32 = v12	
	  
      n = E11/E22
	  m = G12/E22
C
	  delta = (one + v13)*(one - v13 - two*n*v21**2)
C
	  C(1,1) = (E22/delta) * n * (one - n*v21**2)
	  C(1,2) = (E22/delta) * n * v21 * (one + v13)
	  C(1,3) = (E22/delta) * n * (v13 + n* v21**2)
	  C(1,4) = zero
	  C(2,1) = C(1,2)
	  C(2,2) = (E22/delta) * (one - v13**2)
	  C(2,3) = (E22/delta) * n * v21 * (one + v13)
	  C(2,4) = zero
	  C(3,1) = C(1,3)
	  C(3,2) = C(2,3)
	  C(3,3) = (E22/delta) * n * (one - n*v21**2)
	  C(3,4) = zero
	  C(4,1) = zero
	  C(4,2) = zero
	  C(4,3) = zero
	  C(4,4) = G12   
C	  
C START CONTINNUM DAMAGE MODEL
C ============================
C
C INITIAL ELASTIC STEP FOR ABAQUS TESTS
C -------------------------------------	  
C        
      IF (totalTime.eq.zero) THEN
      
          DO i = 1, nblock 
		  
		      stressNew(i,1)= stressOld(i,1)+ C(1,1)*straininc(i,1)
     *           + C(1,2)*straininc(i,2) + C(1,3)*straininc(i,3)
        
		      stressNew(i,2)= stressOld(i,2)+ C(2,1)*straininc(i,1)
     *           + C(2,2)*straininc(i,2) + C(2,3)*straininc(i,3) 

              stressNew(i,3)=stressOld(i,3) + C(3,1)*straininc(i,1)
     *           + C(3,2)*straininc(i,2) + C(3,3)*straininc(i,3)
	 
              stressNew(i,4) = stressOld(i,4) + 
     *			  two*C(4,4)*straininc(i,4)  
              	
			  
          ENDDO
C
C IF NOT INITIAL STEP, CAL. STRESSES ACCRODING TO CDM MODEL
C =========================================================
C      
      ELSE   
C         
          DO i = 1, nblock
		      
			  ! ------------------------------------------------------------ !
			  !     Deformation Gradient Tensor and Right Stretch Tensor     !
			  ! ------------------------------------------------------------ !
			  
			  U = Vec2Matrix(stretchNew(i,:)) ! Right stretch tensor in matrix format
			  F = Vec2Matrix(defgradNew(i,:)) ! Deformation gradient tensor in matrix format
			  F_old = Vec2Matrix(defgradOld(i,:))
			  
			  WRITE(*,*) 'stretchNew=', stretchNew
			  WRITE(*,*) 'defgradNew=', defgradNew
			  WRITE(*,*) 'U=',U
			  WRITE(*,*) 'U(1,1)=', U(1,1)
			  WRITE(*,*) 'U(1,2)=', U(1,2)
			  WRITE(*,*) 'U(1,3)=', U(1,3)
			  WRITE(*,*) 'U(2,1)=', U(2,1)
			  WRITE(*,*) 'U(2,2)=', U(2,2)
			  WRITE(*,*) 'U(2,3)=', U(2,3)
			  WRITE(*,*) 'U(3,1)=', U(3,1)
			  WRITE(*,*) 'U(3,2)=', U(3,2)
			  WRITE(*,*) 'U(3,3)=', U(3,3)
			  WRITE(*,*) 'F=',F
			  WRITE(*,*) 'F(1,1)=', F(1,1)
			  WRITE(*,*) 'F(1,2)=', F(1,2)
			  WRITE(*,*) 'F(1,3)=', F(1,3)
			  WRITE(*,*) 'F(2,1)=', F(2,1)
			  WRITE(*,*) 'F(2,2)=', F(2,2)
			  WRITE(*,*) 'F(2,3)=', F(2,3)
			  WRITE(*,*) 'F(3,1)=', F(3,1)
			  WRITE(*,*) 'F(3,2)=', F(3,2)
			  WRITE(*,*) 'F(3,3)=', F(3,3)
			  
			
			  ! ------------------------------------------------------------ !
			  !     Step-1: The Green-Lagrange strain tensor, GLStrain       !
			  ! ------------------------------------------------------------ !
			  eye = zero
			  Do j = 1,3
				  eye(j,j) = one
			  END DO
			  
			  ! The Green-Lagrange strain
			  GLStrain = (MATMUL(TRANSPOSE(F), F) - eye)*half ! Eq.4
			  WRITE(*,*) 'GLStrain=', GLStrain
			  
			  ! The Green-Lagrange strain components
			  stateNew(i,1) = GLStrain(1,1)
			  stateNew(i,2) = GLStrain(2,2)
			  stateNew(i,3) = GLStrain(3,3)
			  stateNew(i,4) = two*GLStrain(1,2) !Engineering shear strain
			  
			  ! ------------------------------------------------------------ !
			  !     Step-2: The second Piola-Kirchoff stress, S              !
			  ! ------------------------------------------------------------ !
			  
			  ! 2nd Piola-Kichoff stress in matrix format
			  stress = Hooke(C, GLStrain, nshr) ! Eq.5
			  WRITE(*,*) 'stress=', stress
			  
			  ! ------------------------------------------------------------ !
			  !     CDM (INPUT: strain,stress; OUTPUT: degraded stress       !
			  ! ------------------------------------------------------------ !
			  
			  ! ------------------------------------------------------------ !
			  !     Step-3: The Cauchy stress in the current configuration   !
			  ! ------------------------------------------------------------ !
			  
			  ! The cauchy stress
			  Cauchy = convertToCauchy(stress, F) ! Eq.6
			  WRITE(*,*) 'Cauchy=', Cauchy
			  
			  ! ------------------------------------------------------------ !
			  !     Step-4: Determine and store the stress tensor            !
			  ! ------------------------------------------------------------ !
			  
			  ! Rotation tensor
			  Rot = MATMUL(F, MInverse(U))
			  WRITE(*,*) 'Rot=', Rot
			  
			  ! Cauchy stress in the current configuration
			  CauchyABQ = MATMUL(Rot, MATMUL(Cauchy, TRANSPOSE(Rot)))
			  
			  ! Convert to vector format
			  stressNew(i,:) = Matrix2Vec(CauchyABQ, nshr)
			  WRITE(*,*) 'CauchyABQ=', CauchyABQ
			  WRITE(*,*) 'stressNew(i,:)', stressNew(i,:)
			  
	      ENDDO
C	  
	  ENDIF
C	
	  RETURN 	  
	  END



! ----------------------------------------------------------------- !
!   Function GLStrain: to cal the Green-Lagrange strain, GLStrain   !
! ----------------------------------------------------------------- !
	  Function GLStrain(F, ndir)
		! Computes the Green-Lagrange strain from the deformation gradient

		! Input
		Double Precision, intent(IN) :: F(3,3)
		Integer, intent(IN) :: ndir

		! Output
		Double Precision :: GLStrain(ndir,ndir)

		! Locals
		Double Precision :: eye(ndir,ndir)  ! Identity
		Double Precision, parameter :: zero=0.d0, one=1.d0, half=0.5d0
		! -------------------------------------------------------------------- !

		! Initialize identity matrix
		eye = zero; Do I = 1,3; eye(I,I) = one; End Do

		GLStrain = (MATMUL(TRANSPOSE(F), F) - eye)*half

		Return
	  End Function GLStrain

!-------------------------------------------------------------------!
!      Subroutine vumatcoh:                                         !
!-------------------------------------------------------------------!       

#include 'Abaqus_Definitions.f'
#include 'Fatigue_Globals.f'
#include 'Mesh_Utils.f'
#include 'Fatigue_ANN.f'
#include 'Fatigue.f'
#include 'CZM.f'

      !> VUMAT Main entry point from Abaqus
      subroutine vumatcoh(
     *     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew )

        use CZM, only: vumat_cohesive_fatigue
        include 'vaba_param.inc'

        dimension jblock(*), props(nprops),density(*), coordMp(*),
     *     charLength(*), strainInc(*), relSpinInc(*), tempOld(*),
     *     stretchOld(*), defgradOld(*), fieldOld(*), stressOld(*),
     *     stateOld(*), enerInternOld(*),enerInelasOld(*),
     *     tempNew(*), stretchNew(*), defgradNew(*), fieldNew(*),
     *     stressNew(*), stateNew(*), enerInternNew(*), enerInelasNew(*)

        character*80 cmname

        call vumat_cohesive_fatigue ( jblock(1),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(2), jblock(3), jblock(4), jblock(5))

      end subroutine vumatcoh

      
      !> Vexternaldb entry from Abaqus
      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

        use Abaqus_Definitions, only: j_int_StartAnalysis
        use CZM, only: CZM_initialisation
        include 'vaba_param.inc'

        dimension i_Array(niArray), r_Array(nrArray)

        ! Initialisation at start of the analysis
        if (lOp .eq. j_int_StartAnalysis) then
          call CZM_initialisation()
        endif

      end subroutine vexternaldb	

	  
