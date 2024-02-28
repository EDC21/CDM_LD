	Module stress_Mod
	  ! Module for building stiffness tensors, tensor transformations, and converting different definitions of stress

	Contains


	  Pure Function convertToCauchy(stress, F) result(Cauchy)
		! Converts the stress to Cauchy stress

		Use matrixAlgUtil_Mod

		! Input
		Double Precision, intent(IN) :: stress(3,3)
		Double Precision, intent(IN) :: F(3,3)            ! Deformation gradient

		! Output
		Double Precision :: Cauchy(3,3)
		! -------------------------------------------------------------------- !

		Cauchy = MATMUL(F, MATMUL(stress, TRANSPOSE(F)))/MDet(F)

		Return
	  End Function convertToCauchy


	  Pure Function Hooke(C, strain, nshr) result(stress)
		! Calculates the energy conjugate stress based on the input strain and stiffness tensor

		Use matrixAlgUtil_Mod

		! Input
		Double Precision, intent(IN) :: C(3+nshr,3+nshr)         ! Stiffness
		Double Precision, intent(IN) :: strain(3,3)
		Integer, intent(IN) :: nshr

		! Output
		Double Precision :: stress(3,3)

		! Locals
		Double Precision :: strainVec(3+nshr), stressVec(3+nshr)
		Double Precision, parameter :: zero=0.d0, two=2.d0
		! -------------------------------------------------------------------- !

		stress = zero

		! Convert strain to vector format, engineering shear strain
		strainVec = Matrix2Vec(strain, nshr)
		Do I=4,3+nshr; strainVec(I) = two*strainVec(I); End Do

		! Compute stress vector
		stressVec = MATMUL(C,strainVec)

		! Convert back to matrix format
		stress = Vec2Matrix(stressVec)

		Return
	  End Function Hooke

	End Module stress_Mod
