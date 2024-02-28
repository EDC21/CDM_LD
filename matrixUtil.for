	Module matrixAlgUtil_Mod
	  ! Generic utilities for manipulating vectors and matrices

	Contains

	  Pure Function Vec2Matrix(vector)
		! Converts tensors stored in vector format to a matrix format

		! Arguments
		Double Precision, intent(IN) :: vector(:)

		! Output
		Double Precision :: Vec2Matrix(3,3)

		! Locals
		Double Precision, parameter :: zero=0.d0
		! -------------------------------------------------------------------- !

		Vec2Matrix = zero
		Vec2Matrix(1,1) = vector(1)
		Vec2Matrix(2,2) = vector(2)
		Vec2Matrix(3,3) = Vector(3)
		Vec2Matrix(1,2) = vector(4)

		! 2D or 3D
		If (size(vector) > 5) Then ! 3D
		  Vec2Matrix(2,3) = vector(5)

		  ! Symmetric or nonsymmetric
		  If (size(vector) == 6) Then ! 3D, Symmetric
			Vec2Matrix(1,3) = vector(6)

			Vec2Matrix(2,1) = Vec2Matrix(1,2)
			Vec2Matrix(3,2) = Vec2Matrix(2,3)
			Vec2Matrix(3,1) = Vec2Matrix(1,3)

		  Else ! 3D, Nonsymmetric
			Vec2Matrix(3,1) = vector(6)
			Vec2Matrix(2,1) = vector(7)
			Vec2Matrix(3,2) = vector(8)
			Vec2Matrix(1,3) = vector(9)
		  End If

		Else ! 2D
		  If (size(vector) == 5) Then  ! 2D, Nonsymmetric
			Vec2Matrix(2,1) = vector(5)
		  Else                           ! 2D, Symmetric
			Vec2Matrix(2,1) = Vec2Matrix(1,2)
		  End If
		End IF

		Return
	  End Function Vec2Matrix


	  Pure Function Matrix2Vec(mat, nshr, symmetric)
		! Converts a symmetric tensor stored in matrix format (3,3) to a vector

		! Arguments
		Double Precision, intent(IN) :: mat(3,3)
		Integer, intent(IN) :: nshr
		Logical, intent(IN), optional :: symmetric        ! True for symmetric matrix (default=True)

		! Output
		Double Precision, allocatable :: Matrix2Vec(:)

		! Locals
		Integer :: output_length
		Logical :: sym
		Double Precision, parameter :: zero=0.d0
		! -------------------------------------------------------------------- !

		! Optional argument for symmetry
		If (present(symmetric)) Then
		  sym = symmetric
		Else
		  sym = .TRUE.
		End If

		! Initialize output
		If (sym) Then
		  output_length = 3+nshr
		Else
		  If (nshr == 1) Then
			output_length = 5
		  Else
			output_length = 9
		  End If
		End If
		Allocate(Matrix2Vec(output_length))
		Matrix2Vec = zero

		! 2D components
		Do I=1,3; Matrix2Vec(I) = mat(I,I); End Do
		Matrix2Vec(4) = mat(1,2)

		! 3D
		If (sym) Then
		  If (nshr > 1) Then
			Matrix2Vec(5) = mat(2,3)
			Matrix2Vec(6) = mat(3,1)
		  End If
		Else
		  If (nshr > 1) Then
			Matrix2Vec(5) = mat(2,3)
			Matrix2Vec(6) = mat(3,1)
			Matrix2Vec(7) = mat(2,1)
			Matrix2Vec(8) = mat(3,2)
			Matrix2Vec(9) = mat(1,3)
		  Else
			Matrix2Vec(5) = mat(2,1)
		  End If
		End If

		Return
	  End Function Matrix2Vec


	  Pure Function VCmp(vec1, vec2)
		! Checks if vec1 and vec2 have identical components
		! Not optimized for large vectors

		! Arguments
		Double Precision, intent(IN) :: vec1(:), vec2(:)

		! Output
		Logical :: VCmp

		! Locals
		Double Precision, parameter :: eps=1.d-30
		! -------------------------------------------------------------------- !

		If(size(vec1) /= size(vec2)) Then
		  VCmp = .FALSE.
		  Return
		Else
		  Do I=1, size(vec1)
			If(abs(vec1(I) - vec2(I)) > eps) Then
			  VCmp = .FALSE.
			  Return
			End If
		  End Do
		End If

		VCmp = .TRUE.

		Return
	  End Function VCmp


	  Pure Function MInverse(mat)
		! Finds the inverse of a 3x3 matrix

		! Arguments
		Double Precision, intent(IN) :: mat(3,3)

		! Output
		Double Precision :: MInverse(3,3)
		! -------------------------------------------------------------------- !

		MInverse(1,1) = mat(3,3)*mat(2,2) - mat(3,2)*mat(2,3)
		MInverse(2,2) = mat(3,3)*mat(1,1) - mat(3,1)*mat(1,3)
		MInverse(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

		MInverse(1,2) = -(mat(3,3)*mat(1,2) - mat(3,2)*mat(1,3))
		MInverse(2,1) = -(mat(3,3)*mat(2,1) - mat(2,3)*mat(3,1))
		MInverse(1,3) =   mat(2,3)*mat(1,2) - mat(2,2)*mat(1,3)
		MInverse(3,1) =   mat(3,2)*mat(2,1) - mat(2,2)*mat(3,1)
		MInverse(2,3) = -(mat(2,3)*mat(1,1) - mat(2,1)*mat(1,3))
		MInverse(3,2) = -(mat(3,2)*mat(1,1) - mat(1,2)*mat(3,1))

		Do I=1,3
		  Do J=1,3
			MInverse(I,J) = MInverse(I,J) / MDet(mat)
		  End Do
		End Do

		Return
	  End Function MInverse

	  Pure Function MDet(mat)
		! Finds the determinant of a 3x3 matrix

		! Arguments
		Double Precision, intent(IN) :: mat(3,3)

		! Output
		Double Precision :: MDet
		! -------------------------------------------------------------------- !

		MDet = mat(1,1)*(mat(3,3)*mat(2,2) - mat(3,2)*mat(2,3)) - 
     *		mat(2,1)*(mat(3,3)*mat(1,2) - mat(3,2)*mat(1,3)) + 
     *		mat(3,1)*(mat(2,3)*mat(1,2) - mat(2,2)*mat(1,3))

		Return
	  End Function MDet


	  Pure Function OuterProduct(u, v)
		! Calculates the outer product of two vectors

		! Arguments
		Double Precision, intent(IN) :: u(:), v(:)

		! Output
		Double Precision :: OuterProduct(size(u),size(v))
		! -------------------------------------------------------------------- !

		Do I = 1, size(u)
		  Do J = 1, size(v)
			OuterProduct(I,J) = u(I)*v(J)
		  End Do
		End Do

		Return
	  End Function OuterProduct


	  Pure Function CrossProduct(a, b)
		! Calculates the cross product of two vectors with length 3

		! Arguments
		Double Precision, intent(IN) :: a(3), b(3)

		! Output
		Double Precision :: CrossProduct(3)
		! -------------------------------------------------------------------- !

		CrossProduct(1) = a(2)*b(3) - a(3)*b(2)
		CrossProduct(2) = a(3)*b(1) - a(1)*b(3)
		CrossProduct(3) = a(1)*b(2) - a(2)*b(1)

		Return
	  End Function CrossProduct


	  Pure Function Norm(u)
		! Normalizes a vector

		! Arguments
		Double Precision, intent(IN) :: u(:)

		! Output
		Double Precision :: Norm(size(u))

		! Locals
		Double Precision :: u_length
		! -------------------------------------------------------------------- !

		u_length = Length(u)

		Do I = 1, size(u)
		  Norm(I) = u(I)/u_length
		End Do

		Return
	  End Function Norm


	  Pure Function Length(u)
		! Calculates the length of a vector

		! Arguments
		Double Precision, intent(IN) :: u(:)

		! Output
		Double Precision :: Length
		! -------------------------------------------------------------------- !

		Length = 0.d0

		Do I = 1, size(u)
		  Length = Length + u(I)*u(I)
		End Do

		Length = SQRT(Length)

		Return
	  End Function Length




	End Module
