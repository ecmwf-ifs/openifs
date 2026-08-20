! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

!!! lapack and blas interfaces                                  -*- mode: F90; -*-

!!! blas routines used in mcmc-code
!!! dcopy (dsymv) dtrmv (dsyr)
!!! daxpy dswap (if linpack)

!!! integer is now 4 byte, this can be a problem in some implementations

!!! BLAS
! dcopy
INTERFACE 
   SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: DX(*)
     INTEGER(KIND=4) :: INCX
     REAL(KIND=8) :: DY(*)
     INTEGER(KIND=4) :: INCY
   END SUBROUTINE DCOPY
END INTERFACE

INTERFACE 
   SUBROUTINE DSCAL(N,A,X,INCX)
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: A
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
   END SUBROUTINE DSCAL
END INTERFACE

INTERFACE 
   SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: ALPHA
     REAL(KIND=8) :: A(LDA,*)
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
     REAL(KIND=8) :: BETA
     REAL(KIND=8) :: Y(*)
     INTEGER(KIND=4) :: INCY
   END SUBROUTINE DSYMV
END INTERFACE

INTERFACE 
   SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     CHARACTER(LEN=1) :: TRANS
     CHARACTER(LEN=1) :: DIAG
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: A(LDA,*)
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
   END SUBROUTINE DTRMV
END INTERFACE

INTERFACE 
   SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: TRANS
     INTEGER(KIND=4) :: M
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: ALPHA
     REAL(KIND=8) :: A(LDA,*)
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
     REAL(KIND=8) :: BETA
     REAL(KIND=8) :: Y(*)
     INTEGER(KIND=4) :: INCY
   END SUBROUTINE DGEMV
END INTERFACE
 
 
 ! rank one update of symmetric matrix
INTERFACE 
   SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: ALPHA
     REAL(KIND=8) :: X(*)
     INTEGER(KIND=4) :: INCX
     REAL(KIND=8) :: A(LDA,*)
   END SUBROUTINE DSYR
END INTERFACE

!!! LAPACK
! dpotrf, dpotri

INTERFACE 
   SUBROUTINE DPOTRF(UPLO,N,A,LDA,INFO)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: A(LDA,*)
     INTEGER(KIND=4) :: INFO
   END SUBROUTINE DPOTRF
END INTERFACE

INTERFACE 
   SUBROUTINE DPOTRI(UPLO,N,A,LDA,INFO)
     INTEGER(KIND=4) :: LDA
     CHARACTER(LEN=1) :: UPLO
     INTEGER(KIND=4) :: N
     REAL(KIND=8) :: A(LDA,*)
     INTEGER(KIND=4) :: INFO
   END SUBROUTINE DPOTRI
END INTERFACE

interface
   subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
     implicit none
     character(len=1), intent(in) ::  jobu, jobvt
     integer(kind=4), intent(in)  ::  lda, ldu, ldvt, lwork, m, n
     integer(kind=4), intent(out) ::  info
     real(kind=8),intent(inout)   ::  a(lda,*)
     real(kind=8),intent(out)   ::  s(*), u(ldu,*),vt(ldvt,*),work(*)
   end subroutine dgesvd
end interface

!!! symmetric linear solver
!! C = alpha*A*A'+beta*C
interface
   subroutine dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
     implicit none
     character(len=1), intent(in) :: uplo, trans
     integer(kind=4), intent(in) :: n, k, lda, ldc
     real(kind=8), intent(in) :: alpha, beta
     real(kind=8), intent(in) :: a(lda,*)
     real(kind=8), intent(inout) :: c(ldc,*)
   end subroutine dsyrk
end interface
!! B = A\B, A general n*n
interface
   subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
     implicit none
     integer(kind=4), intent(in) :: n, nrhs, lda, ldb
     real(kind=8), intent(inout) :: a(lda,n), b(ldb,n)
     integer(kind=4), intent(out) :: ipiv(n)
     integer(kind=4), intent(out) :: info
   end subroutine dgesv
end interface
!! B = A\B, A symmetric n*n
interface
   subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info) 
     implicit none
     character(len=1), intent(in) :: uplo
     integer(kind=4), intent(in) :: n, nrhs, lda, ldb
     real(kind=8), intent(inout) :: a(lda,n), b(ldb,n)
     integer(kind=4), intent(out) :: info
   end subroutine dposv
end interface

interface
   subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
     implicit none
     character(len=1), intent(in) :: trans
     integer(kind=4), intent(in) :: m, n, nrhs, lda, ldb, lwork
     real(kind=8), intent(inout) :: a(lda,n), b(ldb,nrhs)
     real(kind=8), intent(inout) :: work(lwork)
     integer(kind=4), intent(out) :: info     
   end subroutine dgels
end interface
