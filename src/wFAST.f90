!***************************************************************************
! wFAST is a three-dimensional finite element method program 
! that analyzes the stability of faults (discontinuities) in rock mass.
!***************************************************************************

!!-------------------------------------------------------
!! How to install PETSc with Intel OneAPI.
!!   In $HOME/WFAST directory,
!!     wgets https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.22.3.tar.gz
!!     tar zxvf petsc-3.22.3.tar.gz
!!     cd petsc-3.22.3
!!     ./configure CC=mpiicx CXX=mpiicpx FC=mpiifx --prefix=$HOME/WFAST/petsc-intel
!!     make PETSC_DIR=$HOME/WFAST/petsc-3.22.3-intel PETSC_ARCH=arch-linux-c-debug all
!!     make PETSC_DIR=$HOME/WFAST/petsc-3.22.3-intel PETSC_ARCH=arch-linux-c-debug install
!!   Add this PATH to environment variable LD_LIBRARY_PATH as,
!!     export LD_LIBRARY_PATH=$HOME/WFAST/petsc-intel/lib:$LD_LIBRARY_PATH
!!
!! How to compile wFAST with Intel OneAPI.
!!   In $HOME/WFAST/src_rev directory,
!!     make
!!   wFAST is linked with both PARDISO and PETSc library.
!!
!! How to run wFAST with input file test.cfg.
!!   In $HOME/WFAST/test directory,
!!     ../src_rev/wFAST test.cfg
!!-------------------------------------------------------

!!-------------------------------------------------------
!! How to use PETSc installed in Fujitsu super computer.
!!
!!   Activate spack command first,
!!     . /vol0004/apps/oss/spack/share/spack/setup-env.sh
!!
!!   Activate PETSc, it will take some minites.
!!     spack load petsc
!!   Check PETSc is found.
!!     ls $PETSC_DIR
!!
!! How to compile wFAST with Fujitsu compiler.
!!   In $HOME/WFAST/src_rev directory,
!!   Edit Makefile to use parameters for Fujitsu. Then, 
!!     make
!!   wFAST is linked with PETSc library(not PARDISO).
!!
!!-------------------------------------------------------  

#ifndef ENABLE_PARDISO
#ifndef ENABLE_PETSC
#define ENABLE_PARDISO
#endif
#endif
   
#ifdef ENABLE_PARDISO
#warning "ENABLE PARDISO"
include 'mkl_pardiso.f90'
#endif

#ifdef ENABLE_PETSC
#warning "ENABLE PETSC"
#include <petsc/finclude/petscsnes.h>
#endif

module wFAST
#ifdef ENABLE_PETSC
   use petscsnes
#endif
   implicit none

   !--- constants and limits.
   integer, parameter :: ndim = 3 !< dimension of coordinates.
   integer, parameter :: nst  = 6 !< number of stress (strain) components in Voigt notation.
   real(8), parameter :: M_PI      = 3.141592653589793238462d0 !< Pi, circular constant.
   real(8), parameter :: M_SQRT1_3 = 0.577350269189625764509d0 !< 1/sqrt(3)

   !>
   ! material propertiy for solid element.
   type propertys
      real(8) :: E       !< Young's modulus of bulk rock
      real(8) :: nu      !< Poisson's ratio of bulk rock
      real(8) :: rho     !< density of bulk rock
      real(8) :: Em      !< Young's modulus of matrix
      real(8) :: num     !< Poisson's ratio of matrix
      real(8) :: alpha   !< thermal expansion coeff
      real(8) :: mu      !< static friction coefficient for unmapped faults
      real(8) :: So      !< cohesion for unmapped faults
   end type propertys

   !>
   ! material propertiy for joint element.
   type propertyj
      real(8) :: ksn(3)  !< joint stiffness
      real(8) :: gs, gn  !< viscoplastic coefficients
      real(8) :: multksn !< update factor for joint stiffness in failure
      real(8) :: T       !< tensile strength
      real(8) :: mu      !< static friction coefficient
      real(8) :: So      !< joint cohesion
      real(8) :: dsc     !< threshold stress
      real(8) :: mup     !< dynamic friction coefficient
      real(8) :: dcrit     !< critical slip distance
      real(8) :: shear_mod !< shear modulus
   end type propertyj

   !>
   ! input parameters.
   type Param_type
      integer :: nods    !< number of nodes per solid element.
      integer :: nodb    !< number of nodes per boundary element (surface element).
      integer :: nodj    !< number of nodes per joint element.
      integer :: neles   !< total number of solid elements in the mesh.
      integer :: nelej   !< total number of joint elements in the mesh.
      integer :: nnode   !< total number of nodes in the mesh.
      integer :: nips    !< number of integration points per solid element.
      integer :: nipb    !< number of integration points per boundary element.
      integer :: nprops  !< number of properties for solids.
      integer :: npropj  !< number of properties for joints.
      integer :: nips_stress !< number of integration points per solid element for evaluating stress.
      integer :: nipj_stress !< number of integration points per joint element for evaluating stress.

      real(8) :: convergence     !< convergence limit for nonlinear calculations.
      integer :: use_stress_grad !< stress with gradient.
      integer :: use_press_grad  !< pressure with gradient.
      integer :: use_temp_grad   !< temperature with gradient.
      real(8) :: gravity         !< gravity constant.

      type(propertys), pointer :: props(:) !< material properties for solids (nprops+more).
      type(propertyj), pointer :: propj(:) !< material properties for joints (npropj+more).

      integer :: ntime   !< number of time steps.
      integer :: niter   !< number of nonlinear iterations.

      integer :: nlayer  !< number of layer.
      character(len=20), allocatable :: layer_name(:)  !< layer name (nlayer).
!!$      character(len=20), allocatable :: layer_names(:) !< layer name for solids (nprops).
!!$      character(len=20), allocatable :: layer_namej(:) !< layer name for joints (npropj).
      integer, allocatable :: layer_nele(:)    !< number of elements in each layer (nlayer).
!!$      integer, allocatable :: layer_neles(:)   !< number of elements in each solid layer (nlayer).
!!$      integer, allocatable :: layer_nelej(:)   !< number of elements in each joint layer (nlayer).
      integer, allocatable :: layer_enode(:,:) !< elements in each layer (neles+nelej,nlayer).
   end type Param_type
   type(Param_type) :: param

   !>
   ! FEM geometry.
   type Geometry_type
      real(8), allocatable :: coord(:,:)     !< global nodal coordinates (ndim,nnode).
      real(8), allocatable :: displ(:,:)     !< global nodal displacements (ndim,nnode).
      real(8), allocatable :: displ_inc(:,:) !< increment of nodal displacement, solution of equation (ndim,nnode).

      integer, pointer     :: nodfree(:,:)   !< free coords (ndim,nnode).
      integer, pointer     :: codfree(:)     !< free coords (ndim*nnode).
      real(8), allocatable :: valfix(:,:)    !< fixed value (ndim,nnode).

      integer, allocatable :: enodes(:,:)    !< global element node numbers vector for solids (nods,neles).
      integer, allocatable :: enodej(:,:)    !< global element node numbers vector for joints (nodj,nelej).

      integer, allocatable :: etypes(:)      !< element material index for solids (neles).
      integer, allocatable :: etypej(:)      !< element material index for joints (nelej).

      real(8), allocatable :: ecenters(:,:)  !< element center for solids (ndim,neles).
      real(8), allocatable :: ecenterj(:,:)  !< element center for joints (ndim,nelej).

      character(len=20), allocatable :: elem_layer_names(:) !< layer name for solids (neles).
      character(len=20), allocatable :: elem_layer_namej(:) !< layer name for joints (nelej).

      !--- surface load
      integer, allocatable :: surfj_elem(:)    !< joint elements with surface load.
      integer, allocatable :: surfj_enode(:,:) !< joint nodes with surface load.
      
      !--- coord index in each FEM element.
      integer, allocatable :: esteers(:,:)           !< global element steering vector for solids (ndim*param%nods,param%neles).
      integer, allocatable :: esteerj(:,:)           !< global element steering vector for joints (ndim*param%nodj,param%nelej).

      !--- matrix element index of each coords in FEM element.
      integer, allocatable :: eindexs(:,:,:)         !< index of i-j element in the matrix for solids (ndim*nods,ndim*nods,neles).
      integer, allocatable :: eindexj(:,:,:)         !< index of i-j element in the matrix for joints (ndim*nodj,ndim*nodj,nelej).

      !--- matrix element index.
      integer :: nnz                                 !< total non-zero matrix elements.
      integer, allocatable :: ia(:), ja(:)           !< matrix index for PARDISO (ndim*param%nnode+1), (nnz).
#ifdef ENABLE_PETSC
      PetscInt, allocatable :: irn(:), jcn(:), idx(:) !< matrix index for PETSc (nnz), (ndim*param%nnode).
#endif
   end type Geometry_type
   type(Geometry_type) :: geom

   !--- fields.
   real(8), allocatable :: surfj_load(:,:,:)    !< surface load for joints (ndim,nodb,*).
   real(8), allocatable :: stresss_ini(:,:)     !< initial (residual) stresses for solids (nst,neles).
   real(8), allocatable :: stressj_ini(:,:)     !< initial (residual) stresses for joints (nst,nelej).

   real(8), allocatable :: presss_ini(:)        !< initial pore pressure for solids (neles).
   real(8), allocatable :: pressj_ini(:)        !< initial pore pressure for joints (nelej).
   real(8), allocatable :: presss(:)            !< pore pressure (neles).
   real(8), allocatable :: pressj(:)            !< pore pressure (nelej).

   real(8), allocatable :: strains_ini(:,:)     !< initial (residual) strains (nst,neles).

   logical, allocatable :: failurej(:)          !< joint stability (nelej).

   real(8), allocatable :: tractionj_ini(:,:,:) !< initial joint traction (ndim,nipb,nelej).
   real(8), allocatable :: tractionj(:,:,:)     !< joint traction (ndim,nipb,nelej).

   real(8), allocatable :: coulombj_ini(:,:)    !< initial coulomb adhesion stress (nipb,nelej).
   real(8), allocatable :: coulombj(:,:)        !< coulomb adhesion stress (nipb,nelej).
   real(8), allocatable :: coulombj_dif(:,:)    !< change of coulomb adhesion stress (nipb,nelej).
   
   real(8), allocatable :: cohesionj(:,:)       !< joint cohesion (nipb,nelej).
   real(8), allocatable :: tensilej(:,:)        !< joint tensile strength (nipb,nelej).
   real(8), allocatable :: stiffj(:,:,:)        !< joint stiffness (ndim,nipb,nelej).

   real(8), allocatable :: energy(:,:)          !< strain energy density (nips_stress,neles).

   real(8), allocatable :: robust(:,:)          !< robustness index (nips_stress,neles).
   real(8), allocatable :: robust_ini(:,:)      !< initial robustness index (nips_stress,neles).
   real(8), allocatable :: robust_dif(:,:)      !< change of robustness index (nips_stress,neles).

   real(8), allocatable :: load(:,:)            !< nodal loads (ndim,nnode).
   real(8), allocatable :: stiffmat(:)          !< global stiffness matrix in vector form (nnz).

   real(8), allocatable :: seismic_mom(:)       !< seismic moment (nelej).
   real(8), allocatable :: seismic_int(:)       !< seismic moment intensity (nelej).

   real(8), allocatable :: temp(:)              !< temperature (neles).
   real(8), allocatable :: temp_ini(:)          !< initial temperature (neles).
   real(8), allocatable :: density_dif(:)       !< change of density (neles).

   ! solver for the linear equation.
   character(len=8) :: solver = "PARDISO"
!!$   character(len=8) :: solver = "PETSc"
   logical :: upper_only
   
   ! filenames of input and output files.
   character(len=256) :: filename_gid = ""
   character(len=256) :: filename_dat = ""
   character(len=256) :: filename_err = ""
   character(len=256) :: filename_presss_ini = ""
   character(len=256) :: filename_pressj_ini = ""
   character(len=256) :: filename_point_load = ""
   character(len=256) :: filename_surfj_load = ""
   character(len=256) :: filename_presss_load = ""
   character(len=256) :: filename_pressj_load = ""
   character(len=256) :: filename_temp_load = ""
   character(len=256) :: filename_den_dif = ""
   character(len=256) :: filename_log = ""
   character(len=256) :: filename_res = ""
   character(len=256) :: filename_postmesh = ""
   character(len=256) :: filename_postres = ""


   ! file IOs.
   integer :: iunit_gid
   integer :: ounit_dat
   integer :: ounit_err
   integer :: ounit_presss_ini
   integer :: ounit_pressj_ini
   integer :: ounit_point_load
   integer :: ounit_surfj_load
   integer :: ounit_presss_load
   integer :: ounit_pressj_load
   integer :: ounit_temp_load
   integer :: ounit_den_dif
   integer :: iunit_dat
   integer :: iunit_point_load
   integer :: iunit_surfj_load
   integer :: iunit_presss_load
   integer :: iunit_pressj_load
   integer :: iunit_temp_load
   integer :: iunit_den_dif
   
   integer :: ounit_log
   integer :: ounit_res
   integer :: ounit_postmesh
   integer :: ounit_postres

contains
   real(8) function perpendicular( v, Lmat )
      real(8), intent(in) :: v(3)
      real(8), intent(in), optional :: Lmat(3,3)
      real(8) :: u(3)

      if( present(Lmat) ) then
         u = matmul(v,Lmat)
         perpendicular = u(3)
      else
         perpendicular = v(3)
      end if
   end function perpendicular

   real(8) function parallel( v, Lmat )
      real(8), intent(in) :: v(3)
      real(8), intent(in), optional :: Lmat(3,3)
      real(8) :: u(3)

      if( present(Lmat) ) then
         u = matmul(v,Lmat)
         parallel = sqrt(u(1)**2+u(2)**2)
      else
         parallel = sqrt(v(1)**2+v(2)**2)
      end if
   end function parallel

   real(8) function norm( v )
      real(8), intent(in) :: v(:)

      norm = sqrt(sum(v(:)**2))
   end function norm

   function cross_product( a, b ) result(c)
      real(8), intent(in) :: a(3), b(3)
      real(8) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   end function cross_product

   function outer_product( a, b ) result(c)
      real(8), intent(in) :: a(:), b(:)
      real(8) :: c(size(a),size(b))
      integer :: i, j

      do i=1, size(a)
         do j=1, size(b)
            c(i,j) = a(i)*b(j)
         end do
      end do
   end function outer_product
   
   function determinants( jac ) result(det)
      !-----------------------------------------------------------------------
      ! This function returns the determinant of a 3x3 Jacobian matrix.
      !-----------------------------------------------------------------------
      real(8), intent(in) :: jac(ndim,ndim)
      real(8) :: uni(ndim,ndim)
      real(8) :: normal(ndim)
      real(8) :: det

      uni(:,1) = jac(1,:)
      uni(:,2) = jac(2,:)
      uni(:,3) = jac(3,:)

      normal(:) = cross_product( uni(:,2), uni(:,3) )

      det = sum( uni(:,1)*normal(:) )
   end function determinants

   function determinantj( jac ) result(det)
      !-----------------------------------------------------------------------
      ! This function returns the determinant of a 3x2 Jacobian matrix.
      !-----------------------------------------------------------------------
      real(8), intent(in) :: jac(ndim,2)
      real(8) :: normal(ndim)
      real(8) :: det

      normal(:) = cross_product( jac(:,1), jac(:,2) )

      det = norm(normal)
   end function determinantj

   function invert( mat ) result(inv)
      !---------------------------------------------------------------------------
      ! This subroutine inverts a small square matrix onto itself.
      !---------------------------------------------------------------------------
      real(8), intent(in) :: mat(ndim,ndim) !< matrix containing three vectors.
      real(8) :: uni(ndim,ndim)
      real(8) :: inv(ndim,ndim)
      real(8) :: vol

      uni(:,1) = mat(1,:)
      uni(:,2) = mat(2,:)
      uni(:,3) = mat(3,:)

      inv(:,1) = cross_product( uni(:,2), uni(:,3) )
      inv(:,2) = cross_product( uni(:,3), uni(:,1) )
      inv(:,3) = cross_product( uni(:,1), uni(:,2) )

      vol = sum( uni(:,1)*inv(:,1) )

      ! the inverse matrix calculated from the reciprocal vectors.
      inv(:,:) = inv(:,:)/vol
    end function invert

    subroutine samples( nips, points, weights )
      !---------------------------------------------------------------------------
      ! This subroutine returns the local coordinates and weighting coefficients
      ! of the integrating points.
      !---------------------------------------------------------------------------
      integer, intent(in)  :: nips
      real(8), intent(out) :: points(ndim,nips)
      real(8), intent(out) :: weights(nips)
      real(8) :: a, b

      select case(param%nods)
      case(8,20)
         select case(nips)
         case(1)
            points(1:3,1)=0.0
            weights(1)=8.0
         case(8)   
            points(1,1)=-M_SQRT1_3; points(2,1)=-M_SQRT1_3; points(3,1)=-M_SQRT1_3
            points(1,2)= M_SQRT1_3; points(2,2)=-M_SQRT1_3; points(3,2)=-M_SQRT1_3
            points(1,3)= M_SQRT1_3; points(2,3)= M_SQRT1_3; points(3,3)=-M_SQRT1_3
            points(1,4)=-M_SQRT1_3; points(2,4)= M_SQRT1_3; points(3,4)=-M_SQRT1_3
            points(1,5)=-M_SQRT1_3; points(2,5)=-M_SQRT1_3; points(3,5)= M_SQRT1_3
            points(1,6)= M_SQRT1_3; points(2,6)=-M_SQRT1_3; points(3,6)= M_SQRT1_3
            points(1,7)= M_SQRT1_3; points(2,7)= M_SQRT1_3; points(3,7)= M_SQRT1_3
            points(1,8)=-M_SQRT1_3; points(2,8)= M_SQRT1_3; points(3,8)= M_SQRT1_3
            weights(1:8)=1.0
         case(27)
            a=0.774596669241483d0
            points(1, 1)=-a;    points(2, 1)=-a;    points(3, 1)=-a;    weights( 1)=125.d0/729.d0
            points(1, 2)= a;    points(2, 2)=-a;    points(3, 2)=-a;    weights( 2)=125.d0/729.d0
            points(1, 3)= a;    points(2, 3)= a;    points(3, 3)=-a;    weights( 3)=125.d0/729.d0
            points(1, 4)=-a;    points(2, 4)= a;    points(3, 4)=-a;    weights( 4)=125.d0/729.d0
            points(1, 5)=-a;    points(2, 5)=-a;    points(3, 5)= a;    weights( 5)=125.d0/729.d0
            points(1, 6)= a;    points(2, 6)=-a;    points(3, 6)= a;    weights( 6)=125.d0/729.d0
            points(1, 7)= a;    points(2, 7)= a;    points(3, 7)= a;    weights( 7)=125.d0/729.d0
            points(1, 8)=-a;    points(2, 8)= a;    points(3, 8)= a;    weights( 8)=125.d0/729.d0
            points(1, 9)=0.0;   points(2, 9)=-a;    points(3, 9)=-a;    weights( 9)=200.d0/729.d0
            points(1,10)= a;    points(2,10)=0.0;   points(3,10)=-a;    weights(10)=200.d0/729.d0
            points(1,11)=0.0;   points(2,11)= a;    points(3,11)=-a;    weights(11)=200.d0/729.d0
            points(1,12)=-a;    points(2,12)=0.0;   points(3,12)=-a;    weights(12)=200.d0/729.d0
            points(1,13)=-a;    points(2,13)=-a;    points(3,13)=0.0;   weights(13)=200.d0/729.d0
            points(1,14)=a;     points(2,14)=-a;    points(3,14)=0.0;   weights(14)=200.d0/729.d0
            points(1,15)=a;     points(2,15)= a;    points(3,15)=0.0;   weights(15)=200.d0/729.d0
            points(1,16)=-a;    points(2,16)= a;    points(3,16)=0.0;   weights(16)=200.d0/729.d0
            points(1,17)=0.0;   points(2,17)=-a;    points(3,17)= a;    weights(17)=200.d0/729.d0
            points(1,18)=a;     points(2,18)=0.0;   points(3,18)= a;    weights(18)=200.d0/729.d0
            points(1,19)=0.0;   points(2,19)= a;    points(3,19)= a;    weights(19)=200.d0/729.d0
            points(1,20)=-a;    points(2,20)=0.0;   points(3,20)= a;    weights(20)=200.d0/729.d0
            points(1,21)=0.0;   points(2,21)=0.0;   points(3,21)=-a;    weights(21)=320.d0/729.d0
            points(1,22)=0.0;   points(2,22)=-a;    points(3,22)=0.0;   weights(22)=320.d0/729.d0
            points(1,23)=a;     points(2,23)=0.0;   points(3,23)=0.0;   weights(23)=320.d0/729.d0
            points(1,24)=0.0;   points(2,24)= a;    points(3,24)=0.0;   weights(24)=320.d0/729.d0
            points(1,25)=-a;    points(2,25)=0.0;   points(3,25)=0.0;   weights(25)=320.d0/729.d0
            points(1,26)=0.0;   points(2,26)=0.0;   points(3,26)= a;    weights(26)=320.d0/729.d0
            points(1,27)=0.0;   points(2,27)=0.0;   points(3,27)=0.0;   weights(27)=512.d0/729.d0
         case default
            print*,"SAMPLE: wrong number of integrating points for a quadrilateral"
            stop
         end select
      case(4,10)
         select case(nips)
         case(1)
            points(1,1)=0.25d0; points(2,1)=0.25d0; points(3,1)=0.25d0
            weights(1)=1.0d0/6.0d0   
         case(4)
            a=0.585410196624968d0
            b=0.138196601125010d0
            points(1,1)=b; points(2,1)=b; points(3,1)=b
            points(1,2)=a; points(2,2)=b; points(3,2)=b
            points(1,3)=b; points(2,3)=a; points(3,3)=b
            points(1,4)=b; points(2,4)=b; points(3,4)=a
            weights(1:4)=0.25d0/6.0d0
         case default
            print*,"SAMPLE: wrong number of integrating points for a quadrilateral"
            stop
         end select
      case default
         write(*,*) "SAMPLE: wrong element type"
         stop
      end select
   end subroutine samples
   
   subroutine samplej( nipb, pointj, weightj )
      !---------------------------------------------------------------------------
      ! This subroutine returns the local coordinates and weighting coefficients
      ! of the integrating points.
      !---------------------------------------------------------------------------
      integer, intent(in)  :: nipb
      real(8), intent(out) :: pointj(2,nipb)
      real(8), intent(out) :: weightj(nipb)
      real(8) :: r15, wa, w0

      ! Body of samplej
      select case( param%nodj/2 )
      case(4,8) ! Quadrilateral
         select case(nipb)
         case(1)
            pointj(1:2,1)=0.0
            weightj(1)=4.0
         case(4)
            pointj(1,1)=-M_SQRT1_3; pointj(2,1)=-M_SQRT1_3
            pointj(1,2)= M_SQRT1_3; pointj(2,2)=-M_SQRT1_3
            pointj(1,3)= M_SQRT1_3; pointj(2,3)= M_SQRT1_3
            pointj(1,4)=-M_SQRT1_3; pointj(2,4)= M_SQRT1_3
            weightj(1:4)=1.0
         case(9)
            r15=SQRT(15.0d0)/5.0d0
            wa=5.d0/9.d0
            w0=8.d0/9.d0
            pointj(1,1)=-r15;  pointj(2,1)=-r15;  weightj(1)=wa*wa
            pointj(1,2)= r15;  pointj(2,2)=-r15;  weightj(2)=wa*wa
            pointj(1,3)= r15;  pointj(2,3)= r15;  weightj(3)=wa*wa
            pointj(1,4)=-r15;  pointj(2,4)= r15;  weightj(4)=wa*wa
            pointj(1,5)= 0.0;  pointj(2,5)=-r15;  weightj(5)=w0*wa
            pointj(1,6)= r15;  pointj(2,6)= 0.0;  weightj(6)=wa*w0
            pointj(1,7)= 0.0;  pointj(2,7)= r15;  weightj(7)=w0*wa
            pointj(1,8)=-r15;  pointj(2,8)=0.0;   weightj(8)=wa*w0
            pointj(1,9)= 0.0;  pointj(2,9)=0.0;   weightj(9)=w0*w0
         case default
            print*,"SAMPLE: wrong number of integrating points for a quadrilateral"
            stop
         end select
      case(3,6) ! Triangle ! weights divided by 2
         select case(nipb)
         case(1)
            pointj(1,1)=1.0d0/3.0d0; pointj(2,1)=1.0d0/3.0d0
            weightj(1)=1.0d0/2.0d0
         case(3)
            pointj(1,1)=0.5; pointj(2,1)=0.0; weightj(1)=1.0d0/3.0d0/2.0d0
            pointj(1,2)=0.5; pointj(2,2)=0.5; weightj(2)=1.0d0/3.0d0/2.0d0
            pointj(1,3)=0.0; pointj(2,3)=0.5; weightj(3)=1.0d0/3.0d0/2.0d0
         case default
            print*,"SAMPLE: wrong number of integrating points for a triangle"
            stop
         end select
      case default
         print*,"SAMPLEJ: unknown element type"
         stop
      end select
   end subroutine samplej
   
   subroutine shape_funs(funs,points)
      !------------------------------------------------------------------------
      !   This subroutine computes the values of the shape functions
      !   in local coordinates (3D linear element, 8 nodes)
      !------------------------------------------------------------------------
      real(8), intent(out) :: funs(param%nods)
      real(8), intent(in)  :: points(ndim)

      real(8) :: xi, xim, xip
      real(8) :: et, etm, etp
      real(8) :: ze, zem, zep

      xi = points(1)
      et = points(2)
      ze = points(3)

      etm=1.0-et; xim=1.0-xi; zem=1.0-ze
      etp=et+1.0; xip=xi+1.0; zep=ze+1.0
      select case(param%nods)
      case(8)
         funs(1)=0.125d0*xim*etm*zem
         funs(2)=0.125d0*xip*etm*zem
         funs(3)=0.125d0*xip*etp*zem
         funs(4)=0.125d0*xim*etp*zem
         funs(5)=0.125d0*xim*etm*zep
         funs(6)=0.125d0*xip*etm*zep
         funs(7)=0.125d0*xip*etp*zep
         funs(8)=0.125d0*xim*etp*zep
      case(20) ! Need to check again.
         funs( 1)=0.125d0*(1.0-xi)   *(1.0-et)   *(1.0-ze)   *(-xi-et-ze-2.0)
         funs(13)=0.250d0*(1.0-xi)   *(1.0-et)   *(1.0-ze*ze)
         funs( 5)=0.125d0*(1.0-xi)   *(1.0-et)   *(1.0+ze)   *(-xi-et+ze-2.0)
         funs(17)=0.250d0*(1.0-xi*xi)*(1.0-et)   *(1.0+ze)
         funs( 6)=0.125d0*(1.0+xi)   *(1.0-et)   *(1.0+ze)   *(xi-et+ze-2.0)
         funs(14)=0.250d0*(1.0+xi)   *(1.0-et)   *(1.0-ze*ze)
         funs( 2)=0.125d0*(1.0+xi)   *(1.0-et)   *(1.0-ze)   *(xi-et-ze-2.0)
         funs( 9)=0.250d0*(1.0-xi*xi)*(1.0-et)   *(1.0-ze)
         funs(12)=0.250d0*(1.0-xi)   *(1.0-et*et)*(1.0-ze)
         funs(20)=0.250d0*(1.0-xi)   *(1.0-et*et)*(1.0+ze)
         funs(18)=0.250d0*(1.0+xi)   *(1.0-et*et)*(1.0+ze)
         funs(10)=0.250d0*(1.0+xi)   *(1.0-et*et)*(1.0-ze)
         funs( 4)=0.125d0*(1.0-xi)   *(1.0+et)   *(1.0-ze)   *(-xi+et-ze-2.0)
         funs(16)=0.250d0*(1.0-xi)   *(1.0+et)   *(1.0-ze*ze)
         funs( 8)=0.125d0*(1.0-xi)   *(1.0+et)   *(1.0+ze)   *(-xi+et+ze-2.0)
         funs(19)=0.250d0*(1.0-xi*xi)*(1.0+et)   *(1.0+ze)
         funs( 7)=0.125d0*(1.0+xi)   *(1.0+et)   *(1.0+ze)   *(xi+et+ze-2.0)
         funs(15)=0.250d0*(1.0+xi)   *(1.0+et)   *(1.0-ze*ze)
         funs( 3)=0.125d0*(1.0+xi)   *(1.0+et)   *(1.0-ze)   *(xi+et-ze-2.0)
         funs(11)=0.250d0*(1.0-xi*xi)*(1.0+et)   *(1.0-ze)
      case(4)
         funs(2)=xi
         funs(3)=et
         funs(4)=ze
         funs(1)=1.0-xi-et-ze
      case(10)
         funs( 1)=(1.0-xi-et-ze)*(2*(1.0-xi-et-ze)-1.0)
         funs( 2)= xi*(2*xi-1.0)
         funs( 3)= et*(2*et-1.0)
         funs( 4)= ze*(2*ze-1.0)
         funs( 5)=4*(1.0-xi-et-ze)*xi
         funs( 6)=4*xi*et
         funs( 7)=4*et*(1.0-xi-et-ze)
         funs( 8)=4*(1.0-xi-et-ze)*ze
         funs( 9)=4*xi*ze
         funs(10)=4*et*ze
      case default
         write(*,*) "SHAPE_FUN: SHAPE FUNCTION MUST BE SPECIFIED.", param%nods
         stop
      end select
   end subroutine shape_funs
   
   subroutine shape_funj( funj, pointj )
      !------------------------------------------------------------------------
      !   This subroutine computes the values of the shape functions
      !   in local coordinates (2D joint element)
      !------------------------------------------------------------------------
      real(8), intent(out) :: funj(param%nodj/2)
      real(8), intent(in)  :: pointj(2)

      real(8) :: xi, xim, xip
      real(8) :: et, etm, etp
      real(8) :: c3

      xi = pointj(1)
      et = pointj(2)

      select case(param%nodj/2)
      case(4)
         xim = 0.25d0*(1.0-xi); xip=0.25d0*(1.0+xi)
         etm = 0.25d0*(1.0-et); etp=0.25d0*(1.0+et)
         funj(1) = 4*xim*etm
         funj(4) = 4*xim*etp
         funj(3) = 4*xip*etp
         funj(2) = 4*xip*etm
      case(8)
         funj(1) = 0.25d0*(1.0-xi)*(1.0-et)*(-xi-et-1.0)
         funj(8) = 0.50d0*(1.0-xi)*(1.0-et*et)
         funj(4) = 0.25d0*(1.0-xi)*(1.0+et)*(-xi+et-1.0)
         funj(7) = 0.50d0*(1.0-xi*xi)*(1.0+et)
         funj(3) = 0.25d0*(1.0+xi)*(1.0+et)*(xi+et-1.0)
         funj(6) = 0.50d0*(1.0+xi)*(1.0-et*et)
         funj(2) = 0.25d0*(1.0+xi)*(1.0-et)*(xi-et-1.0)
         funj(5) = 0.50d0*(1.0-xi*xi)*(1.0-et)
      case(3)
         funj(1) = 1.0-xi-et
         funj(2) = xi
         funj(3) = et
      case(6)
         c3 = 1.0-xi-et
         funj(1) = (2*c3-1.0)*c3
         funj(2) = (2*xi-1.0)*xi 
         funj(3) = (2*et-1.0)*et
         funj(4) = 4*c3*xi
         funj(5) = 4*xi*et
         funj(6) = 4*et*c3
      case default
         write(*,*) "SHAPE_FUNJ: SHAPE FUNCTION MUST BE SPECIFIED.",param%nodj/2
         stop    
      end select
   end subroutine shape_funj

   subroutine shape_ders( ders, points )
      !---------------------------------------------------------------------------
      !   This subroutine produces derivatives of shape functions withe respect
      !   to local coordinates.
      !---------------------------------------------------------------------------
      real(8), intent(out) :: ders(param%nods,ndim)
      real(8), intent(in)  :: points(ndim)

      integer :: l
      real(8) :: xi, xim, xip
      real(8) :: et, etm, etp
      real(8) :: ze, zem, zep
      real(8) :: xii(20),eti(20),zei(20),xi0,et0,ze0,c4

      xi = points(1)
      et = points(2)
      ze = points(3)
      
      !--- Body of shape_ders
      xim=1.0-xi; xip=xi+1.0 
      etm=1.0-et; etp=et+1.0 
      zem=1.0-ze; zep=ze+1.0
      
      select case(param%nods)
      case(4)
         ders(1,1)=-1.0; ders(1,2)=-1.0; ders(1,3)=-1.0
         ders(2,1)= 1.0; ders(2,2)= 0.0; ders(2,3)= 0.0
         ders(3,1)= 0.0; ders(3,2)= 1.0; ders(3,3)= 0.0
         ders(4,1)= 0.0; ders(4,2)= 0.0; ders(4,3)= 1.0
      case(10)
         c4=1.0-xi-et-ze
         ders( 1,1)=-4*c4+1.0 ; ders( 1,2)=-4*c4+1.0 ; ders( 1,3)=-4*c4+1.0
         ders( 2,1)= 4*xi-1.0 ; ders( 2,2)= 0.0      ; ders( 2,3)= 0.0
         ders( 3,1)= 0.0      ; ders( 3,2)= 4*et-1.0 ; ders( 3,3)= 0.0
         ders( 4,1)= 0.0      ; ders( 4,2)= 0.0      ; ders( 4,3)= 4*ze-1.0
         ders( 5,1)= 4*(c4-xi); ders( 5,2)=-4*xi     ; ders( 5,3)=-4*xi
         ders( 6,1)= 4*et     ; ders( 6,2)= 4*xi     ; ders( 6,3)= 0.0
         ders( 7,1)=-4*et     ; ders( 7,2)= 4*(c4-et); ders( 7,3)=-4*et
         ders( 8,1)=-4*ze     ; ders( 8,2)=-4*ze     ; ders( 8,3)= 4*(c4-ze)
         ders( 9,1)= 4*ze     ; ders( 9,2)= 0.0      ; ders( 9,3)= 4*xi
         ders(10,1)= 0.0      ; ders(10,2)= 4*ze     ; ders(10,3)= 4*et
      case(8)
         ders(1,1)=-0.125d0*etm*zem; ders(1,2)=-0.125d0*xim*zem; ders(1,3)=-0.125d0*xim*etm 
         ders(2,1)=-0.125d0*etm*zep; ders(2,2)=-0.125d0*xim*zep; ders(2,3)= 0.125d0*xim*etm
         ders(3,1)= 0.125d0*etm*zep; ders(3,2)=-0.125d0*xip*zep; ders(3,3)= 0.125d0*xip*etm   
         ders(4,1)= 0.125d0*etm*zem; ders(4,2)=-0.125d0*xip*zem; ders(4,3)=-0.125d0*xip*etm
         ders(5,1)=-0.125d0*etp*zem; ders(5,2)= 0.125d0*xim*zem; ders(5,3)=-0.125d0*xim*etp   
         ders(6,1)=-0.125d0*etp*zep; ders(6,2)= 0.125d0*xim*zep; ders(6,3)= 0.125d0*xim*etp
         ders(7,1)= 0.125d0*etp*zep; ders(7,2)= 0.125d0*xip*zep; ders(7,3)= 0.125d0*xip*etp    
         ders(8,1)= 0.125d0*etp*zem; ders(8,2)= 0.125d0*xip*zem; ders(8,3)=-0.125d0*xip*etp
      case(20)
         xii= (/-1, 1, 1,-1,-1, 1, 1,-1, 0, 1, 0,-1,-1, 1,1,-1, 0,1,0,-1/)
         eti= (/-1,-1, 1, 1,-1,-1, 1, 1,-1, 0, 1, 0,-1,-1,1, 1,-1,0,1, 0/)
         zei= (/-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0,0, 0, 1,1,1, 1/)
         do l=1, 20
            xi0=xi*xii(l)
            et0=et*eti(l)
            ze0=ze*zei(l)
            select case(l)
            case(  9, 11, 17, 19 )
               ders(l,1)=-0.50d0*xi*(1.0+et0)*(1.0+ze0)
               ders(l,2)= 0.25d0*eti(l)*(1.0-xi*xi)*(1.0+ze0)
               ders(l,3)= 0.25d0*zei(l)*(1.0-xi*xi)*(1.0+et0)
            case( 12, 10, 18, 20 )
               ders(l,1)= 0.25d0*xii(l)*(1.0-et*et)*(1.0+ze0)
               ders(l,2)=-0.50d0*et*(1.0+xi0)*(1.0+ze0)
               ders(l,3)= 0.25d0*zei(l)*(1.0+xi0)*(1.0-et*et)
            case( 13, 14, 15, 16 )
               ders(l,1)= 0.25d0*xii(l)*(1.0+et0)*(1.0-ze*ze)
               ders(l,2)= 0.25d0*eti(l)*(1.0+xi0)*(1.0-ze*ze)
               ders(l,3)=-0.50d0*ze*(1.0+xi0)*(1.0+et0)
            case default
               ders(l,1)=0.125d0*xii(l)*(1.0+et0)*(1.0+ze0) &
                  * (2*xi0+et0+ze0-1.0)
               ders(l,2)=0.125d0*eti(l)*(1.0+xi0)*(1.0+ze0) &
                  * (xi0+2*et0+ze0-1.0)
               ders(l,3)=0.125d0*zei(l)*(1.0+xi0)*(1.0+et0) &
                  * (xi0+et0+2*ze0-1.0)
            end select
         end do
      end select
   end subroutine shape_ders
   
   subroutine shape_derj( derj, pointj )
      !---------------------------------------------------------------------------
      !   This subroutine produces derivatives of shape functions with respect
      !   to local coordinates.
      !---------------------------------------------------------------------------
      real(8),intent(out) :: derj(param%nodj/2,2)
      real(8),intent(in)  :: pointj(2)

      real(8) :: xi, xim, xip
      real(8) :: et, etm, etp
      real(8) :: c3

      xi = pointj(1)
      et = pointj(2)

      !--- Body of shape_der
      xim = 0.25d0*(1.0-xi)
      xip = 0.25d0*(1.0+xi)
      etm = 0.25d0*(1.0-et)
      etp = 0.25d0*(1.0+et)

      select case( param%nodj/2 )
      case(3)
         derj(1,1) =-1.0; derj(1,2) = -1.0
         derj(2,1) = 1.0; derj(2,2) =  0.0
         derj(3,1) = 0.0; derj(3,2) =  1.0
      case(6)
         c3=1.0-xi-et
         derj(1,1) = -4*c3+1.0  ; derj(1,2) = -4*c3+1.0
         derj(2,1) =  4*xi-1.0  ; derj(2,2) =  0.0
         derj(3,1) =  0.0       ; derj(3,2) =  4*et-1.0
         derj(4,1) = -4*xi+4*c3 ; derj(4,2) = -4*xi
         derj(5,1) =  4*et      ; derj(5,2) =  4*xi
         derj(6,1) = -4*et      ; derj(6,2) =  4*c3-4*et
      case(4)                                                              
         derj(1,1) = -etm; derj(1,2)=-xim
         derj(2,1) =  etm; derj(2,2)=-xip
         derj(3,1) =  etp; derj(3,2)= xip
         derj(4,1) = -etp; derj(4,2)= xim
      case(8)
         derj(1,1) =  etm*(2*xi+et); derj(1,2) =  xim*(xi+2*et)
         derj(2,1) =  etm*(2*xi-et); derj(2,2) =  xip*(2*et-xi)
         derj(3,1) =  etp*(2*xi+et); derj(3,2) =  xip*(xi+2*et)
         derj(4,1) =  etp*(2*xi-et); derj(4,2) =  xim*(2*et-xi)
         derj(5,1) = -4*etm*xi     ; derj(5,2) = -8*xim*xip 
         derj(6,1) =  8*etp*etm    ; derj(6,2) = -4*xip*et
         derj(7,1) = -4*etp*xi     ; derj(7,2) =  8*xim*xip
         derj(8,1) = -8*etm*etp    ; derj(8,2) = -4*xim*et   
      case default
         write(*,*) "WARNING: wrong number of nodes in shape_der"
         stop
      end select
   end subroutine shape_derj

   subroutine ellmatj( Lmat, jac )
      !---------------------------------------------------------------------------
      ! Computes normal vectors from Jacobian
      !---------------------------------------------------------------------------
      real(8), intent(out) :: Lmat(ndim,ndim) ! basis vectors at the gauss point on the surface.
      real(8), intent(in)  :: jac(ndim,2)
      real(8) :: det
      real(8), parameter :: tiny=1.0d-10
      
      ! gradient of global coordinates. two tangent vectors on the surface.
      Lmat(:,1:2) = jac(:,1:2)

      ! normal vector on the surface.
      Lmat(:,3) = cross_product( Lmat(:,1), Lmat(:,2) )

      ! Compute Jacobian and unit normal vector
      det = norm(Lmat(:,3))
      if( det < tiny ) then
         print *,"WARNING: small jacobian in ellmatj", det
         stop
      end if

      ! normalization.
      Lmat(:,1) = Lmat(:,1)/norm(Lmat(:,1))
      Lmat(:,2) = Lmat(:,2)/norm(Lmat(:,2))
      Lmat(:,3) = Lmat(:,3)/norm(Lmat(:,3))

      ! another second tangent vector, perpendicular to the first one.
      Lmat(:,2) = cross_product( Lmat(:,3), Lmat(:,1) )
   end subroutine ellmatj

   subroutine beemats( Bmat, gfun )
      !---------------------------------------------------------------------------
      ! This subroutine forms the bee matrix in 3-d
      !---------------------------------------------------------------------------
      real(8), intent(out) :: Bmat(nst,ndim*param%nods)
      real(8), intent(in)  :: gfun(param%nods,ndim)

      integer :: k, l, m, n
!!$      real :: x,y,z ! MIZUHO fixed
      real(8) :: x,y,z ! MIZUHO fixed

      Bmat(:,:) = 0.0d0
      select case(nst)
      case(3,4)
         do m=1, param%nods
            l=2*(m-1)+1
            k=2*(m-1)+2
            x=gfun(m,1)
            y=gfun(m,2)
            Bmat(1,l)=x
            Bmat(3,l)=y
            Bmat(2,k)=y
            Bmat(3,k)=x
         end do
      case(6)
         do m=1, param%nods
            l=3*(m-1)+1
            k=3*(m-1)+2
            n=3*(m-1)+3
            x=gfun(m,1)
            y=gfun(m,2)
            z=gfun(m,3)
            Bmat(1,l)=x
            Bmat(4,l)=y
            Bmat(6,l)=z
            Bmat(4,k)=x
            Bmat(2,k)=y
            Bmat(5,k)=z
            Bmat(6,n)=x
            Bmat(5,n)=y
            Bmat(3,n)=z
         end do
      case default
         write(*,*) "wrong dimension for nst in bee matrix"
      end select
   end subroutine beemats

   subroutine beematj( Bmatj, Lmat, funj )
      real(8), intent(out) :: Bmatj(ndim,ndim*param%nodj)
      real(8), intent(in)  :: Lmat(ndim,ndim), funj(param%nodj/2)
      integer::i,j,k
      real(8) :: nmat(ndim,ndim*param%nodj)
      
      nmat=0.0d0

      do j=1, param%nodj/2
         do i=1, ndim
            k=3*(j-1)+i    
            nmat(i,k) = funj(j)
         end do
      end do

      do j=param%nodj/2+1, param%nodj
         do i=1, ndim
            k=3*(j-1)+i
            nmat(i,k) = -funj(j-param%nodj/2)
         end do
      end do

      Bmatj = matmul(transpose(Lmat),nmat)
   end subroutine beematj

   subroutine deemats( Dmat, e, v )
      !---------------------------------------------------------------------------
      ! This subroutine returns the elastic dee matrix for plane strain
      !---------------------------------------------------------------------------
      real(8), intent(out) :: Dmat(nst,nst)
      real(8), intent(in)  :: e, v

      real(8) :: lam, d11, d12, d44

      lam = e/((1.0+v)*(1.0-2*v))
      d11 = lam*(1.0-v)
      d12 = lam*v
      d44 = lam*(1.0-2*v)*0.5
      
      Dmat(:,:) = 0.0
      Dmat(1,1) = d11
      Dmat(2,2) = d11
      Dmat(3,3) = d11
      Dmat(1,2) = d12
      Dmat(1,3) = d12
      Dmat(2,3) = d12
      Dmat(2,1) = d12
      Dmat(3,1) = d12
      Dmat(3,2) = d12
      Dmat(4,4) = d44
      Dmat(5,5) = d44
      Dmat(6,6) = d44
    end subroutine deemats

    subroutine deeinvs( dinv, e, v )
       real(8), intent(out) :: dinv(nst,nst)
       real(8), intent(in)  :: e, v

       real(8) :: d11, d12, d44

       d11 = 1.0/e
       d12 = -v/e
!!$       d44 = 1*(1+v)/e
       d44 = 2*(1+v)/e ! MIZUHO fixed

       dinv(:,:) = 0.0
       dinv(1,1) = d11
       dinv(2,2) = d11
       dinv(3,3) = d11
       dinv(1,2) = d12
       dinv(1,3) = d12
       dinv(2,3) = d12
       dinv(2,1) = d12
       dinv(3,1) = d12
       dinv(3,2) = d12
       dinv(4,4) = d44
       dinv(5,5) = d44
       dinv(6,6) = d44
    end subroutine deeinvs

    subroutine deematj( Dmatj, k )
       !---------------------------------------------------------------------------
       ! This subroutine returns the elastic matrix for joint element
       !---------------------------------------------------------------------------
       real(8),intent(out) :: Dmatj(ndim,ndim)
       real(8),intent(in)  :: k(3)

       ! Body of deemat
       Dmatj(:,:) = 0.0d0
       Dmatj(1,1) = k(1)
       Dmatj(2,2) = k(2)
       Dmatj(3,3) = k(3)
    end subroutine deematj

   subroutine tensor( sig, sigv )
      real(8), intent(out) :: sig(ndim,ndim)
      real(8), intent(in)  :: sigv(nst)
      
      sig(1,1) = sigv(1)
      sig(2,2) = sigv(2)
      sig(3,3) = sigv(3)
      sig(1,2) = sigv(4)
      sig(2,1) = sigv(4)
      sig(2,3) = sigv(5)
      sig(3,2) = sigv(5)
      sig(1,3) = sigv(6)
      sig(3,1) = sigv(6)
   end subroutine tensor

   subroutine get_principal_stress( eig, sig )
      real(8), intent(out) :: eig(ndim)
      real(8), intent(in)  :: sig(nst)

      real(8) :: mat(ndim,ndim), work(ndim*3), info
      external dsyev

      call tensor(mat,sig)

      call dsyev( "N", "U", ndim, mat, ndim, eig, work, ndim*3, info )
   end subroutine get_principal_stress

   !> 
   ! read filenames for input and output from the config file.
   subroutine loadConfig( filename_cfg )
      character(len=*), intent(in) :: filename_cfg
      integer :: iunit, status, i
      character(len=1024) :: buf, key

      open( newunit=iunit, file=filename_cfg, iostat=status, status="old" )
      if( status/=0 ) then
         write(*,*) "Error: can not open config file: ", trim(filename_cfg)
         call exit(1)
      end if

      do
         read(iunit,'(a)',iostat=status) buf
         if( status/=0 ) exit ! end of file.

         if( buf(1:1)=="#" ) cycle ! skip comment line

         i = index(buf,"#")
         if( i>0 ) buf = buf(1:i-1) ! cut comment words.

         if( len_trim(buf)==0 ) cycle ! skip empty line.

         i = index(buf,"=")
         if( i>0 ) buf(i:i) = " " ! replace "=" by a white space.

         read(buf,*,iostat=status) key ! read the first word as keyword.
         if( status/=0 ) goto 1000

         select case( key )
         case("solver")
            read(buf,*,iostat=status) key, solver
            if( status/=0 ) goto 1000
         case("filename_gid")
            read(buf,*,iostat=status) key, filename_gid
            if( status/=0 ) goto 1000
         case("filename_dat")
            read(buf,*,iostat=status) key, filename_dat
            if( status/=0 ) goto 1000
         case("filename_err")
            read(buf,*,iostat=status) key, filename_err
            if( status/=0 ) goto 1000
         case("filename_presss_ini")
            read(buf,*,iostat=status) key, filename_presss_ini
            if( status/=0 ) goto 1000
         case("filename_pressj_ini")
            read(buf,*,iostat=status) key, filename_pressj_ini
            if( status/=0 ) goto 1000
         case("filename_point_load")
            read(buf,*,iostat=status) key, filename_point_load
            if( status/=0 ) goto 1000
         case("filename_surfj_load")
            read(buf,*,iostat=status) key, filename_surfj_load
            if( status/=0 ) goto 1000
         case("filename_presss_load")
            read(buf,*,iostat=status) key, filename_presss_load
            if( status/=0 ) goto 1000
         case("filename_pressj_load")
            read(buf,*,iostat=status) key, filename_pressj_load
            if( status/=0 ) goto 1000
         case("filename_temp_load")
            read(buf,*,iostat=status) key, filename_temp_load
            if( status/=0 ) goto 1000
         case("filename_den_dif")
            read(buf,*,iostat=status) key, filename_den_dif
            if( status/=0 ) goto 1000
         case("filename_log")
            read(buf,*,iostat=status) key, filename_log
            if( status/=0 ) goto 1000
         case("filename_res")
            read(buf,*,iostat=status) key, filename_res
            if( status/=0 ) goto 1000
         case("filename_postmesh")
            read(buf,*,iostat=status) key, filename_postmesh
            if( status/=0 ) goto 1000
         case("filename_postres")
            read(buf,*,iostat=status) key, filename_postres
            if( status/=0 ) goto 1000
         case default
            write(*,*) "Error: unknown key in config file: ", trim(filename_cfg)
            write(*,*) "Error: unknown key: ", trim(key)
            goto 1000
         end select
      end do

      close(iunit)

      write(*,*) "***************************"
      write(*,*) "configurations"
      write(*,*) " solver:               ", trim(solver)
      write(*,*) " filename_gid:         ", trim(filename_gid)
      write(*,*) " filename_dat:         ", trim(filename_dat)
      write(*,*) " filename_err:         ", trim(filename_err)
      write(*,*) " filename_presss_ini:  ", trim(filename_presss_ini)
      write(*,*) " filename_pressj_ini:  ", trim(filename_pressj_ini)
      write(*,*) " filename_point_load:  ", trim(filename_point_load)
      write(*,*) " filename_surfj_load:  ", trim(filename_surfj_load)
      write(*,*) " filename_presss_load: ", trim(filename_presss_load)
      write(*,*) " filename_pressj_load: ", trim(filename_pressj_load)
      write(*,*) " filename_temp_load:   ", trim(filename_temp_load)
      write(*,*) " filename_den_dif:     ", trim(filename_den_dif)
      write(*,*) " filename_log:         ", trim(filename_log)
      write(*,*) " filename_res:         ", trim(filename_res)
      write(*,*) " filename_postmesh:    ", trim(filename_postmesh)
      write(*,*) " filename_postres:     ", trim(filename_postres)
      write(*,*) "***************************"
      write(*,*)
      
      if( solver == "PARDISO" ) then
         upper_only = .true.
#ifndef ENABLE_PARDISO
         write(*,*) "Error: PARDISO is disabled in this program."
         call exit(1)
#endif
      else if( solver == "PETSc" ) then
         upper_only = .false.
#ifndef ENABLE_PETSC
         write(*,*) "Error: PETSC is disabled in this program."
         call exit(1)
#endif
      else
         write(*,*) "Error: invalid solver: ", trim(solver)
         call exit(1)
      end if

      return

1000  continue
      close(iunit)
      write(*,*) "Error: broken line in config file: ", trim(filename_cfg)
      write(*,*) "Error: broken line: ", trim(buf)
      call exit(1)
   end subroutine loadConfig

   !> 
   ! read FEM and initial and time-series fields from the GiD file.
   subroutine loadGiD
      integer :: nhexahedra, ntetrahedra, nquadrilateral, ntriangle
      integer :: quadratic_type
      integer :: i, k, n, iel, ip, il
      integer :: j
      real(8) :: grad, gradv(nst)
      real(8) :: nodfix(ndim), value(ndim)
      integer :: status
      integer,allocatable :: gidelemnum(:),solidorjoint(:)

      open( newunit=iunit_gid, file=filename_gid, iostat=status, status="old" )
      if( status/=0 ) then
         write(*,*) "Error: can not open input GiD file: ", trim(filename_gid)
         call exit(1)
      end if

      open( newunit=ounit_err, file=filename_err, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create err file: ", trim(filename_err)
         call exit(1)
      end if
      !
      !--- Read data from GiD output data ---
      !
      read(iunit_gid,*) 
      read(iunit_gid,*)
      read(iunit_gid,'(35x,i10)') param%neles
      read(iunit_gid,'(35x,i10)') param%nelej
      read(iunit_gid,'(35x,i10)') nhexahedra
      read(iunit_gid,'(35x,i10)') ntetrahedra
      read(iunit_gid,'(35x,i10)') nquadrilateral
      read(iunit_gid,'(35x,i10)') ntriangle
      read(iunit_gid,'(35x,i10)') param%nnode
      read(iunit_gid,'(35x,i10)') quadratic_type

      if( nhexahedra>0 ) then ! hexahedral elements
         if( quadratic_type==0 ) then ! 8-node hexahedra
            param%nods =  8 ! linear hexahedra element for solid.
            param%nodb =  4 ! linear quadrilateral element for joint.
            param%nodj =  8 ! linear quadrilateral element for joint.
         else if( quadratic_type==1 ) then ! 20-node hexahedra
            param%nods = 20 ! quadratic hexahedra element for solid.
            param%nodb =  8 ! quadratic quadrilateral element for joint.
            param%nodj = 16 ! quadratic quadrilateral element for joint.
         else
            print *, "ERROR: unknown element type"
            stop
         end if
      else ! tetrahedral elements
         if( quadratic_type==0 ) then ! 4-node tetrahedra
            param%nods =  4 ! linear tetrahedra element for solid.
            param%nodb =  3 ! linear triangular element for joint.
            param%nodj =  6 ! linear triangular element for joint.
         else if( quadratic_type==1 ) then ! 10-node tetrahedra
            param%nods = 10 ! quadratic tetrahedra element for solid.
            param%nodb =  6 ! quadratic triangular element for joint.
            param%nodj = 12 ! quadratic triangular element for joint.
         else
            print *, "ERROR: unknown element type"
            stop
         end if
      end if

      read(iunit_gid,*)
      read(iunit_gid,'(50x,i10)')    param%ntime
      read(iunit_gid,'(50x,f25.15)') param%convergence
      read(iunit_gid,'(50x,i10)')    param%niter
      read(iunit_gid,'(50x,i10)')    param%nips

      if( param%nips==0 ) then ! Set the default value
         if( nhexahedra>0 ) then
            param%nips = 8
         else
            param%nips = 4
         end if
      end if

      read(iunit_gid,'(50x,i10)') param%nipb
      if( param%nipb==0 ) then ! Set the default value
         if( nhexahedra>0 ) then
            param%nipb = 4
         else
            param%nipb = 3
         end if
      end if

      read(iunit_gid,'(56x,i10)') param%nips_stress
      read(iunit_gid,'(56x,i10)') param%nipj_stress
      read(iunit_gid,*)
      read(iunit_gid,'(31x,i10)') param%use_stress_grad
      read(iunit_gid,'(31x,i10)') param%use_press_grad ! 0: absolute value, 1: gradient
      read(iunit_gid,'(29x,i10)') param%use_temp_grad
      read(iunit_gid,'(29x,f25.15)') param%gravity

      allocate( geom%coord(ndim,param%nnode) )

      read(iunit_gid,*)
      read(iunit_gid,*)
      do n=1, param%nnode
         read(iunit_gid,*) k, geom%coord(:,k)
      end do

      allocate( geom%displ(3,param%nnode), geom%displ_inc(3,param%nnode) )
      geom%displ(:,:) = 0.0d0 ! initial displacement.

      if( k==param%nnode ) then
         write(ounit_err,*) "* THE COORDINATES OF ALL NODES WERE READ PROPERLY."
      else
         write(ounit_err,*) "ERROR: READING NODAL COORDINATES."
         stop
      end if

      allocate( geom%enodes(param%nods,param%neles), geom%etypes(param%neles) )
      allocate( gidelemnum(param%neles+param%nelej), solidorjoint(param%neles+param%nelej) )
      
      read(iunit_gid,*)
      read(iunit_gid,*)
      read(iunit_gid,*)
      do iel=1, param%neles
         read(iunit_gid,*) k, geom%enodes(:,iel), geom%etypes(iel)
         gidelemnum(k)=iel
         solidorjoint(k)=1
      end do
      ! GiD numbering may start from joint elements
      !if( k==param%neles ) then
      !   write(ounit_err,*) "* THE NODAL CONNECTIVITIES WERE READ PROPERLY."
      !else
      !   write(ounit_err,*) "ERROR: READING NODAL CONNECTIVITY."
      !   stop
      !end if

      read(iunit_gid,*)
      read(iunit_gid,'(30x,i10)') param%nprops
      read(iunit_gid,*)
      allocate( param%props(param%nprops+4))
      do ip=1, param%nprops
         read(iunit_gid,*) k, param%props(k)
      end do

      allocate( geom%enodej(param%nodj,param%nelej), geom%etypej(param%nelej) )

      read(iunit_gid,*)
      read(iunit_gid,*)
      !allocate( param%propj(param%npropj+2))
      do iel=1, param%nelej
         read(iunit_gid,*) k, geom%enodej(1:param%nodb,iel), geom%etypej(iel)
         gidelemnum(k)=iel
         solidorjoint(k)=-1
      end do

      read(iunit_gid,*)
      read(iunit_gid,'(31x,i10)') param%npropj
      allocate( param%propj(param%npropj+2)) ! Furui fixed
      read(iunit_gid,*)
      do ip=1, param%npropj
         read(iunit_gid,*) k, param%propj(k)
      end do

      ! Compute center of each solid element.
      allocate( geom%ecenters(3,param%neles) )
      do iel=1, param%neles
         geom%ecenters(:,iel) = 0.0d0
         do i=1, param%nods
            geom%ecenters(:,iel) = geom%ecenters(:,iel) + geom%coord(:,geom%enodes(i,iel))
         end do
         geom%ecenters(:,iel) = geom%ecenters(:,iel)/dble(param%nods)
      end do

      ! Compute center of each joint element.
      allocate( geom%ecenterj(3,param%nelej) )
      do iel=1, param%nelej
         geom%ecenterj(:,iel) = 0.0d0
         do i=1, param%nodb
            geom%ecenterj(:,iel) = geom%ecenterj(:,iel) + geom%coord(:,geom%enodej(i,iel))
         end do
         geom%ecenterj(:,iel) = geom%ecenterj(:,iel)/dble(param%nodb)
      end do


      read(iunit_gid,*)
      read(iunit_gid,*)
      read(iunit_gid,'(24x,i10)') n
      read(iunit_gid,*)

      allocate( geom%nodfree(ndim,param%nnode), geom%valfix(ndim,param%nnode) ) !, geom%codfree(ndim*param%nnode)
      
      geom%nodfree(:,:) = 1
      geom%valfix(:,:)  = 0.0
      do i=1, n
         !read(iunit_gid,*) k, geom%nodfree(:,k), geom%valfix(:,k)
          read(iunit_gid,*) k, nodfix(:), value(:)
          do j=1,ndim
            if(nodfix(j)==0)then
                geom%nodfree(j,k)=0
                geom%valfix(j,k)=value(j)
            endif
          enddo  
      end do

      read(iunit_gid,*)
      read(iunit_gid,'(46x,i10)') n
      read(iunit_gid,*)
      if( n>0 ) then
         allocate( geom%surfj_elem(n), geom%surfj_enode(param%nodb,n), surfj_load(ndim,param%nodb,n) )
         do iel=1, n
            read(iunit_gid,*) geom%surfj_elem(gidelemnum(iel)), geom%surfj_enode(:,gidelemnum(iel)), surfj_load(:,1,gidelemnum(iel))
            do i=1, param%nodb
               surfj_load(:,i,gidelemnum(iel)) = surfj_load(:,1,gidelemnum(iel))
            end do
         end do
      end if

      allocate( stresss_ini(nst,param%neles) )
      stresss_ini=0! Default values
      read(iunit_gid,*)
      read(iunit_gid,*)
      read(iunit_gid,'(49x,i10)') n
      read(iunit_gid,*)
      if( n>0 ) then
         do k=1, n
            read(iunit_gid,*) iel, stresss_ini(:,gidelemnum(iel)), gradv(:)
            if( param%use_stress_grad==1 ) then
               stresss_ini(:,gidelemnum(iel)) = stresss_ini(:,gidelemnum(iel)) + gradv(:)*geom%ecenters(3,gidelemnum(iel))
            end if
         end do
      end if

      allocate( presss_ini(param%neles) )
      presss_ini=0 ! default values
      read(iunit_gid,*)
      read(iunit_gid,'(55x,i10)') n
      read(iunit_gid,*)
      if( n>0 ) then
         do k=1, n
            read(iunit_gid,*) iel, presss_ini(gidelemnum(iel)), grad
            if( param%use_press_grad==1 ) then
               presss_ini(gidelemnum(iel)) = presss_ini(gidelemnum(iel)) + grad*geom%ecenters(3,gidelemnum(iel))
            end if
         end do
      end if

      allocate( temp_ini(param%neles) )
      temp_ini=0 ! default values
      read(iunit_gid,*)
      read(iunit_gid,'(52x,i10)') n
      read(iunit_gid,*)
      if( n>0 ) then
         do k=1, n
            read(iunit_gid,*) iel, temp_ini(gidelemnum(iel)), grad
            if( param%use_temp_grad==1 ) then
               temp_ini(gidelemnum(iel)) = temp_ini(gidelemnum(iel)) + grad*geom%ecenters(3,gidelemnum(iel))
            end if
         end do
      end if

      read(iunit_gid,*)
      read(iunit_gid,'(49x,i10)') n
      read(iunit_gid,*) 
      if( param%nelej>0 ) then
         allocate( stressj_ini(nst,param%nelej) )
         stressj_ini=0 ! default value
         if( n>0 ) then
            do k=1, n
               read(iunit_gid,*) iel, stressj_ini(:,gidelemnum(iel)), gradv(:)
               if( param%use_stress_grad==1 ) then
                  stressj_ini(:,gidelemnum(iel)) = stressj_ini(:,gidelemnum(iel)) &
                     + gradv(:)*geom%ecenterj(3,gidelemnum(iel))
               end if
            end do
         end if
      end if

      read(iunit_gid,*)
      read(iunit_gid,'(54x,i10)') n
      read(iunit_gid,*)
      if( param%nelej>0 ) then
         allocate( pressj_ini(param%nelej) )
         pressj_ini=0 ! default value
         if( n>0 ) then
            do k=1, n
               read(iunit_gid,*) iel, pressj_ini(gidelemnum(iel)), grad
               if( param%use_press_grad==1 ) then
                  pressj_ini(gidelemnum(iel)) = pressj_ini(gidelemnum(iel)) &
                     + grad*geom%ecenterj(3,gidelemnum(iel))
               end if
            end do
         end if
      end if
      read(iunit_gid,*)
      read(iunit_gid,*)
      read(iunit_gid,'(37x,i10)') n ! number_of_pointload

!!!!! point load is not read.

      allocate( presss(param%neles) )
      presss=0 ! default value
      read(iunit_gid,*)
      read(iunit_gid,'(54x,i10)') n
      read(iunit_gid,*)
      if( n>0 ) then
         do k=1, n
            read(iunit_gid,*) iel, presss(gidelemnum(iel)), grad
            if( param%use_press_grad==1 ) then
               presss(gidelemnum(iel)) = presss(gidelemnum(iel)) + grad*geom%ecenters(3,gidelemnum(iel))
            end if
         end do
      end if

      allocate( temp(param%neles) )
      temp=0 ! default value
      read(iunit_gid,*)
      read(iunit_gid,'(52x,i10)') n
      read(iunit_gid,*)
      if( n>0 ) then
         do k=1, n
            read(iunit_gid,*) iel, temp(gidelemnum(iel)), grad
            if( param%use_temp_grad==1 ) then
               temp(gidelemnum(iel)) = temp(gidelemnum(iel)) + grad*geom%ecenters(3,gidelemnum(iel))
            end if
         end do
      end if

      allocate( density_dif(param%neles) )
      density_dif=0 ! default value
      read(iunit_gid,*)
      read(iunit_gid,'(48x,i10)') n
      read(iunit_gid,*)
      if( n>0 ) then
         do k=1, n
            read(iunit_gid,*) iel, density_dif(gidelemnum(iel))
         end do
      end if

      read(iunit_gid,*)
      read(iunit_gid,'(54x,i10)') n
      read(iunit_gid,*)
      if( param%nelej>0 ) then
         allocate( pressj(param%nelej) )
         pressj=0 ! default value

         if( n>0 ) then
            do k=1, n
               read(iunit_gid,*) iel, pressj(gidelemnum(iel)), grad
               if( param%use_press_grad==1 ) then
                  pressj(gidelemnum(iel)) = pressj(gidelemnum(iel)) &
                     + grad*geom%ecenterj(3,gidelemnum(iel))
               end if
            end do
         end if
      end if

      read(iunit_gid,*)
      read(iunit_gid,'(19x,i10)') param%nlayer
      allocate( param%layer_name(param%nlayer) )
      allocate( param%layer_nele(param%nlayer) )
      allocate( geom%elem_layer_names(param%neles) )
      allocate( geom%elem_layer_namej(param%nelej) )

      allocate( param%layer_enode(param%neles+param%nelej,param%nlayer) )

      do il=1, param%nlayer
         read(iunit_gid,'(13x,a20)') param%layer_name(il)
         read(iunit_gid,'(31x,i10)') param%layer_nele(il)
         read(iunit_gid,*) param%layer_enode(1:param%layer_nele(il),il)
         
         !do k=1,param%layer_nele(il)
         !  param%layer_enode(1:param%layer_nele(il),il)=gidelemnum(param%layer_enode(k,il))
         !enddo

         do i=1, param%layer_nele(il)
            iel = param%layer_enode(i,il)
            if(solidorjoint(iel) == 1 ) then ! solid element.
               geom%elem_layer_names(gidelemnum(iel)) = param%layer_name(il)
            else ! joint element.
               !iel=iel-param%neles ! iel as joint element.
               geom%elem_layer_namej(gidelemnum(iel)) = param%layer_name(il)
            end if
         end do
      end do

      close(iunit_gid)

!!$      call reconstructLayer
      call fixJointFace

   contains

      subroutine fixJointFace
         integer :: ielj, iels
         logical, allocatable :: joint(:)

         allocate( joint(param%nnode) ); joint(:) = .false.

         do ielj=1, param%nelej
            do iels=1, param%neles
               if( isJointOnSolid( geom%enodej(:,ielj), geom%enodes(:,iels) ) ) then
                  if( .not. isValidOrder( geom%enodej(:,ielj), geom%enodes(:,iels) ) ) then
                     call revertOrder( geom%enodej(:,ielj) )
                  end if
                  exit
               end if
            end do

            call addAnotherFace( joint, geom%enodej(:,ielj) )
         end do

         write(ounit_err,'("# of joint nodes =",i10)') count(joint)
         deallocate( joint )
      end subroutine fixJointFace

      logical function isJointOnSolid( enodej, enodes )
         integer, intent(in) :: enodej(param%nodj)
         integer, intent(in) :: enodes(param%nods)

         select case(param%nodb)
         case(3,6)
            isJointOnSolid = isContained( enodej(1:3), enodes(1:4) )
         case(4,8)
            isJointOnSolid = isContained( enodej(1:4), enodes(1:8) )
         end select
      end function isJointOnSolid

      logical function isContained( a, b )
         integer, intent(in) :: a(:), b(:)
         integer :: i

         do i=1, size(a)
            if( .not. any( a(i)==b(:) ) ) then
               isContained = .false.
               return ! all elements are not contained.
            end if
         end do
         isContained = .true. ! all elements are contained.
      end function isContained

      logical function isValidOrder( enodej, enodes )
         integer, intent(in) :: enodej(param%nodj)
         integer, intent(in) :: enodes(param%nods)

         real(8) :: pos(3,4), vec(3,3), nor(3)
         integer :: j

         pos(:,1) = geom%coord(:,enodej(1))
         pos(:,2) = geom%coord(:,enodej(2))
         pos(:,3) = geom%coord(:,enodej(3))

         select case(param%nodb)
         case(3,6)
            j = aValueNotContained( enodej(1:3), enodes(1:4) )
         case(4,8)
            j = aValueNotContained( enodej(1:4), enodes(1:8) )
         end select
         !pos(:,4) = geom%coord(:,j) ! Fixed on July 9, 2026
         pos(:,4) = geom%coord(:,enodes(j))

         vec(:,1) = pos(:,2) - pos(:,1)
         vec(:,2) = pos(:,3) - pos(:,1)
         vec(:,3) = pos(:,4) - pos(:,1)
         nor(:) = cross_product( vec(:,1), vec(:,2) )

         !isValidOrder = dot_product( vec(:,3), nor(:) ) < 0.0 Mizuho (fixed on July 9, 2026)
         isValidOrder = dot_product( vec(:,3), nor(:) ) > 0.0
      end function isValidOrder

      function aValueNotContained( a, b ) result(j)
         integer, intent(in) :: a(:), b(:)
         integer :: j

         do j=1, size(b)
            if( .not. any( b(j)==a(:) ) ) return ! found.
         end do

         j=0 ! not found.
      end function aValueNotContained


      subroutine revertOrder( enodej )
         integer, intent(inout) :: enodej(param%nodj)
         integer :: reverse(param%nodb)

         write(ounit_err,'(10x,"WRONG SURFACE")')

         select case(param%nodb)
         case(3) ! linear triangular element for joint.
            reverse(1) = enodej(1)
            reverse(2) = enodej(3)
            reverse(3) = enodej(2)
         case(6) ! quadratic triangular element for joint.
            reverse(1) = enodej(1)
            reverse(2) = enodej(3)
            reverse(3) = enodej(2)
            reverse(4) = enodej(6)
            reverse(5) = enodej(5)
            reverse(6) = enodej(4)
         case(4) ! linear quadrilateral element for joint.
            reverse(1) = enodej(1)
            reverse(2) = enodej(4)
            reverse(3) = enodej(3)
            reverse(4) = enodej(2)
         case(8) ! quadratic quadrilateral element for joint.
            reverse(1) = enodej(1)
            reverse(2) = enodej(4)
            reverse(3) = enodej(3)
            reverse(4) = enodej(2)
            reverse(5) = enodej(8)
            reverse(6) = enodej(7)
            reverse(7) = enodej(6)
            reverse(8) = enodej(5)
         end select

         enodej(1:param%nodb) = reverse(1:param%nodb)
      end subroutine revertOrder

      subroutine addAnotherFace( joint, enodej )
         logical, intent(inout) :: joint(param%nnode)
         integer, intent(inout) :: enodej(param%nodj)
         integer :: i, j, n

         do j=1, param%nodb
            n = enodej(j)
            i = findAnotherNode( n )
            if( i>0 ) then
               joint(n) = .true.
               enodej(j+param%nodb) = i
            else
               enodej(j+param%nodb) = n
            end if
         end do
      end subroutine addAnotherFace

      function findAnotherNode( j ) result(i)
         integer, intent(in) :: j
         integer :: i

         do i=1, param%nnode
            if( i==j ) cycle ! skip identical node.

            ! check difference of coordinates.
            if( all( abs(geom%coord(:,i) - geom%coord(:,j))<1.0d-10 ) ) return ! found.
         end do

         i=0 ! not found.
      end function findAnotherNode
   end subroutine loadGiD

   !> 
   ! write FEM and initial and time-series fields to some DAT files.
   subroutine saveDAT
      integer :: i, n, iel, ip, il, et, nfix
      integer :: itime
      integer :: status

      open( newunit=ounit_dat, file=filename_dat, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create dat file: ", trim(filename_dat)
         call exit(1)
      end if

      open( newunit=ounit_presss_ini, file=filename_presss_ini, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create initial solid pressure file: ", trim(filename_presss_ini)
         call exit(1)
      end if

      open( newunit=ounit_pressj_ini, file=filename_pressj_ini, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create initial joint pressure file: ", trim(filename_pressj_ini)
         call exit(1)
      end if

      open( newunit=ounit_point_load, file=filename_point_load, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create point load file: ", trim(filename_point_load)
         call exit(1)
      end if

      open( newunit=ounit_surfj_load, file=filename_surfj_load, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create surface load file: ", trim(filename_surfj_load)
         call exit(1)
      end if

      open( newunit=ounit_presss_load, file=filename_presss_load, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create solid pressure file: ", trim(filename_presss_load)
         call exit(1)
      end if

      open( newunit=ounit_pressj_load, file=filename_pressj_load, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create joint pressure file: ", trim(filename_pressj_load)
         call exit(1)
      end if

      open( newunit=ounit_temp_load, file=filename_temp_load, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create temperature file: ", trim(filename_temp_load)
         call exit(1)
      end if

      open( newunit=ounit_den_dif, file=filename_den_dif, iostat=status, status="replace" )
      if( status /= 0 ) then
         write(*,*) "Error: can not create density diff file: ", trim(filename_den_dif)
         call exit(1)
      end if

      !--- CREATE INPUT FILE FOR WFAST ---
      write(ounit_dat,"(a)") "*(1) nod,nodb,nodj,nels,nelsj,nn,nip,nipb,np_types,np_typesj,nstr_solid,nstr_joint,grav"
      write(ounit_dat,10) param%nods, param%nodb, param%nodj, param%neles, param%nelej, param%nnode, &
         param%nips, param%nipb, param%nprops, param%npropj, param%nips_stress, param%nipj_stress, param%gravity
10    format(12(1x,i7),1x,e15.7)    
      write(ounit_dat,"(a)") "*(2) # of time steps, convergence limit for nonlinear calculations, max number of iterations"
      write(ounit_dat,20) param%ntime, param%convergence, param%niter
20    format(1x,i7,e11.3,i7)
      write(ounit_dat,"(a)") "*(3) Material properties, e, v, rho, em, vm, at, fric_umf,co_umf]"
      do ip=1, param%nprops
         write(ounit_dat,30) param%props(ip)
      end do
30    format(1x,8e16.8)

      if( param%nelej>0 ) then
         write(ounit_dat,*) "*(4) ks1,ks2,kn,gs,gn,multksn,tnsile_strength,fric,co,dsc,fric_dynamic"
         do ip=1, param%npropj
            write(ounit_dat,40) param%propj(ip)
         end do
      end if
40    format(1x,13e16.8) 

      write(ounit_dat,*) "*(5) Layer Information, NumLayer,LayerName/#_of_Elems"
      write(ounit_dat,'(1x,i6)') param%nlayer
      do il=1, param%nlayer
         write(ounit_dat,'(4x,a20,1x,i8)') param%layer_name(il), param%layer_nele(il)
      end do

      write(ounit_dat,*) "*(7) Global coordinates [i,g_coord(:,i),in=1,nn]"
      do n=1, param%nnode
         write(ounit_dat,70) n, geom%coord(:,n)
      end do
70    format(1x,i8,3e16.8)

      write(ounit_dat,*) "*(8) Global node numbers for solid elements,etype,layername"
      do iel=1, param%neles
         select case(param%nods)
         case(4)
            write(ounit_dat,80) iel, geom%enodes(:,iel), geom%etypes(iel), geom%elem_layer_names(iel)
!!$            write(ounit_dat,80) iel, geom%enodes(:,iel), geom%etypes(iel), param%layer_names(et)
         case(10)
            write(ounit_dat,81) iel, geom%enodes(:,iel), geom%etypes(iel), geom%elem_layer_names(iel)
         case(8)
            write(ounit_dat,82) iel, geom%enodes(:,iel), geom%etypes(iel), geom%elem_layer_names(iel)
         case(20)
            write(ounit_dat,83) iel, geom%enodes(:,iel), geom%etypes(iel), geom%elem_layer_names(iel)
         end select
      end do
80    format(1x,i8,2x,4i8,1x,i3,1x,a20) 
81    format(1x,i8,2x,10i8,1x,i3,1x,a20) 
82    format(1x,i8,2x,8i8,1x,i3,1x,a20) 
83    format(1x,i8,2x,20i8,1x,i3,1x,a20) 

      if( param%nelej>0 ) then
         write(ounit_dat,*) "*(9) Global node numbers for joint elements,etype,layername"
         do iel=1, param%nelej
            select case(param%nodj)
            case(8)
               write(ounit_dat,90) iel, geom%enodej(:,iel), geom%etypej(iel), geom%elem_layer_namej(iel)
!!$               write(ounit_dat,90) iel, geom%enodej(:,iel), geom%etypej(iel), param%layer_namej(et)
            case(16)
               write(ounit_dat,91) iel, geom%enodej(:,iel), geom%etypej(iel), geom%elem_layer_namej(iel)
            case(6)
               write(ounit_dat,92) iel, geom%enodej(:,iel), geom%etypej(iel), geom%elem_layer_namej(iel)
            case(12)
               write(ounit_dat,93) iel, geom%enodej(:,iel), geom%etypej(iel), geom%elem_layer_namej(iel)
            end select
         end do
      end if
90    format(1x,i8,2x,8i8,1x,i3,1x,a20)
91    format(1x,i8,2x,16i8,1x,i3,1x,a20)
92    format(1x,i8,2x,6i8,1x,i3,1x,a20)
93    format(1x,i8,2x,12i8,1x,i3,1x,a20)

      write(ounit_dat,*) "*(10) Node freedom data [nr/k,nf(:,k),i=1,nr]"
      nfix=0
      do i=1, param%nnode
         if( any( geom%nodfree(:,i)==0 ) ) nfix=nfix+1
      end do
      write(ounit_dat,100) nfix
100   format(3x,i8)

      do n=1, param%nnode
         if( any( geom%nodfree(:,n)==0 ) ) write(ounit_dat,101) n, geom%nodfree(:,n), geom%valfix(:,n)
      end do
      write(ounit_dat,*) "*(11) Initial stress, pore pressure and temperature for solid elements"
101   format(1x,i8,2x,3i2,3e15.7)



      do iel=1, param%neles
         write(ounit_dat,110) iel, stresss_ini(:,iel), presss_ini(iel), temp_ini(iel)
      end do
110   format(1x,i8,8e16.8)

      if( param%nelej>0 ) then
         write(ounit_dat,*) "*(12) Initial tractions and pore pressure for joint elements"
         do iel=1, param%nelej
            write(ounit_dat,120) iel, stressj_ini(:,iel), pressj_ini(iel)
         end do
      end if
120   format(1x,i8,7e16.8)

      write(ounit_dat,*) "!!! Load data is provided in separate files."

      do itime=1, param%ntime
         ! Load step data: Point loading
         write(ounit_point_load, &
            '("*** Load step=",i3,", Point loads data [number of loaded_nodes/k,loads(nf(:,k)),i=1,loaded_nodes]")') itime
         write(ounit_point_load,'(7x,"0")')

!!!!! point load is not writen.

         ! Load step data: Surface loading
         write(ounit_surfj_load,'("*** Load step=",i3,", Loaded surface data [number of loaded surfaces] ")') itime
         write(ounit_surfj_load,140) size(geom%surfj_elem)
140      format(1x,i8)

         if( size(geom%surfj_elem)>0 ) then
            do iel=1, size(geom%surfj_elem)
               select case(param%nods)
               case(4)
                  write(ounit_surfj_load,141) geom%surfj_elem(iel), param%nodb, geom%surfj_enode(:,iel)
                  write(ounit_surfj_load,142) surfj_load(:,:,iel)
               case(10)
                  write(ounit_surfj_load,143) geom%surfj_elem(iel), param%nodb, geom%surfj_enode(:,iel)
                  write(ounit_surfj_load,144) surfj_load(:,:,iel)
               case(8)
                  write(ounit_surfj_load,145) geom%surfj_elem(iel), param%nodb, geom%surfj_enode(:,iel)
                  write(ounit_surfj_load,146) surfj_load(:,:,iel)
               case(20)
                  write(ounit_surfj_load,147) geom%surfj_elem(iel), param%nodb, geom%surfj_enode(:,iel)
                  write(ounit_surfj_load,148) surfj_load(:,:,iel)
               end select
            end do
         end if
141      format(1x,i8,3x,i2,3x,4i8)
142      format(3x,3(3e14.6,2x))

143      format(1x,i8,3x,i2,3x,10i8)
144      format(3x,6(3e14.6,2x))   

145      format(1x,i8,3x,i2,3x,8i8)
146      format(3x,4(3e14.6,2x))   

147      format(1x,i8,3x,i2,3x,20i8)
148      format(3x,8(3e14.6,2x))   

         ! Load step data: Pore pressure (solid elements)
         write(ounit_presss_load,'("*** Load step=",i3,", Pore pressure for solid elements")') itime

         if( itime==1 ) write(ounit_presss_ini,*) " IELEM     X     Y     Z     Pore_Pressure"
         do iel=1, param%neles
            write(ounit_presss_load,150) iel, presss(iel)
150         format(1x,i8,e16.8)    

            if( itime==1 ) then
               write(ounit_presss_ini,162) iel, geom%ecenters(:,iel), presss(iel)
            end if
         end do
162      format(1x,i8,3e16.8,3x,e16.8)

         ! Load step data: Temperature (solid elements)
         write(ounit_temp_load,'("*** Load step=",i3,", Temperature for solid elements")') itime
         do iel=1, param%neles
            write(ounit_temp_load,150) iel, temp(iel)
         end do

         ! Load step data: Density change (solid elements)
         write(ounit_den_dif,'("*** Load step=",i3,", Density change for solid elements")') itime
         do iel=1, param%neles
            write(ounit_den_dif,150) iel, density_dif(iel)
         end do

         ! Load step data: Pore pressure (joint elements)
         if( param%nelej>0 ) then
            write(ounit_pressj_load,'("*** Load step=",i3,", Pore pressure for joint elements")') itime
            if( itime==1 ) write(ounit_pressj_ini,*) " IELEM     X     Y     Z     Po"
            do iel=1, param%nelej
               write(ounit_pressj_load,160) iel, pressj(iel)
160            format(1x,i8,e16.8)    

               if( itime==1 ) then
                  write(ounit_pressj_ini,162) iel, geom%ecenterj(:,iel), pressj(iel)
               end if
            end do
         end if
      end do

      write(ounit_err,*) "wfastinp: Created input files successfully."

      close(ounit_dat)
      close(ounit_err)
      close(ounit_presss_ini)
      close(ounit_pressj_ini)
      close(ounit_point_load)
      close(ounit_surfj_load)
      close(ounit_presss_load)
      close(ounit_pressj_load)
      close(ounit_temp_load)
      close(ounit_den_dif)
   end subroutine saveDAT

   !> 
   ! read FEM and initial fields to the DAT file.
   subroutine loadDAT
      integer :: status
      integer :: i, j, k, l, n, it, iel, ip, ei, ej, et, et_new
      integer, parameter :: cmax = 400
      integer, allocatable :: ja_work(:,:)
      integer :: nnz_ja
      type(propertyj), pointer :: prop

      open( newunit=iunit_dat, file=filename_dat, iostat=status, status="old")
      if( status/=0 ) then
         write(*,*) "Error: can not open input dat file: ", trim(filename_dat)
         call exit(1)
      end if
      
      !
      !--- (1) Read input data and complete the nodal freedom matrix ---
      !
      ! (1) nods,nodb,nodj,neles,nele,nnode,nips,nipb,nprops,npropj,nips_stress,nipj_stress,gravity
      read(iunit_dat,*)
      read(iunit_dat,*) param%nods, param%nodb, param%nodj, param%neles, param%nelej, param%nnode, param%nips, param%nipb, &
         param%nprops, param%npropj, param%nips_stress, param%nipj_stress, param%gravity

      ! (2) # of time steps, convergence limit for nonlinear calculations, max number of iterations
      read(iunit_dat,*)
      read(iunit_dat,*) param%ntime, param%convergence, param%niter

      ! (3) Material properties, e, v, rho, em, vm, at, fric_umf,co_umf]
      read(iunit_dat,*)
      allocate( param%props(param%nprops))
      do it=1, param%nprops
         read(iunit_dat,*) param%props(it)
      end do

      ! (4) ks1,ks2,kn,gs,gn,multksn,tnsile_strength,fric,co,dsc,fric_dynamic
      if( param%nelej>0 ) then
         read(iunit_dat,*)
         allocate( param%propj(param%npropj))
         do it=1, param%npropj
            read(iunit_dat,*) param%propj(it)
         end do
      end if

      ! (5) Layer Information, NumLayer,LayerName/#_of_Elems
      read(iunit_dat,*)
      read(iunit_dat,*) param%nlayer
      allocate( param%layer_name(param%nlayer) )

      do it=1, param%nlayer
         read(iunit_dat,*) param%layer_name(it)
      end do

      ! (7) Global coordinates [i,g_coord(:,i),in=1,nnode]
      allocate( geom%coord(ndim,param%nnode) )
      read(iunit_dat,*)
      do n=1, param%nnode
         read(iunit_dat,*) k, geom%coord(1:3,k)
      end do
      
      allocate( geom%displ(3,param%nnode), geom%displ_inc(3,param%nnode) )
      geom%displ(:,:) = 0.0d0 ! initial displacement.

      ! (8) Global node numbers for solid elements, etype, layer_name
      allocate( geom%enodes(param%nods,param%neles), geom%etypes(param%neles), geom%elem_layer_names(param%neles) )
      read(iunit_dat,*)

      do iel=1, param%neles
         read(iunit_dat,*) k, geom%enodes(:,k), geom%etypes(k), geom%elem_layer_names(k)
      end do

      ! (9) Global node numbers for joint elements, etype, layer_name
      if( param%nelej>0 ) then
         allocate( geom%enodej(param%nodj,param%nelej), geom%etypej(param%nelej), geom%elem_layer_namej(param%nelej) )
         read(iunit_dat,*)

         do iel=1, param%nelej
            read(iunit_dat,*) k, geom%enodej(:,k), geom%etypej(k), geom%elem_layer_namej(k)
         end do
         
      end if

      ! Compute center of each solid element.
      allocate( geom%ecenters(3,param%neles) )
      do iel=1, param%neles
         geom%ecenters(:,iel) = 0.0d0
         do i=1, param%nods
            geom%ecenters(:,iel) = geom%ecenters(:,iel) + geom%coord(:,geom%enodes(i,iel))
         end do
         geom%ecenters(:,iel) = geom%ecenters(:,iel)/dble(param%nods)
      end do

      ! Compute center of each joint element.
      allocate( geom%ecenterj(3,param%nelej) )
      do iel=1, param%nelej
         geom%ecenterj(:,iel) = 0.0d0
         do i=1, param%nodb
            geom%ecenterj(:,iel) = geom%ecenterj(:,iel) + geom%coord(:,geom%enodej(i,iel))
         end do
         geom%ecenterj(:,iel) = geom%ecenterj(:,iel)/dble(param%nodb)
      end do

      
      ! (10) Node freedom data [nr/k,nf(:,k),i=1,nr]
      allocate( geom%nodfree(ndim,param%nnode), geom%valfix(ndim,param%nnode) )! , geom%codfree(ndim*param%nnode)
      geom%nodfree(:,:) = 1
      geom%valfix(:,:)  = 0.0d0
      read(iunit_dat,*)
      read(iunit_dat,*) n
      do i=1, n
         read(iunit_dat,*) k, geom%nodfree(1:3,k), geom%valfix(1:3,k)
      end do
      geom%codfree(1:3*param%nnode) => geom%nodfree

      !
      !--- (2) Obtain global steering vectors ---
      !
      allocate( geom%esteers(ndim*param%nods,param%neles) )
      do iel=1, param%neles ! for solid elements
         do i=1, param%nods
            n=geom%enodes(i,iel)

            geom%esteers( 3*(i-1)+1, iel ) = 3*(n-1)+1
            geom%esteers( 3*(i-1)+2, iel ) = 3*(n-1)+2
            geom%esteers( 3*(i-1)+3, iel ) = 3*(n-1)+3
         end do
      end do

      if( param%nelej>0 ) then
         allocate( geom%esteerj(ndim*param%nodj,param%nelej) )
         do iel=1, param%nelej ! for joint elements
            do i=1, param%nodj
               n=geom%enodej(i,iel)
               
               geom%esteerj( 3*(i-1)+1, iel ) = 3*(n-1)+1
               geom%esteerj( 3*(i-1)+2, iel ) = 3*(n-1)+2
               geom%esteerj( 3*(i-1)+3, iel ) = 3*(n-1)+3
            end do
         end do
      end if

      !--- Initiliaze the internal solver memory pointer. This is only
      ! necessary for the FIRST call of the PARDISO solver.
      if( solver == "PARDISO" ) then
         upper_only = .true.
#ifndef ENABLE_PARDISO
         write(*,*) "Error: PARDISO is disabled in this program."
         call exit(1)
#endif
      end if

      if( solver == "PETSc" ) then
         upper_only = .false.
#ifndef ENABLE_PETSC
         write(*,*) "Error: PETSC is disabled in this program."
         call exit(1)
#endif
      end if
      
      !
      !--- (3) Determine parameters required for Pardiso ---
      !
      allocate( ja_work(cmax,ndim*param%nnode) )
      ja_work(:,:) = 0 ! 0 means no entry.

      do iel=1, param%neles           ! all elements.
         do ei=1, ndim*param%nods         ! all coordinates in this elment.
            i=geom%esteers(ei,iel)     ! row of this coord in the matrix.
            do ej=1, ndim*param%nods      ! all coordinates in this elment.
               j=geom%esteers(ej,iel)  ! column of this coord in the matrix.

               if( upper_only .and. i>j ) cycle ! skip lower triangle part of th matrix.

               if( any( ja_work(:,i) == j ) ) cycle ! already added in the i's list.

               call push_back( ja_work(:,i), j ) ! add coord j at the end of the i's list.
            end do ! ej
         end do ! ei
      end do ! iel

      do iel=1, param%nelej          ! all elements.
         do ei=1, ndim*param%nodj        ! all coordinates in this elment.
            i=geom%esteerj(ei,iel)    ! row of this coord in the matrix.
            do ej=1, ndim*param%nodj     ! all coordinates in this elment.
               j=geom%esteerj(ej,iel) ! column of this coord in the matrix.

               if( upper_only .and. i>j ) cycle ! skip lower triangle part of th matrix.

               if( any( ja_work(:,i) == j ) ) cycle ! already added in the i's list.

               call push_back( ja_work(:,i), j ) ! add coord j at the end of the i's list.
            end do ! ej
         end do ! ei
      end do ! iel

      do i=1, ndim*param%nnode
         call sort( ja_work(:,i) )
      end do ! i

      geom%nnz = count( ja_work(:,:) /= 0 )

      allocate( geom%ia(ndim*param%nnode+1), geom%ja(geom%nnz) )

      l=0
      do i=1, ndim*param%nnode
         geom%ia(i) = l+1

         nnz_ja = count( ja_work(:,i) /= 0 )
         geom%ja(l+1:l+nnz_ja) = ja_work(1:nnz_ja,i)

         l=l+nnz_ja
      end do ! i
      geom%ia(ndim*param%nnode+1) = geom%nnz+1

      deallocate(ja_work)

#ifdef ENABLE_PETSC
      if( solver == "PETSc" ) then
         allocate( geom%irn(geom%nnz), geom%jcn(geom%nnz), geom%idx(ndim*param%nnode) ) ! PETSc

         do i=1, ndim*param%nnode
            do l=geom%ia(i), geom%ia(i+1)-1
               j = geom%ja(l)

               geom%irn(l) = i-1
               geom%jcn(l) = j-1
            end do
            geom%idx(i) = i-1
         end do
      end if
#endif

      allocate( stiffmat(geom%nnz) )

      allocate( geom%eindexs(ndim*param%nods,ndim*param%nods,param%neles) )

      geom%eindexs(:,:,:)=0
      do iel=1, param%neles           ! all elements.
         do ei=1, ndim*param%nods         ! all coordinates in this elment.
            i=geom%esteers(ei,iel)     ! row of this coord in the matrix.
            do ej=1, ndim*param%nods      ! all coordinates in this elment.
               j=geom%esteers(ej,iel)  ! column of this coord in the matrix.

               if( upper_only .and. i>j ) cycle ! skip lower triangle part of th matrix.

               do l=geom%ia(i), geom%ia(i+1)-1
                  if( geom%ja(l) == j ) then
                     geom%eindexs(ei,ej,iel) = l
                     exit
                  end if
               end do

            end do ! ej
         end do ! ei
      end do ! iel

      if( param%nelej>0 ) then
         allocate( geom%eindexj(ndim*param%nodj,ndim*param%nodj,param%nelej) )
         do iel=1, param%nelej          ! all elements.
            do ei=1, ndim*param%nodj        ! all coordinates in this elment.
               i=geom%esteerj(ei,iel)    ! row of this coord in the matrix.
               do ej=1, ndim*param%nodj     ! all coordinates in this elment.
                  j=geom%esteerj(ej,iel) ! column of this coord in the matrix.

                  if( upper_only .and. i>j ) cycle ! skip lower triangle part of th matrix.

                  do l=geom%ia(i), geom%ia(i+1)-1
                     if( geom%ja(l) == j ) then
                        geom%eindexj(ei,ej,iel) = l
                        exit
                     end if
                  end do
               end do ! ej
            end do ! ei
         end do ! iel
      end if

      ! Open input and output files
     
      open(newunit=iunit_point_load, file=filename_point_load, iostat=status, status="old")
      if( status /= 0 ) then
         write(*,*) "Error: can not open point load file: ", trim(filename_point_load)
         call exit(1)
      end if

      open(newunit=iunit_surfj_load, file=filename_surfj_load, iostat=status, status="old")
      if( status /= 0 ) then
         write(*,*) "Error: can not open surface load file: ", trim(filename_surfj_load)
         call exit(1)
      end if

      open(newunit=iunit_presss_load, file=filename_presss_load, iostat=status, status="old")
      if( status /= 0 ) then
         write(*,*) "Error: can not open pore load file: ", trim(filename_presss_load)
         call exit(1)
      end if

      open(newunit=iunit_pressj_load, file=filename_pressj_load, iostat=status, status="old")
      if( status /= 0 ) then
         write(*,*) "Error: can not open joint pore load file: ", trim(filename_pressj_load)
         call exit(1)
      end if

      open(newunit=iunit_temp_load, file=filename_temp_load, iostat=status, status="old")
      if( status /= 0 ) then
         write(*,*) "Error: can not open temperature file: ", trim(filename_temp_load)
         call exit(1)
      end if
      
      open(newunit=iunit_den_dif, file=filename_den_dif, iostat=status, status="old")
      if( status /= 0 ) then
         write(*,*) "Error: can not open density difference file: ", trim(filename_den_dif)
         call exit(1)
      end if
      
      open(newunit=ounit_log, file=filename_log, iostat=status, status="replace")
      if( status /= 0 ) then
         write(*,*) "Error: can not create log file: ", trim(filename_log)
         call exit(1)
      end if

      open(newunit=ounit_res, file=filename_res, iostat=status, status="replace")
      if( status /= 0 ) then
         write(*,*) "Error: can not create result file: ", trim(filename_res)
         call exit(1)
      end if

      open(newunit=ounit_postmesh, file=filename_postmesh, iostat=status, status="replace")
      if( status /= 0 ) then
         write(*,*) "Error: can not create post mesh file: ", trim(filename_postmesh)
         call exit(1)
      end if

      open(newunit=ounit_postres, file=filename_postres, iostat=status, status="replace")
      if( status /= 0 ) then
         write(*,*) "Error: can not create post result file: ", trim(filename_postres)
         call exit(1)
      end if

      ! (11) Initial stress, pore pressure and temperature for solid elements
      allocate( stresss_ini(nst,param%neles) )
      allocate( presss_ini(param%neles), presss(param%neles), temp_ini(param%neles) )

      read(iunit_dat,*)
      do iel=1, param%neles
         read(iunit_dat,*) i, stresss_ini(:,i), presss_ini(i), temp_ini(i)
      end do

      presss(:) = 0.0d0 ! it will be printed prior to being set.

      ! (12) Initial tractions and pore pressure for joint elements
      if( param%nelej>0 ) then
         allocate( stressj_ini(nst,param%nelej) )
         allocate( pressj_ini(param%nelej), pressj(param%nelej) )

         read(iunit_dat,*)
         do iel=1, param%nelej
            read(iunit_dat,*) i, stressj_ini(:,i), pressj_ini(i)
         end do

         pressj(:) = 0.0d0 ! it will be printed prior to being set.
      end if

      close(iunit_dat)

      allocate( temp(param%neles) )
      temp=0 ! default value

      allocate( density_dif(param%neles) )
      density_dif=0 ! default value

      allocate( strains_ini(nst,param%neles) )

      if( param%nelej>0 ) then
         allocate( cohesionj(param%nipb,param%nelej), tensilej(param%nipb,param%nelej), stiffj(ndim,param%nipb,param%nelej) )
      end if

      if( param%nelej>0 ) then
         allocate( tractionj_ini(ndim,param%nipb,param%nelej), tractionj(ndim,param%nipb,param%nelej) )
      end if

      if( param%nelej>0 ) then
         allocate( coulombj_ini(param%nipb,param%nelej), coulombj(param%nipb,param%nelej) )
         allocate( coulombj_dif(param%nipb,param%nelej) )
         coulombj_dif(:,:) = 0.0d0 ! it will be printed prior to being set.
      end if

      if( param%nelej>0 ) then
         allocate( failurej(param%nelej) )
      end if

      allocate( energy(param%nips_stress,param%neles) )
      
      allocate( robust(param%nips_stress,param%neles), robust_ini(param%nips_stress,param%neles) )
      allocate( robust_dif(param%nips_stress,param%neles) )

      allocate( load(ndim,param%nnode) )

      if( param%nelej>0 ) then
         allocate( seismic_mom(param%nelej) )
         seismic_mom(:) = 0.0d0 ! it will be printed prior to being set.
         allocate( seismic_int(param%nelej) )
         seismic_int(:) = 0.0d0 ! it will be printed prior to being set.
      end if

      !--- (6) Set joint stiffness and cohesion for joint elements ---
      !
      if( param%nelej>0 ) then
         do iel=1, param%nelej
            prop => param%propj(geom%etypej(iel))

            do ip=1, param%nipb
               stiffj(1:3,ip,iel) = prop%ksn(1:3)
               tensilej(ip,iel)   = prop%T
               cohesionj(ip,iel)  = prop%So
            end do
         end do
      end if

   contains
      function find( list, value ) result(index)
         character(len=20), intent(in) :: list(:)
         character(len=20), intent(in) :: value
         integer :: index
         
         do index=1, size(list)
            if( list(index) == value ) return
         end do
         index=0
      end function find
      
      subroutine push_back( array, value )
         integer, intent(inout) :: array(:)
         integer, intent(in) :: value
         integer :: len, i

         len = size(array)
         do i=1, len
            if( array(i)==0 ) then
               array(i) = value
               return
            end if
         end do

         write(*,*) "TOO MANY NEIGHBOR NODES. INCREASE cmax ", len
         stop
      end subroutine push_back

      subroutine sort( array )
         integer, intent(inout) :: array(:)
         integer :: len, i, j, work

         len = count( array(:) /= 0 )
         do i=1, len-1
            do j=i+1, len
               if( array(i) > array(j) ) then
                  work     = array(i)
                  array(i) = array(j)
                  array(j) = work
               end if
            end do
         end do
      end subroutine sort
   end subroutine loadDAT

   !> 
   ! calculate initial strain from initial stress and pressure.
   subroutine calcInitStrain
      integer :: iel
      type(propertys), pointer :: prop
      real(8) :: dinv(nst,nst) ! inverse of D matrix (nst,nst)
      real(8) :: sig_eff(nst), eps_eff(nst), eps_neu(nst)

      !
      !--- (5) Compute initial strain for solid elements ---
      !
      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         call deeinvs( dinv, prop%E, prop%nu ) ! inverse of D matrix.
         
         ! Compute initial effective stress
         sig_eff(1:3) = stresss_ini(1:3,iel) + presss_ini(iel)
         sig_eff(4:6) = stresss_ini(4:6,iel)

         ! Compute initial effective strain
         eps_eff(1:6) = matmul(dinv,sig_eff)

         ! Compute initial neutral strain
         eps_neu(1:3) = -(1.0d0-2.0d0*prop%num)/prop%Em * presss_ini(iel) ! MIZUHO fixed
!!$         eps_neu(1:3) = -(1.0-1*prop%num)/prop%Em * presss_ini(iel)
         eps_neu(4:6) = 0.0d0

         ! Compute initial strain
         strains_ini(:,iel) = eps_eff(:) + eps_neu(:)
      end do
   end subroutine calcInitStrain

   !> 
   ! calculate initial traction from initial pressure.
   subroutine calcInitTraction
      integer :: iel, ip
      real(8) :: pointj(2,param%nipb), weightj(param%nipb) ! integrating point local coordinates for joints
      real(8) :: funj(param%nodb,param%nipb), derj(param%nodj/2,2,param%nipb) ! shape function
      type(propertyj), pointer :: prop
      real(8) :: coordj(ndim,param%nodj/2) ! element nodal coordinates (ndim,param%nodj/2)
      real(8) :: sig(ndim,ndim), sig_rot(ndim,ndim)
      real(8) :: jac(ndim,2)   ! Jacobian matrix (ndim,ndim)
      real(8) :: Lmat(ndim,ndim)  ! normal vectors on the surface.

      if( param%nelej==0 ) return

      !
      !--- (7) Compute tractions for joint elements ---
      !
      call samplej( param%nipb, pointj, weightj )
      do ip=1, param%nipb
         call shape_funj( funj(:,ip),   pointj(:,ip) )
         call shape_derj( derj(:,:,ip), pointj(:,ip) )
      end do

      do iel=1, param%nelej
         prop => param%propj(geom%etypej(iel))
         
         coordj(:,:) = geom%coord(:,geom%enodej(1:param%nodj/2,iel))

         call tensor( sig, stressj_ini(:,iel) )

         do ip=1, param%nipb
            jac(:,:) = matmul(coordj,derj(:,:,ip))
            call ellmatj( Lmat, jac ) ! normal vectors on the surface.

            ! Calculate stress tensor in the joint coordinates
            sig_rot = matmul(transpose(Lmat),matmul(sig,Lmat))

            ! tractions in the joint coordinates
            tractionj(1:2,ip,iel) = sig_rot(1:2,3)
            tractionj(  3,ip,iel) = sig_rot(  3,3) + pressj_ini(iel)
         end do
      end do

      tractionj_ini(:,:,:) = tractionj(:,:,:)
   end subroutine calcInitTraction

   !> 
   ! calculate initial coulomb stress from initial traction.
   subroutine evalInitCoulomb
      integer :: iel, ip
      type(propertyj), pointer :: prop
      real(8) :: ts, tn
      
      do iel=1, param%nelej
         prop => param%propj(geom%etypej(iel))

         do ip=1, param%nipb
            ts = parallel(      tractionj_ini(:,ip,iel) )
            tn = perpendicular( tractionj_ini(:,ip,iel) )

            coulombj(ip,iel) = ts + prop%mu*tn
         end do
      end do
            
      coulombj_ini(:,:) = coulombj(:,:)
   end subroutine evalInitCoulomb

   !> 
   ! evaluate initial joint failure from initial traction.
   subroutine evalInitFail
      integer :: iel, ip
      type(propertyj), pointer :: prop
      !
      !---(8) Check the joint stability
      !
      if( param%nelej==0 ) return

      failurej(:) = .false.

      do iel=1, param%nelej
         prop => param%propj(geom%etypej(iel))

         do ip=1, param%nipb
            if( tensilej(ip,iel) <= tractionj_ini(3,ip,iel) ) then
               tensilej(ip,iel) = tractionj_ini(3,ip,iel) + prop%dsc
               failurej(iel) = .true.
            end if

            if( cohesionj(ip,iel) < coulombj_ini(ip,iel) ) then
               cohesionj(ip,iel) = coulombj_ini(ip,iel) + prop%dsc
               failurej(iel) = .true.
            end if
         end do
      end do

      !-----------------------------------------------------
      !---(10) Terminate program if joint failure occurs ---
      !-----------------------------------------------------
      if( any(failurej) ) then
         write(*,*)         " => JOINT ELEMENT SLIPPED UNDER INITIAL STRESS CONDITIONS."
         write(ounit_res,*) " => JOINT ELEMENT SLIPPED UNDER INITIAL STRESS CONDITIONS."

         write(*,*)         " => Parameters related to the joint strength have been reset. Calculation continues..."
         write(ounit_res,*) " => Parameters related to the joint strength have been reset. Calculation continues..."
         !stop
      end if

      failurej(:) = .false. ! reset
   end subroutine evalInitFail

   !> 
   ! evaluate initial robustness from initial displacement.
   subroutine evalInitRobust
      integer :: iel, ip
      real(8) :: points(ndim,param%nips_stress) ! integrating point local coordinates for solids
      real(8) :: weights(param%nips_stress)  ! weighting coefficients for solids
      real(8) :: funs(param%nods,param%nips_stress)   ! shape function ! (nods)
      real(8) :: ders(param%nods,ndim,param%nips_stress)  ! shape function derivatives wrt local coordinates (ndim,nods)
      type(propertys), pointer :: prop
      real(8) :: Dmat(nst,nst) ! D matrix (nst,nst)
      real(8) :: coords(ndim,param%nods) ! element nodal coordinates ()
      real(8) :: displs(ndim*param%nods)     ! element nodal displacements (ndim*nods)
      real(8) :: jac(ndim,ndim), invjac(ndim,ndim)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
      real(8) :: gfun(param%nods,ndim) ! shape function derivatives wrt global coordinates (ndim,nods)
      real(8) :: Bmat(nst,ndim*param%nods)   ! strain-displacement matrix (nst,ndim*nods)
      real(8) :: deps(nst), eps(nst) ! strain terms (nst)
      real(8) :: dsig(nst), sig(nst) ! stress terms

      real(8) :: principal(ndim) ! (ndim)
      real(8) :: radius, center, distance

      call samples( param%nips_stress, points, weights )
      do ip=1, param%nips_stress !nips
         call shape_funs( funs(:,ip),   points(:,ip) ) ! Nm at point
         call shape_ders( ders(:,:,ip), points(:,ip) ) ! dNm/d(xi,et,ze)
      end do

      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         call deemats( Dmat, prop%E, prop%nu ) ! D

         coords(:,:) = geom%coord(:,geom%enodes(:,iel))
         displs(:) = reshape( geom%displ(:,geom%enodes(:,iel)), shape=[ndim*param%nods] )

         do ip=1, param%nips_stress !nips
            jac(:,:) = matmul(coords,ders(:,:,ip))   !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )              ! B=L*Nm

            deps(:) = matmul(Bmat,displs)
            dsig(:) = matmul(Dmat,deps)

            eps(:) = strains_ini(:,iel) + deps(:)
            sig(:) = stresss_ini(:,iel) + dsig(:)

            energy(ip,iel) = &
               + 0.5d0*( sig(1)*eps(1) +   sig(2)*eps(2) +   sig(3)*eps(3) &
               +       2*sig(4)*eps(4) + 2*sig(5)*eps(5) + 2*sig(6)*eps(6) )

            ! Compute formation roboustness index
            call get_principal_stress( principal, sig )
            radius = 0.5d0*(principal(3)-principal(1))
            center = 0.5d0*(principal(3)+principal(1)) + presss_ini(iel) ! use initial pressure
            distance = abs(prop%mu*center-prop%So) / sqrt(prop%mu**2+1.0)
            robust(ip,iel) = distance - radius
         end do
      end do

      robust_ini(:,:) = robust(:,:)
      robust_dif(:,:) = 0.0d0
   end subroutine evalInitRobust

   !> 
   ! calculate equivalent nodal load.
   subroutine setupLoad
      call addPointLoad
      call addSurfaceLoad
      call addPoreLoad
      call addGravLoad
      call addTempLoad
      call addJointPoreLoad
      call fixZeroLoad
   end subroutine setupLoad

   !> 
   ! add point load to equivalent nodal load.
   subroutine addPointLoad
      integer :: n, i, k
      real(8) :: point_load(ndim)

      !----------------------------------
      !---(11) Point loading (File 4) ---
      !----------------------------------
      load(:,:) = 0.0d0 ! initialize load vector

      read(iunit_point_load,*)
      read(iunit_point_load,*) n

      do i=1, n
         read(iunit_point_load,*) k, point_load(:)
         load(:,k) = load(:,k) + point_load(:)
      end do
   end subroutine addPointLoad

   !> 
   ! add surface load to equivalent nodal load.
   subroutine addSurfaceLoad
      integer :: ip, iel, k, m, n
      real(8) :: pointj(2,param%nipb), weightj(param%nipb) ! integrating point local coordinates for joints
      real(8) :: funj(param%nodb,param%nipb), derj(param%nodj/2,2,param%nipb) ! shape function
      integer :: nsnode(param%nodb)
      real(8) :: coordf(ndim,param%nodb)
      real(8) :: jac(ndim,2)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
!!$      real(8) :: Lmat(ndim,ndim), gvect(ndim), surfj_laod(ndim)
      real(8) :: fmj(ndim,param%nodb)

      !----------------------------------
      !---(12) Surface loads (File 5) ---
      !----------------------------------
      call samplej( param%nipb, pointj, weightj )
      do ip=1, param%nipb
         call shape_funj( funj(:,ip),   pointj(:,ip) )
         call shape_derj( derj(:,:,ip), pointj(:,ip) )
      end do

      read(iunit_surfj_load,*)
      read(iunit_surfj_load,*) n

      if( .not. allocated( surfj_load ) ) allocate( surfj_load(ndim,param%nodb,n) )

      do iel=1, n
         read(iunit_surfj_load,*) k, m, nsnode(:)
         read(iunit_surfj_load,*) surfj_load(:,:,iel)

         coordf(:,:) = geom%coord(:,nsnode(:))

         fmj(:,:) = 0.0d0
         do ip=1, param%nipb
            jac(:,:) = matmul(coordf,derj(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)
            
            det = determinantj( jac )               ! |d(x,y,z)/d(xi,et,ze)|
!!$!!!!!
!!$            call ellmatj( Lmat, jac )
!!$            Lmat(:,1)=(/1.0,0.0,0.0/)
!!$            Lmat(:,2)=(/0.0,1.0,0.0/)
!!$            Lmat(:,3)=(/0.0,0.0,1.0/)
!!$            fmj(:,:) = fmj(:,:) &
!!$              + outer_product(matmul(Lmat,matmul(surfj_load,funj(:,ip))),funj(:,ip)) * det*weightj(ip)
!!$!!!!!
            fmj(:,:) = fmj(:,:) &
               + outer_product(matmul(surfj_load(:,:,iel),funj(:,ip)),funj(:,ip)) * det*weightj(ip)
         end do

         load(:,nsnode(:)) = load(:,nsnode(:)) + fmj(:,:)
      end do
   end subroutine addSurfaceLoad

   !> 
   ! add pore pressure to equivalent nodal load.
   subroutine addPoreLoad
      integer :: iel, j, ip
      real(8) :: points(ndim,param%nips) ! integrating point local coordinates for solids
      real(8) :: weights(param%nips)  ! weighting coefficients for solids
      real(8) :: funs(param%nods,param%nips)   ! shape function ! (nods)
      real(8) :: ders(param%nods,ndim,param%nips)  ! shape function derivatives wrt local coordinates (ndim,nods)
      real(8) :: pore_load_inc(nst)
      type(propertys), pointer :: prop
      real(8) :: coords(ndim,param%nods) ! element nodal coordinates ()
      real(8) :: jac(ndim,ndim), invjac(ndim,ndim)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
      real(8) :: gfun(param%nods,ndim) ! shape function derivatives wrt global coordinates (ndim,nods)
      real(8) :: Bmat(nst,ndim*param%nods)   ! strain-displacement matrix (nst,ndim*nods)
      real(8) :: c ! coefficient of Nur's effective stress.
      real(8) :: fm(ndim*param%nods)

      !---------------------------------------------------------------
      !---(13) Loading by pore pressure of solid elements (File 6) ---
      !---------------------------------------------------------------
      call samples(param%nips,points,weights)
      do ip=1, param%nips
         call shape_funs( funs(:,ip),   points(:,ip) ) ! Nm at point
         call shape_ders( ders(:,:,ip), points(:,ip) ) ! dNm/d(xi,et,ze)
      end do

      read(iunit_presss_load,*)
      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         c = 1.0 - prop%E/(1.0-2*prop%nu) * (1.0-2*prop%num)/prop%Em
         
         read(iunit_presss_load,*) j, presss(j) ! j==iel
         pore_load_inc(1:3) = presss(iel) - presss_ini(iel) ! j==iel
         pore_load_inc(4:6) = 0.0d0

         coords(:,:) = geom%coord(:,geom%enodes(:,iel))

         fm(:)=0.0d0
         do ip=1, param%nips
            jac(:,:) = matmul(coords,ders(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )              ! B=L*Nm

            fm(:) = fm(:) + c * matmul(transpose(Bmat),pore_load_inc) * det*weights(ip)
         end do

         load(:,geom%enodes(:,iel)) = load(:,geom%enodes(:,iel)) &
            + reshape( fm(:), shape=[3,param%nods] )
      end do
   end subroutine addPoreLoad
   !>
   ! add temperature to equivalent nodal load.
   subroutine addTempLoad
      integer :: iel, j, ip
      real(8) :: points(ndim,param%nips) ! integrating point local coordinates for solids
      real(8) :: weights(param%nips)  ! weighting coefficients for solids
      real(8) :: funs(param%nods,param%nips)   ! shape function ! (nods)
      real(8) :: ders(param%nods,ndim,param%nips)  ! shape function derivatives wrt local coordinates (ndim,nods)
      real(8) :: temp_load_inc(nst)
      real(8) :: Dmat(nst,nst) ! D matrix (nst,nst)
      type(propertys), pointer :: prop
      real(8) :: coords(ndim,param%nods) ! element nodal coordinates ()
      real(8) :: eps(nst)
      real(8) :: jac(ndim,ndim), invjac(ndim,ndim)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
      real(8) :: gfun(param%nods,ndim) ! shape function derivatives wrt global coordinates (ndim,nods)
      real(8) :: Bmat(nst,ndim*param%nods)   ! strain-displacement matrix (nst,ndim*nods)
      real(8) :: fm(ndim*param%nods)
      !real(8) :: alpha ! thermal expansion coefficient

      !---------------------------------------------------------------
      !---(**) Loading by temperature of solid elements            ---
      !---------------------------------------------------------------
      call samples(param%nips,points,weights)
      do ip=1, param%nips
         call shape_funs( funs(:,ip),   points(:,ip) ) ! Nm at point
         call shape_ders( ders(:,:,ip), points(:,ip) ) ! dNm/d(xi,et,ze)
      end do

      read(iunit_temp_load,*)
      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         call deemats(Dmat,prop%E,prop%nu) ! D
         
         read(iunit_temp_load,*) j, temp(j) ! j==iel
         temp_load_inc(1:3) = temp(iel) - temp_ini(iel) ! j==iel
         temp_load_inc(4:6) = 0.0d0
         eps=prop%alpha*temp_load_inc
         
         coords(:,:) = geom%coord(:,geom%enodes(:,iel))
         fm(:)=0.0d0
         do ip=1, param%nips
            jac(:,:) = matmul(coords,ders(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)
            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|
            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)
            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)
            call beemats( Bmat, gfun )              ! B=L*Nm
            fm(:) = fm(:) + matmul(matmul(transpose(Bmat),Dmat),eps) * det * weights(ip)
         end do

         load(:,geom%enodes(:,iel)) = load(:,geom%enodes(:,iel)) &
            + reshape( fm(:), shape=[3,param%nods] )
      end do
   end subroutine addTempLoad
   !> 
   ! add gravity effect to equivalent nodal load.
   subroutine addGravLoad
      integer :: iel, j, ip
      real(8) :: points(ndim,param%nips) ! integrating point local coordinates for solids
      real(8) :: weights(param%nips)  ! weighting coefficients for solids
      real(8) :: funs(param%nods,param%nips)   ! shape function ! (nods)
      real(8) :: ders(param%nods,ndim,param%nips)  ! shape function derivatives wrt local coordinates (ndim,nods)
      real(8) :: pore_load(nst)
      type(propertys), pointer :: prop
      real(8) :: coord(ndim,param%nods) ! element nodal coordinates ()
      real(8) :: jac(ndim,ndim), invjac(ndim,ndim)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
      real(8) :: gfun(param%nods,ndim) ! shape function derivatives wrt global coordinates (ndim,nods)
      real(8) :: Bmat(nst,ndim*param%nods)   ! strain-displacement matrix (nst,ndim*nods)
      real(8) :: c ! coefficient of Nur's effective stress.
      real(8) :: fm(ndim*param%nods),pzi

      !---------------------------------------------------------------
      !---(**) Loading by gravity change of solid elements ---
      !---------------------------------------------------------------
      call samples(param%nips,points,weights)
      do ip=1, param%nips
         call shape_funs( funs(:,ip),   points(:,ip) ) ! Nm at point
         call shape_ders( ders(:,:,ip), points(:,ip) ) ! dNm/d(xi,et,ze)
      end do

      read(iunit_den_dif,*)
      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         
         read(iunit_den_dif,*) j, density_dif(j) ! j==iel
         
         if(density_dif(j)==0.d0)cycle
         
         coord(:,:) = geom%coord(:,geom%enodes(:,iel))
         do j=1,param%nods
         pzi=0.d0
         do ip=1, param%nips
            jac(:,:) = matmul(coord,ders(:,:,ip))   !  d(x,y,z)/d(xi,et,ze)
            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|
            pzi = pzi - param%gravity * density_dif(iel) * det*weights(ip) * funs(j,ip)
         end do
         load(3,geom%enodes(j,iel)) = load(3,geom%enodes(j,iel)) + pzi
         enddo
         
      end do
   end subroutine addGravLoad

   !> 
   ! add joint pore pressure to equivalent nodal load.
   subroutine addJointPoreLoad
      integer :: iel, ip, j
      real(8) :: pointj(2,param%nipb), weightj(param%nipb) ! integrating point local coordinates for joints
      real(8) :: funj(param%nodb,param%nipb), derj(param%nodj/2,2,param%nipb) ! shape function
      real(8) :: pore_load_inc(ndim)
      real(8) :: coordj(ndim,param%nodb)
      real(8) :: jac(ndim,2)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
      real(8) :: Lmat(ndim,ndim) ! normal vectors on the surface.
      real(8) :: Bmatj(ndim,ndim*param%nodj)
      real(8) :: fmj(ndim*param%nodj)

      if( param%nelej==0 ) return

      !---------------------------------------------------------------
      !---(XX) Loading by pore pressure of joint elements (File 7) ---
      !---------------------------------------------------------------
      call samplej(param%nipb,pointj,weightj) ! Obtain sampling points and weights for surface integration
      do ip=1, param%nipb
         call shape_funj( funj(:,ip),   pointj(:,ip) )
         call shape_derj( derj(:,:,ip), pointj(:,ip) )
      end do

      read(iunit_pressj_load,*)
      do iel=1, param%nelej
         read(iunit_pressj_load,*) j, pressj(j) ! j==iel
         pore_load_inc = [ 0.0d0, 0.0d0, pressj(iel) - pressj_ini(iel) ] ! j==iel

         coordj(:,:) = geom%coord(:,geom%enodej(1:param%nodj/2,iel))

         fmj(:) = 0.0d0
         do ip=1, param%nipb
            jac(:,:) = matmul(coordj,derj(:,:,ip))   !  d(x,y,z)/d(xi,et,ze)
            
            det = determinantj( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            call ellmatj( Lmat, jac ) ! normal vectors on the surface.

            call beematj( Bmatj, Lmat, funj(:,ip) )              ! B=L*Nm

            fmj(:) = fmj(:) + matmul(transpose(Bmatj),pore_load_inc) * det*weightj(ip)  
         end do

         load(:,geom%enodej(:,iel)) = load(:,geom%enodej(:,iel)) &
            + reshape( fmj(:), shape=[3,param%nodj] )
      end do
   end subroutine addJointPoreLoad

   !> 
   ! set load at fixed nodes to zero.
   subroutine fixZeroLoad
      integer :: n, k
      !--- modify force vector for fixed coords.
      do n=1, param%nnode
         do k=1, 3
            if( geom%nodfree(k,n) == 0 ) then ! fixed coord
               load(k,n) = 0.0d0
            end if
         end do
      end do
   end subroutine fixZeroLoad

   !> 
   ! make stiffness matrix.
   subroutine setupStiff
      integer :: ip, iel, ei, ej, i, j, l
      real(8) :: points(ndim,param%nips), weights(param%nips)
      real(8) :: funs(param%nods,param%nips), ders(param%nods,ndim,param%nips)
      type(propertys), pointer :: prop
      real(8) :: Dmat(nst,nst)
      real(8) :: Bmat(nst,ndim*param%nods)
      real(8) :: coords(ndim,param%nods)
      real(8) :: jacs(ndim,ndim), invjac(ndim,ndim), det, gfun(param%nods,ndim)
      real(8), allocatable :: km(:,:)

      real(8) :: pointj(2,param%nipb), weightj(param%nipb)
      real(8) :: funj(param%nodj/2,param%nipb), derj(param%nodj/2,2,param%nipb)
      real(8) :: Dmatj(ndim,ndim)
      real(8) :: Bmatj(ndim,ndim*param%nodj)
      real(8) :: coordj(ndim,param%nodj/2)
      real(8) :: jacj(ndim,2), Lmat(ndim,ndim) ! normal vectors on the surface.
      real(8), allocatable :: kmj(:,:)

      allocate( km(ndim*param%nods,ndim*param%nods), kmj(ndim*param%nodj,ndim*param%nodj) )
      
      !
      !---(XX) Compute element stiffness matrix for solid elements ---
      !
      stiffmat(:) = 0.0d0

      call samples( param%nips, points, weights )
      do ip=1, param%nips
         call shape_funs( funs(:,ip),   points(:,ip) ) ! Nm at point.
         call shape_ders( ders(:,:,ip), points(:,ip) ) ! dNm/d(xi,et,ze)
      end do

      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
      
         call deemats( Dmat, prop%E, prop%nu ) ! D

         coords(:,:) = geom%coord(:,geom%enodes(:,iel))

         km(:,:) = 0.0d0
         do ip=1, param%nips
            jacs(:,:) = matmul(coords,ders(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jacs )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jacs )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac)  ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )               ! B=L*Nm

            km(:,:) = km(:,:) + matmul(matmul(transpose(Bmat),Dmat),Bmat) * det*weights(ip)
         end do

         do ei=1, ndim*param%nods
            i=geom%esteers(ei,iel)
            do ej=1, ndim*param%nods
               j=geom%esteers(ej,iel)

               if( upper_only .and. i>j ) cycle ! skip lower triangle part of th matrix.

               l = geom%eindexs(ei,ej,iel)
               stiffmat(l) = stiffmat(l) + km(ei,ej)
            end do ! ej
         end do ! ei
      end do ! iel

      if( param%nelej>0 ) then
         !
         !---(XX) Compute element stiffness matrix for joint elements ---
         !
         call samplej( param%nipb, pointj, weightj ) ! Obtain sampling points and weights for surface integration
         do ip=1, param%nipb
            call shape_funj( funj(:,ip),   pointj(:,ip) )
            call shape_derj( derj(:,:,ip), pointj(:,ip) )
         end do

         do iel=1, param%nelej
            coordj(:,:) = geom%coord(:,geom%enodej(1:param%nodj/2,iel))

            kmj=0.0d0
            do ip=1, param%nipb
               jacj(:,:) = matmul(coordj,derj(:,:,ip))
               det = determinantj( jacj )               ! |d(x,y,z)/d(xi,et,ze)|

               call ellmatj( Lmat, jacj ) ! normal vectors on the surface.

               call deematj( Dmatj, stiffj(:,ip,iel) ) ! D
               call beematj( Bmatj, Lmat, funj(:,ip) )              ! B=L*Nm

               kmj(:,:) = kmj(:,:) + matmul(matmul(transpose(Bmatj),Dmatj),Bmatj) * det*weightj(ip)
            end do

            do ei=1, ndim*param%nodj
               i=geom%esteerj(ei,iel)
               do ej=1, ndim*param%nodj
                  j=geom%esteerj(ej,iel)

                  if( upper_only .and. i>j ) cycle ! skip lower triangle part of the matrix.

                  l = geom%eindexj(ei,ej,iel)
                  stiffmat(l) = stiffmat(l) + kmj(ei,ej)
               end do ! ej
            end do ! ei
         end do ! iel
      end if
      
      !--- modify stiffness matrix for fixed coords.
      do i=1, ndim*param%nnode
         if( geom%codfree(i) == 0 ) then ! fixed coord
            do l=geom%ia(i), geom%ia(i+1)-1
               j=geom%ja(l)
               if( i==j ) then
                  stiffmat(l) = 1.0
               else
                  stiffmat(l) = 0.0d0
               end if
            end do
         else
            do l=geom%ia(i), geom%ia(i+1)-1
               j=geom%ja(l)
               if( geom%codfree(j) == 0 ) then ! fixed coord
                  stiffmat(l)=0.0d0
               end if
            end do
         end if
      end do

      deallocate( km, kmj )
   end subroutine setupStiff

   !> 
   ! evaluate joint failure from traction.
   subroutine evalFail( residual )
      real(8), intent(out) :: residual

      integer :: iel, ip
      real(8) :: pointj(2,param%nipb), weightj(param%nipb) ! integrating point local coordinates for joints
      real(8) :: funj(param%nodb,param%nipb), derj(param%nodj/2,2,param%nipb) ! shape function
      type(propertyj), pointer :: prop
      real(8) :: coordj(ndim,param%nodb)
      real(8) :: displj_inc(ndim*param%nodj)     ! element nodal displacements (ndim*nodj)
      real(8) :: jac(ndim,2)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
      real(8) :: Lmat(ndim,ndim) ! normal vectors on the surface.
      real(8) :: Dmatj(ndim,ndim)
      real(8) :: Bmatj(ndim,ndim*param%nodj)
      real(8) :: fmj(ndim*param%nodj)
      real(8) :: tractionj_inc(ndim)
      real(8) :: Fn, Fs, ts, tn
      real(8) :: dvp(ndim)
      
      residual = 0.0d0
      if( param%nelej==0 ) return
      
      load(:,:) = 0.0d0

      !-------------------------------------------------
      !---(XX) Check the stability of joint elements ---
      !-------------------------------------------------
      call samplej( param%nipb, pointj, weightj )

      ! Compute shape function matrix
      do ip=1, param%nipb
         call shape_funj( funj(:,ip),   pointj(:,ip) )
         call shape_derj( derj(:,:,ip), pointj(:,ip) )
      end do

      do iel=1, param%nelej
         prop => param%propj(geom%etypej(iel))

         coordj(:,:) = geom%coord(:,geom%enodej(1:param%nodj/2,iel))
         displj_inc(:) = reshape( geom%displ_inc(:,geom%enodej(:,iel)), shape=[3*param%nodj] )

         fmj(:)=0.0d0
         do ip=1, param%nipb
            jac(:,:) = matmul(coordj,derj(:,:,ip))   !  d(x,y,z)/d(xi,et,ze)
            
            det = determinantj( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            call ellmatj( Lmat, jac ) ! normal vectors on the surface.

            ! Compute D matrix for joint element
            call deematj( Dmatj, stiffj(:,ip,iel) ) ! D
            call beematj( Bmatj, Lmat, funj(:,ip) )              ! B=L*Nm

            tractionj_inc(:) = matmul(matmul(Dmatj,Bmatj),displj_inc)

            tractionj(:,ip,iel) = tractionj(:,ip,iel) + tractionj_inc(:)

            if( failurej(iel) ) then
               Fn = tractionj(3,ip,iel)
            else
               Fn = tractionj(3,ip,iel) - tensilej(ip,iel)
            end if

            ts = parallel(      tractionj(:,ip,iel) )
            tn = perpendicular( tractionj(:,ip,iel) )
            if( Fn>=0.0d0 ) then
               Fs = ts
            else
               Fs = ts + prop%mu*tn - cohesionj(ip,iel) 
            end if

            if( Fs>0.0d0 .and. Fn<0.0d0 ) then
               fs = ts + prop%mu*tn
            end if

            if( ts/=0.0d0 ) then
               dvp(1) = phi(Fs)*tractionj(1,ip,iel)/ts/stiffj(1,ip,iel)
               dvp(2) = phi(Fs)*tractionj(2,ip,iel)/ts/stiffj(2,ip,iel)
            else
               dvp(1) = 0.0d0
               dvp(2) = 0.0d0
            end if
            dvp(3) = phi(fn)/stiffj(3,ip,iel)

            tractionj(:,ip,iel) = tractionj(:,ip,iel) - matmul(Dmatj,dvp)

            ! Recalculate Dmatrix with updated ks1,ks2,kn values
            call deematj( Dmatj, stiffj(:,ip,iel) ) ! D

            ! Compute contribution to load vector
            fmj(:) = fmj(:) + matmul(matmul(transpose(Bmatj),Dmatj),dvp) * det*weightj(ip)

            ! Update parameters in stiffness matrix for joint elements
            if( Fs>1.0d-6 ) then
               stiffj(1:2,ip,iel) = prop%ksn(1:2) * prop%multksn 
               failurej(iel) = .true.
            end if

            if( Fn>0.0d0 ) then
               stiffj(1:3,ip,iel) = prop%ksn(1:3) * prop%multksn 
               failurej(iel) = .true.
            else
               stiffj(3,ip,iel) = prop%ksn(3)
            end if
         end do ! ip

         load(:,geom%enodej(:,iel)) = load(:,geom%enodej(:,iel)) &
            + reshape( fmj(:), shape=[3,param%nodj] )
      end do ! iel

      call fixZeroLoad

      residual = sqrt( sum( load(:,:)**2 ) )

   contains
      function phi(F)
         real(8), intent(in) :: F
         real(8) :: phi
         real(8), parameter :: small=1.d-3

         phi = MERGE( 0.0d0, F, F<=small )
      end function phi
   end subroutine evalFail

#ifdef ENABLE_PARDISO
   !>
   ! solve the equilibrium equation of the stiffness matrix by PARDISO.
   subroutine solveDispl_PARDISO
      use mkl_pardiso
      
      !--- PARDISO control parameter
      type(MKL_PARDISO_HANDLE) :: pt(64)
      integer :: maxfct, mnum, mtype, phase, nrhs, iparm(64), msglvl, error, idum(1)
      real(8) :: ddum(1)

!!$      integer :: i, l, j ! DEBUG
!!$      real(16), allocatable :: work(:) ! DEBUG

      !--- Initiliaze internal solver memory pointer.
      pt(:)%DUMMY = 0 
      maxfct    =  1
      mnum      =  1
      mtype     = -2 ! symmetric, indefinite
      nrhs      =  1

      !--- Set up PARDISO control parameter.
      iparm( :) =  0
      iparm( 1) =  1 ! no solver default
      iparm( 2) =  2 ! fill-in reordering from METIS
      iparm( 4) =  0 ! 62 
      iparm( 5) =  0 ! no user fill-in reducing permutation
      iparm( 6) =  0 ! =0 solution on the first n compoments of x
      iparm( 8) =  2 ! numbers of iterative refinement steps
      iparm(10) = 13 ! perturbe the pivot elements with 1E-8
      iparm(11) =  1 ! use nonsymmetric permutation and scaling MPS
      iparm(13) =  0 ! maximum weighted matching algorithm is switched-off
      !              !  (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
      iparm(14) =  0 ! Output: number of perturbed pivots
      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      iparm(19) = -1 ! Output: Mflops for LU factorization
      iparm(20) =  0 ! Output: Numbers of CG Iterations
      msglvl    =  1 ! print statistical information
      error     =  0 ! initialize error flag

      phase = 11 ! only reordering and symbolic factorization
      call pardiso( pt, maxfct, mnum, mtype, phase, ndim*param%nnode, &
         stiffmat, geom%ia, geom%ja, idum, nrhs, iparm, msglvl, ddum, ddum, error )
      if( error/=0 ) then
         write(*,*) 'The following ERROR was detected: ', error
         call exit(1)
      end if
      write(*,*) 'Reordering completed ... '
      write(*,*) 'Number of nonzeros in factors = ',iparm(18)
      write(*,*) 'Number of factorization MFLOPS = ',iparm(19)

      !.. Factorization.
      phase = 22 ! only factorization
      call pardiso( pt, maxfct, mnum, mtype, phase, ndim*param%nnode, &
         stiffmat, geom%ia, geom%ja, idum, nrhs, iparm, msglvl, ddum, ddum, error )
      if( error/=0 ) then
         write(*,*) 'The following ERROR was detected: ', error
         call exit(1)
      end if
      write(*,*) 'Factorization completed ... '

      !.. Back substitution and iterative refinement
      iparm(8) = 2 ! max numbers of iterative refinement steps
      phase = 33 ! only solving
      call pardiso( pt, maxfct, mnum, mtype, phase, ndim*param%nnode, &
         stiffmat, geom%ia, geom%ja, idum, nrhs, iparm, msglvl, load, geom%displ_inc, error )
      if( error/=0 ) then
         write(*,*) 'The following ERROR was detected: ', error
         call exit(1)
      end if
      write(*,*) 'Solve completed ... '

      !.. Termination and release of memory
      phase = -1 ! release internal memory
      call pardiso( pt, maxfct, mnum, mtype, phase, ndim*param%nnode, ddum, &
         idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, error )
      if( error/=0 ) then
         write(*,*) 'The following ERROR on release stage was detected: ',error
         call exit(1)
      end if

      !-----------------------------------------
      !---(XX) Add incremental displacements ---
      !-----------------------------------------
      geom%displ(:,:) = geom%displ(:,:) + geom%displ_inc(:,:)


!!$      !! DEBUG
!!$      write(100,*) "nnode: ", param%nnode
!!$      write(100,*) "neq:", ndim*param%nnode
!!$      write(100,*) "nnz:", geom%nnz
!!$
!!$      do i=1, 3*param%nnode
!!$         write(100,*)
!!$         do l=geom%ia(i), geom%ia(i+1)-1
!!$            write(100,"(2i5,e15.5)") (i-1)/3+1, (geom%ja(l)-1)/3+1, stiffmat(l)
!!$         end do
!!$      end do
!!$
!!$      do i=1, param%nnode
!!$         write(200,"(1i5,3e15.5)") i, load(:,i)
!!$      end do
!!$
!!$      do i=1, param%nnode
!!$         write(300,"(1i5,3e15.5)") i,  geom%displ_inc(:,i)
!!$      end do
!!$
!!$      allocate( work(3*param%nnode) )
!!$      work(:) = 0.0d0
!!$      do i=1, 3*param%nnode
!!$         do l=geom%ia(i), geom%ia(i+1)-1
!!$            j=geom%ja(l)
!!$            ! upper
!!$            work(i) = work(i) + stiffmat(l)*geom%displ_inc(mod(j-1,3)+1,(j-1)/3+1)
!!$            if( i/=j ) then
!!$               ! lower
!!$               work(j) = work(j) + stiffmat(l)*geom%displ_inc(mod(i-1,3)+1,(i-1)/3+1)
!!$            end if
!!$         end do
!!$      end do
!!$      do i=1, 3*param%nnode
!!$         write(400,"(1i5,3e15.5)") i, work(i), load(mod(i-1,3)+1,(i-1)/3+1)
!!$      end do
!!$      deallocate( work )
!!$
!!$      stop "DEBUG"
!!$      !! DEBUG       
   end subroutine solveDispl_PARDISO
#endif
#ifdef ENABLE_PETSC
   !>
   ! solve the equilibrium equation of the stiffness matrix by PETSc.ARDISO.
   subroutine solveDispl_PETSc
      use petscsnes
      use petscksp

      !--- PETSc parameters.
      PetscInt,  parameter :: max_iterations = 1000
      PetscReal, parameter :: tolerance = 1.0d-20
      PetscInt  :: iterations
      PetscReal :: residual
      KSP :: ksp  ! linear solver context
      PC  :: pc   ! preconditioner context
      Mat :: Apetsc
      Vec :: bpetsc, xpetsc
      PetscErrorCode :: error
      integer :: n
      PetscInt :: row(1), col(1)
      PetscScalar :: vala(1)
      PetscScalar, allocatable :: valb(:), valx(:)

!!$      integer :: i, l, j ! DEBUG
!!$      real(16), allocatable :: work(:) ! DEBUG

      allocate( valb(ndim*param%nnode), valx(ndim*param%nnode) )
      
      !   petsc array for matrix data.
      call MatCreate( PETSC_COMM_WORLD, Apetsc, error )
      call MatSetSizes( Apetsc, PETSC_DECIDE, PETSC_DECIDE, ndim*param%nnode, ndim*param%nnode, error )
      call MatSetFromOptions( Apetsc, error )
      call MatSetUP( Apetsc, error )
      call MatSetOption( Apetsc, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, error )
      call MatSetOption( Apetsc, MAT_SYMMETRIC, PETSC_TRUE, error )

      !   petsc array for RHS data.
      call VecCreate( PETSC_COMM_WORLD, bpetsc, error )
      call VecSetSizes( bpetsc, PETSC_DECIDE, ndim*param%nnode, error )
      call VecSetFromOptions( bpetsc, error )
      call VecSetUp( bpetsc, error )

      !   petsc array for solution data.
      call VecCreate( PETSC_COMM_WORLD, xpetsc, error )
      call VecSetSizes( xpetsc, PETSC_DECIDE, ndim*param%nnode, error )
      call VecSetFromOptions( xpetsc, error )
      call VecSetUp( xpetsc, error )

      !   create KSP.
      call KSPCreate( PETSC_COMM_WORLD, ksp, error )
      call KSPSetOperators( ksp, Apetsc, Apetsc, error )

      !   set precondition.
      call KSPGetPC( ksp, pc, error )

      !   SOR.
      call PCSetType( pc, PCSOR, error )
      call PCSORSetIterations( pc, 10, 1, error )
      call PCSORSetOmega( pc, 1.7d0, error )

      !   set KSP and solve.
      call KSPSetType( ksp, KSPBCGS, error )
      call KSPSetTolerances( ksp, tolerance, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, max_iterations, error )
      call KSPSetInitialGuessNonzero( ksp, PETSC_TRUE, error )
      call KSPSetFromOptions( ksp, error )

      !   set matrix data.
      do n=1, geom%nnz
         row(1) = geom%irn(n)
         col(1) = geom%jcn(n)
         vala(1) = stiffmat(n) ! to PetscScalar from real(8).
         call MatSetValues( Apetsc, 1, row, 1, col, vala, INSERT_VALUES, error )
      end do
      call MatAssemblyBegin( Apetsc, MAT_FINAL_ASSEMBLY, error )
      call MatAssemblyEnd( Apetsc, MAT_FINAL_ASSEMBLY, error )

      do n=1, param%nnode
         valb(3*(n-1)+1) = load(1,n) ! to PetscScalar 1D array from real(8) 2D array.
         valb(3*(n-1)+2) = load(2,n) ! to PetscScalar 1D array from real(8) 2D array.
         valb(3*(n-1)+3) = load(3,n) ! to PetscScalar 1D array from real(8) 2D array.
      end do
      call VecSetValues( bpetsc, ndim*param%nnode, geom%idx, valb, INSERT_VALUES, error )
      call VecAssemblyBegin( bpetsc, error )
      call VecAssemblyEnd( bpetsc, error )

      !--- Solve the equation.
      write(*,*) "# petsc start"
      call KSPSolve( ksp, bpetsc, xpetsc, error )
      if( error /= 0 ) write(*,*) "Errror: PETSC"
      write(*,*) "# petsc finished"
      if( error/=0 ) stop "ERROR KSPSolve"

      !   output iterations and residual.
      call KSPGetIterationNumber( ksp, iterations, error )
      call KSPGetResidualNorm( ksp, residual, error )

      write(*,'(a,i5,a,1p,e12.5)') 'BiCGStab  Iterations : ', iterations, ', Residual Norm : ', residual

      call VecGetValues( xpetsc, ndim*param%nnode, geom%idx, valx, error )
      do n=1, param%nnode
         geom%displ_inc(1,n) = valx(3*(n-1)+1) ! to real(8) 2D array from PetscScalar 1D array.
         geom%displ_inc(2,n) = valx(3*(n-1)+2) ! to real(8) 2D array from PetscScalar 1D array.
         geom%displ_inc(3,n) = valx(3*(n-1)+3) ! to real(8) 2D array from PetscScalar 1D array.
      end do

      call MatDestroy( Apetsc, error )
      call VecDestroy( bpetsc, error )
      call VecDestroy( xpetsc, error )
      call KSPDestroy( ksp, error )

      deallocate( valb, valx )

      
      !-----------------------------------------
      !---(XX) Add incremental displacements ---
      !-----------------------------------------
      geom%displ(:,:) = geom%displ(:,:) + geom%displ_inc(:,:)

      
!!$      !! DEBUG
!!$      write(100,*) "nnode: ", param%nnode
!!$      write(100,*) "neq:", ndim*param%nnode
!!$      write(100,*) "nnz:", geom%nnz
!!$
!!$      do i=1, 3*param%nnode
!!$         write(100,*)
!!$         do l=geom%ia(i), geom%ia(i+1)-1
!!$            write(100,"(2i5,e15.5)") (i-1)/3+1, (geom%ja(l)-1)/3+1, stiffmat(l)
!!$         end do
!!$      end do
!!$
!!$      do i=1, param%nnode
!!$         write(200,"(1i5,3e15.5)") i, load(:,i)
!!$      end do
!!$
!!$      do i=1, param%nnode
!!$         write(300,"(1i5,3e15.5)") i,  geom%displ_inc(:,i)
!!$      end do
!!$
!!$      allocate(work(3*param%nnode))
!!$      work(:) = 0.0d0
!!$      do n=1, geom%nnz
!!$         i=geom%irn(n)+1
!!$         j=geom%jcn(n)+1
!!$         work(i) = work(i) + stiffmat(n) * geom%displ_inc(mod(j-1,3)+1,(j-1)/3+1)
!!$      end do
!!$      do i=1, 3*param%nnode
!!$         write(400,"(1i5,3e15.5)") i, work(i), load(mod(i-1,3)+1,(i-1)/3+1)
!!$      end do
!!$      deallocate(work)
!!$
!!$      stop "DEBUG"
!!$      !! DEBUG       
   end subroutine solveDispl_PETSc
#endif
   
   !>
   ! evaluate seismic momemnt from displacement.
   subroutine evalSeismic
      integer :: iel, ip, j
      real(8) :: pointj(2,param%nipb), weightj(param%nipb) ! integrating point local coordinates for joints
      real(8) :: funj(param%nodb,param%nipb), derj(param%nodj/2,2,param%nipb) ! shape function
      type(propertyj), pointer :: prop
      real(8) :: ts, ts_ini
      real(8) :: coordj(ndim,param%nodj/2) ! element nodal coordinates (ndim,param%nodj/2)
      real(8) :: displj(ndim,param%nodj), gap(ndim,param%nodj/2), slip(param%nodj/2)
      real(8) :: jac(ndim,2)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
      real(8) :: Lmat(ndim,ndim) ! normal vectors on the surface.
      real(8) :: eff_shear_modulus, moment, area, mo,magnitude
      
      if( param%nelej==0 ) return
      ! Compute Coulomb stress, Coulomb stress change and seismic moment

      call samplej(param%nipb,pointj,weightj) ! Obtain sampling points and weights for surface integration
      do ip=1, param%nipb
         call shape_funj( funj(:,ip),   pointj(:,ip) )
         call shape_derj( derj(:,:,ip), pointj(:,ip) )
      end do

      do iel=1, param%nelej
         prop => param%propj(geom%etypej(iel))
         eff_shear_modulus = prop%shear_mod
         coordj(:,:) = geom%coord(:,geom%enodej(1:param%nodj/2,iel))
         displj(:,:) = geom%displ(:,geom%enodej(:,iel))
         gap(:,:)    = displj(:,1:param%nodj/2) - displj(:,param%nodj/2+1:param%nodj)
         
         moment = 0.0d0
         area = 0.0d0
         
         do ip=1, param%nipb
            !ts     = parallel( tractionj    (:,ip,iel) )
            !ts_ini = parallel( tractionj_ini(:,ip,iel) )
            !drop   = ts_ini - ts

            jac(:,:) = matmul(coordj,derj(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)

            det = determinantj( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            call ellmatj( Lmat, jac ) ! normal vectors on the surface.

            do j=1, param%nodj/2
               slip(j) = parallel( gap(:,j), Lmat )
            end do
            moment = moment + eff_shear_modulus * sum(slip(:)*funj(:,ip)) * det*weightj(ip)
            !if( drop>0.0d0 ) then
            !   moment = moment + drop * sum(slip(:)*funj(:,ip)) * det*weightj(ip)
            !end if
            area = area + det*weightj(ip)
         end do ! ip

         seismic_mom(iel) = moment
         seismic_int(iel) = moment/area
      end do ! iel

      mo = sum(seismic_mom)
      magnitude = -9999.9999999999d0
      if( mo/=0.0d0 ) magnitude = (log10(mo)-9.05d0)/1.5d0 ! Kanamori model

      write(ounit_log,*)"* Induced Seismicity"
      write(ounit_log,'(1x,"Seismic Moment = ",es20.6," [N-m]")') mo
      write(ounit_log,'(1x,"Magnitude      = ",es20.6)') magnitude
   end subroutine evalSeismic

   !> 
   ! calculate coulomb stress from traction.
   subroutine evalCoulomb
      integer :: iel, ip
      type(propertyj), pointer :: prop
      real(8) :: ts, tn
      
      do iel=1, param%nelej
         prop => param%propj(geom%etypej(iel))

         do ip=1, param%nipb
            ts = parallel(      tractionj(:,ip,iel) )
            tn = perpendicular( tractionj(:,ip,iel) )

            coulombj(ip,iel) = ts + prop%mu*tn
         end do
      end do
            
      coulombj_dif(:,:) = coulombj(:,:) - coulombj_ini(:,:)
   end subroutine evalCoulomb
      
   !> 
   ! evaluate robustness from displacement.
   subroutine evalRobust
      integer :: iel, ip
      real(8) :: points(ndim,param%nips_stress) ! integrating point local coordinates for solids
      real(8) :: weights(param%nips_stress)  ! weighting coefficients for solids
      real(8) :: funs(param%nods,param%nips_stress)   ! shape function ! (nods)
      real(8) :: ders(param%nods,ndim,param%nips_stress)  ! shape function derivatives wrt local coordinates (ndim,nods)
      type(propertys), pointer :: prop
      real(8) :: Dmat(nst,nst) ! D matrix (nst,nst)
      real(8) :: coords(ndim,param%nods) ! element nodal coordinates ()
      real(8) :: displs(ndim*param%nods)     ! element nodal displacements (ndim*nods)
      real(8) :: jac(ndim,ndim), invjac(ndim,ndim)   ! Jacobian matrix (ndim,ndim)
      real(8) :: det      ! determinant of the Jacobian matrix
      real(8) :: gfun(param%nods,ndim) ! shape function derivatives wrt global coordinates (ndim,nods)
      real(8) :: Bmat(nst,ndim*param%nods)   ! strain-displacement matrix (nst,ndim*nods)
      real(8) :: deps(nst), eps(nst) ! strain terms (nst)
      real(8) :: dsig(nst), sig(nst) ! stress terms

      real(8) :: principal(ndim) ! (ndim)
      real(8) :: radius, center, distance

      call samples( param%nips_stress, points, weights )
      do ip=1, param%nips_stress !nips
         call shape_funs( funs(:,ip),   points(:,ip) ) ! Nm at point
         call shape_ders( ders(:,:,ip), points(:,ip) ) ! dNm/d(xi,et,ze)
      end do

      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         call deemats( Dmat, prop%E, prop%nu ) ! D

         coords(:,:) = geom%coord(:,geom%enodes(:,iel))
         displs(:) = reshape( geom%displ(:,geom%enodes(:,iel)), shape=[ndim*param%nods] )

         do ip=1, param%nips_stress !nips
            jac(:,:) = matmul(coords,ders(:,:,ip))   !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )              ! B=L*Nm

            deps(:) = matmul(Bmat,displs)
            dsig(:) = matmul(Dmat,deps)

            eps(:) = strains_ini(:,iel) + deps(:)
            sig(:) = stresss_ini(:,iel) + dsig(:)

            energy(ip,iel) = &
               + 0.5d0*( sig(1)*eps(1) +   sig(2)*eps(2) +   sig(3)*eps(3) &
               +       2*sig(4)*eps(4) + 2*sig(5)*eps(5) + 2*sig(6)*eps(6) )

            ! Compute formation roboustness index
            call get_principal_stress( principal, sig )
            radius = 0.5d0*(principal(3)-principal(1))
            center = 0.5d0*(principal(3)+principal(1)) + presss(iel) ! use current pressure
            distance = abs(prop%mu*center-prop%So) / sqrt(prop%mu**2+1.0)
            robust(ip,iel) = distance - radius
         end do
      end do

      robust_dif(:,:) = robust(:,:) - robust_ini(:,:)
   end subroutine evalRobust

   subroutine initial_report
      integer :: iel

      write(*,'(1x,"# of max time steps, ntime =",i6)') param%ntime

      write(ounit_res,'("wFAST SIMULATION RESULTS (beta 1.0)")')
      write(ounit_res,*)
      write(ounit_res,*)
      write(ounit_res,'("***************************")')
      write(ounit_res,'("*** (1) PRE-CALCULATION ***")')
      write(ounit_res,'("***************************")')
      write(ounit_res,'(1x,"# of MAX ITERATION STEPS, NTIME =",i6)') param%ntime
      write(ounit_res,*) " freedom number for solid elements, g_g"
      do iel=1, param%neles
         write(ounit_res,"(1x,i6,24i8)") iel, geom%esteers(:,iel)
      end do
      write(ounit_res,*) "freedom number for joint elements, g_gj"
      do iel=1, param%nelej
         write(ounit_res,"(1x,i6,24i8)") iel, geom%esteerj(:,iel)
      end do
      write(ounit_res,'((A,I8,A))') " There are ", geom%nnz, &
         " non-zero elements of the coefficient matrix A."

      write(ounit_res,*)
      write(ounit_res,*)
      write(ounit_res,'("******************************")')
      write(ounit_res,'("*** (2) INITIAL CONDITIONS ***")')
      write(ounit_res,'("******************************")')
   end subroutine initial_report

   subroutine progress_report( itime, residual )
      integer, intent(in) :: itime
      real(8), intent(in) :: residual

      ! On screen
      write(*,'(1x,"NONLINEAR CALCULATION: ",i3,"-TH ITERATION")') itime
      write(*,'(1x,"RESIDUAL =",es16.6)') residual

      ! In result file
      write(ounit_res,'(1x,"NONLINEAR CALCULATION: ",i3,"-TH ITERATION")') itime
      write(ounit_log,'(1x,"NONLINEAR CALCULATION: ",i3,"-TH ITERATION")') itime
      write(ounit_res,'(1x,"RESIDUAL =",es16.6)') residual
      write(ounit_log,'(1x,"RESIDUAL =",es16.6)') residual
   end subroutine progress_report

   subroutine post
      implicit none

      integer :: n, iel, i, ip
      real(8) :: displj(ndim,param%nodj), gap(ndim), gapmax(ndim)
      real(8) :: points(ndim,param%nips), weights(param%nips)
      real(8) :: funs(param%nods,param%nips), ders(param%nods,ndim,param%nips)
      type(propertys), pointer :: prop
      real(8) :: Dmat(nst,nst)
      real(8) :: Bmat(nst,ndim*param%nods)
      real(8) :: coords(ndim,param%nods), displs(ndim*param%nods)
      real(8) :: pos(ndim), jac(ndim,ndim), invjac(ndim,ndim), det, gfun(param%nods,ndim)
      real(8) :: deps(nst), dsig(nst), eps(nst), sig(nst)      

      character(len=12), parameter :: line = "------------"
      character(len= 3), parameter :: length = "[L]", dl = "[-]"
      character(len= 8), parameter :: psi = "[M/L/T2]"

      !--- Body of post
      write(ounit_res,'("*** NODAL DISPLACEMENTS ***")')
      write(ounit_res,'(a)') "      NODE        X-DISP        Y-DISP       Z-DISP"
      write(ounit_res,'(6x,"[-]",4x,3(6x,a3,5x))') (length,i=1,3)
      write(ounit_res,'(4(1x,a12,1x))') (line,i=1,4)
      do n=1, param%nnode
         write(ounit_res,'(1x,i8,4x,3es14.6)') n, geom%displ(1:3,n)
      end do

      gapmax(:) = 0.0d0
      do iel=1, param%nelej
         displj(:,:) = geom%displ(:,geom%enodej(:,iel))
         do i=1, param%nodj/2
            gap(:) = displj(:,i) - displj(:,i+param%nodj/2)
            gapmax(:) = max( gapmax(:), abs(gap(:)) )
         end do
      end do
      
      write(*,'(1x,"umax, vmax, wmax =",3es16.6," UNITS IN [L]")') gapmax(:)


      ! recover strains and stresses at nips integrating points
      write(ounit_res,'(/1x,a,i2,a)') "STRESSES EVALUATED AT INTEGRATION POINTS (nip=",param%nips,") ARE:"
      write(ounit_res,'("   Element  GP    X-COORD       Y-COORD       Z-COORD    ' &
         // '    SIG_X         SIG_Y         SIG_Z        SIG_XY        SIG_YZ        SIG_ZX    ' &
         // '    EPSI_X        EPSI_Y        EPSI_Z       EPSI_XY       EPSI_YZ       EPSI_ZX")') 
      write(ounit_res,'(5x,"[-]   [-]",3(6x,a3,5x),6(3x,a8,3x),6(6x,a3,5x))') &
         (length,i=1,3), (psi,i=1,6), (dl,i=1,6)
      write(ounit_res,'(2x,"-------- ----",15(1x,a12,1x))') (line,i=1,15)

      call samples( param%nips, points, weights )
      do ip=1, param%nips
         call shape_funs( funs(:,ip),   points(:,ip) ) ! Nm at point
         call shape_ders( ders(:,:,ip), points(:,ip) ) ! dNm/d(xi,et,ze)
      end do
      
      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         call deemats( Dmat, prop%E, prop%nu ) ! D

         coords(:,:) = geom%coord(:,geom%enodes(:,iel))
         displs(:) = reshape( geom%displ(:,geom%enodes(:,iel)), shape=[ndim*param%nods] )

         do ip=1, param%nips
            pos(:) = matmul(coords,funs(:,ip)  )  ! coord at point
            
            jac(:,:) = matmul(coords,ders(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )              ! B=L*Nm

            deps(:) = matmul(Bmat,displs)
            dsig(:) = matmul(Dmat,deps)

            eps(:) = deps(:) + strains_ini(:,iel)
            sig(:) = dsig(:) + stresss_ini(:,iel)

            write(ounit_res,'(i8,i6,15es14.6)') iel, ip, pos, sig, eps
         end do
      end do
   end subroutine post

   subroutine gidres_header()
      character(len=10) :: element_type_solid
      character(len=13) :: element_type_joint

      select case( param%nods )
      case( 8, 20 )
         element_type_solid = "Hexahedra"
      case( 4, 10 )
         element_type_solid = "Tetrahedra"
      case default
         continue
      end select

      write(ounit_postres,*) "GiD Post Results File 1.0"
      write(ounit_postres,100) element_type_solid
100   format('GaussPoints "Board elements" ElemType',1x,a10)
      write(ounit_postres,*) "Number Of Gauss Points: 1"
      write(ounit_postres,*) "Natural Coordinates: internal"
      write(ounit_postres,*) "end gausspoints"

      write(ounit_postres,110) element_type_solid
110   format('GaussPoints "NST" ElemType',1x,a10)
      write(ounit_postres,*) "Number Of Gauss Points: ", param%nips_stress ! nips
      write(ounit_postres,*) "Natural Coordinates: internal"
      write(ounit_postres,*) "end gausspoints"

      if( param%nelej==0 ) return
      
      select case( param%nodj )
      case( 8, 16 )
         element_type_joint = "Quadrilateral"
      case( 6, 12 )
         element_type_joint = "Triangle"
      case default
         continue
      end select

      write(ounit_postres,120)element_type_joint
120   format('GaussPoints "JGP" ElemType',1x,a13)
      write(ounit_postres,*) "Number Of Gauss Points: ", param%nipb
      write(ounit_postres,*) "Natural Coordinates: Internal"
      write(ounit_postres,*) "end gausspoints"

      write(ounit_postres,140) element_type_joint   ! for seismic moment
140   format('GaussPoints "JJJJ" ElemType',1x,a13)
      write(ounit_postres,*) "Number Of Gauss Points: 1"
      write(ounit_postres,*) "Natural Coordinates: internal"
      write(ounit_postres,*) "end gausspoints"
   end subroutine gidres_header

   subroutine gidmesh
      integer :: i, iel, et, ilayer, layertype(param%nlayer)
      character(len=10) :: element_type_solid
      character(len=13) :: element_type_joint

      select case( param%nods )
      case( 8, 20 )
         element_type_solid = "Hexahedra"
      case( 4, 10 )
         element_type_solid = "Tetrahedra"
      case default
         continue
      end select

      if( param%nelej/=0 ) then
         select case( param%nodj )
         case( 8, 16 )
            element_type_joint = "Quadrilateral"
         case( 6, 12 )
            element_type_joint = "Triangle"
         case default
            stop
         end select
      end if

      layertype(:) = 1 ! initialze as solid layer.
      do ilayer=1, param%nlayer
         ! if any joint element has this layer
         if( any( geom%elem_layer_namej(:) == param%layer_name(ilayer) ) ) then
            layertype(ilayer) = 0 ! this is joint layer.
         end if
      end do

      do ilayer=1, param%nlayer
         if( layertype(ilayer)==1 ) then ! solid layer
            write(ounit_postmesh,10) trim(param%layer_name(ilayer)), element_type_solid, param%nods
10          format('MESH "',a,'" dimension 3 ElemType ',a,' Nnode',i3)
         else ! joint layer
            write(ounit_postmesh,20) trim(param%layer_name(ilayer)), element_type_joint, param%nodj/2
20          format('MESH "',a,'" dimension 3 ElemType ',a,' Nnode',i3)
         end if
         
         if( ilayer==1 ) then
            write(ounit_postmesh,*) "Coordinates"
            do i=1, param%nnode
               write(ounit_postmesh,'(1x,i10,3(e25.15))') i, geom%coord(:,i)
            end do
            write(ounit_postmesh,*) "end coordinates"
         else
            write(ounit_postmesh,*) "Coordinates"
            write(ounit_postmesh,*) "# nodal coordinates already given in the previous table"
            write(ounit_postmesh,*) "end coordinates"
         end if

         write(ounit_postmesh,*) "Elements"
         do iel=1, param%neles
            ! only solid elements of this layer.
            if( geom%elem_layer_names(iel) == param%layer_name(ilayer) ) then
               write(ounit_postmesh,'(1x,i10,21(i10))') iel, geom%enodes(:,iel)  ! Follow GiD node connectivities
            end if
         end do

         do iel=1, param%nelej
            ! only joint elements of this layer.
            if( geom%elem_layer_namej(iel) == param%layer_name(ilayer) ) then
               write(ounit_postmesh,'(1x,i10,21(i10))') iel, geom%enodej(1:param%nodj/2,iel)
            end if
         end do
         write(ounit_postmesh,*) "end elements"
      end do
   end subroutine gidmesh

   subroutine gidres( time )
      implicit none
      real(8), intent(in) :: time

      integer :: n, iel, ip, i
      real(8) :: points(ndim,param%nips_stress), weights(param%nips_stress)
      real(8) :: funs(param%nods,param%nips_stress)
      real(8) :: ders(param%nods,ndim,param%nips_stress)
      type(propertys), pointer :: prop
      real(8) :: Dmat(nst,nst)
      real(8) :: Bmat(nst,ndim*param%nods)
      real(8) :: coords(ndim,param%nods), displs(ndim*param%nods)
      real(8) :: jac(ndim,ndim), invjac(ndim,ndim), det
      real(8) :: gfun(param%nods,ndim)
      real(8) :: deps(nst), dsig(nst), sig(nst), eps(nst)
      real(8) :: pointj(2,param%nipb), weightj(param%nipb) ! MIZUHO fixed
      real(8) :: coordj(ndim,param%nodj/2)
      real(8) :: funj(param%nodj/2,param%nipb)
      real(8) :: gap(ndim), gapmax

      character(len= 3), parameter :: length="[L]"
      character(len=12), parameter :: line="------------"
      character(len= 8), parameter :: psi="[M/L/T2]"

      integer, allocatable :: joint_node(:)

      ! Displacement components
      write(ounit_postres,930) time
930   format('Result "SOLID_displacement" "Elast3D_Analysis " ',f6.2,' Vector OnNodes ')
      write(ounit_postres,*) 'ComponentNames "u","v","w"'
      write(ounit_postres,*) 'Unit "m"'
      write(ounit_postres,*) 'Values'
      do n=1, param%nnode
         write(ounit_postres,'(I9,3e12.4)') n, geom%displ(1:3,n)
      end do
      write(ounit_postres,*) 'End Values'

      ! Pore pressure
      write(ounit_postres,932) time
932   format('Result "SOLID_pore pressure" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "Board elements"')
      write(ounit_postres,*) 'ComponentNames "po"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%neles
         write(ounit_postres,'(I9,e12.4)') iel, presss(iel)
      end do
      write(ounit_postres,*) 'End Values'

      ! Density change
      write(ounit_postres,832) time
832   format('Result "SOLID_density change" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "Board elements"')
      write(ounit_postres,*) 'ComponentNames "drho"'
      write(ounit_postres,*) 'Unit "kg/m3"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%neles
         write(ounit_postres,'(I9,e12.4)') iel, density_dif(iel)
      end do
      write(ounit_postres,*) 'End Values'
      
      call samples( param%nips_stress, points, weights )
      do ip=1, param%nips_stress
         call shape_funs( funs(:,ip),   points(:,ip) ) ! Nm at point
         call shape_ders( ders(:,:,ip), points(:,ip) ) ! dNm/d(xi,et,ze)
      end do
      
      ! Stress components
      write(ounit_postres,933) time
933   format('Result "SOLID_stress" "Elast3D_Analysis " ',f6.2,' Matrix OnGausspoints "NST"')
      write(ounit_postres,*) 'ComponentNames "sx","sy","sz","txy","tyz","tzx"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         call deemats( Dmat, prop%E, prop%nu ) ! D

         coords(:,:) = geom%coord(:,geom%enodes(:,iel))
         displs(:) = reshape( geom%displ(:,geom%enodes(:,iel)), shape=[ndim*param%nods] )

         do ip=1, param%nips_stress
            jac(:,:) = matmul(coords,ders(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )              ! B=L*Nm

            deps(:) = matmul(Bmat,displs)
            dsig(:) = matmul(Dmat,deps)

            sig(:) = dsig(:) + stresss_ini(:,iel)

            if( ip==1 ) then
               write(ounit_postres,'(I9,6e14.6)') iel, sig(:)
            else
               write(ounit_postres,'(9x,6e14.6)')      sig(:)
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Stress changes
      write(ounit_postres,733) time
733   format('Result "SOLID_stress change" "Elast3D_Analysis " ',f6.2,' Matrix OnGausspoints "NST"')
      write(ounit_postres,*) 'ComponentNames "dsx","dsy","dsz","dtxy","dtyz","dtzx"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'

      do iel=1, param%neles
         prop => param%props(geom%etypes(iel))
         call deemats( Dmat, prop%E, prop%nu ) ! D

         coords(:,:) = geom%coord(:,geom%enodes(:,iel))
         displs(:) = reshape( geom%displ(:,geom%enodes(:,iel)), shape=[ndim*param%nods] ) 

         do ip=1, param%nips_stress
            jac(:,:) = matmul(coords,ders(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )              ! B=L*Nm

            deps(:) = matmul(Bmat,displs)
            dsig(:) = matmul(Dmat,deps)

            if( ip==1 ) then
               write(ounit_postres,'(I9,6e14.6)')  iel, dsig(:)
            else
               write(ounit_postres,'(9x,6e14.6)')       dsig(:)
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Strain components
      write(ounit_postres,934) time
934   format('Result "SOLID_strain" "Elast3D_Analysis " ',f6.2,' Matrix OnGausspoints "NST"')
      write(ounit_postres,*) 'ComponentNames "ex","ey","ez","exy","eyz","ezx"'
      write(ounit_postres,*) 'Unit "-"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%neles
         coords(:,:) = geom%coord(:,geom%enodes(:,iel))
         displs(:) = reshape( geom%displ(:,geom%enodes(:,iel)), shape=[ndim*param%nods] )

         do ip=1, param%nips_stress
            jac(:,:) = matmul(coords,ders(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )              ! B=L*Nm

            deps(:) = matmul(Bmat,displs)
            eps(:) = strains_ini(:,iel) + deps(:)

            if( ip==1 ) then
               write(ounit_postres,'(I9,6e14.6)') iel, eps(:)
            else
               write(ounit_postres,'(9x,6e14.6)')      eps(:)
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Strain changes
      write(ounit_postres,734) time
734   format('Result "SOLID_strain change" "Elast3D_Analysis " ',f6.2,' Matrix OnGausspoints "NST"')
      write(ounit_postres,*) 'ComponentNames "dex","dey","dez","dexy","deyz","dezx"'
      write(ounit_postres,*) 'Unit "-"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%neles
         coords(:,:) = geom%coord(:,geom%enodes(:,iel))
         displs(:) = reshape( geom%displ(:,geom%enodes(:,iel)), shape=[ndim*param%nods] )

         do ip=1, param%nips_stress
            jac(:,:) = matmul(coords,ders(:,:,ip))  !  d(x,y,z)/d(xi,et,ze)

            det = determinants( jac )               ! |d(x,y,z)/d(xi,et,ze)|

            invjac = invert( jac )                  ! d(xi,et,ze)/d(x,y,z)

            gfun(:,:) = matmul(ders(:,:,ip),invjac) ! dNm/d(x,y,z)

            call beemats( Bmat, gfun )              ! B=L*Nm

            deps(:) = matmul(Bmat,displs)

            if( ip==1 ) then
               write(ounit_postres,'(I9,6e14.6)') iel, deps(:)
            else
               write(ounit_postres,'(9x,6e14.6)')      deps(:)
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Strain energy density
      write(ounit_postres,334) time
334   format('Result "SOLID_strain_energy_density" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "NST"')
      write(ounit_postres,*) 'ComponentNames "dse"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%neles
         do ip=1, param%nips_stress
            if( ip==1 ) then
               write(ounit_postres,'(I9,E14.6)') iel, energy(ip,iel)
            else
               write(ounit_postres,'(9x,E14.6)')      energy(ip,iel)     
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Formation robustness index
      write(ounit_postres,335) time
335   format('Result "SOLID_robustness_index" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "NST"')
      write(ounit_postres,*) 'ComponentNames "lfr"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%neles
         do ip=1, param%nips_stress
            if( ip==1 ) then
               write(ounit_postres,'(I9,e14.6)') iel, robust(ip,iel)
            else
               write(ounit_postres,'(9x,e14.6)')      robust(ip,iel)     
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Formation robustness index change
      write(ounit_postres,336) time
336   format('Result "SOLID_robustness_index_change" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "NST"')
      write(ounit_postres,*) 'ComponentNames "dlfr"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%neles
         do ip=1, param%nips_stress
            if( ip==1 ) then
               write(ounit_postres,'(I9,e14.6)') iel, robust_dif(ip,iel)
            else
               write(ounit_postres,'(9x,e14.6)')      robust_dif(ip,iel)     
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      if( param%nelej==0 ) return

      ! Traction
      write(ounit_postres,935) time
935   format('Result "JOINT_traction" "Elast3D_Analysis " ',f6.2,' Vector OnGausspoints "JGP"')    
      write(ounit_postres,*) 'ComponentNames "t1","t2","t3"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         do ip=1, param%nipb
            if( ip==1 ) then
               write(ounit_postres,'(I9,3E12.4)') iel, tractionj(:,ip,iel)
            else
               write(ounit_postres,'(9x,3E12.4)')      tractionj(:,ip,iel)
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Relative displacement
      write(ounit_postres,936) time
936   format('Result "JOINT_rel_displacement" "Elast3D_Analysis " ',f6.2,' Vector OnNodes')
      write(ounit_postres,*) 'ComponentNames "d1","d2","d3"'
      write(ounit_postres,*) 'Unit "m"'
      write(ounit_postres,*) 'Values'
      write(ounit_res,'(/"*** JOINT DISPLACEMENT ***")')
      write(ounit_res,'("    JOINT_NODE         X-DISP        Y-DISP        Z-DISP      ABS(DISP)")')
      write(ounit_res,'(7(6x,a3,5x))') (length,i=1,4)
      write(ounit_res,'(1x,7(1x,a12,1x))') (line,i=1,4)
      gapmax=0.0d0

      ! make pair list of jointed nodes.
      if( param%nelej>0 ) then
         allocate(joint_node(param%nnode)); joint_node(:)=0

         do iel=1, param%nelej
            joint_node(geom%enodej(1:param%nodj/2,iel)) = geom%enodej(param%nodj/2+1:param%nodj,iel)
         end do

         do i=1, param%nnode
            if( joint_node(i)==0 ) cycle

            gap(:) = geom%displ(:,i) - geom%displ(:,joint_node(i))
            gapmax = max(gapmax,norm(gap))

            write(ounit_postres,'(1x,i9,3E14.6)') i, gap(:)
            write(ounit_res,'(i9,4es14.6)')       i, gap(:), norm(gap)
         end do

         deallocate(joint_node)
      end if
      
      write(ounit_postres,*) 'End Values'

      write(  *,      '(/"max crack displacement = ",es15.6," [L]")') gapmax
      write(ounit_res,'(/"max crack displacement = ",es15.6," [L]")') gapmax

      ! Pore pressure
      write(ounit_postres,938) time
938   format('Result "JOINT_pore_press" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "JJJJ"')
      write(ounit_postres,*) 'ComponentNames "joint_po"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         write(ounit_postres,'(i9,e12.4)') iel, pressj(iel)
      end do
      write(ounit_postres,*) 'End Values'

      ! Coulomb stress
      write(ounit_postres,941) time
941   format('Result "JOINT_coulomb stress" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "JGP"')
      write(ounit_postres,*) 'ComponentNames "CFS"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         do ip=1, param%nipb
            if( ip==1 ) then
               write(ounit_postres,'(i9,e14.6)') iel, coulombj(ip,iel)
            else
               write(ounit_postres,'(9x,e14.6)')      coulombj(ip,iel)
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Coulomb stress change
      write(ounit_postres,942) time
942   format('Result "JOINT_coulomb stress change" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "JGP"')
      write(ounit_postres,*) 'ComponentNames "delta_CFS"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         do ip=1, param%nipb
            if( ip==1 ) then
               write(ounit_postres,'(i9,e14.6)') iel, coulombj_dif(ip,iel)
            else
               write(ounit_postres,'(9x,e14.6)')      coulombj_dif(ip,iel)
            end if
         end do
      end do
      write(ounit_postres,*)'End Values'

      ! Cohesion
      write(ounit_postres,442) time
442   format('Result "JOINT_cohesion" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "JGP"')
      write(ounit_postres,*) 'ComponentNames "Co"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         do ip=1, param%nipb
            if( ip==1 ) then
               write(ounit_postres,'(i9,e14.6)') iel, cohesionj(ip,iel)
            else
               write(ounit_postres,'(9x,e14.6)')      cohesionj(ip,iel)
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Tensile strength
      write(ounit_postres,443) time
443   format('Result "JOINT_tensile_strength" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "JGP"')
      write(ounit_postres,*) 'ComponentNames "tens"'
      write(ounit_postres,*) 'Unit "Pa"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         do ip=1, param%nipb
            if( ip==1 ) then
               write(ounit_postres,'(i9,e14.6)') iel, tensilej(ip,iel)
            else
               write(ounit_postres,'(9x,e14.6)')      tensilej(ip,iel)
            end if
         end do
      end do
      write(ounit_postres,*) 'End Values'

      ! Seismic moment
      write(ounit_postres,943) time
943   format('Result "JOINT_seismic moment" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "JJJJ"')
      write(ounit_postres,*) 'ComponentNames "Ms"'
      write(ounit_postres,*) 'Unit "N-m"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         write(ounit_postres,'(I9,E12.4)') iel, seismic_mom(iel)
      end do
      write(ounit_postres,*) 'End Values'

      ! Failure index
      write(ounit_postres,944) time
944   format('Result "JOINT_failure_index" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "JJJJ"')
      write(ounit_postres,*) 'ComponentNames "FAILURE"'
      write(ounit_postres,*) 'Unit "-"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         write(ounit_postres,'(I9,I3)') iel, MERGE(1,0,failurej(iel))
      end do
      write(ounit_postres,*) 'End Values'

      ! Seismic moment intensity
      write(ounit_postres,945) time
945   format('Result "JOINT_seismic moment_intensity" "Elast3D_Analysis " ',f6.2,' Scalar OnGausspoints "JJJJ"')
      write(ounit_postres,*) 'ComponentNames "Ms_int"'
      write(ounit_postres,*) 'Unit "N/m"'
      write(ounit_postres,*) 'Values'
      do iel=1, param%nelej
         write(ounit_postres,'(I9,E12.4)') iel, seismic_int(iel)
      end do
      write(ounit_postres,*) 'End Values'

      if( param%nelej>0 ) then
         write(ounit_res,'(/"*** JOINT TRACTION ***")')
         write(ounit_res,'(" ELEMENT  GP         X             Y             Z           TS1           TS2           TN")')
         write(ounit_res,'(14x,3(6x,a3,5x),3(3x,a8,3x))') (length,i=1,3), (psi,i=1,3)
         write(ounit_res,'(" -------  ---  ",6(1x,a12,1x))') (line,i=1,6)

         call samplej( param%nipb, pointj, weightj ) ! MIZUHO fixed
         do ip=1, param%nipb
            call shape_funj( funj(:,ip), pointj(:,ip) )
         end do
         
         do iel=1, param%nelej
            coordj(:,:) = geom%coord(:,geom%enodej(1:param%nodj/2,iel))

            do ip=1, param%nipb
               write(ounit_res,'(i9,4x,i2,2x,6es14.6)') &
                  iel, ip, matmul(coordj,funj(:,ip)), tractionj(:,ip,iel)
            end do
         end do
      end if
   end subroutine gidres

   subroutine freeArray
      deallocate( param%props )
      if( param%nelej>0 ) deallocate( param%propj )
      deallocate( param%layer_name ) !, param%layer_names, param%layer_namej )
      deallocate( geom%coord, geom%displ, geom%displ_inc )
      deallocate( geom%enodes, geom%etypes, geom%elem_layer_names )
      if( param%nelej>0 ) deallocate( geom%enodej, geom%etypej, geom%elem_layer_namej )
      deallocate( geom%ecenters, geom%ecenterj )
      deallocate( geom%nodfree, geom%valfix ) ! , geom%codfree

      deallocate( stresss_ini )
      deallocate( presss_ini, presss )
      if( param%nelej>0 ) deallocate( stressj_ini )
      if( param%nelej>0 ) deallocate( pressj_ini, pressj )
      deallocate( temp_ini, temp, density_dif )

      if( allocated(energy) ) then
         deallocate( strains_ini )
         if( param%nelej>0 ) deallocate( cohesionj, tensilej, stiffj )
         if( param%nelej>0 ) deallocate( tractionj_ini, tractionj )
         if( param%nelej>0 ) deallocate( coulombj_ini, coulombj, coulombj_dif )
         if( param%nelej>0 ) deallocate( failurej )
         deallocate( energy )
         deallocate( robust, robust_ini, robust_dif )
         deallocate( load )
         if( param%nelej>0 ) deallocate( seismic_mom, seismic_int )
         deallocate( geom%esteers )
         if( param%nelej>0 ) deallocate( geom%esteerj )
         deallocate( geom%ia, geom%ja )
#ifdef ENABLE_PETSC
         if( solver == "PETSc" ) deallocate( geom%irn, geom%jcn, geom%idx )
#endif
         deallocate( stiffmat )
         deallocate( geom%eindexs )
         if( param%nelej>0 ) deallocate( geom%eindexj )
      end if
      
      ! Close files
      close(iunit_point_load)
      close(iunit_surfj_load)
      close(iunit_presss_load)
      close(iunit_pressj_load)
      close(ounit_log)
      close(ounit_res)
      close(ounit_postmesh)
      close(ounit_postres)
   end subroutine freeArray
end module wFAST

program main
#ifdef ENABLE_PETSC
   use petscsnes
   use petscksp
#endif
   use wFAST
   implicit none

   character(len=256) :: filename_cfg
   integer :: itime, iteration
   real(8) :: time, residual
#ifdef ENABLE_PETSC
   integer :: error
#endif
   
   !--- Get filename of input or config file.
   if( command_argument_count()/=1 ) then
      write(*,*) "Error: no config file is given."
      call exit(1)
   end if
   call get_command_argument( 1, filename_cfg )

   call loadConfig( filename_cfg )

#ifdef ENABLE_PETSC
   if( solver == "PETSc" ) call PetscInitialize( error )
#endif

   if( filename_gid /= "" ) then
      call loadGiD
      call saveDAT
      call freeArray
      call loadDAT
   else if( filename_dat /= "" ) then
      call loadDAT
   else
      write(*,*) "Error: no GiD or DAT file is given."
      call exit(1)
   end if
   
   write(ounit_log,*) "wFAST: Finished reading initial conditions."

   write(ounit_log,*) "WFAST: Allocate large arrays ok ..."

   call calcInitStrain
   
   !
   !--- (4) Create a pre-calculation files ---
   !
   call initial_report()

   ! GiD files
   call gidmesh()
   call gidres_header()

   write(ounit_log,*) "WFAST: Create GiD mesh file ok ..."

   call calcInitTraction
   call evalInitCoulomb
      
   write(ounit_log,*) "WFAST: Initializing parameters ok ..."

   call evalInitFail

   !
   !---(9) Post-processing for initial conditions
   !
   call evalInitRobust

   write(*,*) "Print initial results ..."
   call post
   write(*,*) "Print GiD results ..."
   time = 0.0d0
   call gidres( time )
   write(ounit_log,*) "WFAST: Print initial results ok ..."
   write(*,*) "Print initial results completed ..."

   write(*,*)
   write(*,*) " => JOINTS ARE STABLE UNDER INITIAL CONDITIONS. CALCULATION CONTINUES..."
   write(*,*)
   write(ounit_res,*)
   write(ounit_res,*) " => JOINTS ARE STABLE UNDER INITIAL CONDITIONS. CALCULATION CONTINUES..."
   write(ounit_res,*)
   write(ounit_res,*)
   write(ounit_res,*)
   write(ounit_res,'("*****************************")')
   write(ounit_res,'("*** (3) LOADED CONDITIONS ***")')
   write(ounit_res,'("*****************************")')

   !-----------------------
   !--- Start load step ---
   !-----------------------

   do itime=1, param%ntime
      call setupLoad
      call setupStiff

      !--------------------------------
      !---(XX) Begin time-step loop ---
      !--------------------------------
      write(  *,'(" => Nonlinear calculaiton has started...")')
      write(ounit_log,'(" => Nonlinear calculaiton has started...")')

      do iteration=1, param%niter
         write(*,'("itime=",i4)') iteration

#ifdef ENABLE_PARDISO
         if( solver == "PARDISO" ) call solveDispl_PARDISO
#endif
#ifdef ENABLE_PETSC
         if( solver == "PETSc"   ) call solveDispl_PETSc
#endif
    
         
         !------------------------------------------------
         !---(XX) Compute residual & check convergence ---
         !------------------------------------------------
         call evalFail( residual )
         call progress_report( iteration, residual )

         if( residual<param%convergence ) exit

         call setupStiff
      end do ! iteration

      if( iteration<=param%niter ) then
         write(  *,'(" => Non-linear calculation successfully converged... ")')
         write(ounit_res,'(" => Non-linear calculation successfully converged... ")')
      else
         write(  *,'(" => Warning: Calculation did not converge... ")')
         write(ounit_res,'(" => Warning: Calculation did not converge... ")')
      end if

      !------------------------------------------------------------
      !---(XX) Post-processing for stress & strain calculations ---
      !------------------------------------------------------------
      call evalSeismic
      call evalCoulomb
      call evalRobust

      !--------------------------------------------
      !---(XX) Print simulation results for GiD ---
      !--------------------------------------------
      call post
      time = dble(itime)
      call gidres( time )
      write(*,*) 'GiD result file is printed...'

   end do ! itime

   call freeArray
   
#ifdef ENABLE_PETSC
   if( solver == "PETSc" ) call PetscFinalize( error )
#endif
end program main