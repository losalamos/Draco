      module CAR_CU_Interface_Class
          implicit none
          integer        :: self
          character * 60 :: file_name = '/home/bta/sun_scalar/bin/top_hat' // ACHAR(0)
          logical        :: verbose = .TRUE.
          integer        :: rtt_format

          private
          public :: construct

          type, public :: CAR_CU_Interface
              private
              integer  :: this
          end type CAR_CU_Interface 

          interface construct
              module procedure CAR_CU_Interface_construct
          end interface

          contains

              subroutine CAR_CU_Interface_construct(self)
                  type(CAR_CU_Interface), intent(inout) :: self
                  integer  :: bool_verbose = 0

                  if (verbose)  bool_verbose = 1
  
                  ! F90 converts names to lower case (have to use the same 
                  ! in C++)
                  call construct_car_cu_interface(self%this, file_name, &
                      bool_verbose, rtt_format)

              end subroutine CAR_CU_Interface_construct

       end module CAR_CU_Interface_Class
