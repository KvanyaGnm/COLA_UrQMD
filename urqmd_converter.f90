
module cola_fortran_converter_impl
  use cola
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public, extends(AbstractFortranConverter) :: URQMDConverter
    character(len=512) :: config_file = ''
    real(8)            :: propagate_time = 0.0 
    logical            :: pass_coll = .true.  
  contains
    procedure :: init => converter_init
    procedure :: convert => converter_convert
    final :: converter_final
  end type URQMDConverter

  interface

    subroutine urqmd_cola_run_propagation(prop_time, parts_in, parts_out, coll)
      real(8), intent(in) :: prop_time
      logical, intent(in) :: coll
      type(EventData), intent(in) :: parts_in
      type(EventData), intent(out) :: parts_out 
    end subroutine

  contains

  subroutine converter_init(self, pmap, err)
    class(URQMDConverter), intent(inout) :: self
    type(ParametersMap), intent(in) :: pmap
    character(len=:), allocatable, intent(out) :: err
    type(ParametersMapItem) :: kv
    character(len=:), allocatable :: key, val
    integer :: i, n

    err = ''
    self%propagate_time = 0.0
    self%pass_coll = .true.
    n = pmap%size()

    do i = 1, n
      kv = pmap%get(i)
      key = kv%get_first()
      val = kv%get_second()

      if (key == 'config_file') then
        self%config_file = trim(val)
      else if (key == 'propagate_time') then
        read(val, *, iostat=i) self%propagate_time
      else if (key == 'coll_to_r0') then
        self%pass_coll = (trim(val) == 'true')
      end if
    end do

    err = ''
  end subroutine converter_init

  function converter_convert(self, data_in, err) result(data_out)
    class(URQMDConverter), intent(in) :: self
    type(EventData), intent(in) :: data_in
    character(len=:), allocatable, intent(out) :: err
    type(EventData) :: data_out
    type(EventParticles) :: parts_in, parts_out
    integer :: n_parts
    character(len=256) :: local_err

    data_out = data_in  ! Копируем начальные данные (метаинформацию события)
    parts_in = data_in%get_particles()
    err = ''    

    call urqmd_cola_run_propagation(self%propagate_time, parts_in, parts_out, self%pass_coll)
    call data_out%set_particles(parts_out)
  end function converter_convert

  subroutine converter_final(self)
    type(URQMDConverter), intent(inout) :: self
  end subroutine converter_final

end module cola_fortran_converter_impl