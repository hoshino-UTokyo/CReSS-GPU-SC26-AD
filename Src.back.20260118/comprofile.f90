!***********************************************************************
! Module for profiling OpenMP parallel sections
! Records execution count, loop length, and timing information
!***********************************************************************

      module m_comprofile

      use omp_lib

      implicit none

!-----7----------------------------------------------------------7-----

! Maximum number of profiled sections
      integer, parameter :: max_sections = 500

! Section information
      character(len=64), save :: prof_file(max_sections)
      character(len=64), save :: prof_subr(max_sections)
      character(len=128), save :: prof_desc(max_sections)

! Profiling data
      integer(8), save :: prof_count(max_sections)      ! Execution count
      integer(8), save :: prof_loops(max_sections)      ! Total loop iterations
      real(8), save :: prof_time(max_sections)          ! Total time (seconds)
      real(8), save :: prof_start(max_sections)         ! Start time (temp)

! Section registration
      integer, save :: prof_nsections = 0
      logical, save :: prof_initialized = .false.

!$omp threadprivate(prof_start)

      contains

!***********************************************************************
! Initialize profiling system
!***********************************************************************

      subroutine profile_init()

      implicit none

      integer :: i

!-----7----------------------------------------------------------7-----

      if (prof_initialized) return

      do i = 1, max_sections
        prof_file(i) = ''
        prof_subr(i) = ''
        prof_desc(i) = ''
        prof_count(i) = 0
        prof_loops(i) = 0
        prof_time(i) = 0.0d0
        prof_start(i) = 0.0d0
      end do

      prof_nsections = 0
      prof_initialized = .true.

      end subroutine profile_init

!***********************************************************************
! Register a new section and return its ID
!***********************************************************************

      function profile_register(filename, subrname, description)        &
     &                          result(section_id)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: subrname
      character(len=*), intent(in) :: description

      integer :: section_id

!-----7----------------------------------------------------------7-----

      if (.not. prof_initialized) call profile_init()

!$omp critical(prof_register_lock)
      prof_nsections = prof_nsections + 1
      section_id = prof_nsections

      if (section_id <= max_sections) then
        prof_file(section_id) = trim(filename)
        prof_subr(section_id) = trim(subrname)
        prof_desc(section_id) = trim(description)
      end if
!$omp end critical(prof_register_lock)

      end function profile_register

!***********************************************************************
! Start timing a section
!***********************************************************************

      subroutine profile_start(section_id)

      implicit none

      integer, intent(in) :: section_id

!-----7----------------------------------------------------------7-----

      if (section_id < 1 .or. section_id > max_sections) return

      prof_start(section_id) = omp_get_wtime()

      end subroutine profile_start

!***********************************************************************
! Stop timing a section and record statistics
!***********************************************************************

      subroutine profile_stop(section_id, loop_length)

      implicit none

      integer, intent(in) :: section_id
      integer(8), intent(in) :: loop_length

      real(8) :: elapsed

!-----7----------------------------------------------------------7-----

      if (section_id < 1 .or. section_id > max_sections) return

      elapsed = omp_get_wtime() - prof_start(section_id)

!$omp atomic
      prof_count(section_id) = prof_count(section_id) + 1
!$omp atomic
      prof_loops(section_id) = prof_loops(section_id) + loop_length
!$omp atomic
      prof_time(section_id) = prof_time(section_id) + elapsed

      end subroutine profile_stop

!***********************************************************************
! Output profiling results
!***********************************************************************

      subroutine profile_output(unit_num)

      implicit none

      integer, intent(in) :: unit_num

      integer :: i
      real(8) :: avg_time, avg_loops
      character(len=256) :: fmt_header, fmt_data

!-----7----------------------------------------------------------7-----

      if (.not. prof_initialized) return
      if (prof_nsections == 0) return

      write(unit_num, '(A)')                                            &
     &  '========================================================'//    &
     &  '========================================================'
      write(unit_num, '(A)') 'OpenMP Parallel Section Profiling Results'
      write(unit_num, '(A)')                                            &
     &  '========================================================'//    &
     &  '========================================================'
      write(unit_num, '(A)')

      write(unit_num, '(A5,2X,A20,2X,A20,2X,A12,2X,A15,2X,A15,2X,A15)')  &
     &  'ID', 'File', 'Subroutine', 'Count',                            &
     &  'Avg Loops', 'Total Time(s)', 'Avg Time(ms)'
      write(unit_num, '(A)')                                            &
     &  '--------------------------------------------------------'//    &
     &  '--------------------------------------------------------'

      do i = 1, min(prof_nsections, max_sections)
        if (prof_count(i) > 0) then
          avg_time = prof_time(i) / dble(prof_count(i)) * 1000.0d0
          avg_loops = dble(prof_loops(i)) / dble(prof_count(i))

          write(unit_num,                                               &
     &      '(I5,2X,A20,2X,A20,2X,I12,2X,F15.1,2X,F15.6,2X,F15.6)')      &
     &      i, prof_file(i)(1:20), prof_subr(i)(1:20),                  &
     &      prof_count(i), avg_loops, prof_time(i), avg_time
        end if
      end do

      write(unit_num, '(A)')
      write(unit_num, '(A)')                                            &
     &  '========================================================'//    &
     &  '========================================================'

      end subroutine profile_output

!***********************************************************************
! Output profiling results to file
!***********************************************************************

      subroutine profile_output_file(filename)

      implicit none

      character(len=*), intent(in) :: filename

      integer :: iunit

!-----7----------------------------------------------------------7-----

      iunit = 99

      open(unit=iunit, file=trim(filename), status='replace',           &
     &     form='formatted')

      call profile_output(iunit)

      close(iunit)

      end subroutine profile_output_file

!***********************************************************************
! Finalize profiling and output results
!***********************************************************************

      subroutine profile_finalize()

      implicit none

!-----7----------------------------------------------------------7-----

      call profile_output(6)
      call profile_output_file('omp_profile.txt')

      end subroutine profile_finalize

      end module m_comprofile
