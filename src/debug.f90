module debug_utils

  implicit none

  logical, public, save :: debug_enableTimeStamps = .false.

contains

subroutine echoTimeStamp(timeString)

  character(len=*), intent(in) :: timeString

  integer, parameter :: nEcho = 20
  character(len=nEcho) :: infoString
  double precision :: currTime, prevTime = 0.0d0
  integer :: nChar

  if(.not. debug_enableTimeStamps) then
    return
  end if

  nChar = len_trim(timeString)
  if(nChar < nEcho) then
    infoString = timeString//repeat(" ",nEcho-nChar)
  else
    infoString = timeString(1:nEcho)
  end if

  call cpu_time(currTime)
  write(*,"(a,2(f11.6,a))") " >>> CPU time at "//infoString//": ",currTime," sec [+",currTime-prevTime," sec]"
  prevTime = currTime

end subroutine echoTimeStamp

end module debug_utils
