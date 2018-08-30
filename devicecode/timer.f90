module timer
    implicit none
    save
    real :: first_call, last_call
    contains
    SUBROUTINE start_timer()
      implicit none
      call cpu_time(first_call)
      last_call=first_call
    end SUBROUTINE start_timer
    SUBROUTINE print_time_step(name)
     use compiletimeconstants
     implicit none
     character(len=*),intent(in) :: name
     real current_time, time1, time2
     if(usetimer==1) then
       call cpu_time(current_time)
       time1=current_time-first_call
       time2=current_time-last_call
       print*,"(",time1,",",time2,")",name
       last_call=current_time
     end if
    end SUBROUTINE print_time_step
end module timer
