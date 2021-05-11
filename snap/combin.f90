program combinations
 
  implicit none
  integer, parameter :: m_max = 4
  integer, parameter :: n_max = 256
  integer:: total,i
  integer, dimension (m_max) :: comb
  integer,dimension(n_max) :: info
  character (*), parameter :: fmt = '(i0' // repeat (', 1x, i0', m_max - 1) // ')'
  total=0
  do i=1,n_max
      info(i)=mod(i,6)
  enddo
 
  call gen (1)
  write(*,*)total
 
contains
 
  recursive subroutine gen (m)
 
    implicit none
    integer, intent (in) :: m
    integer :: n,tp
 
    if (m > m_max) then
        tp=0
        do n=1,m_max
          tp=tp+info(comb(n))
        enddo
        if(tp==6)then
        total=total+1
        endif
    else
      do n = 1, n_max
        if ((m == 1) .or. (n > comb (m - 1))) then
          comb (m) = n
          call gen (m + 1)
        end if
      end do
    end if
 
  end subroutine gen
 
end program combinations
