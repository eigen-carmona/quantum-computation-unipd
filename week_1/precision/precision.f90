program number_precision
implicit none
integer*2:: int_result_2
integer*4:: int_result_4
real*4:: single_result
real*8:: pi = acos(-1.d0), double_result

int_result_2 = 2000000+1
int_result_4 = 2000000+1

write(*,*) int_result_2, int_result_4

single_result = real(pi)*10.0**32+sqrt(2.0)*10.0**21
double_result = pi*10.d0**32+sqrt(2.d0)*10.d0**21

write(*,*) single_result, double_result

end program number_precision
