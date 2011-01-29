#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXF 100
int selected = 0;
FILE * gnuplot[MAXF];
char * command[MAXF];

/********************************************* 
gnuplotselect(channel)

if one wants to use multiple gnuplot windows, then 

  call gnuplot_select(1)
  call gnuplot_open("gnuplot")
  ... commands
  call gnuplot_flush()
  call gnuplot_select(0)
  call gnuplot_open("gnuplot")
  ... commands
  call gnuplot_flush()

default is 0, so if one needs just one window, do not use this command
**********************************************/

void gnuplotselect_(int *n) {
  if (*n >=0 && *n < MAXF) {
    selected = *n;
  }
}

/********************************************* 
gnuplotopen(command)

opens a pipe to program <command> (can be gnuplot or other)
**********************************************/

void gnuplotopen_(char * cmd, int n) {
  command[selected] = malloc((n+1)*sizeof(char));
  strncpy(command[selected], cmd, n);
  command[selected][n]='\0';
  gnuplot[selected] = popen(command[selected], "w");
}


/********************************************* 
gnuplotexecute(command)

send command to the selected pipe. 

notice that in fortran strings are weird, so use it like

  call gnuplot_execute("set xrange [0:1]")
 
but if you need to pass parameters, you need to add a '\0':
  
  character*100 :: str
  ...
  write (str, *) "set plot range [",xmin,",", xmax,"]", char(0)
  call gnuplot_execute(str)

**********************************************/

void gnuplotexecute_(char * cmd, int n) {
  char * str;
  str = malloc((n+1)*sizeof(char));
  strncpy(str, cmd, n);
  str[n]='\0';
/*  printf("%s\n", str); */
  fprintf(gnuplot[selected], "%s\n", str);
  free(str);
}

/********************************************* 
gnuplotwrite(buffer, size, count)

sends binary data from <buffer> (a pointer): 
<count> records of size <size>

examples:
 
- a generated image, 256 levels

  integer, parameter :: L=100
  integer*1 :: x(L, L)
  character*100 :: str
   
  call gnuplot_open("gnuplot") 
  write (str,*) "set xrange [0:",L,"]; set yrange [0:",L,"]", char(0)
  call gnuplot_execute(str)
  
  do t = ...
  
    .. compute x
  
    write (str, *) "plot '-' "// &
    " binary array=",L,"x",L," format='%uchar'"  // &
    " title 't=",t,"'  with image ", char(0)
    call gnuplot_execute(str)
    call gnuplot_write(C, 1,L*L)
    call gnuplot_flush()
  end do
 
 
- a sparse set of arrows (particle velocities)
    
   integer, parameter :: N = 100 
   real :: data (5,N) 

   ...
   
   write (str, *) "plot '-' "// &
   " binary record=",N," format='%float%float%float%float%*float'"  // &
   "' with vectors ", char(0)
   call gnuplot_execute(str)
   call gnuplot_write(data, 4,N*5)
   call gnuplot_flush()


NOTE: 
  array generates the coordinates of pixels or arrows, 
  record reads the coordinates from data
  gnuplot ignores the columns marked with "*" in format

NOTICE THAT ONE CAN ALSO USE POINTERS: 

   integer, parameter :: N = 100 
   real, target :: data (5,N)
   real, pointer :: plot 

   ...
   plot => data(1:4, N)
   
   write (str, *) "plot '-' "// &
   " binary record=",N," format='%float%float%float%float'"  // &
   "' with vectors ", char(0)
   call gnuplot_execute(str)
   call gnuplot_write(plot, 4,N*4)
   call gnuplot_flush()


**********************************************/


void gnuplotwrite_(void * c, int *size, int *n) {
  fwrite(c,*size,*n,gnuplot[selected]);
}

/**********************************************
call gnuplotflush()

empty buffer, so that gnuplot visualizes the data timely
**********************************************/


void gnuplotflush_() {
  fflush(gnuplot[selected]);
}

/**********************************************
call gnuplotclose()

closes channel
**********************************************/

void gnuplotclose_() {
  fclose(gnuplot[selected]);
  free(command[selected]);
}
