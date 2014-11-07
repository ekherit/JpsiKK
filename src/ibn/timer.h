#ifndef IBN_TIMER_H
#define IBN_TIMER_H
/*
 * =====================================================================================
 *
 *       Filename:  timer.h
 *
 *    Description:  Simple timer for measuring time
 *
 *        Version:  1.0
 *        Created:  02/06/2008 07:06:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */

#include <unistd.h>
#include <sys/time.h>


namespace ibn
{
  class timer
  {
    public:
      timer(void)  {restart();};

      void restart(void) 
      {
        gettimeofday(&begin_time,0);
      }


      double elapsed(void)
      {
        timeval current_time;
        gettimeofday(&current_time,0);
        return double(current_time.tv_sec - begin_time.tv_sec)+ double(current_time.tv_usec-begin_time.tv_usec)*1e-6;
      }

      timeval get_begin_time(void) const { return begin_time; }

      void set_begin_time(timeval t) 
      {
        begin_time=t;
      }

    private:
      timeval begin_time;	
  };

  class sleep
  {
    public:
      sleep(double time, bool volatile * flag, double precision=1e-3 /* sec */)
      {
        timeval begin_time;
        gettimeofday(&begin_time,0);
        double dt=0;
        do
        {
          usleep(useconds_t(precision*1e6));

          timeval current_time;
          gettimeofday(&current_time,0);
          dt=double(current_time.tv_sec - begin_time.tv_sec)+ double(current_time.tv_usec-begin_time.tv_usec)*1e-6;
        } while ( dt <= time && flag);
      }
      ~sleep(void){}
  };

};
#endif
