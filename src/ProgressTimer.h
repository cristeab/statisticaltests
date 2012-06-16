/*
This class is used for displaying simulation progress and remaining time
Author: Cristea Bogdan
Revision date: 08.09.2009
Including Progress_Timer class into tr namespace, completely independent from IT++
Revision date: 27.08.2007

Usage example:

#include "Progress_Timer.h"

#define MAX ((int)5e6)

int main(void)
{
  tr::Progress_Timer pt;
  pt.set_max(MAX);
  pt.progress(0);
  for(int n=0;n<MAX;n++)
  {
    pt.progress(n+1);
  }
  pt.toc_print();
}
 */

#ifndef PROGRESS_TIMER_H_
#define PROGRESS_TIMER_H_

#include <iostream>
#include <iomanip>

namespace tr
{

class Progress_Timer //shows simulation progress and remaining time
{
public:
	//! Creates a new timer. Sets max_it to zero.
	Progress_Timer()
	{
		max_it = 0.0;
	};
	//! Shows progress and remaining time using real numbers as input (percent should vary from 0 to 1)
	void progress(const double &percent = 0.0)
	{
	    if(percent==0.0)
	    {
	      //reset and start timer
	      time(&start_time);
	      stop_time = 0;
	      elapsed_time = 0;
	      //display progress
	      std::cout << "0.0%\n";
	      fflush(stdout);//useful when standard output is directed to a file
	    }
	    else
	    {
	      time(&stop_time);
	      if(((stop_time-start_time)>=1.0)||(percent==1.0))//show progress only if the time difference since last display is greater than 1s or we are at the end of our computations
	      {
	    	elapsed_time += stop_time-start_time;//update elapsed time
	    	start_time = stop_time;//update start time
	    	//display progress and remaining time
	        std::cout << std::fixed << std::setw(3) << std::setprecision(1) << 100.0*percent << "%\t";
	    	sec2human(double(elapsed_time)/percent-double(elapsed_time));
	    	std::cout << " remaining\n";
	    	fflush(stdout);
	      }
	    }
	};
	//! Sets maximum number of iterations
	void set_max(const int max)//set maximum number of iterations
	{
		max_it=(double)max;
	};
	//! Shows progress and remaining time using integer numbers as input (iteration should vary from 0 to max_it)
	void progress(const int &iteration = 0)
	{
		progress((double(iteration))/max_it);
	};
	//! Prints the elapsed time in a human readable format
	void toc_print(void)
	{
		std::cout << "Elapsed time = ";
		sec2human(double(elapsed_time));
		std::cout << std::endl;
	};
private:
	double max_it;//maximum number of iterations
	time_t start_time,stop_time,elapsed_time;
	//! Converts seconds to a human readable format and display to standard output
	void sec2human(double sec)//converts seconds to human readable format as string
	{
		int s = (int)(sec+0.5);//round to nearest integer
		int d = s/86400;
		s -= d*86400;
		int h = s/3600;
		s -= h*3600;
		int m = s/60;
		s -= m*60;
		if(d>0)
		{
			std::cout << d << " day, " << h << " hr, " << m << " min, " << s << " sec";
		}
		else if(h>0)
		{
			std::cout << h << " hr, " << m << " min, " << s << " sec";
		}
		else if(m>0)
		{
			std::cout << m << " min, " << s << " sec";
		}
		else
		{
			std::cout << s << " sec";
		}
	};
};

}//namespace tr

#endif
