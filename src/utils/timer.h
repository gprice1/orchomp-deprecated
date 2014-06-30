
#ifndef _TIMER_H_
#define _TIMER_H_

#include <time.h>
#include <vector>
#include <iostream>

extern "C"
{
    #include "os.h"
}

class Timer{
  private:
    class single_timer{
        public:
        std::string name;
        struct timespec tic, toc;
        double total, elapsed;
        
        single_timer(std::string & name) : name( name ), total( 0.0 ){}
    };

    std::vector< single_timer > timers;

    int getIndex( std::string name ){

        for ( int i = 0; i < int(timers.size()); i ++ ){
            if( name == timers[i].name ){ return i; }
        }
        return -1;
    }
    bool timerExists( std::string & name, int index ){
        if (index < 0){
            std::cout << "No timer exists by the name: " 
                      << name << std::endl;
            return false;
        }
        return true;
    }



  public:
    void coutElapsed( std::string name ){
        int index = getIndex( name );
        if ( timerExists( name, index ) ){
            std::cout << "Timer [" << name << "] elapsed time: "
                      << timers[ index ].elapsed << "s\n";
        }
    }

    void coutTotal( std::string name ){
        int index = getIndex( name );
        if ( timerExists( name, index ) ){
            std::cout << "Timer [" << name << "] total time: "
                      << timers[ index ].total << "s\n";
        }
    }


    void start( std::string name ){
        int index = getIndex( name );

        if ( index < 0 ){
            single_timer timer( name );
            timers.push_back( timer );
            index = timers.size() - 1;
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(timers[index].tic));
    }
    double stop( std::string name ){
        
        int index = getIndex( name );

        if ( !timerExists( name, index ) ){ return 0.0;}

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(timers[index].toc));

        CD_OS_TIMESPEC_SUB(&(timers[index].toc), &(timers[index].tic));

        timers[index].elapsed = CD_OS_TIMESPEC_DOUBLE(&(timers[index].toc));
        timers[index].total += timers[index].elapsed;
        return timers[index].elapsed;
    }

    double reset( std::string name){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}

        double temp = timers[index].total;
        timers[index].total = 0;
        return temp;
    }

    double getTotal( std::string name ){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}
        return timers[index].total;
    }

    void wait( double time ){
        //wait for given amount of time
        struct timespec ticks_tic;
        struct timespec ticks_toc;

        /* start timing voxel grid computation */
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ticks_tic);

        while ( true ){
          /* stop timing voxel grid computation */
          clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ticks_toc);
          CD_OS_TIMESPEC_SUB(&ticks_toc, &ticks_tic);
          if ( time < CD_OS_TIMESPEC_DOUBLE(&ticks_toc) ){ break ; }
        }

    }
};

#endif
