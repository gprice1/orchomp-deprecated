
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

  public:
    bool print;
  
  private:
    class single_timer{
        public:
        std::string name;
        struct timespec tic, toc, wall_tic, wall_toc;
        double total, elapsed, wall_total, wall_elapsed;
        unsigned int count;
        bool isStopped;
        
        single_timer(std::string & name) : name( name ), total( 0.0 ),
                                           count( 0 ), isStopped( true ){}
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
    
    Timer() : print(false){}

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
        
        if ( print ){
            std::cout << "Starting " << name << std::endl;
        }
        if ( index < 0 ){
            single_timer timer( name );
            timers.push_back( timer );
            index = timers.size() - 1;
        }
        timers[ index ].isStopped = false;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(timers[index].tic));
        clock_gettime(CLOCK_MONOTONIC, &(timers[index].wall_tic));
    }
    double stop( std::string name ){
        
        int index = getIndex( name );

        if ( print ){
            std::cout << "Stopping " << name << std::endl;
        }

        if ( !timerExists( name, index ) ){ return 0.0;}

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(timers[index].toc));
        CD_OS_TIMESPEC_SUB(&(timers[index].toc), &(timers[index].tic));

        clock_gettime(CLOCK_MONOTONIC, &(timers[index].wall_toc));
        CD_OS_TIMESPEC_SUB(&(timers[index].wall_toc),
                           &(timers[index].wall_tic));

        timers[index].elapsed = 
                CD_OS_TIMESPEC_DOUBLE(&(timers[index].toc));
        timers[index].wall_elapsed =
                CD_OS_TIMESPEC_DOUBLE(&(timers[index].wall_toc));
        timers[index].total += timers[index].elapsed;
        timers[index].wall_total += timers[index].wall_elapsed;

        timers[index].count ++;
        timers[index].isStopped = true;

        return timers[index].elapsed;
    }

    double reset( std::string name){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}

        double temp = timers[index].total;
        timers[index].total = 0;
        timers[index].count = 0;
        timers[index].elapsed = 0;
        return temp;
    }

    double getTotal( std::string name ){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}
        return timers[index].total;
    }
    double getWallTotal( std::string name ){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}
        return timers[index].wall_total;
    }
    double getWallElapsed( std::string name ){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}
        return timers[index].wall_elapsed;
    }
    double getElapsed( std::string name ){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}
        return timers[index].elapsed;
    }

    //returns true if the timer is already in a start state, elsewise,
    //  it is started.
    bool tryStart( std::string name ){
        int index = getIndex( name );

        if ( index < 0 ){
            single_timer timer( name );
            timers.push_back( timer );
            index = timers.size() - 1;
        }

        if ( timers[ index ].isStopped ){
            timers[ index ].isStopped = false;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(timers[index].tic));
            return false;
        }

        return false;
    }

    //returns the elapsed time if if the timer is already in a
    //  stopped state, elsewise, it is started.
    bool tryStop( std::string name ){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}
        if ( timers[ index ].isStopped ){
            return true;
        }

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(timers[index].toc));

        CD_OS_TIMESPEC_SUB(&(timers[index].toc), &(timers[index].tic));

        timers[index].elapsed = CD_OS_TIMESPEC_DOUBLE(&(timers[index].toc));
        timers[index].total += timers[index].elapsed;
        timers[index].count ++;
        timers[index].isStopped = true;

        return false;
    }

    unsigned int getCount( std::string name ){
        int index = getIndex( name );
        if ( !timerExists( name, index ) ){ return 0.0;}
        return timers[index].count;
    }

    void wait( double time ){
        //wait for given amount of time
        struct timespec ticks_tic;
        struct timespec ticks_toc;

        /* start timing voxel grid computation */
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ticks_tic);

        while ( true ){
          /* stop timing voxel grid computation */
          clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ticks_toc);
          CD_OS_TIMESPEC_SUB(&ticks_toc, &ticks_tic);
          if ( time < CD_OS_TIMESPEC_DOUBLE(&ticks_toc) ){ break ; }
        }

    }
};

#endif
