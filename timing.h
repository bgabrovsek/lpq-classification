//  timer.h
//  Created by Boštjan on 4/22/13.
//
//  record time consumption of  a process
//
//  usage:
//  timer->startTimer(1);
//  ...
//  timer->pauseTimer(1);
//  cout << "Time: " << timer->seconds(1);
//  timer->resetTimer(1);
//  ...

#ifndef Lpq_timing_h
#define Lpq_timing_h

#define MAX_TIMERS 10

class Ctimer {
    
public:
    
    double cpu[MAX_TIMERS];
    time_t cpu_start[MAX_TIMERS];
    
    // constructor, initialization
    Ctimer() {for (int i=0;i<MAX_TIMERS;i++) cpu[i] = 0;};
    
    void resetTimer(int i) { cpu[i] = 0; };
    
    void startTimer(int i) { cpu_start[i] = clock();};
    
    void resetAndStartTimer(int i) { cpu[i] = 0; cpu_start[i] = clock(); }
    
    void pauseTimer(int i) { cpu[i] += 1.0*(clock() - cpu_start[i])/ (CLOCKS_PER_SEC);}
    
    // return time for timer i
    double seconds(int i) { return cpu[i]; }
    
};

Ctimer *timer = new Ctimer();


#endif
