#include <fastenvelope/get_mem.h>

extern "C" size_t getCurrentRSS();
extern "C" size_t getPeakRSS();

size_t fastEnvelope::get_mem(){
    return getCurrentRSS()/(1024*1024);
}

size_t fastEnvelope::get_peak_mem(){
    return getPeakRSS()/(1024*1024);
}