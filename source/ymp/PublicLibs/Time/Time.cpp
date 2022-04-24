/* Time.cpp
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 09/17/2014
 * Last Modified    : 09/17/2014
 * 
 */

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Dependencies
#ifdef _WIN32
#include "Time_Windows.ipp"
#else
#include "Time_Posix.ipp"
#endif
#include "Time.h"
namespace ymp{
namespace Time{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
YM_NO_INLINE void print_secs_hrs(double seconds, char color){
    Console::SetColor(color);
    Console::print_fixed(seconds);
    Console::print(" seconds  ( ");
    Console::print_fixed(seconds / 3600);
    Console::print(" hours )");
    if (color != ' ')
        Console::SetColor('w');
}
YM_NO_INLINE void println_secs_hrs(double seconds, char color){
    Console::SetColor(color);
    Console::print_fixed(seconds);
    Console::print(" seconds  ( ");
    Console::print_fixed(seconds / 3600);
    Console::println(" hours )");
    if (color != ' ')
        Console::SetColor('w');
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
}