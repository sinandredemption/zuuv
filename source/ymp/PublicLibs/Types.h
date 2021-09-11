/* Types.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 07/22/2014
 * Last Modified    : 07/22/2014
 * 
 */

#pragma once
#ifndef _ymp_Types_H
#define _ymp_Types_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include <stdint.h>
namespace ymp{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Fixed Sizes
typedef uint16_t    u16_t;
typedef uint32_t    u32_t;
typedef int32_t     s32_t;
typedef uint64_t    u64_t;
typedef int64_t     s64_t;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Logarithmic Size
typedef uint_fast32_t   ukL_t;
typedef int_fast32_t    skL_t;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Pointer Sizes
////////////////////////////////////////////////////////////////////////////////
#if _M_X64 || __x86_64
#define YMP_PTR_MAG     6
const ukL_t PTR_BITS = 64;
typedef u64_t       upL_t;
typedef s64_t       spL_t;
#else
#define YMP_PTR_MAG     5
const ukL_t PTR_BITS = 32;
typedef u32_t       upL_t;
typedef s32_t       spL_t;
#endif
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  File Sizes
typedef u64_t       ufL_t;
typedef s64_t       sfL_t;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Integer Index Size
#define YMP_INDEX_MAG   6
const ukL_t INDEX_BITS = 64;
typedef u64_t uiL_t;
typedef s64_t siL_t;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Float Index Size
typedef double fiL_t;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif
