/* TemplateTools.h
 * 
 * Author           : Alexander J. Yee
 * Date Created     : 11/20/2014
 * Last Modified    : 11/20/2014
 * 
 */

#pragma once
#ifndef _ymp_TemplateTools_H
#define _ymp_TemplateTools_H
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
#include <limits.h>
#include <new>
#include <utility>
#include <type_traits>
#include "CompilerSettings.h"
namespace ymp{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Type ID without RTTI  (Be careful, this is technically undefined behavior.)
template <typename typeA, typename typeB> YM_FORCE_INLINE
bool is_same_type(const typeA& A, const typeB& B){
//    return typeid(a) == typeid(b);
    return *(void**)&A == *(void**)&B;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename T, class... Args>
YM_FORCE_INLINE T* PlaceTrivialObject(void* mem, Args&&... args){
    static_assert(
        std::is_trivially_destructible<T>::value,
        "Object must be trivially destructible."
    );
    return ::new (mem) T(std::forward<Args>(args)...);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Helper for the Flat-Index Micro-Optimization
#if 0
template <typename type> YM_FORCE_INLINE
type& ByteOffset(type* ptr, size_t bytes, size_t objects){
    return *(type*)((char*)ptr + bytes + objects * sizeof(type));
}
#endif
#if 0
#define ByteOffset(ptr, bytes, objects) *(typeof(type)*)((char*)ptr + bytes + objects * sizeof(*type));
#endif
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Is valid word size.
#define YMP_ASSERT_VALID_WORD(wtype) static_assert(WordTraits<wtype>::is_valid, "wtype must be an unsigned 32 or 64-bit integer");
template <typename wtype>
struct WordTraits{
    static const size_t BITS = CHAR_BIT * sizeof(wtype);
    static const size_t MAG  = BITS == 32 ? 5 : 6;
    static const bool is_valid =
        std::is_integral<wtype>::value &&
        std::is_unsigned<wtype>::value &&
        (BITS == 32 || BITS == 64);
    YMP_ASSERT_VALID_WORD(wtype);
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
}
#endif