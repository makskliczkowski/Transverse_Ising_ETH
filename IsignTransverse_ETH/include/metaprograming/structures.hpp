#pragma once
#ifndef __COMMONS
    #include "commons.hpp"
#endif

//<! FORWARD DECLARATIONS
template<size_t idx, typename _ty>
struct GetHelper;
template<typename ... _types>
struct variadic_struct {};

/// @brief Structure to contain symmetry sectors
/// @tparam ..._types remaining template types for sectors
/// @tparam _ty first template type for symmetry sector
template <typename _ty, typename ... _types>
struct variadic_struct<_ty, _types ...>
{
    variadic_struct() = default;
    variadic_struct(const _ty& first, const _types& ... rest)
        : first(first)
        , rest(rest...)
    {}
    
    _ty first;                                
    variadic_struct<_types ... > rest;
    
    template<size_t idx>
    auto get()
        { return GetHelper<idx, variadic_struct<_ty, _types...>>::get(*this); }
};

/// @brief Helper structure to get first element in variadic struct
/// @tparam _ty first type
/// @tparam ..._types remaining template types
template<typename _ty, typename ... _types>
struct GetHelper<0, variadic_struct<_ty, _types ... >>
{ 
    static _ty get(variadic_struct<_ty, _types...>& data)
        { return data.first; }
};

/// @brief General structure to get element in variadic struct by index idx
/// @tparam _ty first type
/// @tparam ..._types remaining template types
/// @tparam idx index of element to get from variadic struct
template<size_t idx, typename _ty, typename ... _types>
struct GetHelper<idx, variadic_struct<_ty, _types ... >>
{ 
    static auto get(variadic_struct<_ty, _types...>& data)
        { return GetHelper<idx - 1, variadic_struct<_types ...>>::get(data.rest); }
};
