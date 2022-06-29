#pragma once
#ifdef HAS_CXX20
	#include <concepts> // has included type_traits, but is available for C++20 or newer
#else
	#include <type_traits>
#endif

// hold non-STL standard traits
namespace traits {

	//<! check if variadic templates have common types (or are inmplicitly convertible) with _baseTy
	template <class _baseTy, class... _otherTy>
	struct is_convertible : std::conjunction<std::is_convertible<_otherTy, _baseTy>...> {};
	
	template <class _baseTy, class... _otherTy>
	inline constexpr bool is_convertible_v = is_convertible<_baseTy, _otherTy...>::value;

	//<! check if variadic templates have common types (or are inmplicitly convertible) with _baseTy
	template <class _baseTy, class... _otherTy>
	struct is_same : std::conjunction<std::is_same<_otherTy, _baseTy>...> {};

	template <class _baseTy, class... _otherTy>
	inline constexpr bool is_same_v = is_same<_baseTy, _otherTy...>::value;

	//<! variadic check if classes are copy constructible
	template <class... Args >
	struct is_copy_constructible : std::conjunction<std::is_copy_constructible<Args>...> {};
	
	template <class... Args >
	inline constexpr bool is_copy_constructible_v = is_copy_constructible<Args...>::value;

	//<! check if typename is among other fixed ones struct (STL only has constexpr)
	template <class _ty, class... fixed>
	struct is_any_of : std::disjunction<std::is_same<_ty, fixed>...> {};

	template <class _ty, class... fixed>
	inline constexpr bool is_any_of_v = is_any_of<_ty, fixed...>::value;

	//<! same as above but checks if whole parameter pack has types among fixed one's
	template <class... input>
	struct amongst {
		template <class... fixed>
		struct is_included : std::disjunction<is_any_of<input, fixed...>...> {};

		template <class... fixed>
		inline static constexpr bool is_included_v = is_included<fixed...>::value;
	};

};
//---------------------NOTES: 
// - std::conjuction:
//		is variadic struct taking in traits
//		and giving compiler error if AT LEAST ONE among args does not satisfy the trait inside.
//		It 'inherits' from bool_constant (here true_type) so we can inherit from indirectly via std::conjunction
// 
// - std::disjunction:
//		similar as std::conjunction but prints compiler error if NONE of the args satisfy the trait.
//		Used to pick one among the input traits. Also inherits from bool_constant (but here false_type)
// --------------------------

