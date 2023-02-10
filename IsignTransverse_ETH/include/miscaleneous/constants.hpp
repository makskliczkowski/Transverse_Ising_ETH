#pragma once

template <typename _type>
struct constants {
	static const _type pi;           //!< pi
	static const _type two_pi;          //!< 2*pi
	static const _type e;            //!< base of the natural logarithm
	static const _type euler;        //!< Euler's constant, aka Euler-Mascheroni constant
	static const _type gratio;       //!< golden ratio
	static const _type sqrt2;        //!< square root of 2
	static const _type sqrt2pi;      //!< square root of 2*pi
	static const _type log_sqrt2pi;  //!< log of square root of 2*pi

};


template<typename _type> const _type constants<_type>::pi			= _type(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679);
template<typename _type> const _type constants<_type>::two_pi	    = _type(6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359);
template<typename _type> const _type constants<_type>::e			= _type(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274);
template<typename _type> const _type constants<_type>::euler		= _type(0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495);
template<typename _type> const _type constants<_type>::gratio		= _type(1.6180339887498948482045868343656381177203091798057628621354486227052604628189024497072072041893911374);
template<typename _type> const _type constants<_type>::sqrt2		= _type(1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727);
template<typename _type> const _type constants<_type>::sqrt2pi		= _type(2.5066282746310005024157652848110452530069867406099383166299235763422936546078419749465958383780572661);
template<typename _type> const _type constants<_type>::log_sqrt2pi  = _type(0.9189385332046727417803297364056176398613974736377834128171515404827656959272603976947432986359541976);