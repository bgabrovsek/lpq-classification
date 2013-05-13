//  knot_preprocessor_directives.h
//  Created by Boštjan on 4/22/13.
//
//  sets the global types to compule, depending on the maximal number of crossings
//  using a lower number of crossings gains on performance

#ifndef knot_preprocessor_directives_h
#define knot_preprocessor_directives_h

#if MAX < 32 // n = 0...31

#define u_number u8 // numbers, number of regions, arc number, region number.
#define s_letter s8 // crossing letter for GW
#define u_sign u32 // sign bits
#define u_region u64 // regions
#define s_small s8 // for orientation, direction.
#define fbf_sign(n) first_bits_full_32(n) // signs
#define fbf_reg(n) first_bits_full_64(n) // regions
#define infinity_number 0xFF // infinity for u_number

#elif MAX < 64 // n = 32...63

#define u_number u8 // 
#define s_letter s8
#define u_sign u64
#define u_region u128
#define s_small s8
#define fbf_sign(n) first_bits_full_64(n)
#define fbf_reg(n) first_bits_full_128(n)
#define infinity_number 0xFF // infinity for u_number


#elif MAX < 128 // n = 64...127

#define u_number u16 // ?
#define s_letter s16 // ?
#define u_sign u128
#define u_region u256
#define s_small s8
#define fbf_sign(n) first_bits_full_128(n)
#define fbf_reg(n) first_bits_full_256(n)
#define infinity_number 0xFFFF // infinity for u_number


#elif MAX < 256 // n = 128...256

#define u_number u16
#define s_letter s16
#define u_sign u256
#define u_region u512
#define s_small s8
#define fbf_sign(n) first_bits_full_256(n)
#define fbf_reg(n) first_bits_full_512(n)
#define infinity_number 0xFFFF // infinity for u_number


#endif




#endif
