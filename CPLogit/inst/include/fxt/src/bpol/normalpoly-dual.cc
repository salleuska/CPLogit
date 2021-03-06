// This file is part of the FXT library.
// Copyright (C) 2010, 2012, 2018 Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the main directory.

#include "bpol/gf2n-trace.h"
#include "bpol/bitpolmod-arith.h"
#include "bpol/bitpolmod-minpoly.h"

#include "fxttypes.h"

//#include "jjassert.h"


ulong
gf2n_xx2k_trace(ulong c, ulong deg)
// Return vector T of traces T[k]=trace(ek),
// where ek = x*x^(2^k), k=0..deg-1, and
//   x is a root of the irreducible polynomial C.
// Must have:  deg == degree(C)
//
// The polynomial C is normal if and only if
//   gcd(x^deg-1, T)==1  where T is taken as a polynomial.
// The dual basis is generated by sum(k=0..deg-1, D[k]*x^(2^k)) where
//   D = T^-1 mod x^deg-1.  The dual basis is also normal.
// If the only nonzero component of T is T[0], then
//   the basis is self-dual.
{
    if ( c==3 )  return 1UL;  // x+1 is self-dual

    const ulong tv = gf2n_trace_vector_x(c, deg);  // traces of x^k
    ulong rt = 2UL;  // root of C
    const ulong h = 1UL << (deg-1);  // aux

    ulong T = 0;
    for (ulong k=0; k<deg; ++k)
    {
        ulong ek = bitpolmod_times_x(rt, c, h);  // == x*x^(2^k)
        ulong tk = gf2n_fast_trace(ek, tv);      // == sum(ek[i]*tk[i])
        T |= (tk<<k);
        rt = bitpolmod_square(rt, c, h);
    }

    return  T;
}
// -------------------------


ulong
gf2n_dual_normal(ulong c, ulong deg, ulong ntc/*=0*/, ulong *ntd/*=0*/)
// Return the minimal polynomial CS for the dual (normal) basis
//  with the irreducible normal polynomial C.
// Return zero if C is not normal.
// Must have:  deg == degree(C).
// If ntc is supplied it must be equal to gf2n_xx2k_trace(c, deg).
// If ntd is nonzero, ntc^-1 (mod x^deg-1) is written to it.
{
    if ( 0==ntc )  ntc = gf2n_xx2k_trace(c, deg);
    const ulong d = bitpolmod_inverse(ntc, 1 | (1UL<<deg) );  // ntc=d^-1 (mod x^deg-1)
    if ( 0==d )  return 0;  // C not normal
    if ( nullptr != ntd )  *ntd = d;

    const ulong h = 1UL << (deg-1);  // aux
    ulong alpha = 2UL;  // 'x', a root of C
    ulong beta = 0;     // root of the dual polynomial
    for (ulong m=d; m!=0; m>>=1)
    {
        if ( m & 1 )  beta ^= alpha;
        alpha = bitpolmod_square(alpha, c, h);
    }

    ulong cs;  // minimal polynomial of beta
    bitpolmod_minpoly(beta, c, deg, cs);
//    csdeg = bitpol_minpoly(beta, c, deg, cs);
//    jjassert( deg == csdeg );

    return cs;
}
// -------------------------
