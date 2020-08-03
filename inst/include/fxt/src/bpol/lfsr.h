#if !defined  HAVE_LFSR_H__
#define       HAVE_LFSR_H__
// This file is part of the FXT library.
// Copyright (C) 2010, 2012, 2014 Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the main directory.

#include "fxttypes.h"
#include "bits/bitsperlong.h"

#include "bpol/poly-tab.h"


class lfsr
// (binary) Linear Feedback Shift Register
// Produces a shift register sequence (SRS)
// generated by  a_k=x^k (mod c) where
// c is a primitive polynomial of degree n.
// The period of SRS is 2^n - 1
// (non-primitive c lead to smaller periods)
{
public:
    ulong a_;  // internal state (polynomial modulo c)
    ulong w_;  // word of the shift_register_sequence (SRS)
    ulong c_;  // (mod 2) poly  e.g. x^4+x+1 == 0x13 == 1..11
    ulong h_;  // highest bit in SRS word  e.g. (above) == 16 = 1...
    ulong mask_;  // mask  e.g. (above) == 15 == 1111
    ulong n_;  // degree of polynomial  e.g. (above) == 4

public:
    explicit lfsr(ulong n, ulong c=0)
    // n: degree of polynomial c
    // c: polynomial (defaults to minimum weight primitive polynomial)
    {
        // must have  1 <= n <= BITS_PER_LONG
        if ( n>BITS_PER_LONG )  n = 1;
        if ( n==0 )  n = 1;

        // no poly supplied: use minimum weight poly
        c_ = c;
        if ( 0==c_ )  c_ = minweight_primpoly[n];

        h_ = (1UL<<(n-1));
        mask_ = (h_ | (h_-1) );
        n_ = n;

//        set_w(1);  // ==> a == c (except highest order)
        set_a(1);
    }

    ~lfsr()  {;}

    ulong next()
    {
        ulong s = a_ & h_;
        a_ <<= 1;
        w_ <<= 1;
        if ( 0!=s )
        {
            a_ ^= c_;
            w_ |= 1;
        }
        w_ &= mask_;
        return w_;
    }

    ulong next_w()
    // produce the next word w
    // 1 <= w <= mask
    {
        for (ulong k=0; k<n_; ++k)  next();
        return w_;
    }

private:
    void prev_a()
    {
        ulong s = a_ & 1;
        a_ >>= 1;
        if ( s )
        {
            a_ ^= (c_>>1);
            a_ |= h_;  // so it works for  n_ == BITS_PER_LONG
        }
    }

public:
    ulong prev()
    {
        prev_a();
        set_a(a_);
        return w_;
    }

    ulong prev_w()
    // produce the previous word w
    // 1 <= w <= mask
    {
        for (ulong k=0; k<n_; ++k)  prev_a();
        set_a(a_);
        return w_;
    }

    ulong get_a()  const  { return a_; }

    void set_a(ulong a)
    {
        a_ = a;
        w_ = 0;
        ulong b = 1;
        for (ulong j=0; j<n_; ++j)
        {
            if ( a & 1 )
            {
                w_ |= b;
                a ^= c_;
            }
            b <<= 1;
            a >>= 1;
        }
    }

    ulong get_w()  const  { return w_; }

    void set_w(ulong w)
    // must have w!=0
    {
        w_ = w;
        a_ = 0;
        ulong c = c_;
        while ( w )
        {
            if ( w & 1 )  a_ ^= c;
            c <<= 1;
            w >>= 1;
        }
        a_ &= mask_;
    }

    ulong max_period()  const
    // period of srs if polynomial is primitive
    {
        return  (~0UL) >> (BITS_PER_LONG - n_);
    }
};
// -------------------------


#endif  // !defined HAVE_LFSR_H__
