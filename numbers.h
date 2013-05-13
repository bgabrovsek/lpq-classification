//  numbers.h
//  Created by Boštjan on 4/22/13.
//
//  128 bit, 256 bit and 512 bit integers
//
//  TODO: generalize the class for all integer operations

#include "common.h"

#ifndef numbers_h
#define numbers_h

#define u8 uint8_t
#define u16 uint16_t
#define u32 uint32_t
#define u64 uint64_t

#define s8 int8_t
#define s16 int16_t
#define s32 int32_t
#define s64 int64_t

class u128;
class u256;

using namespace std;

class u128 {
public:
    u64 hi, lo;
    
    // constructors
    u128() { hi = 0; lo = 0;}
    u128(u64 l) { hi = 0; lo = l;}
    u128(u64 h, u64 l) { hi = h; lo = l;}
    u128(const u128 &i) { hi = i.hi; lo = i.lo;}
    
    // assignment operator
    u128 &operator=(const u128 &rhs) { if (this != &rhs) {lo = rhs.lo; hi = rhs.hi;} return *this; }
    // binary operators
    u128 operator&(const u128 &rhs) { return u128(hi & rhs.hi, lo & rhs.lo); }
    u128 operator&(const int &rhs) { return u128(0, lo & (u64)rhs); }
    u128 operator|(const u128 &rhs) { return u128(hi | rhs.hi, lo | rhs.lo); }
    u128 operator^(const u128 &rhs) { return u128(hi ^ rhs.hi, lo ^ rhs.lo); }
    u128 operator<<(const int &rhs) { if (rhs) return (rhs < 64 ? u128( (hi<<rhs) | (lo>>(64-rhs)), lo<<rhs) : u128(lo<<(rhs-64),0)); else return u128(*this); }
    u128 operator>>(const int &rhs) { if (rhs) return (rhs < 64 ? u128( hi>>rhs, (lo>>rhs) | (hi<<(64-rhs))) : u128(0,hi>>(rhs-64))); else return u128(*this); }
    
    // bit operator assignment
    u128 &operator|=(const u128 &rhs) { hi |= rhs.hi; lo |= rhs.lo; return *this; }
    u128 &operator&=(const u128 &rhs) { hi &= rhs.hi; lo &= rhs.lo; return *this; }
    u128 &operator^=(const u128 &rhs) { hi ^= rhs.hi; lo ^= rhs.lo; return *this; }
    u128 &operator<<=(const int &rhs) { return (*this = *this << rhs);  }
    u128 &operator>>=(const int &rhs) { return (*this = *this >> rhs);  }
    
    // comparison operators
    bool operator==(const u128 &rhs) { return ((lo == rhs.lo) && (hi == rhs.hi)); };
    bool operator!=(const u128 &rhs) { return ((lo != rhs.lo) || (hi != rhs.hi)); };
    bool operator>=(const u128 &rhs) { return ((hi != rhs.hi) ? hi >= rhs.hi : lo >= rhs.lo); };
    bool operator<=(const u128 &rhs) { return ((hi != rhs.hi) ? hi <= rhs.hi : lo <= rhs.lo); };
    bool operator>(const u128 &rhs) { return ((hi != rhs.hi) ? hi > rhs.hi : lo > rhs.lo); };
    bool operator<(const u128 &rhs) { return ((hi != rhs.hi) ? hi < rhs.hi : lo < rhs.lo); };
    // operation operator
    //u128 &operator+=(const u128 &b) { u64 old_lo = lo; lo += b.lo; hi += b.hi; if (lo < old_lo) hi+=1; return *this; }
    // unary operators
    bool operator!() { return (*this == u128(0)); }
    u128 operator~() { return u128(~hi, ~lo); }
    u128 operator++() { ++lo; if (lo == 0) ++hi; return *this;}
    // conversion operators
    operator bool() { return (*this != u128(0));}
    operator int() { if (hi) cout << "ERROR reg->int. " << endl; return (int)lo;}
        
};
        
// output
ostream &operator<<(ostream &out, u128 n) {
    out << "[ ";
    for (int i=0;i<128;i++) if ( i < 64 ? (n.lo & ((u64)1<<i)) : (n.hi & ((u64)1 << (i-64))) ) out << i << " ";
    out << "] ";
    return out;
}


class u256 {
public:
    u128 hi, lo;
    
    // constructors
    u256() { hi = 0; lo = 0;}
    u256(u64 l) { hi = 0; lo = l;}
    u256(u128 h, u128 l) { hi = h; lo = l;}
    u256(const u256 &i) { hi = i.hi; lo = i.lo;}
    
    // assignment operator
    u256 &operator=(const u256 &rhs) { if (this != &rhs) {lo = rhs.lo; hi = rhs.hi;} return *this; }
    
    // binary operators
    u256 operator&(const u256 &rhs) { return u256(hi & rhs.hi, lo & rhs.lo); }
    u256 operator&(const int &rhs) { return u256(0, lo & (u128)(u64)rhs); }
    u256 operator|(const u256 &rhs) { return u256(hi | rhs.hi, lo | rhs.lo); }
    u256 operator^(const u256 &rhs) { return u256(hi ^ rhs.hi, lo ^ rhs.lo); }
    u256 operator<<(const int &rhs) { if (rhs) return (rhs<128 ? u256( (hi<<rhs) | (lo>>(128-rhs)), lo<<rhs):u256(lo<<(rhs-128),0)); else return u256(*this); }
    u256 operator>>(const int &rhs) { if (rhs) return (rhs<128 ? u256( hi>>rhs, (lo>>rhs) | (hi<<(128-rhs))):u256(0,hi>>(rhs-128))); else return u256(*this); }
    
    // bit operator assignment
    u256 &operator|=(const u256 &rhs) { hi |= rhs.hi; lo |= rhs.lo; return *this; }
    u256 &operator&=(const u256 &rhs) { hi &= rhs.hi; lo &= rhs.lo; return *this; }
    u256 &operator^=(const u256 &rhs) { hi ^= rhs.hi; lo ^= rhs.lo; return *this; }
    u256 &operator<<=(const int &rhs) { return (*this = *this << rhs);  }
    u256 &operator>>=(const int &rhs) { return (*this = *this >> rhs);  }
    
    // comparison operators
    bool operator==(const u256 &rhs) { return ((lo == rhs.lo) && (hi == rhs.hi)); };
    bool operator!=(const u256 &rhs) { return ((lo != rhs.lo) || (hi != rhs.hi)); };
    bool operator>=(const u256 &rhs) { return ((hi != rhs.hi) ? hi >= rhs.hi : lo >= rhs.lo); };
    bool operator<=(const u256 &rhs) { return ((hi != rhs.hi) ? hi <= rhs.hi : lo <= rhs.lo); };
    bool operator>(const u256 &rhs) { return ((hi != rhs.hi) ? hi > rhs.hi : lo > rhs.lo); };
    bool operator<(const u256 &rhs) { return ((hi != rhs.hi) ? hi < rhs.hi : lo < rhs.lo); };
    // operation operator
    //u256 &operator+=(const u256 &b) { u128 old_lo = lo; lo += b.lo; hi += b.hi; if (lo < old_lo) hi+=1; return *this; }
    
    // unary operators
    bool operator!() { return (*this == u256(0)); }
    u256 operator~() { return u256(~hi, ~lo); }
    u256 operator++() { ++lo; if (lo == (u128)0) ++hi; return *this;}
    // conversion operators
    operator bool() { return (*this != u256(0));}
    operator int() { if (hi) cout << "ERROR reg->int. " << endl; return (int)lo;}
        
};

class u512 {
public:
    u256 hi, lo;
            
    // constructors
    u512() { hi = 0; lo = 0;}
    u512(u64 l) { hi = 0; lo = l;}
    u512(u256 h, u256 l) { hi = h; lo = l;}
    u512(const u512 &i) { hi = i.hi; lo = i.lo;}
            
    // assignment operator
    u512 &operator=(const u512 &rhs) { if (this != &rhs) {lo = rhs.lo; hi = rhs.hi;} return *this; }
            
    // binary operators
    u512 operator&(const u512 &rhs) { return u512(hi & rhs.hi, lo & rhs.lo); }
    u512 operator&(const int &rhs) { return u512(0, lo & (u256)(u64)rhs); }
    u512 operator|(const u512 &rhs) { return u512(hi | rhs.hi, lo | rhs.lo); }
    u512 operator^(const u512 &rhs) { return u512(hi ^ rhs.hi, lo ^ rhs.lo); }
    u512 operator<<(const int &rhs) { if (rhs) return (rhs<256 ? u512((hi<<rhs) | (lo>>(256-rhs)), lo<<rhs) : u512(lo<<(rhs-256),0)); else return u512(*this); }
    u512 operator>>(const int &rhs) { if (rhs) return (rhs<256 ? u512(hi>>rhs, (lo>>rhs) | (hi<<(256-rhs))) : u512(0,hi>>(rhs-256))); else return u512(*this); }
            
    // bit operator assignment
    u512 &operator|=(const u512 &rhs) { hi |= rhs.hi; lo |= rhs.lo; return *this; }
    u512 &operator&=(const u512 &rhs) { hi &= rhs.hi; lo &= rhs.lo; return *this; }
    u512 &operator^=(const u512 &rhs) { hi ^= rhs.hi; lo ^= rhs.lo; return *this; }
    u512 &operator<<=(const int &rhs) { return (*this = *this << rhs);  }
    u512 &operator>>=(const int &rhs) { return (*this = *this >> rhs);  }
            
    // comparison operators
    bool operator==(const u512 &rhs) { return ((lo == rhs.lo) && (hi == rhs.hi)); };
    bool operator!=(const u512 &rhs) { return ((lo != rhs.lo) || (hi != rhs.hi)); };
    bool operator>=(const u512 &rhs) { return ((hi != rhs.hi) ? hi >= rhs.hi : lo >= rhs.lo); };
    bool operator<=(const u512 &rhs) { return ((hi != rhs.hi) ? hi <= rhs.hi : lo <= rhs.lo); };
    bool operator>(const u512 &rhs) { return ((hi != rhs.hi) ? hi > rhs.hi : lo > rhs.lo); };
    bool operator<(const u512 &rhs) { return ((hi != rhs.hi) ? hi < rhs.hi : lo < rhs.lo); };
    // operation operator
    //u512 &operator+=(const u512 &b) { u256 old_lo = lo; lo += b.lo; hi += b.hi; if (lo < old_lo) hi+=1; return *this; }
            
    // unary operators
    bool operator!() { return (*this == u512(0)); }
    u512 operator~() { return u512(~hi, ~lo); }
    u512 operator++() { ++lo; if (lo == (u256)0) ++hi; return *this;}
    // conversion operators
    operator bool() { return (*this != u512(0));}
    operator int() { if (hi) cout << "ERROR reg->int. " << endl; return (int)lo;}
            
};
        
ostream &operator<<(ostream &out, u256 n) {
    //cout << n.hi << " " << n.lo << " ";
    out << "[ ";
    for (int i=0;i<256;i++) if ( i < 128 ? (n.lo & ((u128)1<<i)) : (n.hi & ((u128)1 << (i-128))) ) out << i << " ";
    out << "] ";
    return out;
}

ostream &operator<<(ostream &out, u512 n) {
    //cout << n.hi << " " << n.lo << " ";
    out << "[ ";
    for (int i=0;i<512;i++) if ( i < 256 ? (n.lo & ((u256)1<<i)) : (n.hi & ((u256)1 << (i-256))) ) out << i << " ";
    out << "] ";
    return out;
}
        
u64 first_bits_full_64_[65] = { 0x0,0x1,0x3,0x7,0xf,0x1f,0x3f,0x7f,
0xff, 0x1ff, 0x3ff, 0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff,
0xffff, 0x1ffff, 0x3ffff, 0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff,
0xffffff, 0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff, 0x3fffffff, 0x7fffffff,
0xffffffff, 0x1ffffffff, 0x3ffffffff, 0x7ffffffff, 0xfffffffff, 0x1fffffffff, 0x3fffffffff, 0x7fffffffff,
0xffffffffff, 0x1ffffffffff, 0x3ffffffffff, 0x7ffffffffff, 0xfffffffffff, 0x1fffffffffff, 0x3fffffffffff, 0x7fffffffffff,
0xffffffffffff, 0x1ffffffffffff, 0x3ffffffffffff, 0x7ffffffffffff, 0xfffffffffffff, 0x1fffffffffffff, 0x3fffffffffffff, 0x7fffffffffffff,
0xffffffffffffff, 0x1ffffffffffffff, 0x3ffffffffffffff, 0x7ffffffffffffff, 0xfffffffffffffff, 0x1fffffffffffffff, 0x3fffffffffffffff, 0x7fffffffffffffff,
        0xffffffffffffffff};

u32 first_bits_full_32(int n) { if (n > 32) throw 20; return (u32)first_bits_full_64_[n];}
u64 first_bits_full_64(int n) { if (n > 64) throw 20; return first_bits_full_64_[n];}
u128 first_bits_full_128(int n) { if (n > 128) throw 20; return (n>64 ? u128(first_bits_full_64(n-64),first_bits_full_64(64)) : u128(0,first_bits_full_64(n)) ); }
u256 first_bits_full_256(int n) { if (n > 256) throw 20; return (n>128 ? u256(first_bits_full_128(n-128),first_bits_full_128(128)) : u256(0,first_bits_full_128(n)) ); }
u512 first_bits_full_512(int n) { return (n>256 ? u512(first_bits_full_256(n-256),first_bits_full_256(256)) : u512(0,first_bits_full_256(n)) ); }


#endif
