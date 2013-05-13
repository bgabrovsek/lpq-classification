//  common.h
//  Created by Boštjan on 4/22/13.
//
//  Various constants, short functions and macros: bit-wise operations, memory functions,
//  macros used for lexicographical comparing, parsing stirngs, HOMFLY varaible names,...

#ifndef Lpq_common_h
#define Lpq_common_h

#include <stdint.h> // for uint64_t
#include <string>
#include <list>
#include <sstream>

#include "numbers.h"

using namespace std;

// GLOBAL CONSTANTS

#define YES true
#define NO false

// maximal knots we can consider
#define MAX_KNOTS 10000

// binary signs
#define PLUS 0
#define MINUS 1

// R-moves flags
#define REMOVE_R_I 0x1
#define CREATE_R_I 0x2
#define REMOVE_R_II 0x4
#define CREATE_R_II 0x8
#define MODIFY_R_III 0x10
#define MODIFY_FLYPE 0x20
#define ROTATE_MERIDIAN 0x40

// simple reduction flags
#define REMOVE_3_KINKS 1
#define REMOVE_AFFINE 2
#define REMOVE_KINK_OUTSIDE_DOT 4
#define REMOVE_ADJACENT_SELECTED 8
#define REMOVE_BST_REMOVABLE 16
#define REMOVE_OBVIOUS_COMPOSITES 32
#define REMOVE_SHRINKABLE 64
#define REMOVE_CONNECTED_SUM 128

// orientation
#define ORIENTED true
#define UNORIENTED false

// directions for looping through arcs
#define RIGHT 0
#define LEFT 1

// knot properties (obsolete?)

// chirality (oboslete?)
#define CHIRAL 1 // different than mirror
#define AMPHICHIRAL 2 // same as mirror

// obsolete
#define FLIPPABLE 1
#define UNFLIPPABLE 2

// obsolete
#define PRIME 1
#define UNPRIME 2

// GLOBAL VARIABLES

// MEMORY MANEGMENT (UNIX SYSTEMS)


#define installed_ram_kb 3300000

int parseLine(char* line){
    int i = (int)strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}

int get_used_memory_kb_linux(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];
    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

int get_used_memory_kb() {
#ifdef __linux__
    return get_used_memory_kb_linux();
#else
    return 0;
#endif
}

// VARIOUS MACROS

#define sign2char(i) ((i==PLUS) ? '+' : '-')

#define SIGN(i) ((i)<0?-1:+1)
#define SIGN2B(i) ((i)<0?1:0)//((i)>>1)&1)
#define B2SIGN(i) ((i)?-1:1)
#define SIGNB(i) ((i)<0?MINUS:PLUS)

#define MOD(a, m) (((a)%(m)+(m))%(m)) // real modulo
#define MOD_EQ(a,b,m) (MOD((a),(m)) == MOD((b),(m)))

#define ABS(a) ((a)>=0?(a):-(a))
#define MINIMUM(a,b) ((a)>=(b) ? (b) : (a))
#define MAXIMUM(a,b) ((a)>=(b) ? (a) : (b))
#define MAXIMUM3(a,b,c) ( MAXIMUM(MAXIMUM(a,b),c) )
#define MINIMUM3(a,b,c) ( MINIMUM(MINIMUM(a,b),c) )

#define ABS_SUCC(a) (SIGN(a) * (ABS(a)+1))
#define ABS_PRED(a) (SIGN(a) * (ABS(a)-1))


// compare two gauss letters, A < B (return -1)
#define COMPARE_LETTER(A,B) ((A) == (B) ? 0 : ( (A) == -(B) ? ( A < 0 ? 1 : -1) : ((ABS(A) > ABS(B) ? 1 : -1))))

// compare signs of crossings A<B (return -1)
#define COMPARE_SIGN(A,B) ((int)(A) - (int)(B))

#define MOD_DISTANCE(a,b,m) (MINIMUM(MAXIMUM(a,b)-MINIMUM(a,b),MINIMUM(a,b)+m-MAXIMUM(a,b))) 

int use_extended_bitset = false;

u_number mod_inverse(u_number a, u_number m) {
    a = MOD(a,m);
    for(u_number x = 1; x < m; x++)
        if ( MOD(a*x,m) == 1) return x;
    return 0;
}

int COMPARE_REGION(u_region r_A, u_region r_B) { // A<B (-1); A>B (1), A==B (0).
    while (r_A | r_B) {
        
        if (!r_A) return -1;
        if (!r_B) return 1;
        if ((r_A & 1) ^ (r_B & 1)) return (int)((r_B & 1) << 1) -1;
        
        r_A >>= 1; r_B >>= 1;
    }
    return 0;
}

// STRING PARSING

#define is_letter(c) ((((c) >='A') && ((c) <='Z')) || (((c) >='a') && ((c) <='z')))
#define is_digit(c) (((c) >='0') && ((c) <='9'))
#define is_capital(c) (((c) >='A') && ((c) <='Z'))

#define hex2int(c)  (int)((c) - (isdigit(c) ? '0' : ('A'-10)))
#define char2sign(c) ((c)=='+' ? PLUS : MINUS)
#define char2int(c) (int)((char)(c)-'0')
#define int2hex(i) (char)((i) <= 9 ? (char)(i) + '0' : (char)(i-10) + 'A')

#define string2int(s) (atoi(s->c_str()))

string itoa(u_region i) {
    stringstream ss;
    ss << (int)i;
    string s = ss.str();
    return s;
}

string itoa(int i) {
    stringstream ss;
    ss << i;
    string s = ss.str();
    return s;
}

string itona(int i, int n) {
    stringstream ss;
    ss << i;
    string s = ss.str();
    while (s.length() < n) s = "0" + s;
    return s;
}

string ERROR(string s) { return string("ERROR: ") + s + string("\n"); }

// gets next int in s, starting at i, i points to next position
int next_int(string s, int *i) {
    int sign = 1, a = 0;
    while (s[*i] == ' ') (*i)++;
    if (s[*i] == '-') {sign=-1; (*i)++; }
    while (s[*i] == ' ') (*i)++;
    if (is_letter(s[*i])) return sign; // ?
    while ((is_digit(s[*i])) && (*i <= s.length())) a = a*10 + char2int(s[(*i)++]);
    return sign*a;
}

// returns the two that are in common, e.g. 3,4,5,4 -> 4
#define common2(a,b,c,d) ((((a)==(b)) || ((a)==(c)) || ((a)==(d))) ? (a) : ((((b)==(c)) || ((b)==(d))) ? (b) : (((c)==(d)) ? (c) : 0xFF)))

#define ON_INTERVAL(x,a,b) ( ((x) >= (a)) && ((x) <= (b)) )

// BIT TWEAKS

#define BIT(i,n) ( ((i) >> (n)) & 1) // Get specific bit from int

#define remove_bits(a, b) (((a)|(b))^(b))

//#define DELETE_BIT(b,i) ((b) = ( ( (b)&firstBitsFull[i] ) | ( (b>>1)&(~firstBitsFull[i]) ) ))
//#define DELETE_BITS(b,i,n) ((b) = ( ( (b)&firstBitsFull[i] ) | ( (b>>n)&(~firstBitsFull[i]) ) ))

u_sign DELETE_BIT(u_sign *b, u_number i) {
    *b = ((*b & fbf_sign(i)) | ((*b >> 1) & (~fbf_sign(i))));
    return *b;
}

u_sign DELETE_BITS(u_sign *b, u_number i, u_number n) {
    *b = ((*b & fbf_sign(i)) | ((*b >> n) & (~fbf_sign(i))));
    return *b;
}

u_region DELETE_BIT(u_region *b, u_number i) {
    *b = ((*b & fbf_reg(i)) | ((*b >> 1) & (~fbf_reg(i))));
    return *b;
}

u_region DELETE_BITS(u_region *b, u_number i, u_number n) {
    *b = ((*b & fbf_reg(i)) | ((*b >> n) & (~fbf_reg(i))));
    return *b;
}


// INSERT(ABCD, 0, X) -> ABCDX
//#define INSERT_BIT(b,p,i) b = (( (b) & firstBitsFull[p]) | ( ((b) & (~firstBitsFull[p]) ) << 1) | ((i) << (p)))
//#define INSERT_BITS(b,p,n,i) b = (( (b) & firstBitsFull[p]) | ( ((b) & (~firstBitsFull[p]) ) << n) | ((i) << (p)))

u_region INSERT_BIT(u_region *b, u_number p, u_region i) {
    *b =  (( (*b) & fbf_reg(p)) | ( ((*b) & (~fbf_reg(p)) ) << 1) | ((i) << (p)));
    return *b;
}

u_region INSERT_BITS(u_region *b, u_number p, u_number n, u_region i) {
    *b =  (( (*b) & fbf_reg(p)) | ( ((*b) & (~fbf_reg(p)) ) << n) | ((i) << (p)));
    return *b;
}

//#define SUBSET(a,b) (((a) & (b)) == (a))

bool SUBSET(u64 a, u64 b) { return (a & b) == a; }
bool SUBSET(u128 a, u128 b) { return (a & b) == a; }
bool SUBSET(u256 a, u256 b) { return (a & b) == a; }
bool SUBSET(u512 a, u512 b) { return (a & b) == a; }

#define SETBIT(a,i,b) a = ( ((a) ^ ( BIT(a,i) << (i) )) | ((b) << (i)) )

string bit_set(u64 b) {// for printing in region style
    string s;
    for (int i=0;i<=63;i++) if (BIT(b,i)) s.append(1,int2hex(i));
    return s;
}

string bit_set(u128 b) {// for printing in region style
    string s;
    for (int i=0;i<=127;i++) if (BIT(b,i)) s.append(1,int2hex(i));
    return s;
}

string bit_set(u256 b) {// for printing in region style
    string s;
    for (int i=0;i<=255;i++) if (BIT(b,i)) s.append(1,int2hex(i));
    return s;
}

string bit_set(u512 b) {// for printing in region style
    string s;
    for (int i=0;i<=512;i++) if (BIT(b,i)) s.append(1,int2hex(i));
    return s;
}

string bitstate(u_region b, u_number n) {// for printing in region style
    string s;
    for (int i=1;i<=n;i++) s.append((BIT(b,i))?"1":"0");
    return s;
}

// for printing in extended region style
string bitset_extended(u_region b) {
    string s;
    s+= "(";
    bool fst = false;
    for (int i=0;i<=(8*sizeof(u_region)-1);i++)
        if (BIT(b,i)) {
            if (fst) s += " ";
            s+= itoa((int)i);
            fst = true;
        }
    s+=")";

    return s;
}

// prints in bit form
//void bitPrint(U_INT_64 i, int start, int n) { for (int j=0; j<n; j++) cout << BIT(i,j+start); }

void permuteBits(u_sign *b, u_number *t, u_number n) { // permute bit 0...n bits in b using renumeration table t
    u_sign b_ = *b;
    *b = 0;
    for (int i=0;i<=n;i++) (*b) |= (BIT(b_,i) << t[i]);
}

void permuteBits(u_sign *b, s_letter *t, u_number n) { // permute bit 0...n bits in b using renumeration table t
    u_sign b_ = *b;
    *b = 0;
    for (int i=0;i<=n;i++) (*b) |= (BIT(b_,i) << t[i]);
}

u_sign permuteBits(u_sign b, u_number *t, u_number n) { // permute bit 0...n bits in b using renumeration table t
    u_sign b_ = 0;
    for (int i=0;i<=n;i++) b_ |= (BIT(b,i) << t[i]);
    return b_;
}

void permuteBits(u_region *b, u_number *t, u_number n) { // permute bit 0...n bits in b using renumeration table t
    u_region b_ = *b;
    *b = (u_region)0;
    for (int i=0;i<=n;i++) *b |= (BIT(b_,i) << t[i]);
}

u_region permuteBits(u_region b, u_number *t, u_number n) { // permute bit 0...n bits in b using renumeration table t
    u_region b_ = (u_region)0;
    for (int i=0;i<=n;i++) b_ |= (BIT(b,i) << t[i]);
    return b_;
}

// count number of bits set to 1
static const unsigned char BitsSetTable256[256] =
{
#define B2(n) n,     n+1,     n+1,     n+2
#define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
};

#define COUNT_BITS_8_(v) (BitsSetTable256[((v) & 0xff)])
#define COUNT_BITS_16_(v) ( COUNT_BITS_8_(v) + COUNT_BITS_8_(((v) >> 8)))
#define COUNT_BITS_32_(v) ( COUNT_BITS_16_(v) + COUNT_BITS_16_(((v) >> 16)))
#define COUNT_BITS_64_(v) ( COUNT_BITS_32_(v) + COUNT_BITS_32_(((v) >> 32)))
#define COUNT_BITS_128_(v) ( COUNT_BITS_64_(v.hi) + COUNT_BITS_64_(v.lo) )
#define COUNT_BITS_256_(v) ( COUNT_BITS_128_(v.hi) + COUNT_BITS_128_(v.lo) )
#define COUNT_BITS_512_(v) ( COUNT_BITS_256_(v.hi) + COUNT_BITS_256_(v.lo) )

u_number count_bits(u64 b) { return COUNT_BITS_64_(b); }
u_number count_bits(u128 b) { return COUNT_BITS_128_(b); }
u_number count_bits(u256 b) { return COUNT_BITS_256_(b); }
u_number count_bits(u512 b) { return COUNT_BITS_512_(b); }

// which first lest significant bit is set to 1
int firstBitSet_t[256] =
{-1,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,
    2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,
    0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,7,0,1,0,2,0,
    1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,
    0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1};


u_number firstBitSet(u_sign b) { // starts counting bits with 0
    int shift = 0;
    do {
        //cout << "b = " << b << ", b&0xff = " << (b & 0xff) << endl;
        if (b & (u_sign)0xff) return firstBitSet_t[(int)(b & (u_sign)0xff)]+shift;
        b >>= 8;
        shift+=8;
    } while (b);
    return -1;
}

u_number firstBitSet(u_region b) { // starts counting bits with 0
    int shift = 0;
    do {
        //cout << "b = " << b << ", b&0xff = " << (b & 0xff) << endl;
        if (b & (u_region)0xff) return firstBitSet_t[(int)(b & (u_region)0xff)]+shift;
        b >>= 8;
        shift+=8;
    } while (b);
    return -1;
}

u_number secondBitSet(u_region b) { // starts counting bits with 0
    b ^= ((u_region)1 << firstBitSet(b));
    return firstBitSet(b);
}

u_number thirdBitSet(u_region b) { // starts counting bits with 0
    b ^= ((u_region)1 << firstBitSet(b));
    b ^= ((u_region)1 << firstBitSet(b));
    return firstBitSet(b);
}

u_number fourthBitSet(u_region b) { // starts counting bits with 0
    b ^= ((u_region)1 << firstBitSet(b));
    b ^= ((u_region)1 << firstBitSet(b));
    b ^= ((u_region)1 << firstBitSet(b));
    return firstBitSet(b);
}


int firstBitUnset(u_region b_) { // starts counting bits with 0
    return firstBitSet(~b_);
}
int firstBitUnset(u_sign b_) { // starts counting bits with 0
    return firstBitSet(~b_);
}
// table 0b, 1b, 11b, 111b, 1111b, 11111b, 111111b,...
u64 firstBitsFull_[65] = // !!!
{0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f,
    0xff, 0x1ff, 0x3ff, 0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff,
    0xffff, 0x1ffff, 0x3ffff, 0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff,
    0xffffff, 0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff, 0x3fffffff, 0x7fffffff,
    0xffffffff, 0x1ffffffff, 0x3ffffffff, 0x7ffffffff, 0xfffffffff, 0x1fffffffff, 0x3fffffffff, 0x7fffffffff,
    0xffffffffff, 0x1ffffffffff, 0x3ffffffffff, 0x7ffffffffff, 0xfffffffffff, 0x1fffffffffff, 0x3fffffffffff, 0x7fffffffffff,
    0xffffffffffff, 0x1ffffffffffff, 0x3ffffffffffff, 0x7ffffffffffff, 0xfffffffffffff, 0x1fffffffffffff, 0x3fffffffffffff, 0x7fffffffffffff,
    0xffffffffffffff, 0x1ffffffffffffff, 0x3ffffffffffffff, 0x7ffffffffffffff, 0xfffffffffffffff, 0x1fffffffffffffff, 0x3fffffffffffffff, 0x7fffffffffffffff,
    0xffffffffffffffff};


u_region circularShiftRight(u_region b, u_number s, u_number n) { // abcdefg -> fgabcde (2,6)
   // if (s > 63) cout << "ERROR CIRCULAR SHIFT" << endl;
    return (b >> s) | ((b & fbf_reg(s)) << (n-s));
}

u_region circularShiftLeft(u_region b, u_number s, u_number n) { // abcdefg -> cdefgab (2,6)
    return (b >> (n-s)) | ((b & fbf_reg(n-s)) << (s)); // or circluarShiftLeft(b,n-s,n)
}

u_region reverseBits(u_region b, u_number n) { // abc -> cba (3)
    u_region a = (u_region)0;
    b &= fbf_reg(n);
    u_number i = 0;
    while (b) { a <<= 1; a |= (b & 1); b >>= 1; i++; }
    return a << (n-i);
}

u_sign reverseBits(u_sign b, u_number n) { // abc -> cba (3)
    u_sign a = (u_sign)0;
    u_number i = 0;
    b &= fbf_sign(n);
    while (b) { a <<= 1; a |= (b & 1); b >>= 1; i++; }
    return a << (n-i);
}

/*void swap(U_INT_8 * a, U_INT_8 *b) { U_INT_8 t = *a; *a = *b; *b = t; }
void swap(U_INT_16 * a, U_INT_16 *b) { U_INT_16 t = *a; *a = *b; *b = t; }
void swap(U_INT_32 * a, U_INT_32 *b) { U_INT_32 t = *a; *a = *b; *b = t; }
void swap(U_INT_64 * a, U_INT_64 *b) { U_INT_64 t = *a; *a = *b; *b = t; }
void swap(INT_8 * a, INT_8 *b) { INT_8 t = *a; *a = *b; *b = t; }
void swap(INT_16 * a, INT_16 *b) { INT_16 t = *a; *a = *b; *b = t; }
void swap(INT_32 * a, INT_32 *b) { INT_32 t = *a; *a = *b; *b = t; }
void swap(INT_64 * a, INT_64 *b) { INT_64 t = *a; *a = *b; *b = t; }
*/

void swap(s_letter *a, s_letter *b) { s_letter t = *a; *a = *b; *b = t;}
void swap(u_number *a, u_number *b) { u_number t = *a; *a = *b; *b = t;}
void swap(int *a, int *b) { int t = *a; *a = *b; *b = t;}

u_number GCD(u_number a, u_number b){
    if (a < b) swap(&a, &b);
    u_number remainder = MOD(a,b);
    if (remainder == 0) return b;
    else return GCD(b,remainder);
}


#define SORT(a,b) if (*(a) > *(b)) swap(a,b)

//swap_signs(INT_8 * a, INT_8 *b) { INT_8 t = *a; *a = ABS(*a)*SIGN(*b); *b = ABS(*b)*SIGN(t);}

void swapBits(u_region * b,u_number i, u_number j) {
    u_region b_ = *b;
    *b &= ~(((u_region)1<<i) | ((u_region)1<<j));
    *b |= ((BIT(b_,i)<<j) | (BIT(b_,j)<<i));
}

void swapBits(u_sign * b, u_number i, u_number j) {
    u_sign b_ = *b;
    (*b) &= ~(((u_sign)1<<i) | ((u_sign)1<<j));
    (*b) |= ((BIT(b_,i)<<j) | (BIT(b_,j)<<i));
}

//{0, 0x1, 0x3 0x}
// NUMBER OF 0 LEAST SIGNIFICANT BITS
//static const int MultiplyDeBruijnBitPosition[32] =
//{ 0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9 };
//#define FIRST_BIT_SET(v) (MultiplyDeBruijnBitPosition[((uint32_t)(((v) & -(v)) * 0x077CB531U)) >> 27] +1)
// FIRST n BITS 1
//#define FIRST_FULL_BITS(n) (uint32_t)((~0 & 0x00FFFFFF) >> (24-(n)))
//#define ONES_EXCEPT_BIT(n) (~0 ^ (1 << (n-1))) //starting from 1

string INT_2_STRING(int i) {
    string s = "";
    string t;
    if (i<0) s += "-";
    if (i) {
        while (i) {
            t = "";
            t.push_back('0'+char(i%10));
            s = t + s;
            i/=10;
        }
    } else s = "0";
    return s;
}

// FOR HOMFLY SKEIN MODULES

#define MID_VAR 'M' // variable that belongs to t_0

#define HOMFLY_MAX_VAR 12

string HOMFLY_BASE_VARIABLES[HOMFLY_MAX_VAR+1] = {"M","LMN","KLMNO","JKLMNOP","IJKLMNOPQ","HIJKLMNOPQR","GHIJKLMNOPQRS","FGHIJKLMNOPQRST","EFGHIJKLMNOPQRSTU","DEFGHIJKLMNOPQRSTUV","CDEFGHIJKLMNOPQRSTUVW","BCDEFGHIJKLMNOPQRSTUVWX","ABCDEFGHIJKLMNOPQRSTUVWXY"};

string HOMFLY_BASE[256] = {
    "","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","", // 0...49
    "","","","","","","","","","","","","","","","t_-12","t_-11","t_-10","t_-9","t_-8","t_-7","t_-6","t_-5","t_-4","t_-3","t_-2","t_-1","t_0",
    "t_1","t_2","t_3","t_4","t_5","t_6","t_7","t_8","t_9","t_10","t_11","t_12","","","","","","","","","","","", // ? ... 99
    "","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","", // ...149
    "","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","", // ...199
    "","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","", // ...249
    "","","","",""};
string HOMFLY_BASE_LATEX[256] = {
    "","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","", // 0...49
    "","","","","","","","","","","","","","","","t_{-12}","t_{-11}","t_{-10}","t_{-9}","t_{-8}","t_{-7}","t_{-6}","t_{-5}","t_{-4}","t_{-3}","t_{-2}","t_{-1}","t_{0}",
    "t_{1}","t_{2}","t_{3}","t_{4}","t_{5}","t_{6}","t_{7}","t_{8}","t_{9}","t_{10}","t_{11}","t_{12}","","","","","","","","","","","", // ? ... 99
    "","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","", // ...149
    "","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","", // ...199
    "","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","", // ...249
    "","","","",""};



string n_space(int i) {
    string s = "";
    cout << "(" << i << ") ";
  //  s += "(" + atoi(i) << ") ";
    for (int k=0;k<i;k++) s += "  ";
    return s;
}


#define b0 0x0
#define b1 0x1
#define b00 0x0
#define b01 0x1
#define b10 0x2
#define b11 0x3
#define b000 0x0
#define b001 0x1
#define b010 0x2
#define b011 0x3
#define b100 0x4
#define b101 0x5
#define b110 0x6
#define b111 0x7
#define b0000 0x0
#define b0001 0x1
#define b0010 0x2
#define b0011 0x3
#define b0100 0x4
#define b0101 0x5
#define b0110 0x6
#define b0111 0x7
#define b1000 0x8
#define b1001 0x9
#define b1010 0xa
#define b1011 0xb
#define b1100 0xc
#define b1101 0xd
#define b1110 0xe
#define b1111 0xf
#define b00000 0x0
#define b00001 0x1
#define b00010 0x2
#define b00011 0x3
#define b00100 0x4
#define b00101 0x5
#define b00110 0x6
#define b00111 0x7
#define b01000 0x8
#define b01001 0x9
#define b01010 0xa
#define b01011 0xb
#define b01100 0xc
#define b01101 0xd
#define b01110 0xe
#define b01111 0xf
#define b10000 0x10
#define b10001 0x11
#define b10010 0x12
#define b10011 0x13
#define b10100 0x14
#define b10101 0x15
#define b10110 0x16
#define b10111 0x17
#define b11000 0x18
#define b11001 0x19
#define b11010 0x1a
#define b11011 0x1b
#define b11100 0x1c
#define b11101 0x1d
#define b11110 0x1e
#define b11111 0x1f
#define b000000 0x0
#define b000001 0x1
#define b000010 0x2
#define b000011 0x3
#define b000100 0x4
#define b000101 0x5
#define b000110 0x6
#define b000111 0x7
#define b001000 0x8
#define b001001 0x9
#define b001010 0xa
#define b001011 0xb
#define b001100 0xc
#define b001101 0xd
#define b001110 0xe
#define b001111 0xf
#define b010000 0x10
#define b010001 0x11
#define b010010 0x12
#define b010011 0x13
#define b010100 0x14
#define b010101 0x15
#define b010110 0x16
#define b010111 0x17
#define b011000 0x18
#define b011001 0x19
#define b011010 0x1a
#define b011011 0x1b
#define b011100 0x1c
#define b011101 0x1d
#define b011110 0x1e
#define b011111 0x1f
#define b100000 0x20
#define b100001 0x21
#define b100010 0x22
#define b100011 0x23
#define b100100 0x24
#define b100101 0x25
#define b100110 0x26
#define b100111 0x27
#define b101000 0x28
#define b101001 0x29
#define b101010 0x2a
#define b101011 0x2b
#define b101100 0x2c
#define b101101 0x2d
#define b101110 0x2e
#define b101111 0x2f
#define b110000 0x30
#define b110001 0x31
#define b110010 0x32
#define b110011 0x33
#define b110100 0x34
#define b110101 0x35
#define b110110 0x36
#define b110111 0x37
#define b111000 0x38
#define b111001 0x39
#define b111010 0x3a
#define b111011 0x3b
#define b111100 0x3c
#define b111101 0x3d
#define b111110 0x3e
#define b111111 0x3f
#define b0000000 0x0
#define b0000001 0x1
#define b0000010 0x2
#define b0000011 0x3
#define b0000100 0x4
#define b0000101 0x5
#define b0000110 0x6
#define b0000111 0x7
#define b0001000 0x8
#define b0001001 0x9
#define b0001010 0xa
#define b0001011 0xb
#define b0001100 0xc
#define b0001101 0xd
#define b0001110 0xe
#define b0001111 0xf
#define b0010000 0x10
#define b0010001 0x11
#define b0010010 0x12
#define b0010011 0x13
#define b0010100 0x14
#define b0010101 0x15
#define b0010110 0x16
#define b0010111 0x17
#define b0011000 0x18
#define b0011001 0x19
#define b0011010 0x1a
#define b0011011 0x1b
#define b0011100 0x1c
#define b0011101 0x1d
#define b0011110 0x1e
#define b0011111 0x1f
#define b0100000 0x20
#define b0100001 0x21
#define b0100010 0x22
#define b0100011 0x23
#define b0100100 0x24
#define b0100101 0x25
#define b0100110 0x26
#define b0100111 0x27
#define b0101000 0x28
#define b0101001 0x29
#define b0101010 0x2a
#define b0101011 0x2b
#define b0101100 0x2c
#define b0101101 0x2d
#define b0101110 0x2e
#define b0101111 0x2f
#define b0110000 0x30
#define b0110001 0x31
#define b0110010 0x32
#define b0110011 0x33
#define b0110100 0x34
#define b0110101 0x35
#define b0110110 0x36
#define b0110111 0x37
#define b0111000 0x38
#define b0111001 0x39
#define b0111010 0x3a
#define b0111011 0x3b
#define b0111100 0x3c
#define b0111101 0x3d
#define b0111110 0x3e
#define b0111111 0x3f
#define b1000000 0x40
#define b1000001 0x41
#define b1000010 0x42
#define b1000011 0x43
#define b1000100 0x44
#define b1000101 0x45
#define b1000110 0x46
#define b1000111 0x47
#define b1001000 0x48
#define b1001001 0x49
#define b1001010 0x4a
#define b1001011 0x4b
#define b1001100 0x4c
#define b1001101 0x4d
#define b1001110 0x4e
#define b1001111 0x4f
#define b1010000 0x50
#define b1010001 0x51
#define b1010010 0x52
#define b1010011 0x53
#define b1010100 0x54
#define b1010101 0x55
#define b1010110 0x56
#define b1010111 0x57
#define b1011000 0x58
#define b1011001 0x59
#define b1011010 0x5a
#define b1011011 0x5b
#define b1011100 0x5c
#define b1011101 0x5d
#define b1011110 0x5e
#define b1011111 0x5f
#define b1100000 0x60
#define b1100001 0x61
#define b1100010 0x62
#define b1100011 0x63
#define b1100100 0x64
#define b1100101 0x65
#define b1100110 0x66
#define b1100111 0x67
#define b1101000 0x68
#define b1101001 0x69
#define b1101010 0x6a
#define b1101011 0x6b
#define b1101100 0x6c
#define b1101101 0x6d
#define b1101110 0x6e
#define b1101111 0x6f
#define b1110000 0x70
#define b1110001 0x71
#define b1110010 0x72
#define b1110011 0x73
#define b1110100 0x74
#define b1110101 0x75
#define b1110110 0x76
#define b1110111 0x77
#define b1111000 0x78
#define b1111001 0x79
#define b1111010 0x7a
#define b1111011 0x7b
#define b1111100 0x7c
#define b1111101 0x7d
#define b1111110 0x7e
#define b1111111 0x7f
#define b00000000 0x0
#define b00000001 0x1
#define b00000010 0x2
#define b00000011 0x3
#define b00000100 0x4
#define b00000101 0x5
#define b00000110 0x6
#define b00000111 0x7
#define b00001000 0x8
#define b00001001 0x9
#define b00001010 0xa
#define b00001011 0xb
#define b00001100 0xc
#define b00001101 0xd
#define b00001110 0xe
#define b00001111 0xf
#define b00010000 0x10
#define b00010001 0x11
#define b00010010 0x12
#define b00010011 0x13
#define b00010100 0x14
#define b00010101 0x15
#define b00010110 0x16
#define b00010111 0x17
#define b00011000 0x18
#define b00011001 0x19
#define b00011010 0x1a
#define b00011011 0x1b
#define b00011100 0x1c
#define b00011101 0x1d
#define b00011110 0x1e
#define b00011111 0x1f
#define b00100000 0x20
#define b00100001 0x21
#define b00100010 0x22
#define b00100011 0x23
#define b00100100 0x24
#define b00100101 0x25
#define b00100110 0x26
#define b00100111 0x27
#define b00101000 0x28
#define b00101001 0x29
#define b00101010 0x2a
#define b00101011 0x2b
#define b00101100 0x2c
#define b00101101 0x2d
#define b00101110 0x2e
#define b00101111 0x2f
#define b00110000 0x30
#define b00110001 0x31
#define b00110010 0x32
#define b00110011 0x33
#define b00110100 0x34
#define b00110101 0x35
#define b00110110 0x36
#define b00110111 0x37
#define b00111000 0x38
#define b00111001 0x39
#define b00111010 0x3a
#define b00111011 0x3b
#define b00111100 0x3c
#define b00111101 0x3d
#define b00111110 0x3e
#define b00111111 0x3f
#define b01000000 0x40
#define b01000001 0x41
#define b01000010 0x42
#define b01000011 0x43
#define b01000100 0x44
#define b01000101 0x45
#define b01000110 0x46
#define b01000111 0x47
#define b01001000 0x48
#define b01001001 0x49
#define b01001010 0x4a
#define b01001011 0x4b
#define b01001100 0x4c
#define b01001101 0x4d
#define b01001110 0x4e
#define b01001111 0x4f
#define b01010000 0x50
#define b01010001 0x51
#define b01010010 0x52
#define b01010011 0x53
#define b01010100 0x54
#define b01010101 0x55
#define b01010110 0x56
#define b01010111 0x57
#define b01011000 0x58
#define b01011001 0x59
#define b01011010 0x5a
#define b01011011 0x5b
#define b01011100 0x5c
#define b01011101 0x5d
#define b01011110 0x5e
#define b01011111 0x5f
#define b01100000 0x60
#define b01100001 0x61
#define b01100010 0x62
#define b01100011 0x63
#define b01100100 0x64
#define b01100101 0x65
#define b01100110 0x66
#define b01100111 0x67
#define b01101000 0x68
#define b01101001 0x69
#define b01101010 0x6a
#define b01101011 0x6b
#define b01101100 0x6c
#define b01101101 0x6d
#define b01101110 0x6e
#define b01101111 0x6f
#define b01110000 0x70
#define b01110001 0x71
#define b01110010 0x72
#define b01110011 0x73
#define b01110100 0x74
#define b01110101 0x75
#define b01110110 0x76
#define b01110111 0x77
#define b01111000 0x78
#define b01111001 0x79
#define b01111010 0x7a
#define b01111011 0x7b
#define b01111100 0x7c
#define b01111101 0x7d
#define b01111110 0x7e
#define b01111111 0x7f
#define b10000000 0x80
#define b10000001 0x81
#define b10000010 0x82
#define b10000011 0x83
#define b10000100 0x84
#define b10000101 0x85
#define b10000110 0x86
#define b10000111 0x87
#define b10001000 0x88
#define b10001001 0x89
#define b10001010 0x8a
#define b10001011 0x8b
#define b10001100 0x8c
#define b10001101 0x8d
#define b10001110 0x8e
#define b10001111 0x8f
#define b10010000 0x90
#define b10010001 0x91
#define b10010010 0x92
#define b10010011 0x93
#define b10010100 0x94
#define b10010101 0x95
#define b10010110 0x96
#define b10010111 0x97
#define b10011000 0x98
#define b10011001 0x99
#define b10011010 0x9a
#define b10011011 0x9b
#define b10011100 0x9c
#define b10011101 0x9d
#define b10011110 0x9e
#define b10011111 0x9f
#define b10100000 0xa0
#define b10100001 0xa1
#define b10100010 0xa2
#define b10100011 0xa3
#define b10100100 0xa4
#define b10100101 0xa5
#define b10100110 0xa6
#define b10100111 0xa7
#define b10101000 0xa8
#define b10101001 0xa9
#define b10101010 0xaa
#define b10101011 0xab
#define b10101100 0xac
#define b10101101 0xad
#define b10101110 0xae
#define b10101111 0xaf
#define b10110000 0xb0
#define b10110001 0xb1
#define b10110010 0xb2
#define b10110011 0xb3
#define b10110100 0xb4
#define b10110101 0xb5
#define b10110110 0xb6
#define b10110111 0xb7
#define b10111000 0xb8
#define b10111001 0xb9
#define b10111010 0xba
#define b10111011 0xbb
#define b10111100 0xbc
#define b10111101 0xbd
#define b10111110 0xbe
#define b10111111 0xbf
#define b11000000 0xc0
#define b11000001 0xc1
#define b11000010 0xc2
#define b11000011 0xc3
#define b11000100 0xc4
#define b11000101 0xc5
#define b11000110 0xc6
#define b11000111 0xc7
#define b11001000 0xc8
#define b11001001 0xc9
#define b11001010 0xca
#define b11001011 0xcb
#define b11001100 0xcc
#define b11001101 0xcd
#define b11001110 0xce
#define b11001111 0xcf
#define b11010000 0xd0
#define b11010001 0xd1
#define b11010010 0xd2
#define b11010011 0xd3
#define b11010100 0xd4
#define b11010101 0xd5
#define b11010110 0xd6
#define b11010111 0xd7
#define b11011000 0xd8
#define b11011001 0xd9
#define b11011010 0xda
#define b11011011 0xdb
#define b11011100 0xdc
#define b11011101 0xdd
#define b11011110 0xde
#define b11011111 0xdf
#define b11100000 0xe0
#define b11100001 0xe1
#define b11100010 0xe2
#define b11100011 0xe3
#define b11100100 0xe4
#define b11100101 0xe5
#define b11100110 0xe6
#define b11100111 0xe7
#define b11101000 0xe8
#define b11101001 0xe9
#define b11101010 0xea
#define b11101011 0xeb
#define b11101100 0xec
#define b11101101 0xed
#define b11101110 0xee
#define b11101111 0xef
#define b11110000 0xf0
#define b11110001 0xf1
#define b11110010 0xf2
#define b11110011 0xf3
#define b11110100 0xf4
#define b11110101 0xf5
#define b11110110 0xf6
#define b11110111 0xf7
#define b11111000 0xf8
#define b11111001 0xf9
#define b11111010 0xfa
#define b11111011 0xfb
#define b11111100 0xfc
#define b11111101 0xfd
#define b11111110 0xfe
#define b11111111 0xff
#define b000000000 0x0
#define b000000001 0x1
#define b000000010 0x2
#define b000000011 0x3
#define b000000100 0x4
#define b000000101 0x5
#define b000000110 0x6
#define b000000111 0x7
#define b000001000 0x8
#define b000001001 0x9
#define b000001010 0xa
#define b000001011 0xb
#define b000001100 0xc
#define b000001101 0xd
#define b000001110 0xe
#define b000001111 0xf
#define b000010000 0x10
#define b000010001 0x11
#define b000010010 0x12
#define b000010011 0x13
#define b000010100 0x14
#define b000010101 0x15
#define b000010110 0x16
#define b000010111 0x17
#define b000011000 0x18
#define b000011001 0x19
#define b000011010 0x1a
#define b000011011 0x1b
#define b000011100 0x1c
#define b000011101 0x1d
#define b000011110 0x1e
#define b000011111 0x1f
#define b000100000 0x20
#define b000100001 0x21
#define b000100010 0x22
#define b000100011 0x23
#define b000100100 0x24
#define b000100101 0x25
#define b000100110 0x26
#define b000100111 0x27
#define b000101000 0x28
#define b000101001 0x29
#define b000101010 0x2a
#define b000101011 0x2b
#define b000101100 0x2c
#define b000101101 0x2d
#define b000101110 0x2e
#define b000101111 0x2f
#define b000110000 0x30
#define b000110001 0x31
#define b000110010 0x32
#define b000110011 0x33
#define b000110100 0x34
#define b000110101 0x35
#define b000110110 0x36
#define b000110111 0x37
#define b000111000 0x38
#define b000111001 0x39
#define b000111010 0x3a
#define b000111011 0x3b
#define b000111100 0x3c
#define b000111101 0x3d
#define b000111110 0x3e
#define b000111111 0x3f
#define b001000000 0x40
#define b001000001 0x41
#define b001000010 0x42
#define b001000011 0x43
#define b001000100 0x44
#define b001000101 0x45
#define b001000110 0x46
#define b001000111 0x47
#define b001001000 0x48
#define b001001001 0x49
#define b001001010 0x4a
#define b001001011 0x4b
#define b001001100 0x4c
#define b001001101 0x4d
#define b001001110 0x4e
#define b001001111 0x4f
#define b001010000 0x50
#define b001010001 0x51
#define b001010010 0x52
#define b001010011 0x53
#define b001010100 0x54
#define b001010101 0x55
#define b001010110 0x56
#define b001010111 0x57
#define b001011000 0x58
#define b001011001 0x59
#define b001011010 0x5a
#define b001011011 0x5b
#define b001011100 0x5c
#define b001011101 0x5d
#define b001011110 0x5e
#define b001011111 0x5f
#define b001100000 0x60
#define b001100001 0x61
#define b001100010 0x62
#define b001100011 0x63
#define b001100100 0x64
#define b001100101 0x65
#define b001100110 0x66
#define b001100111 0x67
#define b001101000 0x68
#define b001101001 0x69
#define b001101010 0x6a
#define b001101011 0x6b
#define b001101100 0x6c
#define b001101101 0x6d
#define b001101110 0x6e
#define b001101111 0x6f
#define b001110000 0x70
#define b001110001 0x71
#define b001110010 0x72
#define b001110011 0x73
#define b001110100 0x74
#define b001110101 0x75
#define b001110110 0x76
#define b001110111 0x77
#define b001111000 0x78
#define b001111001 0x79
#define b001111010 0x7a
#define b001111011 0x7b
#define b001111100 0x7c
#define b001111101 0x7d
#define b001111110 0x7e
#define b001111111 0x7f
#define b010000000 0x80
#define b010000001 0x81
#define b010000010 0x82
#define b010000011 0x83
#define b010000100 0x84
#define b010000101 0x85
#define b010000110 0x86
#define b010000111 0x87
#define b010001000 0x88
#define b010001001 0x89
#define b010001010 0x8a
#define b010001011 0x8b
#define b010001100 0x8c
#define b010001101 0x8d
#define b010001110 0x8e
#define b010001111 0x8f
#define b010010000 0x90
#define b010010001 0x91
#define b010010010 0x92
#define b010010011 0x93
#define b010010100 0x94
#define b010010101 0x95
#define b010010110 0x96
#define b010010111 0x97
#define b010011000 0x98
#define b010011001 0x99
#define b010011010 0x9a
#define b010011011 0x9b
#define b010011100 0x9c
#define b010011101 0x9d
#define b010011110 0x9e
#define b010011111 0x9f
#define b010100000 0xa0
#define b010100001 0xa1
#define b010100010 0xa2
#define b010100011 0xa3
#define b010100100 0xa4
#define b010100101 0xa5
#define b010100110 0xa6
#define b010100111 0xa7
#define b010101000 0xa8
#define b010101001 0xa9
#define b010101010 0xaa
#define b010101011 0xab
#define b010101100 0xac
#define b010101101 0xad
#define b010101110 0xae
#define b010101111 0xaf
#define b010110000 0xb0
#define b010110001 0xb1
#define b010110010 0xb2
#define b010110011 0xb3
#define b010110100 0xb4
#define b010110101 0xb5
#define b010110110 0xb6
#define b010110111 0xb7
#define b010111000 0xb8
#define b010111001 0xb9
#define b010111010 0xba
#define b010111011 0xbb
#define b010111100 0xbc
#define b010111101 0xbd
#define b010111110 0xbe
#define b010111111 0xbf
#define b011000000 0xc0
#define b011000001 0xc1
#define b011000010 0xc2
#define b011000011 0xc3
#define b011000100 0xc4
#define b011000101 0xc5
#define b011000110 0xc6
#define b011000111 0xc7
#define b011001000 0xc8
#define b011001001 0xc9
#define b011001010 0xca
#define b011001011 0xcb
#define b011001100 0xcc
#define b011001101 0xcd
#define b011001110 0xce
#define b011001111 0xcf
#define b011010000 0xd0
#define b011010001 0xd1
#define b011010010 0xd2
#define b011010011 0xd3
#define b011010100 0xd4
#define b011010101 0xd5
#define b011010110 0xd6
#define b011010111 0xd7
#define b011011000 0xd8
#define b011011001 0xd9
#define b011011010 0xda
#define b011011011 0xdb
#define b011011100 0xdc
#define b011011101 0xdd
#define b011011110 0xde
#define b011011111 0xdf
#define b011100000 0xe0
#define b011100001 0xe1
#define b011100010 0xe2
#define b011100011 0xe3
#define b011100100 0xe4
#define b011100101 0xe5
#define b011100110 0xe6
#define b011100111 0xe7
#define b011101000 0xe8
#define b011101001 0xe9
#define b011101010 0xea
#define b011101011 0xeb
#define b011101100 0xec
#define b011101101 0xed
#define b011101110 0xee
#define b011101111 0xef
#define b011110000 0xf0
#define b011110001 0xf1
#define b011110010 0xf2
#define b011110011 0xf3
#define b011110100 0xf4
#define b011110101 0xf5
#define b011110110 0xf6
#define b011110111 0xf7
#define b011111000 0xf8
#define b011111001 0xf9
#define b011111010 0xfa
#define b011111011 0xfb
#define b011111100 0xfc
#define b011111101 0xfd
#define b011111110 0xfe
#define b011111111 0xff
#define b100000000 0x100
#define b100000001 0x101
#define b100000010 0x102
#define b100000011 0x103
#define b100000100 0x104
#define b100000101 0x105
#define b100000110 0x106
#define b100000111 0x107
#define b100001000 0x108
#define b100001001 0x109
#define b100001010 0x10a
#define b100001011 0x10b
#define b100001100 0x10c
#define b100001101 0x10d
#define b100001110 0x10e
#define b100001111 0x10f
#define b100010000 0x110
#define b100010001 0x111
#define b100010010 0x112
#define b100010011 0x113
#define b100010100 0x114
#define b100010101 0x115
#define b100010110 0x116
#define b100010111 0x117
#define b100011000 0x118
#define b100011001 0x119
#define b100011010 0x11a
#define b100011011 0x11b
#define b100011100 0x11c
#define b100011101 0x11d
#define b100011110 0x11e
#define b100011111 0x11f
#define b100100000 0x120
#define b100100001 0x121
#define b100100010 0x122
#define b100100011 0x123
#define b100100100 0x124
#define b100100101 0x125
#define b100100110 0x126
#define b100100111 0x127
#define b100101000 0x128
#define b100101001 0x129
#define b100101010 0x12a
#define b100101011 0x12b
#define b100101100 0x12c
#define b100101101 0x12d
#define b100101110 0x12e
#define b100101111 0x12f
#define b100110000 0x130
#define b100110001 0x131
#define b100110010 0x132
#define b100110011 0x133
#define b100110100 0x134
#define b100110101 0x135
#define b100110110 0x136
#define b100110111 0x137
#define b100111000 0x138
#define b100111001 0x139
#define b100111010 0x13a
#define b100111011 0x13b
#define b100111100 0x13c
#define b100111101 0x13d
#define b100111110 0x13e
#define b100111111 0x13f
#define b101000000 0x140
#define b101000001 0x141
#define b101000010 0x142
#define b101000011 0x143
#define b101000100 0x144
#define b101000101 0x145
#define b101000110 0x146
#define b101000111 0x147
#define b101001000 0x148
#define b101001001 0x149
#define b101001010 0x14a
#define b101001011 0x14b
#define b101001100 0x14c
#define b101001101 0x14d
#define b101001110 0x14e
#define b101001111 0x14f
#define b101010000 0x150
#define b101010001 0x151
#define b101010010 0x152
#define b101010011 0x153
#define b101010100 0x154
#define b101010101 0x155
#define b101010110 0x156
#define b101010111 0x157
#define b101011000 0x158
#define b101011001 0x159
#define b101011010 0x15a
#define b101011011 0x15b
#define b101011100 0x15c
#define b101011101 0x15d
#define b101011110 0x15e
#define b101011111 0x15f
#define b101100000 0x160
#define b101100001 0x161
#define b101100010 0x162
#define b101100011 0x163
#define b101100100 0x164
#define b101100101 0x165
#define b101100110 0x166
#define b101100111 0x167
#define b101101000 0x168
#define b101101001 0x169
#define b101101010 0x16a
#define b101101011 0x16b
#define b101101100 0x16c
#define b101101101 0x16d
#define b101101110 0x16e
#define b101101111 0x16f
#define b101110000 0x170
#define b101110001 0x171
#define b101110010 0x172
#define b101110011 0x173
#define b101110100 0x174
#define b101110101 0x175
#define b101110110 0x176
#define b101110111 0x177
#define b101111000 0x178
#define b101111001 0x179
#define b101111010 0x17a
#define b101111011 0x17b
#define b101111100 0x17c
#define b101111101 0x17d
#define b101111110 0x17e
#define b101111111 0x17f
#define b110000000 0x180
#define b110000001 0x181
#define b110000010 0x182
#define b110000011 0x183
#define b110000100 0x184
#define b110000101 0x185
#define b110000110 0x186
#define b110000111 0x187
#define b110001000 0x188
#define b110001001 0x189
#define b110001010 0x18a
#define b110001011 0x18b
#define b110001100 0x18c
#define b110001101 0x18d
#define b110001110 0x18e
#define b110001111 0x18f
#define b110010000 0x190
#define b110010001 0x191
#define b110010010 0x192
#define b110010011 0x193
#define b110010100 0x194
#define b110010101 0x195
#define b110010110 0x196
#define b110010111 0x197
#define b110011000 0x198
#define b110011001 0x199
#define b110011010 0x19a
#define b110011011 0x19b
#define b110011100 0x19c
#define b110011101 0x19d
#define b110011110 0x19e
#define b110011111 0x19f
#define b110100000 0x1a0
#define b110100001 0x1a1
#define b110100010 0x1a2
#define b110100011 0x1a3
#define b110100100 0x1a4
#define b110100101 0x1a5
#define b110100110 0x1a6
#define b110100111 0x1a7
#define b110101000 0x1a8
#define b110101001 0x1a9
#define b110101010 0x1aa
#define b110101011 0x1ab
#define b110101100 0x1ac
#define b110101101 0x1ad
#define b110101110 0x1ae
#define b110101111 0x1af
#define b110110000 0x1b0
#define b110110001 0x1b1
#define b110110010 0x1b2
#define b110110011 0x1b3
#define b110110100 0x1b4
#define b110110101 0x1b5
#define b110110110 0x1b6
#define b110110111 0x1b7
#define b110111000 0x1b8
#define b110111001 0x1b9
#define b110111010 0x1ba
#define b110111011 0x1bb
#define b110111100 0x1bc
#define b110111101 0x1bd
#define b110111110 0x1be
#define b110111111 0x1bf
#define b111000000 0x1c0
#define b111000001 0x1c1
#define b111000010 0x1c2
#define b111000011 0x1c3
#define b111000100 0x1c4
#define b111000101 0x1c5
#define b111000110 0x1c6
#define b111000111 0x1c7
#define b111001000 0x1c8
#define b111001001 0x1c9
#define b111001010 0x1ca
#define b111001011 0x1cb
#define b111001100 0x1cc
#define b111001101 0x1cd
#define b111001110 0x1ce
#define b111001111 0x1cf
#define b111010000 0x1d0
#define b111010001 0x1d1
#define b111010010 0x1d2
#define b111010011 0x1d3
#define b111010100 0x1d4
#define b111010101 0x1d5
#define b111010110 0x1d6
#define b111010111 0x1d7
#define b111011000 0x1d8
#define b111011001 0x1d9
#define b111011010 0x1da
#define b111011011 0x1db
#define b111011100 0x1dc
#define b111011101 0x1dd
#define b111011110 0x1de
#define b111011111 0x1df
#define b111100000 0x1e0
#define b111100001 0x1e1
#define b111100010 0x1e2
#define b111100011 0x1e3
#define b111100100 0x1e4
#define b111100101 0x1e5
#define b111100110 0x1e6
#define b111100111 0x1e7
#define b111101000 0x1e8
#define b111101001 0x1e9
#define b111101010 0x1ea
#define b111101011 0x1eb
#define b111101100 0x1ec
#define b111101101 0x1ed
#define b111101110 0x1ee
#define b111101111 0x1ef
#define b111110000 0x1f0
#define b111110001 0x1f1
#define b111110010 0x1f2
#define b111110011 0x1f3
#define b111110100 0x1f4
#define b111110101 0x1f5
#define b111110110 0x1f6
#define b111110111 0x1f7
#define b111111000 0x1f8
#define b111111001 0x1f9
#define b111111010 0x1fa
#define b111111011 0x1fb
#define b111111100 0x1fc
#define b111111101 0x1fd
#define b111111110 0x1fe
#define b111111111 0x1ff
#define b0000000000 0x0
#define b0000000001 0x1
#define b0000000010 0x2
#define b0000000011 0x3
#define b0000000100 0x4
#define b0000000101 0x5
#define b0000000110 0x6
#define b0000000111 0x7
#define b0000001000 0x8
#define b0000001001 0x9
#define b0000001010 0xa
#define b0000001011 0xb
#define b0000001100 0xc
#define b0000001101 0xd
#define b0000001110 0xe
#define b0000001111 0xf
#define b0000010000 0x10
#define b0000010001 0x11
#define b0000010010 0x12
#define b0000010011 0x13
#define b0000010100 0x14
#define b0000010101 0x15
#define b0000010110 0x16
#define b0000010111 0x17
#define b0000011000 0x18
#define b0000011001 0x19
#define b0000011010 0x1a
#define b0000011011 0x1b
#define b0000011100 0x1c
#define b0000011101 0x1d
#define b0000011110 0x1e
#define b0000011111 0x1f
#define b0000100000 0x20
#define b0000100001 0x21
#define b0000100010 0x22
#define b0000100011 0x23
#define b0000100100 0x24
#define b0000100101 0x25
#define b0000100110 0x26
#define b0000100111 0x27
#define b0000101000 0x28
#define b0000101001 0x29
#define b0000101010 0x2a
#define b0000101011 0x2b
#define b0000101100 0x2c
#define b0000101101 0x2d
#define b0000101110 0x2e
#define b0000101111 0x2f
#define b0000110000 0x30
#define b0000110001 0x31
#define b0000110010 0x32
#define b0000110011 0x33
#define b0000110100 0x34
#define b0000110101 0x35
#define b0000110110 0x36
#define b0000110111 0x37
#define b0000111000 0x38
#define b0000111001 0x39
#define b0000111010 0x3a
#define b0000111011 0x3b
#define b0000111100 0x3c
#define b0000111101 0x3d
#define b0000111110 0x3e
#define b0000111111 0x3f
#define b0001000000 0x40
#define b0001000001 0x41
#define b0001000010 0x42
#define b0001000011 0x43
#define b0001000100 0x44
#define b0001000101 0x45
#define b0001000110 0x46
#define b0001000111 0x47
#define b0001001000 0x48
#define b0001001001 0x49
#define b0001001010 0x4a
#define b0001001011 0x4b
#define b0001001100 0x4c
#define b0001001101 0x4d
#define b0001001110 0x4e
#define b0001001111 0x4f
#define b0001010000 0x50
#define b0001010001 0x51
#define b0001010010 0x52
#define b0001010011 0x53
#define b0001010100 0x54
#define b0001010101 0x55
#define b0001010110 0x56
#define b0001010111 0x57
#define b0001011000 0x58
#define b0001011001 0x59
#define b0001011010 0x5a
#define b0001011011 0x5b
#define b0001011100 0x5c
#define b0001011101 0x5d
#define b0001011110 0x5e
#define b0001011111 0x5f
#define b0001100000 0x60
#define b0001100001 0x61
#define b0001100010 0x62
#define b0001100011 0x63
#define b0001100100 0x64
#define b0001100101 0x65
#define b0001100110 0x66
#define b0001100111 0x67
#define b0001101000 0x68
#define b0001101001 0x69
#define b0001101010 0x6a
#define b0001101011 0x6b
#define b0001101100 0x6c
#define b0001101101 0x6d
#define b0001101110 0x6e
#define b0001101111 0x6f
#define b0001110000 0x70
#define b0001110001 0x71
#define b0001110010 0x72
#define b0001110011 0x73
#define b0001110100 0x74
#define b0001110101 0x75
#define b0001110110 0x76
#define b0001110111 0x77
#define b0001111000 0x78
#define b0001111001 0x79
#define b0001111010 0x7a
#define b0001111011 0x7b
#define b0001111100 0x7c
#define b0001111101 0x7d
#define b0001111110 0x7e
#define b0001111111 0x7f
#define b0010000000 0x80
#define b0010000001 0x81
#define b0010000010 0x82
#define b0010000011 0x83
#define b0010000100 0x84
#define b0010000101 0x85
#define b0010000110 0x86
#define b0010000111 0x87
#define b0010001000 0x88
#define b0010001001 0x89
#define b0010001010 0x8a
#define b0010001011 0x8b
#define b0010001100 0x8c
#define b0010001101 0x8d
#define b0010001110 0x8e
#define b0010001111 0x8f
#define b0010010000 0x90
#define b0010010001 0x91
#define b0010010010 0x92
#define b0010010011 0x93
#define b0010010100 0x94
#define b0010010101 0x95
#define b0010010110 0x96
#define b0010010111 0x97
#define b0010011000 0x98
#define b0010011001 0x99
#define b0010011010 0x9a
#define b0010011011 0x9b
#define b0010011100 0x9c
#define b0010011101 0x9d
#define b0010011110 0x9e
#define b0010011111 0x9f
#define b0010100000 0xa0
#define b0010100001 0xa1
#define b0010100010 0xa2
#define b0010100011 0xa3
#define b0010100100 0xa4
#define b0010100101 0xa5
#define b0010100110 0xa6
#define b0010100111 0xa7
#define b0010101000 0xa8
#define b0010101001 0xa9
#define b0010101010 0xaa
#define b0010101011 0xab
#define b0010101100 0xac
#define b0010101101 0xad
#define b0010101110 0xae
#define b0010101111 0xaf
#define b0010110000 0xb0
#define b0010110001 0xb1
#define b0010110010 0xb2
#define b0010110011 0xb3
#define b0010110100 0xb4
#define b0010110101 0xb5
#define b0010110110 0xb6
#define b0010110111 0xb7
#define b0010111000 0xb8
#define b0010111001 0xb9
#define b0010111010 0xba
#define b0010111011 0xbb
#define b0010111100 0xbc
#define b0010111101 0xbd
#define b0010111110 0xbe
#define b0010111111 0xbf
#define b0011000000 0xc0
#define b0011000001 0xc1
#define b0011000010 0xc2
#define b0011000011 0xc3
#define b0011000100 0xc4
#define b0011000101 0xc5
#define b0011000110 0xc6
#define b0011000111 0xc7
#define b0011001000 0xc8
#define b0011001001 0xc9
#define b0011001010 0xca
#define b0011001011 0xcb
#define b0011001100 0xcc
#define b0011001101 0xcd
#define b0011001110 0xce
#define b0011001111 0xcf
#define b0011010000 0xd0
#define b0011010001 0xd1
#define b0011010010 0xd2
#define b0011010011 0xd3
#define b0011010100 0xd4
#define b0011010101 0xd5
#define b0011010110 0xd6
#define b0011010111 0xd7
#define b0011011000 0xd8
#define b0011011001 0xd9
#define b0011011010 0xda
#define b0011011011 0xdb
#define b0011011100 0xdc
#define b0011011101 0xdd
#define b0011011110 0xde
#define b0011011111 0xdf
#define b0011100000 0xe0
#define b0011100001 0xe1
#define b0011100010 0xe2
#define b0011100011 0xe3
#define b0011100100 0xe4
#define b0011100101 0xe5
#define b0011100110 0xe6
#define b0011100111 0xe7
#define b0011101000 0xe8
#define b0011101001 0xe9
#define b0011101010 0xea
#define b0011101011 0xeb
#define b0011101100 0xec
#define b0011101101 0xed
#define b0011101110 0xee
#define b0011101111 0xef
#define b0011110000 0xf0
#define b0011110001 0xf1
#define b0011110010 0xf2
#define b0011110011 0xf3
#define b0011110100 0xf4
#define b0011110101 0xf5
#define b0011110110 0xf6
#define b0011110111 0xf7
#define b0011111000 0xf8
#define b0011111001 0xf9
#define b0011111010 0xfa
#define b0011111011 0xfb
#define b0011111100 0xfc
#define b0011111101 0xfd
#define b0011111110 0xfe
#define b0011111111 0xff
#define b0100000000 0x100
#define b0100000001 0x101
#define b0100000010 0x102
#define b0100000011 0x103
#define b0100000100 0x104
#define b0100000101 0x105
#define b0100000110 0x106
#define b0100000111 0x107
#define b0100001000 0x108
#define b0100001001 0x109
#define b0100001010 0x10a
#define b0100001011 0x10b
#define b0100001100 0x10c
#define b0100001101 0x10d
#define b0100001110 0x10e
#define b0100001111 0x10f
#define b0100010000 0x110
#define b0100010001 0x111
#define b0100010010 0x112
#define b0100010011 0x113
#define b0100010100 0x114
#define b0100010101 0x115
#define b0100010110 0x116
#define b0100010111 0x117
#define b0100011000 0x118
#define b0100011001 0x119
#define b0100011010 0x11a
#define b0100011011 0x11b
#define b0100011100 0x11c
#define b0100011101 0x11d
#define b0100011110 0x11e
#define b0100011111 0x11f
#define b0100100000 0x120
#define b0100100001 0x121
#define b0100100010 0x122
#define b0100100011 0x123
#define b0100100100 0x124
#define b0100100101 0x125
#define b0100100110 0x126
#define b0100100111 0x127
#define b0100101000 0x128
#define b0100101001 0x129
#define b0100101010 0x12a
#define b0100101011 0x12b
#define b0100101100 0x12c
#define b0100101101 0x12d
#define b0100101110 0x12e
#define b0100101111 0x12f
#define b0100110000 0x130
#define b0100110001 0x131
#define b0100110010 0x132
#define b0100110011 0x133
#define b0100110100 0x134
#define b0100110101 0x135
#define b0100110110 0x136
#define b0100110111 0x137
#define b0100111000 0x138
#define b0100111001 0x139
#define b0100111010 0x13a
#define b0100111011 0x13b
#define b0100111100 0x13c
#define b0100111101 0x13d
#define b0100111110 0x13e
#define b0100111111 0x13f
#define b0101000000 0x140
#define b0101000001 0x141
#define b0101000010 0x142
#define b0101000011 0x143
#define b0101000100 0x144
#define b0101000101 0x145
#define b0101000110 0x146
#define b0101000111 0x147
#define b0101001000 0x148
#define b0101001001 0x149
#define b0101001010 0x14a
#define b0101001011 0x14b
#define b0101001100 0x14c
#define b0101001101 0x14d
#define b0101001110 0x14e
#define b0101001111 0x14f
#define b0101010000 0x150
#define b0101010001 0x151
#define b0101010010 0x152
#define b0101010011 0x153
#define b0101010100 0x154
#define b0101010101 0x155
#define b0101010110 0x156
#define b0101010111 0x157
#define b0101011000 0x158
#define b0101011001 0x159
#define b0101011010 0x15a
#define b0101011011 0x15b
#define b0101011100 0x15c
#define b0101011101 0x15d
#define b0101011110 0x15e
#define b0101011111 0x15f
#define b0101100000 0x160
#define b0101100001 0x161
#define b0101100010 0x162
#define b0101100011 0x163
#define b0101100100 0x164
#define b0101100101 0x165
#define b0101100110 0x166
#define b0101100111 0x167
#define b0101101000 0x168
#define b0101101001 0x169
#define b0101101010 0x16a
#define b0101101011 0x16b
#define b0101101100 0x16c
#define b0101101101 0x16d
#define b0101101110 0x16e
#define b0101101111 0x16f
#define b0101110000 0x170
#define b0101110001 0x171
#define b0101110010 0x172
#define b0101110011 0x173
#define b0101110100 0x174
#define b0101110101 0x175
#define b0101110110 0x176
#define b0101110111 0x177
#define b0101111000 0x178
#define b0101111001 0x179
#define b0101111010 0x17a
#define b0101111011 0x17b
#define b0101111100 0x17c
#define b0101111101 0x17d
#define b0101111110 0x17e
#define b0101111111 0x17f
#define b0110000000 0x180
#define b0110000001 0x181
#define b0110000010 0x182
#define b0110000011 0x183
#define b0110000100 0x184
#define b0110000101 0x185
#define b0110000110 0x186
#define b0110000111 0x187
#define b0110001000 0x188
#define b0110001001 0x189
#define b0110001010 0x18a
#define b0110001011 0x18b
#define b0110001100 0x18c
#define b0110001101 0x18d
#define b0110001110 0x18e
#define b0110001111 0x18f
#define b0110010000 0x190
#define b0110010001 0x191
#define b0110010010 0x192
#define b0110010011 0x193
#define b0110010100 0x194
#define b0110010101 0x195
#define b0110010110 0x196
#define b0110010111 0x197
#define b0110011000 0x198
#define b0110011001 0x199
#define b0110011010 0x19a
#define b0110011011 0x19b
#define b0110011100 0x19c
#define b0110011101 0x19d
#define b0110011110 0x19e
#define b0110011111 0x19f
#define b0110100000 0x1a0
#define b0110100001 0x1a1
#define b0110100010 0x1a2
#define b0110100011 0x1a3
#define b0110100100 0x1a4
#define b0110100101 0x1a5
#define b0110100110 0x1a6
#define b0110100111 0x1a7
#define b0110101000 0x1a8
#define b0110101001 0x1a9
#define b0110101010 0x1aa
#define b0110101011 0x1ab
#define b0110101100 0x1ac
#define b0110101101 0x1ad
#define b0110101110 0x1ae
#define b0110101111 0x1af
#define b0110110000 0x1b0
#define b0110110001 0x1b1
#define b0110110010 0x1b2
#define b0110110011 0x1b3
#define b0110110100 0x1b4
#define b0110110101 0x1b5
#define b0110110110 0x1b6
#define b0110110111 0x1b7
#define b0110111000 0x1b8
#define b0110111001 0x1b9
#define b0110111010 0x1ba
#define b0110111011 0x1bb
#define b0110111100 0x1bc
#define b0110111101 0x1bd
#define b0110111110 0x1be
#define b0110111111 0x1bf
#define b0111000000 0x1c0
#define b0111000001 0x1c1
#define b0111000010 0x1c2
#define b0111000011 0x1c3
#define b0111000100 0x1c4
#define b0111000101 0x1c5
#define b0111000110 0x1c6
#define b0111000111 0x1c7
#define b0111001000 0x1c8
#define b0111001001 0x1c9
#define b0111001010 0x1ca
#define b0111001011 0x1cb
#define b0111001100 0x1cc
#define b0111001101 0x1cd
#define b0111001110 0x1ce
#define b0111001111 0x1cf
#define b0111010000 0x1d0
#define b0111010001 0x1d1
#define b0111010010 0x1d2
#define b0111010011 0x1d3
#define b0111010100 0x1d4
#define b0111010101 0x1d5
#define b0111010110 0x1d6
#define b0111010111 0x1d7
#define b0111011000 0x1d8
#define b0111011001 0x1d9
#define b0111011010 0x1da
#define b0111011011 0x1db
#define b0111011100 0x1dc
#define b0111011101 0x1dd
#define b0111011110 0x1de
#define b0111011111 0x1df
#define b0111100000 0x1e0
#define b0111100001 0x1e1
#define b0111100010 0x1e2
#define b0111100011 0x1e3
#define b0111100100 0x1e4
#define b0111100101 0x1e5
#define b0111100110 0x1e6
#define b0111100111 0x1e7
#define b0111101000 0x1e8
#define b0111101001 0x1e9
#define b0111101010 0x1ea
#define b0111101011 0x1eb
#define b0111101100 0x1ec
#define b0111101101 0x1ed
#define b0111101110 0x1ee
#define b0111101111 0x1ef
#define b0111110000 0x1f0
#define b0111110001 0x1f1
#define b0111110010 0x1f2
#define b0111110011 0x1f3
#define b0111110100 0x1f4
#define b0111110101 0x1f5
#define b0111110110 0x1f6
#define b0111110111 0x1f7
#define b0111111000 0x1f8
#define b0111111001 0x1f9
#define b0111111010 0x1fa
#define b0111111011 0x1fb
#define b0111111100 0x1fc
#define b0111111101 0x1fd
#define b0111111110 0x1fe
#define b0111111111 0x1ff
#define b1000000000 0x200
#define b1000000001 0x201
#define b1000000010 0x202
#define b1000000011 0x203
#define b1000000100 0x204
#define b1000000101 0x205
#define b1000000110 0x206
#define b1000000111 0x207
#define b1000001000 0x208
#define b1000001001 0x209
#define b1000001010 0x20a
#define b1000001011 0x20b
#define b1000001100 0x20c
#define b1000001101 0x20d
#define b1000001110 0x20e
#define b1000001111 0x20f
#define b1000010000 0x210
#define b1000010001 0x211
#define b1000010010 0x212
#define b1000010011 0x213
#define b1000010100 0x214
#define b1000010101 0x215
#define b1000010110 0x216
#define b1000010111 0x217
#define b1000011000 0x218
#define b1000011001 0x219
#define b1000011010 0x21a
#define b1000011011 0x21b
#define b1000011100 0x21c
#define b1000011101 0x21d
#define b1000011110 0x21e
#define b1000011111 0x21f
#define b1000100000 0x220
#define b1000100001 0x221
#define b1000100010 0x222
#define b1000100011 0x223
#define b1000100100 0x224
#define b1000100101 0x225
#define b1000100110 0x226
#define b1000100111 0x227
#define b1000101000 0x228
#define b1000101001 0x229
#define b1000101010 0x22a
#define b1000101011 0x22b
#define b1000101100 0x22c
#define b1000101101 0x22d
#define b1000101110 0x22e
#define b1000101111 0x22f
#define b1000110000 0x230
#define b1000110001 0x231
#define b1000110010 0x232
#define b1000110011 0x233
#define b1000110100 0x234
#define b1000110101 0x235
#define b1000110110 0x236
#define b1000110111 0x237
#define b1000111000 0x238
#define b1000111001 0x239
#define b1000111010 0x23a
#define b1000111011 0x23b
#define b1000111100 0x23c
#define b1000111101 0x23d
#define b1000111110 0x23e
#define b1000111111 0x23f
#define b1001000000 0x240
#define b1001000001 0x241
#define b1001000010 0x242
#define b1001000011 0x243
#define b1001000100 0x244
#define b1001000101 0x245
#define b1001000110 0x246
#define b1001000111 0x247
#define b1001001000 0x248
#define b1001001001 0x249
#define b1001001010 0x24a
#define b1001001011 0x24b
#define b1001001100 0x24c
#define b1001001101 0x24d
#define b1001001110 0x24e
#define b1001001111 0x24f
#define b1001010000 0x250
#define b1001010001 0x251
#define b1001010010 0x252
#define b1001010011 0x253
#define b1001010100 0x254
#define b1001010101 0x255
#define b1001010110 0x256
#define b1001010111 0x257
#define b1001011000 0x258
#define b1001011001 0x259
#define b1001011010 0x25a
#define b1001011011 0x25b
#define b1001011100 0x25c
#define b1001011101 0x25d
#define b1001011110 0x25e
#define b1001011111 0x25f
#define b1001100000 0x260
#define b1001100001 0x261
#define b1001100010 0x262
#define b1001100011 0x263
#define b1001100100 0x264
#define b1001100101 0x265
#define b1001100110 0x266
#define b1001100111 0x267
#define b1001101000 0x268
#define b1001101001 0x269
#define b1001101010 0x26a
#define b1001101011 0x26b
#define b1001101100 0x26c
#define b1001101101 0x26d
#define b1001101110 0x26e
#define b1001101111 0x26f
#define b1001110000 0x270
#define b1001110001 0x271
#define b1001110010 0x272
#define b1001110011 0x273
#define b1001110100 0x274
#define b1001110101 0x275
#define b1001110110 0x276
#define b1001110111 0x277
#define b1001111000 0x278
#define b1001111001 0x279
#define b1001111010 0x27a
#define b1001111011 0x27b
#define b1001111100 0x27c
#define b1001111101 0x27d
#define b1001111110 0x27e
#define b1001111111 0x27f
#define b1010000000 0x280
#define b1010000001 0x281
#define b1010000010 0x282
#define b1010000011 0x283
#define b1010000100 0x284
#define b1010000101 0x285
#define b1010000110 0x286
#define b1010000111 0x287
#define b1010001000 0x288
#define b1010001001 0x289
#define b1010001010 0x28a
#define b1010001011 0x28b
#define b1010001100 0x28c
#define b1010001101 0x28d
#define b1010001110 0x28e
#define b1010001111 0x28f
#define b1010010000 0x290
#define b1010010001 0x291
#define b1010010010 0x292
#define b1010010011 0x293
#define b1010010100 0x294
#define b1010010101 0x295
#define b1010010110 0x296
#define b1010010111 0x297
#define b1010011000 0x298
#define b1010011001 0x299
#define b1010011010 0x29a
#define b1010011011 0x29b
#define b1010011100 0x29c
#define b1010011101 0x29d
#define b1010011110 0x29e
#define b1010011111 0x29f
#define b1010100000 0x2a0
#define b1010100001 0x2a1
#define b1010100010 0x2a2
#define b1010100011 0x2a3
#define b1010100100 0x2a4
#define b1010100101 0x2a5
#define b1010100110 0x2a6
#define b1010100111 0x2a7
#define b1010101000 0x2a8
#define b1010101001 0x2a9
#define b1010101010 0x2aa
#define b1010101011 0x2ab
#define b1010101100 0x2ac
#define b1010101101 0x2ad
#define b1010101110 0x2ae
#define b1010101111 0x2af
#define b1010110000 0x2b0
#define b1010110001 0x2b1
#define b1010110010 0x2b2
#define b1010110011 0x2b3
#define b1010110100 0x2b4
#define b1010110101 0x2b5
#define b1010110110 0x2b6
#define b1010110111 0x2b7
#define b1010111000 0x2b8
#define b1010111001 0x2b9
#define b1010111010 0x2ba
#define b1010111011 0x2bb
#define b1010111100 0x2bc
#define b1010111101 0x2bd
#define b1010111110 0x2be
#define b1010111111 0x2bf
#define b1011000000 0x2c0
#define b1011000001 0x2c1
#define b1011000010 0x2c2
#define b1011000011 0x2c3
#define b1011000100 0x2c4
#define b1011000101 0x2c5
#define b1011000110 0x2c6
#define b1011000111 0x2c7
#define b1011001000 0x2c8
#define b1011001001 0x2c9
#define b1011001010 0x2ca
#define b1011001011 0x2cb
#define b1011001100 0x2cc
#define b1011001101 0x2cd
#define b1011001110 0x2ce
#define b1011001111 0x2cf
#define b1011010000 0x2d0
#define b1011010001 0x2d1
#define b1011010010 0x2d2
#define b1011010011 0x2d3
#define b1011010100 0x2d4
#define b1011010101 0x2d5
#define b1011010110 0x2d6
#define b1011010111 0x2d7
#define b1011011000 0x2d8
#define b1011011001 0x2d9
#define b1011011010 0x2da
#define b1011011011 0x2db
#define b1011011100 0x2dc
#define b1011011101 0x2dd
#define b1011011110 0x2de
#define b1011011111 0x2df
#define b1011100000 0x2e0
#define b1011100001 0x2e1
#define b1011100010 0x2e2
#define b1011100011 0x2e3
#define b1011100100 0x2e4
#define b1011100101 0x2e5
#define b1011100110 0x2e6
#define b1011100111 0x2e7
#define b1011101000 0x2e8
#define b1011101001 0x2e9
#define b1011101010 0x2ea
#define b1011101011 0x2eb
#define b1011101100 0x2ec
#define b1011101101 0x2ed
#define b1011101110 0x2ee
#define b1011101111 0x2ef
#define b1011110000 0x2f0
#define b1011110001 0x2f1
#define b1011110010 0x2f2
#define b1011110011 0x2f3
#define b1011110100 0x2f4
#define b1011110101 0x2f5
#define b1011110110 0x2f6
#define b1011110111 0x2f7
#define b1011111000 0x2f8
#define b1011111001 0x2f9
#define b1011111010 0x2fa
#define b1011111011 0x2fb
#define b1011111100 0x2fc
#define b1011111101 0x2fd
#define b1011111110 0x2fe
#define b1011111111 0x2ff
#define b1100000000 0x300
#define b1100000001 0x301
#define b1100000010 0x302
#define b1100000011 0x303
#define b1100000100 0x304
#define b1100000101 0x305
#define b1100000110 0x306
#define b1100000111 0x307
#define b1100001000 0x308
#define b1100001001 0x309
#define b1100001010 0x30a
#define b1100001011 0x30b
#define b1100001100 0x30c
#define b1100001101 0x30d
#define b1100001110 0x30e
#define b1100001111 0x30f
#define b1100010000 0x310
#define b1100010001 0x311
#define b1100010010 0x312
#define b1100010011 0x313
#define b1100010100 0x314
#define b1100010101 0x315
#define b1100010110 0x316
#define b1100010111 0x317
#define b1100011000 0x318
#define b1100011001 0x319
#define b1100011010 0x31a
#define b1100011011 0x31b
#define b1100011100 0x31c
#define b1100011101 0x31d
#define b1100011110 0x31e
#define b1100011111 0x31f
#define b1100100000 0x320
#define b1100100001 0x321
#define b1100100010 0x322
#define b1100100011 0x323
#define b1100100100 0x324
#define b1100100101 0x325
#define b1100100110 0x326
#define b1100100111 0x327
#define b1100101000 0x328
#define b1100101001 0x329
#define b1100101010 0x32a
#define b1100101011 0x32b
#define b1100101100 0x32c
#define b1100101101 0x32d
#define b1100101110 0x32e
#define b1100101111 0x32f
#define b1100110000 0x330
#define b1100110001 0x331
#define b1100110010 0x332
#define b1100110011 0x333
#define b1100110100 0x334
#define b1100110101 0x335
#define b1100110110 0x336
#define b1100110111 0x337
#define b1100111000 0x338
#define b1100111001 0x339
#define b1100111010 0x33a
#define b1100111011 0x33b
#define b1100111100 0x33c
#define b1100111101 0x33d
#define b1100111110 0x33e
#define b1100111111 0x33f
#define b1101000000 0x340
#define b1101000001 0x341
#define b1101000010 0x342
#define b1101000011 0x343
#define b1101000100 0x344
#define b1101000101 0x345
#define b1101000110 0x346
#define b1101000111 0x347
#define b1101001000 0x348
#define b1101001001 0x349
#define b1101001010 0x34a
#define b1101001011 0x34b
#define b1101001100 0x34c
#define b1101001101 0x34d
#define b1101001110 0x34e
#define b1101001111 0x34f
#define b1101010000 0x350
#define b1101010001 0x351
#define b1101010010 0x352
#define b1101010011 0x353
#define b1101010100 0x354
#define b1101010101 0x355
#define b1101010110 0x356
#define b1101010111 0x357
#define b1101011000 0x358
#define b1101011001 0x359
#define b1101011010 0x35a
#define b1101011011 0x35b
#define b1101011100 0x35c
#define b1101011101 0x35d
#define b1101011110 0x35e
#define b1101011111 0x35f
#define b1101100000 0x360
#define b1101100001 0x361
#define b1101100010 0x362
#define b1101100011 0x363
#define b1101100100 0x364
#define b1101100101 0x365
#define b1101100110 0x366
#define b1101100111 0x367
#define b1101101000 0x368
#define b1101101001 0x369
#define b1101101010 0x36a
#define b1101101011 0x36b
#define b1101101100 0x36c
#define b1101101101 0x36d
#define b1101101110 0x36e
#define b1101101111 0x36f
#define b1101110000 0x370
#define b1101110001 0x371
#define b1101110010 0x372
#define b1101110011 0x373
#define b1101110100 0x374
#define b1101110101 0x375
#define b1101110110 0x376
#define b1101110111 0x377
#define b1101111000 0x378
#define b1101111001 0x379
#define b1101111010 0x37a
#define b1101111011 0x37b
#define b1101111100 0x37c
#define b1101111101 0x37d
#define b1101111110 0x37e
#define b1101111111 0x37f
#define b1110000000 0x380
#define b1110000001 0x381
#define b1110000010 0x382
#define b1110000011 0x383
#define b1110000100 0x384
#define b1110000101 0x385
#define b1110000110 0x386
#define b1110000111 0x387
#define b1110001000 0x388
#define b1110001001 0x389
#define b1110001010 0x38a
#define b1110001011 0x38b
#define b1110001100 0x38c
#define b1110001101 0x38d
#define b1110001110 0x38e
#define b1110001111 0x38f
#define b1110010000 0x390
#define b1110010001 0x391
#define b1110010010 0x392
#define b1110010011 0x393
#define b1110010100 0x394
#define b1110010101 0x395
#define b1110010110 0x396
#define b1110010111 0x397
#define b1110011000 0x398
#define b1110011001 0x399
#define b1110011010 0x39a
#define b1110011011 0x39b
#define b1110011100 0x39c
#define b1110011101 0x39d
#define b1110011110 0x39e
#define b1110011111 0x39f
#define b1110100000 0x3a0
#define b1110100001 0x3a1
#define b1110100010 0x3a2
#define b1110100011 0x3a3
#define b1110100100 0x3a4
#define b1110100101 0x3a5
#define b1110100110 0x3a6
#define b1110100111 0x3a7
#define b1110101000 0x3a8
#define b1110101001 0x3a9
#define b1110101010 0x3aa
#define b1110101011 0x3ab
#define b1110101100 0x3ac
#define b1110101101 0x3ad
#define b1110101110 0x3ae
#define b1110101111 0x3af
#define b1110110000 0x3b0
#define b1110110001 0x3b1
#define b1110110010 0x3b2
#define b1110110011 0x3b3
#define b1110110100 0x3b4
#define b1110110101 0x3b5
#define b1110110110 0x3b6
#define b1110110111 0x3b7
#define b1110111000 0x3b8
#define b1110111001 0x3b9
#define b1110111010 0x3ba
#define b1110111011 0x3bb
#define b1110111100 0x3bc
#define b1110111101 0x3bd
#define b1110111110 0x3be
#define b1110111111 0x3bf
#define b1111000000 0x3c0
#define b1111000001 0x3c1
#define b1111000010 0x3c2
#define b1111000011 0x3c3
#define b1111000100 0x3c4
#define b1111000101 0x3c5
#define b1111000110 0x3c6
#define b1111000111 0x3c7
#define b1111001000 0x3c8
#define b1111001001 0x3c9
#define b1111001010 0x3ca
#define b1111001011 0x3cb
#define b1111001100 0x3cc
#define b1111001101 0x3cd
#define b1111001110 0x3ce
#define b1111001111 0x3cf
#define b1111010000 0x3d0
#define b1111010001 0x3d1
#define b1111010010 0x3d2
#define b1111010011 0x3d3
#define b1111010100 0x3d4
#define b1111010101 0x3d5
#define b1111010110 0x3d6
#define b1111010111 0x3d7
#define b1111011000 0x3d8
#define b1111011001 0x3d9
#define b1111011010 0x3da
#define b1111011011 0x3db
#define b1111011100 0x3dc
#define b1111011101 0x3dd
#define b1111011110 0x3de
#define b1111011111 0x3df
#define b1111100000 0x3e0
#define b1111100001 0x3e1
#define b1111100010 0x3e2
#define b1111100011 0x3e3
#define b1111100100 0x3e4
#define b1111100101 0x3e5
#define b1111100110 0x3e6
#define b1111100111 0x3e7
#define b1111101000 0x3e8
#define b1111101001 0x3e9
#define b1111101010 0x3ea
#define b1111101011 0x3eb
#define b1111101100 0x3ec
#define b1111101101 0x3ed
#define b1111101110 0x3ee
#define b1111101111 0x3ef
#define b1111110000 0x3f0
#define b1111110001 0x3f1
#define b1111110010 0x3f2
#define b1111110011 0x3f3
#define b1111110100 0x3f4
#define b1111110101 0x3f5
#define b1111110110 0x3f6
#define b1111110111 0x3f7
#define b1111111000 0x3f8
#define b1111111001 0x3f9
#define b1111111010 0x3fa
#define b1111111011 0x3fb
#define b1111111100 0x3fc
#define b1111111101 0x3fd
#define b1111111110 0x3fe
#define b1111111111 0x3ff

#define b1001001001001 0x1249
#define b1001001001001001 0x9249
#define b1001001001001001001 0x49249

#endif

