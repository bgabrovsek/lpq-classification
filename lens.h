//  links.h
//  Created by Boštjan on 4/22/13.
//
//  Manipulation of knots in the lens space L(p,q) (slide moves, winding number, HOMLFY skein module in L(p,q), KBSM of L(p,q), etc.)
//
//  L(p,q), GCD(p,q) = 1, 0 <= q < p
//  Homotopy equivalence: q*q' = ±n^2 (mod p) n ∈ N
//  Homeomorphic equivalence: q = ±q'^±1 (mod p)

#ifndef Lpq_lens_h
#define Lpq_lens_h

#include <algorithm>
#include "global_vars.h"
#include "numbers.h"
#include "knot.h"
#include "link.h"
#include "common.h"
#include "reidemeister_moves.h"
#include "in_out.h"

#define CW 0 // clockwise
#define CCW 1 // counterclockwise
#define wise(a) (a ? "ccw " : "cw ")

#define DEBUG_LENS false

// HOMLFY skein module globals

#define ELEMENTARY_N 13 //15 // number of elementary generators: t_0, t_1, t_-1,...

int s_r[ELEMENTARY_N] = {0,2,1,4,3,6,5,8,7,10,9,12,11}; // sorting of elementary generators
int s_g[ELEMENTARY_N] = {0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6}; // sorted elementary generators
int torus_HSM_generator_table[MAX_HSM_GEN][ELEMENTARY_N]; // torus generators: 1, -1, 11, -1-1,....
Cterm * torus_HSM_generator_terms[MAX_HSM_GEN]; // torus generators as Cterm-s
int torus_HSM_generator_n; // number of torus generators

Cmultivariate * generator[N_HOMEO][MAX_HSM_GEN]; // the HSM of a torus generator in L(p,q)

// KBSM globals

Cterm *torus_KBSM_generator_terms[MAX_KBSM_GEN]; // torus generator terms (x^0, x^1,...)
int torus_KBSM_generator_n = MAX_KBSM_GEN; // number of KBSM generaotrs
Cmultivariate * generator_KBSM[N_HOMEO][MAX_HSM_GEN]; // the KBSM torus generator in L(p,q)

// HSM functions

bool is_lpq_generator(int p, int tgi) { // is torus_HSM_generator_table[tgi] generator of L(p,*)
    int max_g;
    for (max_g = ELEMENTARY_N-1; max_g > 0; max_g--) if (torus_HSM_generator_table[tgi][max_g]) break; // maximal genererator
    return ((0.5*p > -1.0*s_g[max_g]) && (0.5*p >= 1.0*s_g[max_g])); // -p/2 < s_g[max_g] <= p/2
}

int reverse(int g_index) { // get index of reverse generator
    bool b;
    int i;
    for (int g = 0; g < torus_HSM_generator_n; g++) {
        b = true;
        for (i=ELEMENTARY_N-1; i>0; i--) if (torus_HSM_generator_table[g][i] != torus_HSM_generator_table[g_index][s_r[i]]) {b = false; break;}
        if (b) return g;
    }
    throw 20;
}

bool is_not_generator_but_reverse_is(int p, int g_index) { // g_index is not a genererator of L(p,-), but the reverse of g_index is, obsolete
    int g_reverse = reverse(g_index);
    return ((!is_lpq_generator(p, g_index)) && (is_lpq_generator(p, g_reverse)));
}

int get_lpq_index(int p, int q) { // get index of L(p,q)
    for (int i=0; i < N_HOMEO; i++)
        if ((lens_homeo[i][0] == p) && (lens_homeo[i][1] == q)) return i;
    throw 20;
}

bool reverse_is_smaller(int g_index) { // true if the reverse of generator is previously in the list of generators
    for (int i=ELEMENTARY_N-1; i>0; i-=2) {
        if (torus_HSM_generator_table[g_index][i] > torus_HSM_generator_table[g_index][i-1]) return true;
        if (torus_HSM_generator_table[g_index][i] < torus_HSM_generator_table[g_index][i-1]) return false;
    }
    
    return false; // if equal, it is not smaller
}

// calculates the generators of the HSM of the solid torus,
// must be called before generating lens space generators

void generate_torus_HSM_generator_table() {
    int num[ELEMENTARY_N],i, sum;
    bool go = true;
    Cterm *term;
    
    memset(torus_HSM_generator_table,0,sizeof(torus_HSM_generator_table));
    torus_HSM_generator_n = 0;
    for (i=0;i<ELEMENTARY_N;i++) num[i]=0;
    
    while (1) {
        for (num[1]++,i=1;i<ELEMENTARY_N;i++) 
            if (num[i]*ABS(s_g[i]) > 6) {
                if (i>=ELEMENTARY_N-1) { go = false; break;}
                num[i] = 0;
                num[i+1]++;
            }
        if (!go) break;
        for (sum=0,i=1;i<ELEMENTARY_N;i++) sum += num[i]*ABS(s_g[i]);
        if (sum > 6) continue;
        for (i=1;i<ELEMENTARY_N;i++) torus_HSM_generator_table[torus_HSM_generator_n][i] = num[i];
        
        term = torus_HSM_generator_terms[torus_HSM_generator_n] = new Cterm(1);
        for (i=1;i<ELEMENTARY_N;i++) term->power[ reverse_variable[MID_VAR + s_g[i]] ] = torus_HSM_generator_table[torus_HSM_generator_n][i];
        
        torus_HSM_generator_n++;
    }

}

// is x^power a generator of KBSM(l)?
bool is_KBSM_generator(int l, int power) { return (power <= (lens_homeo[l][0]/2)); }

// string "L(p,q)" to index of L(p,q)
int stolpq(string s) {
    for (int l=0;l<N_HOMEO;l++) if (lpq_s[l] == s) return l;
    throw 20;
}

// KBSM functions

// calculate KBSM torus generators in terms of L(p,q) torus generators (polynomials obtained by Mathematica)
void calculate_KBSM_generators() { 

    int g, l;
    
    for (g=0;g<torus_KBSM_generator_n;g++)
        torus_KBSM_generator_terms[g] = new Cterm("x^" + itoa(g));
    
    for (l=0;l<N_HOMEO;l++) {
        for (g=0;g<torus_HSM_generator_n;g++)
            if (is_KBSM_generator(l, g))
                generator_KBSM[l][g] = new Cmultivariate("x^" + itoa(g));
    }
    
// Other relations calculated via given relations (double checked):
    
    generator_KBSM[stolpq("L(2,1)")][2] = mv("1 + a^-6 + a^-4 + a^-2");
    generator_KBSM[stolpq("L(2,1)")][3] = mv("2x + xa^-8 + xa^-4");
    generator_KBSM[stolpq("L(2,1)")][4] = mv("2 + a^-16 + a^-14 + a^-12 + a^-10 + a^-8 + 3a^-6 + 3a^-4 + 3a^-2");
    generator_KBSM[stolpq("L(2,1)")][5] = mv("5x + xa^-20 + xa^-16 + xa^-12 + 4xa^-8 + 4xa^-4");
    generator_KBSM[stolpq("L(2,1)")][6] = mv("5 + a^-30 + a^-28 + a^-26 + a^-24 + a^-22 + a^-20 + a^-18 + 5a^-16 + 5a^-14 + 5a^-12 + 5a^-10 + 5a^-8 + 9a^-6 + 9a^-4 + 9a^-2");
    generator_KBSM[stolpq("L(3,1)")][2] = mv("1 + a^-4 + xa^-1");
    generator_KBSM[stolpq("L(3,1)")][3] = mv("a^-7 + a^-3 + 2x + xa^-4");
    generator_KBSM[stolpq("L(3,1)")][4] = mv("2 + a^-8 + 3a^-4 + xa^-9 + xa^-5 + 3xa^-1");
    generator_KBSM[stolpq("L(3,1)")][5] = mv("a^-15 + a^-11 + 4a^-7 + 4a^-3 + 5x + xa^-12 + xa^-8 + 4xa^-4");
    generator_KBSM[stolpq("L(3,1)")][6] = mv("5 + a^-20 + a^-16 + a^-12 + 5a^-8 + 9a^-4 + xa^-17 + xa^-13 + 5xa^-9 + 5xa^-5 + 9xa^-1");
    generator_KBSM[stolpq("L(4,1)")][3] = mv("2x + xa^-4 + xa^-2");
    generator_KBSM[stolpq("L(4,1)")][4] = mv("-1 + a^-8 + 3x^2 + x^2a^-4");
    generator_KBSM[stolpq("L(4,1)")][5] = mv("5x + xa^-10 + xa^-8 + xa^-6 + 4xa^-4 + 4xa^-2");
    generator_KBSM[stolpq("L(4,1)")][6] = mv("-4 + 4a^-8 + 9x^2 + x^2a^-12 + x^2a^-8 + 5x^2a^-4");
    generator_KBSM[stolpq("L(5,1)")][3] = mv("-a^-5 - a^-1 + 2x + xa^-4 + x^2a^-1");
    generator_KBSM[stolpq("L(5,1)")][4] = mv("-1 - a^-4 + xa^-3 + 3x^2 + x^2a^-4");
    generator_KBSM[stolpq("L(5,1)")][5] = mv("-4a^-5 - 4a^-1 + 5x + xa^-8 + 4xa^-4 + x^2a^-5 + 4x^2a^-1");
    generator_KBSM[stolpq("L(5,1)")][6] = mv("-4 - a^-8 - 5a^-4 + xa^-11 + xa^-7 + 5xa^-3 + 9x^2 + x^2a^-8 + 5x^2a^-4");
    generator_KBSM[stolpq("L(5,2)")][3] = mv("a^-8 + a^-4 + 2x + xa^-4 - x^2a^-4");
    generator_KBSM[stolpq("L(5,2)")][4] = mv("-1 - a^-12 - xa^-8 + 3x^2 + x^2a^-8");
    generator_KBSM[stolpq("L(5,2)")][5] = mv("4a^-8 + 4a^-4 + 5x + xa^-12 + 4xa^-4 - x^2a^-12 - 4x^2a^-4");
    generator_KBSM[stolpq("L(5,2)")][6] = mv("-4 - a^-16 - 5a^-12 - xa^-20 - xa^-16 - 5xa^-8 + 9x^2 + x^2a^-16 + 5x^2a^-8");
    generator_KBSM[stolpq("L(6,1)")][4] = mv("-1 - a^-6 - a^-4 - a^-2 + 3x^2 + x^2a^-4 + x^2a^-2");
    generator_KBSM[stolpq("L(6,1)")][5] = mv("-3x - xa^-4 + 4x^3 + x^3a^-4");
    generator_KBSM[stolpq("L(6,1)")][6] = mv("-4 - a^-8 - 5a^-6 - 5a^-4 - 5a^-2 + 9x^2 + x^2a^-8 + x^2a^-6 + 5x^2a^-4 + 5x^2a^-2");
    generator_KBSM[stolpq("L(7,1)")][4] = mv("-1 - a^-4 - xa^-5 - 2xa^-1 + 3x^2 + x^2a^-4 + x^3a^-1");
    generator_KBSM[stolpq("L(7,1)")][5] = mv("-a^-7 - a^-3 - 3x - 2xa^-4 + x^2a^-3 + 4x^3 + x^3a^-4");
    generator_KBSM[stolpq("L(7,1)")][6] = mv("-4 - a^-8 - 5a^-4 - xa^-9 - 6xa^-5 - 10xa^-1 + 9x^2 + x^2a^-8 + 5x^2a^-4 + x^3a^-5 + 5x^3a^-1");
    generator_KBSM[stolpq("L(7,2)")][4] = mv("-1 - a^-4 + xa^-8 + 2xa^-4 + 3x^2 + x^2a^-4 - x^3a^-4");
    generator_KBSM[stolpq("L(7,2)")][5] = mv("a^-12 + a^-8 - 3x - xa^-12 - xa^-8 - x^2a^-8 + 4x^3 + x^3a^-8");
    generator_KBSM[stolpq("L(7,2)")][6] = mv("-4 - a^-16 - 5a^-4 + xa^-16 + xa^-12 + 5xa^-8 + 10xa^-4 + 9x^2 + x^2a^-12 + 5x^2a^-4 - x^3a^-12 - 5x^3a^-4");
    generator_KBSM[stolpq("L(8,1)")][5] = mv("-3x - xa^-6 - 2xa^-4 - 2xa^-2 + 4x^3 + x^3a^-4 + x^3a^-2");
    generator_KBSM[stolpq("L(8,1)")][6] = mv("1 - a^-8 - 6x^2 - 2x^2a^-4 + 5x^4 + x^4a^-4");
    generator_KBSM[stolpq("L(8,3)")][5] = mv("-3x - xa^-10 - 2xa^-8 - 2xa^-6 + 4x^3 + x^3a^-8 + x^3a^-6");
    generator_KBSM[stolpq("L(8,3)")][6] = mv("1 - 2a^-16 + a^-8 - 6x^2 + x^2a^-16 - 3x^2a^-8 + 5x^4 + x^4a^-8");
    generator_KBSM[stolpq("L(9,1)")][5] = mv("a^-5 + a^-1 - 3x - 2xa^-4 - x^2a^-5 - 3x^2a^-1 + 4x^3 + x^3a^-4 + x^4a^-1");
    generator_KBSM[stolpq("L(9,1)")][6] = mv("1 + a^-4 - xa^-7 - 2xa^-3 - 6x^2 - 3x^2a^-4 + x^3a^-3 + 5x^4 + x^4a^-4");
    generator_KBSM[stolpq("L(9,2)")][5] = mv("-a^-8 - a^-4 - 3x - 2xa^-4 + x^2a^-8 + 3x^2a^-4 + 4x^3 + x^3a^-4 - x^4a^-4");
    generator_KBSM[stolpq("L(9,2)")][6] = mv("1 + a^-12 + xa^-12 + 2xa^-8 - 6x^2 - x^2a^-12 - 2x^2a^-8 - x^3a^-8 + 5x^4 + x^4a^-8");
    generator_KBSM[stolpq("L(10,1)")][6] = mv("1 + a^-6 + a^-4 + a^-2 - 6x^2 - x^2a^-6 - 3x^2a^-4 - 3x^2a^-2 + 5x^4 + x^4a^-4 + x^4a^-2");
    generator_KBSM[stolpq("L(10,3)")][6] = mv("1 - a^-16 + a^-12 + a^-10 + a^-8 + a^-6 - 6x^2 - x^2a^-10 - 3x^2a^-8 - 3x^2a^-6 + 5x^4 + x^4a^-8 + x^4a^-6");
    generator_KBSM[stolpq("L(11,1)")][6] = mv("1 + a^-4 + 2xa^-5 + 3xa^-1 - 6x^2 - 3x^2a^-4 - x^3a^-5 - 4x^3a^-1 + 5x^4 + x^4a^-4 + x^5a^-1");
    generator_KBSM[stolpq("L(11,2)")][6] = mv("1 + a^-4 - 2xa^-8 - 3xa^-4 - 6x^2 - 3x^2a^-4 + x^3a^-8 + 4x^3a^-4 + 5x^4 + x^4a^-4 - x^5a^-4");
    generator_KBSM[stolpq("L(11,3)")][6] = mv("1 - a^-16 + a^-12 + a^-8 + xa^-15 + xa^-11 + 3xa^-3 - 6x^2 - 3x^2a^-8 - x^3a^-11 - 4x^3a^-3 + 5x^4 + x^4a^-8 + x^5a^-3");
}

// general lens space functions

// global variables used by calculation of winding strands of a torus knot

u_region w_path; // bit winding path
u_number w_arc[MAX+1]; // winding arcs sequence (in order from center)
s_small w_orient[MAX+1]; // orientation of arcs (CW/CCW) (in order from center)
u_number w_sort[MAX+1]; // sorting of w_arc (in order from appearance in gw)
u_number w_n; // number of winding arcs

// get winding strands of a knot
void get_winding_strands(Cknot *K) {
    
    u_region s_path, c_path;
    u_number reg = K->reg_0;
    u_number arc;
    s_small orient, i, pos;
    
    // get winding paths
    s_path = shortestPath(K);
        
    c_path = s_path;
    w_n = 0;

    do {
        arc = firstBitSet(K->regions[reg] & c_path);
        c_path ^= ((u_region)1<< arc);
        
        orient = K->regionSide(reg, arc);
        
        w_arc[w_n] = arc;
        w_orient[w_n] = orient;
        w_n++;
        reg = K->adjacentRegion(reg, arc);
    } while ( reg != K->reg_1);

    w_path = s_path;
    
    // sort winding paths
    pos = 0;
    while (s_path) {
        arc = firstBitSet(s_path);
        s_path ^= ((u_region)1<<arc);
        for (i=0; i<w_n; i++) if (w_arc[i] == arc) { w_sort[pos++] = i; break; };
    }
}

// get the winding number of K
int winding_number(Cknot *K) {
    int w = 0;

    if (K->n == 0) {
        if (K->reg_0 == K->reg_1) return 0; else
        if (K->reg_0 < K->reg_1) return 1; else
            return -1;
    }
    get_winding_strands(K);
    for (int i=0; i<w_n;i++)
        w += (int) (w_orient[i] == CCW ? 1 : -1);
    return w;
}

// calculate the winding number of K after an slide move would be performed
int winding_number_after_omega_4(int l, Cknot *K, s_small extra) {
    int w = winding_number(K);
    return w - (w_orient[w_n-1] == CCW ? 1 : -1)*lens_homeo[l][0];
}

// include the over-strands for the lens-space-move? if so, how much to shift each part?

// generate a (P,Q) torus knot
// default arguments: generate standard (P,Q) knot
// otherwise leave space and crossings for strands appearing in the arraching of (P,Q) to another not (for the slide move) 
Cknot * get_torus_knot(int P, int Q, bool over_strands = false, int knot_n = 0, int winding_n = 0) {
    if ((P-1)*Q  > MAX-1) return NULL; // knot too bigs
    int i, j, x, M = (P-1)*Q, pos = 0;
    u_region b_0 = (u_region)0, b_1 = (u_region)0;
    Cknot *K = new Cknot();
    K->n = (P-1)*Q + winding_n*Q;
    K->n2 = K->n << 1;
    x = 0;
    for (j = 0; j < Q; j++) {
        // diagonal strands
        for (i = 1; i < P; i++) {
            K->gw_[pos++] = -(x+1 + knot_n + winding_n * 2 * Q);
            x = MOD(x+1,M);
        }
        b_0 |= (u_region)1<<(pos-1);

        // over knot strands        
        for (i=0; i<winding_n; i++) // going up
            K->gw_[pos++] = -(i + knot_n + 1 + MOD((j*MOD(P,Q)),Q)*winding_n*2);
        for (i=0; i<winding_n; i++) // going down
            K->gw_[pos++] = (i + knot_n + winding_n + 1 + MOD((j*MOD(P,Q)),Q)*winding_n*2);
        
        // wave interlace strands
        for (i = 1; i < P; i++) {
            x = MOD(x + P - 2, M);
            K->gw_[pos++] = (x+1 +  + knot_n + winding_n * 2 * Q);
        }
        x = MOD(x + P - 1, M);
        b_1 |= (u_region)1<<(pos-1);
    }
    
    if (over_strands) return K;
    
    K->signs = 0;
    K->ID = P * 100 + Q;
    K->generateRegions();
    K->reg_0 = K->bin2reg(b_0);
    K->reg_1 = K->bin2reg(b_1);
    return K;
}

Cknot * get_elementary_generator(int i) {
    Cknot *Q;
    if (i == 0) { Q = new Cknot(); Q->n_reg = 2; Q->reg_0 = Q->reg_1 = 0; return Q; }
    if (ABS(i) == 1) { Q = new Cknot(); Q->n_reg = 2; Q->reg_0 = (i>0?0:1); Q->reg_1 = (i>0?1:0); return Q;}
    Q = get_torus_knot(ABS(i),1);
    if (i < 0) Q->reverse();
    return Q;
}

// attach knot B to knot A via crossing x, without signs
Cknot * attach(Cknot *A, Cknot *B, s_letter x) {
    Cknot *Q = new Cknot();
    u_number i, j, pos = 0, xpos;
    s_letter l;
    
    Q->n = (A->n) + (B->n) - 1;
    Q->n2 = Q->n*2;
    
    for (i=0; i<A->n2; i++) {
        l = A->gw(i);
        if (ABS(l) != x) {
            Q->gw_[pos++] = (ABS(l) > x ? ABS_PRED(l) : l); // insert next gauss letter from A
        } else {
            // insert knot B
            xpos = B->gw_pos(x);
            for (j=1; j<B->n2; j++) {
                l = B->gw(xpos+j);
                Q->gw_[pos++] = (ABS(l) > x ? ABS_PRED(l) : l);
            }
        }
    }
    return Q;
}



void insert_crossings_on_arc(Clink *L, u_number arc, s_letter x1, s_letter x2) { // inserts letters x1x2 on arc
    
    int i, arc_0;
    u8 c;
    
    for (arc_0 = 0, c=0;c<L->n_components;arc_0 += L->n_arcs[c++]) // find component
        if (arc < arc_0 + L->n_arcs[c]) break;
    
    arc -= arc_0;
    for (i=L->n_arcs[c]-1;i>arc;i--) L->gw_[c][i+2] = L->gw_[c][i];
    
    L->gw_[c][arc+1] = x1;
    L->gw_[c][arc+2] = x2;
    L->n++;
    L->n_arcs[c]+=2;
}

// perform slide move (extra: the phi slide move, extra = 0 performs the theta move, extra = 1 performs the psi move (in Hoste: KBSM of L(p,q)
// lpqi: select lens space
Cknot * Omega_4(int lpqi, Cknot *K_, s_small extra) {
    
    int l_p = lens_homeo[lpqi][0];
    int l_q = lens_homeo[lpqi][1];

    Cknot *T, *K, *sum; // torus knot, knot K, their sum (i.e. attaching T to K)
    int i, q, a, a_, last_strand, curr_strand, torus_knot_orientation;
    u_number k_n, curr_sign, attach_crossing;
    
    if (K_->n > 0)  K = K_->deepCopy();
    else {
        if (K_->reg_0 == K_->reg_1) throw 20; // affine
        K = new Cknot( (K_->reg_0 == 0 ? "-1 1 + *0 &01 1" : "1 -1 + *0 &01 1") ); // non-affine trivial knot
    }
    
    k_n = K->n;
    
    get_winding_strands(K);
    
    // insert crossings with the knot
    for (i = w_n-1; i>=0; i--) {  // go through each winding strand
        a = k_n + w_n - w_sort[i]; //w_n - w_sort[i] + k_n; // crossing pair of new crossings the winding strand
        a_ = k_n + w_n + 1 + w_sort[i]; //2*w_n+1 + 2*k_n - (w_n - w_sort[i] + k_n) ;
        
        if (w_orient[w_sort[i]] == CW) { // winding strand goes counter-clockwise
            for (q = 0; q < l_q; q++) // repeat q-times
                insertCrossingToGaussWord(K, w_arc[w_sort[i]]+1, w_arc[w_sort[i]]+1, - (a_ + q*w_n*2), (a + q*w_n*2));
        } else { // clockwise
            for (q = l_q-1; q >= 0; q--)
                insertCrossingToGaussWord(K, w_arc[w_sort[i]]+1, w_arc[w_sort[i]]+1, (a + q*w_n*2), -(a_ + q*w_n*2));
        }
        
    }
    
    // is last "attaching" strand cw or ccw
    last_strand = w_orient[w_n-1];
    
    T = get_torus_knot(l_p, l_q, true, k_n, w_n);

    // is the torus knot direction reversed?
    if ((last_strand == CCW) ^ (extra)) T->reverse();
    
    attach_crossing = k_n+w_n*l_q*2;
    
    sum = attach(K,T,attach_crossing);

    // make signs
    sum->signs |= (K->signs & fbf_sign(k_n+1)); // copy signs of original knot crossings

    for (i=w_n-1; i >= 0; i--) { // from outer strand to inner
        a = k_n + w_n - i; //w_n - i + k_n; // crossing pair of new crossings the winding strand
        a_ = k_n + w_n + 1 + i; //w_n+1 + k_n + i) ;
        curr_strand = w_orient[i];
        curr_sign = 1 ^ curr_strand ^ last_strand ^ extra;
        for (q=0; q < l_q; q++)
            sum->signs |= ((u_sign)curr_sign << (a+q*w_n*2)) | ((u_sign)curr_sign << (a_+q*w_n*2));
    }
    
    SETBIT(sum->signs,(u_sign)attach_crossing, (u_sign)0); // we deleted the last attaching crossing
    
    // all other bits are positive (=0)
    
    delete K;
    delete T;
    
    sum->generateRegions();
    //sum->check_ok();
    
    torus_knot_orientation = last_strand ^ extra ^ 1;
    
    // get new punctured regions
    if (torus_knot_orientation == CCW) {
        sum->reg_0 = sum->regionGet(sum->gw_pos(-(int)(k_n+w_n)), LEFT);
        sum->reg_1 = sum->regionGet(sum->gw_pos((int)attach_crossing), RIGHT);
    } else {
        sum->reg_0 = sum->regionGet(sum->gw_pos((int)k_n+w_n+1), RIGHT);
        sum->reg_1 = sum->regionGet(sum->gw_pos(-(int)attach_crossing), LEFT);
    }
    
    return sum;
}

// delete letter i on component c of L
void delete_gauss_letter(Clink *L, u8 c, u_number i) {
    for (int j=i+1; j<L->n_arcs[c];j++) L->gw_[c][j-1] = L->gw_[c][j];
    L->n_arcs[c]--;
}

void delete_zero_letters(Clink *L) {
    int i;
    for (u8 c=0;c<L->n_components;c++)
        for (i=L->n_arcs[c]-1;i>=0;i--)
            if (L->gw_[c][i] == 0) 
                delete_gauss_letter(L, c, i); // shift left
}

// find letter x in L (return component c and position p)
void find_gauss_letter(Clink *L, s_letter x, u8 *c, u_number *p) {
    for (*c = 0; *c < L->n_components; (*c)++)
        for (*p=0;*p<L->n_arcs[*c];(*p)++)
            if (L->gw_[*c][*p] == x) return;
}

u_number wl_arc[MAX+1]; // winding arcs sequence (in order from center)
s_small wl_orient[MAX+1]; // orientation of arcs (CW/CCW) (in order from center)
u_number wl_sort[MAX+1]; // sorting of w_arc (in order from appearance in gw)
u_number wl_n; // number of winding arcs

// perform the slide move on a link (similar than the function of performing a slide move on a knot)
Clink * Omega_4(int lpqi, Clink *L) {

    int l_p = lens_homeo[lpqi][0], l_q = lens_homeo[lpqi][1];
    Cknot *T;
    int i,a, a_, l_n, q;
    u_number attach_crossing;
    
    l_n = L->n;
    
    // 1st stage: adding crossings to components
    
    for (i = wl_n-1; i>=0; i--) {  // go through each winding strand
        
        a = l_n + wl_n - wl_sort[i]; //w_n - w_sort[i] + k_n; // crossing pair of new crossings the winding strand
        a_ = l_n + wl_n + 1 + wl_sort[i]; //2*w_n+1 + 2*k_n - (w_n - w_sort[i] + k_n) ;

        if (wl_orient[wl_sort[i]] == -1/*CW*/) { // winding strand goes counter-clockwise
            for (q = 0; q < l_q; q++) // repeat q-times
                insert_crossings_on_arc(L, wl_arc[wl_sort[i]], -(a_ + q*wl_n*2), (a + q*wl_n*2));
        } else { // clockwise
            for (q = l_q-1; q >= 0; q--)
                insert_crossings_on_arc(L, wl_arc[wl_sort[i]], (a + q*wl_n*2), -(a_ + q*wl_n*2));
        }
    }

    delete_zero_letters(L);
    
    // 2nd stage: generate torus knot
    
    T = get_torus_knot(l_p, l_q, true, l_n, wl_n);
    
    // 3rd stage: attach torus knot to link
    
    attach_crossing = l_n+wl_n*l_q*2;
    
    u_number pos, pos_;
    u8 c = L->n_components-1;
    
    if (wl_orient[wl_n-1] == 1) T->reverse(); // attach to last component, adjust orientation
    
    find_gauss_letter(L, -attach_crossing, &c, &pos);
    pos_ = T->gw_pos(attach_crossing);
    
    // shift gauss word
    for (i = L->n_arcs[c]-1; i > pos; i--)
        L->gw_[c][ i + T->n*2-2 ] = L->gw_[c][ i ];
    
    L->n += T->n;
    L->n--;
    L->n_arcs[c] += T->n*2;
    L->n_arcs[c] -= 2;
    
    // insert T
    
    s_letter l;
    for (i=0; i<T->n*2-1; i++) {
        l = T->gw_[ MOD( 1+ i + pos_ ,T->n*2) ];
        L->gw_[c][ pos + i] = (ABS(l) > attach_crossing ? ABS_PRED(l) : l);
    }
    
    // signs
    u_number curr_sign;
    L->signs = 0;
    
    for (i=wl_n-1; i >= 0; i--) { // from outer strand to inner
        a = l_n + wl_n - i; //w_n - i + k_n; // crossing pair of new crossings the winding strand
        a_ = l_n + wl_n + 1 + i; //w_n+1 + k_n + i) ;
        curr_sign = 1 ^ SIGN2B(wl_orient[i]) ^ SIGN2B(wl_orient[wl_n-1]) ^ 0;
        for (q=0; q < l_q; q++)
            L->signs |= ((u_sign)curr_sign << (a+q*wl_n*2)) | ((u_sign)curr_sign << (a_+q*wl_n*2));
    }
    
    SETBIT(L->signs,(u_sign)attach_crossing, (u_sign)0); // we deleted the last attaching crossing
    
    // regions
    L->generate_regions();
    
    if ( wl_orient[wl_n-1] == -1) {// torus knot CCW
        L->reg_0 = L->get_region(L->arc_after(-(int)(l_n+wl_n)), LEFT);
        L->reg_1 = L->get_region(L->arc_after((int)attach_crossing), RIGHT);
    } else {
        L->reg_0 = L->get_region(L->arc_after((int)l_n+wl_n+1), RIGHT);
        L->reg_1 = L->get_region(L->arc_after(-(int)attach_crossing), LEFT);
    }

    return L;
}

// return a L(p,q) generator of HSM as a link class
Clink * generator_link(int gi) {
    
    int i, i_, j, k, s, as;
    Cknot *K;
    Clink *L;
    u_number arc_0;

    L = new Clink();
    
    wl_n = 0; // winding strands
    arc_0 = 0; // 1st arc of component
    
    for (i_=1; i_<ELEMENTARY_N; i_++) { // loop through all elementary generators
       
        i = i_; // obsolete
        for (j=0; j<torus_HSM_generator_table[gi][i];j++) {
            s = s_g[i];
            as = ABS(s_g[i]);
            if (as == 1) {
                
                wl_arc[wl_n] = arc_0;
                wl_orient[wl_n] = (s > 0 ? 1 : -1); // +-1
                wl_n++;
                K = new Cknot();
                
            } else {
                
                K = get_elementary_generator(s_g[i]);
                if (s > 0) { // positive orientation
                    for (k=0; k<as; k++) {
                        wl_arc[wl_n] = k + arc_0 + as - 2;
                        wl_orient[wl_n] = 1;
                        wl_n++;
                    }
                } else { // negative orientation
                    for (k=as-2; k>=0; k--) {
                        wl_arc[wl_n] = k + arc_0;
                        wl_orient[wl_n] = -1;
                        wl_n++;
                    }
                    wl_arc[wl_n] = 2*(as-1)-1 + arc_0;
                    wl_orient[wl_n] = -1;
                    wl_n++;
                }
            }

            L->add_disjunct_component(K);
            arc_0 += (as != 1 ? K->n*2 : 1);
            delete K;
        }
    }

    // sort winding paths
    u_region sl_path = 0;
    u_number arc;
    for (i=0;i<wl_n;i++) sl_path |= (u_region)1 << wl_arc[i];
    int pos = 0;
    while (sl_path) {
        arc = firstBitSet(sl_path);
        sl_path ^= ((u_region)1<<arc);
        for (i=0; i<wl_n; i++) if (wl_arc[i] == arc) { wl_sort[pos++] = i; break; };
    }

    return L;
}

#define CALC_HOMFLY true
#define SHOW_CROSSING_NUMBERS true

void load_generators(string);
void save_generators_to_file(string, bool);
void print_generators();

// t to index of t
int get_generator_index(Cterm *t) {
    for (int g=0; g<torus_HSM_generator_n; g++ )
        if (t->eq_term(torus_HSM_generator_terms[g])) return g;
    throw 20;
}

void print_generators() {
    int l,g;
    for (l=0;l< N_HOMEO;l++) {
        cout << "L(" << lens_homeo[l][0] << "," << lens_homeo[l][1] << ")" << endl << endl;
        for (g=0; g<torus_HSM_generator_n; g++ ) {
            cout << torus_HSM_generator_terms[g] << " =" << endl;
            if (generator[l][g] != NULL) cout << generator[l][g] << endl;
            else cout << "NULL" << endl;
            cout << endl;
        }
    }
    
}

void print_generators_binary() {
    int i = 0;
    for (int l=0;l< N_HOMEO;l++) {
        cout << "L(" << lens_homeo[l][0] << "," << lens_homeo[l][1] << ") ";
        for (int g=0; g<torus_HSM_generator_n; g++ )
            if ((i++%2) && (generator[l][g] != NULL)) cout << "*"; else cout << "-";
        cout << endl;
    }
}

// polynomial manipulation

// converts the HOMLFY SM of a solid torus knot to the HOMLFY SM of lens space l

Cmultivariate *HSM_lpq(int l, Cmultivariate * poly_) { 
    bool change;
    Cterm * t;
    int gi, p, q;
    Cmultivariate * poly = poly_->deepCopy();
    
    p = lens_homeo[l][0];
    q = lens_homeo[l][1];
    
    do {
        change = false;
        
        for (list<Cterm*>::iterator it=poly->terms.begin(); it!=poly->terms.end(); it++) { // loop terms
        
            if (!(*it)->has_upper()) continue;
        
            t = (*it)->upper();
            gi = get_generator_index(t);
            
            if (!is_lpq_generator(p, gi)) {
                poly->replace_exact_subfactor_with_poly(t, generator[l][gi]);
                change = true;
                break;
            }
            
        }
    } while (change);
    poly->simplify();
    return poly;
}

// converts the KBSM SM of a solid torus knot to the KBSM SM of lens space l

Cmultivariate *KBSM_lpq(int l, Cmultivariate * poly_) { 
    bool change;
    Cterm * t;
    int gi, p, q;
    Cmultivariate * poly = poly_->deepCopy();
    
    p = lens_homeo[l][0];
    q = lens_homeo[l][1];
    
    do {
        change = false;
        
        for (list<Cterm*>::iterator it=poly->terms.begin(); it!=poly->terms.end(); it++) { // loop through terms
            
            if (!(*it)->has_x()) continue; // if no x-term, skip 
            t = (*it)->x(); // get power of x in term *it
          
            gi = t->power[reverse_variable['x']];

            if (!is_KBSM_generator(l, gi)) { // if generator in in the basis of KBSM(L(p,q)), do a replacement
                poly->replace_exact_subfactor_with_poly(t, generator_KBSM[l][gi]);
                change = true;
                break;
            }
            
        }
    } while (change);
    poly->simplify();
    return poly;
}

Cmultivariate * HSM_lpq(int l, Cknot *K) { // returns the HSM(K) in l=L(p,q)
    Cmultivariate * poly = HOMFLY(K);
    Cmultivariate * b = HSM_lpq(l,poly);
    delete poly;
    return b;
}

Cknot * get_knot_by_id(int);

// calculates the HSMs of the knots in the set KNOTS, stores the results in HOMLFYS_lpq (and the reverse knot in HOMFLYS_lpq_reverse) 
void calculate_HSMs_lpq(t_knot_set * KNOTS) {
    
    int l, ID;
    Cmultivariate * poly;
    
    for (l=0;l< N_HOMEO;l++) {
        
        if (lens_homeo[l][1] <= 2) {
            
            HOMFLYS_lpq[l].clear();
            HOMFLYS_lpq_reverse[l].clear();
            
            // include the unknot in the list
            HOMFLYS_lpq[l].push_back(new Cmultivariate("- vz^-1 + v^-1z^-1")); // unknot
            HOMFLYS_lpq_reverse[l].push_back(new Cmultivariate("- vz^-1 + v^-1z^-1")); // unknot
            
            // loop throuh all knots and insert the HOMFLYs
            for (ID=1; ID <= number_of_knots; ID++) {
                poly = HSM_lpq(l, HOMFLYS[ID]);
                HOMFLYS_lpq[l].push_back(poly->deepCopy());
                delete poly;
                poly =HSM_lpq(l, HOMFLYS_reverse[ID]);
                HOMFLYS_lpq_reverse[l].push_back(poly->deepCopy());
                delete poly;
            }

            
        } else {
            // we do not know the HSM of L(p,q)
            HOMFLYS_lpq[l].push_back(new Cmultivariate("1")); // unknot
            HOMFLYS_lpq_reverse[l].push_back(new Cmultivariate("1")); // unknot
            
            for (ID=1; ID <= number_of_knots; ID++) {
                HOMFLYS_lpq[l].push_back(new Cmultivariate("1")); // unknot
                HOMFLYS_lpq_reverse[l].push_back(new Cmultivariate("1")); // unknot
            }
        }
        
    }
    
    
}

// calculates the KBSMs of the knots in the set KNOTS, stores the results in KBSM_lpq
void calculate_KBSMs_lpq(t_knot_set * KNOTS) {
    
    int l, ID;
    Cmultivariate * poly;
 
    for (l=0;l< N_HOMEO;l++) {
        
        // inlcude the unknot
        KBSMS_lpq[l].clear();
        KBSMS_lpq[l].push_back(new Cmultivariate("-a^2 -a^-2"));
        
        for (ID=1; ID <= number_of_knots; ID++) {
            poly = KBSM_lpq(l, KBSMS[ID]); // get the KBSM
            normalize_framing(poly); // normalize the poly, since we do not distinguish between framings
            KBSMS_lpq[l].push_back(poly->deepCopy());
            delete poly;
        }
    }
}

// estimate the number of crossings after a lpq move
int estimate_number_of_crossings(int l, Cknot * K) {
    get_winding_strands(K);
    return K->n + 2*lens_homeo[l][1]*w_n + (lens_homeo[l][0]-1)*lens_homeo[l][1] - 1;
}

// completely expresses generators of the torus in terms of the lens space, once they are partly expressed after one recursive step
void process_generators_iteratively() { 
    
    int l, g, g_, p, q;
    Cterm * t;
    Cmultivariate *poly;
    bool change;
    
    for (l=0;l< N_HOMEO;l++) if (lens_homeo[l][1] <= 2) { // select lens space
        
        // 1st stage, express "reverse is smaller" generators in terms of calculated reversed generators
  
        p = lens_homeo[l][0];
        q = lens_homeo[l][1];
        
        for (g=0; g<torus_HSM_generator_n; g++ ) {
            
            if (is_lpq_generator(p, g)) continue;
            
            poly = generator[l][g];
           
            do { // repeat expressing the terms in lower terms until the polynomial does not change anymore
            
                change = false;
                for (list<Cterm*>::iterator it=poly->terms.begin(); it!=poly->terms.end(); it++) { // go through all terms
                    if (!(*it)->has_upper()) continue;
                    t = (*it)->upper();
                    
                    g_ = get_generator_index(t);
                
                    if (!is_lpq_generator(p, g_)) { // found a term that is not a generator
                        poly->replace_exact_subfactor_with_poly(t, generator[l][g_]);
                        change = true;
                        break;
                    }
                }
            } while (change); 
            poly->simplify();
        
        } // g
    
    } // l
}

// free memory of generators
void free_generators() {
    int l, g;
    for (g=0; g<torus_HSM_generator_n; g++ ) delete torus_HSM_generator_terms[g];
    for (l=0;l< N_HOMEO;l++)
        for (g=0; g<torus_HSM_generator_n; g++ )
            if (generator[l][g] != NULL) delete generator[l][g];
}

// calculates how torus generators are expressed in lower terms of lens space generators
// first performs one step of the recursion process and saves results in GENERATORS_PATH_INTERMEDIATE on each step
// then the rest is done by process_generators_iteratively(), saves results in GENERATORS_PATH
void lens_space_generators_in_terms_of_torus_generators() {
    
    int l, g;
    Cmultivariate *poly;
    Clink *L;
    
    for (l=0;l<N_HOMEO;l++)
        for (g=0; g<torus_HSM_generator_n; g++) generator[l][g] = NULL;
    
    
    for (l=0;l< N_HOMEO;l++) { // select lens space
        
        cout << lpq_s[l] << ": " << flush;
    
        for (g = 0; g < torus_HSM_generator_n; g++ ) { // select torus generator
            
            if (!is_lpq_generator(lens_homeo[l][0], g)) { // is it a generator of L(p,q)?
                if (lens_homeo[l][1] <= 2) { 
                
                    L = generator_link(g);
                    Omega_4(l, L);
                
                    poly = HOMFLY(L);
                    poly->simplify(); // needed?
                    
                    cout << torus_HSM_generator_terms[g] << ", " <<flush;
                    
                    generator[l][g] = poly;
                
                    // save generators on each step
                    save_generators_to_file(GENERATORS_PATH_INTERMEDIATE, false);
                } else
                if (is_lpq_generator(lens_homeo[l][0], g)) {
                    poly = new Cmultivariate;
                    poly->zerofy();
                    *poly += torus_HSM_generator_terms[g]->copy();
                    generator[l][g] = poly;
                }
            } // is lpq generator
        } // g
        cout << endl;
    } // l
    
    process_generators_iteratively();
    save_generators_to_file(GENERATORS_PATH, false);
    free_generators();
 }

// classification of lens spaces

#define MAX_LPQ 100

class clpq { public: int p, q; };

clpq lpq_homotopy[MAX_LPQ];
int n_homotopy = 0;
clpq lpq_homeo[MAX_LPQ];
int n_homeo = 0;

// is L(p1, q1) homotopic to L(p2,q2) ?
bool homotopic_eq(int p1, int q1, int p2, int q2) {
    if (p1 != p2) return false;
    for (int i = 0; i <= p1*p2; i++) // p1*p2 enough?
        if ((MOD_EQ(q1*q2, i*i, p1)) || (MOD_EQ(q1*q2, -i*i, p1))) return true;
    return false;
}

// is L(p1, q1) homeomorphic to L(p2,q2) ?
bool homeo_eq(int p1, int q1, int p2, int q2) {
    if (p1 != p2) return false;
  //  cout << "test L(" << p1 << "," << q1 << ") = L(" << p2 << "," << q2 << ")" << endl;
    if ((MOD_EQ(q1, q2, p1)) || (MOD_EQ(q1, -q2, p1))) return true;
    int q2_inv = mod_inverse(q2, p1);
    if (q2_inv == 0) return false;
    if ((MOD_EQ(q1, q2_inv, p1)) || (MOD_EQ(q1, -q2_inv, p1))) return true;
    return false;
}

// generate all homotopy equivalent L(p,q)'s up to  0 <= q < p <= n
void generate_homotopy_equvalent_lpqs(int n) {
    int p, q, i;
    bool b;
    for (p = 1; p <= n; p++)
        for (q = 1; q < p; q++)
            if (GCD(p,q) == 1) { // (p,q) candidates
                for (b=true,i=0;i<n_homotopy;i++) b &= !homotopic_eq(p, q, lpq_homotopy[i].p, lpq_homotopy[i].q);
                if (b) {
                    cout << "{" <<p << "," <<q << "}, ";
                    lpq_homotopy[n_homotopy].p = p;
                    lpq_homotopy[n_homotopy].q = q;
                    n_homotopy++;
                }
            }
}

// generate all homeo equivalent L(p,q)'s up to  0 <= q < p <= n
void generate_homeo_equvalent_lpqs(int n) { // 0 <= q < p <= n
    int p, q, i;
    bool b;
    for (p = 1; p <= n; p++)
        for (q = 1; q < p; q++)
            if (GCD(p,q) == 1) { // (p,q) candidates
                for (b=false,i=0;i<n_homeo;i++) b |= homeo_eq(p, q, lpq_homeo[i].p, lpq_homeo[i].q);
                if (!b) {
                    cout << "{" <<p << "," <<q << "}, ";
                    lpq_homeo[n_homeo].p = p;
                    lpq_homeo[n_homeo].q = q;
                    n_homeo++;
                }
            }
}

#endif
