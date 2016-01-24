//  links.h
//  Created by Boštjan on 4/22/13.
//
//  A link class that represents links in the solid torus, used for calculating the HOMLFY skein module or manipulating links

#ifndef Lpq_links_h
#define Lpq_links_h

#include <list>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <climits>

#include "common.h"
#include "numbers.h"
#include "knot.h"

#define MAX_COMPONENTS 7

// main class

class Clink {
public:
    u_number n_arcs[MAX_COMPONENTS]; // arcs of component
    u_number n; // overall crossings
    s_letter gw_[MAX_COMPONENTS][MAX<<1]; // Gauss word e.g. 1,-1,2,-2
    
    u_sign signs; // in binary, starting with 1st bit (not 0th)
    u_region regions[MAX+2]; // region in 64 bit binary (arcs in binary starting from 0th bit)
    u_number n_reg , reg_0, reg_1; // number of regions, index of 0 region, index of 1 region
    s16 ID; // knot ID
    u8 n_components;
    
    s_letter arc_s(int arc); // arc start
    s_letter arc_e(int arc); // arc end
    
    u_number arc_before(s_letter l); // arc left of l
    u_number arc_after(s_letter l); // arc right of l
    u_number get_region(u_number arc, s_small side); // get region adjacent to arc on side side
    u_number adjacent_region(u_number r_, u_number arc); // get region adjacent ot r_ wrt arc
    
    u_number region_index(u_region b); // get index of binary region b
    
    void swap_components(u8 c1, u8 c2); // swap link components
    
    bool is_disjunct_component(u8 c); // is components c not linked to other components
    
    void letter_indices(s_letter x, u8 *c1, u_number *i1, u8 *c2, u_number *i2); // finds occurences of x
    
    bool self_crossing(u8 c, s_letter x); // is x a self crossing in component c
    bool self_crossing(s_letter x); // is x a self crossing
    
    u_number turn(u_number arc, s_small *o, s_small direction); // turn to direction,o on arc, returns new arc
    
    void generate_regions(); // generate the regions from the gauss word
    
    u_region generate_region(u_number arc, s_small side); // generate the region on side adjacen to arc
    
    Cknot * knot(u8 c = 0); // convert the component c to a knot
    
    Cknot * seperate_component(u8 c);
    
    Cknot * seperate_disjunct_component(u8 c); 

    void add_disjunct_component(Cknot *K); // add knot K to the link

    void zerofyLink(); // initialization
    
    Clink * deep_copy();
    
    // constructors
    
    Clink() { zerofyLink(); } // create empty link
    Clink (Cknot *K); // create link from knot
    Clink (string s); // main constructor
    

};

ostream &operator<<(ostream &out, Clink L);

bool in_interval_(u_number a, u_number b, u_number x) { return ((x >= a) && (x <= b)); } // TODO: use existing interval function in common

void Clink::swap_components(u8 c1, u8 c2) { // swap link components
        
    u_number *arc_renum = new u_number[n*2];
    u8 c;
    u_number arc1, arc2, r, a;
    
    if (c1 > c2) { c=c1; c1=c2; c2=c;} // sort
    
    arc2 = arc1 = 0;
    for (c=0;c<c1;c++) arc1 += n_arcs[c];
    for (c=0;c<c2;c++) arc2 += n_arcs[c];
    
    // create region permutation table
    for (a=0; a<n*2; a++) {
        if (in_interval_(arc1,arc1+n_arcs[c1]-1,a))
            arc_renum[a] = a + arc2 - arc1  + n_arcs[c2] - n_arcs[c1];
        else
        if (in_interval_(arc2,arc2+n_arcs[c2]-1,a))
            arc_renum[a] = a + arc1 - arc2;
        else
        if (in_interval_(arc1+n_arcs[c1],arc2-1,a))
            arc_renum[a] = a + n_arcs[c2] - n_arcs[c1];
        else arc_renum[a] = a;
    }
    
    // swap regions
    for (r=0;r<n_reg;r++) permuteBits(&(regions[r]),arc_renum,n*2);
    
    // swap Gauss words
    for (a=0;a<MAXIMUM(n_arcs[c1],n_arcs[c2]);a++)
        swap(&gw_[c1][a],&gw_[c2][a]);
    
    a = n_arcs[c1]; n_arcs[c1] = n_arcs[c2]; n_arcs[c2] = a;
    
    delete arc_renum;
    
}

void Clink::add_disjunct_component(Cknot *K) { // ignore regions
    u_number i;
    for (i=0;i<K->n*2;i++) gw_[n_components][i] = (ABS(K->gw_[i]) + n) * SIGN(K->gw_[i]);
    signs |= (K->signs << K->n);
    n_arcs[n_components] = (K->n == 0 ? 1 : K->n*2);
    n += K->n;
    n_components++;
}

// main constructor, e.g. " 1 -2, 2 1 +-", components are seperated by commas

Clink::Clink(string s) {
    
    int i, p, ps;
    bool long_form;
    
    s = s + string(" ");
    
    
    long_form = (s.find("(") !=  string::npos); // are the gauss code regions in long form "(1 2 3)" instead of "123"
    
    zerofyLink();
    
    // find knot length
    p = (int)s.find(" +");
    if (p==string::npos) p = (int)s.find(" --");
    if (p==string::npos) p = (int)s.find(" -+");
    if (p==string::npos) p = (int)s.find(" - ");
    ps = p;
    
    n = (int)s.find(" ",p+1) - p - 1;

    p = (s.find("Link ")==string::npos ? 0 : 4);
    if (p) { ID = next_int(s,&p); p++; } // knot ID
    
    n_components = 0;
    
    // parse gauss word
    while (p < ps) {
        i = 0;
        while ((s[p] != ',') && (p < ps))
            gw_[n_components][i++]  = next_int(s,&p);
        n_arcs[n_components] = i;
        p+=2; n_components++;
    }
    p-=2;
    
    while (s[p] == ' ') p++;
    for (i=1;i<=n;i++) signs |= ((u_sign)char2sign(s[p++]) << i); // signs
    while (s[p] == ' ') p++;
    
    reg_0 = reg_1 = 0; // affine
   
}

// 1 component link
Clink::Clink(Cknot *K) {
    int i;
    zerofyLink();
    n = K->n;
    signs = K->signs;
    n_arcs[0] = K->n*2;
    n_components = 1;
    for (i=0;i<K->n_reg;i++) regions[i] = K->regions[i];
    for (i=0;i<K->n*2;i++) gw_[0][i] = K->gw_[i];
    n_reg = K->n_reg;
    reg_0 = K->reg_0;
    reg_1 = K->reg_1;
    ID = K->ID;
}


// get index of binary region b
u_number Clink::region_index(u_region b) {
    u_number r;
    for (r=0; r<n_reg ; r++ )
        if (regions[r]==b) return r;
    throw 20;//no region found
}

// finds ajacent region to r wrt arc
u_number Clink::adjacent_region(u_number r_, u_number arc) {
    u_number r;
    for (r=0; r<n_reg ; r++ )
        if (( BIT(regions[r],arc) ) && (r != r_)) return r;
    throw 20;//no region found
}

// finds occurences of x in the GW
void Clink::letter_indices(s_letter x, u8 *c1, u_number *i1, u8 *c2, u_number *i2) {
    u_number i = 0;
    u8 c = 0;
    *c1 = *c2 = 0xFF;

    x = ABS(x);
    
    for (c = 0; c < n_components; c++) {
        for (i = 0; i< n_arcs[c];i++) {
            if (ABS(gw_[c][i]) == x) {
                if (*c1 == 0xFF) {*i1 = i; *c1 = c; }
                else {
                    if (*c2 == 0xFF) {*i2 = i; *c2 = c;}
                    else throw 20;
                }
                
                if ((*c1 != 0xFF) && (*c2 != 0xFF)) return;
            }
        }
    }
    
    // not found
    throw 20;
}

Clink * Clink::deep_copy() {
    Clink *L = new Clink();
    L->n = n;
    memcpy(L->n_arcs,n_arcs,sizeof(n_arcs));
    L->n_components = n_components;
    memcpy(L->gw_,gw_,sizeof(gw_));
    L->signs = signs;
    memcpy(L->regions,regions,sizeof(regions));
    L->n_reg = n_reg;
    L->reg_0 = reg_0; L->reg_1 = reg_1;
    L->ID = ID;
    return L;
}

// generate region on side adjacent to arc
u_region Clink::generate_region(u_number arc, s_small side) {
    u_region b = (u_region)0;
    s_small o = side;
    do { // travel through knot
        b |= ((u_region)1 << arc); // add arc
        arc = turn(arc,&o,RIGHT); // find next arc
    } while (!BIT(b,arc));
    return b;
}

// get region on side adjacent to arc
u_number Clink::get_region(u_number arc, s_small side) {
    u_region b = generate_region(arc,side);
    for (u_number r=0;r<n_reg;r++) if (regions[r] == b) return r;
    throw 20;
}

// merge regions adjacent to arc
void merge_regions(u_region *regions, u_region arcs, u_number *n_reg, u_number *reg_0, u_number *reg_1) {
    u_region b;
    u_number arc;
    int r, r_, r__;
    
    // merge
    while (arcs) {
        arc = firstBitSet(arcs); // first non-processed arc
        b = (u_region)0;
        for (r=0;r<*n_reg;r++) if (BIT(regions[r],arc)) b |= regions[r];
        for (r=0;r<*n_reg;r++) if (BIT(regions[r],arc)) regions[r] = b;
        arcs ^= (u_region)1 << arc; // remove arc from set
    }
    
    // delete duplicate regions
    for (r=*n_reg-1; r >= 0; r--) {
        for (r_=0;r_<r;r_++) if (regions[r] == regions[r_]) break;
        
        if (r_ < r) { // same region found in lower regions, remove r, keep r_
            if (*reg_0 == r) *reg_0 = r_;
            if (*reg_1 == r) *reg_1 = r_;
            // delete region
            for (r__=r+1;r__<*n_reg;r__++) regions[r__-1] = regions[r__];
            if (*reg_1 > r) (*reg_1)--;
            if (*reg_0 > r) (*reg_0)--;
            (*n_reg)--;
            
        }
    }
    
}

// get maximal value of integer
template<class T>
inline bool max_value(const T t) { return 0x7f; /*std::numeric_limits<T>::max();*/}

void shift_arcs(Cknot *K, u_number shift) {
    for (u_number r = 0; r < K->n_reg; r++) K->regions[r] >>= shift;
}

void zerofy_arcs_above(Clink *L, u_number shift) {
    for (u_number r = 0; r < L->n_reg; r++) L->regions[r] &= fbf_reg(shift);
}

// normalize GC of a knot
void normalize_gauss(Cknot *K) {
    
    if (K->n == 0) return; // nothing to do.
    u_number i, next = 0, renum[MAX];
    s_letter l;
    u8 b_renum[MAX];
    u_sign signs;

    memset(b_renum, 0, sizeof(b_renum));
    for (i=0; i < K->n2; i++) {
        if (!b_renum[l = ABS(K->gw_[i])]) {renum[l] = ++next; b_renum[l] = 1;}
        K->gw_[i] = SIGN(K->gw_[i]) * renum[l];
    }

    // renumerate signs
    signs = K->signs;
    K->signs = 0;
    for (i=0;(int)i<MAX;i++)
        if (b_renum[i]) K->signs |= (BIT(signs,i) << renum[i]);
}

// normalize GC of a link
void normalize_gauss(Clink *L) {
    
    if (L->n == 0) return; // nothing to do.
    u_number i, next = 0, renum[MAX];
    s_letter l;
    u8 c;
    u8 b_renum[MAX];
    u_sign signs;
    memset(b_renum, 0, sizeof(b_renum));
    for (c=0;c<L->n_components;c++) {
        for (i=0;i<L->n_arcs[c];i++) {
            l = ABS(L->gw_[c][i]);
            if (!b_renum[l]) {renum[l] = ++next; b_renum[l] = 1;}
            L->gw_[c][i] = SIGN(L->gw_[c][i]) * renum[l];
        }
    }
    
    signs = L->signs;
    L->signs = 0;
    for (i=0;(int)i<MAX;i++) {
        if (b_renum[i]) L->signs |= (BIT(signs,i) << renum[i]);
    }
}


void delete_crossings_and_renumerate_arcs(Cknot *K, u_sign crossings) {
    
    u_number r,i, pos;
    u_number *arc_renum = new u_number[K->n*2]; // ?
    u_number dummy = 0;
    u_number shift, c_arcs;
    shift = 0;
    
    for (i=0;i<K->n*2;i++) arc_renum[i] = max_value(dummy);
    
    for (i=0,c_arcs=0;i<K->n*2;i++) c_arcs += ( BIT(crossings, ABS(K->gw_[i])) ? 0 : 1);
    
    pos = 0;
    
    for (i=0;i<K->n*2;i++) { // gw of component
        if (!BIT(crossings, ABS(K->gw_[i]))) {// letter not in deleting crossings
            K->gw_[pos] = K->gw_[i]; // put letter to position pos_
            pos++;
        }
        arc_renum[i] = MOD(pos-1,K->n*2);
    }
    
    for (i=0;i<K->n*2;i++) if (arc_renum[i] >= c_arcs) arc_renum[i] = c_arcs-1;
    
    for (r=0;r<K->n_reg;r++) {
        permuteBits(&(K->regions[r]),arc_renum,K->n*2);
        K->regions[r] &= fbf_reg((K->n*2+1));
    }
    
    
    K->n = c_arcs>>1; //MOD(pos_+1,L->n_arcs[c]);
    K->n2 = K->n*2;
    
    normalize_gauss(K);
    
    delete arc_renum;
    
}


void delete_crossings_and_renumerate_arcs(Clink *L, u_sign crossings) {
    
    u_number r,i, pos;
    u_number *arc_renum = new u_number[L->n*2]; // ?
    u_number shift, c_arcs, shift_;
    u8 c;
    shift = 0;
    shift_ = 0;
    
    for (i=0;i<L->n*2;i++) arc_renum[i] = -1; //max_value(dummy);
    
    for (c=0;c<L->n_components;c++) { // all components
        
        for (i=0,c_arcs=0;i<L->n_arcs[c];i++) c_arcs += ( BIT(crossings, ABS(L->gw_[c][i])) ? 0 : 1);
        
        pos = 0;
        
        for (i=0;i<L->n_arcs[c];i++) { // gw of component
            if (!BIT(crossings, ABS(L->gw_[c][i]))) {// letter not in deleting crossings
                L->gw_[c][pos] = L->gw_[c][i]; // put letter to position pos_
                pos++;
            }
            if (c_arcs >0 ) arc_renum[i+shift_] = MOD(pos-1,c_arcs) + shift; else
                arc_renum[i+shift_] = pos + shift -1;
        }
        shift_ += L->n_arcs[c];
        L->n_arcs[c] = c_arcs; //MOD(pos_+1,L->n_arcs[c]);
        shift += L->n_arcs[c];
    }
    
    for (r=0;r<L->n_reg;r++) {
        permuteBits(&(L->regions[r]),arc_renum,L->n*2);
        L->regions[r] &= fbf_reg((L->n*2+1));
    }
    
    u_number n_all_arcs = 0;
    for (c=0;c<L->n_components;c++) n_all_arcs += L->n_arcs[c];
    L->n = n_all_arcs >> 1;
    
    normalize_gauss(L);
    
    delete arc_renum;
}

#define DEBUG_WIND false

// return the sum of winding directions (1, -1) of arcs
int winding_direction(Clink *L, u_region arcs) {
    
    u_number arc, r;
    u_region i_arcs = 0, i_arcs_;
    int direction;
    
    i_arcs = arcs; //interior arcs
    
    for (arc=0; arc < L->n*2; arc++) if (BIT(arcs,arc)) i_arcs |= L->regions[ L->get_region(arc,RIGHT)]; // arcs on the left
    

    do { 
        i_arcs_ = i_arcs;
        for (r=0;r<L->n_reg;r++)
            if (!(L->regions[r] & arcs) && (L->regions[r] & i_arcs)) { i_arcs |= L->regions[r];} //cout << "r" << (int)r;}
    } while (i_arcs != i_arcs_);

    
    if (i_arcs == arcs) {
        
        r = L->get_region(firstBitSet(arcs),RIGHT);
        if ((r==L->reg_0) && (r==L->reg_1)) return 0;
        if (r==L->reg_0) return -1;
        if (r==L->reg_1) return 1;
        
        if ((r!=L->reg_0) && (r!=L->reg_1)) return 0; // this shouldn't happen
        
        r = L->get_region(firstBitSet(arcs),LEFT);
        if ((r==L->reg_0) && (r==L->reg_1)) return 0;
        if (r==L->reg_0) return 1;
        if (r==L->reg_1) return -1;

        throw 20;
        
        
    } else {
        
        i_arcs_ = (i_arcs | arcs) ^ arcs;
        
        direction = 0;
        if (i_arcs_ & L->regions[L->reg_0]) direction--;
        if (i_arcs_ & L->regions[L->reg_1]) direction++;
        
        return direction;
        
    }
    
}

#define DEBUG_S false

// to be added
Cknot * Clink::seperate_disjunct_component(u8 c) {
    return NULL;
}

bool Clink::is_disjunct_component(u8 c) {
    bool b  = true;
    for (u_number i = 0; i < n_arcs[c]; i++) b &= self_crossing(c, ABS(gw_[c][i]));
    return b;
}

Cknot * Clink::seperate_component(u8 c) {
    
    Cknot * K; // convert last component to knot
    int i;
    u_number arc_0; // start arc of region c
    u_region arcs_k, arcs_l;//, arcs_c;
    u_sign crossings; // linked crossings
    int direction_k, direction_l;
    bool disj_c;
    
    for (arc_0=0,i=0;i<c;i++) arc_0 += n_arcs[i];
    
    // get arcs
    arcs_k = fbf_reg(n_arcs[c]) << arc_0; // arcs of the knot
    arcs_l = fbf_reg((n*2)) ^ arcs_k; // arcs of the link, all expect those of the knot
    
    // get directions in case we get trivial components
    
    disj_c = is_disjunct_component(c);
    if (!disj_c) {
        direction_k = winding_direction(this, arcs_k);
        direction_l = winding_direction(this, arcs_l);
    }

    K = knot(c); // last component to knot

    //get linked crossings
    crossings = 0;
    for (i=0;i<n_arcs[c];i++) if (!self_crossing(c, gw_[c][i])) crossings |= (u_sign)1 << ABS(gw_[c][i]);
    
    merge_regions(K->regions, arcs_l, &K->n_reg, &K->reg_0, &K->reg_1);
    merge_regions(regions, arcs_k, &n_reg, &reg_0, &reg_1);
    
    shift_arcs(K, arc_0);
    
    delete_crossings_and_renumerate_arcs(K, crossings);
    
    zerofy_arcs_above(this, arc_0);

    delete_crossings_and_renumerate_arcs(this, crossings);

    // trivial components
    if (K->n == 0) {
        if (disj_c) throw 20;
        K->regions[0] = K->regions[1] = 1;
        if (K->reg_0 != K->reg_1) {
            K->reg_0 = (direction_k == 1 ? 0 : 1);
            K->reg_1 = (direction_k == 1 ? 1 : 0);
        } else K->reg_0 = K->reg_1 = 0;
    }
    
    if (n == 0) {
        if (disj_c) throw 20;
        regions[0] = regions[1] = 1;
        if (n_components != 1) { cout << *this << endl; throw 20; }
        if (reg_0 != reg_1) {
            reg_0 = (direction_l == 1 ? 0 : 1);
            reg_1 = (direction_l == 1 ? 1 : 0);
        } else reg_0 = reg_1 = 0;
    }
    
    return K;
    
}

// convert component c to a knot
Cknot * Clink::knot(u8 c) {
    
    if ((c == 0) && (n_components != 1)) throw 20; // can convert if only one component
    if ((c != 0) && (c != n_components-1)) throw 20; // or can extract only the last component;
    
    Cknot *K = new Cknot();
    K->ID = ID;
    K->n = n_arcs[c]/2;
    K->n2 = n_arcs[c];
    memcpy(K->gw_, gw_[c], sizeof(K->gw_));
    K->n_reg = n_reg;
    K->signs = signs;
    if (n_reg >= MAX+2) throw 20; // too much regions
    for (int i = 0; i<n_reg; i++) K->regions[i] = regions[i];
    K->reg_0 = reg_0;
    K->reg_1 = reg_1;
    
    n_components--;
    
    return K;
}

// is x a self crossing on c?
bool Clink::self_crossing(u8 c, s_letter x) {
    int i, count = 0;
    for (i=0;i<n_arcs[c]; i++) if (ABS(x) == ABS(gw_[c][i])) count++;
    if (count == 2) return true;
    if (count == 1) return false;
    throw 20; // crossing not found
}

// is x a self crossing in the link?
bool Clink::self_crossing(s_letter x) {
    int i, c, count;
    for (c=0;c<n_components;c++) {
        for (count=0,i=0;i<n_arcs[c]; i++) if (ABS(x) == ABS(gw_[c][i])) count++;
        if (count == 2) return true;
        if (count == 1) return false;
    }
    throw 20; // crossing not found
}

u_number Clink::arc_before(s_letter l) {
    int i, arc_add = 0;
    u8 c;
    for (c=0;c<n_components;c++) {
        for (i=0; i< n_arcs[c];i++)
            if (gw_[c][i] == l) return MOD(i+n_arcs[c]-1,n_arcs[c])+arc_add;
        
        arc_add += n_arcs[c];
        
    }
    throw 20; // letter not found
}

u_number Clink::arc_after(s_letter l) {
    u_number i, arc_add = 0;
    u8 c;
    for (c=0;c<n_components;c++) {
        for (i=0; i< n_arcs[c];i++)
            if (gw_[c][i] == l) return i+arc_add;
        arc_add += n_arcs[c];
    }
    throw 20; // letter l not found
}

s_letter Clink::arc_s(int arc) {
    int c = 0;
    arc = MOD(arc, n*2); // real mod?
    while (arc >= n_arcs[c]) { arc -= n_arcs[c]; c++; } // find component of arc
    return gw_[c][arc];
}

s_letter Clink::arc_e(int arc) {
    int c = 0;
    arc = MOD(arc, n*2); // real mod?
    while (arc >= n_arcs[c]) { arc -= n_arcs[c]; c++; } // find component of arc
    return gw_[c][MOD((arc+1),n_arcs[c])];
}

u_number Clink::turn(u_number arc, s_small *o, s_small direction) { // *o=0 - towards orientation, *o=1 - backwards
    s_letter x = ((*o) ? arc_s(arc) : arc_e(arc));
    *o = (s_small)SIGNB(x) ^ (s_small)(*o) ^ (s_small)(int)BIT(signs,(u_sign)ABS(x)) ^ (s_small)1 ^ (s_small)direction; // new direction
    return ((*o) ? arc_before(-x) : arc_after(-x));
}

// generate link regions from Gauss code
void Clink::generate_regions() {
    
    u_region usedArcsR = fbf_reg(n<<1), usedArcsL = fbf_reg(n<<1); // arcs that are not used
    u_number arc;
    s_small o;
    
    memset(regions,0,sizeof(regions));
    n_reg = 0;
    
    while (usedArcsR | usedArcsL) { // if unused arcs exist
        
        arc = firstBitSet(usedArcsR | usedArcsL);// find 1st available arc
        
        o = !BIT(usedArcsR,(u_region)arc); // 0 positive (towards orientation), - negative (opposite of orientation)

        while (1) { // travel through knot
            
            (o ? usedArcsL : usedArcsR) ^= ((u_region)1 << arc); // remove arc
            
            regions[n_reg] |= ((u_region)1 << arc); // add arc to region
            
            arc = turn(arc,&o,RIGHT); // find next arc
            
            if (BIT(regions[n_reg],arc)) { // close loop
                n_reg++;
                if (n_reg > n*2+2 + n_components) throw 20; // too many regions
                break;
            }
        }
    }
}

void Clink::zerofyLink() {
    n  = 0;
    memset(n_arcs,0,sizeof(n_arcs));
    memset(gw_,0,sizeof(gw_));
    signs = 0;
    memset(regions,0,sizeof(regions));
    n_reg = 0;
    reg_0 = reg_1 = -1;
    ID = -1;
    n_components = 0;
}

ostream &operator<<(ostream &out, Clink * L) {// print link
    int i, c;
    
    out << "Link " << (int)L->ID << ":";
    
    for (c=0;c<L->n_components;c++) {
        for (i=0;i<L->n_arcs[c];i++) out << ' '<<(int)L->gw_[c][i];
        if (c < L->n_components-1) out << ',';
    }
    cout << ' ';
    for (i=1;i<=L->n;i++) out << sign2char((int)BIT(L->signs,i));
    cout << ' ';
    
    for (i=0;i<L->n_reg;i++) {
        if ((i==L->reg_0) | (i==L->reg_1)) out << ((i==L->reg_0)?"*":"") << ((i==L->reg_1)?"&":"");
        
        if (L->n > 14) {
        cout << bitset_extended(L->regions[i]) << " ";
        } else {
        cout << bit_set(L->regions[i]) << ' ';
        }
    }
    
    return out;
}
    
ostream &operator<<(ostream &out, Clink L) {// print knot
    out << &L;
    return out;
}

#endif
