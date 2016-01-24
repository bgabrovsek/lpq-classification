//  HOMFLY.h
//  Created by Boštjan on 4/22/13.
//
//  Calculation of the HOMLFY skein module of a knot in the solid torus
//
//  using HOMLFY relation:
//  v^-1 H(L+) - v H(L-) - z H(L0) = 0
//  L- = v^-2 L+ - zv^-1 L0
//  L+ = v^2 L- + vz L0
//  L0 = v^-1 z^-1 L+ - v z^-1 L-
//  unknot = v^-1z^-1 - vz^-1
//
//  for the trivial non-affine knot we use the convention:
//  *0 &1 - counterclockwise direction (HSM = x)
//  &1 *0 - clockwise direcion

#ifndef Lpq_HOMFLY_h
#define Lpq_HOMFLY_h

#include "multivariate_laurent.h"
#include "link.h"
#include "knot.h"
#include "reidemeister_moves.h"

#define GEN_SIGN 1 // positive crossings in generators (we use the choice of basis that a (p,1)-knot with positive crossings is a generator)

// all crossing signs in K are the same?
bool all_signs_eq(Cknot *K, s_small sign) {
    int b = true;
    for (int i=1; i<=K->n; i++) if (BIT(K->signs, i) != (u_sign)sign) b = false;
    return b;
}

Cmultivariate * HOMFLY(Clink *L);
Cmultivariate * HOMFLY(Cknot *K);

// calculates the HSM of a knot with at most 1 crossing by hand (for optimization, since these knots often appear in the calculation of the HSM)
Cmultivariate * small_HOMFLY(Cknot *K) {
    
    if (K->n == 0) { // HSM of a trivial knot
        if (K->reg_0 == K->reg_1) return new Cmultivariate("- vz^-1 + v^-1z^-1"); // trivial affine
        if ((K->reg_0 == 0) && (K->reg_1 == 1)) { /*cout << "&";*/ return new Cmultivariate("N"); } // t_1
        if ((K->reg_0 == 1) && (K->reg_1 == 0)) return new Cmultivariate("L"); // t_-1
        throw 20;
        
    }
    
    if (K->n == 1) {
        if (K->reg_0 == K->reg_1) return new Cmultivariate("- vz^-1 + v^-1z^-1"); // trivial affine
        if ((count_bits(K->regions[K->reg_0]) == 1) && (count_bits(K->regions[K->reg_1]) == 1)) { // t_2, t_-2
            
            if (K->signs & 2) { // negative
                if (K->gw_[firstBitSet(K->regions[K->reg_0])] > 0) { /*cout << "%";*/ return new Cmultivariate("v^-2O - zv^-1N^2"); } // v^-2 t_2 - zv^-1 t_1^2
                if (K->gw_[firstBitSet(K->regions[K->reg_0])] < 0) return new Cmultivariate("v^-2K - zv^-1L^2"); // v^-2 t^-2 - zv^-1t^-1^2
                throw 20;
            } else { // positive
                if (K->gw_[firstBitSet(K->regions[K->reg_0])] > 0) return new Cmultivariate("K"); // t_-2
                if (K->gw_[firstBitSet(K->regions[K->reg_0])] < 0) return new Cmultivariate("O"); // t_2
                throw 20;
            }
        }
        
        if ((count_bits(K->regions[K->reg_0]) == 1) && (count_bits(K->regions[K->reg_1]) == 2)) { // t_1, t_-1
            if (K->regionSide(K->reg_0, firstBitSet(K->regions[K->reg_0])) == RIGHT) { /*cout << "*";*/ return new Cmultivariate("L"); }
            else return new Cmultivariate("N");
        }
        
        if ((count_bits(K->regions[K->reg_0]) == 2) && (count_bits(K->regions[K->reg_1]) == 1)) { // t_1, t_-1
            if (K->regionSide(K->reg_1, firstBitSet(K->regions[K->reg_1])) == LEFT) return new Cmultivariate("L");
            else return new Cmultivariate("N");
        }
        throw 20;
    }
    throw 20;
}

// is K represents a HSM generator, return the HSM
Cmultivariate * return_generator(Cknot *K) {
    
    Cmultivariate *poly = NULL;
    u_number i, shift;
    Cterm *t;
    
    poly = new Cmultivariate();
    
    if (K->n == 0) { // unknot
        
        if (K->reg_0 == K->reg_1) { // affine unknot, "v^-1z^-1 - vz^-1"
            t = new Cterm(1); t->power[reverse_variable['v']] = -1; t->power[reverse_variable['z']] = -1;
            *poly += *t;
            delete t;
            t = new Cterm(1); t->power[reverse_variable['v']] = 1; t->power[reverse_variable['z']] = -1;
            *poly += *t;
            delete t;
            return poly;
        }
        
        if ((K->reg_0 == 0) && ((K->reg_1 == 1))) { // nonaffine positive unknot, "t_1"
            *poly += *(new Cterm(1,MID_VAR+1,1));
            return poly;
        }
        
        if ((K->reg_0 == 0) && ((K->reg_1 == 1))) { // nonaffine positive unknot, "t_-1"
            *poly += *(new Cterm(1,MID_VAR-1,1));
            return poly;
        }
        
        throw 20;
    }
    
    // generator of higher degree?
    if (K->reg_0 == K->reg_1) return NULL; // affine
    if (!all_signs_eq(K,SIGN2B(GEN_SIGN))) return NULL; // all signs equal?
    if ((count_bits(K->regions[K->reg_0]) != 1) || (count_bits(K->regions[K->reg_1]) != 1)) return NULL; // selected regions have 1 arc?
    
    // is the gauss code a palindrom with reg_0 as the center?
	shift = firstBitSet(K->regions[K->reg_0]);
	for (i=1;i<=K->n;i++) {
		if (K->gw(shift)*K->gw(shift+i) >= 0) return NULL;
		if (K->gw(shift+i) + K->gw(shift-i+1) != 0) return NULL;
	}
    
    char ch = MID_VAR + ((K->gw(shift) > 0) ? -1 : 1)*(char)(K->n+1);
    t = new Cterm(1,ch,1);
    *poly += *t;
   // delete t;
    return poly;
}

// does K represent a HSM generator?
bool is_generator(Cknot *K) {

    u_number i, shift;

    if (K->n == 0) { // unknot
        if (K->reg_0 == K->reg_1) return YES; // affine unknot, "v^-1z^-1 - vz^-1"
        if ((K->reg_0 == 0) && ((K->reg_1 == 1))) return YES; // nonaffine positive unknot, "t_1"
        if ((K->reg_0 == 0) && ((K->reg_1 == 1))) return YES; // nonaffine positive unknot, "t_-1"
        throw 20;
    }
    
    // generator of higher degree?
    if (K->reg_0 == K->reg_1) return NO; // affine
    if (!all_signs_eq(K,SIGN2B(GEN_SIGN))) return NO; // all signs equal?
    if ((count_bits(K->regions[K->reg_0]) != 1) || (count_bits(K->regions[K->reg_1]) != 1)) return NO; // selected regions have 1 arc?
    
    // is the gauss code a palindrom with reg_0 as the center?
	shift = firstBitSet(K->regions[K->reg_0]);
	for (i=1;i<=K->n;i++) {
		if (K->gw(shift)*K->gw(shift+i) >= 0) return NO; //return NULL;
		if (K->gw(shift+i) + K->gw(shift-i+1) != 0) return NO; //return NULL;
	}

    return YES;
}

// manipulation of links

// returns -1 if last component of L should be saperated by lowering it with undercrossings, otherwise 1, 0 if already seperated
int last_component_height(Clink *L) { 
    
    int i, overcrossing = 0, undercrossings = 0;
    s_letter x;
    u8 c = L->n_components-1;
    for (i=0;i<L->n_arcs[c];i++) {
        x = L->gw_[c][i]; // letter
        if (!L->self_crossing(c, x)) (x>0 ? overcrossing : undercrossings)++;
    }
    if ((overcrossing == 0) || (undercrossings == 0)) return 0;
    return (overcrossing >= undercrossings ? 1 : -1);
}

// if height = 1, find first undercrossing and vice versa if height = -1
s_letter first_linked_crossing_on_last_component(Clink *L, int height) { 
    u8 c = L->n_components-1;
    s_letter x;
    for (u_number i=0;i<L->n_arcs[c];i++) {
        x = L->gw_[c][i];
        if ((!L->self_crossing(c, x)) && (height * x < 0)) return ABS(x);
    }
    throw 20;
}

// joins components by making a 0-skein join, leaves L as it is, returns the new link
Clink *join_smooth(Clink *L, s_letter x) {
    Clink *L_;
    u8 c, c1, c2; // components to join
    u_number i, i1, i2, p;
    u_number *arc_renum = new u_number[L->n*2];
    u_number arc_1_0, arc_2_0, arc_0, r_a, r_b, r;
    
    L_ = L->deep_copy();
    
    L->letter_indices(x, &c1, &i1, &c2, &i2); // arcs and components to join
    
    arc_1_0 = arc_2_0 = 0; // find starting arcs
    for (c=0;c<c1;c++) arc_1_0 += L->n_arcs[c];
    for (c=0;c<c2;c++) arc_2_0 += L->n_arcs[c];
    
    for (i=0;i<L->n*2;i++) arc_renum[i] = -1;
    
    p = 0;
    for (i=0;i<i1;i++) { // 1st part of 1st component
        L_->gw_[c1][p] = L->gw_[c1][i];
        arc_renum[i + arc_1_0] = p + arc_1_0;
        p++;
    }
    
    for (i=1;i<L->n_arcs[c2];i++) { // 2nd component
        L_->gw_[c1][p] = L->gw_[c2][MOD(i+i2,L->n_arcs[c2])];
        arc_renum[MOD(i+i2,L->n_arcs[c2])  + arc_2_0] = p  + arc_1_0;
        p++;
    }
    
    for (i=i1+1;i<L->n_arcs[c1];i++) { // 2nd part of 1st component
        L_->gw_[c1][p] = L->gw_[c1][i];
        arc_renum[i + arc_1_0] = p  + arc_1_0;
        p++;
    }

    arc_0 = 0;
    u_number arc_0_ = 0;
    
    for (c=0;c<L->n_components;c++) {
        if ((c != c1) && (c != c2)) {
            for (i=0;i<L->n_arcs[c];i++) {
                arc_renum[i+arc_0] = (i + arc_0_) + ((i+arc_0 > arc_1_0) ? -1 : 0) + ((i+arc_0 > arc_2_0) ? -1 : 0);
            }
        }
        arc_0 += L->n_arcs[c];
        arc_0_ += L->n_arcs[c];
        if (c == c1) arc_0_ += L->n_arcs[c2] -1;
    }
    
     L_->n_arcs[c1] = L->n_arcs[c1] + L->n_arcs[c2] -2;
    
    // join regions
    if (L->gw_[c1][i1] == ABS(x)) {
        if (BIT(L->signs,ABS(x))) {
            r_a = L->get_region(MOD(i1-1,L->n_arcs[c1])+arc_1_0, LEFT);
            r_b = L->get_region(i1+arc_1_0, RIGHT);
        } else {
            r_a = L->get_region(MOD(i1-1,L->n_arcs[c1])+arc_1_0, RIGHT);
            r_b = L->get_region(i1+arc_1_0, LEFT);
        }
    } else
        if (L->gw_[c2][i2] == ABS(x)) {
            if (BIT(L->signs,ABS(x))) {
                r_a = L->get_region(MOD(i2-1,L->n_arcs[c2])+arc_2_0, LEFT);
                r_b = L->get_region(i2+arc_2_0, RIGHT);
            } else {
                r_a = L->get_region(MOD(i2-1,L->n_arcs[c2])+arc_2_0, RIGHT);
                r_b = L->get_region(i2+arc_2_0, LEFT);
            }
        } else throw 20;
    
    // join regions
    
    if (r_b < r_a) {r = r_b; r_b = r_a; r_a = r;} // sort
    
    L_->regions[r_a] |= L_->regions[r_b];
    
    for (r=r_b;r<L->n_reg-1;r++) L_->regions[r] = L_->regions[r+1];

    if (L_->reg_0 == r_b) L_->reg_0 = r_a; else
    if (L_->reg_0 > r_b) L_->reg_0--;
    if (L_->reg_1 == r_b) L_->reg_1 = r_a; else
    if (L_->reg_1 > r_b) L_->reg_1--;
    
    L_->n_reg--;
    
    // remove component
    for (c=0;c<L->n_components;c++) {
        if (c > c2) {
            memcpy(&L_->gw_[c-1], &L_->gw_[c], sizeof(L_->gw_[c]));
            L_->n_arcs[c-1] = L_->n_arcs[c];
        }
    }
    
    L_->n_components--;
    
    for (r=0;r<L_->n_reg;r++) {
        permuteBits(&L_->regions[r], arc_renum, L->n*2);
        L_->regions[r] &= fbf_reg(L_->n*2+1);
    }
    
    L_->n--;
    normalize_gauss(L_);
    delete arc_renum;
    return L_;
}

// changes the sign of c in K
void change_crossing_signature(Cknot *K, s_letter c) {
    K->signs ^= ((u_sign)1 << ABS(c));
    swap( &(K->gw_[K->gw_pos(c)]), &(K->gw_[K->gw_pos(-c)]));
}

// does crossing x represnt a kink?
bool is_kink_at_crossing(Clink *L, s_letter x) {
    u8 c1, c2; 
    u_number i1, i2;
    L->letter_indices(x, &c1, &i1, &c2, &i2);
    return ((c1 == c2) && ( (i1+1==i2) || ( i1 + L->n_arcs[c1] == i2+1 ) ));
}

// smoothens a punctured kink in L, return the unknot in K and the link with removed component in A
// TODO: if not puctured kink, then adjust K
void kink_smooth(Clink *L, s_letter x, Clink **A, Cknot **K) {
    
    u8 c1, c2; // components to join
    u_number i1, i2, reg, outer_reg;
    u_number arc;
    int direction;
    bool affine;
    
    L->letter_indices(x, &c1, &i1, &c2, &i2);
    
    arc = (((i1 == 0) && (i2 == L->n_arcs[c1]-1)) ? (i2) : (i1)); // kink arc
    
    direction = winding_direction(L, (u_region)1<<arc); // positive or negative generator
    
    reg = L->region_index((u_region)1<<arc); // get region of kink
    
    outer_reg = L->adjacent_region(reg,arc); // get region adjacent to kink
    
    affine = !((reg == L->reg_0) ^ (reg == L->reg_1)); // xnor, is kink affine?
    
    if (c1 != c2) throw 20; // kink should be on same component
    
    *A = L->deep_copy();
    
    if (!affine) (reg==L->reg_0 ? (*A)->reg_0 : (*A)->reg_1) = outer_reg; // change dotted region
    
    removeRI(*A, c1, reg); // remove the kink
    
    // new unknot, get the direction (CCW or CW)
    (*K) = new Cknot();
    (*K)->n = 0; (*K)->n2 = 0;
    (*K)->n_reg = 2;
    if (!affine) {
        (*K)->reg_0 = (direction > 0 ? 0:1);
        (*K)->reg_1 = (direction > 0 ? 1:0);
    } else (*K)->reg_0 = (*K)->reg_1 = 0;
}

// join components by making a 0-skein join at x, leaves L as it is
Clink *split_smooth(Clink *L, s_letter x) {
    
    Clink *L_;
    u8 c, c1, c2; // components to join
    u_number i, i1, i2, p,p_;
    u_number *arc_renum = new u_number[L->n*2];
    u_number arc_1_0, arc_2_0, arc_3_0, arc_0, r_a, r_b, r;
    
    L_ = L->deep_copy();
    L->letter_indices(x, &c1, &i1, &c2, &i2); // arcs and components to join
    
    arc_1_0 = arc_2_0 = 0; // find starting arcs
    for (c=0;c<c1;c++) arc_1_0 += L->n_arcs[c];
    for (c=0;c<L->n_components;c++) arc_2_0 += L->n_arcs[c];
    c2 = L->n_components;
    arc_2_0 -= L->n_arcs[c1] - (i2-i1) +1; // elliminate
    arc_3_0 = L->n_arcs[c1] - (i2-i1) +1;
    
    for (i=0;i<L->n*2;i++) arc_renum[i] = -1;
    
    p = 0; p_ = 0; // positions of arcs
    for (i=i1+1; i<i2;i++) {
        L_->gw_[c1][p] = L->gw_[c1][i];
        arc_renum[i + arc_1_0] = p + arc_1_0;
        p++;
    }
    
    for (i=MOD(i2+1,L->n_arcs[c1]); i!= i1; i = MOD(i+1,L->n_arcs[c1])) {
        L_->gw_[c2][p_] = L->gw_[c1][i];
        arc_renum[i + arc_1_0] = p_ + arc_2_0;
        p_++;
    }
    
    // renumerate ramaining of arcs
    arc_0 = 0;
    for (c=0;c<L->n_components;c++) {
        if ((c != c1))
            for (i=0;i<L->n_arcs[c];i++)
                arc_renum[i+arc_0] = (i + arc_0) + ((i+arc_0 > arc_1_0) ? -arc_3_0 : 0);
        arc_0 += L->n_arcs[c];
    }
    
    L_->n_arcs[c2] = L_->n_arcs[c1] - (i2 - i1 + 1);
    L_->n_arcs[c1] = i2 - i1 - 1;
    
    // join regions
    if (L->gw_[c1][i1] == ABS(x)) {
        if (BIT(L->signs,ABS(x))) {
            r_a = L->get_region(MOD(i1-1,L->n_arcs[c1])+arc_1_0, LEFT);
            r_b = L->get_region(i1+arc_1_0, RIGHT);
        } else {
            r_a = L->get_region(MOD(i1-1,L->n_arcs[c1])+arc_1_0, RIGHT);
            r_b = L->get_region(i1+arc_1_0, LEFT);
        }
    } else
        if (L->gw_[c1][i2] == ABS(x)) {
            if (BIT(L->signs,ABS(x))) {
                r_a = L->get_region(MOD(i2-1,L->n_arcs[c1])+arc_1_0, LEFT);
                r_b = L->get_region(i2+arc_1_0, RIGHT);
            } else {
                r_a = L->get_region(MOD(i2-1,L->n_arcs[c1])+arc_1_0, RIGHT);
                r_b = L->get_region(i2+arc_1_0, LEFT);
            }
        } else throw 20;
    
    if (r_b < r_a) {r = r_b; r_b = r_a; r_a = r;} // sort
    if (r_b != r_a) {
        
        L_->regions[r_a] |= L_->regions[r_b];
        
        for (r=r_b;r<L->n_reg-1;r++)
            L_->regions[r] = L_->regions[r+1];
        
        if (L_->reg_0 == r_b) L_->reg_0 = r_a; else
        if (L_->reg_0 > r_b) L_->reg_0--;
        if (L_->reg_1 == r_b) L_->reg_1 = r_a; else
        if (L_->reg_1 > r_b) L_->reg_1--;

        L_->n_reg--;
    }
    
    //add component
    
    L_->n_components++;
    
    for (r=0;r<L_->n_reg;r++) {
        permuteBits(&L_->regions[r], arc_renum, L->n*2);
        L_->regions[r] &= fbf_reg(L_->n*2+1);
    }
    
    L_->n--;
    normalize_gauss(L_);
    
    delete arc_renum;
    return L_;
    
}

// performs a skein smoothening, decides either to perform a
// kink smooth, split smooth (a smoothening that splits the link) or join smooth ( a smoothening that joints the link)
Clink * skein_smooth(Clink *L, s_letter x, Cknot **K ) {
    *K = NULL;
    if (L->self_crossing(x)) {
        if (is_kink_at_crossing(L, x)) {
            Clink *L_smooth;
            kink_smooth(L, x, &L_smooth, K);
            return L_smooth;
        } else {
            return split_smooth(L, x);
        }
    } else {
        return join_smooth(L, x);
    }
}

// performs a skein crossing change relation
Clink * skein_pm(Clink *L_, s_letter x) { 
    int c, i, count = 0;
    x = ABS(x);
    
    Clink *L = L_->deep_copy();
    L->signs ^= ((u_sign)1 << x);
    
    for (c=L->n_components-1; c>=0; c--) // reverse, since last component has x
        for (i=0;i<L->n_arcs[c];i++)
            if (ABS(L->gw_[c][i]) == x) { L->gw_[c][i] = -L->gw_[c][i]; if (++count == 2) return L;}
    
    throw 20;
}

// is K an affine unknot?
bool affine_unknot(Cknot *K) { return ((K->n == 0) && (K->reg_0 == K->reg_1)); }

// the direction of the generator in case K is a generator
int generator_direction(Cknot *K) {
    if (K->n > 0)
        return ( K->gw(firstBitSet(K->regions[K->reg_0])) > 0 ? -1 : 1);
    else {
        if ((K->reg_0 == 0) && (K->reg_1 == 1)) return 1;
        if ((K->reg_0 == 1) && (K->reg_1 == 0)) return -1;
        throw 20;
    }
}

// make a simple guess of the crossing to perform a HOMLFY relation that simplifies the knot
s_letter simple_crossing2simlify(Cknot *K) {
    u_number r;
    for (r=0;r<K->n_reg;r++)
        if ((count_bits(K->regions[r]) == 2) && (r != K->reg_0) && (r != K->reg_1) &&
            (K->gw(firstBitSet(K->regions[r])) * K->gw(firstBitSet(K->regions[r])+1) < 0 )) { // one negative, one positive
            if (BIT(K->signs, ABS(K->gw(firstBitSet(K->regions[r])))))
                return ABS(K->gw(firstBitSet(K->regions[r])));
            else return ABS(K->gw(firstBitSet(K->regions[r])+1));
            
        }
    return -1;
}

// make a guess of the crossing to perform a HOMLFY relation that simplifies the knot
// changes each crossing and sees which of the crossings reduces the knot the most
s_letter crossing2simlify(Cknot *K) {
    Cknot * K_min, *Q;
    s_letter c, c_min = -1;
    K_min = K;
    
    for (c = 1; c <= K->n; c++) {
        Q = K->deepCopy();
        change_crossing_signature(Q, c);
        BST_shrink(&Q, ORIENTED,(int)(Q->n),1); // false = also increase
        if ( Q->n < K_min->n ) {
            delete Q;
            return c;
            K_min = Q;
            c_min = c;
            
        } else delete Q;
    } // crossing

    return c_min;
}

// span of the powers of the HOMLFY generators in p
void HOMFLY_span(Cmultivariate *p, int *a, int *b) {
    
    *a = 100;
    *b = -100;
    int c;
    p->simplify();
    
    for (list<Cterm*>::iterator it= p->terms.begin(); it != p->terms.end(); it++) {
        for(int i=0;i<n_v;i++)
            if ((is_capital(variables[i])) && ((*it)->power[i]!=0)) {
                c = (int)variables[i] - (int)MID_VAR;
                if (c < *a) *a = c;
                if (c > *b) *b = c;
            }
    }
}

// returns the HOMFLY skein module of K
Cmultivariate * HOMFLY(Cknot *K) {
    
    Cmultivariate *poly, *poly_pm, *poly_smooth; // the final poly, the poly after a crossing change, the poly after a smoothening
    s_letter crossing;
    s_small crossing_sign;
    
    Cknot *K_pm;
    Clink *L_smooth;
    
    // if (K->n > 18) BST_shrink(&K, ORIENTED, K->n*2,5);
    // poly = new Cmultivariate();
    
    if (K->n <= 1) { poly = small_HOMFLY(K); return poly; }
    
    if (is_generator(K)) {
        
        poly = new Cmultivariate();
        *poly += new Cterm(1,MID_VAR + (generator_direction(K))*(char)(K->n+1),1);
		
        return poly;
	}

    K->canonical(ORIENTED);

    BST_shrink(&K, ORIENTED, K->n, 6 ); // minimize the knot
    
    K->canonical(ORIENTED);
    
    if (K->n <= 1) { poly = small_HOMFLY(K); return poly;}
    
    if (is_generator(K)) {
        poly = new Cmultivariate();
        *poly += new Cterm(1,MID_VAR + (generator_direction(K))*(char)(K->n+1),1);
		return poly;
	}
    
    // find a good crossing for performing the skein crossing change relation
    
    crossing = -1; // no crossing choosen
    
    crossing = simple_crossing2simlify(K); // make simple guess
    
    if (crossing == -1) crossing = crossing2simlify(K); // make an advanced guess
    
    // make a guess by making all crossing either positive or negative (depending on the basis of the HSM choosen)
    if ( GEN_SIGN == -1) {
        if ((crossing == -1) && (K->n - count_bits(K->signs) )) {
            crossing = firstBitUnset(K->signs | (u_sign)1);
        }
    } else {
        if ((crossing == -1) && (count_bits(K->signs) )) {
            crossing = firstBitSet(K->signs); //ssss
        }
    }
    
    // if no sucessful guess made, use a random crossing
    if (crossing == -1) crossing = (rand()%K->n)+1;
    
    crossing_sign = B2SIGN(BIT(K->signs,(u_sign)crossing));
    
    // make the crossing change
    K_pm = K->deepCopy();
    change_crossing_signature(K_pm, crossing);
    L_smooth = new Clink(K);
    poly_pm = HOMFLY(K_pm);

    // make a smoothening
    Cknot *K_smooth;
    Clink *L_smooth_ = skein_smooth(L_smooth,ABS(crossing),&K_smooth);
    
    poly_smooth = HOMFLY(L_smooth_);
    
    if (K_smooth != NULL) {
        
        Cmultivariate *poly_K;
        
        poly_K = HOMFLY(K_smooth);
        *poly_smooth *= *poly_K;
        
        delete poly_K;
        delete K_smooth;
    }
    
    poly_smooth->multiplyByVariable(crossing_sign, 'z', 1);
    poly_smooth->multiplyByVariable(1, 'v', crossing_sign);
    poly_pm->multiplyByVariable(1, 'v', 2*crossing_sign);
    
    poly = *poly_smooth + *poly_pm;
 
    delete K_pm;
    delete L_smooth;
    delete L_smooth_;
    delete poly_pm;
    delete poly_smooth;
    
    return poly;
    
}

// are components c1 and c2 disjunct in L?
bool disjunct_components(Clink *L, u8 c1, u8 c2) {
    u_sign xb1, xb2;
    xb1 = xb2 = 0;
    int i;
    for (i=0;i<L->n_arcs[c1];i++) xb1 |= (u_sign)1 << ABS(L->gw_[c1][i]);
    for (i=0;i<L->n_arcs[c2];i++) xb2 |= (u_sign)1 << ABS(L->gw_[c2][i]);
    return (bool)( xb1 & xb2 );
}

// is the last component of L such, that after removing it, it will disconnect the other components?
s8 last_component_is_disconnecting(Clink *L) {
    s8 c;
    u_sign b, bl;
    int i;
    
    for (bl=0,i=0;i<L->n_arcs[L->n_components-1];i++) bl |= (u_sign)1 << ABS(L->gw_[L->n_components-1][i]);
    
    for (c=0;c<L->n_components-1;c++) {
        for (b=0,i=0;i<L->n_arcs[c];i++) b |= (u_sign)1 << ABS(L->gw_[c][i]);
        if (SUBSET(b,bl)) return c;
    }
    
    return -1;
}

// adjust the order of the components in such a way, that the last components is not disconnecting
void adjust_component_order(Clink *L) {
    s8 c;
    if (L->n_components > 2) {
        do {
            c = last_component_is_disconnecting(L);
            if (c != -1) L->swap_components(c,L->n_components-1);
        } while (c != -1);
    }
}

// returns the HOMLFY of the link L
Cmultivariate * HOMFLY(Clink *L) {
    
    if (L->n_components == 0) throw 20;
    
    Cmultivariate * poly_pm, *poly_smooth, *poly, *poly_L, *poly_K;
    Cknot * K_smooth, * K_last_component;
    Clink * L_tmp;
    Cknot * K_tmp;
    int height;
    int sign;
    s_letter x;
    
    // if only one component, calculate the HOMFLY of the knot
    if (L->n_components == 1) {
        K_tmp = L->knot();
        poly = HOMFLY(K_tmp);
        delete K_tmp;
        return poly;
    }

    adjust_component_order(L);
    height = last_component_height(L);

    if (height == 0) { // component is seperated

        K_last_component = L->seperate_component(L->n_components-1);
        
        poly_L = HOMFLY(L);
        
        poly_K = HOMFLY(K_last_component);
        
        poly = *poly_L * *poly_K;
    
        delete poly_L;
        delete poly_K;
        
        delete K_last_component;
        
        return poly;
        
        
    } else { // make skein to seperate components
        
        x = first_linked_crossing_on_last_component(L, height);
        sign = (BIT(L->signs,ABS(x)) ? -1 : 1);
        
        // smoothening skein term
        L_tmp = skein_smooth(L,x,&K_smooth);
        poly_smooth = HOMFLY(L_tmp);
        delete L_tmp;
        
        // change of crossing skein term
        L_tmp = skein_pm(L,x);
        poly_pm = HOMFLY(L_tmp);
        delete L_tmp;

        // skein relation variable shift
        
        if (K_smooth != NULL) {
            
            poly_K = HOMFLY(K_smooth);
            *poly_smooth *= *poly_K;
            
            delete poly_K;
            delete K_smooth;
            
        }
     
        poly_smooth->multiplyByVariable(sign, 'z', 1);
        poly_smooth->multiplyByVariable(1, 'v', sign);
        poly_pm->multiplyByVariable(1, 'v', 2*sign);
        
        poly = *poly_smooth + *poly_pm;
        poly->simplify();
        
        delete poly_smooth;
        delete poly_pm;
       // delete L; // is this safe?
        return poly;
        
    }
    
}

#endif
