//  reidemeister_moves.h
//  Created by Boštjan on 4/22/13.
//
//  Performing and searching for Reidemeister moves and flypes

#ifndef Lpq_reidemeister_moves_h
#define Lpq_reidemeister_moves_h

#include <list>

#include "knot.h"
#include "link.h"
#include "common.h"


// data structure of a Reidemeister move
class CRSite {
public:
    s_small type;
    u_number data[5]; // which regions is involved, which arcs are involved, type of move, type of reid. move
    CRSite(s_small t, u_number d0, u_number d1 = 0, u_number d2 = 0, u_number d3 = 0, u_number d4 = 0)
    { type = t; data[0] = d0; data[1] = d1; data[2] = d2; data[3] = d3; data[4] = d4;}
// type = REMOVE_R_I, data[0] = region
// type = CREATE_R_I, data[0] = region, data[1] = arc, data[2] = sign
// type = REMOVE_R_II, data[0] = region
// type = CREATE_R_II, data[0] = region, data[1] = over_arc, data[2] = under_arc
// type = MODIFY_R_III, data[0] = region, data[1] = over_arc, dara[2] = under_arc, data[3] = mixed arc
// type = FLYPE, data[0] = left arc, data[1] = left arc, data[2] = right arc, data[3] = right arc, data[4] = flype ouside or inside
};

// list of R-moves
typedef list<CRSite> t_reidemeister_list;

void quadrants(Cknot *K, s_letter c, u_number *r0, u_number *r1, u_number *r2, u_number *r3 ); // returns the quadrants of K at c
void deleteRegion(Cknot *K, u_number reg); // deletes region from region list
void deleteRegion(Clink *K, u_number reg); // deletes region from region list
void deleteCrossingFromGaussWord(Cknot *K, s_letter x);
void deleteCrossingFromGaussWord(Clink *L, u8 c, s_letter x);
void insertCrossingToGaussWord(Cknot *K, u_number pos1, u_number pos2, s_letter x1, s_letter x2); // if pos1 = pos2, inserts x1x2, otherwise inserts x1 at p1 and x2 at p2
void deleteArcsFromRegions(Cknot *K, u_number arc1, u_number arc2); // deletes arc1 and arc2 from all regions
void swapCrossingInGaussWord(Cknot *K, u_number pos_a, u_number pos_b, u_number pos_a_, u_number pos_b_, s_small turn_around); // swaps letters from pos_a -> pos_a_, pos_b -> pos_b_


// print all possible Reidemeister moves, not needed
ostream &operator<<(ostream &out, CRSite site) {
    switch (site.type) {
        case REMOVE_R_I  : out << "Remove RI: " << "region " << (int)site.data[0]; break;
        case CREATE_R_I  : out << "Create RI: " << "region " << (int)site.data[0] << ", arc " << (int)site.data[1] << ", sign " << (int)site.data[2]; break;
        case REMOVE_R_II : out << "Remove RII: " << "region " << (int)site.data[0]; break;
        case CREATE_R_II : out << "Create RII: " << "region " << (int)site.data[0] << ", over arc " << (int)site.data[1] << ", under arc " << (int)site.data[2] << " ds " << (int)site.data[3]; break;
        case MODIFY_R_III: out << "Modify RIII: "<<"region "<<(int)site.data[0]<<", arcs "<<(int)site.data[1]<<", "<<(int)site.data[2]<<", "<<(int)site.data[3]; break;
        case MODIFY_FLYPE: out << "Modify FLYPE: "<<"a "<<(int)site.data[0]<<", b "<<(int)site.data[1]<<", a_ "<<(int)site.data[2]<<", b_ "<<(int)site.data[3] << ((site.data[4])?" complement":""); break;
    }
    return out;
}

// print all possible Reidemeister moves
ostream &operator<<(ostream &out, list<CRSite> sites) {
    for (list<CRSite>::iterator it = sites.begin(); it != sites.end(); it++) out << (*it) << endl;
    return out;
}

// find and return arcs involved in a R-III moves
bool find_RIII_arcs(Cknot *K, u_number reg, u_number*over, u_number *under, u_number *mixed) {
    u_number arc;
    *over = *under = *mixed = infinity_number;
    
    arc = firstBitSet(K->regions[reg]);
    if ((K->gw(arc) > 0) && (K->gw(arc+1) > 0)) *over = arc; else if ((K->gw(arc) < 0) && (K->gw(arc+1) < 0)) *under = arc; else *mixed = arc;
    arc = secondBitSet(K->regions[reg]);
    if ((K->gw(arc) > 0) && (K->gw(arc+1) > 0)) *over = arc; else if ((K->gw(arc) < 0) && (K->gw(arc+1) < 0)) *under = arc; else *mixed = arc;
    arc = thirdBitSet(K->regions[reg]);
    if ((K->gw(arc) > 0) && (K->gw(arc+1) > 0)) *over = arc; else if ((K->gw(arc) < 0) && (K->gw(arc+1) < 0)) *under = arc; else *mixed = arc;
   
    if ((*over == infinity_number) || (*over == infinity_number) || (*over == infinity_number)) return false;
    
    if ((K->gw(*over) == -K->gw(*over+1)) || (K->gw(*under) == -K->gw(*under+1)) || (K->gw(*mixed) == -K->gw(*mixed+1))) return false;
    
    return true;
    
    /*if (( K->gw(firstBitSet(K->regions[reg])) > 0 ) && ( K->gw(firstBitSet(K->regions[reg])+1) > 0 )) over_arc = firstBitSet(K->regions[reg]);
    if (( K->gw(firstBitSet(K->regions[reg])) < 0 ) && ( K->gw(firstBitSet(K->regions[reg])+1) < 0 )) under_arc = firstBitSet(K->regions[reg]);
    if ( K->gw(firstBitSet(K->regions[reg])) * K->gw(firstBitSet(K->regions[reg])+1) < 0 ) mixed_arc = firstBitSet(K->regions[reg]);
    if (( K->gw(secondBitSet(K->regions[reg])) > 0 ) && ( K->gw(secondBitSet(K->regions[reg])+1) > 0 )) over_arc = secondBitSet(K->regions[reg]);
    if (( K->gw(secondBitSet(K->regions[reg])) < 0 ) && ( K->gw(secondBitSet(K->regions[reg])+1) < 0 )) under_arc = secondBitSet(K->regions[reg]);
    if ( K->gw(secondBitSet(K->regions[reg])) * K->gw(secondBitSet(K->regions[reg])+1) < 0 ) mixed_arc = secondBitSet(K->regions[reg]);
    if (( K->gw(thirdBitSet(K->regions[reg])) > 0 ) && ( K->gw(thirdBitSet(K->regions[reg])+1) > 0 )) over_arc = thirdBitSet(K->regions[reg]);
    if (( K->gw(thirdBitSet(K->regions[reg])) < 0 ) && ( K->gw(thirdBitSet(K->regions[reg])+1) < 0 )) under_arc = thirdBitSet(K->regions[reg]);
    if ( K->gw(thirdBitSet(K->regions[reg])) * K->gw(thirdBitSet(K->regions[reg])+1) < 0 ) mixed_arc = thirdBitSet(K->regions[reg]);
    */
    
}

// find possible sites to make Reidemeister moves and stores them into rList
// K: knot
// RT: flags which R-moves to find
// if rList = NULL only return true/false
// return: true if R-moves are found, false otherwise
bool findReidemeisterMoveSites(Cknot *K, u_number RT, t_reidemeister_list * rList = NULL) {
    
    int r, a, b, a_,b_, i, j;;
    u_number d[4];
    s_letter c;
    u_region cr_or, cr_xor, cr_or_, cr_xor_, regs, regs_;
    u_number reg_A, reg_B; // for flype
    s_small orient, reg_0_inside_flype, reg_1_inside_flype;

    if (rList != NULL) rList->clear();
    
    // crossing-decreasing R_I
    if (RT & REMOVE_R_I) {
        for (r=0;r<K->n_reg;r++)
            if ((count_bits(K->regions[r]) == 1) && (r != K->reg_0) && (r != K->reg_1)) {
                if (rList != NULL) rList->push_back(CRSite(REMOVE_R_I, r)); else return true;
            }
    } // remove R I

    // crossing increasing R_I
    if (RT & CREATE_R_I) {
        for (r=0;r<K->n_reg;r++)
            for (a=0;a<K->n2;a++)
                if (BIT(K->regions[r],a)) {
                    if (rList != NULL) {
                        rList->push_back(CRSite(CREATE_R_I,r,a,0));
                        rList->push_back(CRSite(CREATE_R_I,r,a,1));
                    } else return true;
                } 
    } // create R I
    
    // crossing decreasing R-II
    if (RT & REMOVE_R_II) {
        // if (R_DEBUG) cout << "r2" << endl;
        for (r=0;r<K->n_reg;r++)
            if ((count_bits(K->regions[r]) == 2) && (r != K->reg_0) && (r != K->reg_1) &&
                (K->gw(firstBitSet(K->regions[r])) * K->gw(firstBitSet(K->regions[r])+1) > 0 )) { // both positive or both negative
                    if (rList != NULL) rList->push_back(CRSite(REMOVE_R_II, r)); else return true;
                }
    } // remove R II
    
    // crossing-increasing R-II
    if (RT & CREATE_R_II) {
        for (r=0;r<K->n_reg;r++)
            for (a=0;a<K->n2;a++)
                for (b=0;b<K->n2;b++)
                    if ((BIT(K->regions[r],a)) && (BIT(K->regions[r],b))) { // two arcs in same region
                        if (rList != NULL) {
                            if (a != b) { // dont perform RII on self-arc, this is taken care of by RI
                                if ((r != K->reg_0) && (r != K->reg_1)) // RII not in dotted region
                                    rList->push_back(CRSite(CREATE_R_II,r,a,b));
                                else {
                                    rList->push_back(CRSite(CREATE_R_II,r,a,b,0));
                                    rList->push_back(CRSite(CREATE_R_II,r,a,b,1));
                                }
                            }
                            
                            
                        } else return true;
                    }
    } // create R II
    
    // R-III
    if (RT & MODIFY_R_III) {
        for (r=0;r<K->n_reg;r++)
            if ((count_bits(K->regions[r]) == 3) && (r != K->reg_0) && (r != K->reg_1)) {
                if (find_RIII_arcs(K, r, &d[0], &d[1], &d[2])) {
                    if (rList != NULL) rList->push_back(CRSite(MODIFY_R_III, r, d[0],d[1],d[2]));
                    else return true;
                }
            }
    } // modify R III
    
    
    // flype
    if (RT & MODIFY_FLYPE) {
        for (c=1; c <= K->n; c++) { // all crossings
            
            a = K->gw_pos(c); b = K->gw_pos(-c); // positions of crossing c in Gauss word
            SORT(&a,&b);
            if (b-a > K->n) { swap(&a,&b); b += K->n2; }
            
            if ((b-a == K->n) || (b-a <= 2)) continue; // trivial flype
            
            cr_or = cr_xor = (u_region)0; // get crossigns between a and b
            a_ = 0xFFF; b_ = 0; // min and max to search for 2nd block
            for (i=a+1;i<b;i++) { // find 1st block of crossings
                cr_or |= ((u_region)1<<ABS(K->gw(i)));
                cr_xor ^= ((u_region)1<<ABS(K->gw(i)));
            }

            // find 2nd corresponding block            
            for (i=a+1;i<b;i++) // find where the block can lie on
                if (BIT(cr_xor,ABS(K->gw(i)))) {
                    j = K->gw_pos(-K->gw(i));
                    if (j < a) j += K->n2;
                    a_ = MINIMUM(a_,j);
                    b_ = MAXIMUM(b_,j);
                }

            cr_or_ = cr_xor_ = (u_region)0;
            for (i=a_; i<=b_; i++) { // find corresponding 2nd block
                cr_or_ |= ((u_region)1<<ABS(K->gw(i)));
                cr_xor_ ^= ((u_region)1<<ABS(K->gw(i)));
            }
            
            if (a_ == 0xFFF) continue;

            if ((cr_xor ^ cr_xor_)) { // if there are extra crossing, search for them left and right
                i = a_-1;
                while ((MOD(i-a,K->n2)) && (MOD(i-b,K->n2)) && BIT(cr_xor_,ABS(K->gw(i)))) {
                    cr_or_ |= ((u_region)1<<ABS(K->gw(i)));
                    cr_xor_ ^= ((u_region)1<<ABS(K->gw(i)));
                    i--; a_--;
                }
                i = b_+1;
                while ((MOD(i-a,K->n2)) && (MOD(i-b,K->n2)) && BIT(cr_xor_,ABS(K->gw(i)))) {
                    cr_or_ |= ((u_region)1<<ABS(K->gw(i)));
                    cr_xor_ ^= ((u_region)1<<ABS(K->gw(i)));
                    i++; b_++;
                }
                
            }
            // get region A & B
            orient = K->regionGet(MOD(b-1,K->n2),RIGHT) == K->regionGet(MOD(a-1,K->n2),LEFT); // how are the flype incomming arcs oriented?
            reg_A = K->regionGet(MOD(a-1,K->n2), orient);
            reg_B = K->regionGet(MOD(b,K->n2), orient);
            
            if (!(cr_xor ^ cr_xor_) && (count_bits(cr_or_ | cr_or) < K->n-1)) {
                
                regs = (u_region)0; // arcs inside the flype
                regs_ = (u_region)0;
      
                for (i=a; i<b; i++) regs |= ((u_region)1 << MOD(i,K->n2));
                for (i=a_; i<b_; i++) regs_ |= ((u_region)1 << MOD(i,K->n2));
                
                if (regs & regs_) continue;
                
                reg_0_inside_flype = SUBSET(K->regions[K->reg_0],(regs | regs_));
                reg_1_inside_flype = SUBSET(K->regions[K->reg_1],(regs | regs_));
                
                if ((!reg_0_inside_flype) & (!reg_1_inside_flype)) { // can perform flype
                    if (rList != NULL)
                        rList->push_back(CRSite(MODIFY_FLYPE, MOD(a,K->n2),MOD(b,K->n2),MOD(a_,K->n2),MOD(b_,K->n2),0));
                    else return true;
                } else
                if ((reg_0_inside_flype & reg_1_inside_flype) ||
                    (reg_0_inside_flype & ((K->reg_1 == reg_A) || (K->reg_1 == reg_B))) ||
                    (reg_1_inside_flype & ((K->reg_0 == reg_A) || (K->reg_0 == reg_B)))) { // can perform flype in the complement
                    //cout << "FLYPE: " <<(int)MOD(a,K->n2)<<(int)MOD(b,K->n2)<<(int)MOD(a_,K->n2)<<(int)MOD(b_,K->n2)<<endl;
                    if (rList != NULL) 
                        rList->push_back(CRSite(MODIFY_FLYPE, MOD(a,K->n2),MOD(b,K->n2),MOD(a_,K->n2),MOD(b_,K->n2),1));
                    else return true;
                } 
            } else {
               // check, return error
            }
        }
    } // modify FLYPE
    
    if (rList == NULL) return false;
    
    return (rList->size() != 0);
}

// TODO: optimize for oriented projective trivial knots

// removes crossing (either by R-I or R-II), if K is the an obvious unknot
// we use a seperate function since we define the projective trivial knot by the indicies of reg_0 and reg_1
// TODO: check if this function can be avoided
bool remove_unknot(Cknot *K, int omega) {
    
    int g = 0;
    bool r0_a, r0_b, r1_a, r1_b, b_s;
    
    if (K->n != omega) return false;
    
    K->canonical(ORIENTED);
    
    if (K->n == 1) { // remove by R-I
        
        if (K->reg_0 == K->reg_1) { K->zerofyKnot(); K->n_reg = 2; K->reg_0 = 0; K->reg_1 = 0; return true;} // affine unknot
        
        // TODO: simplify with XOR
        
        b_s = BIT(K->signs,1);
        
        r0_a = (K->regions[K->reg_0] == (u_region)(1<<0));
        r0_b = (K->regions[K->reg_0] == (u_region)(1<<1));
        r1_a = (K->regions[K->reg_1] == (u_region)(1<<0));
        r1_b = (K->regions[K->reg_1] == (u_region)(1<<1));
        
        if ((b_s == PLUS) && (!r0_a) && (r0_b) && (!r1_a) && (!r1_b)) g = 1; else
        if ((b_s == PLUS) && (!r0_a) && (!r0_b) && (r1_a) && (!r1_b)) g = 1; else
        if ((b_s ==MINUS) && (r0_a) && (!r0_b) && (!r1_a) && (!r1_b)) g = 1; else
        if ((b_s ==MINUS) && (!r0_a) && (!r0_b) && (!r1_a) && (r1_b)) g = 1; else
        
        if ((b_s == PLUS) && (r0_a) && (!r0_b) && (!r1_a) && (!r1_b)) g = -1; else
        if ((b_s == PLUS) && (!r0_a) && (!r0_b) && (!r1_a) && (r1_b)) g = -1; else
        if ((b_s ==MINUS) && (!r0_a) && (r0_b) && (!r1_a) && (!r1_b)) g = -1; else
        if ((b_s ==MINUS) && (!r0_a) && (!r0_b) && (r1_a) && (!r1_b)) g = -1;
        
    } else
    if (K->n == 2) { // remove by R-II
        
        if ((K->reg_0 == K->reg_1) || ((count_bits(K->regions[K->reg_0]) == 1) && (count_bits(K->regions[K->reg_1]) == 1)))
        { K->zerofyKnot(); K->n_reg = 2; K->reg_0 = 0; K->reg_1 = 0; return true;} // affine unknot
        
        // TODO: simplify with XOR
        
        b_s = BIT(K->signs,1);
        
        r0_a = (K->regions[K->reg_0] == (u_region)(1<<0));
        r0_b = (K->regions[K->reg_0] == (u_region)(1<<2));
        r1_a = (K->regions[K->reg_1] == (u_region)(1<<0));
        r1_b = (K->regions[K->reg_1] == (u_region)(1<<2));
        
        if ((b_s == PLUS) && (!r0_a) && (!r0_b) && (r1_a) && (!r1_b)) g = 1; else
        if ((b_s == PLUS) && (!r0_a) && (!r0_b) && (!r1_a) && (r1_b)) g = 1; else
        if ((b_s ==MINUS) && (r0_a) && (!r0_b) && (!r1_a) && (!r1_b)) g = 1; else
        if ((b_s ==MINUS) && (!r0_a) && (r0_b) && (!r1_a) && (!r1_b)) g = 1; else
                        
        if ((b_s == PLUS) && (r0_a) && (!r0_b) && (!r1_a) && (!r1_b)) g = -1; else
        if ((b_s == PLUS) && (!r0_a) && (r0_b) && (!r1_a) && (!r1_b)) g = -1; else
        if ((b_s ==MINUS) && (!r0_a) && (!r0_b) && (r1_a) && (!r1_b)) g = -1; else
        if ((b_s ==MINUS) && (!r0_a) && (!r0_b) && (!r1_a) && (r1_b)) g = -1;
        
    }
    
    if (g == 0) throw 20;
    
    K->zerofyKnot();
    K->n_reg = 2;
    K->reg_0 = ((g == 1) ? 0 : 1); // (0,1) - t_1
    K->reg_1 = ((g == 1) ? 1 : 0); // (1,0) - t_-1
    
    return true;
}

// Performing R-moves

// Crossing decreasing R-I
void removeRI(Cknot *K, u_number reg) { /* OPTIMIZED */
    
    if (remove_unknot(K,1)) return;
    u_number r, arc;
    s_letter x = ABS(K->gw_[arc = firstBitSet(K->regions[reg])]); // get crossing
    deleteRegion(K,reg); // delete from regions
    for (r=0; r<K->n_reg; r++) DELETE_BITS(&(K->regions[r]), arc,2); // delete arc from regions
    if (arc == K->n2-1) for (r=0; r<K->n_reg; r++) K->regions[r] >>= 1; // if last arc deleted, remove 0th arc
    deleteCrossingFromGaussWord(K,x);
}

// Crossing decreasing R-II (link)
void removeRI(Clink *L, u8 c, u_number reg) { /* OPTIMIZED */

    u_number r, arc;
    s_letter x = ABS(L->gw_[c][arc = firstBitSet(L->regions[reg])]); // get crossing
    
    deleteRegion(L,reg); // delete from regions
    
    for (r=0; r<L->n_reg; r++) DELETE_BITS(&(L->regions[r]), arc,2); // delete arc from regions
    if (arc == L->n*2-1) for (r=0; r<L->n_reg; r++) L->regions[r] >>= 1; // if last arc deleted, remove 0th arc
    
    deleteCrossingFromGaussWord(L,c,x);
}

// crossing increasing R-I
void createRI(Cknot *K, u_number reg, u_number arc, int sgn) { // sgn = +1, -1  /* OPTIMIZED */
    
    u_number r, reg_ = K->adjacentRegion(reg, arc);
    int first_sign;
    
    first_sign = B2SIGN(K->regionSide(reg, arc)^SIGN2B(sgn)); // what sign should the first letter in the GW be?
    
    for (r=0;r<K->n_reg;r++) INSERT_BITS(&K->regions[r], arc, 2, (u_region)0); // shift arcs
    
    K->regions[reg] |= ((u_region)0x3 << arc); // insert created arcs
    K->regions[reg_] |= ((u_region)0x1 << arc);
    K->regions[K->n_reg++] = (u_region)1 << (arc+1); // add new region
    
    insertCrossingToGaussWord(K, arc+1, arc+1, first_sign*(s_letter)(K->n+1), -1*first_sign*(s_letter)(K->n+1));

    K->signs |= SIGN2B(sgn) << K->n; // insert sign
}

// crossing decreasing R-II
// TODO: optimize
void removeRII(Cknot *K, u_number reg) {  /* NOT OPTIMIZED */
    
    if (remove_unknot(K,2)) return;
    
    int r,g1,g2,g3,g4,r1,r2;
    
    u_number arca = firstBitSet(K->regions[reg]); // find the two arcs
    u_number arcb = secondBitSet(K->regions[reg]);
    
    s_letter x = ABS(K->gw(arca)); // 1st crossing
    s_letter y = ABS(K->gw(arca+1)); // 2nd crossing
    if (x > y) { s_letter z=x;x=y;y=z; } // sort them
    
    if ((ABS(K->gw_[arca])) == ABS(K->gw_[arcb])) { // same if the arc some in from the same directions
        g1 = MOD(arca-1,K->n2); g2 = MOD(arcb-1,K->n2); g3 = MOD(arca+1,K->n2); g4 = MOD(arcb+1,K->n2); }
    else {
        g1 = MOD(arca-1,K->n2); g4 = MOD(arcb-1,K->n2); g3 = MOD(arca+1,K->n2); g2 = MOD(arcb+1,K->n2); }

    // find two regions to join
    r1 = r2 = -1;
    for (r=0;r<K->n_reg;r++) {
        if ((BIT(K->regions[r],g1)) && (BIT(K->regions[r],g2)) && (!BIT(K->regions[r],arca)) && (!BIT(K->regions[r],arcb)) && (r1 == -1)) r1 = r;
        if ((BIT(K->regions[r],g3)) && (BIT(K->regions[r],g4)) && (!BIT(K->regions[r],arca)) && (!BIT(K->regions[r],arcb)) && (r2 == -1)) r2 = r;
    }

    if (r1+r2 != -2) { // join the two regions, if not joined already
        K->regions[r1] |= K->regions[r2]; // union of arcs
        // delete the region
        if (K->reg_0 == r2) K->reg_0 = r1;
        if (K->reg_1 == r2) K->reg_1 = r1;
        deleteRegion(K,r2);
        if (r2 < reg) reg--;
    }
    
    deleteRegion(K,reg);
    deleteArcsFromRegions(K,arca,arcb);
   
    // change the gauss word
    if (x > y) { deleteCrossingFromGaussWord(K,x); deleteCrossingFromGaussWord(K,y); }
    else { deleteCrossingFromGaussWord(K,y); deleteCrossingFromGaussWord(K,x); }
}

// crossing increasing R-II
// deletes the region that splits, and puts the new two ones at the end of the region list
// if the split region contains a dot, we can manualy put the dot into the last or the second last reigon
void createRII(Cknot *K, u_number reg, u_number over_arc, u_number under_arc, u_number dotted_select) {  /* OPTIMIZED */
    u_number r;
    u_number min_arc = MINIMUM(over_arc,under_arc), max_arc = MAXIMUM(over_arc,under_arc);
    s_small over_o = K->regionSide(reg, over_arc), under_o = K->regionSide(reg, under_arc);
    s_small max_o = (max_arc == over_arc ? over_o : under_o), min_o = (min_arc == over_arc ? over_o : under_o);
    u_number r_adj_max = K->adjacentRegion(reg, max_arc), r_adj_min = K->adjacentRegion(reg, min_arc);
    s_small same_direction = (under_o != over_o);
    
    // expand regions
    for (r=0;r<K->n_reg;r++) {
        INSERT_BITS(&K->regions[r], max_arc, 2, (u_region)0); // shift arcs
        INSERT_BITS(&K->regions[r], min_arc, 2, (u_region)0);
    }
    
    // gauss word
    int max_sign = B2SIGN(max_arc != over_arc); // sign of gauss letter of maximal arc
    insertCrossingToGaussWord(K, max_arc+1, max_arc+1, max_sign*(s_letter)(K->n+1), max_sign*(s_letter)(K->n+2));
    insertCrossingToGaussWord(K, min_arc+1, min_arc+1, -max_sign*(s_letter)(K->n+(same_direction?0:1)), -max_sign*(s_letter)(K->n+(same_direction?1:0)));
    K->signs |= ((u_sign)((max_o^same_direction^(max_arc==over_arc)) ? 1 : 2)) << (K->n-1);

    // fix the regions
    K->regions[r_adj_max] |= ((u_region)1 << (min_arc+1)) | ((u_region)1 << (max_arc+2));
    K->regions[r_adj_min] |= ((u_region)1 << (min_arc+0)) | ((u_region)1 << (max_arc+3));
    
    if (reg == K->reg_0) K->reg_0 = K->n_reg + dotted_select + 1;
    if (reg == K->reg_1) K->reg_1 = K->n_reg + dotted_select + 1;
    
    deleteRegion(K, reg); // delete split region
    K->regions[K->n_reg++] = ((u_region)1<< (min_arc+1)) | ((u_region)1 << (max_arc+3)); // add newly created region
    K->regions[K->n_reg++] = K->gen_reg(min_arc,min_o); // split region 1
    K->regions[K->n_reg++] = K->gen_reg((same_direction ? max_arc+4 : max_arc+2),max_o); // split region 2

}

// functions used in R-III
// TODO: optimize search for adjacent arcs

bool crossing_at_region(Cknot *K, u_number reg, s_letter x) {
    return (BIT(K->regions[reg],K->gw_pos(x)) || BIT(K->regions[reg],MOD(K->gw_pos(x)-1,K->n2)));
}

s_letter common_crossing(Cknot *K, u_number arc_a, u_number arc_b) {
    if (K->gw(arc_a) == -K->gw(arc_b)) return ABS(K->gw(arc_a));
    if (K->gw(arc_a) == -K->gw(arc_b+1)) return ABS(K->gw(arc_a));
    if (K->gw(arc_a+1) == -K->gw(arc_b)) return ABS(K->gw(arc_a+1));
    if (K->gw(arc_a+1) == -K->gw(arc_b+1)) return ABS(K->gw(arc_a+1));
    return 0xFF;
}

// get region of opposite corner of crossing at arc_a & arc_b
u_number antipodal_region(Cknot *K, u_number reg, u_number arc_a, u_number arc_b) {
    u_number r[4] = {
        K->regionGet(MOD(arc_a+1,K->n2), !K->regionSide(reg, arc_a)),
        K->regionGet(MOD(arc_a-1,K->n2), !K->regionSide(reg, arc_a)),
        K->regionGet(MOD(arc_b+1,K->n2), !K->regionSide(reg, arc_b)),
        K->regionGet(MOD(arc_b-1,K->n2), !K->regionSide(reg, arc_b))};
    
    if (!crossing_at_region(K, r[0], common_crossing(K,arc_a,arc_b))) r[0] = 0x7F;
    if (!crossing_at_region(K, r[1], common_crossing(K,arc_a,arc_b))) r[1] = 0x7E;
    if (!crossing_at_region(K, r[2], common_crossing(K,arc_a,arc_b))) r[2] = 0x7D;
    if (!crossing_at_region(K, r[3], common_crossing(K,arc_a,arc_b))) r[3] = 0x7C;
    
    return common2(r[0],r[1],r[2],r[3]);
}

// Apply the R-III move
void modifyRIII(Cknot *K, u_number reg, u_number over_arc = 0, u_number under_arc = 0, u_number mixed_arc = 0 ) {
    
    if (!(over_arc | under_arc | mixed_arc)) // if arcs not provided, find them
        find_RIII_arcs(K, reg, &over_arc, &under_arc, &mixed_arc);
    
    u_number over_outside = K->adjacentRegion(reg,over_arc);
    u_number under_outside = K->adjacentRegion(reg,under_arc);
    u_number mixed_outside = K->adjacentRegion(reg,mixed_arc);
    
    u_number anti_over, anti_under, anti_mixed; // regions on the other side of the arc
    
    anti_over = antipodal_region(K, reg, under_arc, mixed_arc);
    anti_under = antipodal_region(K, reg, over_arc, mixed_arc);
    anti_mixed = antipodal_region(K, reg, over_arc, under_arc);
    
    K->regions[over_outside] &= ~((u_region)1<<over_arc); // delete arcs
    K->regions[under_outside] &= ~((u_region)1<<under_arc);
    K->regions[mixed_outside] &= ~((u_region)1<<mixed_arc);

    K->regions[anti_over] |= ((u_region)1<<over_arc); // add arcs
    K->regions[anti_under] |= ((u_region)1<<under_arc);
    K->regions[anti_mixed] |= ((u_region)1<<mixed_arc);
    
    //change gauss word
    swap( &(K->gw_[over_arc]), &(K->gw_[(over_arc+1)%K->n2]) );
    swap( &(K->gw_[under_arc]), &(K->gw_[(under_arc+1)%K->n2]) );
    swap( &(K->gw_[mixed_arc]), &(K->gw_[(mixed_arc+1)%K->n2]) );
}

// Apply the flype
void flype(Cknot *K, u_number a, u_number b, u_number a_, u_number b_, u_number complement) {
    
    u_number arc_in, arc_out, ar_in, ar_out; // dotted arcs
    u_number r, reg_O, reg_X, reg_A, reg_B, reg_Y;
    s_letter i, orient;
    u_region TMP, INNER_ARCS = (u_region)0, REG_A = (u_region)0, REG_B = (u_region)0, REG_O;
    arc_in = MOD(a-1,K->n2); // get arcs involved in the flype crossing
    arc_out = MOD(b,K->n2);
    ar_in = MOD(a,K->n2);
    ar_out = MOD(b-1,K->n2);
    
    orient = (K->regionGet(ar_out,RIGHT) == K->regionGet(arc_in,LEFT)); // how are the flype incomming arcs oriented?
    
    reg_O = K->regionGet(ar_out, orient); // get regions involved on the edges of the flype
    reg_X = K->regionGet(arc_in, !orient);
    reg_A = K->regionGet(arc_in, orient);
    reg_B = K->regionGet(arc_out, orient);
    reg_Y = 000;

    /* // used for debugging
    bool er = true;
    
    if (reg_O >= K->n_reg) cout << "ERRRRR 0" << endl; else
    if (reg_X >= K->n_reg) cout << "ERRRRR x" << endl; else
    if (reg_A >= K->n_reg) cout << "ERRRRR a" << endl; else
    if (reg_B >= K->n_reg) cout << "ERRRRR b" << endl; else
    if (a_ >= K->n2) cout << "ERRRRR a_" << endl; else
    if (b_ >= K->n2) cout << "ERRRRR b_" << endl; else
        
    er= false;
    
    if (er) {
        cout << "n = " << (int)K->n << " n2= " << (int)K->n2 << endl;
        cout << "arc_in "<<(int)arc_in << ", arc_out " << (int)arc_out << ", ar_in " << (int)ar_in << ", ar_out " << (int)ar_out << endl;
        cout << "reg O " << (int)reg_O << ", reg_X " << (int)reg_X << ", reg_A " << (int)reg_A << ", reg_B " << (int)reg_B << " ,reg_Y " << (int)reg_Y << endl;
        
        cout << endl;
        
    }
    cout << "reg O " << (int)reg_O << ", reg_X " << (int)reg_X << ", reg_A " << (int)reg_A << ", reg_B " << (int)reg_B << " ,reg_Y " << (int)reg_Y << endl;
     */
    
    INNER_ARCS = (u_region)0; // get inner arcs
    for (i=MOD(a+1,K->n2); i!=MOD(b-1,K->n2); i = MOD(i+1,K->n2)) INNER_ARCS |= ((u_region)1 << i);
    for (i=a_; i!=b_; i = MOD(i+1,K->n2)) INNER_ARCS |= ((u_region)1 << i);
    
    // get Y arcs info
    u_number arc_in_ = MOD(a_-1,K->n2);
    u_number arc_out_ = b_;
    u_number arc_A_ = BIT(K->regions[reg_A],arc_in_) ? arc_in_ : arc_out_;
    u_number arc_B_ = BIT(K->regions[reg_B],arc_in_) ? arc_in_ : arc_out_;
    s_small turn_around = (s_small)(int)BIT(K->regions[reg_A],arc_in_); // optimized if put to previous 2 lines
    
    reg_Y = K->adjacentRegion(reg_A, arc_A_);
    
    // modify new regions
    REG_A = (K->regions[reg_A] & ~(INNER_ARCS | ((u_region)1<<ar_in) | ((u_region)1<<ar_out))) | (K->regions[reg_B] & INNER_ARCS);
    REG_B = (K->regions[reg_B] & ~(INNER_ARCS | ((u_region)1<<ar_in) | ((u_region)1<<ar_out))) | (K->regions[reg_A] & INNER_ARCS);
    
    REG_O = K->regions[reg_O];
    
    K->regions[reg_O] = K->regions[reg_Y] & INNER_ARCS; // replace REG_Y1 with REG_O
    K->regions[reg_Y] &= ~(INNER_ARCS);
    K->regions[reg_X] |= REG_O & INNER_ARCS;
    
    
    K->regions[reg_A] = REG_A;
    K->regions[reg_B] = REG_B;
    
    // renumerate all regions
    for (r=0;r<K->n_reg;r++) {
        TMP = (u_region)0;
        for (i=0;i<K->n2;i++)
            if ((BIT(K->regions[r],i)) && (i != a) && (i != MOD(b-1,K->n2)))
                TMP |= (u_region)1 << MOD( i - (i>a?1:0) - (i>=b?1:0) + (i>=a_?1:0) + (i>=b_?1:0), K->n2);
        K->regions[r] = TMP;
    }
    
    // add newly created arcs
    
    arc_A_ = MOD(arc_A_ - (arc_A_>a?1:0) - (arc_A_>=b?1:0) + (arc_A_>=a_?1:0) + (arc_A_>=b_?1:0),K->n2);
    arc_B_ = MOD(arc_B_ - (arc_B_>a?1:0) - (arc_B_>=b?1:0) + (arc_B_>=a_?1:0) + (arc_B_>=b_?1:0),K->n2);
    
    if (turn_around) {
        K->regions[reg_A] |= (u_region)1 << MOD(arc_B_-1,K->n2);
        K->regions[reg_B] |= (u_region)1 << MOD(arc_A_+1, K->n2);
        K->regions[reg_O] |= ((u_region)1 << MOD(arc_B_-1,K->n2)) | ((u_region)1 << MOD(arc_A_+1,K->n2));
    } else {
        K->regions[reg_A] |= (u_region) 1 << MOD(arc_B_+1,K->n2);
        K->regions[reg_B] |= (u_region)1 << MOD(arc_A_-1, K->n2);
        K->regions[reg_O] |= ((u_region)1 << MOD(arc_B_+1,K->n2)) | ((u_region)1 << MOD(arc_A_-1,K->n2));
    }
    
    
    // Gauss word
    
    // mirror the flyped part
    for (i=a+1; MOD(b-i,K->n2); i++) K->gw_[MOD(i,K->n2)] *= -1;
    for (i=a_; MOD(b_+1-i,K->n2); i++) K->gw_[MOD(i,K->n2)] *= -1;
    swapCrossingInGaussWord(K, a, b, a_, MOD(b_+1,K->n2),turn_around);
    if (complement) { // if flype performed "on the outside of the loop"
        for (i=0; i<K->n2; i++) K->gw_[i] *= -1;
        
        bool A0, A1, B0, B1, Y0, Y1, O0, O1;
        
        A0 = (reg_A == K->reg_0);
        A1 = (reg_A == K->reg_1);
        B0 = (reg_B == K->reg_0);
        B1 = (reg_B == K->reg_1);
        Y0 = (reg_Y == K->reg_0);
        Y1 = (reg_Y == K->reg_1);
        O0 = (reg_O == K->reg_0);
        O1 = (reg_O == K->reg_1);
        
        if (A0) K->reg_0 = reg_B;
        if (A1) K->reg_1 = reg_B;
        if (B0) K->reg_0 = reg_A;
        if (B1) K->reg_1 = reg_A;
        if (Y0) K->reg_0 = reg_O;
        if (Y1) K->reg_1 = reg_O;
        if (O0) K->reg_0 = reg_X;
        if (O1) K->reg_1 = reg_X;
    }
}

// functions used by R-moves

// delete reg from knot K
void deleteRegion(Cknot *K, u_number reg) {
    for (int i=reg+1;i<K->n_reg;i++) K->regions[i-1] = K->regions[i];
    if (K->reg_1 > reg) K->reg_1--;
    if (K->reg_0 > reg) K->reg_0--;
    K->n_reg--;
}

// delete reg from link K
void deleteRegion(Clink *L, u_number reg) {
    for (int i=reg+1;i<L->n_reg;i++) L->regions[i-1] = L->regions[i];
    if (L->reg_1 > reg) L->reg_1--;
    if (L->reg_0 > reg) L->reg_0--;
    L->n_reg--;
}

// delete crossing x in K
void deleteCrossingFromGaussWord(Cknot *K, s_letter x) {
    int z,i,j,m = K->n2;
    // change the Gauss word 
    DELETE_BIT(&K->signs,x);
    for (i=0, m = K->n2; i<m;i++)
        if (x == ABS(K->gw_[i])) { for (j=i+1;j<m;j++) K->gw_[j-1] = K->gw_[j]; i--; m--; }
    K->n--; K->n2 -= 2;
    // renumerate gauss word
    for (i=0;i<K->n2;i++)
        if (ABS(z=K->gw_[i]) > x) K->gw_[i] = (ABS(z)-1) * SIGN(z);
}

// delete crossing x from component c in link L
void deleteCrossingFromGaussWord(Clink *L, u8 c, s_letter x) {
    // change the gauss word
    int z,i,j,m = L->n_arcs[c];
    DELETE_BIT(&L->signs,x);
    for (i=0, m = L->n_arcs[c]; i<m;i++)
        if (x == ABS(L->gw_[c][i])) { for (j=i+1;j<m;j++) L->gw_[c][j-1] = L->gw_[c][j]; i--; m--; }
    L->n--; L->n_arcs[c] -= 2;
    // renumerate gauss word
    for (c=0;c<L->n_components;c++)
        for (i=0;i<L->n_arcs[c];i++)
        if (ABS(z=L->gw_[c][i]) > x) L->gw_[c][i] = (ABS(z)-1) * SIGN(z);
}


// insert the crossing x to the position pos1 and pos2 in the Gauss code
// if pos1 = pos2, inserts x1x2
void insertCrossingToGaussWord(Cknot *K, u_number pos1, u_number pos2, s_letter x1, s_letter x2) { 
    int i;
    if (pos2 >= pos1) for (i=K->n2-1;i>=pos1;i--) K->gw_[i + (i<pos2 ? 1 : 2) ] = K->gw_[i];
    else for (i=K->n2-1;i>=pos2;i--) K->gw_[i + (i<pos1 ? 1 : 2) ] = K->gw_[i];
    if (pos1 == pos2) { K->gw_[pos1]  = x1; K->gw_[pos1+1] = x2; }
    else { K->gw_[pos1]  = x1; K->gw_[pos2] = x2; }
    K->n++;
    K->n2+=2;
}

// deletes arc1 and arc2 from all regions of K
// TODO: use only one argument arc
void deleteArcsFromRegions(Cknot *K, u_number arc1, u_number arc2) {
    int r, i, c;
    // delete the arcs
    for (r=0;r<K->n_reg;r++) K->regions[r] &= ~(((u_region)1<<arc1) | ((u_region)1<<MOD(arc1+1,K->n2)) |
                                                ((u_region)1<<arc2) | ((u_region)1<<MOD(arc2+1,K->n2)));
    // renumerate the arcs in regions
    int m = K->n2;
    for (i=0;i<m;i++) {
        c = 0;
        for (r=0;r<K->n_reg;r++) if (BIT(K->regions[r],i)) c++;
        if (!c) {
            for (r=0;r<K->n_reg;r++)
                DELETE_BIT(&K->regions[r], i);
            i--; m--;
        }
    }
}

// swap two crossings in K: pos_a <-> pos_a_ and pos_b <-> pos_b_ (first swap a, then b; or vice-versa if turn_around = true)
void swapCrossingInGaussWord(Cknot *K, u_number pos_a, u_number pos_b, u_number pos_a_, u_number pos_b_, s_small turn_around) {
    int i;
    int pos_a_t = pos_a, pos_b_t = pos_b;
    if (pos_b_ == 0) pos_b_ = K->n2;
    
    s_letter a = K->gw_[pos_a];
    s_letter b = K->gw_[pos_b];
    if (turn_around) swap(a,b);
    
    if (pos_a <= pos_b) pos_b--; else pos_a--;
    
    for (i=MINIMUM(pos_a,pos_b);i<K->n2-2;i++)
        K->gw_[i] = K->gw( i + ( (i>=pos_a)?1:0 ) + ( (i>=pos_b)?1:0 ) ); // remove crossing
    
    pos_a_ -= ((pos_a_>pos_a_t)?1:0 ) + ( (pos_a_>pos_b_t)?1:0 ); // TODO: check if equaliy 
    pos_b_ -= ((pos_b_>pos_a_t)?1:0 ) + ( (pos_b_>pos_b_t)?1:0 );
    
    if (pos_a_ <= pos_b_) pos_b_++; else pos_a_++;
    
    for (i=K->n2-1; i > MINIMUM(pos_a_,pos_b_); i--)
        K->gw_[i] = K->gw( i - ((i>=pos_a_)?1:0) - ((i>=pos_b_)?1:0) );
    
    K->gw_[pos_a_] = a;
    K->gw_[pos_b_] = b;
    
}

// returns the region quadrants of crossing c in K, enumerated as:
//  ^
// 3|0
// -|->
// 2|1
void quadrants(Cknot *K, s_letter c, u_number *r0, u_number *r1, u_number *r2, u_number *r3 ) {
    u_region r0_ = (u_region)0, r1_ = (u_region)0, r2_ = (u_region)0, r3_ = (u_region)0;
    
    r0_ = K->gen_reg(K->gw_pos(c), 0);
    r1_ = K->gen_reg(K->gw_pos(-c), 0);
    r2_ = K->gen_reg(MOD(K->gw_pos(c)-1,K->n2), 1);
    r3_ = K->gen_reg(MOD(K->gw_pos(-c)-1,K->n2), 1);
    
    for (u_number r = 0; r < K->n_reg; r++)
        if (K->regions[r] == r0_) *r0 = r; else
            if (K->regions[r] == r1_) *r1 = r; else
                if (K->regions[r] == r2_) *r2 = r; else
                    if (K->regions[r] == r3_) *r3 = r;
    
    if (!BIT(K->signs,c)) { swap(r0,r1); swap(r2,r3); }
}

// performs a R-move on K given by rSite
void rMove(Cknot *K, CRSite * rSite) {

    switch(rSite->type) {
            
        case MODIFY_R_III:
            modifyRIII(K,rSite->data[0],rSite->data[1],rSite->data[2],rSite->data[3]);
            break;
            
        case MODIFY_FLYPE:
            flype(K,rSite->data[0],rSite->data[1],rSite->data[2],rSite->data[3],rSite->data[4]);
            break;
            
        case REMOVE_R_I:
            removeRI(K, rSite->data[0]);
            break;
            
        case REMOVE_R_II:
            removeRII(K,rSite->data[0]);
            break;
            
        case CREATE_R_I:
            createRI(K, rSite->data[0],rSite->data[1],rSite->data[2]);
            break;
            
        case CREATE_R_II:
            createRII(K, rSite->data[0],rSite->data[1],rSite->data[2],rSite->data[3]);
            break;
    }
}

#endif
