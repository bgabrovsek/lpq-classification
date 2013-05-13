//  algorithms.h
//  Created by Boštjan on 4/22/13.
//
//  Various algorithms: partitioning knots by invariant, reducing knots, classification algorithms, ... 

#ifndef Lpq_algorithms_h
#define Lpq_algorithms_h

#include <list>
#include <vector>
#include <set>
#include "global_vars.h"
#include "common.h"
#include "knot.h"
#include "reidemeister_moves.h"
#include "HOMFLY.h"
#include "multivariate_laurent.h"
#include "KBSM.h"
#include "classification_data.h"
#include "in_out.h"
#include "lens.h"


/*void log_(string s, int i = -1, int j = -1, int k = -1) {
    
    ofstream file (LOG_FILE, ios::app);
    if (file.is_open()) {
        file << s << ": ";
        if (i!=-1) file << i << " ";
        if (j!=-1) file << j << " ";
        if (k!=-1) file << k << " ";
        file << endl;
        file.close();
    }
    
}
*/

// global variables

// store the KBSMs of affine knots up to n = 5
Cmultivariate *prime_KBSM[7]; // KBSM of prime affine knots (for checking conneced sums)
Cmultivariate *prime_KBSM_mod[8]; // modified KBSM of prime affine knots (for checking affinity)

// maximal memory that can be consumed by BST_reducible
u64 constant_memory = sizeof(Cknot)*6000;
u64 memory_factor = sizeof(Cknot) + 32;


// group knots that have the same polynomial
class CknotGroup { 
public:
    s16 group_ID;
    Cmultivariate * poly;
    t_knot_set * KNOTS;
    CknotGroup() { KNOTS = new t_knot_set(); group_ID = -1; }
};

// get knot by its ID
Cknot * get_knot_by_id(int ID) {
    for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++)
        if ((*iK)->ID == ID) return (*iK);
    return NULL;
}


// set first knots in group to prime (torus)
void set_single_knots_to_classified() {
    for (t_group_list::iterator iG = G_HOMFLY->begin(); iG != G_HOMFLY->end(); iG++)
        classification[ (*((*iG)->KNOTS->begin()))->ID ].prime = true;
}

// set first knots in group to prime (L(p,q))
void set_first_knots_to_classified_lpq() {
    for (int l = 0; l < N_HOMEO; l++) /*if (consider_lpq[l])*/ {
        t_group_list * G = GROUP_OF_CHOICE(l)[l]; // select either the group of KBSMs or HSMs
            for (t_group_list::iterator iG = G->begin(); iG != G->end(); iG++)
                classification_lpq[l][ (*((*iG)->KNOTS->begin()))->ID ].prime = true;
    }
    
}

void renumerate_prime_knot_table(t_knot_set *);

// generates knot candidates, reduces those that are simplifiable by removing crossings by R-moves
void generate_knot_candidates(t_knot_set *KNOTS) {
    
    generateAllKnots(1, 5, KNOTS, REMOVE_3_KINKS|REMOVE_AFFINE|REMOVE_KINK_OUTSIDE_DOT|REMOVE_ADJACENT_SELECTED|REMOVE_BST_REMOVABLE|REMOVE_OBVIOUS_COMPOSITES/*|REMOVE_CONNECTED_SUM*/);
    
    renumerate_prime_knot_table(KNOTS);

}

// mark primes as duplicates in the classification table (composites checked by hand)
void mark_composites(t_knot_set *KNOTS) {
    for (t_knot_set::iterator it = KNOTS->begin(); it != KNOTS->end(); it++)
        if (get_knot_property((*it)->deepCopy()) & PROPERTY_SUM) {
            classification[(*it)->ID].prime = NO;
            classification[(*it)->ID].duplicate = YES;
        }
}

// free memory
void deleteKnots(t_knot_set *KNOTS) {
    for (t_knot_set::iterator it = KNOTS->begin(); it != KNOTS->end(); it++) delete *it;
    KNOTS->clear();
}

// free memory
void deleteGroups(t_group_list *G) {
    for (t_group_list::iterator it = G->begin(); it != G->end(); it++) {
        delete (*it)->poly;
        (*it)->KNOTS->clear();
    }
    G->clear();
}

// get index of group from the polynomial of the group (torus)
int group_index_from_poly(t_group_list *KNOT_GROUP, Cmultivariate *poly) {
    poly->simplify();
    for (t_group_list::iterator itg = KNOT_GROUP->begin(); itg != KNOT_GROUP->end(); itg++) {
        if (poly_equal_orientation((*itg)->poly, poly, false, false)) {
            return (*itg)->group_ID;
        }
    }
    throw 20;
}

// get interator of group from the polynomila of the group (torus)
t_group_list::iterator group_from_poly(t_group_list *KNOT_GROUP, Cmultivariate *poly) {
    t_group_list::iterator itg;
    poly->simplify();
    for (itg = KNOT_GROUP->begin(); itg != KNOT_GROUP->end(); itg++) {
        if (poly_equal_orientation((*itg)->poly, poly, false, false)) {
            return itg;
        }
    }
    return KNOT_GROUP->end();
}

// get index of group from the polynomial of the group (L(p,q))
int group_index_from_poly_lpq(int l, t_group_list *KNOT_GROUP, Cmultivariate *poly) {
    t_group_list::iterator itg;
    poly->simplify();
    for (itg = KNOT_GROUP->begin(); itg != KNOT_GROUP->end(); itg++) {
        if (poly_equal_orientation_lpq(l,(*itg)->poly, poly, true, true)) {
            return (*itg)->group_ID;
        }
    }
    throw 20;
}


// calculate candidates of chiral knots, store to array chirality
void get_chirality(t_knot_set *KNOTS, t_group_list *GROUP) {
    
    t_knot_set::iterator iK;
    t_group_list::iterator iG, iGB;
    Cknot *K;
    Cmultivariate *poly;
    int IDA, IDB, up_id;
    
    for (iG = GROUP->begin(); iG != GROUP->end(); iG++) {
        K = (*((*iG)->KNOTS->begin()))->deepCopy();
        K->mirror();
        poly = HOMFLY(K); poly->simplify();
        IDA = (*iG)->group_ID; // knot
        IDB = group_index_from_poly(GROUP, poly); // knot of mirror
        up_id = -1;
        if (IDA <= IDB) {
            for (iK = (*iG)->KNOTS->begin(); iK != (*iG)->KNOTS->end(); iK++)
                if (!classification[(*iK)->ID].duplicate) {
                    if (up_id == -1) up_id =(*iK)->ID;
                    chirality[(*iK)->ID] = ((IDA == IDB)? 3 : 1); 
                }
            iGB = group_from_poly(GROUP, poly);
            
            for (iK = (*iGB)->KNOTS->begin(); iK != (*iGB)->KNOTS->end(); iK++)
                if (!classification[(*iK)->ID].duplicate) {
                    chirality[(*iK)->ID] = ((IDA == IDB)? 3 : 2);
                    chirality_id[(*iK)->ID] = up_id; //(*(*group_from_poly(GROUP, poly))->KNOTS->begin())->ID;

                        
                }
        }
    }
}

// print knots that are mphichiral candidates
void get_amphichiral_candidates(t_knot_set *KNOTS, t_group_list *GROUP) {
    
    t_knot_set::iterator iK;
    t_group_list::iterator iG, iGB;
    Cknot *K;
    Cmultivariate *poly;
    int IDA, IDB;

    for (iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
        
        if (!classification[(*iK)->ID].duplicate) { // not a duplicate
            K = (*iK)->deepCopy();
            K->mirror();
            poly = HOMFLY(K); poly->simplify();
            if ( poly_equal_orientation((HOMFLYS)[(*iK)->ID], poly) ) {
                if (!classification[(*iK)->ID].prime) cout << "?"; 
                cout << (*iK)->ID << ", " << flush;
            }
            
        }
        
    }
    
    cout << endl << endl << "Unclassified: " << endl;
    for (iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
        if ((!classification[(*iK)->ID].duplicate) && (!classification[(*iK)->ID].prime))
            cout << (*iK)->ID << ", ";
    }
    
    cout << endl << endl << "Mirrors: " << endl << endl;
    for (iG = GROUP->begin(); iG != GROUP->end(); iG++) {
        K = (*((*iG)->KNOTS->begin()))->deepCopy();
        K->mirror();
        poly = HOMFLY(K); poly->simplify();
        
        IDA = (*iG)->group_ID;
        IDB = group_index_from_poly(GROUP, poly);
        
        if (IDA <= IDB) {
            
            cout << " (";
            for (iK = (*iG)->KNOTS->begin(); iK != (*iG)->KNOTS->end(); iK++)
                if (!classification[(*iK)->ID].duplicate) 
                    cout << ((*iK)->ID) << " ";
            
            iGB = group_from_poly(GROUP, poly);
            
            for (iK = (*iGB)->KNOTS->begin(); iK != (*iGB)->KNOTS->end(); iK++)
                if (!classification[(*iK)->ID].duplicate)
                    cout <<((*iK)->ID) << "";
         
            cout << "),";
        }
    }
 
    
    cout << endl << endl;
    
}

// remove knots in the groups that are duplicates (non-prime)
void remove_group_duplicates(t_group_list *G) {
    bool b;
    int ID;
    do {
        b = false;
        for (t_group_list::iterator iG = G->begin(); iG != G->end(); iG++) { // loop through all polynomial groups
            for (t_knot_set::iterator iK = (*iG)->KNOTS->begin(); iK != (*iG)->KNOTS->end(); iK++) { // loop through all knots in group
                ID = (*iK)->ID;
                if ((classification[ID].duplicate) || (classification[ID].composite)) {
                    (*iG)->KNOTS->erase(iK);
                    b = true;
                    break;
                } // if
            } // knots
            if (b) break;
        } // group
    } while (b); // do
}

// are punctured regions adjacent in K?
bool adjacentRegionsSelected(Cknot *K) {
    return (bool)(K->regions[K->reg_0] & K->regions[K->reg_1]);
}

// is K obviously a connected sum?
bool obious_composite(Cknot *K) { 
    u_number dot, r, reg, /*reg_,*/ arc_0, arc_1;
    for (dot = 0; dot <= 1; dot++) {
        reg = (dot?K->reg_0:K->reg_1); // observing region
        if (count_bits(K->regions[reg]) > 1) { // has more than 1 arc
            for (arc_0=0;arc_0<K->n2-1;arc_0++)
                if (BIT(K->regions[reg],arc_0))
                    for (arc_1=arc_0+1;arc_1<K->n2;arc_1++)
                        if (BIT(K->regions[reg],arc_1))
                            for (r=0;r<K->n_reg;r++)
                                if ((r != reg) && (BIT(K->regions[r],arc_0) && BIT(K->regions[r],arc_1) )) return true;
        }
    }
    return false;
}

// HSM(A) = HSM(B) ?
bool HOMFLY_equal(Cknot *A, Cknot *B) {
    bool b;
    Cknot *A_ = A->deepCopy();
    Cknot *B_ = B->deepCopy();
    Cmultivariate * poly_A = HOMFLY(A_); poly_A->simplify();
    Cmultivariate * poly_B = HOMFLY(B_); poly_B->simplify();
    b = poly_equal_orientation(poly_A, poly_B);
    delete poly_A; delete poly_B;
    delete A_; delete B_;
    return b;
}

// caclulate the KBSMs of affine knots up to n = 5, must be called before is_connected_sum() or get_knot_property
void init_prime_affine() {
    // Kauffman bracket skein modules
    prime_KBSM[0] = mv("a^4 + a^12 - a^16"); // left trefoil
    prime_KBSM[1] = mv("-a^-16 + a^-12 + a^-4"); // right trefoil
    prime_KBSM[2] = mv("1 + a^-8 - a^-4 - a^4 + a^8"); // figure eight knot
    prime_KBSM[3] = mv("a^8 + a^16 - a^20 + a^24 - a^28"); // left pentafoil
    prime_KBSM[4] = mv("-a^-28 + a^-24 - a^-20 + a^-16 + a^-8"); // right pentafoil
    prime_KBSM[5] = mv("a^4 - a^8 + 2*a^12 - a^16 + a^20 - a^24"); // left three-twist knot
    prime_KBSM[6] = mv("-a^-24 + a^-20 - a^-16 + 2a^-12 - a^-8 + a^-4"); // right three-twist knot
    
    // Modified Kauffman bracket skein modules
    prime_KBSM_mod[0] = mv("a^18 - a^10 - a^6 - a^2"); // left trefoil
    prime_KBSM_mod[1] = mv("-a^-2 - a^-6 - a^-10 + a^-18"); // right trefoil
    prime_KBSM_mod[2] = mv("-a^10 - a^-10"); // figure eight knot
    prime_KBSM_mod[3] = mv("a^30 - a^14 - a^10 - a^6"); // left pentafoil
    prime_KBSM_mod[4] = mv("-a^-6 - a^-10 - a^-14 + a^-30"); // right pentafoil
    prime_KBSM_mod[5] = mv("a^26 - a^14 - a^10 - a^2"); // left three-twist knot
    prime_KBSM_mod[6] = mv("-a^-2 - a^-10 - a^-14 + a^-26"); // right three-twist knot
    prime_KBSM_mod[7] = mv("-a^2 - a^-2"); // unknot in modified
    
}

// obsolete
/*void init_prime_affine() {
    Cknot *primes[8];
    primes[0] = new Cknot("Knot 0: 1 -2 3 -1 2 -3 --- 03 *&024 14 135 25"); // left trefoil
    primes[1] = primes[0]->deepCopy(); primes[1]->mirror(); // right trefoil
    primes[2] = new Cknot("Knot 0: 1 -2 3 -1 4 -3 2 -4 +--+ 025 036 15 *&146 247 37"); // figure 8 knot
    primes[3] = new Cknot("Knot 0: 1 -2 3 -4 5 -1 2 -3 4 -5 ----- 05 *&02468 13579 16 27 38 49"); // left pentafoil
    primes[4] = primes[3]->deepCopy(); primes[4]->mirror(); // right pentafoil
    primes[5] = new Cknot("Knot 0: 1 -2 3 -4 5 -1 4 -3 2 -5 ----- 0257 048 17 *&1368 26 359 49"); // left three-twist knot
    primes[6] = primes[5]->deepCopy(); primes[6]->mirror(); // right three-twist knot
    primes[7] = new Cknot("Knot 0: 1 -1 + *&0 1 01"); // left trefoil
    for (int i=0;i<8;i++)
        prime_KBSM[i] = KBSM(primes[i]);
    for (int i=0;i<8;i++) delete primes[i];
}*/

// is the knot with the KBSM p a candidate for a connected sum with a trivial projective knot?
// obsolete, replaced with get_knot_property()
/*bool connected_sum_with_projective_trivial(Cmultivariate * poly) {
    bool bp1 = true, bm1 = true;
    for (list<Cterm*>::iterator it = poly->terms.begin(); it != poly->terms.end(); it++) {
        for(int i=0;i<n_v;i++)
            if (is_capital(variables[i])) {
                if ((variables[i] != 'N') && (variables[i] != 'L'))
                    if ((*it)->power[i] != 0) {bp1 = bm1 = false;}
            }
        if ( (*it)->power[reverse_variable['N']] != 1) {bp1 = false; }
        if ( (*it)->power[reverse_variable['L']] != 1) {bm1 = false; }
    }
    return (bp1 || bm1);
}*/

// is K a connected sum (with an affine knot) candidate?
// obsolete, replaced with get_knot_property()
/*bool is_connected_sum(Cknot *K) {
    int i;
    bool b = false;
    Cmultivariate * poly = KBSM(K);
    for (i=0;i<7;i++) b |= divisible(poly, prime_KBSM[i]);
    return b;
    
}*/

// get properties of the knot K (connected sum, affine) in solid torus (l = -1) or L(p,q)
u16 get_knot_property(Cknot *K, int l) {
    u16 property = 0;
    
    int i;
    //Cmultivariate * hsm; // obsolete
    Cmultivariate * kbsm_mod; // modified KBSM
    Cmultivariate * kbsm; // KBSM
    Cmultivariate * poly;
    
    kbsm_mod = (l == -1 ? KBSM(K) : KBSM_lpq(l, KBSM(K))); // is the knot in the solid torus or L(p,q) ?
    //kbsm = (*kbsm_mod) * (*(mv("-a^2 - a^-2")));
    
    kbsm = (l != -1 ? kbsm_mod : (*kbsm_mod * (*(mv("-a^2 - a^-2")))));
    
    for (i=0;i<7;i++) if (divisible(kbsm, prime_KBSM_mod[i])) property |= PROPERTY_SUM;
    
    normalize_framing(kbsm);
    for (i=0;i<8;i++) {
        poly = prime_KBSM_mod[i]->deepCopy();
        normalize_framing(poly);
        if (*poly == *kbsm) property |= PROPERTY_AFFINE;
        delete poly;
    }
    
    
   /* Cmultivariate * p_t = new Cmultivariate("-a^2 - a^-2");
    normalize_framing(p_t);
    p_t->simplify();
    
    for (int l = 0; l < N_HOMEO; l++) {
        cout << lpq_s[l] << (lens_homeo[l][0] > 9 ? "" :" ")<<": ";
        for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
            if (classification[(*iK)->ID].duplicate) continue;
            Cmultivariate *p = KBSM_lpq(l,KBSM(*iK));
            normalize_framing(p);
            p->simplify();
            if (*p == *p_t) cout << output_knot_name((*iK)->ID) << ", ";
        }
        cout << endl;
     */
        
    // obsolete
    /*if ( l != -1)
        if (lens_homeo[l][1] <= 2) {
            hsm = HSM_lpq(l, HOMFLY(K));
            if (is_affine(hsm)) property |= PROPERTY_AFFINE;
        }
    */
    return property;
}



// is K a candidate for a connected sum in L(p,q)?
// obsolete, replaced with get_knot_property()
/*int is_connected_sum_candidate_lpq(int l, Cknot *K) {
    bool b = false;
    Cmultivariate * q;
    Cmultivariate * poly = KBSM_lpq(l, KBSM(K));
    if (lens_homeo[l][1] <= 2) q = HSM_lpq(l, HOMFLY(K));
    if (lens_homeo[l][1] <= 2) if (connected_sum_with_projective_trivial(q)) return 3;
    if (lens_homeo[l][1] <= 2) if (is_affine(q)) return 2;
    for (int i=0;i<7;i++) b |= divisible(poly, prime_KBSM[i]);
    return (b ? 1 : 0);
}*/

#define MAX_SETS 100 // maximal depth of the breadth-first search BST algorithm

// class to manipulate the binary search trees
class CBST {
public:
    t_knot_set *BST[MAX_SETS]; // the BST at certain depth
    t_knot_set *BST_ALL; // stores all knots
  
    int max_depth;
    bool over_limit, critical_limit, affine_knot; // has the BST search reached certain memory limits?
    
    bool compare_with_third_knot; // are we comparing the knot with itself or another knot?
    Cknot *compare_knot; // if comparing to another knot, the knot is stored here
    
    bool insert(int depth, Cknot *Q); // insert a new knot at depth
    bool find(Cknot *Q); // find knot in the BST
    
    t_knot_set::iterator find_and_return_knot(Cknot *Q); // return the iterator of Q in BST
    
    bool smaller_knot_found(); // knot has been reduced
    void set_compared_knot(Cknot *Q); // init knot we are comparing to 
    
    Cknot *minimal(); // gets the minimal knot in the BST
    
    // constructor
    CBST(int max_depth);
    
    // destructor
    ~CBST();

};

// get the minimal knot of the BST (exclude compare knot, which we add artificially to the BST)
Cknot *CBST::minimal() {
    for (t_knot_set::iterator iK = BST_ALL->begin(); iK != BST_ALL->end(); iK++)
        if (*iK != compare_knot) return *iK;
    throw 20;
}

// knot has been reduced?
bool CBST::smaller_knot_found() {
    return ((affine_knot) || (*BST_ALL->begin() != ( compare_with_third_knot ? compare_knot : *BST[0]->begin() )));
}

// insert new knot
bool CBST::insert(int depth, Cknot *Q) {
    pair<t_knot_set::iterator, bool> ret;
    
    affine_knot |= (Q->reg_0 == Q->reg_1);
    affine_knot |= adjacentRegionsSelected(Q);
    
    ret = BST_ALL->insert(Q); // insert the new knot into BST
    
    if (ret.second) { // insertation successful, check for memory
        
        BST_last_knot = Q;
     
        BST[depth]->insert(Q); // if inserted, insert also to BST at current depth

        // check if memory has reached its limits
        if (BST_ALL->size() % 5000 == 0) {
            if (!over_limit) over_limit = (get_used_memory_kb() > installed_ram_kb*3/5);
            if (!critical_limit) critical_limit = (get_used_memory_kb() > installed_ram_kb);
        }
    }
    
    return ret.second;
}

// Q in BST?
bool CBST::find(Cknot *Q) {
    return (BST_ALL->find(Q) != BST_ALL->end());
}

// return iterator of Q in BST
t_knot_set::iterator CBST::find_and_return_knot(Cknot *Q) {
    return BST_ALL->find(Q);
}

// insert the knot we are comparing to
void CBST::set_compared_knot(Cknot *Q) {
    compare_with_third_knot = true;
    compare_knot = Q;
    insert(0,Q);
}

// constructor, initialization
CBST::CBST(int max_depth_) {
    max_depth = max_depth_;
    over_limit = false;
    critical_limit = false;
    compare_with_third_knot = false;
    affine_knot = false;
    compare_knot = NULL;
    BST_ALL = new t_knot_set(); // create BS tree
    BST_last_knot = NULL;
    for (int d=0; d <= max_depth; d++) BST[d] = new t_knot_set();
}

// destructor, free memory
CBST::~CBST() {
    if (BST_last_knot == NULL) cout << "null" << flush;
    BST_last_knot = ((BST_last_knot == NULL) ? NULL : BST_last_knot->deepCopy());
    BST_minimal_knot = minimal()->deepCopy();
    for (int d=0; d<=max_depth;d++) { BST[d]->clear(); delete BST[d]; }
    for (t_knot_set::iterator iQ = BST_ALL->begin(); iQ != BST_ALL->end(); iQ++) { /*if (c++ < 10) cout << *iQ << endl; */delete (*iQ); }
   // BST_ALL->clear(); delete BST_ALL;
}

// main breadth-first search algorithem using BSTs
// compare_knot: use NULL if we are not finding a reduction of a particular knot, but only the knot K
// K_: knot we are reducing
// oriented: true/false (consider knot up to orientation?)
// max_depth: depht of the breadth-first search
// b_increase: set of bits at which we pefrorm a Reidemeister increasing move (e.g. 0b1001001 - perfoms a R. incr. move at depth 1, 4 and 7
// max_additional_crossings: maximal crossings we wish the knot to grow
bool BST_reducible(Cknot *compare_knot_, Cknot *K_, bool oriented, int max_depth, u64 b_increase, int max_additional_crossings) {
    
    t_reidemeister_list rList; // list of all Reidemeister moves of current knot
    t_reidemeister_list::iterator ir;
    t_knot_set::iterator iK;
    CBST BST(max_depth);
    Cknot *Q; // working knot
    
    int depth, repeat;
    u64 r_filter;
    bool force_increase;
    
    // set compared knot, if exists
    if (compare_knot_ != NULL) {
        
        Q = compare_knot_->deepCopy();
        Q->canonical_meridional(oriented);
        BST.set_compared_knot(Q);
    }
    
    Q = K_->deepCopy(); // make a copy of knot and check if its in its canonical form
    Q->canonical_meridional(oriented); // TODO: check speed
    
    BST.insert(0,Q); // insert root node
    
    if ((compare_knot_ == NULL) && (smaller_knot(Q,K_))) return YES; // check if knot is smaller after a canonical meridional operation
    
    // main loop
    for (depth = 1; depth <= max_depth; depth++) {

        force_increase = NO; 
        
        for (repeat=0; repeat <= 1; repeat++) { // if no Reidemeister moves found, force a crossing-increasing R-move, repeat = 1 sets force_increase to YES, otherwise loop breaks
            
            for (iK = BST.BST[depth-1]->begin(); iK != BST.BST[depth-1]->end(); iK++) // scan previous depth level
                
                if (*iK != BST.compare_knot) { // ignore knot we are comparing to
                    
                    r_filter = REMOVE_R_II | REMOVE_R_I | MODIFY_R_III | MODIFY_FLYPE; // non-crossing-increasing moves
                    
                    if (!BST.over_limit) { // if memory nearly critical, do not perform crossing-increasing moves
                        // should we perform a crossing-increasing move?
                        if (( (*iK)->n + 1 <= K_->n+max_additional_crossings) && ((BIT(b_increase, depth-1))/* || (force_increase)*/)) r_filter |= CREATE_R_I;
                        if (( (*iK)->n + 2 <= K_->n+max_additional_crossings) && ((BIT(b_increase, depth-1))|| (force_increase))) r_filter |= CREATE_R_II;
                    }
                    
                    // find possible R-moves
                    findReidemeisterMoveSites((*iK), r_filter, &rList);
                    
                    for (ir = rList.begin(); ir != rList.end(); ir++) { // loop through R-moves
                        
                        Q = (*iK)->deepCopy();
                        
                        rMove(Q,&(*ir)); // perform the Reidemeister move
                        
                        Q->canonical_meridional(oriented); // knot into its canonical form
                        
                        if (BST.insert(depth,Q)) { // knot inserted?
                            
                            if (BST.smaller_knot_found()) return YES; // reduction found
                            
                        } else delete Q; // if no insertation, we can free memory
                
                        if (BST.critical_limit) break; // break loop if critical memory reached
                    } // reidemeister
                    
                    if (BST.critical_limit) break; // break loop if critical memory reached
                } // K
            
                if ((repeat==0) && (BST.BST[depth]->size()>0)) break; // new knots found
                if ((repeat==0) && (BST.BST[depth]->size()==0)) force_increase = YES; // no new knots found, force increasing Reidemeister moves
                if (BST.critical_limit) break;
        } // repeat
        
        if (BST.critical_limit) break;
    } // depth

    return NO; // reduction not found
}

// TODO: join with BST_reducible
bool BST_reducible_up_to_n_crossings(Cknot *K_, bool oriented, int max_depth, u64 b_increase, int max_additional_crossings, int n) {
    
    t_reidemeister_list rList; // list of all Reidemeister moves of current knot
    t_reidemeister_list::iterator ir;
    t_knot_set::iterator iK;
    CBST BST(max_depth);
    Cknot *Q; // working knot
    
    int depth, repeat;
    u64 r_filter;
    bool force_increase;
    
    Q = K_->deepCopy(); // make a copy of knot and check if its in its canonical form
    Q->canonical_meridional(oriented);
    
    if (BST.insert(0,Q)) // insert root
        if (BST.smaller_knot_found()) return YES;
    
    for (depth = 1; depth <= max_depth; depth++) {
        
        // cout <<  "." << flush;
        
        // cout << depth << " "<< BST.BST_ALL->size()/1000 << endl;
        force_increase = NO;
        
        for (repeat=0; repeat <= 1; repeat++) { // no new knots?
            for (iK = BST.BST[depth-1]->begin(); iK != BST.BST[depth-1]->end(); iK++) // previous level
                if (*iK != BST.compare_knot) { // ignore knot we are comparing to
                    if (depth % 3 == 1) r_filter |= CREATE_R_I + CREATE_R_II;
                    r_filter = REMOVE_R_II | REMOVE_R_I | MODIFY_R_III | MODIFY_FLYPE;
                    
                    if (!BST.over_limit) { // should we perform crossing-increasing moves?
                        if (( (*iK)->n + 1 <= K_->n+max_additional_crossings) && ((BIT(b_increase, depth-1))/* || (force_increase)*/)) r_filter |= CREATE_R_I;
                        if (( (*iK)->n + 2 <= K_->n+max_additional_crossings) && ((BIT(b_increase, depth-1))|| (force_increase))) r_filter |= CREATE_R_II;
                    }
                    
                    findReidemeisterMoveSites((*iK), r_filter, &rList);
                    
                    for (ir = rList.begin(); ir != rList.end(); ir++) { // loop through R-moves
                        
                        Q = (*iK)->deepCopy();
                        
                        rMove(Q,&(*ir));
                        
                        Q->canonical_meridional(oriented);
                        
                        if (BST.insert(depth,Q)) { // knot inserted
                            
                            if (Q->n <= n) return YES;
                           // if (BST.smaller_knot_found()) return YES;
                            
                        } else delete Q;
                        
                        if (BST.critical_limit) break;
                    } // reidemeister
                    
                    if (BST.critical_limit) break;
                } // K
            
            if ((repeat==0) && (BST.BST[depth]->size()>0)) break; // new knots found
            if ((repeat==0) && (BST.BST[depth]->size()==0)) force_increase = YES; // no new knots found, force increasing Reidemeister moves
            if (BST.critical_limit) break;
        } // repeat
        
        if (BST.critical_limit) break;
    } // depth
    
    //cout << BST.minimal() << endl;
    return NO;
    
}

// transform A and B and see if there are intersections, OBSOLETE
bool BST_two_sided_test(Cknot *A_, Cknot *B_, bool oriented, int max_depth, u64 b_increase, int max_additional_crossings) {
    
    t_reidemeister_list rList; // list of all Reidemeister moves of current knot
    t_reidemeister_list::iterator ir;
    t_knot_set::iterator iK;
    CBST A_BST(max_depth);
    CBST B_BST(max_depth);
    Cknot *Q; // working knot
    
    int depth, repeat;
    u64 r_filter;
    bool force_increase;
    

    Q = A_->deepCopy();
    Q->canonical_meridional(oriented);
    A_BST.insert(0,Q);
    
    Q = B_->deepCopy(); // make a copy of knot and check if its in its canonical form
    Q->canonical_meridional(oriented);
    
    if (B_BST.insert(0,Q)) if (B_BST.smaller_knot_found()) return YES;
    
    for (depth = 1; depth <= max_depth; depth++) {
        // A
        force_increase = NO;
        for (repeat=0; repeat <= 1; repeat++) { // no new knots?
            
            for (iK = A_BST.BST[depth-1]->begin(); iK != A_BST.BST[depth-1]->end(); iK++) {
         
                r_filter = REMOVE_R_II | REMOVE_R_I | MODIFY_R_III | MODIFY_FLYPE;
                if (!A_BST.over_limit) { // should we perform crossing-increasing moves?
                    if (( (*iK)->n + 1 <= A_->n+max_additional_crossings) && ((BIT(b_increase, depth-1))/* || (force_increase)*/)) r_filter |= CREATE_R_I;
                    if (( (*iK)->n + 2 <= A_->n+max_additional_crossings) && ((BIT(b_increase, depth-1))|| (force_increase))) r_filter |= CREATE_R_II;
                }
                findReidemeisterMoveSites((*iK), r_filter, &rList);
                for (ir = rList.begin(); ir != rList.end(); ir++) { // loop through R-moves
                    Q = (*iK)->deepCopy();
                    rMove(Q,&(*ir));
                    Q->canonical_meridional(oriented);
                    if (!A_BST.insert(depth,Q)) delete Q;
                    if (A_BST.critical_limit) break;
                } // r
                if (A_BST.critical_limit) break;
            } // k
            if ((repeat==0) && (A_BST.BST[depth]->size()>0)) break; // new knots found
            if ((repeat==0) && (A_BST.BST[depth]->size()==0)) force_increase = YES; // no new knots found, force increasing Reidemeister moves
            if (A_BST.critical_limit) break;
        } // repeat
        
        if (A_BST.critical_limit) break;
        
        // B
        force_increase = NO;
        for (repeat=0; repeat <= 1; repeat++) { // no new knots?
            
            for (iK = B_BST.BST[depth-1]->begin(); iK != B_BST.BST[depth-1]->end(); iK++) {
                
                r_filter = REMOVE_R_II | REMOVE_R_I | MODIFY_R_III | MODIFY_FLYPE;
                if (!B_BST.over_limit) { // should we perform crossing-increasing moves?
                    if (( (*iK)->n + 1 <= A_->n+max_additional_crossings) && ((BIT(b_increase, depth-1))/* || (force_increase)*/)) r_filter |= CREATE_R_I;
                    if (( (*iK)->n + 2 <= A_->n+max_additional_crossings) && ((BIT(b_increase, depth-1))|| (force_increase))) r_filter |= CREATE_R_II;
                }
                findReidemeisterMoveSites((*iK), r_filter, &rList);
                for (ir = rList.begin(); ir != rList.end(); ir++) { // loop through R-moves
                    Q = (*iK)->deepCopy();
                    rMove(Q,&(*ir));
                    Q->canonical_meridional(oriented);
                    if (B_BST.insert(depth,Q)) {
                        //if (B_BST.smaller_knot_found()) return YES;
                        if (A_BST.find(Q)) {
                            cout << "B: " << Q << endl;
                            cout << "A: " << *A_BST.find_and_return_knot(Q) << endl;
                            return YES;
                        }
                        
                    } else delete Q;
                    if (B_BST.critical_limit) break;
                } // r
                if (B_BST.critical_limit) break;
            } // k
            if ((repeat==0) && (A_BST.BST[depth]->size()>0)) break; // new knots found
            if ((repeat==0) && (A_BST.BST[depth]->size()==0)) force_increase = YES; // no new knots found, force increasing Reidemeister moves
            if (B_BST.critical_limit) break;
        } // repeat
        if (A_BST.critical_limit) break;
    } // depth
    return NO;
    
}

// clears a list of knots
void clear_list(list<Cknot *> * k_list) {
    for (list<Cknot *>::iterator ik = k_list->begin(); ik != k_list->end(); ik++)
        delete *ik;
    k_list->clear();
}

// simple reduction: elliminates crossings by R-I and R-I while possible (seperator function to gain performance)
bool elliminate_crossings_while_possible(Cknot **Q) {
    bool b = false;
    t_reidemeister_list rList;
    do {
        // rList.clear();
        findReidemeisterMoveSites(*Q, REMOVE_R_II | REMOVE_R_I, &rList);
        if (rList.size() == 0) return b;
        rMove(*Q,&(*(rList.begin())));
        b = true;
    } while (1);
    throw 20;
}

// reduces K by not searching the entire BST tree, but only performing non-crossing increasing Reidemeister moves
// i.e. performs all R-II moves at a given depth and chooses the minimal knot to manipulate it at the next stage
// used in the calulation of the HOMLFY, to minimize knot at each skein step and when minimizing the knot after the slide moves
// much faster then BST_reducible
// stops when no crossing-decreasing move is found max_row times in a row
// TODO: use CBST class
Cknot *BST_shrink(Cknot *K_, bool oriented, int max_depth, int max_row) {

    max_depth = MINIMUM(max_depth, 100); // ?
    
    int depth; // current depth
    t_reidemeister_list rList; // list of all Reidemeister moves of current knot
    t_reidemeister_list::iterator ir;
    t_knot_set *BST[max_depth+1]; // binary search tree (depth by depth), only references (originals in BST_ALL)
    t_knot_set *BST_ALL; // BST of all knots processed so far
    t_knot_set::iterator iK, iQ, iK_;
    pair<t_knot_set::iterator, bool> ret;
    Cknot *Q;
    Cknot *K_minimized;
    bool could_elliminate;
    
    int no_crossing_reduction_in_a_row;
    
    Q = K_->deepCopy(); // make a copy of knot and check if its in its canonical form
    Q->canonical(oriented); // TODO: meridional?
    
    // insert root knot to BST
    BST_ALL = new t_knot_set(); // create BS tree
    BST[(depth = 0)] = new t_knot_set();
    BST_ALL->insert(Q);
    BST[0]->insert(Q);
    
    no_crossing_reduction_in_a_row = 0;
    
    for (depth = 1; ((depth <= max_depth) && (BST[depth-1]->size() > 0) && (no_crossing_reduction_in_a_row <= max_row)); depth++) { // go throug all depths
        
        BST[depth] = new t_knot_set();
        could_elliminate = false;
        
        for (iK = BST[depth-1]->begin(); iK != BST[depth-1]->end(); iK++) {
            
            // first try to elliminate crossings
            
            Q = (*iK)->deepCopy();
            if (elliminate_crossings_while_possible(&Q)) {
                Q->canonical(oriented);
                ret = BST_ALL->insert(Q);
                if (!ret.second) throw 20; // crossings removed but knot found in the BST
                
                BST[depth]->clear(); // put only new knot in binary search tree
                BST[depth]->insert(Q);
                could_elliminate = true;
                break;
            }
            delete Q;
            
            // if no crossings were elliminated, try all R-III moves
            
            findReidemeisterMoveSites((*iK), MODIFY_R_III, &rList);
            
            for (ir = rList.begin(); ir != rList.end(); ir++) {
                Q = (*iK)->deepCopy();
                rMove(Q,&(*ir));
                Q->canonical(oriented);
                ret = BST_ALL->insert(Q);
                if (ret.second) { // insertation OK?
                    BST[depth]->insert(Q);
                } else delete Q;
            }
        } // previous BST level
        
        if (could_elliminate) no_crossing_reduction_in_a_row = 0; else no_crossing_reduction_in_a_row++;
        
    } // depth
    
    
    // clean up
    K_minimized = (*BST_ALL->begin())->deepCopy();
    for (int i=0; i<depth;i++) { BST[i]->clear(); delete BST[i]; }
    for (iQ = BST_ALL->begin(); iQ != BST_ALL->end(); iQ++) delete (*iQ);
    BST_ALL->clear(); delete BST_ALL;
    
    return K_minimized;
   
}

// BST_shrink, but stops when knot is reduced to n crossings
// TODO: join with BST_shrink
Cknot *BST_shrink_up_to_n_crossings(Cknot *K_, bool oriented, int max_depth, int max_row, int n) {
    max_depth = MINIMUM(max_depth, 100); // ?
    int depth; // current depth
    t_reidemeister_list rList; // list of all Reidemeister moves of current knot
    t_reidemeister_list::iterator ir;
    t_knot_set *BST[max_depth+1]; // binary search tree (depth by depth), only references (originals in BST_ALL)
    t_knot_set *BST_ALL; // BST of all knots processed so far
    t_knot_set::iterator iK, iQ, iK_;
    pair<t_knot_set::iterator, bool> ret;
    Cknot *Q;
    Cknot *K_minimized;
    bool could_elliminate;
    
    bool quit_shrink;
    
    int no_crossing_reduction_in_a_row;
    
    Q = K_->deepCopy(); // make a copy of knot and check if its in its canonical form
    Q->canonical(oriented); // meridional?
    
    // insert root knot to BST
    BST_ALL = new t_knot_set(); // create BS tree
    BST[(depth = 0)] = new t_knot_set();
    BST_ALL->insert(Q);
    BST[0]->insert(Q);
    
    no_crossing_reduction_in_a_row = 0;
    quit_shrink = false;
    
    for (depth = 1; ((depth <= max_depth) && (BST[depth-1]->size() > 0) && (no_crossing_reduction_in_a_row <= max_row)); depth++) { // go throug all depths
        
        BST[depth] = new t_knot_set();
        could_elliminate = false;

        for (iK = BST[depth-1]->begin(); iK != BST[depth-1]->end(); iK++) {
            
            // first try to elliminate crossings
            
            Q = (*iK)->deepCopy();
            
            if (elliminate_crossings_while_possible(&Q)) {
                
                Q->canonical(oriented);
                ret = BST_ALL->insert(Q);
                if (!ret.second) throw 20;
          
                BST[depth]->clear(); // put only new knot in binary search tree
                BST[depth]->insert(Q);
                could_elliminate = true;
                
                if (Q->n <= n) quit_shrink = YES;
                
                break;
            }
            delete Q;
            
            if (quit_shrink) break;
            
            // if no crossings were elliminated, try R-III
            
            findReidemeisterMoveSites((*iK), MODIFY_R_III, &rList);
            
            for (ir = rList.begin(); ir != rList.end(); ir++) {
                Q = (*iK)->deepCopy();
                rMove(Q,&(*ir));
                Q->canonical(oriented);
                ret = BST_ALL->insert(Q);
                
                if (ret.second) { // insertation OK?
                    BST[depth]->insert(Q);
                    if (Q->n <= n) { quit_shrink = YES; break; }
                } else delete Q;
                
            }
            if (quit_shrink) break;
            
        } // previous BST level
        
        if (could_elliminate) no_crossing_reduction_in_a_row = 0; else no_crossing_reduction_in_a_row++;
        if (quit_shrink) break;
    } // depth
    
    K_minimized = (*BST_ALL->begin())->deepCopy();
    for (int i=0; i<depth;i++) { BST[i]->clear(); delete BST[i]; }
    for (iQ = BST_ALL->begin(); iQ != BST_ALL->end(); iQ++) delete (*iQ);
    BST_ALL->clear(); delete BST_ALL;
    
    return K_minimized;
}

// BST_shrink, but replaces the knot with the minimized one (i.e. does not return the reduced knot)
void BST_shrink(Cknot **K_, bool oriented, int max_depth, int max_row) {
    Cknot *Q = BST_shrink(*K_, oriented, max_depth, max_row);
    (*K_)->copyKnotData(Q);
    delete Q;
    
}

// TODO: join with BST_shrink
void BST_shrink_up_to_n_crossings(Cknot **K_, bool oriented, int max_depth, int max_row, int n) {
    Cknot *Q = BST_shrink_up_to_n_crossings(*K_, oriented, max_depth, max_row,n);
    (*K_)->copyKnotData(Q);
    delete Q;
}

// returns if K can be "BST_shrink"-able
bool BST_shrinkable(Cknot *K_, bool oriented, int max_depth, int max_row) {
    Cknot *Q = BST_shrink(K_, oriented, max_depth, max_row);
    bool shrinkable = smaller_knot(K_, Q);
    delete Q;
    return shrinkable;
}

// calculates the HSMs of all the knots (in the solid torus) in the set KNOTS, stores to HOMFLYS and HOMFLYS_reverse 
void calculate_HSMs(t_knot_set * KNOTS) {
    Cmultivariate *poly, *poly_reverse;
    Cknot *K;
    
    HOMFLYS.clear();
    HOMFLYS_reverse.clear();
    int i = 0;
    
    // insert the unknot
    HOMFLYS.push_back(new Cmultivariate("- vz^-1 + v^-1z^-1")); // unknot
    HOMFLYS_reverse.push_back(new Cmultivariate("- vz^-1 + v^-1z^-1")); // unknot
    for (t_knot_set::iterator it = KNOTS->begin(); it != KNOTS->end(); it++) {
       // cout << (i+1) << "," << (*it)->ID << endl;
        if ((++i) != (*it)->ID) throw 20; // double check if the knot ID is the next knot in the set
        K = (*it)->deepCopy();
        
        poly = HOMFLY(K);
        delete K;
        poly->simplify();
        
        poly_reverse = reverse_orientation(poly);
        poly_reverse->simplify();
        
        HOMFLYS.push_back(poly);
        HOMFLYS_reverse.push_back(poly_reverse);
    }
}

// calculates the KBSMs of all the knots (in the solid torus) in the set KNOTS, stores to KBSMS 
void calculate_KBSMs(t_knot_set * KNOTS) {
    Cmultivariate *poly;
    KBSMS.clear();
    int i = 0;
    
    KBSMS.push_back(new Cmultivariate("-a^2 - a^-2")); // unknot

    for (t_knot_set::iterator it = KNOTS->begin(); it != KNOTS->end(); it++) {
        if ((++i) != (*it)->ID) throw 20; // double check
        poly = KBSM(*it);
        poly->simplify();
        KBSMS.push_back(poly);
    }
}

// obsolete
/*
Cknot * minimize_wrapping_number(int l, Cknot *K) {
    Cknot *Q, *Q_;
    int wind, wind_;
    int goal = lens_homeo[l][0]>>1; //(lens_homeo[l][0]-1)>>1; // floor((p-1)/2)
    
    cout << " L(" <<lens_homeo[l][0]<<"," << lens_homeo[l][1] << ") " << K->n <<": ";
    Q = K->deepCopy();
    wind = winding_number(Q);
    if (wind < 0) { Q->reverse(); wind *= -1; }
    
    while (abs(wind) > goal) {
        cout << "(" <<wind << ")" << flush;
        wind_ = winding_number_after_omega_4(l, Q, 0);
        Q_ = Omega_4(l, Q, (abs(wind_) < abs(wind)) ? 0 : 1);
        delete Q;
        Q = Q_;
        cout << "."<<Q->n << flush;
        BST_shrink(&Q, YES, 20, 6);
        cout << ","<<Q->n << flush;
        wind = winding_number(Q);
        if (wind < 0) { Q->reverse(); wind *= -1; }
    }
    //cout << wind <<endl;
    
    if (wind < 0) Q->reverse();
    return Q;
}*/

// groups (partitions) the set of KNOTS into seperate partitions KNOT_GROUP by the knots HSMs
void groupKnotsByHSM(t_group_list * KNOT_GROUP, t_knot_set * KNOTS) {
    
    t_knot_set::iterator it;
    t_group_list::iterator itg;
    bool inserted;
    CknotGroup * KG;
    Cmultivariate * poly, * poly_reverse;
    
    int gid = 1; // group iD
    
    KNOT_GROUP->clear();
    
    for (it = KNOTS->begin(); it != KNOTS->end(); it++) {
        
        inserted = false;
        poly = (HOMFLYS)[(*it)->ID];
        poly_reverse = (HOMFLYS_reverse)[(*it)->ID];
        
        for (itg = KNOT_GROUP->begin(); itg != KNOT_GROUP->end(); itg++) {
            if ((*((*itg)->poly) == *poly) || (*((*itg)->poly) == *poly_reverse)) {
                (*itg)->KNOTS->insert(*it);
                inserted = true;
                break;
            }
            
        }
        
        if (!inserted) { // create new group of knots
            KG = new CknotGroup();
            KG->group_ID = (gid++);
            KG->poly = poly;
            KG->KNOTS->insert(*it);
            KNOT_GROUP->push_back(KG);
        } //else delete poly;
    }
}

// groups (partitions) the set of KNOTS into seperate partitions KNOT_GROUP by the knots KBSMs
void groupKnotsByKBSM(t_group_list * KNOT_GROUP, t_knot_set * KNOTS) {
    
    t_knot_set::iterator it;
    t_group_list::iterator itg;
    bool inserted;
    CknotGroup * KG;
    Cmultivariate * poly;
    
    int gid = 1; // group ID
    KNOT_GROUP->clear();
    
    for (it = KNOTS->begin(); it != KNOTS->end(); it++) {
        
        inserted = false;
        poly = (KBSMS)[(*it)->ID];
        
        for (itg = KNOT_GROUP->begin(); itg != KNOT_GROUP->end(); itg++) {
            if ((*((*itg)->poly) == *poly)) {
                (*itg)->KNOTS->insert(*it);
                inserted = true;
                break;
            }
        }
        
        if (!inserted) { // create new group of knots
            KG = new CknotGroup();
            KG->group_ID = (gid++);
            KG->poly = poly;
            KG->KNOTS->insert(*it);
            KNOT_GROUP->push_back(KG);
        } //else delete poly;
    }
}

// groups (partitions) the set of KNOTS in all L(p,q) into seperate partitions KNOT_GROUP[-] by the knots HSMs
void groupKnotsByHSM_lpq(t_group_list * KNOT_GROUP[N_HOMEO], t_knot_set * KNOTS) {
    
    t_knot_set::iterator it;
    t_group_list::iterator itg;
    bool inserted;
    CknotGroup * KG;
    Cmultivariate * poly, * poly_reverse;
    int l;
    
    int gid = 1; // group ID, TODO: should be initialized after l-loop
    
    for (l=0;l<N_HOMEO; l++) { // loop through all lens spaces
        KNOT_GROUP[l]->clear();
    
        for (it = KNOTS->begin(); it != KNOTS->end(); it++) {
            inserted = false;
            poly = (HOMFLYS_lpq[l])[(*it)->ID];
            poly_reverse = (HOMFLYS_lpq_reverse[l])[(*it)->ID];
        
            for (itg = KNOT_GROUP[l]->begin(); itg != KNOT_GROUP[l]->end(); itg++) {
                if ((*((*itg)->poly) == *poly) || (*((*itg)->poly) == *poly_reverse)) {
                    (*itg)->KNOTS->insert(*it);
                    inserted = true;
                    break;
                }
            }
        
            if (!inserted) { // create new group of knots
                KG = new CknotGroup();
                KG->group_ID = (gid++);
                KG->poly = poly;
                KG->KNOTS->insert(*it);
                KNOT_GROUP[l]->push_back(KG);
            } //else delete poly;
        }
    }
}

// groups (partitions) the set of KNOTS in all L(p,q) into seperate partitions KNOT_GROUP[-] by the knots KBSMSs
void groupKnotsByKBSM_lpq(t_group_list * KNOT_GROUP[N_HOMEO], t_knot_set * KNOTS) {
    
    t_knot_set::iterator it;
    t_group_list::iterator itg;
    bool inserted;
    CknotGroup * KG;
    Cmultivariate * poly;
    int l;
    
    int gid = 1; // knot ID, TODO: should be initialized after l-loop
    
    for (l=0;l<N_HOMEO; l++) {
        KNOT_GROUP[l]->clear();
        
        for (it = KNOTS->begin(); it != KNOTS->end(); it++) {
            inserted = false;
            poly = (KBSMS_lpq[l])[(*it)->ID];
            
            for (itg = KNOT_GROUP[l]->begin(); itg != KNOT_GROUP[l]->end(); itg++) {
                if (*((*itg)->poly) == *poly) {
                    (*itg)->KNOTS->insert(*it);
                    inserted = true;
                    break;
                }
            }
            
            if (!inserted) { // create new group of knots
                KG = new CknotGroup();
                KG->group_ID = (gid++);
                KG->poly = poly;
                KG->KNOTS->insert(*it);
                KNOT_GROUP[l]->push_back(KG);
            } //else delete poly;
            
        }
    }
}

// gets iterator to group by the group ID, TODO: check if function duplicated
t_group_list::iterator get_group_by_id(t_group_list * G, int ID) {
    t_group_list::iterator it;
    for (it=G->begin(); it != G->end(); it++)
        if ( (*it)->group_ID == ID) return it;
    throw 20;
}

// renumerates the ID of KNOTS, used if knots are removed (i.e. after deleting connected sums)
// TODO: check if function is duplicated
void renumerate_prime_knot_table(t_knot_set *KNOTS) {
    t_knot_set::iterator it;
    int kid;
    for (it = KNOTS->begin(), kid = 1; it != KNOTS->end(); kid++, it++)
        if ((*it)->ID != kid) {
            cout << "ERROR: " << (*it) << endl;
            throw 20;
        }
}

// Gauss code generation

// replaces K with the next Gauss code of K wrt to the lexicographical ordering
bool nextRawGaussCode(Cknot* K) { // doesn't check if Gauss word makes sense
    int i, j, k, next;
    u64 count[2];
    count[0] = count[1] = first_bits_full_64(K->n+1) ^ 1;
    
    for (i=K->n2-1;i>=0;i--) {
        count[K->gw_[i]<0] &= ~((u64)1<<ABS(K->gw_[i])); // letter used
        
        next = 0;
        // could there be a different letter on i-th place?
        if (BIT(count[0] | count[1],ABS(K->gw_[i]))) {
            for (j=ABS(K->gw_[i])+1;j<=K->n;j++)
                if (!(BIT(count[0] & count[1],j))) { next = (BIT(count[0],j) ? -j : j); break; }
        } else
            if (K->gw_[i] > 0) next = -K->gw_[i];
        if (next) break;
    }
    //fill up the rest of the GW
    if ((i > 0) && (next)) {
        K->gw_[i] = next;
        count[next<0] |= 1<<ABS(next);
        for (j=i+1;j<2*K->n;j++)
            for (k=1;k<=K->n;k++)
                if (!BIT(count[0],k)) {K->gw_[j] = k; count[0] |= 1<<k; break;}
                else
                    if (!BIT(count[1],k)) {K->gw_[j] = -k; count[1] |= 1<<k; break;}
    } else
        return false;
    return true;
}

// true if letters in even distance (1221), simple check for realizability
bool evenPairing(Cknot *K) { 
    for (int i=1; i<=K->n; i++)
        if (!MOD(K->gw_pos(i)-K->gw_pos(-i),2)) return false;
    return true;
}

// number of kinks in the Gauss word
int numberOfKinks(Cknot *K) { 
    int c = 0;
    for (int i=0; i<2*K->n; i++)
        c += (K->gw(i) == -K->gw(i+1));
    return c;
}

// is there a kink outside the two dots?
bool KinkOutsideDot(Cknot *K) { 
    for (u_number r = 0; r < K->n_reg; r++)
        if ((r != K->reg_0) && (r != K->reg_1) && (count_bits(K->regions[r]) == 1)) return true;
    return false;
}

// order the regions in the lexicographical ordering, does not replace reg_0, reg_1 (i.e. ignores the punctures)
// TODO: rewrite to quicksort
void sortRegions(Cknot *K) { 
    bool change;
    u_number r;
    u_region tmp;
    do {
        change = false;
        for (r=0;r<K->n_reg-1;r++) {
            if (COMPARE_REGION(K->regions[r], K->regions[r+1]) > 0)
            {
                tmp = K->regions[r];
                K->regions[r] = K->regions[r+1];
                K->regions[r+1] = tmp;
                change = true;
            }
        }
    } while (change);
}

// order the regions in the lexicographical ordering, also makes sure the punctures are sorted
// TODO: rewrite to quicksort, join with sort_regions
void sortRegionsDots(Cknot *K) { 
    bool change;
    u_number r;
    u_region tmp;
    do {
        change = false;
        for (r=0;r<K->n_reg-1;r++) {
            if (COMPARE_REGION(K->regions[r], K->regions[r+1]) > 0)
            {
                tmp = K->regions[r];
                K->regions[r] = K->regions[r+1];
                K->regions[r+1] = tmp;
                if (K->reg_0 == r) K->reg_0 = r+1; else
                    if (K->reg_0 == r+1) K->reg_0 = r;
                if (K->reg_1 == r) K->reg_1 = r+1; else
                    if (K->reg_1 == r+1) K->reg_1 = r;
                
                change = true;
            }
        }
    } while (change);
}

// global variables for couting diagrams and Gauss codes
long int count_realizable, count_canonical, count_realizable_extended, count_canonical_extended;

// generates all Gauss codes up to N crossings, counts all (realizable) diagrams, extended diagrams,...
void count_gauss_words(int n) {
    Cknot *K = new Cknot(), *Q;
    int i;
    u_sign sign_bits;
    count_realizable = count_canonical = count_realizable_extended = count_canonical_extended = 0;
    
    K->n = n;
    K->n2 = 2*n;
        
    for (i=0;i<n;i++) {K->gw_[2*i] = i+1; K->gw_[2*i+1] = -(i+1); } // generate lowest possible knot
        do {
            if (!evenPairing(K)) continue;
            
            for (sign_bits = (u_sign)0; sign_bits < ((u_sign)1 << n); ++sign_bits) { // loop through all possible placement of crossings signs
                
                memset(K->regions,0,sizeof(K->regions));
                K->signs = reverseBits(sign_bits,n+1);
                if (!K->realizable()) continue;
                Q = K->deepCopy();
      
                if (!Q->canonical(UNORIENTED)) { delete Q; continue; }
                count_canonical++;
                count_realizable += f(n)*2; // permute letters
                
                K->generateRegions();
                sortRegions(K);
                
                for (K->reg_0=0; K->reg_0 < K->n_reg; K->reg_0++) // put dots in all possible regions
                    for (K->reg_1=0; K->reg_1 < K->n_reg; K->reg_1++) {
                        Q = K->deepCopy();
                        count_realizable_extended += f(n)*2; // permute letters
                        if (!Q->canonical(UNORIENTED)) { delete Q; continue; }
                        count_canonical_extended++;
                    } // reg 1
            } //sign
            //  cout << "ok" <<endl;
        }
        while (nextRawGaussCode(K));
}

// generates all diagrams for N_MIN <= n <= N_MAX crossings, also does a simple reduction check (up to operations in flags)
// returns the set of Gauss codes in KNOTS
void generateAllKnots(int N_MIN, int N_MAX, t_knot_set *KNOTS, u64 flags) {
    
    Cknot *K = new Cknot(), *Q;
    int n,i;
    u_sign sign_bits;
    int count_primes = 0;
    
    KNOTS->clear();
    
    for (n=N_MIN;n<=N_MAX;n++) {
        K->n = n;
        K->n2 = 2*n;
        
        for (i=0;i<n;i++) {K->gw_[2*i] = i+1; K->gw_[2*i+1] = -(i+1); } // generate lowest possible knot
        do {
            
            if (!evenPairing(K)) continue; // "121" ?
            
            if ((flags & REMOVE_3_KINKS) && (numberOfKinks(K) >= 3)) continue; // if any kinks outside punctures, ignore
            
            for (sign_bits = (u_sign)0; sign_bits < ((u_sign)1 << n); ++sign_bits) { // loop through all possible placement of crossings signs
            
                K->signs = reverseBits(sign_bits,n+1);
            
                if (!K->realizable()) continue;
                
                K->generateRegions();
                
                sortRegions(K);
                
                for (K->reg_0=0; K->reg_0 < K->n_reg; K->reg_0++) // put dots in all possible regions
                    for (K->reg_1=0; K->reg_1 < K->n_reg; K->reg_1++) {
                    
                        if ((flags & REMOVE_AFFINE) && (K->reg_0 == K->reg_1)) continue;
                        
                        if ((flags & REMOVE_KINK_OUTSIDE_DOT) && (KinkOutsideDot(K))) continue;
                        
                        if ((flags & REMOVE_ADJACENT_SELECTED) && (adjacentRegionsSelected(K))) continue;
                
                        //if ((remove_equivalent) && (findReidemeisterMoveSites(K, REMOVE_R_II))) continue;
                        
                        if ((flags & REMOVE_OBVIOUS_COMPOSITES) && (obious_composite(K))) continue;
                       
                        if ((flags & REMOVE_BST_REMOVABLE) && (BST_reducible(NULL, K, UNORIENTED, 10, 0, 0))) continue;
                                               
                        if ((flags & REMOVE_SHRINKABLE) && (BST_shrinkable(K, UNORIENTED, K->n *2,6))) continue;
                        
                        //if ((flags & REMOVE_CONNECTED_SUM) && (is_connected_sum(K))) continue;
 
                        if ((flags & REMOVE_CONNECTED_SUM) && (get_knot_property(K) & PROPERTY_SUM)) continue;
                        
                        Q = K->deepCopy();
                        
                        if (!Q->canonical(UNORIENTED)) { delete Q; continue; }
                        
                        Q->ID = ++count_primes;
                                            
                        if (!KNOTS->insert(Q).second) throw 20; // knot duplicated?
                        
                    } // reg 1

            } // signs
                        
        } while (nextRawGaussCode(K));
    } // n
   
}

// prints the candidates for affine knots
// obsolete ?
void get_affine_candidates_lpq(t_group_list *G[N_HOMEO]) {
    t_group_list::iterator iG;
    t_knot_set::iterator iK;
    for (int l=0;l<N_HOMEO; l++) {
        
        for (iG = G[l]->begin(); iG != G[l]->end(); iG++) {
            
            if (is_affine((*iG)->poly)) cout << (*iG)->poly  << endl << "[" << (*iG)->group_ID << "] "; else continue;
            
            for (iK = (*iG)->KNOTS->begin(); iK != (*iG)->KNOTS->end(); iK++) {
                if (classification_lpq[l][(*iK)->ID].duplicate) continue;
                cout << (*iK)->ID << " ";
            }
            cout << endl;
        }
    }
}

// main classification algorithm for knots in the splid torus
// depth: to what depth we perform R-moves
// pattern: at what positions do we perform crossing-increasing R-moves
// start_group, end_group: which groups do we consider checking (in case we wish to seperate the search into subparts)
// stores results in classification array, saves the classification at each step
void classify(int depth, u64 pattern, int additional_crossings, int start_group = 0, int end_group = 0x7FFF) {
    
    int ID; // knot id
    t_knot_set::iterator iK;
    
    for (t_group_list::iterator iG = G_HOMFLY->begin(); iG != G_HOMFLY->end(); iG++)
        if (((*iG)->group_ID >= start_group) && ((*iG)->group_ID <= end_group)) { // loop through all polynomial groups

            if (unclassified_in_group((*iG)->KNOTS, classification) >= 1) {

                //cout << (*iG)->group_ID <<" "<< flush; // print
                
                for (iK = (*iG)->KNOTS->begin(), iK++; iK != (*iG)->KNOTS->end(); iK++) { // loop through all knots in group
                    
                    ID = (*iK)->ID;
                    
                    // do not continue if knot already classified
                    if (classification[ID].duplicate) continue;
                    if (classification[ID].prime) continue;
                    if (classification[ID].composite) continue;
                    
                    cout << (*iK)->ID << " " << flush; // prints
                    
                    classification[ID].duplicate |= BST_reducible(NULL,*iK, UNORIENTED, depth, pattern, additional_crossings);
                    
                    if (classification[ID].prime) cout << "* " << flush; // knot has been classified
                    if (classification[ID].duplicate) cout << "* " << flush; // knot as been classified
                    if ((!classification[ID].prime) && (!classification[ID].duplicate)) cout << "- " << flush; // we were not able to classify the knot

                    OR_loadClassification(CLASSIFICATION_PATH); // load the prevous classification and add join with the new classification
                    saveClassification(CLASSIFICATION_PATH); // save the classification
                }
            } // size > 1
        }
    OR_loadClassification(CLASSIFICATION_PATH);
    saveClassification(CLASSIFICATION_PATH);
}

// classification function for torus, two sided test
// obolete
void classify_two_sided(int depth, u64 pattern, int additional_crossings) {
    
    int ID; // knot id
    t_knot_set::iterator iK;
    
    for (t_group_list::iterator iG = G_HOMFLY->begin(); iG != G_HOMFLY->end(); iG++)
        if (1) { // loop through all polynomial groups
            
            if (unclassified_in_group((*iG)->KNOTS,classification) >= 1) {
            
                cout << (*iG)->group_ID <<" "<< flush;
                
                for (iK = (*iG)->KNOTS->begin(), iK++; iK != (*iG)->KNOTS->end(); iK++) { // loop through all knots in group
                    
                    ID = (*iK)->ID;
                    
                    if (classification[ID].duplicate) continue;
                    if (classification[ID].prime) continue;
                    if (classification[ID].composite) continue;
                    cout << endl << "compare " << (*(*iG)->KNOTS->begin()) << endl << " and    " << *iK << endl << endl;
                    classification[ID].duplicate |= BST_two_sided_test((*(*iG)->KNOTS->begin()), *iK, UNORIENTED, depth, pattern, additional_crossings);
                    
                    if (classification[ID].prime) cout << "* " << flush;
                    if (classification[ID].duplicate) cout << "* " << flush;
                    if ((!classification[ID].prime) && (!classification[ID].duplicate)) cout << "- " << flush;;
                }
                
            } // size > 1
        }
    // TODO: saving?
}

// make an omega move in l depending on parameters extra and reverse
Cknot * omega_move(int l, Cknot *K_, s_small extra, bool reverse, bool merid_rot) {
    Cknot *K, *Q;
    K = K_->deepCopy();
    if (reverse) K->reverse();
    if (merid_rot) K->meridional_rotation();
    
    Q = Omega_4(l, K, extra);
    delete K;
    if (abs(winding_number(Q)) > 6) { delete Q; return NULL; }
    return Q;
}

char ch_;

bool lpq_try_to_reduce(int l, int ID, Cknot *Q_min, Cknot *Q, int depth, u64 pattern, int additional_crossings, int shrink_depth, int shrink_row) {
    
    if (Q == NULL) return classification_lpq[l][ID].duplicate;
    
    BST_shrink(&Q, UNORIENTED, shrink_depth, shrink_row);
    classification_lpq[l][ID].duplicate |= BST_reducible(Q_min, Q, UNORIENTED, depth, pattern,additional_crossings);
    
    if (classification_lpq[l][ID].duplicate) cout << ch_ << ", " << flush;
    
    //cout << (classification_lpq[l][ID].duplicate ? c_ : "") << flush;
    
    return classification_lpq[l][ID].duplicate;
    
}

// main classification algorithm for knots in the splid torus
// depth: to what depth we perform R-moves
// pattern: at what positions do we perform crossing-increasing R-moves
// shrink_depth, chrink_row: the arguments used to call BST_shrink when reducing knot after slide move
// stores results in classification_lpq array, saves the classification at each step
void classify_lpq(int depth, u64 pattern, int additional_crossings, int shrink_depth, int shrink_row) {
    int l, ID; // knot id
    t_knot_set::iterator iK;
    t_group_list * G_KNOTS;
    Cknot *Q;
    
    for (l = 0; l < N_HOMEO; l++) { // loop through all lens spaces
        
        G_KNOTS = GROUP_OF_CHOICE(l)[l]; // select partition of KBSM or HSM
        cout << lpq_s[l] << ": " << flush;
        
        for (t_group_list::iterator iG = G_KNOTS->begin(); iG != G_KNOTS->end(); iG++) { // loop through all groups
            
            for (iK = (*iG)->KNOTS->begin(), iK++; iK != (*iG)->KNOTS->end(); iK++) { // loop through all knots in group
                    
                ID = (*iK)->ID;
                
                if (ID == 168) continue;
                if (ID == 171) continue;
                
                // if knot already classified, skip
                if ((classification_lpq[l][ID].duplicate) || (classification_lpq[l][ID].prime) || (classification_lpq[l][ID].composite)) continue;
                
                cout << ID << " " << flush;
                
                // if a slide move will have more crossings than MAX, continue
                if (estimate_number_of_crossings(l, *iK)+additional_crossings + 1 >= MAX) { cout << "/" << flush; continue;}
                
            
                
                ch_ = 'A';
                
                Q = omega_move(l, *iK, 0, false, false);
                if (lpq_try_to_reduce(l, ID, *iK, Q, depth,pattern,additional_crossings,shrink_depth,shrink_row)) continue;
                
                ch_ = 'B';
                Q = omega_move(l, *iK, 0, true, false);
                if (lpq_try_to_reduce(l, ID, *iK, Q, depth,pattern,additional_crossings,shrink_depth,shrink_row)) continue;
                ch_ = 'C';
                Q = omega_move(l, *iK, 0, false, true);
                if (lpq_try_to_reduce(l, ID, *iK, Q, depth,pattern,additional_crossings,shrink_depth,shrink_row)) continue;
            //    ch_ = 'D';
              //  Q = omega_move(l, *iK, 0, true, true);
                //if (lpq_try_to_reduce(l, ID, *iK, Q, depth,pattern,additional_crossings,shrink_depth,shrink_row)) continue;
          
                
                ch_ = 'X';
                Q = omega_move(l, *iK, 1, false, false);
                if (lpq_try_to_reduce(l, ID, *iK, Q, depth,pattern,additional_crossings,shrink_depth,shrink_row)) continue;
              //  ch_ = 'Y';
              //  Q = omega_move(l, *iK, 1, true, false);
              //  if (lpq_try_to_reduce(l, ID, *iK, Q, depth,pattern,additional_crossings,shrink_depth,shrink_row)) continue;
              //  ch_ = 'Z';
              //  Q = omega_move(l, *iK, 1, false, true);
               /// if (lpq_try_to_reduce(l, ID, *iK, Q, depth,pattern,additional_crossings,shrink_depth,shrink_row)) continue;
               // ch_ = 'W';
               // Q = omega_move(l, *iK, 1, true, true);
               // if (lpq_try_to_reduce(l, ID, *iK, Q, depth,pattern,additional_crossings,shrink_depth,shrink_row)) continue;
                
                cout << "- " << flush;
            } // knots
        } // group
        cout << endl;
    } // lens space 
    save_classification_lpq(CLASSIFICATION_PATH_LPQ);
}


#define OMEGA_MOVES_N 1
void classify_lpq_(int depth, u64 pattern, int additional_crossings, int shrink_depth, int shrink_row) {
    int l, i;
    int ID; // knot id
    t_knot_set::iterator iK;
    t_group_list * G_KNOTS;
    Cknot *K_temp, *K;
    int omega_0, omega_1, rotation;
    
    for (l = 0; l < N_HOMEO; l++) { // loop through all lens spaces
        if (lens_homeo[l][0] < 7) continue;
            
        G_KNOTS = GROUP_OF_CHOICE(l)[l]; // select partition of KBSM or HSM
        cout << lpq_s[l] << ": " << flush;
        
        for (t_group_list::iterator iG = G_KNOTS->begin(); iG != G_KNOTS->end(); iG++) {
            
            if (unclassified_in_group_lpq(l,(*iG)->KNOTS) >= 1) {
                
                //cout << "[" <<(*iG)->group_ID <<"] "<< flush; // print
                
                for (iK = (*iG)->KNOTS->begin(), iK++; iK != (*iG)->KNOTS->end(); iK++) { // loop through all knots in group
                    
                    ID = (*iK)->ID;
                    
                    if ((classification_lpq[l][ID].duplicate) || (classification_lpq[l][ID].prime) || (classification_lpq[l][ID].composite)) continue;
                
                    cout << ID << " " << flush;
                    for (omega_1=0; omega_1<=OMEGA_MOVES_N; omega_1++) {
                        for (omega_0=0; omega_0<=OMEGA_MOVES_N; omega_0++) {
                            for (rotation=0; rotation<=1; rotation++) {
                                
                                if ((omega_0 == 0) && (omega_1 == 0)) continue;
                                if ((omega_0 == 1) && (omega_1 == 1)) continue;
                                
                                if (estimate_number_of_crossings(l, *iK)+additional_crossings + 1 >= MAX) { cout << "/" << flush; continue;}
                                
                                K = (*iK)->deepCopy();
                                
                                if (rotation) K->meridional_rotation();
                                
                                for (i=0;i<omega_0;i++) { K_temp = Omega_4(l,K,0); delete K; K = K_temp; BST_shrink(&K, UNORIENTED, shrink_depth, shrink_row);}
                                for (i=0;i<omega_1;i++) { K_temp = Omega_4(l,K,1); delete K; K = K_temp; BST_shrink(&K, UNORIENTED, shrink_depth, shrink_row);}
                                
                                if (abs(winding_number(K)) > 6) continue;
                                
                                //cout << endl <<K << endl;
                                
                                classification_lpq[l][ID].duplicate |= BST_reducible(*iK, K, UNORIENTED, depth, pattern, additional_crossings);
                                
                                delete K;
                                
                                //cout << endl << BST_minimal_knot << endl;
                                
                                if (classification_lpq[l][ID].duplicate) break;
                            } // rotation
                            if (classification_lpq[l][ID].duplicate) break;
                        } // omega_0
                        if (classification_lpq[l][ID].duplicate) break;
                    } // omega_1
                    
                    cout << (classification_lpq[l][ID].duplicate ? "* " : "- ") << flush;
                    if (classification_lpq[l][ID].duplicate) save_classification_lpq(CLASSIFICATION_PATH_LPQ);
                } // iK
            }
        } // iG
        cout << endl;
    } // l
    save_classification_lpq(CLASSIFICATION_PATH_LPQ);
}

// obsolete functions

// obolete
#define HB(p,a) ((-(float)(p)/2.0 < ((float)a)) && ((float)(p)/2.0 >= ((float)a)))

// obsolete
bool is_HOMFLY_in_lpq_basis(int l, Cmultivariate *p) {
    int gmin, gmax;
    HOMFLY_span(p, &gmin, &gmax);
    cout << "(" << gmin << "," << gmax << ") p = " << lens_homeo[l][0] << " ";
    return ( (HB(lens_homeo[l][0],gmin) && HB(lens_homeo[l][0],gmax)) || (HB(lens_homeo[l][0],-gmax) && HB(lens_homeo[l][0],-gmin)));
}

// obsolete
void calculate_HOMLFY_KBSM_lpq(t_knot_set *KNOTS) {
    
    t_knot_set::iterator iK;
    Cknot *Q;
    
    for (iK = KNOTS->begin(); iK != KNOTS->end(); iK++)
        if (!classification[(*iK)->ID].duplicate) {
            
            if (((*iK)->ID == 484) && ((*iK)->ID >= 0)) { // ?

                cout << (*iK)->ID << endl;
                
                Q = (*iK)->deepCopy();
                cout << (is_HOMFLY_in_lpq_basis(3, HOMFLY(Q)) ? " YES " : " NO ") << winding_number(Q) << endl << endl;
                
                Q = (*iK)->deepCopy();
                Q->reverse();
                cout << (is_HOMFLY_in_lpq_basis(3, HOMFLY(Q)) ? " YES " : " NO ") << winding_number(Q) << endl << endl;
                
                Q = (*iK)->deepCopy();
                Q = Omega_4(3, Q, 0);
                BST_shrink(Q, YES, 20, 6);
                cout << (is_HOMFLY_in_lpq_basis(3, HOMFLY(Q)) ? " YES " : " NO ") << winding_number(Q) << endl << endl;
                
                Q = (*iK)->deepCopy();
                Q = Omega_4(3, Q, 1);
                BST_shrink(Q, YES, 20, 6);
                cout << (is_HOMFLY_in_lpq_basis(3, HOMFLY(Q)) ? " YES " : " NO ") << winding_number(Q) << endl << endl;
  
                Q = (*iK)->deepCopy();
                Q = Omega_4(3, Q, 0); cout << "." << flush;
                BST_shrink(Q, YES, 20, 6);
                Q->reverse();cout << "." << flush;
                Q = Omega_4(3, Q, 1);
                BST_shrink(Q, YES, 30, 10);
                cout << "." << flush;
                cout << (is_HOMFLY_in_lpq_basis(3, HOMFLY(Q)) ? " YES " : " NO ") << winding_number(Q) << endl << endl;
                
                Q = (*iK)->deepCopy();
                Q->reverse();
                Q = Omega_4(3, Q, 1);
                BST_shrink(Q, YES, 20, 6);
                cout << (is_HOMFLY_in_lpq_basis(3, HOMFLY(Q)) ? " YES " : " NO ") << winding_number(Q) << endl << endl;

            }
        }    
}

// manual classification

#define T_MERIDIAN 1 // do a meridional rotation
#define T_REVERSE 2 // reverse the knot
#define T_1 4 // do a slide move with parameter 1

// transform a knot K in lens space l with flags w1
// if simple = false, also BST_reduce the knot
Cknot * transform1(Cknot *K, int l, u8 w1, bool simple = false) {
    
    Cknot *Q = K->deepCopy();
    if (w1 & T_MERIDIAN) Q->meridional_rotation();
    if (w1 & T_REVERSE) Q->reverse();
    Q = Omega_4(l, Q, (w1 & T_1 ? 1 : 0) );
    BST_shrink(&Q, UNORIENTED, 80, 30);
    
    if (!simple) {
        BST_reducible(NULL, Q, ORIENTED, 30, 0x1001001, 3);
        Q = BST_minimal_knot;
    }
    
    return Q;
}

// transform a knot K in lens space l with flags w (do a transformation twice)
// if simple = false, also BST_reduce the knot
// TODO: clean up, so only one transformation function is used
Cknot * transform2(Cknot *K, int l, u8 w1, u8 w2, bool simple = false) {
    
    Cknot *Q = K->deepCopy();
    if (w1 & T_MERIDIAN) Q->meridional_rotation();
    if (w1 & T_REVERSE) Q->reverse();
    Q = Omega_4(l, Q, (w1 & T_1 ? 1 : 0) );
    BST_shrink(&Q, UNORIENTED, 80, 30);
    
    if (!simple) {
        BST_reducible(NULL, Q, ORIENTED, 30, 0x1001001, 3);
        Q = BST_minimal_knot;
    }
    
    if (w2 & T_MERIDIAN) Q->meridional_rotation();
    if (w2 & T_REVERSE) Q->reverse();
    Q = Omega_4(l, Q, (w2 & T_1 ? 1 : 0) );
    BST_shrink(&Q, UNORIENTED, 80, 30);
    
    if (!simple) {
        BST_reducible(NULL, Q, ORIENTED, 30, 0x1001001, 3);
        Q = BST_minimal_knot;
    }
    
    if (!simple) {
        BST_reducible(NULL, Q, ORIENTED, 30, 0x1001001, 3);
        Q = BST_minimal_knot;
    }
    return Q;
}

// transform a knot K in lens space l with flags w (do a transformation three times)
// if simple = false, also BST_reduce the knot
// TODO: clean up, so only one transformation fucntion is used
Cknot * transform3(Cknot *K, int l, u8 w1, u8 w2, u8 w3, bool simple = false) {
    
    Cknot *Q = K->deepCopy();
    if (w1 & T_MERIDIAN) Q->meridional_rotation();
    if (w1 & T_REVERSE) Q->reverse();
    Q = Omega_4(l, Q, (w1 & T_1 ? 1 : 0) );
    BST_shrink(&Q, UNORIENTED, 80, 30);
    
    if (!simple) {
        BST_reducible(NULL, Q, ORIENTED, 30, 0x1001001, 3);
        Q = BST_minimal_knot;
    }
    
    if (w2 & T_MERIDIAN) Q->meridional_rotation();
    if (w2 & T_REVERSE) Q->reverse();
    Q = Omega_4(l, Q, (w2 & T_1 ? 1 : 0) );
    BST_shrink(&Q, UNORIENTED, 80, 30);
    
    if (simple) {
        BST_reducible(NULL, Q, ORIENTED, 30, 0x1001001, 3);
        Q = BST_minimal_knot;
    }
    
    if (w3 & T_MERIDIAN) Q->meridional_rotation();
    if (w3 & T_REVERSE) Q->reverse();
    Q = Omega_4(l, Q, (w3 & T_1 ? 1 : 0) );
    BST_shrink(&Q, UNORIENTED, 80, 30);
    
    if (!simple) {
        BST_reducible(NULL, Q, ORIENTED, 30, 0x1001001, 3);
        Q = BST_minimal_knot;
    }
    return Q;
}

// reduce the knot Q to K by any means (repeat until the equality is found)
// use only when we are sure the knot Q reduces to K
bool find_equality_by_any_cost(Cknot *K, Cknot *Q) {
    
    Cknot *K_  = K->deepCopy(), *Q_ = Q->deepCopy();
    Cknot *K__ = new Cknot(), *Q__ = new Cknot();
    
    while (!equal_knots(Q_, K_)) {
        if ((equal_knots(Q_, Q__)) && (equal_knots(K_, K__))) throw 20;
        if (!equal_knots(Q_, Q__)) {
            Q__ = Q_->deepCopy();
            BST_reducible(NULL, Q_, UNORIENTED, 30, 0x1001001001, 3);
            Q_ = BST_minimal_knot;
            Q_->canonical(UNORIENTED);
        }
        cout << "." << flush;
        if (!equal_knots(K_, K__)) {
            K__ = K_->deepCopy();
            BST_reducible(NULL, K_, UNORIENTED, 30, 0x1001001001, 3);
            K_ = BST_minimal_knot;
            K_->canonical(UNORIENTED);
        }
    }
    return YES;
}

// try determining eqaliies of knots k_id_1 = k_id_2 in lens space l
// m - index of method used (each method corresponds to different combination of slide moves)
bool try_classification_method(int m, int l, int k_id_1, int k_id_2) {
    
#define ph //cout << endl << (int)Q1->n << "," << (int)Q2->n << "; " << HOMFLY(Q1->deepCopy()) << endl << HOMFLY(Q2->deepCopy()) << endl;
    
    Cknot *Q1, *Q2;
    bool b_equality = NO;
    
    if (k_id_2 < k_id_1) throw 20; // should be ordered
    
    cout << lpq_s[l] << ": " <<  k_id_2 <<flush;
    
    switch (m) {
        case 1: // classification method
            Q1 = transform2(get_knot_by_id(k_id_1), l, 0, 0, YES);
            Q2 = get_knot_by_id(k_id_2);
            b_equality = equal_knots(Q1, Q2);
            break;
            
        case 2:
            Q1 = transform1(get_knot_by_id(k_id_1), l, 0, YES); //Q1->reverse();
            Q2 = transform1(get_knot_by_id(k_id_2), l, T_MERIDIAN, YES);
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 21:
            Q1 = transform1(get_knot_by_id(k_id_1), l, 0, YES); Q1->reverse();
            Q2 = transform1(get_knot_by_id(k_id_2), l, T_MERIDIAN, YES);
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 3:
            Q1 = transform2(get_knot_by_id(k_id_1), l, 0, T_MERIDIAN, YES);
            Q2 = get_knot_by_id(k_id_2);
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 4:
            Q1 = transform1(get_knot_by_id(k_id_1), l, 0,  YES);
            Q2 = transform1(get_knot_by_id(k_id_2), l, 0, YES); 
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 5:
            Q1 = transform1(get_knot_by_id(k_id_1), l, T_MERIDIAN, YES); //Q1->reverse();
            Q2 = transform1(get_knot_by_id(k_id_2), l, T_MERIDIAN, YES); 
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 51:
            Q1 = transform3(get_knot_by_id(k_id_1), l, 0, T_MERIDIAN, 0, YES);
            Q2 = transform3(get_knot_by_id(k_id_2), l, 0, T_MERIDIAN, 0, YES); 
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
        case 6:
            Q1 = transform2(get_knot_by_id(k_id_1), l, 0,T_REVERSE, NO);
            Q2 = transform2(get_knot_by_id(k_id_2), l, T_MERIDIAN, 0, NO); 
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 7:
            Q1 = transform2(get_knot_by_id(k_id_1), l, T_MERIDIAN, T_MERIDIAN, YES); Q1->meridional_rotation();
            Q2 = transform2(get_knot_by_id(k_id_2), l, T_MERIDIAN , T_MERIDIAN, YES); Q2->meridional_rotation(); 
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 8:
            Q1 = transform1(get_knot_by_id(k_id_1), l, T_MERIDIAN, YES); Q1->meridional_rotation();
            Q2 = get_knot_by_id(k_id_2);
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 9:
            Q1 = transform1(get_knot_by_id(k_id_1), l, T_MERIDIAN, YES);
            Q2 = get_knot_by_id(k_id_2);
            b_equality = find_equality_by_any_cost(Q1, Q2);
            break;
            
        case 10:
            Q1 = get_knot_by_id(k_id_1);
            Q2 = transform1(get_knot_by_id(k_id_2), l , T_MERIDIAN + T_REVERSE, YES);
            b_equality = equal_knots(Q1, Q2);
            break;
            
    }
    
    classification_lpq[l][k_id_2].duplicate = b_equality;
    cout << (b_equality ? " *" : " -") << endl;
    return b_equality;
}

// manualy find the R-moves of knots in L(p,q)
// used for knots where the ordinary classification fails

void manual_classification_lpq() {
    
    u16 prop;
   
    
    // L(2,1)
    
    try_classification_method(2, stolpq("L(2,1)"), 6, 63);
    try_classification_method(2, stolpq("L(2,1)"), 6, 307);
    try_classification_method(1, stolpq("L(2,1)"), 14, 367);
    try_classification_method(1, stolpq("L(2,1)"), 60, 352);
    
    // L(3,1)
    try_classification_method(21, stolpq("L(3,1)"), 14, 37);
    try_classification_method(2, stolpq("L(3,1)"), 34, 60);
    try_classification_method(3, stolpq("L(3,1)"), 47, 190);
    try_classification_method(4, stolpq("L(3,1)"), 160, 195);
    try_classification_method(5, stolpq("L(3,1)"), 160, 328);
    try_classification_method(51, stolpq("L(3,1)"), 197, 247);
    
    // L(4,1)
    try_classification_method(6, stolpq("L(4,1)"), 54, 213);
    try_classification_method(7, stolpq("L(4,1)"), 216, 294);
    
    // L(5,2)
    try_classification_method(8, stolpq("L(5,2)"), 299, 302);

    // L(6,1)
    try_classification_method(9, stolpq("L(6,1)"), 23, 120);
    
    // L(7,1)
    try_classification_method(9, stolpq("L(7,1)"), 65, 120);
    
    // L(7,2)
    try_classification_method(9, stolpq("L(7,2)"), 6, 18);
    try_classification_method(10, stolpq("L(7,2)"), 69, 299);
    
    // L(9,2)
    try_classification_method(10, stolpq("L(9,2)"), 18, 69);
    
    // L(11,2)
    try_classification_method(10, stolpq("L(11,2)"), 69, 338);
    
    cout << endl;
    
    // remove affine knots and connected sums, checked by hand
    
    cout << endl << endl << "Also removing knots: " << endl;
    for (int l = 0; l < N_HOMEO; l++) {
        
        cout << lpq_s[l] << (lens_homeo[l][0] > 9 ? "" :" ")<<": ";
        
        for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
        
            if (classification_lpq[l][(*iK)->ID].duplicate) continue;
            
            prop = get_knot_property(*iK, l);
            
            if (prop & PROPERTY_AFFINE) {
                cout << (int)(*iK)->ID << " (affine), ";
                classification_lpq[l][(*iK)->ID].duplicate = true;
            }
            
            if (prop & PROPERTY_SUM) {
                cout << (int)(*iK)->ID << " (composite), ";
                classification_lpq[l][(*iK)->ID].duplicate = true;
            }
        }
        cout << endl;
    }
    
    
    /* obsolete
    
    cout << lpq_s[l = stolpq("L(2,1)")] << endl; // equal knots: [20, 560] ok, [95, 541] ok, [230, 231], [233, 234]
    
    // knots 20, 560
    K = transform2(get_knot_by_id(20), l, 0, 0, YES);
    Q = get_knot_by_id(560);
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 95, 541
    K = transform2(get_knot_by_id(95), l, 0, 0, YES);
    Q = get_knot_by_id(541);
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // L(3,1)
    cout << lpq_s[l = stolpq("L(3,1)")] << endl; // equal knots: [20, 52], [47, 95] ok, [72, 253], [221, 262, 507], [230, 231], [233, 234], [264, 362]
    
    // knots 20, 52
    
    K = transform1(get_knot_by_id(20), l, 0, YES); K->reverse();
    Q = transform1(get_knot_by_id(52), l, T_MERIDIAN, YES); Q->ID = 52;
    if (find_equality_by_any_cost(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 47, 95
    K = transform1(get_knot_by_id(47), l, 0, YES);
    Q = transform1(get_knot_by_id(95), l, T_MERIDIAN, YES); Q->ID = 95;
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 72, 253

    K = transform2(get_knot_by_id(72), l, 0, T_MERIDIAN, YES);
    Q = get_knot_by_id(253);
    if (find_equality_by_any_cost(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 221, 262,
    
    K = transform1(get_knot_by_id(221), l, 0,  YES);
    Q = transform1(get_knot_by_id(262), l,  0, YES); Q->ID = 262;
    if (find_equality_by_any_cost(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 221, 507
    
    K = transform1(get_knot_by_id(221), l, T_MERIDIAN, YES); K->reverse();
    Q = transform1(get_knot_by_id(507), l, T_MERIDIAN, YES); Q->ID = 507;
    //cout << endl << "(" << winding_number(K) <<") " <<HOMFLY(K)<< endl <<  "(" << winding_number(Q) <<") " << HOMFLY(Q) << endl;
    if (find_equality_by_any_cost(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 264, 362
    K = transform3(get_knot_by_id(264), l, 0, T_MERIDIAN, 0, YES);
    Q = transform3(get_knot_by_id(362), l, 0, T_MERIDIAN, 0, YES); Q->ID = 362;
    
    if (find_equality_by_any_cost(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    
    cout << lpq_s[l = stolpq("L(4,1)")] << endl; // equal knots: [83, 295], [230, 231], [233, 234], [301, 447]
    
    // knots 83, 295
    
    K = transform2(get_knot_by_id(83), l, 0,T_REVERSE, NO);
    Q = transform2(get_knot_by_id(295), l, T_MERIDIAN, 0, NO); Q->ID = 295;
    if (find_equality_by_any_cost(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 301, 447
    
    K = transform2(get_knot_by_id(301), l, T_MERIDIAN, T_MERIDIAN, YES); K->meridional_rotation(); K->reverse();
    Q = transform2(get_knot_by_id(447), l, T_MERIDIAN , T_MERIDIAN, YES); Q->meridional_rotation(); Q->ID = 447;
    if (find_equality_by_any_cost(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    //cout << lpq_s[l = stolpq("L(5,1)")] << endl; // equal knots: [230, 231], [233, 234]
    
    
    cout << lpq_s[l = stolpq("L(5,2)")] << endl; // equal knots: [25, 101] ok, [230, 231], [233, 234], [456, 459] ok
    
    // knots 25, 101
    
    K = transform1(get_knot_by_id(25), l, T_MERIDIAN, YES);
    K->meridional_rotation(); BST_shrink(&K, UNORIENTED, 50, 20);
    Q = get_knot_by_id(101);
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 456, 459
    
    K = transform1(get_knot_by_id(456), l, T_MERIDIAN, YES);
    K->meridional_rotation(); K->reverse(); BST_shrink(&K, UNORIENTED, 50, 20);
    Q = get_knot_by_id(459);
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    cout << lpq_s[l = stolpq("L(6,1)")] << endl; // equal knots: [30, 145] ok, [99, 141] ok, [230, 231], [233, 234]
    
    // knots 30, 145
    
    K = transform1(get_knot_by_id(30), l, T_MERIDIAN, YES);
    Q = get_knot_by_id(145);
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 99, 141
    
    K = transform1(get_knot_by_id(99), l, T_MERIDIAN, YES);
    Q = get_knot_by_id(141);
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    
    cout << lpq_s[l = stolpq("L(7,1)")] << endl; // equal knots: [101, 145] ok, [230, 231], [233, 234]
    
    // knots 101, 145
    
    K = transform1(get_knot_by_id(101), l, T_MERIDIAN, YES);
    Q = get_knot_by_id(145);
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    
    cout << lpq_s[l = stolpq("L(7,2)")] << endl; // equal knots: [7, 25] ok, [109, 456] ok, [230, 231], [233, 234]
    
    // knots 7, 25
    
    K = transform1(get_knot_by_id(7), l, T_MERIDIAN, YES);
    Q = get_knot_by_id(25);
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    // knots 109, 456
    
    K = get_knot_by_id(109);
    Q = transform1(get_knot_by_id(456), l , T_MERIDIAN + T_REVERSE, YES); Q->ID = 456;
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    //cout << lpq_s[l = stolpq("L(8,1)")] << endl; // equal knots: [230, 231], [233, 234]
    
    //cout << lpq_s[l = stolpq("L(8,3)")] << endl; // equal knots KBSM: [101, 139], [105, 125], [230, 231], [233, 234], [459, 475]
    
    //cout << lpq_s[l = stolpq("L(9,1)")] << endl; // equal knots: [230, 231], [233, 234]
    
    cout << lpq_s[l = stolpq("L(9,2)")] << endl; // equal knots: [25, 109] ok, [230, 231], [233, 234]
    
    // knots 25, 109
    
    K = get_knot_by_id(25);
    Q = transform1(get_knot_by_id(109), l , T_MERIDIAN + T_REVERSE, YES); Q->ID = 109;
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    //cout << lpq_s[l = stolpq("L(10,1)")] << endl; // equal knots: [230, 231], [233, 234]
    
    //cout << lpq_s[l = stolpq("L(10,3)")] << endl; // equal knots KBSM: [101, 139], [105, 125], [230, 231], [233, 234], [459, 475]
    
    //cout << lpq_s[l = stolpq("L(11,1)")] << endl; // equal knots: [230, 231], [233, 234]
    
    cout << lpq_s[l = stolpq("L(11,2)")] << endl; // equal knots: [109, 527] ok, [230, 231], [233, 234]
    
    // knots 109, 527
    
    K = get_knot_by_id(109);
    Q = transform1(get_knot_by_id(527), l , T_MERIDIAN + T_REVERSE, YES); Q->ID = 527;
    if (equal_knots(K, Q)) { classification_lpq[l][Q->ID].duplicate = YES; cout << "Knot " << Q->ID << " reducible." << endl; }
    
    //cout << lpq_s[l = stolpq("L(11,3)")] << endl; // equal knots KBSM: [101, 139], [105, 125], [230, 231], [233, 234], [459, 475]
    
    //cout << lpq_s[l = stolpq("L(12,1)")] << endl; // equal knots: [230, 231], [233, 234]
    
    //cout << lpq_s[l = stolpq("L(12,5)")] << endl; // equal knots KBSM: [101, 139], [105, 125], [230, 231], [233, 234], [459, 475]
   
    // remove affine knots and connected-sum knots from the classification (all of the candidates are indeed affines/composites, checked by hand)
    u16 property;
    init_prime_affine();
    
    // affines
    cout << endl <<"Affine knots: ";
    for (int l = 0; l < N_HOMEO; l++) {
        cout << endl <<lpq_s[l] << ": " << flush;
        for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
            K = (*iK)->deepCopy();
            if (classification_lpq[l][K->ID].duplicate) continue;
            property = get_knot_property_lpq(l, K);
            if (property & PROPERTY_AFFINE) {
                cout << K->ID << ", ";
                classification_lpq[l][K->ID].duplicate = true;
            }
        }
    }
    
    cout << endl << endl << "Connected sums: ";
    for (int l = 0; l < N_HOMEO; l++) {
        cout << endl <<lpq_s[l] << ": " << flush;
        for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
            K = (*iK)->deepCopy();
            if (classification_lpq[l][K->ID].duplicate) continue;
            property = get_knot_property_lpq(l, K);
            if (property & PROPERTY_SUM) {
                cout << K->ID << ", ";
                classification_lpq[l][K->ID].duplicate = true;
            }
        }
    }
    
    cout << endl;
     */
}

#endif
