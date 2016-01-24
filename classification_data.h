//  classification_data.h
//  Created by Boštjan on 4/22/13.
//
//  classification class used in the classification

#ifndef Lpq_classification_data_h
#define Lpq_classification_data_h

#include "knot.h"
#include "global_vars.h"

#define MAX_GROUPS 4000 // maximal number of groups (partitions)
#define MAX_GROUP_SIZE 100 // maximal knots in each group

// stores info of a knot
class Cknot_info {
public:
    bool prime;
    bool duplicate;
    bool composite; // obsolete
    
    Cknot_info() { /*ID = -1;*/ prime = false; duplicate = false; composite = false;  }
};

// classification of knots in the solid torus
Cknot_info classification[MAX_KNOTS];

// classification of knots in L(p,q)
Cknot_info classification_lpq[N_HOMEO][MAX_KNOTS];


// minimal representative of the knot class (for each knot store the ID of the first knot in the class
int classification_equality_lpq[N_HOMEO][MAX_KNOTS];

// flags
#define PRIME_KNOT 1
#define MIRROR_KNOT 2
#define AMPHICHIRAL_KNOT 3


int n_prime; // number of prime knots
string knot_name[1000]; // names of knots "0_1", "1_1",...
int knot_name_ID[1000]; // ID of the ith knot in the classification table
int chirality[1000]; // 1 - prime, 2 mirror of prime, 3 - amphichiral
int chirality_id[1000]; // id of the mirror knot of a knot


void get_chirality(t_knot_set *KNOTS, t_group_list *GROUP);

// counts non-duplicate knots from KNOTS wrt the classification table ct with n crossings
int count_prime_knots(t_knot_set * KNOTS, Cknot_info ct[MAX_KNOTS], int n) {
    int count = 0;
    for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
        if ((*iK)->n != n) continue; // select only knots witn n crossings
        if (ct[ (*iK)->ID].duplicate) continue; // select only non-duplicate knots
        count++;
    }
    return count;
}

// generate the knots names "0_1", "1_1",...
// TODO: optimize
void generate_knot_names(bool ignore_mirror = true) {

    t_knot_set::iterator iK;
    int n, k;

    
    get_chirality(KNOTS, G_HOMFLY);
    
    // unknot
    knot_name[0] = "0_1";
    knot_name_ID[0] = 0;
    n_prime = 1;
    n = 0;
    k = 2;
    
    for (iK = KNOTS->begin(); iK != KNOTS->end(); iK++) { // loop through all knots
        
        if (!classification[(*iK)->ID].duplicate) { // only consider knots that are not duplicates
            
            if ((chirality[(*iK)->ID] != 2)||(!ignore_mirror)) { // ignore mirrors?
                
                if (ignore_mirror)
                    if (( chirality[(*iK)->ID] != 1 ) && (chirality[(*iK)->ID] != 3)) throw 20; // if not chiral or amphichiral, throw exception
                
                if ((*iK)->n > n) {n = (*iK)->n; k = 1;}
                
                knot_name_ID[n_prime] = (*iK)->ID;
                
                if (chirality[(*iK)->ID] != 2)
                    knot_name[n_prime] = string(itoa(n)) + "_{" + string(itoa(k)) + "}";
                
                if (chirality[(*iK)->ID] != 2) n_prime++;
                if (chirality[(*iK)->ID] != 2) k++;
            }
        }
    }
}

// get knot name index from ID
int name_from_id(int ID) {
    int i = 0;
    while (knot_name_ID[i] != ID) i++;
    return i;
}

// print the names of the knotss
void print_knot_names() {
    for (int i = 0; i < n_prime; i++) {
        if (i % 10 == 0) cout << endl;
        cout << knot_name[i] << " = " << knot_name_ID[i] << ", " << flush;
    }
    cout << endl << endl;
}

// generate knot names by respecting chiral and amphichiral knots (naming the mirror paris by a bar)
void generate_knot_names_bar() {
    
    t_knot_set::iterator iK;
    int n, k;
    
    get_chirality(KNOTS, G_HOMFLY);
    
    knot_name[0] = "0_1"; // unknot
    knot_name_ID[0] = 0;
    n_prime = 1;
    n = 0;
    k = 2;
    
    for (iK = KNOTS->begin(); iK != KNOTS->end(); iK++) { // loop through all knows
        if (!classification[(*iK)->ID].duplicate) { // prime knots
            
            knot_name_ID[n_prime] = (*iK)->ID; // save knot ID
            
            if (chirality[(*iK)->ID] != MIRROR_KNOT) {
                if ((*iK)->n > n) {n = (*iK)->n; k = 1;} // number of crossings
                knot_name[n_prime] = string(itoa(n)) + "_{" + string(itoa(k)) + "}"; // knot name "n_k"
                if (chirality[(*iK)->ID] == 3) knot_name[n_prime] = "*" + knot_name[n_prime]; // amphichiral?
            }
            else { // mirror of a lower index prime knot
                if ((*iK)->ID == 171) knot_name[n_prime] = "\\bar{5_{27}}"; else //manual
                knot_name[n_prime] = "\\bar{" + knot_name[ name_from_id(chirality_id[(*iK)->ID  ]) ] + "}";
            }
            n_prime++;
            if (chirality[(*iK)->ID] != 2) k++; // only increase index, if not a mirror
        }
    }

}

// obsolete
int mirror_groups[MAX_KNOTS];
int mirror_groups_lpq[N_HOMEO][MAX_KNOTS];

// clears the classification for the solid torus
void clear_classification_data() {
    for (int i=0; i < MAX_KNOTS; i++)
        classification[i].prime = classification[i].duplicate = classification[i].composite = false;
}

// clears the classification for L(p,q)
void clear_classification_data_lpq() {
    for (int l = 0; l < N_HOMEO; l++)
        for (int i=0; i < MAX_KNOTS; i++)
            classification_lpq[l][i].prime = classification_lpq[l][i].duplicate = classification_lpq[l][i].composite = false;
}

// transfer the "duplicate" flag from knots in the torus to knots in L(p,q)
void or_duplicates_from_torus() {
    for (int l = 0; l < N_HOMEO; l++)
        for (int i=0; i < MAX_KNOTS; i++)
            if (classification[i].duplicate) classification_lpq[l][i].duplicate = true;
}

// number of unclassified knots in the knot group KS (torus)
int unclassified_in_group(t_knot_set * KS, Cknot_info * classification) {
    int unclassified = 0;
    for (t_knot_set::iterator iK = KS->begin(); iK != KS->end(); iK++) {
        if ((!classification[(*iK)->ID].prime) && (!classification[(*iK)->ID].duplicate) && (!classification[(*iK)->ID].composite))
            unclassified++;
    }
    return unclassified;
}

// number of classified knots in the knot group KS (torus)
int classified_in_group(t_knot_set * KS, Cknot_info * classification) {
    int classified = 0;
    for (t_knot_set::iterator iK = KS->begin(); iK != KS->end(); iK++) {
        if ((!classification[(*iK)->ID].duplicate)) classified++;
    }
    return classified;
}

// number of unclassified knots in the knot group KS  (in lens space l)
int unclassified_in_group_lpq(int l, t_knot_set * KS) {
    int unclassified = 0;
    for (t_knot_set::iterator iK = KS->begin(); iK != KS->end(); iK++) {
        if ((!classification_lpq[l][(*iK)->ID].prime) && (!classification_lpq[l][(*iK)->ID].duplicate) && (!classification_lpq[l][(*iK)->ID].composite))
            unclassified++;
    }
    return unclassified;
}

// print stats, obsolete
ostream &operator<<(ostream &out, Cknot_info * G) {
    int p, d, c, pd, pc, dc, pdc, n, a;
    p = d = c = pd = pc = dc = pdc = n = a = 0;

    // primes
    for (int i=1;i<=number_of_knots;i++) {
        a++;
        if ((!classification[i].prime) && (!classification[i].duplicate) && (!classification[i].composite)) n++;
        if ((classification[i].prime) && (!classification[i].duplicate) && (!classification[i].composite)) p++;
        if ((!classification[i].prime) && (classification[i].duplicate) && (!classification[i].composite)) d++;
        if ((!classification[i].prime) && (!classification[i].duplicate) && (classification[i].composite)) c++;
        if ((classification[i].prime) && (classification[i].duplicate) && (!classification[i].composite)) pd++;
        if ((classification[i].prime) && (!classification[i].duplicate) && (classification[i].composite)) pc++;
        if ((!classification[i].prime) && (classification[i].duplicate) && (classification[i].composite)) dc++;
        if ((classification[i].prime) && (classification[i].duplicate) && (classification[i].composite)) pdc++;
    }
    
    out << "Knots: " << a << " = " << p << " primes + " << d << " duplicates + " << c << " composites. " << endl;
    out << "Unknown: " << n << endl;
    out << "Prime+duplicate: " << pd << " Prime+comsposite: " << pc << " Duplicate+composite: " << dc << " Prime+duplicate+composite: " << pdc << endl;
    out << (100.0*(a - n)/a) << "% classified.";
    out << endl << endl;
    return out;
}
            
#endif
