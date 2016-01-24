//  global_vars.h
//  Created by Boštjan on 4/22/13.
//
//  global variables and definitions used throughout the project space
//
//  TODO: seperate types and variables in common.h and global_vats.h to make more sense

#ifndef Lpq_global_vars_h
#define Lpq_global_vars_h

#include <set>
#include <list>
#include <vector>
#include "common.h"

// path to files (loading/saving)

#define LOG_FILE "./lpq.log"
#define GENERATORS_PATH "./lpq-homfly-generators.txt"
#define GENERATORS_PATH_INTERMEDIATE "./lpq-homfly-intermediate-generators.txt"
#define HOMFLY_PATH "./torus-HOMFLY.txt"

#define SAVED_KNOTS_PATH "./torus-gauss-codes.txt"
#define CLASSIFICATION_PATH "./torus-classification.txt"
#define CLASSIFICATION_BACKUP_PATH "./backup/torus-classification"

#define CLASSIFICATION_PATH_LPQ "./lpq-classification.txt"
#define CLASSIFICATION_BACKUP_PATH_LPQ "./backup/lpq-classification"

// knot property flags
#define PROPERTY_SUM 1
#define PROPERTY_AFFINE 2
#define PROPERTY_PROJECTIVE 4

// loading flags
#define LOAD_GC 1 // load Gauss codes
#define LOAD_TORUS 2 // load solid torus classification
#define CALC_HSM 4 // calculate HSM of solid torus
#define CALC_KBSM 8 // calculate KBSM of solid torus
#define CALC_HSM_LPQ 16 // load and calculate HSM generators
#define CALC_KBSM_LPQ 32 // calculate KBSM of L(p,q)
#define LOAD_LPQ 64 // load L(p,q) classification

#define SAVE_GC 128 // save Gauss codes
#define SAVE_TORUS 256 // save solid torus classification
#define SAVE_LPQ 512 // save L(p,q) classification
#define CLEAN_UP 1024 // free memory
#define QUIT 2048 // quit program
#define CLEAN_QUIT (CLEAN_UP + QUIT) // free memory and quit

#define OVER_LIMIT     (2300000000) // memory usage to stop complicating the knot (~2.3 GB)
#define CRITICAL_LIMIT (3500000000) // knot search shutdown memory usage (~3.5 GB)

#define MAX_HSM_GEN 250 // number of generators of the HSM of the solid torus for up to 6 crossings, obsolete
#define MAX_KBSM_GEN 7 // number of generators of the KBSM of the solid torus for up to 6 crossings, obsolete

#define N_HOMEO 19 // number of lens spaces up to p = 12

// non-homeomorphic lens spaces up to p = 12
int lens_homeo[N_HOMEO][2] = {{2,1}, {3,1}, {4,1}, {5,1}, {5,2}, {6,1}, {7,1}, {7,2}, {8,1}, {8,3}, {9,1}, {9,2}, {10,1}, {10,3}, {11,1}, {11,2}, {11,3}, {12,1}, {12,5}};

// non-homeomorphic lens spaces up to p = 12 (strings)
string lpq_s[N_HOMEO] = {"L(2,1)", "L(3,1)", "L(4,1)", "L(5,1)", "L(5,2)", "L(6,1)", "L(7,1)", "L(7,2)", "L(8,1)", "L(8,3)", "L(9,1)",
    "L(9,2)", "L(10,1)", "L(10,3)", "L(11,1)", "L(11,2)", "L(11,3)", "L(12,1)", "L(12,5)"};//, "L(13,1)", "L(13,2)", "L(13,3)", "L(13,5)", "L(14,1)", "L(14,3)"};

// which lens spaces do we consider? obsolete
//bool consider_lpq[N_HOMEO] = {YES,YES,YES,YES,YES,YES,YES,YES,YES,YES, YES,YES,YES,YES,YES,YES,YES,YES,YES};

class CknotGroup;
class Cknot;
struct compareKnots;
class Cknot_info;
class Cmultivariate;

// list of knots
typedef std::list<CknotGroup *> t_group_list;
// oriented set of knots (BST)
typedef std::set<Cknot *,compareKnots> t_knot_set;
// vectors of polynomials
typedef std::vector<Cmultivariate *> t_poly_vector;

void generateAllKnots(int N_MIN, int N_MAX, t_knot_set *KNOTS, u64 flags);
bool BST_reducible(Cknot *compare_knot_, Cknot *K_, bool oriented, int max_depth, u64 b_increase, int max_additional_crossings);
bool BST_shrinkable(Cknot *, bool, int, int);
void BST_shrink(Cknot **, bool, int, int);
Cknot * BST_shrink(Cknot *, bool, int, int);
void saveClassification(string s);
void OR_loadClassification(string s);
void save_classification_lpq(string s);
Cknot * Omega_4(int lpqi, Cknot *K_, s_small extra = 0);
int estimate_number_of_crossings(int l, Cknot * K);
Cmultivariate *HSM_lpq(int l, Cmultivariate * poly_);
u16 get_knot_property(Cknot *K, int l = -1);

Cknot *BST_last_knot;
Cknot *BST_minimal_knot;

#define number_of_knots (int)(KNOTS->size())
#define number_of_groups (int)(G_HOMFLY->size())

#define mv(x) new Cmultivariate(x)

// which invaraint do we consider using for partitionig knots
#define GROUP_OF_CHOICE(l) (lens_homeo[l][1] < 3 ? G_HOMFLY_LPQ : G_KBSM_LPQ) // <= 1

// global variables and global lists

// list of all knots
t_knot_set * KNOTS = new t_knot_set();
// list of HOMLFY groups
t_group_list * G_HOMFLY = new t_group_list();
// list of HOMLFY groups
t_group_list * G_KBSM = new t_group_list();

#define tgl new t_group_list()
#define tgl19 tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl, tgl

// list of HOMLFY groups
t_group_list * G_HOMFLY_LPQ[N_HOMEO] = {tgl19};
// list of KBSM groups
t_group_list * G_KBSM_LPQ[N_HOMEO] = {tgl19}; 

// store HSMs of torus
t_poly_vector HOMFLYS;
t_poly_vector HOMFLYS_reverse;
// store KBSM of torus
t_poly_vector KBSMS;

// store HSMs of lpq
t_poly_vector HOMFLYS_lpq[N_HOMEO];
t_poly_vector HOMFLYS_lpq_reverse[N_HOMEO];
// store KBSMs of Lpq
t_poly_vector KBSMS_lpq[N_HOMEO];

// global variable used by BST_reduce (stores the minimal knot when freeing memory after BST_reduction())
Cknot *reduced_knot;



#endif
