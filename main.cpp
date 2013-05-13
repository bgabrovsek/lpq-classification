//  main.cpp
//
//  Created by Bostjan Gabrovsek 31/1/13, bostjan.gabrovsek@fmf.uni-lj.si
// 
//  Classification of knots in lens spaces

/*
Maximum number of crossings
If MAX < 31, 64 bit architecture will be used,
if MAX >= 32, 128 bit architecture will be used,
if MAX >= 64, 256 bit architecture will be used,
if MAX >= 128, 512 bit architecture will be used.

if not sure, use the upper limit provided (e.g. 63 for 128 architecture) or use MAX = 220 (this may take more time).

Command line options:

-g
 
Generates all ralizable Gauss codes up to N crossings, saves the GC to file,
calculates the HOMLFY skein modules partitions,
sets first knot in each partition as prime
stores the classification to file
(for this option code should be compiled with MAX >= 12).
 
-c
 
Tries to reduce each non-prime knot using R-moves, flypes, and meridional rotations
(for this option code should be compiled with MAX >= 12).
 
-h
Calculates HOMLFY skein module generators in terms of solid torus HOMLFY skein module generators
(for this option code should be compiled with MAX >= 220).
 
-i
Initialized classification of knots in L(p,q)
(for this option code should be compiled with MAX >= 12).
 
-l
Tries to reduce non-prime knots in L(p,q) using R-moves, flypes, and meridional rotations
(for this option code should be compiled with MAX >= 127).
 
-m
Manualy reduce L(p,q) knots by R-moves, flypes, meridional rotations, and slides
(for this option code should be compiled with MAX => 63).
 
-pt
Prints the classification table(s).

-pf
Prints knots that share HSM or KBSM
 
-ph
Prints HSM of knots

-pm
Prints KBSM of knots
*/

#define MAX 31

#include "knot_preprocessor_directives.h"

#include <iostream>
#include <string.h> 
#include <set>
#include <map>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "knot.h"
#include "global_vars.h"
#include "classification_data.h"
#include "multivariate_laurent.h"
#include "timing.h"
#include "in_out.h"
#include "KBSM.h"
#include "reidemeister_moves.h"
#include "algorithms.h"
#include "HOMFLY.h"
#include "timing.h"
#include "lens.h"
#include "link.h"
#include "tree.h"
#include "randomized.h"
//#include "numbers.h"

#define done cout << " DONE." << endl
#define dot cout << "." << flush

// global variables

char argument = ' '; // primary argument from command line
char sub_argument = ' '; // secondary argument in command line

// print help when run without arguments
void print_help() {
    cout << "usage: Lpq [-gcihklmwp(-)r]" << endl << endl;
    cout << "  -g : generates extended Gauss codes in the solid torus and save, clear classification data (64 bit)" << endl;
    cout << "  -c : classifies knots in the solid torus (64 bit)" << endl;
    cout << "  -i : initializes lpq classification, marks duplicates from torus, set first knot in each group as prime (64 bit)" << endl;
    cout << "  -h : calculates the lens space HOMFLY skein module generators (512 bit)" << endl;
    cout << "  -k : print KBSM relations in L(p,q) in terms of generators of the solid torus (64 bit)" << endl;
    cout << "  -l : classify knots in lens spaces (256 bit)" << endl;
    cout << "  -m : manual classification of knots in lens spaces, i.e. adjusted methods for unclassified knots (128 bit)" << endl;
    // cout << "  LN: classify lens spaces with strength N by alternative method" << endl;
    cout << "  -w : counts Gauss words (64 bit)" << endl;
    cout << "  -pf : print knots that share HSM or KBSM (64 bit)" << endl;
    cout << "  -pt : print classification table(s)" << endl;
    cout << "  -ph : print HSM of knots" << endl;
    cout << "  -pm : print KBSM of knots" << endl;
    //cout << "  a : print amphichiral candidates" << endl;
    cout << "  -r : uses a randomized algorithm to try finding equivalences of unclassified knot(s) (128 bit)" << endl;
    cout << endl;
    cout << "The parenthesis indicates the minimal architecture the code must be compiled in for that specific option." << endl;
    cout << "Only one input argument can be used at a time." << endl;
    cout << endl;
    cout << "For a complete classification of knots in the solid torus and L(p,q), run the program with the following sequence of arguments:" << endl;
    cout << "\"-g\", \"-c\", \"-h\",  \"-i\", \"-k\", \"-l\",  \"-m\", \"-pt\"" << "." << endl;
}

void make_operation(u32 flags) { // load/save Guass codes; calculate KBSM, HSM; load/save/clear classification
    
    if (flags & LOAD_GC) {
        cout << "Loading Gauss codes..." << flush;
        loadKnotsFromFile(KNOTS, SAVED_KNOTS_PATH); done;
    }
    
    if (flags & LOAD_TORUS) {
        cout << "Loading classification of solid torus..." << flush;
        loadClassification(CLASSIFICATION_PATH); done;
    }
    
    if (flags & CALC_HSM) {
        cout << "Calculating HOMFLY skein modules." << flush;
        calculate_HSMs(KNOTS); dot;
        groupKnotsByHSM(G_HOMFLY, KNOTS); dot; done;
    }
    
    if (flags & CALC_KBSM) {
        cout << "Calculating Kauffman bracket skein modules." << flush;
        calculate_KBSMs(KNOTS); dot;
        groupKnotsByKBSM(G_KBSM, KNOTS); done;
    }
    
    if (flags & CALC_HSM_LPQ) {
        cout << "Loading and calculating L(p,q) HOMLFY skein modules" << flush;
        load_generators(GENERATORS_PATH); dot;
        //process_generators_iteratively();
        calculate_HSMs_lpq(KNOTS); dot;
        groupKnotsByHSM_lpq(G_HOMFLY_LPQ, KNOTS); dot; done;
    }
    
    if (flags & CALC_KBSM_LPQ) {
        cout << "Calculating L(p,q) Kauffman bracket skein modules" << flush;
        calculate_KBSM_generators(); dot;
        calculate_KBSMs_lpq(KNOTS); dot;
        groupKnotsByKBSM_lpq(G_KBSM_LPQ, KNOTS); dot; done;
    }
    
    if (flags & LOAD_LPQ) {
        cout << "Loading L(p,q) classification..." << flush;
        load_classification_lpq(CLASSIFICATION_PATH_LPQ); done;
    }
    
    if (flags & SAVE_GC) {
        saveKnotsToFile(KNOTS, SAVED_KNOTS_PATH);
        cout << "Gauss codes saved to " << SAVED_KNOTS_PATH << endl;
    }

    if (flags & SAVE_TORUS) {
        saveClassification(CLASSIFICATION_PATH);
        cout << "Solid torus classification saved to " << CLASSIFICATION_PATH << endl;
    }
    
    if (flags & SAVE_LPQ) {
        save_classification_lpq(CLASSIFICATION_PATH_LPQ);
        cout << "L(p,q) classification saved to " << CLASSIFICATION_PATH_LPQ << endl;
    }
    
    if (flags & CLEAN_UP) {
        deleteGroups(G_HOMFLY);
        deleteKnots(KNOTS);
        delete KNOTS;
    }
    
    if (flags & QUIT) exit(0);

}

int main(int argc, const char * argv[]) {
    
    // COMMAND LINE ARGUMENTS
    
    if (argc <= 1) { print_help(); exit(0); } // if no command line arguments, print help
    if (argc > 2) {cout << "Too many arguments." << endl; exit(0); } // if more then one command line argument
    argument = argv[1][1]; // primary argument
    if (argv[1][2]) sub_argument = argv[1][2]; // secondary argument
    
    // BACKUP
    
    // backup the classificatin file
    // backup_file(CLASSIFICATION_PATH, CLASSIFICATION_BACKUP_PATH);
    
    // INITAILIZATION
    
    // init the multivariable laurent variables
    initMultivariateLaurent(HOMFLY_BASE_VARIABLES[HOMFLY_MAX_VAR] + (string)"vzax"); // v,z for HOMFLY, a,x for KBSM
    
    GROUP_OUTPUT_BY_FIRST_VARIABLE = false; // do not group common terms
    HOMFLY_STYLE = true; // print HSM variables as "t_1", "t_2",...

    generate_torus_HSM_generator_table(); // generate the table of HSM generators
    init_prime_affine(); // calculate KBSM of prime affine knots (for checking composites)
    
    // ARGUMENT DEPENDANT
    
    // count Gauss codes, realizable Gauss codes, etc.
    
    if (argument == 'w') {
        
        for (int n = 1; n<=7; n++) {
            cout << "N = " << n << endl;
            count_gauss_words(n);
            cout << "gauss:      " << f(2*n)* (1<< n) << endl;
            cout << "realizable: " << count_realizable << endl;
            cout << "canonical:  " << count_canonical << endl;
            cout << "real ext.:  " << count_realizable_extended << endl;
            cout << "canon ext.: " << count_canonical_extended << endl;
        }
        exit(0);
    }
    
    // initialize the solid torus classification
    
    if (argument == 'g') {
        
        cout << "Generating realizable Gauss codes..." << flush;
        generate_knot_candidates(KNOTS); done;
        
        make_operation(SAVE_GC + CALC_HSM);
        
        cout << "Clearing solid torus classification" << flush;
        clear_classification_data(); dot;
        set_single_knots_to_classified(); dot;
        mark_composites(KNOTS); dot; done;
        
        make_operation(SAVE_TORUS + CLEAN_QUIT);
    }
    
    // initialize the L(p,q) classification
    
    if (argument == 'i') {
        
        make_operation(LOAD_GC + LOAD_TORUS + CALC_HSM + CALC_KBSM + CALC_HSM_LPQ + CALC_KBSM_LPQ);
        
        cout << "Clearing L(p,q) classificarion" << flush;
        clear_classification_data_lpq(); dot;
        or_duplicates_from_torus(); dot;
        set_first_knots_to_classified_lpq(); dot; done;
        
        make_operation(SAVE_LPQ + CLEAN_QUIT);
    }
    
    // print the KBSM relations of the generators in L(p,q), Mathematica format
    
    if (argument == 'k') {
        print_KBSM_relations();
        make_operation(QUIT);
    }
    
    // caclulate the HSM generetos of L(p,q) in terms of solid torus HSM generators
    
    if (argument == 'h') {
        cout << "Calculating HOMFLY skein module generators" << endl;
        lens_space_generators_in_terms_of_torus_generators(); done;
        make_operation(QUIT);
    }
    
    // classify knots in the solid torus
    
    if (argument == 'c') {
        
        make_operation(LOAD_GC + CALC_HSM + LOAD_TORUS);

        cout << "Reducing knots, reducible (*), unknown (-):" << endl;
        
        classify(15, b1, 2); // weaker
        
        cout << endl;
        
        classify(30, b1001001, 3); // stronger
        
        make_operation(CLEAN_QUIT); // classification is saved by classify()
    }
    
    // classify knots in L(p,q)
    
    if (argument == 'l') { // classify lens spaces
        
        make_operation(LOAD_GC + LOAD_TORUS + CALC_HSM + CALC_KBSM + CALC_HSM_LPQ + CALC_KBSM_LPQ + LOAD_LPQ);
        
        cout << "Reducing knots in L(p,q), reducible (*), unknown (-):" << endl;
        
        //classify_lpq(10, b0, 2, 5, 2);
        classify_lpq(30, b0, 2, 10, 3);
        
        //classify_lpq(40, b0000001001, 2, 25, 10);
        
        // the obsolete options below do not produce new results:
        // classify_lpq(35, b0001001001, 2, 20, 7);
        // classify_lpq(40, b0000001001, 2, 25, 10);
        // classify_lpq(40, b1001001001, 5, 25,10);
        
        make_operation(SAVE_LPQ + CLEAN_QUIT);
    }
    
    // manual classification of knots in L(p,q)
    // manualy deal with the exceptions that the above classification method fails
    
    if (argument == 'm') {
        
        make_operation(LOAD_GC + LOAD_TORUS + CALC_HSM + CALC_KBSM + CALC_HSM_LPQ + CALC_KBSM_LPQ + LOAD_LPQ);
        
        manual_classification_lpq();
        
        make_operation(SAVE_LPQ + CLEAN_QUIT);
    }
    
    // print results
    
    if (argument=='p') {
        
        make_operation(LOAD_GC + LOAD_TORUS + CALC_HSM + CALC_KBSM + CALC_HSM_LPQ + CALC_KBSM_LPQ + LOAD_LPQ);
        
        // print classification table
        if (sub_argument == 't') {  classification_table();}
        
        // print torus knot groups
        if (sub_argument == 'g') print_knots(G_HOMFLY, GROUP_POLY | PRINT_GROUP_POLY);
        
        // print torus knots
        if (sub_argument == 'k') print_knots(G_HOMFLY, PRINT_GROUPED_BY_LINE | PRINT_KNOT_ID | PRINT_GROUP_ID | NO_DUPLICATES);
        
        // print latex knot table
        if (sub_argument == 'l') {  latex_tabulate_knot_table(); }
        
        // print latex torus HOMLFY table
        if (sub_argument == 'h') { latex_tabulate_HOMFLY_table(); }
        
        // print latex torus HOMLFY table
        if (sub_argument == 'm') {  latex_tabulate_KBSM_table(); }
    
        // full print knot groups
        if (sub_argument == 'f') { print_knot_groups(); }
        
        make_operation(CLEAN_QUIT);
    }
    
    // random search
    
    if (argument == 'r') {
        
        make_operation(LOAD_GC);
        cout << "Performing random walks..." << endl;
        
        random_walk(get_knot_by_id(167), get_knot_by_id(168));
        
        make_operation(QUIT);
    }
    
    return 0;
}
