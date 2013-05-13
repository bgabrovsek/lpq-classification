
// randomized.h
// make randomized Reidemeister move walks, simplify the knot and compare to another knot

#ifndef Lpq_randomized_h
#define Lpq_randomized_h

#include "knot.h"
#include "reidemeister_moves.h"
#include <stdlib.h> 
#include <time.h>

#define N_TRIES 100000 // number of tries
#define WALK_LENGTH 30000 // length of the walk
#define MAX_CROSSINGS_WALK 20 // number of allowed crossings in a walk

// make a random Reidemeister move b
void random_move(Cknot **K, u16 b) {
    
    t_reidemeister_list rList; // list of Reidemeister moves
    t_reidemeister_list::iterator ir;
    
    findReidemeisterMoveSites(*K, b , &rList);
    
    if (rList.size() <= 1) return;
    
    int r = rand() % (rList.size());
    ir = rList.begin();
    while (r-- > 0) ir++; // get i-th Reidemeister move
    
    rMove(*K,&(*ir));
}


// make a random walk on knot K and print the result

void random_walk(Cknot *K, Cknot *K_desired) {

    srand((int)time(NULL));
    
    for (int i = 0; i < N_TRIES; i++) {
        
        Cknot *Q = K->deepCopy();
        
        // perform the walk
        for (int j = 0; j < WALK_LENGTH; j++) {
            
            if (Q->n <= MAX_CROSSINGS_WALK-2) random_move(&Q, MODIFY_FLYPE | MODIFY_R_III | REMOVE_R_I | REMOVE_R_II | CREATE_R_II | CREATE_R_I); else
            if (Q->n <= MAX_CROSSINGS_WALK-1) random_move(&Q, MODIFY_FLYPE | MODIFY_R_III | REMOVE_R_I | REMOVE_R_II | CREATE_R_I); else
            random_move(&Q, MODIFY_FLYPE | MODIFY_R_III | REMOVE_R_I | REMOVE_R_II);
        }
        
        cout << "(" << i <<") ";
        
        // reduce the number of crossings to those of the original knot K
    
        while (Q->n > K->n) {
            
           BST_shrink_up_to_n_crossings(&Q, UNORIENTED, 30, 10, K->n);
    
            if (Q->n > K->n) {
                BST_reducible_up_to_n_crossings(Q, UNORIENTED, 30, 0x1001001001, 2, K->n);
                Q = BST_minimal_knot;
            }
        }
        
        sortRegionsDots(Q); // canonical regions
        
        cout << Q << endl;
        
        // does the new knot equal the desired one?
        if (equal_knots(Q,K_desired)) cout << "Equivalence found" << endl;
        
    }
}

// candidates: Knot 231: 1 -1 -2 3 -4 5 -3 6 7 2 -6 4 -5 -7 +-+++-+ *0 269 258C 36A 018D 47AC 4B &35B 179D
//Knot 231: 1 2 3 -4 5 -1 6 7 -2 -3 4 -5 -7 -6 +-++++- &6C 0248A 3A 057C 1379B 46BD *18 5D 29
//Knot 231: 1 2 3 -4 5 -1 6 7 -2 -3 4 -5 -7 -6 +-++++- &6C 0248A 3A 057C 1379B 46BD *18 5D 29
//Knot 231: 1 2 3 -4 5 -1 6 7 -2 -3 4 -5 -7 -6 +-++++- &6C 0248A 3A 057C 1379B 46BD *18 5D 29


#endif
