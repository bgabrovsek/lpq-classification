//  KBSM.h
//  Created by Boštjan on 4/22/13.
//
//  Calculation of the KBSM of a knot in the solid torus

#ifndef Lpq_KBSM_h
#define Lpq_KBSM_h

#include <list>

#include "knot.h"
#include "common.h"
#include "tree.h"

// Global variables (list of arcs in each circle of the states)
list<u_region> arcs;

// join arcs a and b
void joinArcs(u_region a, u_region b) {
    list<u_region>::iterator it_a = arcs.end(), it_b = arcs.end();
    for (list<u_region>::iterator it=arcs.begin(); it!=arcs.end(); it++) {
        if ((*it) & ((u_region)1<<a)) it_a = it;
        if ((*it) & ((u_region)1<<b)) it_b = it;
    }
    (*it_a) |= (*it_b);
    if (it_a != it_b) arcs.erase(it_b);
}

// get the shortest path of punctured regions from region 0 to region 1
// TODO: optimization
u_region shortestPath(Cknot *K) {
    
    if (K->reg_0 == K->reg_1) throw 20; // affine
    
    list<u_region> white;
    list< pair<u_region,u_region> > black, grey;
    list<u_region>::iterator itw;
    list< pair<u_region,u_region> >::iterator itb;

    for (int i=0; i<K->n_reg; i++) // push available regions to white, reg_0 to grey
        if (i == K->reg_0) grey.push_back(make_pair(K->regions[K->reg_0],0)); else white.push_back(K->regions[i]);

    while (!white.empty()) {
        for (itw = white.begin(); itw != white.end(); itw++) // loop through unused regions
            for (itb = grey.begin(); itb != grey.end(); itb++)
                if (*itw & (*itb).first) { // region between white and grey intersect, add new path to black
                    black.push_front(make_pair(*itw|(*itb).first,(*itb).second|((u_region)1<<firstBitSet((*itw&(*itb).first))))); // add new path
                    if (*itw == K->regions[K->reg_1]) return (*black.begin()).second; // path found
                    itw = white.erase(itw); itw--; // delete from white
                    break;
                }
        grey.clear();
        grey.splice(grey.begin(), black);
    }
    throw 20;
    // return u_region(0xffffffffffffffff);
}

// calculate the winding number using the shortest path from the punctured riogns
// arcs of circle and arcs of shortest path
int windingNumber(u_region circle, u_region path) {
    return count_bits(circle & path) & 1;
}

// writhe of a diagram
int writhe(Cknot *K) { return K->n - 2*count_bits(K->signs);}

// calculate the KBSM of K in the solid torus
Cmultivariate * KBSM(Cknot *K) {

    int i, triv_n, proj_n;
    u_sign state;
    u_region sPath; // shortest path between reg_0 and reg_1
    
    // get the shortest path from the punctured regions
    if (K->reg_0 == K->reg_1) sPath = (u_region)0; else sPath = shortestPath(K);
    
    Cmultivariate * poly = new Cmultivariate(); // the KBSM
    Cmultivariate * poly_triv = new Cmultivariate("-a^2 -a^-2"); // for substituting trivial affine circles
    Cmultivariate * poly_proj = new Cmultivariate("x"); // for substituting trivial non-affine circles
    Cmultivariate * tmp = new Cmultivariate();

    for (state=0; state < ((u_sign)2<<K->n); ++state, ++state)  { // cube of states
        
        arcs.clear();
        for (i=0; i<K->n<<1;i++) arcs.push_back((u_region)1 << i); // push all arcs into a list

        for (i=1; i<=K->n; i++) { // loop through crossings
            if ( BIT(state,i) ^ BIT(K->signs,i)) { // positive smoothening and positive crossings or vice versa
                joinArcs((u_region)K->gw_pos(i),(u_region)K->gw_pos(-i));
                joinArcs((u_region)MOD(K->gw_pos(i)-1,K->n2),(u_region)MOD(K->gw_pos(-i)-1,K->n2));
            } else {
                joinArcs((u_region)K->gw_pos(i),(u_region)MOD(K->gw_pos(-i)-1,K->n2));
                joinArcs((u_region)K->gw_pos(-i),(u_region)MOD(K->gw_pos(i)-1,K->n2));
            }
        }

        triv_n = proj_n = 0; // trivial & projective components
        for (list<u_region>::iterator it=arcs.begin(); it!=arcs.end(); it++)
            if (windingNumber(*it, sPath)) proj_n++; else triv_n++;
        
        // calculate the polynomial of the current state
        tmp->zerofy();
        *tmp += *( *(*poly_triv ^ triv_n) * *(*poly_proj ^ proj_n) );
        tmp->multiplyByVariable('a',K->n - 2*count_bits(state));
        
        *poly += *tmp; // add current state to the KBSM
    }
    
  
    poly->multiplyByVariable((writhe(K)&1) ? -1 : 1, 'a',-3*writhe(K)); // unnormalized
    poly->simplify();

    return poly;

}

#define ABS_SUB(a, x) ((a)>=0 ? ((a)-(x)) : ((a)+(x)))

long int f(long int n) { return (n == 0 ? 1 : n*f(n-1));} // n!

long int choose(long int k, long int n) { if (k>n) throw 20; return f(n) / f(k) / f(n-k); } // choose(k,n)

int writhe(Ctree *t) { // writhe of a arrow diagram knot presented with a tree
    if (!t->has_grandsons()) {
        int w = 0;
        for (t_c_it it = t->child.begin(); it != t->child.end(); it++) w += (*it)->val;
        return w;
    } else {
        int depth = t->depth();
        Ctree * t_ = t->deepest_child();
        int val = t_->val;
        t_->parent->child.remove(t_);
        return (depth-1)*2*abs(val) + writhe(t) + val;
    }
}

// calculates the KBSM of an arrow diagram without crossings (for calculating the KBSM basis relations)
Cmultivariate * KBSM_arrow_(Ctree *t) { 
    
    if (!t->has_grandsons()) { // flat
        Ctree * t_ = t->find_first_son_interval(0, 1, false); // find a son with val != {0,1}
        if (t_ == NULL) { // no such sun
            int c_0 = t->count_vals(0);
            int c_1 = t->count_vals(1);
            return *(Cmultivariate("-a^2 - a^-2")^c_0) * *(Cmultivariate("-a^3x")^c_1); // framed "x", unframed "-a^3"
        } else {
            int v = t_->val;
            t_->val = ABS_SUB(v, 2);
            Ctree * t__ = t->deep_copy();
            t_->val = ABS_SUB(v,1);
            t_->parent->child.push_back(new Ctree(1, t_->parent->parent));
            if (v > 0) return *(Cmultivariate("-a^-2") * *KBSM_arrow_(t)) + *(Cmultivariate("-a^2") * *KBSM_arrow_(t__));
            else return *(Cmultivariate("-a^-4") * *KBSM_arrow_(t)) + *(Cmultivariate("-a^-2") * *KBSM_arrow_(t__));
        }
    } else {
        
        Ctree * t_ = t->deepest_child();
        if (!t_->is_leaf()) throw 20;
        if (t_->val == 0) { // unknot
            t_->parent->child.remove(t_); // remove child
            return (Cmultivariate("-a^2 -a^-2") * *KBSM_arrow_(t));
        
        } else
        if (t_->val == 1) {
            
            t_->parent->child.remove(t_);
            t_->parent->val++;
            Ctree *t__ = t->deep_copy();
            t_->parent->val--;
            t_->create_uncle(1);
            return *(Cmultivariate("1 - a^4") * *KBSM_arrow_(t__)) + *(Cmultivariate("a^-2") * *KBSM_arrow_(t));
            
        } else { // value != {0,1}
            
            int v = t_->val;
            
            t_->val = ABS_SUB(v, 2);
            Ctree * t__ = t->deep_copy();
            t_->val = ABS_SUB(v,1);
            t_->create_sibling(1);
            
            if (v > 0) return *(Cmultivariate("-a^-2") * *KBSM_arrow_(t)) + *(Cmultivariate("-a^2") * *KBSM_arrow_(t__));
            else return *(Cmultivariate("-a^-4") * *KBSM_arrow_(t)) + *(Cmultivariate("-a^-2") * *KBSM_arrow_(t__));
        }
    }
}

// returns the KBSM od a knot arrow diagram (without crossings) presented by a tree structre
Cmultivariate * KBSM_arrow(Ctree *t) {
    Ctree *t_ = t->deep_copy();
    return KBSM_arrow_(t_);
}

// returns the KBSM od a knot arrow diagram (without crossings) presented by a string of the tree structre
Cmultivariate * KBSM_arrow(string s) {
    Ctree *T = new Ctree(s);
    return KBSM_arrow(T);
}

// performs a L(p,1) slide move on a knot arrow diagram represented by a string of the tree structure and returns the KBSM
Cmultivariate * slidep1(string s, int p) {
    string s1;
    int i = 2; // position of the number of arrows of the outer-most circle
    int g = next_i(s,&i); i++;
    s.erase(0,i);
    if (s.length() > 2) // perform the slide
        s1 = "[[" + itoa(p-g) + "[" + s + "]";
    else
        s1 = "[[" + itoa(p-g) + "]]";
    
    return KBSM_arrow(new Ctree(s1)); // poly->multiplyByVariable(-1, 'a',3); (not necessary)
}

// performs a L(p,2) slide move on a knot arrow diagram represented by a string of the tree structure and returns the KBSM
Cmultivariate * slidep2(string s, int p) {
    
    string s1, s2;
    int b = p/2 + 1; // a nd b new arrows (outer-most and second outer-most circle)
    int a = p/2;
    int i = 2;
    int g = next_i(s,&i); i++;
    s.erase(0,i);
    
    if (s.length() > 2) {
        s1 = "[[" + itoa(b) + "[" + itoa(a-g) +"[" + s +"]]"; //?
        s2 = "[[" + itoa(g+1) + "]" + s + "]"; // ?
    } else {
        s1 = "[[" + itoa(b) + "[" + itoa(a-g) + "]]]";
        s2 = "[[" + itoa(g+1) + "]]";
    }

    // calculate the combination of both new arrow diagrams by the skein relation
    Cmultivariate * poly = *(Cmultivariate("a") * *KBSM_arrow(new Ctree(s1))) +
    *(Cmultivariate("a^-1") * *KBSM_arrow(new Ctree(s2)));

    return poly;
}

// performs a L(p,3) slide move on a knot arrow diagram represented by a string of the tree structure and returns the KBSM
Cmultivariate * slidep3(string s, int p) {
    
    string s1, s2, s3, s4;
    
    int c = p/3 + ((p%3)>=1?1:0); // a,b,c - arrows on the outer-most circles
    int b = p/3 + ((p%3)>=2?1:0);
    int a = p/3;
    int i = 2;
    int g = next_i(s,&i); i++;
    s.erase(0,i);
    
    if (s.length() > 2) {
        s1 = "[[" + itoa(c) + "[" + itoa(b) + "[" + itoa(a-g) + "]" + s + "]]";
        s2 = "[[" + itoa(c-b) + "][" + itoa(a-g) + "]" + s + "]";
        s3 = "[[" + itoa(c) + "[" + itoa(b-(a-g)) + "[" + s + "]]";
        s4 = "[[" + itoa(c-b+a-g) + "]" + s + "]";
    } else {
        s1 = "[[" + itoa(c) + "[" + itoa(b) +"[" + itoa(a-g) + "]]]]";
        s2 = "[[" + itoa(c-b)+ "][" + itoa(a-g) + "]]]]";
        s3 = "[[" + itoa(c) + "[" + itoa(b-(a-g)) + "]]]";
        s4 = "[[" + itoa(c-b+a-g) + "]]";
    }
    
    // calculate the KBSM by the skein realation of the arrow diagrams of si
    Cmultivariate * poly =  *(*(*(Cmultivariate("a^2") * *KBSM_arrow(new Ctree(s1)))
    + *KBSM_arrow(new Ctree(s2)))
    + *KBSM_arrow(new Ctree(s3)))
    + *(Cmultivariate("a^-2") * *KBSM_arrow(new Ctree(s4)));

    // poly->multiplyByVariable(-1, 'a',9);
    return poly;
}

// performs a L(p,3) slide move on a knot arrow diagram represented by a string of the tree structure and returns the KBSM
Cmultivariate * slidep5(string s, int p) {
    
    if (p!=12) throw 20;
    
    string s1, s2, s3, s4, s5, s6, s7, s8_, s9;
    int i = 2;
    int g = next_i(s,&i); i++;
    s.erase(0,i);
    
    if (s.length() > 2) {
    } else {
        s1 = "[[3[3[2[2[" + itoa(2-g) + "]]]]]]"; //a^4
        s2 = "[[3[3[2[" + itoa(g) + "]]]]]"; // a^2
        s3 = "[[3[3[" + itoa(2-g) + "]]]]"; //-a^4+1
        s4 = "[[3[1][2[" + itoa(2-g) +"]]]]"; //a^2
    
        s5 = "[[2[2[" + itoa(2-g)+ "]]]"; //-a^4
        s6 = "[[" + itoa(2-g) + "]]"; //a^4
        s7 = "[[2[" + itoa(g)+ "]]]"; //-a^2
        s8_= "[[3[1][" + itoa(g)+"]]]"; //1
        s9= "[[3[" + itoa(1+g)+"]]]"; // a^-2
    }

    Cmultivariate *poly;
    
    // calculate the KBSM by the skein realation of the arrow diagrams of si
    poly = Cmultivariate("a^4") * *KBSM_arrow(new Ctree(s1));
    *poly += *(Cmultivariate("a^2") * *KBSM_arrow(new Ctree(s2)));
    *poly += *(Cmultivariate("-a^4 + 1") * *KBSM_arrow(new Ctree(s3)));
    *poly += *(Cmultivariate("a^2") * *KBSM_arrow(new Ctree(s4)));
    *poly += *(Cmultivariate("-a^4") * *KBSM_arrow(new Ctree(s5)));
    *poly += *(Cmultivariate("a^4") * *KBSM_arrow(new Ctree(s6)));
    *poly += *(Cmultivariate("-a^2") * *KBSM_arrow(new Ctree(s7)));
    *poly += *KBSM_arrow(new Ctree(s8_));
    *poly += *(Cmultivariate("a^-2") * *KBSM_arrow(new Ctree(s9)));

    return poly;
}

// calculate the framing requred to normalize the KBSM with the minimal term of power m
int required_framing(int m) {
    int i = 0;
    if (m > 0) {
        do {
            m-=3;
            i--;
        } while (m >= 0);
            return i+1;
    } else
    if (m < 0)
    {
        do {
            m += 3;
            i++;
        } while (m < 0);
        return i;
    } else return 0;
}

// normalize the KBSM
void normalize_framing(Cmultivariate * poly) {
    int fr = required_framing( poly->min_power('a'));
    poly->multiplyByVariable((fr&1) ? -1 : 1, 'a',3*fr);
}

// perform a slide move in the lens space l to a arrow diagram repreesnted by the tree stored in string s
Cmultivariate * slide(int l, string s) {
    switch (lens_homeo[l][1]) {
        case 1:
            return slidep1(s,lens_homeo[l][0]);
        case 2:
            return slidep2(s,lens_homeo[l][0]);
        case 3:
            return slidep3(s,lens_homeo[l][0]);
        case 5:
            return slidep5(s,lens_homeo[l][0]);
        default: throw 20;
    }
}

// output the L(p,q) relations of the KBSM generators of the solid torus, output Mathematica format
void print_KBSM_relations() {

    for (int l= 0; l<N_HOMEO-2; l++) { // loop through all lens apces
        
        cout << "(* " <<lpq_s[l] << " *)" << endl;

        cout << "A" << lens_homeo[l][0] << lens_homeo[l][1] << " = " << "Solve[{";
        switch (lens_homeo[l][1]) {
            case 1: case 2: case 3: case 5:
                for (int i=0;i<(13-lens_homeo[l][0])/2;i++) {
                    cout << KBSM_arrow(new Ctree("[[" + itoa((lens_homeo[l][0]-1)/2 - i )+"]]")) << " == ";
                    cout << slide(l,"[[" + itoa( (lens_homeo[l][0]-1)/2 - i ) + "]]");
                    if (i+1<(13-lens_homeo[l][0])/2) cout << ", ";
                }
                
                if (l == 16) cout << ", "<<KBSM_arrow(new Ctree("[[3]]")) << " == " << slide(l,"[[3]]"); // extra relation
                break;
                
            default: throw 20;
        }
        
        // output variables we wish to solve
        cout << "}, {";
        if (lens_homeo[l][0] <= 3) cout << "x^2, x^3, x^4, x^5, x^6"; else
        if (lens_homeo[l][0] <= 5) cout << "x^3, x^4, x^5, x^6"; else
        if (lens_homeo[l][0] <= 7) cout << "x^4, x^5, x^6"; else
        if (lens_homeo[l][0] <= 9) cout << "x^5, x^6"; else
        if (lens_homeo[l][0] <= 10) cout << "x^6"; else
        if ((lens_homeo[l][0] <= 11) && (lens_homeo[l][1] <= 2)) cout << "x^6"; else
        cout << "x^6, x^8";
        cout << "}]";
        
        cout << endl << endl;
        
    }
}

/*
 
 Print["generator_KBSM[stol(\"L(2,1)\")][2] = mv(\"",InputForm[Expand[ReplaceAll[x2, A21[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(2,1)\")][3] = mv(\"",InputForm[Expand[ReplaceAll[x3, A21[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(2,1)\")][4] = mv(\"",InputForm[Expand[ReplaceAll[x4, A21[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(2,1)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A21[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(2,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A21[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(3,1)\")][2] = mv(\"",InputForm[Expand[ReplaceAll[x2, A31[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(3,1)\")][3] = mv(\"",InputForm[Expand[ReplaceAll[x3, A31[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(3,1)\")][4] = mv(\"",InputForm[Expand[ReplaceAll[x4, A31[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(3,1)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A31[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(3,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A31[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(4,1)\")][3] = mv(\"",InputForm[Expand[ReplaceAll[x3, A41[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(4,1)\")][4] = mv(\"",InputForm[Expand[ReplaceAll[x4, A41[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(4,1)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A41[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(4,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A41[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(5,1)\")][3] = mv(\"",InputForm[Expand[ReplaceAll[x3, A51[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(5,1)\")][4] = mv(\"",InputForm[Expand[ReplaceAll[x4, A51[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(5,1)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A51[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(5,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A51[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(5,2)\")][3] = mv(\"",InputForm[Expand[ReplaceAll[x3, A52[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(5,2)\")][4] = mv(\"",InputForm[Expand[ReplaceAll[x4, A52[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(5,2)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A52[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(5,2)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A52[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(6,1)\")][4] = mv(\"",InputForm[Expand[ReplaceAll[x4, A61[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(6,1)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A61[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(6,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A61[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(7,1)\")][4] = mv(\"",InputForm[Expand[ReplaceAll[x4, A71[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(7,1)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A71[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(7,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A71[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(7,2)\")][4] = mv(\"",InputForm[Expand[ReplaceAll[x4, A72[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(7,2)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A72[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(7,2)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A72[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(8,1)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A81[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(8,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A81[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(8,3)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A83[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(8,3)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A83[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(9,1)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A91[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(9,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A91[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(9,2)\")][5] = mv(\"",InputForm[Expand[ReplaceAll[x5, A92[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(9,2)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A92[[1]]]]], "\");"];
 
 Print["generator_KBSM[stol(\"L(10,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A101[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(10,3)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A103[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(11,1)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A111[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(11,2)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A112[[1]]]]], "\");"];
 Print["generator_KBSM[stol(\"L(11,3)\")][6] = mv(\"",InputForm[Expand[ReplaceAll[x6, A113[[1]]]]], "\");"];
 
 */

/* 
 
 KBSM:
 
 L(2"1)
 
 x^2 "1 + a^2 + a^4 + a^6"
 x^3 "x*a^-2 + a^2x + 2a^6x"
 x^4 "1 + a^-4 + a^-2 + a^2 + a^4 + 3a^6 + 3a^8 + 3a^10 + 2a^12"
 x^5 "x + x*a^-8 + x*a^-4 + 4a^4x + 4a^8x + 5a^12x"
 x^6 "1 + a^-12 + a^-10 + a^-8 + a^-6 + a^-4 + a^-2 + 5a^2 + 5a^4 + 5a^6 + 5a^8 + 5a^10 + 9a^12 + 9a^14 + 9a^16 + 5a^18"
 
 L(3"1)
 
 x^2 "a^2 + a^6 - a^2x"
 x^3 "-a^2 - a^6 + a^2x + 2a^6x"
 x^4 "a^4 + 3a^8 + 2a^12 - x - a^4x - 3a^8x"
 x^5 "-1 - a^4 - 4a^8 - 4a^12 + x + a^4x + 4a^8x + 5a^12x"
 x^6 "a^-2 + a^2 + a^6 + 5a^10 + 9a^14 + 5a^18 - x*a^-2 - a^2x - 5a^6x - 5a^10x - 9a^14x"
 
 L(4"1)
 
 x^3 "a^2x + a^4x + 2a^6x"
 x^4 "a^4 - a^12 + a^2x^2 + 3a^6x^2"
 x^5 "a^2x + a^4x + a^6x + 4a^8x + 4a^10x + 5a^12x"
 x^6 "4 a^10 - 4a^18 + x^2 + a^4x^2 + 5a^8x^2 + 9a^12x^2"
 
 L(5"1)
 
 x^3 "a^4 + a^8 + a^2x + 2a^6x - a^2x^2"
 x^4 "-a^8 - a^12 - a^6x + a^2x^2 + 3a^6x^2"
 x^5 "4 a^10 + 4a^14 + a^4x + 4a^8x + 5a^12x - a^4x^2 - 4a^8x^2"
 x^6 "-a^10 - 5a^14 - 4a^18 - a^4x - a^8x - 5a^12x + a^4x^2 + 5a^8x^2 + 9a^12x^2"
 
 L(5"2)
 
 x^3 "-a - a^5 + a ax + 2a^6x + x^2*a^-1"
 x^4 "-1 - a^12 + ax + x^2*a^-2 + 3a^6x^2"
 x^5 "-4 a^7 - 4a^11 + ax*a^-1 + 4a^7 ax + 5a^12x + x^2*a^-3 + 4a^5x^2"
 x^6 "-a^2 - 5a^6 - 4a^18 + x*a^-5 + x*a^-1 + 5a^7x + x^2*a^-4 + 5a^4x^2 + 9a^12x^2"
 
 L(6"1)
 
 x^4 "-a^6 - a^8 - a^10 - a^12 + a^2x^2 + a^4x^2 + 3a^6x^2"
 x^5 "-a^8x - 3a^12x + a^2x^3 + 4a^6x^3"
 x^6 "-a^10 - 5a^12 - 5a^14 - 5a^16 - 4a^18 + a^4x^2 + a^6x^2 + 5a^8x^2 + 5a^10x^2 + 9a^12x^2"
 
 L(7"1)
 
 x^4 "-a^8 - a^12 + a^4x + 2a^8x + a^2x^2 + 3a^6x^2 - a^2x^3"
 x^5 "a^8 + a^12 - 2a^8x - 3a^12x - a^6x^2 + a^2x^3 + 4a^6x^3"
 x^6 "-a^10 - 5a^14 - 4a^18 + a^6x + 6a^10x + 10a^14x + a^4x^2 + 5a^8x^2 + 9a^12x^2 - a^4x^3 - 5a^8x^3"
 
 
 L(7"2)
 
 x^4 "-a^8 - a^12 - ax - 2a^5x + a^2x^2 + 3a^6x^2 + x^3*a^-1"
 x^5 "-a^3 - a^7 - a^3 ax - x - 3a^12x + ax^2 + x^3*a^-2 + 4a^6x^3"
 x^6 "-a^2 - 5a^14 - 4a^18 - x*a^-1 - a^3x - 5a^7x - 10a^11x + x^2 + 5a^8x^2 + 9a^12x^2 + x^3*a^-3 + 5a^5x^3"
 
 L(8"1)
 
 x^5 "-a^6x - 2a^8x - 2a^10x - 3a^12x + a^2x^3 + a^4x^3 + 4a^6x^3"
 x^6 "-a^10 + a^18 - 2a^8x^2 - 6a^12x^2 + a^2x^4 + 5a^6x^4"
 
 L(8"3)
 
 x^5 "-a^2x - 2a^4x - 2a^6x - 3a^12x + x^3 + x^3*a^-2 + 4a^6x^3"
 x^6 "-2 a^2 + a^10 + a^18 + x^2*a^-4 - 3a^4x^2 - 6a^12x^2 + x^4*a^-2 + 5a^6x^4"
 
 L(9"1)
 
 x^5 "-a^10 - a^14 - 2a^8x - 3a^12x + a^4x^2 + 3a^8x^2 + a^2x^3 + 4a^6x^3 - a^2x^4"
 x^6 "a^14 + a^18 + a^8x + 2a^12x - 3a^8x^2 - 6a^12x^2 - a^6x^3 + a^2x^4 + 5a^6x^4"
 
 L(9"2)
 
 x^5 "a^7 + a^11 - 2a^8x - 3a^12x - ax^2 - 3a^5x^2 + a^2x^3 +  4 a^6x^3 + x^4*a^-1"
 x^6 "a^6 + a^18 - a^3x - 2a^7x - x^2 - 2a^4x^2 - 6a^12x^2 + ax^3 + x^4*a^-2 + 5a^6x^4"
 
 L(10"1)
 
 x^6 "a^12 + a^14 + a^16 + a^18 - a^6x^2 - 3a^8x^2 - 3a^10x^2 - 6a^12x^2 + a^2x^4 + a^4x^4 + 5a^6x^4"
 
 L(10"3)
 
 x^6 "-a^2 + a^6 + a^8 + a^10 + a^12 + a^18 - a^2x^2 - 3a^4x^2 - 3a^6x^2 - 6a^12x^2 + x^4 + x^4*a^-2 + 5a^6x^4"
 
 L(11"1)
 
 x^6 "a^14 + a^18 - 2a^10x - 3a^14x - 3a^8x^2 - 6a^12x^2 + a^4x^3 + 4a^8x^3 + a^2x^4 + 5a^6x^4 - a^2x^5"
 
 L(11"2)
 
 x^6 "a^14 + a^18 + 2a^7x + 3a^11x - 3a^8x^2 - 6a^12x^2 - ax^3 - 4a^5x^3 + a^2x^4 + 5a^6x^4 + x^5*a^-1"
 
 L(11"3)
 
 x^6 "-a^2 + a^6 + a^10 + a^18 - x - a^4x - 3a^12x - 3a^4x^2 - 6a^12x^2 + x^3*a^-2 + 4a^6x^3 + x^4*a^-2 + 5a^6x^4 - x^5"

 
 */


#endif
