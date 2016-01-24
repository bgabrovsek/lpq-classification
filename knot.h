//  knots.h
//  Created by Boštjan on 4/22/13.
//
//  Knot class and knot manipulation functions
//
//  KNOT DATA STRUCTURE:
//  Gauss word: 1  -1   2  -3   4  -2   5  -4   3  -5
//  Arcs:         0   1   2   3   4   5   6   7   8   9
//  Region 0: binary "0000000001" (= arc 0)
//  Region 1: binary "0100100100" (= arc 2, 5, 8)
//  Signs ++--+ = binary "01100" (+ = 0, - = 1, least significant bit = crossing #1, most significant bit = crossing #n)

#ifndef Lpq_knot_h
#define Lpq_knot_h

#include <list>
#include <vector>
#include <string>
#include <set>
#include <iostream>

#include "common.h"
#include "numbers.h"

class Cknot;
ostream &operator<<(ostream &out, Cknot * K);

class Cknot {
public:
    u_number n, n2; // number of crossings, twice
    s_letter gw_[MAX<<1]; // Gauss word e.g. 1,-1,2,-2
    u_sign signs; // in binary, starting with 1st bit (not 0th)
    u_region regions[MAX+2]; // region in 64 bit binary (arcs in binary starting from 0th bit)
    u_number n_reg , reg_0, reg_1; // number of regions, index of 0 region, index of 1 region
    s16 ID; // knot ID
    
    // position of letter l in the GW
    u_number gw_pos(s_letter l) { for (int i=0; i<n<<1;i++) if (gw_[i] == l) return i;
        cout << this << endl << "ERROR finding " << (int)l << " n= " << (int)n << " n<<1= "<< (int)(n<<1)<<endl; throw 20;}
    
    s_letter gw(int i) { return gw_[MOD(i,n2)]; } // gets the gauss letter at position i (i from Z)
    
    void zerofyKnot(); // set all data to 0

    bool canonical(bool oriented); // put Gauss word into its canonical minimal form, returns true if canonical, false if it has been changed
    bool canonical_meridional(bool oriented); // put Gauss word into its canonical minimal form, also including a meridional rotation
    
    bool realizable(); // must first call generate_regions
    void generateRegions(); // generate regions from scrach
    void check_ok(bool also_regions = true); // checks if letters, signs and regions present a knot (used for debugging)
    string check_ok_reason(); // does a check and returns a string of what went wrong
    
    u_number turn(u_number arc, s_small *orient, s_small direction);
    u_region gen_reg(u_number arc, s_small side); // generates region next to arc
    u_number turn2(u_number arc, s_small *orient, s_small direction);
    u_region gen_reg2(u_number arc, s_small side); // generates region next to arc
    u_number regionGet(u_number arc, s_small side); // get the region that is right or left of arc
    s_small regionSide(u_number r, u_number arc) { return (r==regionGet(arc,RIGHT))?RIGHT:LEFT; } // 0 if r right of arc, else 1
    u_number adjacentRegion(u_number reg, u_number arc); // get adjacent region to reg in terms of arc
    
    Cknot * deepCopy(); // deep copy whole knot, except regions
    
    void copyKnotData(Cknot *K); // deep copy whole knot, except regions
    
    void mirror() { for(u_number i=0;i<n2;i++) gw_[i] *= -1; signs ^= fbf_sign(n+1) ^ (u_sign)1;}
    void flip() {u_number tmp = reg_0; reg_0 = reg_1; reg_1 = tmp;}
    void reverse() {if (n==0){swap(reg_0,reg_1);return;} u_number i; for (i=0;i<n;i++) swap(&(gw_[i]),&(gw_[n2-i-1])); for (i=0;i<n_reg;i++) regions[i]=reverseBits(regions[i],(u_number)(2*n-1)) | ( regions[i] & ((u_region)1 << (2*n-1))  );}
    
    void meridional_rotation() { for(u_number i=0;i<n2;i++) gw_[i] *= -1; u_number tmp = reg_0; reg_0 = reg_1; reg_1 = tmp; }
        
    int bin2reg(u_region b) { for (u_number i=0;i<n_reg;i++) if (regions[i] == b) return i; return -1;}

    // constructors
    Cknot() { zerofyKnot(); }
    Cknot(string s); // "Knot 139: 1 -1 2 -3 4 -5 5 -4 3 -2 +-++- *1 12A 239A 345789 567 *6 48"
    
    //destructor
    ~Cknot() { /*cout << "deleting knot " << ID << endl; */ }
};

// A = B ?
bool equal_knots(Cknot *A, Cknot *B) {
    
    if ((A->n != B->n) || // number of crossings
        (A->n_reg != B->n_reg) || // number of regions
        (((A->signs)&fbf_sign(A->n+1)^(u_sign)1) != ((B->signs)&fbf_sign(B->n+1)^(u_sign)1)) || // signs
        (A->regions[A->reg_0] != B->regions[B->reg_0]) || // dotted regions
        (A->regions[A->reg_1] != B->regions[B->reg_1]))
        return false;
    
    for (u_number i = 0; i<A->n2; i++) // gauss word
        if (A->gw_[i] != B->gw_[i]) return false;
    
    return true;
}

// return region adjacent to arc and reg
u_number Cknot::adjacentRegion(u_number reg, u_number arc) {
    for (u_number r = 0; r<n_reg; r++)
        if ((r != reg) && (((regions[r]) >> (arc)) & 1)) return r;
    throw 20;
}

int winding_number(Cknot *K);

// prints the knot
ostream &operator<<(ostream &out, Cknot * K) {// print knot
    
    int i;
    
    use_extended_bitset = (K->n2 > 20); // if a large number of arcs, print the regions in extended form
    // print knot ID
    out << "Knot " << K->ID << ": ";
    // print GW
    for (i=0;i<2*K->n;i++) out << (int)K->gw_[i] << ' ';
    // print signs
    for (i=1;i<=K->n;i++) out << sign2char((int)BIT(K->signs,i));
    out << " ";
    // print regions
    for (i=0;i<K->n_reg;i++)
        if ((!use_extended_bitset | 1) | (i==K->reg_0) | (i==K->reg_1))
        out << ((i==K->reg_0)?"*":"") << ((i==K->reg_1)?"&":"") <<
        (use_extended_bitset ? bitset_extended(K->regions[i]) : bit_set(K->regions[i])) << (i!=K->n_reg-1?" ":"");

    return out;
}

// print knot
ostream &operator<<(ostream &out, Cknot K) {// print knot
    out << &K;
    return out;
}

// true if A < B, false otherwise (in the lexicographical ordering)
bool smaller_knot(Cknot *A, Cknot *B) { 
    s_small compare;
    u_number i;
    
    // compare number of crossings
    if (A->n != B->n) return (A->n < B->n); 
    
    // compare Gauss word
    for (i=0; i < A->n2; i++) {
        compare = COMPARE_LETTER(A->gw_[i], B->gw_[i]);
        if (compare) return (compare < 0);
    }
    
    // compare signs
    for (i=1; i<= A->n; i++) {
        compare = COMPARE_SIGN(BIT(A->signs,i), BIT(B->signs,i));
        if (compare) return (compare < 0);
    }
    
    // compare punctured regions
    compare = COMPARE_REGION(A->regions[A->reg_0],B->regions[B->reg_0]);
    if (compare) return (compare < 0);
    compare = COMPARE_REGION(A->regions[A->reg_1],B->regions[B->reg_1]);
    if (compare) return (compare < 0);
    
    return false; // knots are equal
}

// obsolete
bool smaller_knot_up_to_flip_mirror(Cknot *A, Cknot *B) {
    s_small compare;
    u_number i;
    
    u_region r0, r1, s0, s1;

    if (A->n != B->n) { return (A->n < B->n); }
    else {
        for (i=0; i < A->n2; i++) {
            compare = COMPARE_LETTER(A->gw_[i], B->gw_[i]);
            if (compare) return (compare < 0);
        }
        
        for (i=1; i<= A->n; i++) {
            compare = COMPARE_SIGN(BIT(A->signs,i), BIT(B->signs,i));
            if (compare) return (compare < 0);
        }
        
        r0 = MINIMUM(A->regions[A->reg_0],A->regions[A->reg_1]);
        r1 = MAXIMUM(A->regions[A->reg_0],A->regions[A->reg_1]);
        s0 = MINIMUM(B->regions[B->reg_0],B->regions[B->reg_1]);
        s1 = MAXIMUM(B->regions[B->reg_0],B->regions[B->reg_1]);
        
        compare = COMPARE_REGION(r0,s0);
        if (compare) return (compare < 0);
        compare = COMPARE_REGION(r1,s1);
        if (compare) return (compare < 0);
        
    }
    return false; // they are equal
    
}

// compare structure for <map>
struct compareKnots
{
    bool operator()(Cknot * K1, Cknot * K2) const // K1 < K2 (true)
    {
        return smaller_knot(K1, K2); // (*K1) < (*K2);
    }
};


// make a deep copy of the knot
Cknot * Cknot::deepCopy() {
    Cknot *L = new Cknot();
    L->n = n; L->n2 = n2;
    memcpy(L->gw_,gw_,sizeof(gw_));
    L->signs = signs;
    memcpy(L->regions,regions,sizeof(regions));
    L->n_reg = n_reg;
    L->reg_0 = reg_0; L->reg_1 = reg_1;
    L->ID = ID;
    return L;
}


// copy the knot into K
void Cknot::copyKnotData(Cknot *K) {
    n = K->n; n2 = K->n2;
    memcpy(gw_,K->gw_,sizeof(gw_));
    signs = K->signs;
    memcpy(regions,K->regions,sizeof(regions));
    n_reg = K->n_reg;
    reg_0 = K->reg_0; reg_1 = K->reg_1;
    ID = K->ID;
}

// start traveling through used arcs and fill the region
void fillReg(Cknot *K, u_number arc, s_small orient, s_small rl, u_region *usedArcs, u_region *regRight, u_region *regLeft) {
    
    *regLeft = *regRight = (u_region)0;
    u_number arc_, n_arc_;
    s_small o_,  n_o_;
    o_ = orient;
    arc_ = arc;
    
    do {
        ((rl^o_) ? *regLeft : *regRight) |= (u_region)1 << arc_;
        n_o_ = o_;
        n_arc_ = K->turn(arc_, &n_o_, rl);
        
        if (BIT(usedArcs[rl ^ n_o_],n_arc_)) {
            o_ = n_o_;
            arc_ = n_arc_;
        } else {
            arc_ = MOD(arc_+(o_?-1:1),K->n2);
            if (!BIT(usedArcs[rl^n_o_],arc_)) {
                o_ = !o_;
                arc_ = MOD(arc_+(o_?-1:1),K->n2);
            }
        }
        
    } while ((arc != arc_) || (orient != o_));
}

// check if a Gauss code is realizable
bool Cknot::realizable() {
    u_number i, j;
    u_region usedArcs[2]; // used right/left crossings
    u_number left, right;
    s_small o_l, o_r;
    u_number c_reg, n_reg;
    u_number crl, crr;
    s_small crlo, crro;
    
    u_region regionsLeft[MAX+2],  regionsRight[MAX+2];
    
    u_region regR, regL;
    
    for (i=0;i<MAX+1;i++) regionsLeft[i] = 0;
    for (i=0;i<MAX+1;i++) regionsRight[i] = 0;
    
    usedArcs[0] = usedArcs[1] = (u_region)0; // 0 ... right, 1 ... left;

    c_reg = 0;
    n_reg = 1;
    
    for (i=0;i<n2;i++) { // go through Gauss word
        
        usedArcs[0] |= (u_region)1 << i;
        usedArcs[1] |= (u_region)1 << i;
        
        regionsRight[c_reg] |= (u_region)1 << i;
        regionsLeft[c_reg] |= (u_region)1 << i;
        
        o_l = o_r = 0;
        left = turn(i,&o_l,LEFT);
        right = turn(i,&o_r,RIGHT);
  
        
        if (( BIT(usedArcs[o_r],right) ) && ( !BIT( (o_r ? regionsLeft[c_reg] :  regionsRight[c_reg]), right ) ))
        {/*cout << "NON RALIZABLE RIGHT." << endl; */return false;}

        if (( BIT(usedArcs[!o_l],left) ) && ( !BIT( (o_l ? regionsRight[c_reg] :  regionsLeft[c_reg]), left ) ))
        {/*cout << "NON RALIZABLE LEFT." << endl;*/ return false;}

        regR = regL = (u_region)0;
        if ( BIT(usedArcs[o_r],right) ) { // right hit
            fillReg(this,right,o_r,0,usedArcs,&regR,&regL);
        } else
        if ( BIT(usedArcs[!o_l],left) ) { // right hit
            fillReg(this,left,o_l,1,usedArcs,&regR,&regL);
        }
        
        if (regR | regL) {
            regionsRight[c_reg] &= ~regR;
            regionsLeft[c_reg] &= ~regL;
            regionsRight[n_reg] = regR;
            regionsLeft[n_reg] = regL;
            n_reg++;
            
            // change the current region
            crlo = crro = 1;
            crr = turn(MOD(i+1,n2), &crro, 0); // one enough
            crl = turn(MOD(i+1,n2), &crlo, 1);
            
            for (j=0;j<n_reg;j++) {
                if (BIT((crro ? regionsLeft[j] : regionsRight[j]),crl)) { c_reg = j; /*cout << "current " << (int)j << "    REG: "; */break; }
                if (BIT((crro ? regionsLeft[j] : regionsRight[j]),crr)) { c_reg = j; /*cout << "current " << (int)j << "    REG: "; */break; }
            }

        }
    }
    
    return true; // process has ended, no errors found
}

// compares the knot before and after its meridional rotation, returns the smallest
bool Cknot::canonical_meridional(bool oriented) {
    Cknot *K = deepCopy();
    bool b, b_;
    b = canonical(oriented);
    K->meridional_rotation();
    b_ = K->canonical(oriented);
    
    if (smaller_knot(K, this)) {
        copyKnotData(K);
        delete K;
        return true;
    }
    delete K;
    return b;
}

// global macros
#define SHIFT_REGIONS if (reverse==1) for (r=0;r<n_reg;r++) regions[r] = circularShiftRight(regions[r],shift,n2); else for (r=0;r<n_reg;r++) regions[r] = circularShiftRight(reverseBits(regions[r],n2),-shift,n2)
#define SHIFT_SIGNS permuteBits(&signs, renum, n)
#define ADJUST_VARS shift = 0; reverse = 1; changed = true
#define RENUM
        

// puts the knot into its canonical form (up to orientation of oriented = true)
bool Cknot::canonical(bool oriented) {

    u_number rnu, renum[n+1], reverse_renum[n+1];
    s_letter tmp, i, shift, new_x, compare, gw_tmp[2*MAX], r, reverse_counter, counter, reverse, sgn_comp;
    u_region reg_0t, reg_1t;
    u_sign r_signs;

    bool changed = false;

    for (reverse = reverse_counter = 1; reverse_counter >= -1; reverse_counter -= 2, reverse -= 2) { // reverse is 1 or -1
        for (counter = shift = 0; counter<n2; shift+=reverse, counter++) { // loop through all shifts

            r_signs = 0; rnu = 0; // clear renumeration table
            
            for (i=0;i<n2;i++) { // loop through all crossings
                
                tmp = ABS(gw(i*reverse+shift)); // crossing letter
                
                if (!(BIT(r_signs, tmp))) { // new letter not yet been renumerated
                    renum[tmp] = ++rnu;
                    reverse_renum[renum[tmp]] = tmp;
                    r_signs |= (u_sign)1 << tmp;
                }
                
                new_x = (s_letter)renum[tmp] * SIGN(gw(i*reverse+shift));
                
                if ((compare = COMPARE_LETTER(new_x,gw(i)))) break;
            }
            
            if (compare > 0) continue; // new gauss word bigger, continue
            if (compare < 0) { // gauss word smaller, we have a new canonical form
            
                memcpy(gw_tmp,gw_,sizeof(gw_tmp));
                
                for (;i<n2;i++) {
                    tmp = ABS(gw_tmp[MOD(i*reverse+shift,n2)]);
                    
                    if (!(BIT(r_signs, tmp))) { // new letter not yet been renumerated
                        renum[tmp] = ++rnu;
                        reverse_renum[renum[tmp]] = tmp;
                        r_signs |= (u_sign)1 << tmp;
                    }
                    
                    gw_[i] = (s_letter)renum[tmp] * SIGN(gw_tmp[MOD(i*reverse+shift,n2)]);
                    
                    // gw_[i] = (s_letter)(renum[tmp] ? renum[tmp] : (renum[tmp] = ++rnu)) * SIGN(gw_tmp[MOD(i*reverse+shift,n2)]);
                }
                
                SHIFT_REGIONS;
                SHIFT_SIGNS;
                ADJUST_VARS;
                continue;

            }
            
            // CHECK SIGNS, since gauss word equal
            
            for (sgn_comp = 0, i=1;i<=n;i++) if ( (sgn_comp = COMPARE_SIGN(BIT(signs,i), BIT(signs,reverse_renum[i]))) ) break;
            
            if (sgn_comp < 0) continue; // new signs bigger, continue
            
            if (sgn_comp > 0) { // new shift makes signs smaller
                SHIFT_REGIONS;
                SHIFT_SIGNS;
                ADJUST_VARS;
                continue;
            }
            
            // CHECK REGIONS; since signs and gauss word equal
            
            reg_0t = ((reverse==1) ? circularShiftRight(regions[reg_0],shift,n2) : circularShiftRight(reverseBits(regions[reg_0],n2),-shift,n2));
            reg_1t = ((reverse==1) ? circularShiftRight(regions[reg_1],shift,n2) : circularShiftRight(reverseBits(regions[reg_1],n2),-shift,n2));
            
            if ((COMPARE_REGION(reg_0t,regions[reg_0])<0) ||
                ((COMPARE_REGION(reg_1t,regions[reg_1])<0) && (COMPARE_REGION(reg_0t,regions[reg_0])==0))) {
                SHIFT_REGIONS;
                SHIFT_SIGNS;
                ADJUST_VARS;
                continue;
            }

        } // shift
        if (oriented) break; // if considering orientations, do not continue with reverse orientation
    } // reverse

    return (!changed);
}

void Cknot::zerofyKnot() {
    n = n2 = 0;
    memset(gw_,0,sizeof(gw_));
    signs = 0;
    memset(regions,0,sizeof(regions));
    n_reg = 0;
    reg_0 = reg_1 = -1;
    ID = -1;
}

// string to knot
Cknot::Cknot(string s) {
    int i, p;
    bool long_form;
    s = s + string(" ");
    
    long_form = (s.find("(") !=  string::npos);
    
    zerofyKnot();
    p = (int)s.find(" +");
    if (p==string::npos) p = (int)s.find(" --");
    if (p==string::npos) p = (int)s.find(" -+"); 
    if (p==string::npos) p = (int)s.find(" - "); // knot length
    
    n = (int)s.find(" ",p+1) - p - 1;
    
    n2 = n << 1;
    p = (s.find("Knot ")==string::npos ? 0 : 4);
    if (p) { ID = next_int(s,&p); p++; } // knot ID
    for (i=0;i<n2;i++) { int w; gw_[i] = w = next_int(s,&p); } // gauss word

    while (s[p] == ' ') p++;
    for (i=1;i<=n;i++) signs |= ((u_sign)char2sign(s[p++]) << i); // signs
    while (s[p] == ' ') p++;
    
    if (p >= s.length()) return; // no regions
    
    if (!long_form) {
        while (p <= s.length()-1) { //read all regions
            if (s[p] == '*') {reg_0 = n_reg; p++;}
            if (s[p] == '&') {reg_1 = n_reg; p++;}
            regions[n_reg] = (u_region)0;
            while ((s[p] != ' ') && (p < s.length())) { regions[n_reg] |= ((u_region)1<<(hex2int(s[p]))); /*cout << hex2int(s[p]) << " "; */ p++;}
            n_reg++;
            if (!regions[n_reg-1]) n_reg--;
            p++;
        }
    } else {
        
        cout << s << endl << endl;
        int r;
        string ss;
        while (p <= s.length()-1) { //read all regions
            if (s[p] == '*') {reg_0 = n_reg; p++;}
            if (s[p] == '&') {reg_1 = n_reg; p++;}
            if (s[p] == '(') p++;
            regions[n_reg] = (u_region)0;

            while ((s[p] != ')') && (p < s.length())) {
                r = 0; ss = "";
                while (is_digit(s[p])) {ss = ss + s[p]; p++;}
                while (s[p] == ' ') p++;

                r = atoi(ss.c_str());
                regions[n_reg] |= ((u_region)1<<r); /*cout << hex2int(s[p]) << " "; */
                //p++;
            }
            p++;
            n_reg++;
            if (!regions[n_reg-1]) n_reg--;
            p++;
        }
    }
    //cout << endl;
}


// turn in a knot region on arc, current arc, orient - against orien. (1) or towards (0), direction = right/left
u_number Cknot::turn(u_number arc, s_small *o, s_small direction) {
   // if ((arc == 126) || (arc == 135))cout << "arc 126 " << "orientation = " << (*o ? "against" : "toward") << (direction==RIGHT ? " right" : " left");
    s_letter x = gw(arc+1-(*o)); // csrossing
    *o = (s_small)SIGNB(x) ^ (s_small)(*o) ^ (s_small)(int)BIT(signs,(u_sign)ABS(x)) ^ (s_small)1 ^ (s_small)direction; // new direction
    //if ((arc == 126) || (arc == 135))  cout << " ---> arc " << MOD((int)gw_pos(-x)-(*o),n2) << " " << (*o? "against" : "toward") << endl;
    return MOD((int)gw_pos(-x)-(*o),n2);
}
    
        
string Cknot::check_ok_reason() {

    string s = "";
    
    u_region rx = u_region(0), ro = u_region(0);
    u_sign cx = u_sign(0), co = u_sign(0);
    int i;
    
    // crossings
    for (i = 0; i < n2; i++) {
        co |= (u_sign)1 << ABS(gw_[i]);
        cx ^= (u_sign)1 << ABS(gw_[i]);
    }

    // regions
    for (i = 0; i < n_reg; i++) {
        ro |= regions[i];
        rx ^= regions[i];
    }
    
    if (co != (fbf_sign(n+1) ^(u_sign)1)) s += "no crossing ";
    if (cx != u_sign(0)) s += "missing crossing " + itoa(firstBitSet(rx)) + " ";
    
    if (ro != fbf_reg(n2)) { s += "no arc "; cout << ro << " " << fbf_reg(n2) << endl; }
    if (rx != u_region(0)) { s += "parity arc " + itoa(firstBitSet(ro)) + " "; cout << "arc parity: "<<rx << endl;}
    
    return s;
    
}
        
void Cknot::check_ok(bool also_regions) {
  
    u_region rx = u_region(0), ro = u_region(0);
    u_sign cx = u_sign(0), co = u_sign(0);
    int i;
    
    // crossings
    for (i = 0; i < n2; i++) {
        co |= (u_sign)1 << ABS(gw_[i]);
        cx ^= (u_sign)1 << ABS(gw_[i]);
    }
    
    // regions
    for (i = 0; i < n_reg; i++) {
        ro |= regions[i];
        rx ^= regions[i];
    }
    
    if (co != (fbf_sign(n+1) ^ (u_sign)1)) throw 20;
    if (cx != u_sign(0)) throw 20;
    if (ro != fbf_reg(n2)) throw 20;
    if (rx != u_region(0)) throw 20;
    
    if (!also_regions) return;
    
    for (i=0; i<n2;i++) {
        if (bin2reg(gen_reg(i,RIGHT)) == -1) throw 20;
        if (bin2reg(gen_reg(i,LEFT)) == -1) throw 20;
    }
    
}

// generates the regions of a knot
void Cknot::generateRegions() {

    s_small o;
    u_number arc; // orientation, current arc
    u_region usedArcsR = fbf_reg(n<<1), usedArcsL = fbf_reg(n<<1); // all possible arcs set to 1

    memset(regions,0,sizeof(regions)); // init
    n_reg = 0;
    
    while (usedArcsR | usedArcsL) {

        arc = firstBitSet(usedArcsR | usedArcsL);// find 1st available arc
        
        o = !BIT(usedArcsR,(u_region)arc); // 0 positive (towards orientation), - negative (opposite of orientation)
        
        while (1) { // travel through knot
            
            (o ? usedArcsL : usedArcsR) ^= ((u_region)1 << arc); // remove arc

            regions[n_reg] |= ((u_region)1 << arc); // add arc // problrms

            arc = turn(arc,&o,RIGHT); // find next arc

            if (BIT(regions[n_reg],arc)) {
                
                n_reg++;
                
                if (n_reg > n2+2) throw 20; // to many regions
                break;
            
            } // closed loop
        }
    }
    
}
    
// generates a region adjacent to arc on side side
u_region Cknot::gen_reg(u_number arc, s_small side) {
    u_region b = (u_region)0;
    s_small o = side;
    do { // travel through knot
        b |= ((u_region)1 << arc); // add arc
        arc = turn(arc,&o,RIGHT); // find next arc
    } while (!BIT(b,arc));
    return b;
}

// current arc, orient - travelling against orientation (1) or towards (0), direction = right/left
// TODO: use new funcitons instead of this one
u_number Cknot::turn2(u_number arc, s_small *o, s_small direction) {
    s_letter x = gw(arc+1-(*o)); // crossing
    *o = SIGNB(x) ^ (*o) ^ (int)BIT(signs,ABS(x)) ^ 1 ^ direction; // new direction
    return MOD((int)gw_pos(-x)-(*o),n2);
}

// TODO: use new functions instead of this one
u_region Cknot::gen_reg2(u_number arc, s_small side) {
    u_region b = (u_region)0;
    s_small o = side;
    do { // travel through knot
        b |= (u_region)((u_region)1 << arc); // add arc
        arc = turn2(arc,&o,RIGHT); // find next arc
    } while (!BIT(b,arc));
    return b;
}

// get region adjacent to arc on side side
// TODO: use new functions instead of this one
u_number Cknot::regionGet(u_number arc, s_small side) {
    u_region b = gen_reg(arc,side);
    for (u_number r=0;r<n_reg;r++) if (regions[r] == b) return r;
    //cout << this << endl <<"region of ARC "<<(int)arc << " side " << (side==LEFT ? "left " : "right ") << "should have region " << b << endl;
    throw 20;
}


#endif
