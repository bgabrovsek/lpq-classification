//  tree.h
//  Created by Boštjan on 4/22/13.
//
//  A tree class with integers nodes used for manipulating the KBSM of arrow diagrams
//
//  Example : Ctree T = Ctree("[[ 3 [4 5] [5 [6 0]]]]");

#ifndef Lpq_tree_h
#define Lpq_tree_h

#include <set>
#include <list>

class Ctree;

using namespace std;

typedef list<Ctree *> t_children; // list of children of a node
typedef list<Ctree *>::iterator t_c_it;

// TODO: use existing next_i funcion
int next_i(string s, int *i) { // gets next int in s, starting at i, i points to next position
    int sign = 1, a = 0;
    while (s[*i] == ' ') (*i)++;
    if (s[*i] == '-') {sign=-1; (*i)++; }
    while (s[*i] == ' ') (*i)++;
    if (is_letter(s[*i])) return sign; // ?
    while ((is_digit(s[*i])) && (*i <= s.length())) a = a*10 + char2int(s[(*i)++]);
    while (s[*i] == ' ') (*i)++;
    return sign*a;
}

// TODO: use existing replace string function
void replace_s(std::string& subject, const std::string& search, const std::string& replace) { // search & replace all
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) { subject.replace(pos, search.length(), replace); pos += replace.length(); }
}


// Tree class

class Ctree {
    
public:
    int val;
    
    t_children child;
    
    Ctree * parent;
    
    bool is_leaf() { return child.empty(); }
    
    // depth of the tree
    int depth() { int d_, d=0; for (t_c_it it=child.begin();it!=child.end();it++) if ((d_=(*it)->depth()) > d-1) d = d_+1; return d;}
    
    // count the number of times v occours in the tree
    int count_vals(int v) {int c=0; for (t_c_it it=child.begin();it!=child.end();it++) if ((*it)->val == v) c++; return c; }
    
    // get the deepest child
    Ctree * deepest_child() {int d=depth(); if (d==0) return this; for (t_c_it it=child.begin();it!=child.end();it++) if ((*it)->depth()==d-1) return (*it)->deepest_child(); throw 20; } // not optimized
    
    // get the first leaf of the tree
    Ctree * a_leaf() { return (is_leaf() ? this : (*(child.begin()))->a_leaf());}
    
    void create_sibling(int v) { parent->child.push_back(new Ctree(v, parent)); }
    void create_uncle(int v) { parent->parent->child.push_back(new Ctree(v, parent->parent)); }
    
    bool has_grandsons() { for (t_c_it it = child.begin(); it != child.end(); it++) if (!(*it)->is_leaf()) return true; return false; }
    
    // finds son with val = v
    Ctree * find_first_son(int v) { for (t_c_it it = child.begin(); it != child.end(); it++) if ((*it)->val == v) return *it; return NULL; }
 
    // finds son with (!) lo <= val <= hi
    Ctree * find_first_son_interval(int lo, int hi, bool inside = true) {for (t_c_it it=child.begin();it!=child.end();it++) if ((((*it)->val>=lo)&&((*it)->val<=hi))^(!inside)) return *it; return NULL; }

    // constructors
    Ctree(Ctree *dad = NULL) { val = 0; parent = dad;}
    Ctree(int v, Ctree *dad = NULL) { val = v; parent = dad;}
    Ctree(string s, int *i_ = NULL, Ctree *dad = NULL); 
    
    Ctree * deep_copy();
    
    //~Ctree() { if (parent != NULL) { parent->child.erase(&(*this)); }}
};

// main constructor, e.g. s = "[[ 3 [4 5] [5 [6 0]]]]"
Ctree::Ctree(string s, int *i_, Ctree *dad) {
    replace_s(s," ","");
    cout << s << endl;
    int i = (i_ == NULL ? 0 : *i_);
    if (s[i++] != '[') throw 20;
    val = next_i(s,&i);
    parent = dad;
    while (s[i] == '[') child.push_back(new Ctree(s,&i,this));
    if (s[i++] != ']') throw 20;
    if (i_ != NULL) *i_ = i;
}

// copy the tree structure
Ctree * Ctree::deep_copy() {
    Ctree *t_, *t = new Ctree(val);
    for (t_c_it it=child.begin();it!=child.end();it++) {
        t_ = (*it)->deep_copy();
        t_->parent = t;
        t->child.push_back(t_);
    }
    return t;
}

// print tree
ostream &operator<<(ostream &out, Ctree * T) { // print tree
    out << "[" << T->val;
    for (list<Ctree *>::iterator it = T->child.begin(); it != T->child.end(); it++)
        out << " "<<*it;
    out << "]";
    return out;
}

#endif
