//  multivariate_laurent.h
//  Created by Boštjan on 4/22/13.
//
//  Multivariable laurent polyomial library for up to 256 variables

#ifndef Lpq_multivariate_laurent_h
#define Lpq_multivariate_laurent_h

// TODO: modify to a stand-alone library without common.h
#include "common.h"

#include <list>
#include <string>
#include <string.h>
#include <cmath>

#define SIMPLIFY_ON_ADDITION true // simplify polynomial after adding polynomials?
#define SIMPLIFY_ON_MULTIPLICATION true // simplify polynomial after multiplying polynomials?

bool LATEX = false; // print polynomial in latex mode
bool GROUP_OUTPUT_BY_FIRST_VARIABLE = true; // group output by first variable? e.g. (1+z)x + x^3 vs. x + xz + x^3
bool HOMFLY_STYLE = false; // use variables "t_3", "t_2",... instead of capitals, TODO: generalize this functionality for a stand-alone library

using namespace std;

// global variables
// TODO: do not use global variables

int n_v; // number of variables, global
char * variables; // names of variables (x,y,z,...), global
int reverse_variable[256]; // stores position of variable a, b, c,...

// class of a single term in the polynomial
// TODO: avoid public variables
class Cterm {
public:
    int coeff; // coedfficient of the term
    int * power; // list of powers of the variables

    //constructors
    // zero term
    Cterm() {power = new int[n_v]; coeff = 0; memset(power,0,n_v*sizeof(int)); /*  for (int i=0;i<n_v;i++) power[i] = 0; */ }
    
    // constant term
    Cterm(int cf) {power = new int[n_v]; coeff = cf; memset(power,0,n_v*sizeof(int));}
    
    // term from string, "-3x^3y^2"
    Cterm(string s);
    
    // term cf*varaibl^power
    Cterm(int cf, char variable, int pwr) {power = new int[n_v]; coeff = cf; memset(power,0,n_v*sizeof(int)); power[reverse_variable[variable]] = pwr;}

    // compares if powers are equal (disregards coefficient)
    bool eq_pow(Cterm *q) { for(int i=0;i<n_v;i++) if (power[i] != q->power[i]) return false; return true;}
    
    // compare if terms are equal
    bool eq_term(Cterm *q) { if (coeff != q->coeff) return false; for(int i=0;i<n_v;i++) if (power[i] != q->power[i]) return false; return true;}
    
    // compares if capital terms "ABC,..." are equal (for HOMFLY), TODO: delete or generalize this
    bool eq_capital_term(Cterm *q) { for(int i=0;i<n_v;i++) if ((power[i] != q->power[i]) && (is_capital(variables[i]))) return false; return true;}
    
    // TODO: rename to deep_copy
    Cterm * copy() { Cterm *t = new Cterm(); for(int i=0;i<n_v;i++) t->power[i] = power[i]; t->coeff = coeff; return t;}

    // get only upper factors of the term, TODO: delete or generalize
    Cterm *upper() {Cterm *t=copy(); t->coeff=1; for(int i=0;i<n_v;i++) if (!is_capital(variables[i])) t->power[i]=0; return t;}
    
    // get factor "x", TODO: delete or generalize
    Cterm *x() {Cterm *t=copy(); t->coeff=1; for(int i=0;i<n_v;i++) if (variables[i]!='x') t->power[i]=0; return t;}
    
    // does the term have upper factors "ABC...", TODO: delete or generalize
    bool has_upper() {for(int i=0;i<n_v;i++) if ((power[i] > 0) && is_capital(variables[i])) return true; return false;}
    
    // does the term have a "x" factor? // TODO: delete or generalize
    bool has_x() {for(int i=0;i<n_v;i++) if ((power[i] > 0) && (variables[i]=='x')) return true; return false;}
    
    // operator overloads
    
    // assignment operator
    Cterm * operator = (const Cterm & q); // creates new Cterm and copies data
    
    // compound assignment operators
    void operator *= (const int & q) {coeff *= q;} // scalar mul
    void operator *= (const Cterm & q) {coeff *= q.coeff; for(int i=0;i<n_v;i++) power[i] += q.power[i];} // mul, keeps q
    void operator /= (const Cterm & q) {coeff /= q.coeff; for(int i=0;i<n_v;i++) power[i] -= q.power[i];} // mul, keeps q

    // binary operators
    Cterm * operator *(const Cterm& q) { Cterm *t = new Cterm; t->coeff = coeff*q.coeff; for(int i=0;i<n_v;i++) t->power[i] = power[i] + q.power[i]; return t;}

    //destructor
    ~Cterm() {delete(power); }
};

// multivariable laurent polynomial class
// TODO: avoid public variables
class Cmultivariate {
public:

    //list of terms
    list<Cterm*> terms;

    //constructor
    Cmultivariate() {} // empty or zero polynomial
    Cmultivariate(string s); // string to polynomial "5xy -4yz^2x^-1 +5"

    // simplifies the polynomial (joins terms and sorts)
    void simplify();

    // deep copy the entire polynomial
    Cmultivariate * deepCopy();
    
    // operator overloads

    // compound assignment operators
    void operator += (Cmultivariate &); // addition, duplicates q
    void operator -= (Cmultivariate &); // substraction, duplicates q
    void operator *= (Cmultivariate &); // substraction, duplicates q
    void operator *= (const int &); // scalar multiplication
    void operator ^= (int P) {if(P) {Cmultivariate *s = deepCopy(); for (int i=1;i<P;i++) *this *= *s; delete s;} else {clear_terms(); terms.push_front(new Cterm(1));}}
    void operator *= (Cterm & t) {for (list<Cterm*>::iterator it=terms.begin(); it!=terms.end(); it++) *(*it) *= t;}
    void operator += (Cterm & t) {terms.push_back(&t); simplify();}
    void operator += (Cterm * t) {terms.push_back(t); simplify();}
    
    // binary operators
    Cmultivariate * operator * (Cmultivariate & q) {Cmultivariate * r = deepCopy(); (*r) *= q; return r;}
    Cmultivariate * operator + (Cmultivariate & q) {Cmultivariate * r = deepCopy(); (*r) += q; return r;}
    Cmultivariate * operator ^ (int & p) {Cmultivariate * r = deepCopy(); (*r) ^= p; return r;}
    Cmultivariate * operator * (Cterm & t) {Cmultivariate * r = deepCopy(); (*r) *= t; return r;}

    // max power of variable c in polynomial
    int max_power(char v) {int max=-100; for(list<Cterm*>::iterator it=terms.begin();it!=terms.end();it++)
        if ((*it)->power[reverse_variable[v]]>max) max=(*it)->power[reverse_variable[v]]; return max; }

    // min power of variable c in polynomial
    int min_power(char v) {int min=100; for(list<Cterm*>::iterator it=terms.begin();it!=terms.end();it++)
        if ((*it)->power[reverse_variable[v]]<min) min=(*it)->power[reverse_variable[v]]; return min; }

    Cmultivariate * extract(char c, int power); // extract "6+y" from "(6+y)x^2 + 65y" (c=x, power=2)

    // varaible multiplications
    // TODO: use operators, not seperate functions
    Cmultivariate *  multiplyByVariable(int c, char v, int p); // e.g. (-4,'x', 5) multiplies the polynomial by -4x^5
    Cmultivariate *  multiplyByVariable(char v, int p) { return multiplyByVariable(1,v,p); }
    Cmultivariate *  multiplyByVariable(char v) { return multiplyByVariable(1,v,1); }
    
    // replaces evert occurrence of variable v with poly
    void replace_factor_with_poly(char v, Cmultivariate * poly);
    // replace all terms that are divisible by t with poly
    void replace_exact_subfactor_with_poly(Cterm *t, Cmultivariate * poly);
    
    // clear poly to empty poly 
    void clear_terms() { for (list<Cterm*>::iterator it=terms.begin(); it!=terms.end(); it++) delete (*it);  terms.clear();}
    
    void zerofy() { clear_terms(); }

    // destructor
    ~Cmultivariate() { clear_terms(); /*terms.clear();*/}

private:

};

// does poly have the capitals terms equal to those in t
// TODO: delete or generalize
bool has_capital_term(Cmultivariate *poly, Cterm *t) {
    
    for (list<Cterm*>::iterator it=poly->terms.begin(); it!=poly->terms.end(); it++)
        if ( (*it)->eq_capital_term(t) ) return true;
    
    return false;

}

// true if term has no capital factors, false otherwise
// TODO: delete or generalize
bool is_affine(Cterm *t) {
    for(int i=0;i<n_v;i++) if ((t->power[i] != 0) && (is_capital(variables[i]))) return false;
    return true;
}

// true if poly has no terms with capital factors, false otherwise
// TODO: delete or generalize
bool is_affine(Cmultivariate * poly) {
    poly->simplify();
    for (list<Cterm*>::iterator it=poly->terms.begin(); it!=poly->terms.end(); it++)
        if ( !is_affine(*it) ) return false;
    return true;
}

// TODO: check for memory leaks
void Cmultivariate::replace_factor_with_poly(char v_factor, Cmultivariate * poly) {
    for (list<Cterm*>::iterator it=terms.begin(); it!=terms.end(); it++) {
        if ((*it)->power[reverse_variable[v_factor]] > 0) {
            Cmultivariate * q = poly->deepCopy();
            *q ^= (*it)->power[reverse_variable[v_factor]];
            (*it)->power[reverse_variable[v_factor]] = 0;
            *q *= *(*it);
            terms.splice(terms.end(), q->terms);
            delete q;
            (*it)->coeff = 0;
        }
    }
    simplify();
}

void Cmultivariate::replace_exact_subfactor_with_poly(Cterm *t, Cmultivariate * poly) { // leaks?
    
    bool b;
    Cmultivariate * q;
    
    if (t->coeff != 1) throw 20;
    for (list<Cterm*>::iterator it=terms.begin(); it!=terms.end(); it++) {
        
        b = true;
        for(int i=0;i<n_v;i++)
            if ((t->power[i] > 0) && ( t->power[i] != (*it)->power[i])) b = false;
        if (b) { // replacement should be done
            q = poly->deepCopy(); 
            
            for(int i=0;i<n_v;i++)
                if (t->power[i] > 0) (*it)->power[i] = 0; // in current factor, put variable powers to zero
            *q *= *(*it); // multiply poly by remaining factors
            terms.splice(terms.end(), q->terms); // add to poly
            
            delete q;
            
            (*it)->coeff = 0; // delete current term
        }
    }
    simplify();
}


void clear_Cterms(list<Cterm*> *terms) {
    for (list<Cterm*>::iterator it=terms->begin(); it!=terms->end(); it++) delete (*it);  terms->clear();
}

// comparison operators

// terms A = B?
bool operator==(Cterm &A, Cterm &B) {
    if (A.coeff != B.coeff) return false;
    for(int i=0;i<n_v;i++) if (A.power[i] != B.power[i]) return false;
    return true;
}

// terms A != B?
bool operator!=(Cterm &A, Cterm &B) {
    //cout << "term ";
    if (A.coeff != B.coeff) return true;
    for(int i=0;i<n_v;i++) if (A.power[i] != B.power[i]) return true;
    return false;
}

Cmultivariate * Cmultivariate::deepCopy() {
    Cmultivariate * p = new Cmultivariate();
	for (list<Cterm*>::iterator it=terms.begin(); it != terms.end(); it++) {
        p->terms.push_back((*it)->copy());
    }
    return p;
}


// initialization function, must be called before initializing the Cmultivatiable class
// argument s: string of variables, e.g. "xyz"
void initMultivariateLaurent(string s) {
    n_v = (int)s.length();
    variables = new char[n_v];
    memset(reverse_variable, 0, 256 * sizeof reverse_variable[0]);
    for (int i=0; i<n_v; i++) {variables[i] = s[i]; reverse_variable[s[i]] = i;};
}

// returns first possible integer and remove it, ex. stoi_r("123abc") = 123 (abc)
int stoi_r(string *s)
{
	int i = 0, sign = 1;
	if ((*s)[0] == '-') {
		sign = -1; s->erase(0,1);
	} else
        if ((*s)[0] == '+') s->erase(0,1);

	while (isdigit((*s)[0])) {
		i = i*10 + ((*s)[0]- '0');
		s->erase(0,1);
	}
    return i/(sign);
}


// compare two terms, for sorting
bool compareTerm(Cterm *a, Cterm *b) {
    for (int i=0; i<n_v; i++) // order by first variable
        if (a->power[i] != b->power[i]) return a->power[i] > b->power[i];
	return (a->coeff > b->coeff);

}

// constructor, parse string term
Cterm::Cterm(string s) {
    char c;
    power = new int[n_v];
    for (int i=0;i<n_v;i++) power[i] = 0;

    if (s[0] == '+') s.erase(0,1);
    if (is_letter(s[0])) coeff = 1; // if starts with letter, coeff = 1
    else
        if ((s[0] == '-') && (is_letter(s[1]))) {coeff = -1; s.erase(0,1);}
    else
    coeff = stoi_r(&s);

    while (s.length() > 0) {
        c = s[0]; s.erase(0,1);
        if (s[0] == '^') {
            s.erase(0,1);
            power[reverse_variable[c]] = stoi_r(&s);
        } else power[reverse_variable[c]] = 1;
    }
}

Cmultivariate * Cmultivariate::multiplyByVariable(int c, char v, int p) {
    for (list<Cterm*>::iterator it=terms.begin(); it != terms.end(); it++) {
        (*it)->power[reverse_variable[v]] += p;
        (*it)->coeff *= c;
    }
    return this;
}


Cmultivariate * Cmultivariate::extract(char c, int pw) {
    Cmultivariate *p = new Cmultivariate();

    for (list<Cterm*>::iterator it=terms.begin(); it != terms.end(); it++)
        if ((*it)->power[reverse_variable[c]] == pw)
            p->terms.push_back((*it)->copy());
    for (list<Cterm*>::iterator it=p->terms.begin(); it != p->terms.end(); it++) ((*it)->power)[reverse_variable[c]] = 0;
    return p;
}

// constructor, parse string polynomial, e.g. "2xy + 3x^-1 -8y^3z^-4"
Cmultivariate::Cmultivariate(string s) {
    string x;
    Cterm *t;
    int p; // position of ' '

    s.append(" ");

    for (int i=0;i<s.length();i++)
        if ((s[i] == '-') || (s[i] == '+') || (s[i] == ' ') || (s[i] == '^'))
            while (s[i+1] == ' ') s.erase(i+1,1);

    while (s.length() > 0) { // subdivide by spaces
        p = (int)s.find(' ');
		x = s.substr(0,p);
		s.erase(0,p+1);

        if ((x.length() > 0) && (x[0] != ' ')) {
            t = new Cterm(x);
            terms.push_back(t);
        }
    }
}

ostream &operator<<(ostream &out, Cmultivariate * p);

void Cmultivariate::simplify() {

	list<Cterm*>::iterator it, it_;
    
    // sort terms
    terms.sort(compareTerm);
    
	// group similar terms
	for (it=terms.begin(); it != terms.end(); it++) {
        it_ = it; it_++;
        while ((it_ != terms.end()) && ((*it)->eq_pow(*it_))) { (*it)->coeff += (*it_)->coeff; delete (*it_); it_ = terms.erase(it_); }
	}

	// delete zero terms
	for (it=terms.begin(); it != terms.end();)
        if ((*it)->coeff == 0) {delete *it; it = terms.erase(it);} else it++;
    
    
}

// Operator overloads

Cmultivariate * operator * (Cmultivariate p, int q) { // multiply polynomial with constant, creates new poly
    Cmultivariate *w = p.deepCopy();
    *w *= q;
    return w;
}

void Cmultivariate::operator += (Cmultivariate &q) {
    Cmultivariate *p = q.deepCopy();
    terms.splice(terms.end(), p->terms);
	delete p;
    //terms.splice(terms.end(), q.terms);
    if (SIMPLIFY_ON_ADDITION) simplify();
}

void Cmultivariate::operator -= (Cmultivariate &q) {
    
    terms.splice(terms.end(), (q*(-1))->terms);
	if (SIMPLIFY_ON_ADDITION) simplify();
}

void Cmultivariate::operator *= (const int & q) {
    for (list<Cterm*>::iterator it=terms.begin(); it != terms.end(); it++)
        (**it) *= q;
}

void Cmultivariate::operator *= (Cmultivariate &q) // poly multiplication
{
    list<Cterm*> a;
	list<Cterm*>::iterator ita, itb;

    a.splice(a.end(),terms); // copy terms to a

	for (ita=a.begin(); ita!=a.end(); ita++)
		for (itb=q.terms.begin(); itb!=q.terms.end(); itb++)
            terms.push_back((**ita) * (**itb));
    clear_Cterms(&a);
	if (SIMPLIFY_ON_MULTIPLICATION) simplify();

}

// Global overloads

bool CAPITALS_FIRST = true;

// term to string stream
ostream &operator<<(ostream &out, Cterm * t) {
    bool z = true;
    if (t->coeff < 0) out << "- ";
    // are all powers 0 ?
    for(int i=0;i<n_v;i++) z &= (t->power[i] == 0);
    //print coefficient
    if ((abs(t->coeff) != 1) || z) out << abs(t->coeff);
    if (t->coeff) {
        if (CAPITALS_FIRST) {
            for(int i=0;i<n_v;i++)
            if (t->power[i]!=0) {
                if ((!HOMFLY_STYLE) || (!is_capital(variables[i])))
                    out << variables[i];
                else out << (LATEX ? HOMFLY_BASE_LATEX[variables[i]] : HOMFLY_BASE[variables[i]]);
                    
                if (t->power[i] != 1)
                out << (LATEX?"^{":"^") << t->power[i] << (LATEX?"}":"");
               // out << "*";
            }
        } else {
            
            for(int i=0;i<n_v;i++)
                if (t->power[i]!=0) {
                    if (is_capital(variables[i])) continue;
                    out << variables[i];
                    if (t->power[i] != 1) out << (LATEX?"^{":"^") << t->power[i] << (LATEX?"}":"");
                }
            
            for(int i=0;i<n_v;i++)
                if (t->power[i]!=0) {
                    if (!is_capital(variables[i])) continue;
                    out << (LATEX ? HOMFLY_BASE_LATEX[variables[i]] : HOMFLY_BASE[variables[i]]);
                    if (t->power[i] != 1) out << (LATEX?"^{":"^") << t->power[i] << (LATEX?"}":"");
                }
        }
    }
   // out << "1";
    return out;
}

// pointer to Cmultivariate to string stream
ostream &operator<<(ostream &out, Cmultivariate * p) { // print pointer to Cmultivariate
    if (p->terms.size() == 0) out << 0;
    else
    if (GROUP_OUTPUT_BY_FIRST_VARIABLE) {

        int power = 9999; // ?? TODO: check why this is set to 9999
        Cmultivariate * temp;

        p->terms.sort(compareTerm); // order by 1st variable
        for (list<Cterm*>::iterator it= p->terms.begin(); it != p->terms.end(); it++)
        {
            if ((*it)->power[0] != power) { // new degree of 1st variable?
                if (it != p->terms.begin()) out << " + ";
                temp = p->extract(variables[0], power = (*it)->power[0]);
                GROUP_OUTPUT_BY_FIRST_VARIABLE = false;
                if (power) out << "(" << temp << ")" << variables[0] << ((power!=1)?"^":""); else out << temp;
                if ((power) && (power!=1)) out << power;
                  //  cout << "*";
                GROUP_OUTPUT_BY_FIRST_VARIABLE = true;
                delete temp;
            }
        }
    }

    else // just print, don't group

    for (list<Cterm*>::iterator it= p->terms.begin(); it != p->terms.end(); it++) {
        if (it != p->terms.begin()) out << (((*it)->coeff >= 0) ? " + " : " ");
        out << *it;
    }

    return out;
}

// true if b divides a, false otherwise
// TODO: optimize, generalize
bool divisible(Cmultivariate *a, Cmultivariate *b) {
    Cmultivariate * numerator;
    Cmultivariate * divisor;
    Cterm * numerator_factor;
    Cterm * divisor_factor;
    
    numerator = a->deepCopy();
    
    unsigned long max = MAXIMUM(a->terms.size(),b->terms.size());
    max = a->terms.size() + b->terms.size();
    for (int m = 0; m < max; m++) {
        
        divisor = b->deepCopy();
        
        numerator_factor = (*(numerator->terms.begin()))->copy();
        divisor_factor = (*(divisor->terms.begin()))->copy();;
        
        *numerator_factor /= *divisor_factor;
        *divisor *= *numerator_factor;
        *numerator -= *divisor;
        numerator->simplify();

       // delete divisor;
        delete numerator_factor;
        delete divisor_factor;
        
        if (numerator->terms.size() <= 0) return true;
    }
    return false;
}

// Cmultivariate to string stream
ostream &operator<<(ostream &out, Cmultivariate & p) {
    out << &p; return out;
}
          
// in the case t is a term of the HSM, reverse it t^n -> t^-n
// TODO: delete or generalize
Cterm * reverse_orientation(Cterm * t) {

    Cterm *t_ = t->copy();
    int i;
    
    for(i=0;i<n_v;i++) if (is_capital(variables[i])) t_->power[i] = 0;
    
    for(i=0;i<n_v;i++)
        if (is_capital(variables[i])) {
            (t_)->power[ reverse_variable[ 2*MID_VAR - variables[i]] ] = (t)->power[i];
        }
    return t_;
}
 
// in the case t is a polynomial of the HSM, reverse terms:  t^n -> t^-n
// TODO: delete or generalize
Cmultivariate * reverse_orientation(Cmultivariate *A) {

    Cmultivariate *B = A->deepCopy();
    
    list<Cterm*>::iterator ita, itb;
    
    for (ita=A->terms.begin(), itb=B->terms.begin(); ita!=A->terms.end(); ita++, itb++) {
        
        for(int i=0;i<n_v;i++) if (is_capital(variables[i])) (*itb)->power[i] = 0;
        
        for(int i=0;i<n_v;i++)
            if (is_capital(variables[i]))
                (*itb)->power[ reverse_variable[ 2*MID_VAR - variables[i]] ] = (*ita)->power[i];
            
    }
    return B;
}

// polynmial A = polynomial B ?
// assume polynomial are simplified (i.e. in their canonical forms)
bool operator==(Cmultivariate &A, Cmultivariate &B) {
    list<Cterm*>::iterator ita, itb;
    //A.simplify();
    //B.simplify();
    if (A.terms.size() != B.terms.size()) return false;
    for (ita=A.terms.begin(), itb=B.terms.begin(); ita!=A.terms.end(); ita++, itb++)
        if ((**ita) != (**itb)) return false;
    return true;
}
                    
// A = B up to change of variables t_i <-> t_-i
// TODO: generalize or delete
bool poly_equal_orientation(Cmultivariate *A, Cmultivariate *B, bool simplify_A = true, bool simplify_B = true) { 
    if (simplify_A) A->simplify();
    if (simplify_B) B->simplify();
    if (*A == *B) return true;
    Cmultivariate *C = reverse_orientation(B);
    C->simplify();
    if (*A == *C) {delete C; return true; }
    delete C;
    return false;
}
                    
// A = B up to change of variables t_i <-> t_-i in L(p,q)
// TODO: delete
bool poly_equal_orientation_lpq(int l, Cmultivariate *poly_a, Cmultivariate *poly_b, bool simplify_A = true, bool simplify_B = true) {
    Cmultivariate *b_a, *b_b, *b_b_, *poly;
    b_a = HSM_lpq(l,poly_a);
    b_b = HSM_lpq(l,poly_b);
    b_a->simplify();
    b_b->simplify();
    if (*b_a == *b_b) return true;
    poly = reverse_orientation(poly_b);
    b_b_ = HSM_lpq(l, poly);
    b_b_->simplify();
    if (*b_a == *b_b_) return true;
    return false;
}
 
// comparison operator
bool operator!=(Cmultivariate &A, Cmultivariate &B) {
    
    list<Cterm*>::iterator ita, itb;
    
    A.simplify();
    B.simplify();
    
    if (A.terms.size() != B.terms.size()) return true;
    
    for (ita=A.terms.begin(), itb=B.terms.begin(); ita!=A.terms.end(); ita++, itb++)
        if ((**ita) != (**itb)) return true;
    return false;
}


#endif
