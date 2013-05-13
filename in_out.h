//  in_out.h
//  Created by Boštjan on 4/22/13.
//
//  input/output functions: saving/loading data from/to files, printing on screen

#ifndef Lpq_in_out_h
#define Lpq_in_out_h

#include <string>
#include <fstream>
#include <vector>
#include <time.h>

#include "multivariate_laurent.h"
#include "classification_data.h"
#include "KBSM.h"
#include "algorithms.h"
#include "knot.h"
#include "HOMFLY.h"
#include "global_vars.h"
#include "lens.h"

using namespace std;

// GENERAL FUNCTIONS

// read file by lines and place lines into a vector
vector<string> readStringsFromFile(string f_name) {
    vector<string> vs;
    ifstream ifs(f_name.c_str(), ifstream::in);
    string temp;
    while (getline(ifs, temp)) vs.push_back(temp);
    return vs;
}

// "search & replace all"
void replace_strings(std::string& subject, const std::string& search, const std::string& replace) { 
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) { subject.replace(pos, search.length(), replace); pos += replace.length(); }
}

// convert HSM generators from human to computer readable form
void human_to_computer(std::string &s) {
    replace_strings(s, "t_10", "W"); replace_strings(s, "t_11", "X"); replace_strings(s, "t_12", "Y");
    replace_strings(s, "t_0", "M"); replace_strings(s, "t_1", "N"); replace_strings(s, "t_2", "O");
    replace_strings(s, "t_3", "P"); replace_strings(s, "t_4", "Q"); replace_strings(s, "t_5", "R");
    replace_strings(s, "t_6", "S"); replace_strings(s, "t_7", "T"); replace_strings(s, "t_8", "U");
    replace_strings(s, "t_9", "V");
    
    replace_strings(s, "t_-10", "C"); replace_strings(s, "t_-11", "B"); replace_strings(s, "t_-12", "A");
    replace_strings(s, "t_-1", "L"); replace_strings(s, "t_-2", "K"); replace_strings(s, "t_-3", "J");
    replace_strings(s, "t_-4", "I"); replace_strings(s, "t_-5", "H"); replace_strings(s, "t_-6", "G");
    replace_strings(s, "t_-7", "F"); replace_strings(s, "t_-8", "E"); replace_strings(s, "t_-9", "D");
}

// KNOT FUNCTIONS

// save all knots to file
void saveKnotsToFile(t_knot_set *KNOTS, string s) {
    ofstream file;
    file.open (s.c_str());
    bool save_bitset_format = use_extended_bitset; // if printing in long form
    use_extended_bitset = false;
    for (t_knot_set::iterator it=KNOTS->begin(); it!=KNOTS->end(); it++) file << (*it) << endl;
    use_extended_bitset = save_bitset_format;
    file.close();
}

// load all knot gauss words from file
void loadKnotsFromFile(t_knot_set *KNOTS, string s) {
    vector<string> sKnots = readStringsFromFile(s);
    Cknot *K;
    for (vector<string>::iterator it=sKnots.begin(); it!=sKnots.end(); it++) {
        K = new Cknot(*it);
        KNOTS->insert(K);
    }
}

// HOMLFY SKEIN MODULES

// save HSM generators to file
void save_generators_to_file(string s, bool write_if_null) {
    
    ofstream file;
    file.open (s.c_str());
    int l,g;
    
    for (l = 0; l<N_HOMEO; l++)
        for (g=0; g<torus_HSM_generator_n; g++)
        if ( (generator[l][g] != NULL) || (write_if_null))
            //if ((write_if_null)^(generator[l][g] == NULL))
        {
            file << "L(" << lens_homeo[l][0] << "," << lens_homeo[l][1] << "): ";
            file << torus_HSM_generator_terms[g] << " = ";
            if (generator[l][g] != NULL) {
                file << generator[l][g] << endl;
            } else file << "NULL" << endl;
        }
    file.close();
    
}

// load intermediate HSM generators
void load_generators(string s) { 
    
    vector<string> sKnots = readStringsFromFile(s);
    
    int pos, p, q, pos_;
    string ss;
    string ss_;
    Cterm *term;
    Cmultivariate *poly;
    
    for (vector<string>::iterator it=sKnots.begin(); it!=sKnots.end(); it++) {
        
        pos = 2;
        p = next_int(*it, &pos); pos++;
        q = next_int(*it, &pos); pos+= 3;
        
        pos_ = (int)(*it).find(" = ");
        
        ss = (*it).substr(pos,pos_-pos); pos_+=3; // generator
        ss_ = (*it).substr(pos_); // polynomial
        
        human_to_computer(ss);
        human_to_computer(ss_);
        
        term = new Cterm(ss);
        poly = new Cmultivariate(ss_);
        
        generator[get_lpq_index(p, q)][get_generator_index(term)] = poly;
        
        delete term;
    }
}

// saves a vector of HOMLFY skein modules to file
void save_HOMFLYs(t_poly_vector *HOMLFYS, string s) {
    ofstream file;
    file.open (s.c_str());
    int i = 0;
    for (t_poly_vector::iterator it=HOMFLYS.begin(); it!=HOMFLYS.end(); it++)
        file << "P(" << i++ << ") = " << (*it) << endl;
    file.close();
}


// CLASSIFICATION FUNCTIONS

// saves which knots are marked as primes and which are marked as duplicates (torus)
void saveClassification(string s) {
    ofstream file;
    file.open (s.c_str());
    for (int i = 1; i <= number_of_knots; i++) {
        file << i << " ";
        if (classification[i].prime) file << "prime ";
        if (classification[i].duplicate) file << "duplicate ";
        if (classification[i].composite) file << "composite "; // obsolete
        file << endl;
    }
    file.close();
}

// saves which knots are marked as primes and which are marked as duplicates (L(p,q))
void save_classification_lpq(string s) {
    ofstream file;
    file.open (s.c_str());
    
    for (int l = 0; l<N_HOMEO; l++) {
    
        for (int i = 1; i <= number_of_knots; i++) {
            file << lpq_s[l] << ": ";
            file << i << " ";
            if (classification_lpq[l][i].prime) file << "prime ";
            if (classification_lpq[l][i].duplicate) file << "duplicate ";
            if (classification_lpq[l][i].composite) file << "composite ";
        
            file << endl;
        }
    }
    
    file.close();
    
    
    // obsolete:
    /*
    generate_knot_names(); // 0_1, 1_1, 1_2,...
    file.open ("./visualizaton-lpq.txt");
    
    for (int k=1;k<n_prime;k++) {
        file << knot_name_ID[k] << ": ";
        if (classification[knot_name_ID[k]].prime) file << " "; else
        if (classification[knot_name_ID[k]].duplicate) file << " "; else
        file << "-" << " ";
        for (int l = 0; l<N_HOMEO; l++) if (consider_lpq[l]) {
            if (classification_lpq[l][knot_name_ID[k]].prime) file << "*"; else
                if (classification_lpq[l][knot_name_ID[k]].duplicate) file << "/"; else
                    file << "-";
            file << " ";
        }
        file << endl;
    }
    file.close();*/
}

// loads the classification of torus knots
void loadClassification(string s) {
    
    int p;
    int ID;
    
    clear_classification_data();
    vector<string> sKnots = readStringsFromFile(s);
    
    for (vector<string>::iterator it=sKnots.begin(); it!=sKnots.end(); it++) {
        p = 0;
        ID = next_int(*it, &p);
        classification[ ID ].prime = ((*it).find("prime") != string::npos);
        classification[ ID ].duplicate = ((*it).find("duplicate") != string::npos);
        classification[ ID ].composite = ((*it).find("composite") != string::npos); // obsolete
    }
}

// loads the classification of L(p,q) knots
void load_classification_lpq(string s) {
    int pos, p, q, l;
    int ID;
    string ss;
    
    clear_classification_data_lpq();
    vector<string> sKnots = readStringsFromFile(s);
    
    for (vector<string>::iterator it=sKnots.begin(); it!=sKnots.end(); it++) {
        //pos = 0;
        pos = 2;
        p = next_int(*it, &pos); pos++;
        q = next_int(*it, &pos); pos+= 3;
        ss = (*it).substr(pos); // regular

        l = get_lpq_index(p, q);
        
        pos = 0;
        ID = next_int(ss, &pos);

        classification_lpq[l][ ID ].prime = (ss.find("prime") != string::npos);
        classification_lpq[l][ ID ].duplicate = (ss.find("duplicate") != string::npos);
        classification_lpq[l][ ID ].composite = (ss.find("composite") != string::npos); // obsolete
    }
}

// loads and adds the classification of torus knots to the existing classification
void OR_loadClassification(string s) {
    int p, ID;
    vector<string> sKnots = readStringsFromFile(s);
    for (vector<string>::iterator it=sKnots.begin(); it!=sKnots.end(); it++) {
        p = 0;
        ID = next_int(*it, &p);
        classification[ ID ].prime |= ((*it).find("prime") != string::npos);
        classification[ ID ].duplicate |= ((*it).find("duplicate") != string::npos);
        classification[ ID ].composite |= ((*it).find("composite") != string::npos); // obsolete
    }
}

// print all knots in the set KNOTS, obsolete
ostream &operator<<(ostream &out, t_knot_set KNOTS) {
    for (t_knot_set::iterator it = KNOTS.begin(); it != KNOTS.end(); it++)
        out << *it << endl;
    return out;
}
        
ostream &operator<<(ostream &out, t_knot_set * KNOTS) {
    cout << *KNOTS;
    return out;
}
        
// print stats of knots in the group G, oboslete
ostream &operator<<(ostream &out, t_group_list * G) {
    int num_knots = 0, num_polys = 0, num_dupli = 0;
    for (t_group_list::iterator ig = G->begin(); ig != G->end(); ig++) {
        num_polys++;
        if ((*ig)->KNOTS->size() > 1) num_dupli++;
        for (t_knot_set::iterator iK = (*ig)->KNOTS->begin(); iK != (*ig)->KNOTS->end(); iK++) num_knots++;
    }
    out << "Knots: " << num_knots << endl;
    out << "Polys: " << num_polys;
    out << " (" << (num_polys-num_dupli) << " single + " << num_dupli << " duplicates)" << endl;
    return out;
}
        
// copy the file source to dest and append the current time to the file name (for making backup of files)
void backup_file(string source, string dest) {
    time_t rawtime;
    char * c_time;
    time ( &rawtime );
    c_time = ctime(&rawtime);
    dest += " - " + string(c_time) + ".txt";
    ifstream f1(source.c_str(), fstream::binary);
    ofstream f2(dest.c_str(), fstream::trunc|fstream::binary);
    f2 << f1.rdbuf();
}
        
// printing groups/sets of knots
        
// flags
#define GROUP_POLY (1<<0)
#define PRINT_KNOT_ID (1<<1)
#define PRINT_GAUSS (1<<2)
#define PRINT_KNOT_POLY (1<<3)
#define PRINT_GROUP_ID (1<<4)
#define PRINT_GROUP_POLY (1<<5)
#define PRINT_CLASSIFICATION (1<<6)
#define PRINT_GROUPED_BY_LINE (1<<7)
#define NO_SINGLE_GROUPS (1<<8)
#define NO_DUPLICATES (1<<9)
#define NO_COMPOSITES (1<<10)
#define NO_PRIMES (1<<11)

// general print with multiple options
void print_knots(t_group_list *G_HOMFLY, u64 flag) {

    if ((flag & GROUP_POLY) ||(flag & PRINT_GROUPED_BY_LINE) ) { // list groups
        
        for (t_group_list::iterator iG = G_HOMFLY->begin(); iG != G_HOMFLY->end(); iG++) { // loop through all polynomial groups
            
            if (flag & NO_SINGLE_GROUPS)
                if (unclassified_in_group((*iG)->KNOTS,classification) == 0) continue;
            
            if (flag & PRINT_GROUP_ID) cout << "["<<(*iG)->group_ID << "] "; //<<endl;
            if (flag & PRINT_GROUP_POLY) cout << (*iG)->poly << endl;
            
            
            for (t_knot_set::iterator iK = (*iG)->KNOTS->begin(); iK != (*iG)->KNOTS->end(); iK++) { // loop through all knots in group
                
                if ((flag & NO_PRIMES) && (classification[(*iK)->ID].prime)) continue;
                if ((flag & NO_DUPLICATES) && (classification[(*iK)->ID].duplicate)) continue;
                if ((flag & NO_COMPOSITES) && (classification[(*iK)->ID].composite)) continue;
                
                
                if (flag & PRINT_KNOT_ID) cout << (*iK)->ID << " ";
                if (flag & PRINT_GAUSS) cout << (*iK);
                if (flag & PRINT_CLASSIFICATION) {
                    if (classification[(*iK)->ID].prime) cout << " prime";
                    if (classification[(*iK)->ID].duplicate) cout << " duplicate";
                    if (classification[(*iK)->ID].composite) cout << " composite"; // obsolete
                    
                }
                if (!PRINT_GROUPED_BY_LINE) cout << endl;
            }
            cout << endl;
        }
        
                
    } else { // list all knots
        
        for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
            
            if ((flag & NO_PRIMES) && (classification[(*iK)->ID].prime)) continue;
            if ((flag & NO_DUPLICATES) && (classification[(*iK)->ID].duplicate)) continue;
            if ((flag & NO_COMPOSITES) && (classification[(*iK)->ID].composite)) continue;
            
            
            if (flag & PRINT_KNOT_ID) cout << (*iK)->ID << endl;
            if (flag & PRINT_GAUSS) cout << (*iK) << endl;
            if (flag & PRINT_KNOT_POLY) cout << HOMFLY((*iK)) << endl;
            if (flag & PRINT_CLASSIFICATION) {
                if (classification[(*iK)->ID].prime) cout << " prime";
                if (classification[(*iK)->ID].duplicate) cout << " duplicate";
                if (classification[(*iK)->ID].composite) cout << " composite"; // obsolete
                cout << endl;
            }
        }
        
    }
}
        
// for each unclassified knot in G calculate which knot in K it could be equal (=first knot in each group)
// TODO: use other methods for this
void calculate_representatives_lpq(int l, t_group_list *G) {
    int ID;
    for (t_group_list::iterator iG = G->begin(); iG != G->end(); iG++) {
        ID = (*(*iG)->KNOTS->begin())->ID;
        for (t_knot_set::iterator iK = (*iG)->KNOTS->begin(); iK != (*iG)->KNOTS->end(); iK++) {
            classification_equality_lpq[l][(*iK)->ID] = ID;
        }
    }
}

string output_knot_name(int kid);
        
// print knots that are not marked as duplicates (primes + unknowns)
void print_non_duplicates(t_group_list *G, bool single, Cknot_info classification[MAX_KNOTS], bool use_names = false) {
    bool first, first_g;
    first_g = true;
    for (t_group_list::iterator iG = G->begin(); iG != G->end(); iG++) {
        if (!single) if (classified_in_group((*iG)->KNOTS,classification) <= 1) continue;
        first = true;
        cout << (first_g ? "[" : ", ["); first_g = false;
        for (t_knot_set::iterator iK = (*iG)->KNOTS->begin(); iK != (*iG)->KNOTS->end(); iK++) {
            if (!single) if (classification[(*iK)->ID].duplicate) continue;
            cout << (first ? "" : ", ");
            if (use_names) cout << output_knot_name(((*iK)->ID)); else cout << (*iK)->ID;
            first = false;
        }
        cout << "]";
    }
}

// print integers i,j,k to log file, obsolete
void log(string s, int i = -1, int j = -1, int k = -1) {
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

// return string of knot id (in case we want to switch output by printing IDs or knots in the form n_k)
string output_knot_name(int kid) {
    return knot_name[name_from_id(kid)];
}
        
// print the main classification table, both torus knots and L(p,q) knots
// prints in latex format
void classification_table() {
    int i, j, k;
    //u16 property; // obsolete
    t_knot_set::iterator iK;
    t_group_list:: iterator iG;
    
    // generate_knot_names(false);
    generate_knot_names_bar(); // generate knot names 0_1, 1_1, 2_1, 2_2, ...
    print_knot_names(); // print the conversion table between knot names and knot IDs
    // init_prime_affine();
    
    // torus classification
    
    cout << endl << "% ---------------------------------"<< endl<< "Solid torus" << endl << endl;;

    cout << endl <<"Unclassified torus HSM:  ";
    print_non_duplicates(G_HOMFLY, false, classification, YES);
    
    cout << endl;
    cout << endl<< "Unclassified torus KBSM: ";
    print_non_duplicates(G_KBSM, false, classification, YES);
    
    // lens space classification
    
    for (int l=0; l < N_HOMEO; l++) {
        
        calculate_representatives_lpq(l, lens_homeo[l][1] <= 1 ? G_HOMFLY_LPQ[l] : G_KBSM_LPQ[l]);
        
        cout << endl << "% ---------------------------------"<< endl;
        cout << "\\section{$"<<lpq_s[l] <<"$}"<< endl << endl;;
        
        // duplicates
        
        cout << "Duplicates: " << endl;
      
      /*  for (int i_ = 0; i_ < n_prime; i_++)  {
            i = knot_name_ID[i_];
            if (classification_lpq[l][i].duplicate) {
                cout << output_knot_name(i)<< " & " <<output_knot_name(classification_equality_lpq[l][i])  << "";
                cout << endl;
               // if ((++cb) % 6) cout << "& "; else cout << "\\\\" << endl;
            }
        }*/
        for (int i_ = 0; i_ < n_prime; i_++)  {
            i = knot_name_ID[i_];
            
            if (classification_equality_lpq[l][i] != i) continue; // must be the smallest of its class
            
            k = 0; // count duplicates
            for (int j_ = 0; j_ < n_prime; j_++) {
                j = knot_name_ID[j_];
                if ((classification_lpq[l][j].duplicate) &&
                    (classification_equality_lpq[l][i] == classification_equality_lpq[l][j])) k++;
            }
            if (k > 0) {
                cout << " & ";
                for (int j_ = 0; j_ < n_prime; j_++) {
                    j = knot_name_ID[j_];
                    if ((classification_lpq[l][j].duplicate) &&
                        (classification_equality_lpq[l][i] == classification_equality_lpq[l][j])) cout << "$" <<output_knot_name(j) << "$, ";
                }
                
                cout << " & $" << output_knot_name(i) << "$"<< endl;
            }
           /* if (classification_lpq[l][i].duplicate) {
                cout << output_knot_name(i)<< " & " <<output_knot_name(classification_equality_lpq[l][i])  << "";
                cout << endl;
                // if ((++cb) % 6) cout << "& "; else cout << "\\\\" << endl;
            }*/
        }
        
        cout << endl << endl;

        cout << "Unknown status: " << flush;
        for (int i_ = 0; i_ < n_prime; i_++)  {
            i = knot_name_ID[i_];
            if (i_ != 0)
            if ((!classification_lpq[l][i].duplicate) && (!classification_lpq[l][i].prime)) cout << output_knot_name(i) << ", ";
        }
        
        // print (un)classified knots that share the HSM and/or KBSM
    
        if (lens_homeo[l][1] <= 2) {
            cout << endl<<endl<< "Unclassified by HSM:  ";
            print_non_duplicates(G_HOMFLY_LPQ[l], false, classification_lpq[l], YES);
        }
        
        cout << endl<<endl<< "Unclassified KBSM: ";
        print_non_duplicates(G_KBSM_LPQ[l], false, classification_lpq[l], YES);
        
        cout << endl;
    }
    
    cout << endl << endl << "NON-HOMEOMORPHIC LENS SPACES: " << endl;
    for (int l=0;l<N_HOMEO;l++) cout << "$" <<lpq_s[l] << "$" <<", ";
    
    cout << endl << endl << "NUMBER OF PRIME KNOTS: " << endl;
    cout << "Solid torus: ";
    for (i=0;i<=5;i++) cout << "(" << i << ") " <<( (i==0 ? 1 : count_prime_knots(KNOTS, classification, i)) ) << ", ";
    cout << endl << endl;

    for (int l = 0; l < N_HOMEO; l++) {
        cout << lpq_s[l] << (lens_homeo[l][0] > 9 ? "" :" ")<<": ";
        for (i=0;i<=5;i++) cout << "(" << i << ") " <<( (i==0 ? 1 : count_prime_knots(KNOTS, classification_lpq[l], i)) ) << ", ";
        cout << endl;

    }
    
    cout << endl << endl << "AFFINE TRIVIAL KNOTS: " << endl;
    
    Cmultivariate * p_t = new Cmultivariate("-a^2 - a^-2");
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
        
    }
    
    cout << endl << endl << "AFFINE KNOTS: " << endl;
    for (int l = 0; l < N_HOMEO; l++) {
        cout << lpq_s[l] << (lens_homeo[l][0] > 9 ? "" :" ")<<": ";
        for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
            if (classification[(*iK)->ID].duplicate) continue;
            if (get_knot_property(*iK, l) & PROPERTY_AFFINE) cout << output_knot_name((*iK)->ID) << ", ";
        }
        cout << endl;
    }
    
    cout << endl << endl << "CONNECTED SUMS: " << endl;
    for (int l = 0; l < N_HOMEO; l++) {
        cout << lpq_s[l] << (lens_homeo[l][0] > 9 ? "" :" ")<<": ";
        for (t_knot_set::iterator iK = KNOTS->begin(); iK != KNOTS->end(); iK++) {
            if (classification[(*iK)->ID].duplicate) continue;
            if (get_knot_property(*iK, l) & PROPERTY_SUM) cout << output_knot_name((*iK)->ID) << ", ";
        }
        cout << endl;
    }
    
}

// generates the latex knot table (in terms of PDFs)
void latex_tabulate_knot_table() {
    
    #define PGN 5 // number of pages
    
    int page[PGN][2] = {{4,6},{6,5},{6,5},{6,5},{5,5}}; // how many rows and column each page has
    
    int i, pg;
    int c,r,k;
    
    generate_knot_names(); // 0_1, 1_1, 1_2,...
    
    cout << endl;

    k = 0;
    for (pg=0; pg<PGN; pg++) {
        cout << "\\begin{center}";
        cout << "\\begin{tabular}{";
        for (i=0;i<page[pg][1];i++) cout << "c";
        cout << "}" << endl;
        
        for (r=0;r<page[pg][0];r++) {
            
            for (c=0;c<page[pg][1];c++)
                if (k+c < n_prime)
                    cout << "\\knot{" << knot_name_ID[k+c]+1 << "}" << ((c < page[pg][1]-1) ? " & " : ""); // +1 since unknot in 1st place
                else cout << ((c < page[pg][1]-1) ? " & " : "");
            cout << "\\\\" << endl;
            
            for (c=0;c<page[pg][1];c++)
                if (k+c < n_prime) {
                    if (chirality[knot_name_ID[k+c]] == 3) cout << "${^\\ast}"; else cout << "$";
                    cout << knot_name[k+c] << "$" << ((c < page[pg][1]-1) ? " & " : "");
                }
                else cout << ((c < page[pg][1]-1) ? " & " : "");
            
            k+=page[pg][1];
            if (r<page[pg][0]-1) cout << "\\\\[0.5cm]";
            cout << endl;
        }
        cout << "\\end{tabular}" << "\\end{center}" << endl << endl; //\\clearpage
    }
    cout << "Legend: $\\ast$ -- amphichiral knot." << endl << endl;
}
        

// extract the knot name from string, b = string includes a "bar"
string trim_string(string s, bool *b) {
    size_t f1 = s.find('(');
    size_t f2 = s.find(')');
    string sa = s.substr(f1+1, f2-f1-1); // trim
    *b = (sa.find("bar") != sa.npos); // does the string have a "bar"?
    if (*b) { // if it has a bar, trim it off
        sa.erase(0,5);
        sa.erase(sa.end()-1);
    }
    if (sa.find("*") != sa.npos) sa.erase(0,1);
    return sa;
}
        
// sort two strings by knot name, i.e. 
bool name_sort(string a, string b) {
            //    return (i<j);
    string sa, sb;
    bool bar_a, bar_b;
    int crs_a, crs_b, ind_a, ind_b, ia, ib;
    // get knot name and bar info
    sa = trim_string(a, &bar_a);
    sb = trim_string(b, &bar_b);
    
    // convert number of crossings and index to integers
    ia = ib = 0;
    crs_a = next_int(sa,&ia); ia+=2;
    ind_a = next_int(sa,&ia);
    crs_b = next_int(sb,&ib); ib+=2;
    ind_b = next_int(sb,&ib);
    // return
    if (crs_a != crs_b) return (crs_a < crs_b);
    if (ind_a != ind_b) return (ind_a < ind_b);
    if (bar_a == bar_b) return false;
    return !bar_a;
    //return (bar_a == bar_b ? false : (!bar_a));
}
        
struct less_than_name {
    inline bool operator() (const string& s1, const string& s2)
    {
        return name_sort(s1, s2);
    }
    };

        
// tabulate the HOMLFY skein modules in latex form
void latex_tabulate_HOMFLY_table() {
    Cknot *K;
    
    ostringstream ss;
    string s, *ps;
    vector<string> vs;
    
    generate_knot_names_bar(); // 0_1, 1_1, 1_2,...
    LATEX = YES; //
    CAPITALS_FIRST = false; //
    
    cout << endl << endl;
    cout << "\\chapter{The HOMFLY skein modules}" << endl << endl;

    cout << "\\section*{Knots in the solid torus}\\markright{Knots in the solid torus}" << endl << endl;
    
    vs.clear();
    
    // print HSMs of the knots in the solid torus
    for (int k=0;k<n_prime;k++) {
        
        // print output to string stream ss
        // empty string
        ss.str(""); ss.clear();
        
        ss <<  "{\\small $\\hsm(" <<knot_name[k] <<  ") = ";
        K = (k == 0 ? new Cknot("Knot 0: 1 -1 + 0 *1 &01") : get_knot_by_id(knot_name_ID[k]));
        if (winding_number(K) < 0) K->reverse();
        ss << HOMFLY(K);
        ss << "$}\\\\" <<  endl;
        s = ss.str();
        ps = new string(s);
        // place the string into the vector vs
        vs.push_back(*ps);
    }
    
    // sort strings by knot name
    sort (vs.begin(), vs.end(), less_than_name());
    // print HSM
    for (vector<string>::iterator is = vs.begin(); is != vs.end(); is++) cout << (*is);
    
    //generate_knot_names(false); // 0_1, 1_1, 1_2,...
    cout << endl << endl;
    
    Cmultivariate * pl, *pt;
    
    cout << "\\newpage" << endl;
    


    for (int l=0;l<N_HOMEO;l++) if (lens_homeo[l][1] <= 1){ // loop through all lens spaces

        cout << endl << endl;
        cout << "\\section*{Knots in $" << lpq_s[l] <<"$}\\markright{Knots in $" << lpq_s[l] <<"$}" << endl << endl;
        
        vs.clear();

        for (int k=0;k<n_prime;k++) { // loop through knots
        
            if (classification_lpq[l][knot_name_ID[k]].duplicate) continue;
        
            K = (k == 0) ? new Cknot("Knot 0: 1 -1 + 0 *1 &01") : K = get_knot_by_id(knot_name_ID[k]);
        
            if (winding_number(K) < 0) K->reverse();
            pt = HOMFLY(K); pt->simplify();
            pl = HSM_lpq(l,HOMFLY(K)); pl->simplify();
        
            if (*pt == *pl) continue; // do not print knots if the HOMLFY is equal to that of the solid torus

            // print output to string stream ss
            // empty string
            ss.str(""); ss.clear();
            
            ss <<  "{\\small $\\hsm(" <<knot_name[k] <<  ") = ";
    
            if (writhe(K) < 0) K->reverse(); // orient the knot so the writhe is non-negative (so the leading term is a "positive" generator)
        
            ss << HSM_lpq(l,HOMFLY(K));
            ss << "$}\\\\" <<  endl;
            
            s = ss.str();
            ps = new string(s);
    
            // place the string into the vector vs
            vs.push_back(*ps);
        }
  
        // sort strings by knot name
        sort (vs.begin(), vs.end(), less_than_name());
        // print HSM
        for (vector<string>::iterator is = vs.begin(); is != vs.end(); is++)
            cout << (*is);
    }
    return;
}
        
// tabulate the KBSMs in latex form
void latex_tabulate_KBSM_table() {
    Cknot *K;
    
    ostringstream ss;
    string s;
    vector<string> vs;
    
    LATEX = YES; // ?
    CAPITALS_FIRST = NO; // ?
            
    generate_knot_names_bar();
    
    cout << endl << endl;
    cout << "\\chapter{The Kauffman bracket skein modules}" << endl << endl;
            
    cout << "\\section*{Knots in the solid torus}\\markright{Knots in the solid torus}" << endl << endl;
    
    vs.clear();
    // knots in the torus
    for (int k=0;k<n_prime;k++) {
        
        // print output to string stream ss
        // empty string
        ss.str(""); ss.clear();
        
        ss <<  "{\\small $\\kbsm(" <<knot_name[k] <<  ") = ";
        K = (k == 0) ? new Cknot("Knot 0: 1 -1 + 0 *1 &01") : get_knot_by_id(knot_name_ID[k]);
        ss << KBSM(K);
        ss << "$}\\\\" <<  endl;
        s = ss.str();
        
        // place the string into the vector vs
        vs.push_back(s);
    }
    
    // sort strings by knot name
    sort (vs.begin(), vs.end(), less_than_name());
    // print HSM
    for (vector<string>::iterator is = vs.begin(); is != vs.end(); is++) cout << (*is);
    
    cout << endl << endl;
            
    Cmultivariate * pl, *pl_, *pt, *pt_;
    
    for (int l=0;l<N_HOMEO;l++) { // loop through all lens spaces
        cout << endl << endl;
        cout << "\\section*{Knots in $" << lpq_s[l] <<"$}\\markright{Knots in $" << lpq_s[l] <<"$}" << endl << endl;
        
        vs.clear();
        for (int k=0;k<n_prime;k++) {
                    
            if (classification_lpq[l][knot_name_ID[k]].duplicate) continue;
            
            K = (k == 0) ? new Cknot("Knot 0: 1 -1 + 0 *1 &01") : get_knot_by_id(knot_name_ID[k]);
            if (writhe(K) < 0) K->reverse();  // not needed in KBSM
                    
            pt = KBSM(K);
            pl = KBSM_lpq(l, KBSM(K));
                    
            pt_ = pt->deepCopy();
            pl_ = pl->deepCopy();
                    
            normalize_framing(pt_);
            normalize_framing(pl_);
                    
            if (*pt_ == *pl_) continue;
            
            ss.str(""); ss.clear();
            
            ss <<  "{\\small $\\kbsm(" <<knot_name[k] <<  ") = ";
            ss << KBSM_lpq(l,KBSM(K));
            ss << "$}\\\\" <<  endl;
            s = ss.str();
            
            // place the string into the vector vs
            vs.push_back(s);
        }
        // sort strings by knot name
        sort (vs.begin(), vs.end(), less_than_name());
        // print HSM
        for (vector<string>::iterator is = vs.begin(); is != vs.end(); is++) cout << (*is);
    }
            
    return;
}
        
// print groups of knots that share HSM and KBSM
void print_knot_groups() {
    cout << endl <<"Unclassified torus HSM:  ";
    print_non_duplicates(G_HOMFLY, false, classification);
    cout << endl;
    cout << endl<< "Unclassified torus KBSM: ";
    print_non_duplicates(G_KBSM, false, classification);
            
            
    for (int l=0; l<N_HOMEO; l++) {
                
        if (lens_homeo[l][1] <= 2) {
                    
            cout << endl<<endl<< "All          " << lpq_s[l] <<" HSM:  ";
            print_non_duplicates(G_HOMFLY_LPQ[l], true, classification_lpq[l]);
                    
            cout << endl<<endl<< "Unclassified " << lpq_s[l] <<" HSM:  ";
            print_non_duplicates(G_HOMFLY_LPQ[l], false, classification_lpq[l]);
        }
                
        cout << endl<<endl<< "All          " << lpq_s[l] <<" KBSM: ";
        print_non_duplicates(G_KBSM_LPQ[l], true, classification_lpq[l]);
                
        cout << endl<<endl<< "Unclassified " << lpq_s[l] <<" KBSM: ";
        print_non_duplicates(G_KBSM_LPQ[l], false, classification_lpq[l]);
                
    }
            
}
        
#endif
