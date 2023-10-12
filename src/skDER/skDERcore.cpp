/*
AUTHOR - Rauf Salamzade
DATE - 07/23/2023
PROGRAM NAME - skDERcore.cpp
DESCRIPTION - Core program of skDER to dereplicate genomes.
*/

#include <sys/resource.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <math.h>
using namespace std;

string delim = "\t";

/*
Taken from Arafat Hasan's response on StackOverflow:
https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
*/
// for string delimiter
vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

int main (int argc, char* argv[]) {
    if ( argv[1]==NULL || (argv[1][0]=='-' && argv[1][1]=='h') || (argv[1][0]=='-' && argv[1][1]=='-' && argv[1][2]=='h') ) {
	    cout << "Usage:" << endl;
	    cout << "skDERcore <skani triangle edge listing file> <n50 listing file> <min ANI> <min AF> <max AF difference for disqualification>" << endl;
	    return 0;
    }
    else {
         /*
         Parse cutoff as float
         */

        double min_ani;
        sscanf(argv[3],"%lf",&min_ani);
        double min_af;
        sscanf(argv[4],"%lf",&min_af);
        double max_af_difference;
        sscanf(argv[5],"%lf",&max_af_difference);
        
        /*
        First pass through skani edge file to assess connectivity of each genome.
        */
        string line, query, subject;
        int line_count_1 = 0;
        double ani, af_query, af_subject, delta_af;
        int split_counter;
        map<string, int> connectivity_dictionary;
        vector<string> v;
        ifstream input_file;
        input_file.open (argv[1]);
        if (input_file.is_open()) {
            while (input_file.good()) {
                getline (input_file,line);
                if (!line.empty() and line_count_1 > 0) {
                    split_counter = 0;
                    v = split (line, delim);
                    for (auto i : v) {
                        if (split_counter == 0) {
                            query = i;
                        }
                        else if (split_counter == 1) {
                            subject = i;
                        }
                        else if (split_counter == 2) {
                            ani = stod(i);
                        }
                        else if (split_counter == 3) {
                            af_query = stod(i);
                        }
                        else if (split_counter == 4) {
                            af_subject = stod(i);
                        }
                        split_counter++;
                    }
                    if ((ani >= min_ani) && (af_query >= min_af || af_subject >= min_af)) {
                        connectivity_dictionary[query]++;
                        connectivity_dictionary[subject]++;
                    }
                }
                line_count_1++;
            }
        } else {
  	        cout << "ERROR: Unable to open file " + (string)argv[1] << endl;
        }
        input_file.close();

        /*
        Parse N50 statistics per assembly.
        */
        map<string, int> n50_dictionary;
        int n50_val;
	    string sample;
        input_file.open (argv[2]);
        if (input_file.is_open()) {
            while (input_file.good()) {
                getline (input_file,line);
                if (!line.empty()) {
                    split_counter = 0;
                    v = split (line, delim);
                    for (auto i : v) {
                        if (split_counter == 0) {
                            sample = i;
                        }
                        else if (split_counter == 1) {
                            n50_val = stoi(i);
                        }
                        split_counter++;
                    }
                    n50_dictionary[sample] = n50_val;
                }
            }
        } else {
  	        cout << "ERROR: Unable to open file " + (string)argv[1] << endl;
        }
        input_file.close();

        /*
        Second pass to actually perform dereplication.
        */
        input_file.open (argv[1]);
        set<string> redundancy_set;
        int line_count_2 = 0;
        map<string, double> reference_score;
        double query_score, subject_score;
        if (input_file.is_open()) {
            while (input_file.good()) {
                getline (input_file,line);
                if (!line.empty() && line_count_2 > 0) {
                    split_counter = 0;
                    v = split (line, delim);
                    for (auto i : v) {
                        if (split_counter == 0) {
                            query = i;
                        }
                        else if (split_counter == 1) {
                            subject = i;
                        }
                        else if (split_counter == 2) {
                            ani = stod(i);
                        }
                        else if (split_counter == 3) {
                            af_query = stod(i);
                        }
                        else if (split_counter == 4) {
                            af_subject = stod(i);
                        }
                        split_counter++;
                    }
                    if ((ani >= min_ani) && (af_query >= min_af || af_subject >=  min_af)) {
                        delta_af = af_query - af_subject;
                        if (delta_af <= max_af_difference) {
                            if (af_query > af_subject) {
                                redundancy_set.insert(query);
                            } else {
                                redundancy_set.insert(subject);
                            }
                        } else {
                            query_score = (double)n50_dictionary[query]*(double)connectivity_dictionary[query];
                            subject_score = (double)n50_dictionary[subject]*(double)connectivity_dictionary[subject];
                            if (query_score >= subject_score) {
                                redundancy_set.insert(subject);
                            } else {
                                redundancy_set.insert(query);
                            }
                        }
                    }
                }
                line_count_2++;
            }
        } else {
  	        cout << "ERROR: Unable to open file " + (string)argv[1] << endl;
        }
        input_file.close();

        /*
        Parse N50 file to see if they are redundancy set.
        */
        input_file.open (argv[2]);
        if (input_file.is_open()) {
            while (input_file.good()) {
                getline (input_file,line);
                if (!line.empty()) {
                    split_counter = 0;
                    v = split (line, delim);
                    for (auto i : v) {
                        if (split_counter == 0) {
                            sample = i;
                        }
                    	split_counter++;
		    }
                    if (!redundancy_set.count(sample)) {
                        cout << sample << endl;
                    }
		  
                }
            }
        } else {
  	        cout << "ERROR: Unable to open file " + (string)argv[1] << endl;
        }
        input_file.close();

        return 0;
    }
}
