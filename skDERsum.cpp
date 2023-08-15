/*
AUTHOR - Rauf Salamzade
DATE - 07/23/2023
PROGRAM NAME - skDERsum.cpp
DESCRIPTION - Program in skDER to get information for greedy clustering.
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


string joinVectorString (vector<string> lst) {
    string ret;
    for(const auto &s : lst) {
        if(!ret.empty())
            ret += "; ";
        ret += s;
    }
    return ret;
}

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

string doubleToString(double val) {
	string result;
	ostringstream convertd;
	convertd << val;
	result = convertd.str();
	return result;
}

int main (int argc, char* argv[]) {
    if ( argv[1]==NULL || (argv[1][0]=='-' && argv[1][1]=='h') || (argv[1][0]=='-' && argv[1][1]=='-' && argv[1][2]=='h') ) {
	    cout << "Usage:" << endl;
	    cout << "skDERsum <skani triangle edge listing file> <n50 listing file> <alignment fraction cutoff>" << endl;
	    return 0;
    }
    else {
	double min_af;
	sscanf(argv[3],"%lf",&min_af);

        /*
        First pass through skani edge file to assess connectivity of each genome.
        */
        string line, query, subject;
        int line_count_1 = 0;
        double af_query, af_subject;
        int split_counter;
	double query_af, subject_af; 
        set<string> all_samples;
        map<string, int> connectivity_dictionary;
        map<string, vector<string>> connected_members;
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
			else if (split_counter == 3) {
			    query_af = stod(i);
			} 
			else if (split_counter == 4) {
			    subject_af = stod(i);
			}
                        split_counter++;
                    }
		    // this is the converse of what I had originally 
		    // because the larger genomes were getting disregarded
		    if (subject_af >= min_af) {
			all_samples.insert(query);
                        connectivity_dictionary[query]++;
			connected_members[query].push_back(subject);
		    }
		    if (query_af >= min_af) {
			all_samples.insert(subject);
                        connectivity_dictionary[subject]++;
		        connected_members[subject].push_back(query);
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
        int n50_val;
        double score;
	string sample, member_list;
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
		    if (all_samples.count(sample) != 0) {
                        score = (double)n50_val*(double)connectivity_dictionary[sample];
			member_list = joinVectorString(connected_members[sample]);
                        cout << sample + '\t' + doubleToString(score) + '\t' + member_list << endl;
                    } else {
			cout << sample + "\t0.0\t" << endl;
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
