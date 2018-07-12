//This program runs two jet clustering algorithms on the same data,
//ConstituentSubtractor and IteratedConstituentSubtractor.
//Algorithm parameters are read off of FILL IN LATER, and the
//deviances of the jet parameters are written to FILL IN LATER.

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "ConstituentSubtractor.hh"
#include "IterativeConstituentSubtractor.hh"

#include "functions.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;
using namespace fastjet;

vector< vector<double> > ConstituentSubtract(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, double max_distance, double alpha, double ghost_area, double grid_spacing, int common_bge, int cutoff);
vector< vector<double> > IteratedConstituentSubtract(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, vector<double> max_distances, vector<double> alphas, double ghost_area, double grid_spacing, int common_bge, int cutoff);
vector<string> get_pieces(string line);
vector< vector<int> > match_jets(vector <vector<double> > hard_coors, vector< vector<double> > sub_coors);
double R_dist(vector<double> point1, vector<double> point2);
double phi_dist(double phi1, double phi2);

int main() {

	//Load in all the algorithm parameters from a file
	ifstream paramfile;
	paramfile.open("subtractor_parameters.txt");
	string data;
	string data2;

	//First, set up eta cutoffs
	getline(paramfile, data);
	double max_eta = atof(get_pieces(data)[2].c_str());
	getline(paramfile, data);
	double max_eta_jet = atof(get_pieces(data)[2].c_str());
	getline(paramfile, data);

	//Determine the jet algorithm
	int algorithm number = atoi(get_pieces(data)[2].c_str());
	getline(paramfile, data2);
	double R = atof(get_pieces(data2)[2].c_str());
	JetDefinition jet_def;
	if (algorithm_number == 0) {
		jet_def = JetDefinition(kt_algorithm, R);
	}
	if (algorithm_number == 1) {
		jet_def = JetDefinition(cambridge_algorithm, R);
	}
	if (algorithm_number == 2) {
		jet_def = JetDefinition(antikt_algorithm, R);
	}

	//Read in an extra line to discard the distance type.
	//I may just leave that permanently on deltaR since I
	//don't know what the other options are.
	getline(paramfile, data);

	//Get the max_distances and alphas for the CS and ICS algorithms
	getline(paramfile, data);
	double max_distance_CS = atof(get_pieces(data)[2].c_str());
	getline(paramfile, data);
	vector<double> max_distances_ICS;
	vector<string> temp = get_pieces(data);
	for (int i = 2; i < temp.size(); i++) {
		max_distances_ICS.push_back(atof(temp[i].c_str()));
	}

	getline(paramfile, data);
	double alpha_CS = atof(get_pieces(data)[2].c_str());
	getline(paramfile, data);
	vector<double> alpha_ICS;
	temp = get_pieces(data);
	for (int i = 2; i < temp.size(); i++) {
		alphas_ICS.push_back(atof(temp[i].c_str()));
	}

	//Now get stuff that has to do with background estimation
	getline(paramfile, data);
	double ghost_area = atof(get_pieces(data)[2].c_str());
	getline(paramfile, data);
	double bge_grid_spacing = atof(get_pieces(data)[2].c_str());
	int common_bge = 0;
	getline(paramfile, data);
	if (get_pieces(data)[2] == "Y") {
		common_bge = 1;
	}

	paramfile.close();

	//read in the input particles. THIS WILL BE MODIFIED FOR ROOT LATER
	vector<PseudoJet> hard_event, full_event;
	read_event(hard_event, full_event);

	hard_event = SelectorAbsEtaMax(max_eta)(hard_event);
	full_event = SelectorAbsEtaMax(max_eta)(full_event);

	//Make a jet selector and a cutoff of 20 GeV. MODIFY BASED ON DATA
	int cutoff = 20;

	//Run the clustering on the hard jets and the full jets
	ClusterSequence clust_seq_hard(hard_event, jet_def);
	ClusterSequence clust_seq_full(full_event, jet_def);
	vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets(cutoff));
	vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets(cutoff));

	//Now pull in the subtracted jets
	vector< vector<double> > corrected_CS_jet_params = ConstituentSubtract(full_event, max_eta, sel_jets, jet_def, max_distance_CS, alpha_CS, ghost_area, grid_spacing, common_bge, cutoff);
	vector< vector<double> > corrected_ICS_jet_params = IteratedConstituentSubtract(full_event, max_eta, sel_jets, jet_def, max_distances_ICS, alphas_ICS, ghost_area, grid_spacing, common_bge, cutoff);

	//Make an "answer key" of which hard jets match to which subtracted jets
	//First, get the (eta,phi) coordinates of all the jets
	vector< vector<double> > hard_coors;
	for (int i = 0; i < hard_jets.size(); i++) {
		vector<double> coors;
		coors.push_back(hard_jets[i].eta());
		coors.push_back(hard_jets[i].phi());
		hard_jet_coors.push_back(coors);
	}

	vector< vector<double> > CS_sub_coors;
	for (int i = 0; i < corrected_CS_jet_params.size(); i++) {
		vector<double> coors;
		coors.push_back(corrected_CS_jet_params[i][1]);
		coors.push_back(corrected_CS_jet_params[i][4]);
		CS_sub_coors.push_back(coors);
	}

	vector< vector<double> > ICS_sub_coors;
	for (int i = 0; i < corrected_ICS_jet_params.size(); i++) {
		vector<double> coors;
		coors.push_back(corrected_ICS_jet_params[i][1]);
		coors.push_back(corrected_ICS_jet_params[i][4]);
		ICS_sub_coors.push_back(coors);
	}

	//Now, make the matching key
	vector< vector<int> > CS_key;
	CS_key = match_jets(hard_coors, CS_sub_coors);
	vector< vector<int> > ICS_key;
	ICS_key = match_jets(hard_coors, ICS_sub_coors);

	JetWidth width;

	//Get the deviation on all the main jet parameters
	ofstream results;
	results.open("jet_differences.txt");
	results << hard_jets.size() << endl;
	for (int i = 0; i < hard_jets.size(); i++) {
		PseudoJet jet = hard_jets[i];
		int CS_tag = CS_key[i][1];
		vector<double> CSjet = corrected_CS_jet_params[CS_tag];
		results << abs(jet.pt() - CSjet[0]) / jet.pt() << " " << abs(jet.rap() - CSjet[1]) / jet.rap() << " " << abs(jet.m() - CSjet[2]) / jet.m() << " " << abs(width(jet) - CSjet[3]) / width(jet) << endl;
	}
	for (int i = 0; i < hard_jets.size(); i++) {
		PseudoJet jet = hard_jets[i];
		int ICS_tag = ICS_key[i][1];
		vector<double> ICSjet = corrected_ICS_jet_params[ICS_tag];
		results << abs(jet.pt() - ICSjet[0]) / jetpt() << " " << abs(jet.rap() - ICSjet[1]) / jet.rap() << " " << abs(jet.m() - ICSjet[2]) / jet.m() << " " << abs(width(jet) - ICSjet[3]) / width(jet) << endl;
	}

	return 0;
}

//this function does the CS algorithm and returns the relevant jet parameters
vector< vector<double> > ConstituentSubtract(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, double max_distance, double alpha, double ghost_area, double grid_spacing, int common_bge, int cutoff) {
	contrib::ConstituentSubtractor subtractor;
	subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR);
	subtractor.set_max_distance(max_distance);
	subtractor.set_alpha(alpha);
	subtractor.set_ghost_area(ghost_area);
	subtractor.set_max_eta(max_eta);
	subtractor.initialize();

	//Set up the background estimator
	GridMedianBackgroundEstimator bge_rho(max_eta, grid_spacing);
	bge_rho.set_particles(full_event);
	subtractor.set_background_estimator(&bge_rho);
	if (common_bge == 1) {
		subtractor.set_common_bge_for_rho_and_rhom(true);
	}
	if (common_bge == 0) {
		subtractor.set_common_bge_for_rho_and_rhom(false);
	}

	//Run the subtractor and cluster the new event
	vector<PseudoJet> corrected_event = subtractor.subtract_event(full_event);
	ClusterSequence clust_seq_corr(corrected_event, jet_def);
	vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets(cutoff));

	//Load in the jet parameters from the corrected jet
	JetWidth width;

	vector< vector<double> > jet_params;
	for (unsigned int i = 0; i < corrected_jets.size(); i++) {
		const PseudoJet &jet = corrected_jets[i];
		vector<double> temp_params;
		temp_params.push_back(jet.pt());
		temp_params.push_back(jet.rap());
		temp_params.push_back(jet.m());
		temp_params.push_back(width(jet));
		temp_params.push_back(jet.phi());
		jet_params.push_back(temp_params);
	}

	return jet_params;
}

//this function does the ICS algorithm and returns the relevant jet parameters
vector< vector<double> > IteratedConstituentSubtract(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, vector<double> max_distances, vector<double> alphas, double ghost_area, double grid_spacing, int common_bge, int cutoff) {
	contrib::IteratedConstituentSubtractor subtractor;
	subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR);
	subtractor.set_parameters(max_distances, alphas);
	subtractor.set_remove_remaining_proxies(true);
	subtractor.set_ghost_area(ghost_area);
	subtractor.set_max_eta(max_eta);
	subtractor.initialize();

	//Set up the background estimator
	GridMedianBackgroundEstimator bge_rho(max_eta, grid_spacing);
	bge_rho.set_particles(full_event);
	subtractor.set_background_estimator(&bge_rho);
	if (common_bge == 1) {
		subtractor.set_common_bge_for_rho_and_rhom(true);
	}
	if (common_bge == 0) {
		subtractor.set_common_bge_for_rho_and_rhom(false);
	}

	//Run the subtractor and cluster the new event
	vector<PseudoJet> corrected_event = subtractor.subtract_event(full_event);
	ClusterSequence clust_seq_corr(corrected_event, jet_def);
	vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets(cutoff));

	//Load in the jet parameters from the corrected jet
	JetWidth width;

	vector< vector<double> > jet_params;
	for (unsigned int i = 0; i < corrected_jets.size(); i++) {
		const PseudoJet &jet = corrected_jets[i];
		vector<double> temp_params;
		temp_params.push_back(jet.pt());
		temp_params.push_back(jet.rap());
		temp_params.push_back(jet.m());
		temp_params.push_back(width(jet));
		temp_params.push_back(jet.phi());
		jet_params.push_back(temp_params);
	}

	return jet_params;
}

//This function takes a string and returns a vector of substrings, using space as a delimiter
vector<string> get_pieces(string line) {

	int pos = -1;
	int backmark = 0;
	int length = 0;
	string temp;
	vector<string> tokens;
	int pos1;
	tokens.push_back(line.substr(0, line.find(' ')));
	while (line.find(' ', pos + 1) != string::npos) {
		pos = line.find(' ', backmark);
		length = line.find(' ', pos + 1) - pos;
		temp = line.substr(pos + 1, length);
		backmark = pos + 1;
		tokens.push_back(temp);
	}

	return tokens;
}

vector< vector<int> > match_jets(vector< vector<double> > hard_coors, vector< vector<double> > sub_coors) {
	vector< vector<int> > key;
	for (int i = 0; i < hard_coors.size(); i++) {
		vector<double> jet = hard_coors[i];
		int closest_sub = 0;
		double min_dist = 1000;
		for (int j = 0; i < sub_coors.size(); j++) {
			if (R_dist(jet, sub_coors[j]) < min_dist) {
				min_dist = R_dist(jet, sub_coors[j]);
				closest_sub = j;
			}
		}
		vector<int> temp;
		temp.push_back(i);
		temp.push_back(closest_sub);
		key.push_back(temp);
	}
	return key;
}

double R_dist(vector<double> point1, vector<double> point2) {
	double dist = pow(pow(point1[0] - point2[0], 2) + pow(phi_dist(point1[1], point2[1]), 2), 0.5);

	return dist;
}

double phi_dist(double phi1, double phi2) {
	double dist;
	double temp_dist = abs(phi1 - phi2);
	double opt_1 = (atan(1) * 8 - phi2) + phi1;
	double opt_2 = (atan(1) * 8 - phi1) + phi2;
	if (temp_dist <= opt_1 && temp_dist <= opt_2) {
		dist = temp_dist;
	}
	else if (opt_1 < opt_2) {
		dist = opt_1;
	}
	else {
		dist = opt_2;
	}

	return dist;
}