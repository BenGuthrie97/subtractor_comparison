//this program runs two jet clustering algorithms on the same data,
//identical except that one uses ConstituentSubtractor, and the other
//uses IteratedConstituentSubtractor. The results of these two different
//algorithms are then written to a text file... will be implemeneted later.

//run this with ./subtractor-comparison < ../data/Pythia and just tab-complete to get the rest.

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "ConstituentSubtractor.hh"
#include "IterativeConstituentSubtractor.hh"
#include "SoftKiller/SoftKiller.hh"

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

vector< vector<double> > CS_SK(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, double max_distance, double alpha, double ghost_area, double grid_spacing, int common_bge, int cutoff, double SK_grid_size);

vector< vector<double> > IteratedConstituentSubtract(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, vector<double> max_distances, vector<double> alphas, double ghost_area, double grid_spacing, int common_bge, int cutoff, bool proxies);

vector< vector<double> > ICS_SK(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, vector<double> max_distances, vector<double> alphas, double ghost_area, double grid_spacing, int common_bge, int cutoff, double SK_grid_spacing, bool proxies);

vector< vector<double> > get_coors(vector< vector<double> > params);

vector<int> adjust_excess(vector< vector<int> > key, vector<int> excesses);

vector<string> write_results(vector< vector<int> > key, vector< vector<double> > sigs, vector< vector<double> > params, int type);

int findNHardest(vector<PseudoJet> full_event, JetDefinition jet_def, double max_eta_jet);
vector<string> get_pieces(string line);
vector< vector<int> > match_jets(vector< vector<double> > hard_coors, vector< vector<double> > sub_coors);
vector< vector<int> > match_all_jets(vector< vector<double> > signal_coors, vector< vector<double> > CS_coors, vector< vector<double> > ICS_coors);
vector< vector<int> > match_some_jets(vector< vector<double> > signal_coors, vector< vector<double> > full_coors);
void print_coors(vector< vector<double> > coors);
void print_key(vector< vector<int> > key);
double R_dist(vector<double> point1, vector<double> point2);
double phi_dist(double phi1, double phi2);
double phi_dist_alt(double phi1, double phi2);

int main(){

  cout << "Running." << endl;

  //Load in all the parameters from a file
  ifstream paramfile;
  paramfile.open("subtractor_parameters.txt");
  string data;
  string data2;

  //First, set up eta cutoffs
  getline(paramfile,data);
  double max_eta = atof(get_pieces(data)[2].c_str());
  getline(paramfile,data);
  double max_eta_jet = atof(get_pieces(data)[2].c_str());
  getline(paramfile,data);

  //Determine the jet algorithm
  int algorithm_number = atoi(get_pieces(data)[2].c_str());
  getline(paramfile,data2);
  double R = atof(get_pieces(data2)[2].c_str());
  JetDefinition jet_def;
  if (algorithm_number == 0){
    jet_def = JetDefinition(kt_algorithm,R);
  }
  if (algorithm_number == 1){
    jet_def = JetDefinition(cambridge_algorithm,R); 
  }
  if (algorithm_number == 2){
    jet_def = JetDefinition(antikt_algorithm,R);
  }

  //Read in an extra line to discard the distance type. I may just leave that permanently on deltaR, since I don't even know what the other options are
  getline(paramfile,data);

  //Get the max_distances and alphas for the CS and ICS algorithms
  getline(paramfile,data);
  double max_distance_CS = atof(get_pieces(data)[2].c_str());
  getline(paramfile,data);
  vector<double> max_distances_ICS;
  vector<string> temp = get_pieces(data);
  for (int i = 2; i < temp.size(); i++){
    max_distances_ICS.push_back(atof(temp[i].c_str()));
  }
  
  getline(paramfile,data);
  double alpha_CS = atof(get_pieces(data)[2].c_str());
  getline(paramfile,data);
  vector<double> alphas_ICS;
  temp = get_pieces(data);
  for (int i = 2; i < temp.size(); i++){
    alphas_ICS.push_back(atof(temp[i].c_str()));
  }

  //Get the ICS parameter about removing proxies
  bool proxies;
  getline(paramfile,data);
  if (get_pieces(data)[2] == "0"){
    proxies = false;
  }
  if (get_pieces(data)[2] == "1"){
    proxies = true;
  }

  //Get the SK parameters here.
  getline(paramfile,data);
  double SK_grid_spacing = atof(get_pieces(data)[2].c_str());

  //Now get stuff that has to do with background estimation
  getline(paramfile,data);
  double ghost_area = atof(get_pieces(data)[2].c_str());
  getline(paramfile,data);
  double bge_grid_spacing = atof(get_pieces(data)[2].c_str());
  int common_bge = 0;
  getline(paramfile,data);
  if (get_pieces(data)[2] == "Y"){
    common_bge = 1;
  }

  paramfile.close();

  //Go ahead and make the selectors before I start working on the root data.
  int cutoff = 20;
  Selector sel_jets = SelectorAbsEtaMax(max_eta_jet);
  //Pretty sure that cutoff is in units of GeV. Check on that for ttbar data

  //also, this has to go somewhere 
  JetWidth width;

  //read in input particles
  vector<PseudoJet> hard_event, full_event;
  read_event(hard_event, full_event);
  
  //Read in all the ttbar root events in their own vector of PseudoJets
  vector< vector<PseudoJet> > root_events;
  ifstream root_data;
  root_data.open("ttbar_event_output.txt"); 
  double pt, y, phi, m;
  string cluster;
  vector<string> cluster_data;
  int ncl;

  ofstream ratios;
  ratios.open("jet_ratios.txt");

  ofstream checks;
  checks.open("sanity_checks.txt");


  //There seem to be an unusual number of jets with zero mass. Check these out
  ofstream mzero;
  mzero.open("mzero_jets.txt");

  cout << "Loaded all the parameters." << endl;

  //Get cluster data
  while (getline(root_data,cluster)){
    cluster_data = get_pieces(cluster);
    //If there's only a single number, that's the counter indicating all the next
    //however-many lines are the same event.
    if (cluster_data.size() == 1){
      ncl = atoi(cluster_data[0].c_str());
      vector<PseudoJet> root_event;
      for (int i = 0; i < ncl; i++){
        getline(root_data,cluster);
        cluster_data = get_pieces(cluster);
        pt = atof(cluster_data[0].c_str());
        y = atof(cluster_data[1].c_str());
        phi = atof(cluster_data[2].c_str());
        m = atof(cluster_data[3].c_str());

        PseudoJet temp_jet = PseudoJet();
	temp_jet.reset_PtYPhiM(pt,y,phi,m);
        root_event.push_back(temp_jet);

      }
      root_events.push_back(root_event);
    } else{
      cout << "Crap, the line-reading from ttbar_cluster_output screwed up." << endl;
    }
  }

  root_data.close();

  //Do this again, for signal jets
  vector< vector< vector<double> > > all_signal_jets;
  ifstream sj_data;
  sj_data.open("ttbar_signal_jet_output.txt");
  double j_E, j_pt, j_phi, j_y, j_m;
  string s_jet;
  vector<string> s_jet_data;
  int n_sjs;

  while (getline(sj_data,s_jet)){
    s_jet_data = get_pieces(s_jet);
    if (s_jet_data.size() == 1){
      n_sjs = atoi(s_jet_data[0].c_str());
      vector< vector<double> > jets_in_event;
      for (int i = 0; i < n_sjs; i++){
        getline(sj_data,s_jet);
        s_jet_data = get_pieces(s_jet);
        j_pt = atof(s_jet_data[0].c_str());
        j_y = atof(s_jet_data[1].c_str());
        j_phi = atof(s_jet_data[2].c_str());
        j_E = atof(s_jet_data[3].c_str());
	j_m = atof(s_jet_data[4].c_str());

        vector<double> temp_s_jet;
        temp_s_jet.push_back(j_pt);
        temp_s_jet.push_back(j_y);
        temp_s_jet.push_back(j_phi);
        temp_s_jet.push_back(j_E);
	temp_s_jet.push_back(j_m);

        jets_in_event.push_back(temp_s_jet);
      }
      all_signal_jets.push_back(jets_in_event);
    } else{
      cout << "Crap, the line-reading from ttbar_signal_jet_output screwed up." << endl;
    }
  }

  sj_data.close();

  //Now, if everything is correct, the events should be fully loaded
  cout << "Number of events " << root_events.size() << endl;
  ofstream results;
  results.open("jet_differences.txt");
  
  for (int j = 0; j < root_events.size(); j++){

    //since this takes a freaking long time, print out progress reports
    if (j % 100 == 0){
      cout << "Checking event number " << j << endl;
    }

    vector<PseudoJet> root_event = root_events[j];
    root_event = SelectorAbsEtaMax(max_eta)(root_event);

    ClusterSequence clust_seq_root(root_event, jet_def);
    vector<PseudoJet> root_jets = sel_jets(clust_seq_root.inclusive_jets(cutoff));

    //Now run the subtraction algorithms.
    vector < vector<double> > corrected_CS_root_jet_params = ConstituentSubtract(root_event, max_eta, sel_jets, jet_def, max_distance_CS, alpha_CS, ghost_area, bge_grid_spacing, common_bge, cutoff);
    vector< vector<double> > corrected_ICS_root_jet_params = IteratedConstituentSubtract(root_event, max_eta, sel_jets, jet_def, max_distances_ICS, alphas_ICS, ghost_area, bge_grid_spacing, common_bge, cutoff, proxies); 

    vector< vector<double> > corrected_CS_SK_root_jet_params = CS_SK(root_event, max_eta, sel_jets, jet_def, max_distance_CS, alpha_CS, ghost_area, bge_grid_spacing, common_bge, cutoff, SK_grid_spacing);
    vector< vector<double> > corrected_ICS_SK_root_jet_params = ICS_SK(root_event, max_eta, sel_jets, jet_def, max_distances_ICS, alphas_ICS, ghost_area, bge_grid_spacing, common_bge, cutoff, SK_grid_spacing, proxies);

    const vector< vector<double> > &jets_in_event = all_signal_jets[j];

    //Now that I have all the jet info, check out the jets with no mass (zeroth entry)
    //Write them to the mzero file
    
    //First, write the event number
    mzero << j << endl;

    mzero << "CS" << endl;
    for (int i = 0; i < corrected_CS_root_jet_params.size(); i++){
      const vector<double> &jet = corrected_CS_root_jet_params[i];
      if (jet[5] == 0){
        mzero << jet[0] << " " << jet[1] << " " << jet[4] << " " << jet[2] << " " << jet[5] << " " << jet[6] << endl;
      }
    }
    mzero << "ICS" << endl;
    for (int i = 0; i < corrected_ICS_root_jet_params.size(); i++){
      const vector<double> &jet = corrected_ICS_root_jet_params[i];
      if (jet[5] == 0){
        mzero << jet[0] << " " << jet[1] << " " << jet[4] << " " << jet[2] << " " << jet[5] << " " << jet[6] << endl;
      }
    }
    mzero << "CSSK" << endl;
    for (int i = 0; i < corrected_CS_SK_root_jet_params.size(); i++){
      const vector<double> &jet = corrected_CS_SK_root_jet_params[i];
      if (jet[5] == 0){
        mzero << jet[0] << " " << jet[1] << " " << jet[4] << " " << jet[2] << " " << jet[5] << " " << jet[6] << endl;
      }
    }
    mzero << "ICSSK" << endl;
    for (int i = 0; i < corrected_ICS_SK_root_jet_params.size(); i++){
      const vector<double> &jet = corrected_ICS_SK_root_jet_params[i];
      if (jet[5] == 0){
        mzero << jet[0] << " " << jet[1] << " " << jet[4] << " " << jet[2] << " " << jet[5] << " " << jet[6] << endl;
      }
    }
    //This just lets the program know that the event is done.
    mzero << "Pass" << endl;

/*
    cout << "There are " << root_jets.size() << " jets in the " << j << " event of the root data." << endl;
    for (unsigned int i = 0; i < root_jets.size(); i++){
      const PseudoJet &jet = root_jets[i];
      cout << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << endl;
    }

    cout << endl;

    cout << "There are " << corrected_CS_root_jet_params.size() << " jets in the " << j << " event of the root data after CS-subtraction." << endl;
    for (unsigned int i = 0; i < corrected_CS_root_jet_params.size(); i++){
      const vector<double> &params = corrected_CS_root_jet_params[i];
      cout << "pt = " << params[0] << ", rap = " << params[1] << ", E = " << params[2] << ", width = " << params[3] << ", phi = " << params[4] << endl;
    }

    cout << endl;

    cout << "There are " << corrected_CS_SK_root_jet_params.size() << " jets in the " << j << " event of the root data after CS_SK-subtraction." << endl;
    for (unsigned int i = 0; i < corrected_CS_SK_root_jet_params.size(); i++){
      const vector<double> &params = corrected_CS_SK_root_jet_params[i];
      cout << "pt = " << params[0] << ", rap = " << params[1] << ", E = " << params[2] << ", width = " << params[3] << ", phi = " << params[4] << endl;
    }

    cout << endl;


    cout << "There are " << corrected_ICS_root_jet_params.size() << " jets in the " << j << " event of the root data after ICS-subtraction." << endl;
    for (unsigned int i = 0; i < corrected_ICS_root_jet_params.size(); i++){
      const vector<double> &params = corrected_ICS_root_jet_params[i];
      cout << "pt = " << params[0] << ", rap = " << params[1] << ", E = " << params[2] << ", width = " << params[3] << ", phi = " << params[4] << endl;
    }

    cout << endl;

    cout << "There are " << corrected_ICS_SK_root_jet_params.size() << " jets in the " << j << " event of the root data after CS_SK-subtraction." << endl;
    for (unsigned int i = 0; i < corrected_ICS_SK_root_jet_params.size(); i++){
      const vector<double> &params = corrected_ICS_SK_root_jet_params[i];
      cout << "pt = " << params[0] << ", rap = " << params[1] << ", E = " << params[2] << ", width = " << params[3] << ", phi = " << params[4] << endl;
    }
    
    cout << endl;

    cout << "There are " << jets_in_event.size() << " signal jets in the event." << endl;
    for (unsigned int i = 0; i < jets_in_event.size(); i++){
      const vector<double> &params = jets_in_event[i];
      cout << "pt = " << params[0] << ", rap = " << params[1] << ", phi = " << params[2] << ", E = " << params[3] << ", m = " << params[4] << endl;
    }

    cout << endl;

    //I have to make these outside the lines when I print the results because subtracting negative integers in those lines screws up.
    int full_excess;
    full_excess = root_jets.size() - jets_in_event.size();
    int CS_excess;
    CS_excess = corrected_CS_root_jet_params.size() - jets_in_event.size();
    int ICS_excess;
    ICS_excess = corrected_ICS_root_jet_params.size() - jets_in_event.size();
    int CS_SK_excess; 
    CS_SK_excess = corrected_CS_SK_root_jet_params.size() - jets_in_event.size();
    int ICS_SK_excess;
    ICS_SK_excess = corrected_ICS_SK_root_jet_params.size() - jets_in_event.size();
    //I believe I can now handle the cases when the algorithms find too few jets. I'll just have to throw out the zero-jet cases
    //For now, I'll quit if either has too few, but I ought to make a case for when only one or the other has too few.
    if (corrected_CS_root_jet_params.size() == 0 || corrected_ICS_root_jet_params.size() == 0 || corrected_CS_SK_root_jet_params.size() == 0 || corrected_ICS_SK_root_jet_params.size() == 0){
      //Write -1 to tell the program just to read the jet excesses, then stop
      results << -1 << endl;
      results << full_excess << " " << CS_excess << " " << ICS_excess << " " << CS_SK_excess << " " << ICS_SK_excess << endl;
      continue;
    }

    //Get all the coordinates so I can do delta-R matching.
    vector< vector<double> > signal_coors;
    for (int i = 0; i < jets_in_event.size(); i++){
      vector<double> temp_coors;
      temp_coors.push_back(jets_in_event[i][1]);
      temp_coors.push_back(jets_in_event[i][2]);
      signal_coors.push_back(temp_coors);
    }

    vector< vector<double> > full_coors;
    for (int i = 0; i < root_jets.size(); i++){
      vector<double> temp_coors;
      temp_coors.push_back(root_jets[i].rap());
      temp_coors.push_back(root_jets[i].phi_std());
      full_coors.push_back(temp_coors);
    }
 
    vector< vector<double> > CS_coors = get_coors(corrected_CS_root_jet_params);
    vector< vector<double> > ICS_coors = get_coors(corrected_ICS_root_jet_params);
    vector< vector<double> > CS_SK_coors = get_coors(corrected_CS_SK_root_jet_params);
    vector< vector<double> > ICS_SK_coors = get_coors(corrected_ICS_SK_root_jet_params);

    //I'll always assume that the signal/hard event has the fewest number of jets.
    //THAT WAS A BAD ASSUMPTION.
    //The key has the indices for matching jets in the order (signal, CS, ICS).
    vector< vector<int> > matching_key;
    matching_key = match_all_jets(signal_coors,CS_coors,ICS_coors);
    vector< vector<int> > matching_key_SK = match_all_jets(signal_coors,CS_SK_coors,ICS_SK_coors);
    vector< vector<int> > matching_key_full = match_some_jets(signal_coors,full_coors);

    //With the matching key, write the deviances to the file.
    //Also write the number of excessive jets before the deviances
    //Parameters are written in the following order: pt, y, phi, E.
    //Do I need to do jet width?

    //I need to figure out how many jets will be written to the file for this event.
    //Figure out many of each kind of jet will be written to the results file.
    int numFull = 0;
    int numCS = 0;
    int numICS = 0;
    int numCSSK = 0;
    int numICSSK = 0;
    for (int i = 0; i < matching_key_full.size(); i++){
      if(matching_key_full[i][1] != -1){
	numFull++;
      }
    }
    for (int i = 0; i < matching_key.size(); i++){
      if (matching_key[i][1] != -1){
	numCS++;
      }
    }
    for (int i = 0; i < matching_key.size(); i++){
      if (matching_key[i][2] != -1){
	numICS++;
      }
    }
    for (int i = 0; i < matching_key_SK.size(); i++){
      if (matching_key_SK[i][2] != -1){
	numCSSK++;
      }
    }
    for (int i = 0; i < matching_key_SK.size(); i++){
      if (matching_key_SK[i][2] != -1){
	numICSSK++;
      }
    }
    
    //Now write those as instructions for the graphing file for how many of each to read
    results << numFull << " " << numCS << " " << numICS << " " << numCSSK << " " << numICSSK << endl;

    //Use the keys to adjust excesses for when jets failed to match
    vector<int> normal_excesses;
    normal_excesses.push_back(CS_excess);
    normal_excesses.push_back(ICS_excess);
    vector<int> SK_excesses;
    SK_excesses.push_back(CS_SK_excess);
    SK_excesses.push_back(ICS_SK_excess);
    vector<int> full_excesses;
    full_excesses.push_back(full_excess);
    normal_excesses = adjust_excess(matching_key,normal_excesses);
    SK_excesses = adjust_excess(matching_key_SK,SK_excesses);
    full_excesses = adjust_excess(matching_key_full,full_excesses);
    CS_excess = normal_excesses[0];
    ICS_excess = normal_excesses[1];
    CS_SK_excess = SK_excesses[0];
    ICS_SK_excess = SK_excesses[1];

    results << full_excess << " " <<  CS_excess << " " << ICS_excess << " " << CS_SK_excess << " " << ICS_SK_excess << endl;

    //Now write the signed deviances to the results file, in the order pt, rap, phi, E, m.

    //Write the full ones first
    for (int i = 0; i < matching_key_full.size(); i++){
      int Sigtag = matching_key_full[i][0];
      const vector<double> &sigjet = jets_in_event[Sigtag];
      int fulltag = matching_key_full[i][1];
      if (fulltag != -1){
	const PseudoJet &fulljet = root_jets[fulltag];
	results << (fulljet.pt() - sigjet[0])/sigjet[0] << " " << (fulljet.rap() - sigjet[1])/sigjet[1] << " " << phi_dist_alt(sigjet[2],fulljet.phi_std())/sigjet[2] << " " << (fulljet.E() - sigjet[3])/sigjet[3] << " " << (fulljet.m() - sigjet[4])/sigjet[4] << endl;
      }
    }

    //cout << "Wrote full jet results." << endl;

    //Then CS
    for (int i = 0; i < matching_key.size(); i++){
      int Sigtag = matching_key[i][0];
      const vector<double> &sigjet = jets_in_event[Sigtag];
      int jettag = matching_key[i][1];
      if (jettag != -1){
        const vector<double> &CSjet = corrected_CS_root_jet_params[jettag];
        results << (CSjet[0] - sigjet[0])/sigjet[0] << " " << (CSjet[1] - sigjet[1])/sigjet[1] << " " << phi_dist_alt(sigjet[2], CSjet[4])/sigjet[2] << " " << (CSjet[2] - sigjet[3])/sigjet[3] << " " << (CSjet[5] - sigjet[4])/sigjet[4] << endl;
      }
    }

    //cout << "Wrote CS results." << endl;

    //Then ICS
    for (int i = 0; i < matching_key.size(); i++){
      int Sigtag = matching_key[i][0];
      const vector<double> &sigjet = jets_in_event[Sigtag];
      int jettag = matching_key[i][2];
      if (jettag != -1){
        const vector<double> &ICSjet = corrected_ICS_root_jet_params[jettag];
        results << (ICSjet[0] - sigjet[0])/sigjet[0] << " " << (ICSjet[1] - sigjet[1])/sigjet[1] << " " << phi_dist_alt(sigjet[2], ICSjet[4])/sigjet[2] << " " << (ICSjet[2] - sigjet[3])/sigjet[3] << " " << (ICSjet[5] - sigjet[4])/sigjet[4] << endl;
      }
    }

    //cout << "Wrote ICS results." << endl;

    //Then CS_SK
    for (int i = 0; i < matching_key_SK.size(); i++){
      int Sigtag = matching_key_SK[i][0];
      const vector<double> &sigjet = jets_in_event[Sigtag];
      int jettag = matching_key_SK[i][1];
      if (jettag != -1){
        const vector<double> &CSSKjet = corrected_CS_SK_root_jet_params[jettag];
        results << (CSSKjet[0] - sigjet[0])/sigjet[0] << " " << (CSSKjet[1] - sigjet[1])/sigjet[1] << " " << phi_dist_alt(sigjet[2], CSSKjet[4])/sigjet[2] << " " << (CSSKjet[2] - sigjet[3])/sigjet[3] << " " << (CSSKjet[5] - sigjet[4])/sigjet[4] << endl;
      }
    }

    //cout << "Wrote CS_SK results." << endl;

    //Then ICS_SK
    for (int i = 0; i < matching_key_SK.size(); i++){
      int Sigtag = matching_key_SK[i][0];
      const vector<double> &sigjet = jets_in_event[Sigtag];
      int jettag = matching_key_SK[i][2];
      if (jettag != -1){
        const vector<double> &ICSSKjet = corrected_ICS_SK_root_jet_params[jettag];
        results << (ICSSKjet[0] - sigjet[0])/sigjet[0] << " " << (ICSSKjet[1] - sigjet[1])/sigjet[1] << " " << phi_dist_alt(sigjet[2], ICSSKjet[4])/sigjet[2] << " " << (ICSSKjet[2] - sigjet[3])/sigjet[3] << " " << (ICSSKjet[5] - sigjet[4])/sigjet[4] << endl;
      }
    }

    //this is an entirely different part of the program, used to make output for jet_ratio_plotter.cxx

    //Now, make some keys for each kind of jet, compared to the FULL jets, not the signal!
    vector< vector<int> > key_CS = match_some_jets(full_coors, CS_coors);
    vector< vector<int> > key_ICS = match_some_jets(full_coors, ICS_coors);
    vector< vector<int> > key_CSSK = match_some_jets(full_coors, CS_SK_coors);
    vector< vector<int> > key_ICSSK = match_some_jets(full_coors, ICS_SK_coors);


    ratios << key_CS.size() << endl;
    for (int i = 0; i < key_CS.size(); i++){
      int fulltag = key_CS[i][0];
      int jettag = key_CS[i][1];
      if (jettag != -1){
        const PseudoJet &fulljet = root_jets[fulltag];
        const vector<double> &CSjet = corrected_CS_root_jet_params[jettag];
        ratios << CSjet[0]/fulljet.pt() << " " << CSjet[5]/fulljet.m() << endl;
      }
    }


    ratios << key_ICS.size() << endl;
    for (int i = 0; i < key_ICS.size(); i++){
      int fulltag = key_ICS[i][0];
      int jettag = key_ICS[i][1];
      if (jettag != -1){
        const PseudoJet &fulljet = root_jets[fulltag];
        const vector<double> &ICSjet = corrected_ICS_root_jet_params[jettag];
        ratios << ICSjet[0]/fulljet.pt() << " " << ICSjet[5]/fulljet.m() << endl;
      }
    }


    ratios << key_CSSK.size() << endl;
    for (int i = 0; i < key_CSSK.size(); i++){
      int fulltag = key_CSSK[i][0];
      int jettag = key_CSSK[i][1];
      if (jettag != -1){
        const PseudoJet &fulljet = root_jets[fulltag];
        const vector<double> &CSSKjet = corrected_CS_SK_root_jet_params[jettag];
        ratios << CSSKjet[0]/fulljet.pt() << " " << CSSKjet[5]/fulljet.m() << endl;
      }
    }


    ratios << key_ICSSK.size() << endl;
    for (int i = 0; i < key_ICSSK.size(); i++){
      int fulltag = key_ICSSK[i][0];
      int jettag = key_ICSSK[i][1];
      if (jettag != -1){
        const PseudoJet &fulljet = root_jets[fulltag];
        const vector<double> &ICSSKjet = corrected_ICS_SK_root_jet_params[jettag];
        ratios << ICSSKjet[0]/fulljet.pt() << " " << ICSSKjet[5]/fulljet.m() << endl;
      }
    }

    //And now we'll do this section for sanity checks.
    //Go full, CS, ICS, CSSK, ICSSK.
    //Write the number of jets, then on each line, pt, eta, mass
    checks << root_jets.size() << endl;
    for (int i = 0; i < root_jets.size(); i++){
      checks << root_jets[i].pt() << " " << root_jets[i].eta() << " " << root_jets[i].m() << endl;
    }
    checks << corrected_CS_root_jet_params.size() << endl;
    for (int i = 0; i < corrected_CS_root_jet_params.size(); i++){
      checks << corrected_CS_root_jet_params[i][0] << " " << corrected_CS_root_jet_params[i][6] << " " << corrected_CS_root_jet_params[i][5] << endl;
    }
    checks << corrected_ICS_root_jet_params.size() << endl;
    for (int i = 0; i < corrected_ICS_root_jet_params.size(); i++){
      checks << corrected_ICS_root_jet_params[i][0] << " " << corrected_ICS_root_jet_params[i][6] << " " << corrected_ICS_root_jet_params[i][5] << endl;
    }
    checks << corrected_CS_SK_root_jet_params.size() << endl;
    for (int i = 0; i < corrected_CS_SK_root_jet_params.size(); i++){
      checks << corrected_CS_SK_root_jet_params[i][0] << " " << corrected_CS_SK_root_jet_params[i][6] << " " << corrected_CS_SK_root_jet_params[i][5] << endl;
    }
    checks << corrected_ICS_SK_root_jet_params.size() << endl;
    for (int i = 0; i < corrected_ICS_SK_root_jet_params.size(); i++){
      checks << corrected_ICS_SK_root_jet_params[i][0] << " " << corrected_ICS_SK_root_jet_params[i][6] << " " << corrected_ICS_SK_root_jet_params[i][5] << endl;
    }


    //cout << "Wrote ICS_SK results." << endl;

//    for (int i = 0; i < jets_in_event.size(); i++){
//      const vector<double> &jet = jets_in_event[i];
//      int CS_tag = matching_key[i][1];
//      const vector<double> &CSjet = corrected_CS_root_jet_params[CS_tag];

//      results << (CSjet[0] - jet[0])/jet[0] << " " << (CSjet[1] - jet[1])/abs(jet[1]) << " " << phi_dist_alt(jet[2],CSjet[4])/abs(jet[2]) << " " << (CSjet[2] - jet[3])/jet[3] << " " << (CSjet[5]-jet[4])/jet[4] << endl;
//    }

//    for (int i = 0; i < jets_in_event.size(); i++){
//      const vector<double> &jet = jets_in_event[i];
//      int ICS_tag = ICS_root_key[i][1];
//      const vector<double> &ICSjet = corrected_ICS_root_jet_params[ICS_tag]; 
//      results << (ICSjet[0] - jet[0])/jet[0] << " " << (ICSjet[1] - jet[1])/abs(jet[1]) << " " << phi_dist_alt(jet[2],ICSjet[4])/abs(jet[2]) << " " << (ICSjet[2] - jet[3])/jet[3] << " " << (ICSjet[5]-jet[4])/jet[4] << endl;
//    }
    //cout << "Wrote ICS_SK results." << endl;

//    for (int i = 0; i < jets_in_event.size(); i++){
//      const vector<double> &jet = jets_in_event[i];
//      int CS_tag = matching_key[i][1];
//      const vector<double> &CSjet = corrected_CS_root_jet_params[CS_tag];

//      results << (CSjet[0] - jet[0])/jet[0] << " " << (CSjet[1] - jet[1])/abs(jet[1]) << " " << phi_dist_alt(jet[2],CSjet[4])/abs(jet[2]) << " " << (CSjet[2] - jet[3])/jet[3] << " " << (CSjet[5]-jet[4])/jet[4] << endl;
//    }

//    for (int i = 0; i < jets_in_event.size(); i++){
//      const vector<double> &jet = jets_in_event[i];
//      int ICS_tag = ICS_root_key[i][1];
//      const vector<double> &ICSjet = corrected_ICS_root_jet_params[ICS_tag]; 
//      results << (ICSjet[0] - jet[0])/jet[0] << " " << (ICSjet[1] - jet[1])/abs(jet[1]) << " " << phi_dist_alt(jet[2],ICSjet[4])/abs(jet[2]) << " " << (ICSjet[2] - jet[3])/jet[3] << " " << (ICSjet[5]-jet[4])/jet[4] << endl;
//    }
    
    //cout << "Here are the differences between the parameters: " << endl;
    //for (int i = 0; i < root_jets.size(); i++){
    //  const PseudoJet &jet = root_jets[i];
    //  const vector<double> &params = corrected_ICS_root_jet_params[i];
    //  cout << "ptdiff = " << jet.pt() - params[0] << ", rapdiff = " << jet.rap() - params[1] << ", mass = " << jet.m() - params[2] << endl;
    //  sketchy_diffs << abs(jet.pt() - params[0])/jet.pt() << " " << abs(jet.rap() - params[1])/jet.rap() << " " << abs(jet.m() - params[2])/jet.m() << " " << abs(width(jet) - params[3])/width(jet) << endl;
    //}

    //cout << endl;
*/    
  }

 // results.close();
 // ratios.close();
 // checks.close();


  mzero.close();
  return 0;
}

vector< vector<double> > get_coors(vector< vector<double> > params){
  vector< vector<double> > all_coors;

  for (int i = 0; i < params.size(); i++){
    vector<double> coors;
    coors.push_back(params[i][1]);
    coors.push_back(params[i][4]);
    all_coors.push_back(coors);
  }

  return all_coors;
}

//Based on the -1 in the matching keys, if no good match was made,
//this adjusts the information on the jet excesses
vector<int> adjust_excess(vector< vector<int> > key, vector<int> excesses){
  for (int i = 0; i < key.size(); i++){
    for (int j = 0; j < excesses.size(); j++){
      int excess = excesses[j];
      if (key[i][j+1] == -1){
        if (excess > 0){
	  excess = excess + 1;
        }
        if (excess < 0){
	  excess = excess - 2;
        }
        if (excess == 0){
	  int decide = rand() % 2;
	  if (decide == 0){
	    excess = excess + 2;
	  }
	  if (decide == 1){
	    excess = excess - 2;
	  }
        }
      }
      excesses[j] = excess;
    }
  }

  return excesses;
}

//vector<string> write_results(vector< vector<int> > key, vector< vector<double> > sigs, vector< vector<double> > params, int type){
//  vector<string> results;
//  for (int i = 0; i < key.size(); i++){
//    int Sigtag = key[i][0];
//    const vector<double> &sigjet = sigs[Sigtag];
//    int tag = key[i][type];
//    if (tag != -1){
//      const vector<double> &jet = params[tag];
//      string line = (jet[0] - sigjet[0])/sigjet[0] + " " + (jet[1] - sigjet[1])/sigjet[1] + " " + phi_dist_alt(sigjet[2],jet[4])/sigjet[2] + " " + (jet[2] - sigjet[3])/sigjet[3] + " " + (jet[5] - sigjet[4])/sigjet[4] + "\n";
//      results.push_back(line);
//    }
//  }

//  return results;
//}

vector< vector<double> > ConstituentSubtract(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, double max_distance, double alpha, double ghost_area, double grid_spacing, int common_bge, int cutoff){
  //Set up the subtractor with various parameters
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
  if (common_bge == 1){
    subtractor.set_common_bge_for_rho_and_rhom(true);
  }
  if (common_bge == 0){
    subtractor.set_common_bge_for_rho_and_rhom(false);
  }
  
  //Run the subtractor and cluster the new event
  //First, let's see if the subtractor is doing its job.
  //cout << "Number of particles in the full event, pre-subtraction: " << full_event.size() << endl;
  vector<PseudoJet> corrected_event = subtractor.subtract_event(full_event);
  //cout << "Number of particles in the full event, post-subtraction: " << corrected_event.size() << endl;
  ClusterSequence clust_seq_corr(corrected_event, jet_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets(cutoff));

  //cout << "Witihn ConstituentSubtract, here's the jet pt's:" << endl;
  //for (int i = 0; i < corrected_jets.size(); i++){
  //  cout << corrected_jets[i].pt() << endl;
  //}

  //Print the information from the corrected jet
  
  JetWidth width;

 // myfile << "Subtracted fulljets with CS algorithm: \n" <<endl;
  vector< vector<double> > jet_params;
  for (unsigned int i = 0; i < corrected_jets.size(); i++){
    const PseudoJet &jet = corrected_jets[i];
    vector<double> temp_params;
   // myfile << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << endl;
   temp_params.push_back(jet.pt());
   temp_params.push_back(jet.rap());
   temp_params.push_back(jet.E());
   temp_params.push_back(width(jet));
   temp_params.push_back(jet.phi_std());
   temp_params.push_back(jet.m());
   temp_params.push_back(jet.eta());
   jet_params.push_back(temp_params);
  }

  return jet_params;
}

//This method does everything the last one did, except it applies
//SoftKiller after CS.
vector< vector<double> > CS_SK(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, double max_distance, double alpha, double ghost_area, double grid_spacing, int common_bge, int cutoff, double SK_grid_size){
  //Set up the subtractor with various parameters
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
  if (common_bge == 1){
    subtractor.set_common_bge_for_rho_and_rhom(true);
  }
  if (common_bge == 0){
    subtractor.set_common_bge_for_rho_and_rhom(false);
  }
  
  //Run the subtractor and cluster the new event
  //First, let's see if the subtractor is doing its job.
  //cout << "Number of particles in the full event, pre-subtraction: " << full_event.size() << endl;
  vector<PseudoJet> corrected_event = subtractor.subtract_event(full_event);
 
  //For the first few particles, print out info for the event before and after soft_killing, for debugging
  //cout << "CS-subtracted event: " << endl;
  //for (int i = 0; i < 5; i++){
    //cout << "pt = " << corrected_event[i].pt() << ", rap = " << corrected_event[i].rap() << endl;
  //}

  //Now, I need to apply SoftKiller to the event.
  contrib::SoftKiller soft_killer(max_eta, SK_grid_size);
  double pt_threshold;
  vector<PseudoJet> soft_killed_event;
  soft_killer.apply(corrected_event, soft_killed_event, pt_threshold);

  //Okay, did soft_killer.apply do anything?
  //cout << "CS-SK subtracted event: " << endl;
  //for (int i = 0; i < 5; i++){
    //cout << "pt = " << soft_killed_event[i].pt() << ", rap = " << soft_killed_event[i].rap() << endl;
  //}

  //cout << "Number of particles in the full event, post-subtraction: " << corrected_event.size() << endl;
  ClusterSequence clust_seq_corr_kill(soft_killed_event, jet_def);
  vector<PseudoJet> corrected_killed_jets = sel_jets(clust_seq_corr_kill.inclusive_jets(cutoff));

  //cout << "Witihn ConstituentSubtract, here's the jet pt's:" << endl;
  //for (int i = 0; i < corrected_jets.size(); i++){
  //  cout << corrected_jets[i].pt() << endl;
  //}

  //Print the information from the corrected jet
  
  JetWidth width;

 // myfile << "Subtracted fulljets with CS algorithm: \n" <<endl;
  vector< vector<double> > jet_params;
  for (unsigned int i = 0; i < corrected_killed_jets.size(); i++){
    const PseudoJet &jet = corrected_killed_jets[i];
    vector<double> temp_params;
   // myfile << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << endl;
   temp_params.push_back(jet.pt());
   temp_params.push_back(jet.rap());
   temp_params.push_back(jet.E());
   temp_params.push_back(width(jet));
   temp_params.push_back(jet.phi_std());
   temp_params.push_back(jet.m());
   temp_params.push_back(jet.eta());
   jet_params.push_back(temp_params);
  }
  //cout << endl;
  //cout << "Inside the method, CS-SK jets are: " << endl;
  //for (int i = 0; i < corrected_killed_jets.size(); i++){
    //const PseudoJet &jet =  corrected_killed_jets[i];
    //cout << "pt: " << jet.pt() << ", rap: " << jet.rap() << ", phi: " << jet.phi_std() << ", E: " << jet.E() << endl;
  //}

  //cout << endl;

  return jet_params;
}


vector< vector<double> > IteratedConstituentSubtract(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, vector<double> max_distances, vector<double> alphas, double ghost_area, double grid_spacing, int common_bge, int cutoff, bool proxies){
 //Set up the subtractor with various parameters
  contrib::IterativeConstituentSubtractor subtractor;
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR);
  subtractor.set_parameters(max_distances, alphas);
  subtractor.set_remove_remaining_proxies(proxies);
  subtractor.set_ghost_area(ghost_area);
  subtractor.set_max_eta(max_eta);
  subtractor.initialize();

  //Set up the background estimator
  GridMedianBackgroundEstimator bge_rho(max_eta, grid_spacing);
  bge_rho.set_particles(full_event);
  subtractor.set_background_estimator(&bge_rho);
  if (common_bge == 1){
    subtractor.set_common_bge_for_rho_and_rhom(true);
  }
  if (common_bge == 0){
    subtractor.set_common_bge_for_rho_and_rhom(false);
  }

  //Run the subtractor and cluster the new event
  //cout << "Number of particles before ICS subtraction: " << full_event.size() << endl;
  //cout << "Here's data on the first ten particles before subtraction: " << endl;
  //for (int i  = 0; i < 10; i++){
  //  const PseudoJet &particle = full_event[i];
  //  cout << particle.pt() << endl;
  //}
  vector<PseudoJet> corrected_event = subtractor.subtract_event(full_event);
  //cout << "Number of particles after ICS subtraction: " << corrected_event.size() << endl;
  //cout << "Here's the data on the first tne particles after subtraction: " << endl;
  //for (int i  = 0; i < 10; i++){
  //  const PseudoJet &particle = corrected_event[i];
  //  cout << particle.pt() << endl;
  //}
  ClusterSequence clust_seq_corr(corrected_event, jet_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets(cutoff));
  //cout << "Here are the jets, right after subtraction." << endl;
  //for (int i = 0; i < corrected_jets.size(); i++){
  //  const PseudoJet &newjet = corrected_jets[i];
  //  cout << "pt = " << newjet.pt() << ", eta = " << newjet.eta() << ", mass = " << newjet.m() << endl;
  //}

  //Print the information from the corrected jet
  
  JetWidth width;
 
  vector< vector<double> > jet_params;

  for (unsigned int i = 0; i < corrected_jets.size(); i++){
    const PseudoJet &jet = corrected_jets[i];
    vector<double> temp_params;
    temp_params.push_back(jet.pt());
    temp_params.push_back(jet.rap());
    temp_params.push_back(jet.E());
    temp_params.push_back(width(jet));
    temp_params.push_back(jet.phi_std());
    temp_params.push_back(jet.m());
    temp_params.push_back(jet.eta());
    jet_params.push_back(temp_params);
   // myfile << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << endl;
  }

  //Return the corrected jets
  return jet_params;

}

vector< vector<double> > ICS_SK(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, vector<double> max_distances, vector<double> alphas, double ghost_area, double grid_spacing, int common_bge, int cutoff, double SK_grid_size, bool proxies){
 //Set up the subtractor with various parameters
  contrib::IterativeConstituentSubtractor subtractor;
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR);
  subtractor.set_parameters(max_distances, alphas);
  subtractor.set_remove_remaining_proxies(proxies);
  subtractor.set_ghost_area(ghost_area);
  subtractor.set_max_eta(max_eta);
  subtractor.initialize();

  //Set up the background estimator
  GridMedianBackgroundEstimator bge_rho(max_eta, grid_spacing);
  bge_rho.set_particles(full_event);
  subtractor.set_background_estimator(&bge_rho);
  if (common_bge == 1){
    subtractor.set_common_bge_for_rho_and_rhom(true);
  }
  if (common_bge == 0){
    subtractor.set_common_bge_for_rho_and_rhom(false);
  }

  //Run the subtractor and cluster the new event
  //cout << "Number of particles before ICS subtraction: " << full_event.size() << endl;
  //cout << "Here's data on the first ten particles before subtraction: " << endl;
  //for (int i  = 0; i < 10; i++){
  //  const PseudoJet &particle = full_event[i];
  //  cout << particle.pt() << endl;
  //}
  vector<PseudoJet> corrected_event = subtractor.subtract_event(full_event);
  //cout << "Number of particles after ICS subtraction: " << corrected_event.size() << endl;
  //cout << "Here's the data on the first tne particles after subtraction: " << endl;
  //for (int i  = 0; i < 10; i++){
  //  const PseudoJet &particle = corrected_event[i];
  //  cout << particle.pt() << endl;
  //}
  ClusterSequence clust_seq_corr(corrected_event, jet_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets(cutoff));
  //cout << "Here are the jets, right after subtraction." << endl;
  //for (int i = 0; i < corrected_jets.size(); i++){
  //  const PseudoJet &newjet = corrected_jets[i];
  //  cout << "pt = " << newjet.pt() << ", eta = " << newjet.eta() << ", mass = " << newjet.m() << endl;
  //}

  //Now, I need to apply SoftKiller to the event.
  contrib::SoftKiller soft_killer(max_eta, SK_grid_size);
  double pt_threshold;
  vector<PseudoJet> soft_killed_event;
  soft_killer.apply(corrected_event, soft_killed_event, pt_threshold);

  ClusterSequence clust_seq_corr_kill(soft_killed_event, jet_def);
  vector<PseudoJet> corrected_killed_jets = sel_jets(clust_seq_corr_kill.inclusive_jets(cutoff)); 

  //Print the information from the corrected jet
  
  JetWidth width;
 
  vector< vector<double> > jet_params;

  for (unsigned int i = 0; i < corrected_killed_jets.size(); i++){
    const PseudoJet &jet = corrected_killed_jets[i];
    vector<double> temp_params;
    temp_params.push_back(jet.pt());
    temp_params.push_back(jet.rap());
    temp_params.push_back(jet.E());
    temp_params.push_back(width(jet));
    temp_params.push_back(jet.phi_std());
    temp_params.push_back(jet.m());
    temp_params.push_back(jet.eta());
    jet_params.push_back(temp_params);
  }

  //Return the corrected jets
  return jet_params;

}



//This function clusters the full event without selecting the n hardest jets and judges where the dropoff in pt is to impose the NHardest cutoff
//Actually, making a relative cutoff isn't all that physically useful. Just ignore this.
int findNHardest(vector<PseudoJet> full_event, JetDefinition jet_def, double max_eta_jet){

  Selector sel_jets = SelectorAbsEtaMax(max_eta_jet);

  ClusterSequence clust_seq_full(full_event, jet_def);
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());
  int size = full_jets.size();
  double flipped_pts[size];
  for (int i = 0; i < full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    flipped_pts[i] = jet.pt();
  }
  sort(flipped_pts, flipped_pts + size);

  //The problem is, sort goes from least to greatest. Flip it
  double pts[size];
  for (int i = 0; i < size; i++){
    pts[i] = flipped_pts[size-i-1];
  }
 
  int cutoff = 1;
  for (int i = 0; i < size - 1; i++){
    if (pts[i]/pts[i+1] < 2){
      cutoff++;
    }
    else{
      break;
    }
  }

  return cutoff;
}

//This takes a string and returns a vector of substrings, using space as delimieter
vector<string> get_pieces(string line){

  int pos = -1;
  int backmark = 0;
  int length = 0;
  string temp;
  vector<string> tokens;
  int pos1;
  tokens.push_back(line.substr(0,line.find(' ')));
  while (line.find(' ',pos+1) != string::npos){
    pos = line.find(' ',backmark);
    length = line.find(' ',pos+1)-pos;
    temp = line.substr(pos+1,length);
    backmark = pos + 1;
    tokens.push_back(temp);
  }

  return tokens;
}

vector< vector<int> > match_jets(vector< vector<double> > hard_coors, vector< vector< double> > sub_coors){
  vector< vector<int> > key;

  for (int i = 0; i < hard_coors.size(); i++){
    vector<double> jet = hard_coors[i];
    int closest_sub = 0;
    double min_dist = 1000;
    for (int j = 0; j < sub_coors.size(); j++){
      if (R_dist(jet,sub_coors[j]) < min_dist){
        min_dist = R_dist(jet,sub_coors[j]);
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

vector< vector<int> > match_all_jets(vector< vector<double> > signal_coors, vector< vector<double> > CS_coors, vector< vector<double> > ICS_coors){

  vector < vector<int> > key;

  //Figure out which set of coordinates has the fewest number of jets. All our matching will be centered around those.
  vector < vector<double> > min_coors;
  if (signal_coors.size() <= CS_coors.size() && signal_coors.size() <= ICS_coors.size()){
    min_coors = signal_coors;
  } else if(CS_coors.size() <= ICS_coors.size()){
    min_coors = CS_coors;
  } else{
    min_coors = ICS_coors;
  }

  //Cycle through the smallest one, make matches.
  for (int i = 0; i < min_coors.size(); i++){
    vector<double> jet = min_coors[i];
    int closest_signal = 0;
    int closest_CS = 0;
    int closest_ICS = 0;
    double minDistSig = 1000;
    double  minDistCS = 1000;
    double  minDistICS = 1000;
    for (int j = 0; j < signal_coors.size(); j++){
      if (R_dist(jet, signal_coors[j]) < minDistSig){
	minDistSig = R_dist(jet, signal_coors[j]);
	closest_signal = j;
      }
    }
    for (int j = 0; j < CS_coors.size(); j++){
      if (R_dist(jet, CS_coors[j]) < minDistCS){
	minDistCS = R_dist(jet,CS_coors[j]);
	closest_CS = j;
      }
    }
    for (int j = 0; j < ICS_coors.size(); j++){
      if (R_dist(jet,ICS_coors[j]) < minDistICS){
	minDistICS = R_dist(jet, ICS_coors[j]);
	closest_ICS = j;
      }
    }
    vector<int> temp;

    temp.push_back(closest_signal);
    temp.push_back(closest_CS);
    temp.push_back(closest_ICS);

    key.push_back(temp);
  }
  //Now that the key is done, I need to prune it of matches which are farther
  //apart than a cutoff
  double cutoff = 0.4;
  vector<int> matches;
  vector<double> sig,CS,ICS;
  for (int i = 0; i < key.size(); i++){
    matches = key[i];
    sig = signal_coors[matches[0]];
    CS = CS_coors[matches[1]];
    ICS = ICS_coors[matches[2]];

    //If any of the matched jets are too far apart, change out their tag for -1.
    if (R_dist(sig,CS) > cutoff){
      key[i][1] = -1;
    }
    if (R_dist(sig,ICS) > cutoff){
      key[i][2] = -1;
    }
  }

  return key;
 
}

//This is a copy for when I'm just matching the full jets to the signal jets.
vector< vector<int> > match_some_jets(vector< vector<double> > signal_coors, vector< vector<double> > full_coors){

  vector < vector<int> > key;

  //Figure out which set of coordinates has the fewest number of jets. All our matching will be centered around those.
  vector < vector<double> > min_coors;
  if (signal_coors.size() <= full_coors.size()){
    min_coors = signal_coors;
  } else {
    min_coors = full_coors;
  }

  //Cycle through the smallest one, make matches.
  for (int i = 0; i < min_coors.size(); i++){
    vector<double> jet = min_coors[i];
    int closest_signal = 0;
    int closest_full = 0;
    double minDistSig = 1000;
    double  minDistFull = 1000;
    for (int j = 0; j < signal_coors.size(); j++){
      if (R_dist(jet, signal_coors[j]) < minDistSig){
	minDistSig = R_dist(jet,signal_coors[j]);
	closest_signal = j;
      }
    }
    for (int j = 0; j < full_coors.size(); j++){
      if (R_dist(jet, full_coors[j]) < minDistFull){
	minDistFull = R_dist(jet,full_coors[j]);
        closest_full = j;
      }
    }
    vector<int> temp;

    temp.push_back(closest_signal);
    temp.push_back(closest_full);

    key.push_back(temp);
  }

  //Now that the key is done, I need to prune it of matches jets which are farther
  //apart than a cutoff
  double cutoff = 0.4;
  vector<int> matches;
  vector<double> sig,full;
  for (int i = 0; i < key.size(); i++){
    matches = key[i];
    sig = signal_coors[matches[0]];
    full = full_coors[matches[1]];

    //If any of the matched jets are too far apart, change out their tag for -1.
    if (R_dist(sig,full) > cutoff){
      key[i][1] = -1;
    }
  }

  return key;
 
}


//Just prints out jet coordinates
void print_coors(vector< vector<double> > coors){
  for (int i = 0; i < coors.size(); i++){
    cout << coors[i][0] << " " << coors[i][1] << endl;
  }
  cout << endl;
}

//Prints out keys, but only the two-entry kind.
void print_key(vector< vector<int> > key){
  for (int i = 0; i < key.size(); i++){
    cout << key[i][0] << " " << key[i][1] << endl;
  }
  cout << endl;
}


//Right now, this is set up for phi to be in [-pi,pi].
double R_dist(vector<double> point1, vector<double> point2){
  double dist = pow(pow(point1[0] - point2[0],2) + pow(phi_dist_alt(point1[1],point2[1]),2),0.5);

  return dist;
}

//Note that this assumes phi in [0,2pi].
double phi_dist(double phi1, double phi2){
  double dist;
  double temp_dist = abs(phi1-phi2);
  double opt_1 = (atan(1)*8 - phi2) + (phi1);
  double opt_2 = (atan(1)*8 - phi1) + (phi2);
  if (temp_dist <= opt_1 && temp_dist <= opt_2){
    dist = temp_dist;
  } else if(opt_1 < opt_2) {
    dist = opt_1;
  } else {
    dist = opt_2;
  }

  return dist;
}

//Note that this assumes phi in [-pi,pi].
double phi_dist_alt(double phi1, double phi2){
  double dist;
  double temp;
  if (phi2 < phi1){
    temp = phi1;
    phi1 = phi2;
    phi2 = temp;
  }
  double opt_1 = phi2 - phi1;
  double opt_2 = (atan(1)*4 - phi2) + (phi1 + atan(1)*4);
  if (opt_1 < opt_2){
    dist = opt_1;
  } else {
    dist = opt_2;
  }

  return dist;
}
