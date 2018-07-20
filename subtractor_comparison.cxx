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
int findNHardest(vector<PseudoJet> full_event, JetDefinition jet_def, double max_eta_jet);
vector<string> get_pieces(string line);
vector< vector<int> > match_jets(vector< vector<double> > hard_coors, vector< vector<double> > sub_coors);
double R_dist(vector<double> point1, vector<double> point2);
double phi_dist(double phi1, double phi2);
double phi_dist_alt(double phi1, double phi2);

int main(){

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
  double px, py, pz, E;
  string cluster;
  vector<string> cluster_data;
  int ncl;
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
        px = atof(cluster_data[0].c_str());
        py = atof(cluster_data[1].c_str());
        pz = atof(cluster_data[2].c_str());
        E = atof(cluster_data[3].c_str());
        PseudoJet temp_jet = PseudoJet(px,py,pz,E);
        root_event.push_back(temp_jet);
      }
      root_events.push_back(root_event);
    } else{
      cout << "Crap, the line-reading from ttbar_event_output screwed up." << endl;
    }
  }

  //Now, if everything is correct, the events should be fully loaded
  //For brevity, I'm loading the differences from the first fifty events.
  cout << "Number of events " << root_events.size() << endl;
  ofstream sketchy_diffs("sketchy_diffs.txt");
  for (int i = 0; i < 1000; i++){
    vector<PseudoJet> root_event = root_events[i];
    root_event = SelectorAbsEtaMax(max_eta)(root_event);
    cout << "There are " << root_event.size() << " particles in ttbar event number " << i << endl;

    ClusterSequence clust_seq_root(root_event, jet_def);
    vector<PseudoJet> root_jets = sel_jets(clust_seq_root.inclusive_jets(cutoff));

    //Now run the clustering algorithms.
//    vector < vector<double> > corrected_CS_root_jet_params = ConstituentSubtract(root_event, max_eta, sel_jets, jet_def, max_distance_CS, alpha_CS, ghost_area, bge_grid_spacing, common_bge, cutoff);
    vector< vector<double> > corrected_ICS_root_jet_params = IteratedConstituentSubtract(root_event, max_eta, sel_jets, jet_def, max_distances_ICS, alphas_ICS, ghost_area, bge_grid_spacing, common_bge, cutoff); 
 
    //cout << "There are " << root_jets.size() << " jets in the first event of the root data." << endl;
    //for (unsigned int i = 0; i < root_jets.size(); i++){
    //  const PseudoJet &jet = root_jets[i];
    //  cout << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << endl;
    //}

    //cout << endl;

//    cout << "There are " << corrected_CS_root_jet_params.size() << " jets in the first event of the root data after CS-subtraction." << endl;
//    for (unsigned int i = 0; i < corrected_CS_root_jet_params.size(); i++){
//      const vector<double> &params = corrected_CS_root_jet_params[i];
//      cout << "pt = " << params[0] << ", rap = " << params[1] << ", mass = " << params[2] << ", width = " << params[3] << endl;
//    }

    //cout << endl;

    //cout << "There are " << corrected_ICS_root_jet_params.size() << " jets in the first event of the root data after ICS-subtraction." << endl;
    //for (unsigned int i = 0; i < corrected_ICS_root_jet_params.size(); i++){
    //  const vector<double> &params = corrected_ICS_root_jet_params[i];
    //  cout << "pt = " << params[0] << ", rap = " << params[1] << ", mass = " << params[2] << ", width =" << params[3] << endl;
    //}


    //cout << "Here are the differences between the parameters: " << endl;
    for (int i = 0; i < root_jets.size(); i++){
      const PseudoJet &jet = root_jets[i];
      const vector<double> &params = corrected_ICS_root_jet_params[i];
      //cout << "ptdiff = " << jet.pt() - params[0] << ", rapdiff = " << jet.rap() - params[1] << ", mass = " << jet.m() - params[2] << endl;
      sketchy_diffs << abs(jet.pt() - params[0])/jet.pt() << " " << abs(jet.rap() - params[1])/jet.rap() << " " << abs(jet.m() - params[2])/jet.m() << " " << abs(width(jet) - params[3])/width(jet) << endl;
    }

    cout << endl;
    
  }

  hard_event = SelectorAbsEtaMax(max_eta)(hard_event);
  full_event = SelectorAbsEtaMax(max_eta)(full_event);

  //cout << "There are " << hard_event.size() << " particles in the hard event." << endl;
  //cout << "There are " << full_event.size() << " particles in the full event." << endl;

  //Make the jet definition and the selector for the jets.
  //int cutoff = findNHardest(full_event, jet_def, max_eta_jet);
  //cout << "Studying the  " << cutoff << " hardest jets." << endl;
  //Selector sel_jets = SelectorNHardest(cutoff) * SelectorAbsEtaMax(max_eta_jet);
  
  //Peter recommends not actually using the cutoff. I'll just keep the AbsEtaMax selector for now and see what happens.
  //I'll also insert an artificial cutoff in units of... I think GeV?

  //run the clustering on the hard jets and the full jets
  ClusterSequence clust_seq_hard(hard_event, jet_def);
  ClusterSequence clust_seq_full(full_event, jet_def);
  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets(cutoff));
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets(cutoff));

  //Now pull in the subtracted jets
  vector< vector<double> > corrected_CS_jet_params = ConstituentSubtract(full_event, max_eta, sel_jets, jet_def, max_distance_CS, alpha_CS, ghost_area, bge_grid_spacing, common_bge, cutoff);
  vector< vector<double> > corrected_ICS_jet_params = IteratedConstituentSubtract(full_event, max_eta, sel_jets, jet_def, max_distances_ICS, alphas_ICS, ghost_area, bge_grid_spacing, common_bge, cutoff);

  vector< vector <double> > hard_jet_coors;
  for (int i = 0; i < hard_jets.size(); i++){
    vector <double> coors;
    coors.push_back(hard_jets[i].eta());
    coors.push_back(hard_jets[i].phi());
    hard_jet_coors.push_back(coors);
  }

  for (int i = 0; i < hard_jet_coors.size(); i++){
    cout << hard_jet_coors[i][0] << " " << hard_jet_coors[i][1] << endl;
  }

  //Find all the comparable jets using deltaR matching.
  //First, get their coordinates.
  vector < vector<double> > hard_coors;
  for (int i = 0; i < hard_jets.size(); i++){
    vector <double> temp;
    temp.push_back(hard_jets[i].eta());
    temp.push_back(hard_jets[i].phi());
    hard_coors.push_back(temp);
  }

  vector < vector<double> > CS_sub_coors;
  for (int i = 0; i < corrected_CS_jet_params.size(); i++){
    vector <double> temp;
    temp.push_back(corrected_CS_jet_params[i][1]);
    temp.push_back(corrected_CS_jet_params[i][4]);
    CS_sub_coors.push_back(temp);
  }

  vector < vector<double> > ICS_sub_coors;
  for (int i = 0; i < corrected_ICS_jet_params.size(); i++){
    vector <double> temp;
    temp.push_back(corrected_ICS_jet_params[i][1]);
    temp.push_back(corrected_ICS_jet_params[i][4]);
    ICS_sub_coors.push_back(temp);
  }

  //Now make the matching key.
  vector< vector<int> > CS_key;
  CS_key = match_jets(hard_coors,CS_sub_coors);
  vector < vector<int> > ICS_key;
  ICS_key = match_jets(hard_coors, ICS_sub_coors);


  //Write the results of the jets to a text file, but only the ones
  //which can be matched to a hard jet with deltaR matching.
  //Print everything for good measure.
  cout << setprecision(8);
  //ofstream myfile;
  //myfile.open("subtractor_comparison_results.txt");

  //myfile << hard_jets.size() << endl;

  //cout << "# original hard jets" << endl;
  //cout << "Number of hard jets: " << hard_jets.size() << endl;
 // myfile << "Original hard jets \n" << endl;
  //for (unsigned int i = 0; i < hard_jets.size(); i++){
  //  const PseudoJet &jet = hard_jets[i];
  //  cout << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << endl;
   // myfile << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << "\n" << endl;
  //  myfile << jet.pt() << " "  << jet.rap() << " " << jet.m() << " " << width(jet) << endl;
  //}
  //cout << endl;

  //cout << "# unsubtracted full jets" << endl;
  //cout << "Number of full jets: " << full_jets.size() << endl;
 // myfile << "Unsubtracted full jets \n" << endl;
  //for (unsigned int i = 0; i < full_jets.size(); i++){
  //  const PseudoJet &jet = full_jets[i];
  //  cout << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << endl;
   // myfile << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << "\n" << endl;
  //  myfile << jet.pt() << " " << jet.rap() << " " << jet.m() << " " << width(jet) << endl;
  //}
  //cout << endl;

  //cout << "# subtracted full jets with CS algorithm" << endl;
 // myfile << "Subtracted full jets with CS algorithm \n" << endl;
  //cout << "Number of CS-subtracted jets: " << corrected_CS_jet_params.size() << endl;
  //for (unsigned int i = 0; i < corrected_CS_jet_params.size(); i++){
    //const vector<double> &params = corrected_CS_jet_params[i];
    //cout << "pt = " << params[0] << ", rap = " << params[1] << ", mass = " << params[2] << ", width = " << params[3] << endl;
   // myfile << "pt = " << params[0] << ", rap = " << params[1] << ", mass = " << params[2] << ", width = " << params[3] << "\n" << endl;
   // myfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] << endl;
  //}

  //cout << endl;

  //cout << "# subtracted full jets with ICS algorithm" << endl;
 // myfile << "Subtracted full jets with ICS algorithm \n" << endl;
  //cout << "Number of ICS-subtracted jets: " << corrected_ICS_jet_params.size() << endl;
  //for (unsigned int i = 0; i < corrected_ICS_jet_params.size(); i++){
  //  const vector<double> &params = corrected_ICS_jet_params[i];
  //  cout << "pt = " << params[0] << ", rap = " << params[1] << ", mass = " << params[2] << ", width = " << params[3] << endl;
   // myfile << "pt = " << params[0] << ", rap = " << params[1] << ", mass = " << params[2] << ", width = " << params[3] << "\n" << endl;
  //  myfile << params[0] << " " << params[1] << " " << params[2] << " " << params[3] << endl;
  //}

  //myfile.close();

  //Now, we need to get the percent error for the subtraction algorithms.
  ofstream results;
  results.open("jet_differences.txt");
  results << hard_jets.size() << endl;
  for (int i = 0; i < hard_jets.size(); i++){
    PseudoJet jet = hard_jets[i];
    int CS_tag = CS_key[i][1];
    vector<double> CSjet = corrected_CS_jet_params[CS_tag];
    results << abs(jet.pt() - CSjet[0])/jet.pt() << " " << abs(jet.rap() - CSjet[1])/jet.rap() << " " << abs(jet.m() - CSjet[2])/jet.m() << " " << abs(width(jet) - CSjet[3])/width(jet) << endl;
  }
  for (int i = 0; i < hard_jets.size(); i++){
    PseudoJet jet = hard_jets[i];
    int ICS_tag = ICS_key[i][1];
    vector<double> ICSjet = corrected_ICS_jet_params[ICS_tag]; 
    results << abs(jet.pt() - ICSjet[0])/jet.pt() << " " << abs(jet.rap() - ICSjet[1])/jet.rap() << " " << abs(jet.m() - ICSjet[2])/jet.m() << " " << abs(width(jet) - ICSjet[3])/width(jet) << endl;
  }

  return 0;
}

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
   temp_params.push_back(jet.m());
   temp_params.push_back(width(jet));
   temp_params.push_back(jet.phi());
   jet_params.push_back(temp_params);
  }
  cout << endl;

  return jet_params;
}

vector< vector<double> > IteratedConstituentSubtract(vector<PseudoJet> full_event, double max_eta, Selector sel_jets, JetDefinition jet_def, vector<double> max_distances, vector<double> alphas, double ghost_area, double grid_spacing, int common_bge, int cutoff){
 //Set up the subtractor with various parameters
  contrib::IterativeConstituentSubtractor subtractor;
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

 // myfile << "Subtracted full jets with ICS algorithm \n" << endl;
 
  vector< vector<double> > jet_params;

  for (unsigned int i = 0; i < corrected_jets.size(); i++){
    const PseudoJet &jet = corrected_jets[i];
    vector<double> temp_params;
    temp_params.push_back(jet.pt());
    temp_params.push_back(jet.rap());
    temp_params.push_back(jet.m());
    temp_params.push_back(width(jet));
    temp_params.push_back(jet.phi());
    jet_params.push_back(temp_params);
   // myfile << "pt = " << jet.pt() << ", rap = " << jet.rap() << ", mass = " << jet.m() << ", width = " << width(jet) << endl;
  }
  cout << endl;

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
