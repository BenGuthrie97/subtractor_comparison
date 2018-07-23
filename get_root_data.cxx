#include<iostream>

void get_root_data(){
  ofstream myfile;
  myfile.open("ttbar_event_output.txt");
  
  TFile *f = TFile::Open("ttbar.root");
  f->cd("PlaintiveAmericanbobtail");
  TTree *t1 = (TTree*)f->Get("PlaintiveAmericanbobtail/nominal");
  vector<float> *cl_pt;
  vector<float> *cl_eta;
  vector<float> *cl_m;
  vector<float> *cl_phi;
  Int_t ncl;

  t1->SetBranchAddress("cl_pt",&cl_pt);
  t1->SetBranchAddress("cl_eta", &cl_eta);
  t1->SetBranchAddress("cl_m",&cl_m);
  t1->SetBranchAddress("cl_phi",&cl_phi);
  t1->SetBranchAddress("ncl", &ncl);

  Int_t nentries = (Int_t)t1->GetEntries();

  double pt, eta, phi, m, px, py , pz, E;
  for (int i = 0; i < nentries; i++){
    t1->GetEntry(i);
    
    //First, write the number of particles in this entry to the file so the reader
    //knows how long it has to read for this event
    myfile << ncl << endl;

    //then go through, particle by particle, and convert the data to what PseudoJet needs
    for (int j = 0; j < ncl; j++){
      pt = cl_pt->at(j);
      eta = cl_eta->at(j);
      phi = cl_phi->at(j);
      m = cl_m->at(j);

      px = pt * sin(phi);
      py = pt * cos(phi);
      pz = 0.5*exp(-1*eta)*(exp(2*eta) - 1)*pt;
      E = sqrt(pow(pt,2) + pow(pz,2) + pow(m,2));

      //Correct sign ambiguity
      if (eta < 0) {
        pz = pz * (-1);
      }
      //and write it to the file
      myfile << px << " " << py << " " << pz << " " << E << endl;
      }
    }

  myfile.close();
}
