//this program reads off the deviations of the CS and ICS algorithms
//from the file jet_differences.txt

#include <iostream>
#include "TH1D.h"
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;

void jet_deviance_plotter(){

	//Make the histograms for CS deviances
	TH1F* pt_CS_histo = new TH1F("pt CS errors", "Deviances from the signal jets' pt and CS-subtracted jets' pt; Deviation; # occurrences",20,0,1);
	TH1F* rap_CS_histo = new TH1F("rapidity CS errors", "Deviances from the signal jets' rap and CS-subtracted jets' rap; Deviation; # occurrences",20,0,1);
	TH1F* phi_CS_histo = new TH1F("phi CS errors", "Deviances from the signal jets' phi and CS-subtracted jets' phi; Deviation; # occurrences",20,0,1);
	TH1F* E_CS_histo = new TH1F("E CS errors", "Deviances from the signal jets' energy and CS-subtracted jets' energy; Deviation; # occurrences",20,0,1);

	//Make the same histograms, but for ICS deviances
	TH1F* pt_ICS_histo = new TH1F("pt ICS errors", "Deviances from the signal jets' pt and ICS-subtracted jets' pt; Deviation; # occurrences",20,0,1);
	TH1F* rap_ICS_histo = new TH1F("rapidity ICS errors", "Deviances from the signal jets' rap and ICS-subtracted jets' rap; Deviation; # occurrences",20,0,1);
	TH1F* phi_ICS_histo = new TH1F("phi ICS errors", "Deviances from the signal jets' phi and ICS-subtracted jets' phi; Deviation; # occurrences",20,0,1);
	TH1F* E_ICS_histo = new TH1F("E ICS errors", "Deviances from the signal jets' energy and ICS-subtracted jets' energy; Deviation; # occurrences",20,0,1);	

	//Don't forget the ones for jet excesses
	TH1F* excess_CS = new TH1F("number of excess CS jets","Number of jets present in CS-subtracted event not present in hard event; # excess jets; # occurrences",16,0,15);
	TH1F* excess_ICS = new TH1F("number of excess ICS jets","Number of jets present in ICS-subtracted event not present in hard event; # excess jets; # occurrences",16,0,15);
	

	ifstream myfile;
	myfile.open("jet_differences.txt");
	string line;
	vector<string> numbers;
	//while (getline(myfile,line)){
	//	numbers = get_pieces(line);
	//	for (int i = 0; i < numbers.size(); i++){
	//		histo->Fill(atof(numbers[i].c_str()));
	//	}
	//}
	string jets;
	while (getline(myfile,jets)){
                //Get the number of jets in this event
		int numjets = atoi(jets.c_str());

                //The next line is how many jets the subtractions have in excess of the "true" value.
                getline(myfile, line);
                numbers = get_pieces(line);
                excess_CS->Fill(atoi(numbers[0].c_str()));
                excess_ICS->Fill(atoi(numbers[1].c_str()));

                //Now get all the CS differences. Each line is for one jet
		for(int i = 0; i < numjets; i++){
			getline(myfile,line);
			numbers = get_pieces(line);
			pt_CS_histo->Fill(atof(numbers[0].c_str()));
			rap_CS_histo->Fill(atof(numbers[1].c_str()));
			phi_CS_histo->Fill(atof(numbers[2].c_str()));
			E_CS_histo->Fill(atof(numbers[3].c_str()));
		}

		//Do the same for ICS.
		for(int i = 0; i < numjets; i++){
			getline(myfile,line);
			numbers = get_pieces(line);
			pt_ICS_histo->Fill(atof(numbers[0].c_str()));
			rap_ICS_histo->Fill(atof(numbers[1].c_str()));
			phi_ICS_histo->Fill(atof(numbers[2].c_str()));
			E_ICS_histo->Fill(atof(numbers[3].c_str()));
			
		}
	}

	TFile f1("jet_deviances.root","RECREATE");
	f1.mkdir("CS");
	f1.mkdir("ICS");
	f1.cd("/CS");
	pt_CS_histo->Write();
	rap_CS_histo->Write();
	phi_CS_histo->Write();
	E_CS_histo->Write();
	excess_CS->Write();
	f1.cd("../ICS");
	pt_ICS_histo->Write();
	rap_ICS_histo->Write();
	phi_ICS_histo->Write();
	E_ICS_histo->Write();
	excess_ICS->Write();
}

vector<string> get_pieces(string line){

	int pos = -1;
	int backmark = 0;
	int length = 0;
	string temp;
	vector<string> numbers;
	int pos1;
	numbers.push_back(line.substr(0,line.find(' ')));
	while(line.find(' ',pos+1)!=string::npos){
		pos = line.find(' ',backmark);
		length = line.find(' ',pos+1)-pos;
		temp = line.substr(pos+1,length);
		backmark = pos+1;
		numbers.push_back(temp);
	}

	return numbers;

}
