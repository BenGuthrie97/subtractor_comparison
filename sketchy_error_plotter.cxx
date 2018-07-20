//this program reads off the deviations of the CS and ICS algorithms
//from the file jet_differences.txt

#include <iostream>
#include "TH1D.h"
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;

void sketchy_error_plotter(){
	TH1F* pt_histo = new TH1F("pt errors","Differences from the full event pt and ICS-subtracted event pt; Deviation from full event pt;# occurrences",20,0,pow(10,-15));
	TH1F* rap_histo = new TH1F("rap errors", "Differences from the full event rap and ICS-subtracted event rap; Deviation from full event rap;# occurences", 20, 0, 0.5*pow(10,-4));
	TH1F* m_histo = new TH1F("m errors", "Differences from the full event mass and ICS-subtracted event mass; Deviation from full event mass;# occurences", 20, 0, 0.001);
	TH1F* width_histo = new TH1F("width errors", "Differences from the full event width and ICS-subtracted event width; Deviation from full event width;# occurences", 20, 0, 0.5*pow(10,-4));
	
	ifstream myfile;
	myfile.open("sketchy_diffs.txt");
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
		int numjets = atoi(jets.c_str());
		for(int i = 0; i < numjets; i++){
			getline(myfile,line);
			numbers = get_pieces(line);
			pt_histo->Fill(atof(numbers[0].c_str()));
			rap_histo->Fill(atof(numbers[1].c_str()));
			m_histo->Fill(atof(numbers[2].c_str()));
			width_histo->Fill(atof(numbers[3].c_str()));
		}
	}

	TFile f1("sketchy_deviances.root","RECREATE");
	f1.mkdir("CS");
	f1.mkdir("ICS");
	f1.cd("/ICS");
	pt_histo->Write();
	rap_histo->Write();
	m_histo->Write();
	width_histo->Write();
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
