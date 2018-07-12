//This program reads off jet parameter deviations of the CS- and
//ICS-subtracted jets' parameters from the hard jets' parameters,
//then plots them on ROOT histograms.

#include<iostream>
#include "TH1D.h"
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;

void error_plotter() {
	TH1F* CS_histo = new TH1F("errors", "Errors from the CS algorithm; Deviation from the hard jet;# occurrences", 20, 0, 0.5);

	ifstream myfile;
	myfile.open("jet_differences.txt");
	string line;
	vector<string> numbers;
	string jets;
	getline(myfile, jets);
	int numjets = atoi(jet.c_str());
	for (int i = 0; i < numjets; i++) {
		getline(myfile, line);
		numbers = get_pieces(line);
		for (int j = 0; j < numbers.size(); j++) {
			CS_histo->Fill(atof(numbers[j].c_str()));
		}
	}
	CS_histo->Draw();

	//Now, make another one but for ICS errors
	TH1F* ICS_histo = new TH1F("ICS_errors", "Errors from the ICS algorithm; Deviation from the hard jet; # occurrences", 20, 0, 0.5);
	for (int i = 0; i < numjets; i++) {
		getline(myfile, line);
		numbers = get_pieces(line);
		for (int j = 0; j < numbers.size(); j++) {
			ICS_histo->Fill(atof(numbers[j].c_str()));
		}
	}
	ICS_hist->Draw();
}

vector<string> get_Pieces(string line) {

	int pos = -1;
	int backmark = 0;
	int length = 0;
	string temp;
	vector<string> numbers;
	int pos1;
	numbers.push_back(line.substr(0, line.find(' ')));
	while (line.find(' ', pos + 1) != string::npos) {
		pos = line.find(' ', backmark);
		length = line.find(' ', pos + 1) - pos;
		temp = line.substr(pos + 1, length);
		backmark = pos + 1;
		numbers.push_back(temp);
	}

	return numbers;
}