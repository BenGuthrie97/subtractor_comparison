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

	//Histograms for deviances of the full jets
	TH1F* pt_full_histo = new TH1F("pt full errors", "Deviances from the signal jets' pt and full jets' pt; Deviation; # occurrences",40,-1,1);
	TH1F* rap_full_histo = new TH1F("rapidity full errors", "Deviances from the signal jets' rap and full jets' rap; Deviation; # occurrences",40,-1,1);
	TH1F* phi_full_histo = new TH1F("phi full errors", "Deviances from the signal jets' phi and full jets' phi; Deviation; # occurrences",40,-1,1);
	TH1F* E_full_histo = new TH1F("E full errors", "Deviances from the signal jets' energy and full jets' energy; Deviation; # occurrences",40,-1,1);
	TH1F* m_full_histo = new TH1F("m full errors", "Deviances from the signal jets' mass and full jets' mass; Deviation; # occurrences", 40,-1,1);


	//Make the histograms for CS deviances
	TH1F* pt_CS_histo = new TH1F("pt CS errors", "Deviances from the signal jets' pt and CS-subtracted jets' pt; Deviation; # occurrences",40,-1,1);
	TH1F* rap_CS_histo = new TH1F("rapidity CS errors", "Deviances from the signal jets' rap and CS-subtracted jets' rap; Deviation; # occurrences",40,-1,1);
	TH1F* phi_CS_histo = new TH1F("phi CS errors", "Deviances from the signal jets' phi and CS-subtracted jets' phi; Deviation; # occurrences",40,-1,1);
	TH1F* E_CS_histo = new TH1F("E CS errors", "Deviances from the signal jets' energy and CS-subtracted jets' energy; Deviation; # occurrences",40,-1,1);
	TH1F* m_CS_histo = new TH1F("m CS errors", "Deviances from the signal jets' mass and CS-subtracted jets' mass; Deviation; # occurrences", 40,-1,1);


	//Make the histograms for CS_SK deviances
	TH1F* pt_CSSK_histo = new TH1F("pt CSSK errors", "Deviances from the signal jets' pt and CSSK-subtracted jets' pt; Deviation; # occurrences",40,-1,1);
	TH1F* rap_CSSK_histo = new TH1F("rapidity CSSK errors", "Deviances from the signal jets' rap and CSSK-subtracted jets' rap; Deviation; # occurrences",40,-1,1);
	TH1F* phi_CSSK_histo = new TH1F("phi CSSK errors", "Deviances from the signal jets' phi and CSSK-subtracted jets' phi; Deviation; # occurrences",40,-1,1);
	TH1F* E_CSSK_histo = new TH1F("E CSSK errors", "Deviances from the signal jets' energy and CSSK-subtracted jets' energy; Deviation; # occurrences",40,-1,1);
	TH1F* m_CSSK_histo = new TH1F("m CSSK errors", "Deviances from the signal jets' mass and CSSK-subtracted jets' mass; Deviation; # occurrences", 40,-1,1);

	//Make the same histograms, but for ICS deviances
	TH1F* pt_ICS_histo = new TH1F("pt ICS errors", "Deviances from the signal jets' pt and ICS-subtracted jets' pt; Deviation; # occurrences",40,-1,1);
	TH1F* rap_ICS_histo = new TH1F("rapidity ICS errors", "Deviances from the signal jets' rap and ICS-subtracted jets' rap; Deviation; # occurrences",40,-1,1);
	TH1F* phi_ICS_histo = new TH1F("phi ICS errors", "Deviances from the signal jets' phi and ICS-subtracted jets' phi; Deviation; # occurrences",40,-1,1);
	TH1F* E_ICS_histo = new TH1F("E ICS errors", "Deviances from the signal jets' energy and ICS-subtracted jets' energy; Deviation; # occurrences",40,-1,1);	
	TH1F* m_ICS_histo = new TH1F("m ICS errors", "Deviances from the signal jets' mass and ICS-subtracted jets' mass; Deviation; # occurrences",40,-1,1);

	//Make the same histograms, but for ICSSK deviances
	TH1F* pt_ICSSK_histo = new TH1F("pt ICSSK errors", "Deviances from the signal jets' pt and ICSSK-subtracted jets' pt; Deviation; # occurrences",40,-1,1);
	TH1F* rap_ICSSK_histo = new TH1F("rapidity ICSSK errors", "Deviances from the signal jets' rap and ICSSK-subtracted jets' rap; Deviation; # occurrences",40,-1,1);
	TH1F* phi_ICSSK_histo = new TH1F("phi ICSSK errors", "Deviances from the signal jets' phi and ICSSK-subtracted jets' phi; Deviation; # occurrences",40,-1,1);
	TH1F* E_ICSSK_histo = new TH1F("E ICSSK errors", "Deviances from the signal jets' energy and ICSSK-subtracted jets' energy; Deviation; # occurrences",40,-1,1);	
	TH1F* m_ICSSK_histo = new TH1F("m ICSSK errors", "Deviances from the signal jets' mass and ICSSK-subtracted jets' mass; Deviation; # occurrences",40,-1,1);

	//Don't forget the ones for jet excesses
	TH1F* excess_CS = new TH1F("number of excess CS jets","Number of jets present in CS-subtracted event not present in hard event; # excess jets; # occurrences",21,-5,15);
	TH1F* excess_ICS = new TH1F("number of excess ICS jets","Number of jets present in ICS-subtracted event not present in hard event; # excess jets; # occurrences",21,-5,15);
	TH1F* excess_full = new TH1F("number of excess full jets","Number of jets present in full event not present in hard event; # excess jets; # occurrences",21,-5,15);
	TH1F* excess_CSSK = new TH1F("number of excess CS_SK jets","Number of jets present in CS_SK-subtracted event not present in hard event; # excess jets; # occurrences",21,-5,15);
	TH1F* excess_ICSSK = new TH1F("number of excess ICS_SK jets","Number of jets present in ICS_SK-subtracted event not present in hard event; # excess jets; # occurrences",21,-5,15);


	ifstream myfile;
	myfile.open("jet_differences.txt");
	string line;
	vector<string> numbers;
	int numFull, numCS, numICS, numCSSK, numICSSK;
	//while (getline(myfile,line)){
	//	numbers = get_pieces(line);
	//	for (int i = 0; i < numbers.size(); i++){
	//		histo->Fill(atof(numbers[i].c_str()));
	//	}
	//}
	string jets;
	while (getline(myfile,jets)){

                //Get the number of each king of jet in this event
		vector<string> numjets = get_pieces(jets);


                //The next line is how many jets the subtractions have in excess of the "true" value.
                getline(myfile, line);
                numbers = get_pieces(line);
		excess_full->Fill(atoi(numbers[0].c_str()));
                excess_CS->Fill(atoi(numbers[1].c_str()));
                excess_ICS->Fill(atoi(numbers[2].c_str()));
		excess_CSSK->Fill(atoi(numbers[3].c_str()));
		excess_ICSSK->Fill(atoi(numbers[4].c_str()));


		//Quit if there were too few jets for proper analysis. May be adjusted later.
		if (atoi(numjets[0].c_str()) == -1){
			continue;
		}

		//Now that we know there's actual data to take, figure out how many of each to take
		numFull = atoi(numjets[0].c_str());
		numCS = atoi(numjets[1].c_str());
		numICS = atoi(numjets[2].c_str());
		numCSSK = atoi(numjets[3].c_str());
		numICSSK = atoi(numjets[4].c_str());


                //Now get all the full differences. Each line is for one jet
		for(int i = 0; i < numFull; i++){
			getline(myfile,line);
			numbers = get_pieces(line);
			pt_full_histo->Fill(atof(numbers[0].c_str()));
			rap_full_histo->Fill(atof(numbers[1].c_str()));
			phi_full_histo->Fill(atof(numbers[2].c_str()));
			E_full_histo->Fill(atof(numbers[3].c_str()));
			m_full_histo->Fill(atof(numbers[4].c_str()));
		}
		

                //Now get all the CS differences. Each line is for one jet
		for(int i = 0; i < numCS; i++){
			getline(myfile,line);
			numbers = get_pieces(line);
			pt_CS_histo->Fill(atof(numbers[0].c_str()));
			rap_CS_histo->Fill(atof(numbers[1].c_str()));
			phi_CS_histo->Fill(atof(numbers[2].c_str()));
			E_CS_histo->Fill(atof(numbers[3].c_str()));
			m_CS_histo->Fill(atof(numbers[4].c_str()));
		}

		//Do the same for ICS.
		for(int i = 0; i < numICS; i++){
			getline(myfile,line);
			numbers = get_pieces(line);
			pt_ICS_histo->Fill(atof(numbers[0].c_str()));
			rap_ICS_histo->Fill(atof(numbers[1].c_str()));
			phi_ICS_histo->Fill(atof(numbers[2].c_str()));
			E_ICS_histo->Fill(atof(numbers[3].c_str()));
			m_ICS_histo->Fill(atof(numbers[4].c_str()));
		}

                //Same for CSSK
		for(int i = 0; i < numCSSK; i++){
			getline(myfile,line);
			numbers = get_pieces(line);
			pt_CSSK_histo->Fill(atof(numbers[0].c_str()));
			rap_CSSK_histo->Fill(atof(numbers[1].c_str()));
			phi_CSSK_histo->Fill(atof(numbers[2].c_str()));
			E_CSSK_histo->Fill(atof(numbers[3].c_str()));
			m_CSSK_histo->Fill(atof(numbers[4].c_str()));
		}

                //Same for ICSSK
		for(int i = 0; i < numICSSK; i++){
			getline(myfile,line);
			numbers = get_pieces(line);
			pt_ICSSK_histo->Fill(atof(numbers[0].c_str()));
			rap_ICSSK_histo->Fill(atof(numbers[1].c_str()));
			phi_ICSSK_histo->Fill(atof(numbers[2].c_str()));
			E_ICSSK_histo->Fill(atof(numbers[3].c_str()));
			m_ICSSK_histo->Fill(atof(numbers[4].c_str()));
		}


	}

	//It would be nice to put these into folders, but I won't mess with that now
	TFile f1("jet_deviances.root","RECREATE");
	f1.mkdir("full");
	f1.cd("/full");
	pt_full_histo->Write();
	rap_full_histo->Write();
	phi_full_histo->Write();
	E_CS_histo->Write();
	m_CS_histo->Write();
	excess_full->Write();
	f1.cd("..");
	f1.mkdir("CS");
	f1.cd("/CS");
	pt_CS_histo->Write();
	rap_CS_histo->Write();
	phi_CS_histo->Write();
	E_CS_histo->Write();
	m_CS_histo->Write();
	excess_CS->Write();
	f1.cd("..");
	f1.mkdir("ICS");
	f1.cd("/ICS");
	pt_ICS_histo->Write();
	rap_ICS_histo->Write();
	phi_ICS_histo->Write();
	E_ICS_histo->Write();
	m_ICS_histo->Write();
	excess_ICS->Write();
	f1.cd("..");
	f1.mkdir("CSSK");
	f1.cd("/CSSK");	
	pt_CSSK_histo->Write();
	rap_CSSK_histo->Write();
	phi_CSSK_histo->Write();
	E_CSSK_histo->Write();
	m_CSSK_histo->Write();
	excess_CSSK->Write();
	f1.cd("..");
	f1.mkdir("ICSSK");
	f1.cd("/ICSSK");
	pt_ICSSK_histo->Write();
	rap_ICSSK_histo->Write();
	phi_ICSSK_histo->Write();
	E_ICSSK_histo->Write();
	m_ICSSK_histo->Write();
	excess_ICSSK->Write();


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
