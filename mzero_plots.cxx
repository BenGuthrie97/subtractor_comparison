#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace std;

void mzero_plots(){

	TH1F* CS_pt = new TH1F("pt CS jets","p_t distribution of CS-subtracted jets; p_t; # occurrences",60,0,300);
	TH1F* CS_eta = new TH1F("eta CS jets","eta distribution of CS-subtracted jets; eta; # occurrences",40,-2.5,2.5);
	TH1F* CS_m = new TH1F("m CS jets", "mass distribution of CS-subtracted jets; mass; # occurrences",100,0,300);
	TH1F* CS_E = new TH1F("E CS jets", "energy distribution of CS-subtracted jets; energy; # occurrences",100,0,300);
	TH1F* ICS_pt = new TH1F("pt ICS jets","p_t distribution of ICS-subtracted jets; p_t; # occurrences",60,0,300);
	TH1F* ICS_eta = new TH1F("eta ICS jets","eta distribution of ICS-subtracted jets; eta; # occurrences",40,-2.5,2.5);
	TH1F* ICS_m = new TH1F("m ICS jets", "mass distribution of ICS-subtracted jets; mass; # occurrences",100,0,300);
	TH1F* ICS_E = new TH1F("E ICS jets", "energy distribution of ICS-subtracted jets; energy; # occurrences",100,0,300);
	TH1F* CSSK_pt = new TH1F("pt CSSK jets","p_t distribution of CSSK-subtracted jets; p_t; # occurrences",60,0,300);
	TH1F* CSSK_eta = new TH1F("eta CSSK jets","eta distribution of CSSK-subtracted jets; eta; # occurrences",40,-2.5,2.5);
	TH1F* CSSK_m = new TH1F("m CSSK jets", "mass distribution of CSSK-subtracted jets; mass; # occurrences",100,0,300);
	TH1F* CSSK_E = new TH1F("E CSSK jets", "energy distribution of CSSK-subtracted jets; energy; # occurrences",100,0,300);
	TH1F* ICSSK_pt = new TH1F("pt ICSSK jets","p_t distribution of ICSSK-subtracted jets; p_t; # occurrences",60,0,300);
	TH1F* ICSSK_eta = new TH1F("eta ICSSK jets","eta distribution of ICSSK-subtracted jets; eta; # occurrences",40,-2.5,2.5);
	TH1F* ICSSK_m = new TH1F("m ICSSK jets", "mass distribution of ICSSK-subtracted jets; mass; # occurrences",100,0,300);
	TH1F* ICSSK_E = new TH1F("energy ICSSK jets", "energy distribution of ICSSK-subtracted jets; energy; # occurrences",100,0,300);

	ifstream data;
	data.open("mzero_jets.txt");
	string line;
	string j;
	vector<string> numbers;

	while(getline(data,line)){
		//Right now, the line holds the event number
		j = line;
		getline(data,line);
		//Now it holds "CS". Check that
		if (line != "CS"){
			cout << "Whoops, something messed up. Couldn't find the CS in event " << j << endl;
		}
		getline(data,line);
		while (line != "ICS" ){
			numbers = get_pieces(line);
			CS_pt->Fill(atof(numbers[0].c_str()));
			CS_eta->Fill(atof(numbers[5].c_str()));
			CS_m->Fill(atof(numbers[4].c_str()));
			CS_E->Fill(atof(numbers[3].c_str()));
			getline(data,line);
		}
		if (line != "ICS"){
			cout << "Whoops, something messed up. Couldn't find the ICS in event " << j << endl;
		}
		getline(data,line);
		while (line != "CSSK" ){
			numbers = get_pieces(line);
			ICS_pt->Fill(atof(numbers[0].c_str()));
			ICS_eta->Fill(atof(numbers[5].c_str()));
			ICS_m->Fill(atof(numbers[4].c_str()));
			ICS_E->Fill(atof(numbers[3].c_str()));
			getline(data,line);
		}
		if (line != "CSSK"){
			cout << "Whoops, something messed up. Couldn't find the CSSK in event " << j << endl;
		}
		getline(data,line);
		while (line != "ICSSK" ){
			numbers = get_pieces(line);
			CSSK_pt->Fill(atof(numbers[0].c_str()));
			CSSK_eta->Fill(atof(numbers[5].c_str()));
			CSSK_m->Fill(atof(numbers[4].c_str()));
			CSSK_E->Fill(atof(numbers[3].c_str()));
			getline(data,line);
		}
		if (line != "ICSSK"){
			cout << "Whoops, something messed up. Couldn't find the ICSSK in event " << j << endl;
		}
		getline(data,line);
		while (line != "Pass"){
			numbers = get_pieces(line);
			ICSSK_pt->Fill(atof(numbers[0].c_str()));
			ICSSK_eta->Fill(atof(numbers[5].c_str()));
			ICSSK_m->Fill(atof(numbers[4].c_str()));
			ICSSK_E->Fill(atof(numbers[3].c_str()));
			getline(data,line);
		}
	}

	TFile f1("mzero_plots.root","RECREATE");
	f1.mkdir("CS");
	f1.cd("/CS");
	CS_pt->Write();
	CS_eta->Write();
	CS_m->Write();
	CS_E->Write();
	f1.cd("..");
	f1.mkdir("ICS");
	f1.cd("/ICS");
	ICS_pt->Write();
	ICS_eta->Write();
	ICS_m->Write();
	ICS_E->Write();
	f1.cd("..");
	f1.mkdir("CSSK");
	f1.cd("/CSSK");
	CSSK_pt->Write();
	CSSK_eta->Write();
	CSSK_m->Write();
	CSSK_E->Write();
	f1.cd("..");
	f1.mkdir("ICSSK");
	f1.cd("/ICSSK");
	ICSSK_pt->Write();
	ICSSK_eta->Write();
	ICSSK_m->Write();
	ICSSK_E->Write();
}

vector<string> get_pieces(string line){

	int pos = -1;
	int backmark = 0;
	int length = 0;
	string temp;
	vector<string> numbers;
	int pos1;
	numbers.push_back(line.substr(0,line.find(' ')));
	while(line.find(' ',pos+1) != string::npos){
		pos = line.find(' ',backmark);
		length = line.find(' ',pos+1)-pos;
		temp = line.substr(pos+1,length);
		backmark = pos+1;
		numbers.push_back(temp);
	}

	return numbers;

}
