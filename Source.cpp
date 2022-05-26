#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <random>
using namespace std;

void writeTCLContact(const string& tclFileName, const string& trajectoryFileName, 
	const string& contactFileName, const float& cutoff, const float& cutoffLo);
int addCrosslink(const string& dataFileName, const string& tempFileName, const string& contactFileName, const double& percent );
int addCrosslink2(const string& dataFileName, const string& tempFileName, const string& contactFileName, const float& health, int& initialBonds, int& numberOfMolecules);
int measureBeadsPerMolecule(const string& dataFileName, int& numberOfCollagenAtoms, int& numberOfMolecules);
void modifyTheNumberOfBonds(const string& tempFileName, const string& outputDataFileName, const int& realNumberOfBonds);

void main(int argc, char* argv[])
{
	float cutoff = 16.70f;
	float cutoffLo = 0.0f;
	double percent = 0.0;
	int initialBonds, numberOfMolecules;
	float AGEPerMolecule = 500.0f;
	float health;
	string healthStr;
	cout << "Every how many collagen molecules there is one AGE: (from 2 for healthy to 5 for diabetes" << endl;
	cin >> health;
	string cutoffStr, cutoffLoStr, percentStr;
	stringstream ss;
	ss << health;
	ss >> healthStr;
	ss.str("");
	string dataFileName("4-MT-Mineralized-11.5249");
	string trajectoryFileName("dump.equilibrate.lammpstrj");
	string outputDataFileName = dataFileName + "-Crosslinked-Every" + healthStr + "Collagens.dat";
	string tempFileName("temp");
	string tclFileName("contact.tcl");
	string contactFileName("contacts.dat");
	string vmdCommand("vmd -dispdev text -e ");
	dataFileName += ".dat";
	int realNumberOfBonds;
	writeTCLContact(tclFileName, trajectoryFileName, contactFileName, cutoff, cutoffLo);
	vmdCommand = vmdCommand + contactFileName;
	//system(vmdCommand.c_str());

	float healthTarget = health;
	float lastAGEPerMol = 0.0f;
	int counter = 0;
	while (abs(AGEPerMolecule - health) > 0.000001f)
	{
		counter++;
		if (AGEPerMolecule > health) healthTarget -= 0.05f;
		else healthTarget += 0.05f;
		int addedAGE = 0;
		realNumberOfBonds = addCrosslink2(dataFileName, tempFileName, contactFileName, healthTarget, initialBonds, numberOfMolecules);
		modifyTheNumberOfBonds(tempFileName, outputDataFileName, realNumberOfBonds);
		addedAGE = realNumberOfBonds - initialBonds;
		AGEPerMolecule = ((float)numberOfMolecules)/((float)addedAGE);
		cout << "Iteration " << counter << ": Molecule Per AGE is " << AGEPerMolecule << " Health " << healthTarget  << endl;
		if (abs((AGEPerMolecule - lastAGEPerMol)) < 0.01f)
		{
			break;
		}
		if (counter == 6) break;
		else lastAGEPerMol = AGEPerMolecule;
	}


	system("PAUSE");
}

void writeTCLContact(const string& tclFileName, const string& trajectoryFileName, 
	const string& contactFileName, const float& cutoff, const float& cutoffLo)
{
	ofstream tclFile(tclFileName);
	tclFile << "mol load lammpstrj " << trajectoryFileName << endl;
	tclFile << "set contactFile [open " << contactFileName << " w]" << endl;
	tclFile << "set type2 [atomselect top \"type 2\"]" << endl;
	tclFile << "set type3 [atomselect top \"type 3\"]" << endl;
	tclFile << "set contact2 [lindex [measure contacts " << cutoff << " $type2 $type3] 0]" << endl;
	tclFile << "set contact3 [lindex [measure contacts " << cutoff << " $type2 $type3] 1]" << endl;
	tclFile << "set contactSize [llength $contact2]" << endl;
	tclFile << "for {set i 0} {$i < $contactSize} {incr i} {" << endl;
	tclFile << "    set bondLength [measure bond [list [lindex $contact2 $i] [lindex $contact3 $i]]]" << endl;
	tclFile << "    if {$bondLength > " << cutoffLo << " } {" << endl;
	tclFile << "        puts $contactFile \"[lindex $contact2 $i] [lindex $contact3 $i]\"" << endl;
	tclFile << "    }" << endl;
	tclFile << "}" << endl;
	tclFile << "close $contactFile" << endl;
	tclFile << "exit" << endl;
}

int addCrosslink(const string& dataFileName, const string& outputDataFileName, const string& contactFileName, const double& percent)
{
    ifstream contactFile(contactFileName);
	ifstream dataFile(dataFileName);
	ofstream outputDataFile(outputDataFileName);
	string line;
	int numberOfCrosslinks = 0;
	int realNumberOfCrosslinks = 0;
	int numberOfAtoms;
	int numberOfBonds;
	vector<int> first, second;
	vector<bool> atomFlags;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(0, 100);
	vector<double> randomNumbers;
	while (getline(contactFile, line))
	{
		stringstream iss(line);
		int one, two;
		iss >> one;
		iss >> two;
		first.push_back(one);
		second.push_back(two);
		numberOfCrosslinks++;
	}

	for (int i = 0; i < numberOfCrosslinks; i++)
	{
		randomNumbers.push_back(dis(gen));
	}

	int lineNumber = 0;
	int flag = true;
	while (flag)
	{
		getline(dataFile, line);
		lineNumber++;
		int lineSize = line.size();
		if (lineNumber == 3)
		{
			stringstream iss(line);
			int number;
			iss >> number;
			numberOfAtoms = number;
			outputDataFile << number << " atoms" << endl;
			for (int i = 0; i < numberOfAtoms; i++)
				atomFlags.push_back(true);
		}
		else if (lineNumber == 4)
		{
			stringstream iss(line);
			int number;
			iss >> number;
			numberOfBonds = number;
			number += numberOfCrosslinks;
			outputDataFile << number << " bonds" << endl;
		}
		else if (lineNumber == 10)
		{
			stringstream iss(line);
			int number;
			iss >> number;
			number++;
			outputDataFile << number << " bond types" << endl;
		}
		else if (lineSize == 5)
		{
			if (strcmp(line.c_str(), "Bonds") == 0)
			{
				outputDataFile << line << endl;
				getline(dataFile, line); // Skipping the empty line
				outputDataFile << line << endl;
				for (int i = 0; i < numberOfBonds; i++)
				{
					getline(dataFile, line);
					outputDataFile << line << endl;
				}
				for (int i = 0; i < numberOfCrosslinks; i++)
				{

					if (atomFlags[first[i]] && atomFlags[second[i]] && 
						(randomNumbers[i] < 0.86*percent))
					{
						outputDataFile << "  " << numberOfBonds + i + 1 << " 2 " << first[i] + 1 << " " << second[i] + 1 << endl;
						atomFlags[first[i]] = false;
						atomFlags[second[i]] = false;
						realNumberOfCrosslinks++;
					}
				}
				while (getline(dataFile, line))
					outputDataFile << line << endl;
				flag = false;
			}
			else
				outputDataFile << line << endl;
		}
		else
			outputDataFile << line << endl;
	}
	return realNumberOfCrosslinks + numberOfBonds;
}

int addCrosslink2(const string& dataFileName, const string& outputDataFileName, const string& contactFileName, const float& health, int& initialBonds, int& numberOfMolecules)
{
	ifstream contactFile(contactFileName);
	ifstream dataFile(dataFileName);
	int numberOfCollagenAtoms;
	int beadsPerMolecule;
	int healthInt = (int)(100.0f*health);
	beadsPerMolecule = measureBeadsPerMolecule(dataFileName, numberOfCollagenAtoms, numberOfMolecules);
	ofstream outputDataFile(outputDataFileName);
	string line;
	int numberOfCrosslinks = 0;
	int realNumberOfCrosslinks = 0;
	int numberOfBonds;
	float zLow, zHi;
    int numberOfZSections;
	vector<int> first, second;
	vector<int> atomSections;
	vector<float> zSectionLow;
	vector<bool> atomFlags;
	vector<bool> zSectionFlags;
	vector<bool> moleculeFlags;
	for (int i = 0; i < numberOfMolecules; i++)
	{
		if ((100*i) % healthInt < 200)
			moleculeFlags.push_back(true);
		else
			moleculeFlags.push_back(false);
	}

	while (getline(contactFile, line))
	{
		stringstream iss(line);
		int one, two;
		iss >> one;
		iss >> two;
		first.push_back(one);
		second.push_back(two);
		numberOfCrosslinks++;
	}

	int lineNumber = 0;
	int flag = true;
	while (flag)
	{
		getline(dataFile, line);
		lineNumber++;
		int lineSize = line.size();
		if (lineNumber == 10)
		{
			stringstream iss(line);
			int number;
			iss >> number;
			number++;
			outputDataFile << number << " bond types" << endl;
		}
		else if (lineNumber == 4)
		{
			stringstream iss(line);
			int number;
			iss >> number;
			numberOfBonds = number;
			initialBonds = number;
			number += numberOfCrosslinks;
			outputDataFile << number << " bonds" << endl;
		}
		else if (lineNumber == 15)
		{
			numberOfZSections = (int) (numberOfMolecules / health);
			outputDataFile << line << endl;
			stringstream iss("");
			iss << line;
			iss >> zLow;
			iss >> zHi;
			for (int i = 0; i < numberOfZSections; i++)
			{
				float sectionLow = zLow + i*(zHi - zLow) / numberOfZSections;
				zSectionFlags.push_back(true);
				zSectionLow.push_back(sectionLow);
			}
		}
		else if (lineSize == 5)
		{
			if (strcmp(line.c_str(), "Bonds") == 0)
			{
				outputDataFile << line << endl;
				getline(dataFile, line); // Skipping the empty line
				outputDataFile << line << endl;
				for (int i = 0; i < numberOfBonds; i++)
				{
					getline(dataFile, line);
					outputDataFile << line << endl;
				}
				for (int i = 0; i < numberOfCrosslinks; i++)
				{
					int firstMoleculeID  = first[i] / beadsPerMolecule;
					int secondMoleculeID = second[i] / beadsPerMolecule;
					int firstSection  = atomSections[first[i]  - 1];
					int secondSection = atomSections[second[i] - 1];
					/*if (moleculeFlags[firstMoleculeID] &&
						moleculeFlags[secondMoleculeID] &&
						(firstMoleculeID != secondMoleculeID) &&
						zSectionFlags[atomSections[first[i]-1]] &&
						zSectionFlags[atomSections[second[i]-1]] )*/
						if (moleculeFlags[firstMoleculeID] &&
							moleculeFlags[secondMoleculeID] &&
							zSectionFlags[firstSection] &&
							zSectionFlags[secondSection])
					{
						realNumberOfCrosslinks++;
						outputDataFile << "  " << numberOfBonds + realNumberOfCrosslinks << " 2 " << first[i] + 1 << " " << second[i] + 1 << endl;
						moleculeFlags[firstMoleculeID]  = false;
						moleculeFlags[secondMoleculeID] = false;
						zSectionFlags[atomSections[first[i] - 1]] = false;
						zSectionFlags[atomSections[second[i] - 1]] = false;
					}
				}
				while (getline(dataFile, line))
					outputDataFile << line << endl;
				flag = false;
			}
			else if (strcmp(line.c_str(), "Atoms") == 0)
			{
				outputDataFile << line << endl;
				getline(dataFile, line); // Skipping the empty line
				outputDataFile << line << endl;
				for (int i = 0; i < numberOfCollagenAtoms; i++)
				{
					int dumInt;
					float dumFloat;
					float zVal;
					int sectionNumber;
					getline(dataFile, line);
					outputDataFile << line << endl;
					stringstream iss("");
					iss << line;
					iss >> dumInt;
					iss >> dumInt;
					iss >> dumInt;
					iss >> dumFloat;
					iss >> dumFloat;
					iss >> dumFloat;
					iss >> zVal;
					sectionNumber = (int) ((zVal - zLow) / ((zHi - zLow) / numberOfZSections));
					atomSections.push_back(sectionNumber);
				}
			}
			else
				outputDataFile << line << endl;
		}
		else
			outputDataFile << line << endl;
	}
	return realNumberOfCrosslinks + numberOfBonds;
}
int measureBeadsPerMolecule(const string& dataFileName, int& numberOfCollagenAtoms, int& numberOfMolecules)
{
	int beadsPerMolecule = 0;
	int lineNumber = 0;
	int numberOfAtomTypes = 0;
	bool atomsSection = false;
	string line;
	ifstream dataFile(dataFileName);
	numberOfCollagenAtoms = 0;
	while (getline(dataFile, line))
	{
		lineNumber++;
		if (atomsSection)
		{
			stringstream iss(line);
			int atomID, moleculeID, atomType;
			iss >> atomID;
			iss >> moleculeID;
			iss >> atomType;
			iss.str("");
			numberOfCollagenAtoms++;
			if (moleculeID == 1)
			{
				beadsPerMolecule++;
			}
			else if (atomType == numberOfAtomTypes)
			{
				numberOfCollagenAtoms--;
				break;
			}
			else
			{
				numberOfMolecules = moleculeID;
			}
		}
		else if (lineNumber == 9)
		{
			stringstream iss3(line);
			iss3 >> numberOfAtomTypes;
		}
		else if (line.length() >= 5)
			if (!strcmp(line.c_str(), "Atoms"))
			{
				getline(dataFile, line);
				atomsSection = true;
			}
	}
	return beadsPerMolecule;
}

void modifyTheNumberOfBonds(const string& tempFileName, const string& outputDataFileName, const int& realNumberOfBonds)
{
	ifstream tempFile(tempFileName);
	ofstream outputDataFile(outputDataFileName);
	string line;
	float percent = (float)(realNumberOfBonds - 33790)*100.0f / (float)(38108-33790);
	int lineNumber = 0;
	while (getline(tempFile, line))
	{
		int length = line.length();
		lineNumber++;
		if (lineNumber == 4)
			outputDataFile << realNumberOfBonds << " bonds" << endl;
		else if (length >= 11)
		{
			if (line.substr(0, 11) == "Bond Coeffs")
			{
				outputDataFile << line << endl;
				getline(tempFile, line); // Empty line
				lineNumber++;
				outputDataFile << line << endl;
				getline(tempFile, line);
				lineNumber++;
				outputDataFile << "  1 8.565 14.0 62.78 18.20 0      21.00 21.20" << endl
					           << "  2 0.318 16.0  4.28 21.91 10.876 23.71 25.11" << endl;
			}
			else outputDataFile << line << endl;
		}
		else
			outputDataFile << line << endl;
	}
}
