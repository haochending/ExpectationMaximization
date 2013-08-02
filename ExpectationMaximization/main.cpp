//
//  main.cpp
//  ExpectationMaximization
//  CS 486 A4P1
//
//  Created by Haochen Ding on 2013-07-24.
//  Copyright (c) 2013 Haochen Ding. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <time.h>
#include <sys/time.h>


using namespace std;

// structure of the patient record
struct Patient {
    int s;
    int f;
    int deg;
    int t;
    int dun;
};

// structure of CPTs
struct CPT {
    string name;
    int *condition1;
    int *condition2;
    double *probability;
};

double WHITE_NOISE;

// vector of patients data read from trainData.txt
vector<Patient> patients;
// vector of test data read from testData.txt
vector<Patient> tests;
// probability table for all the situations
pair<Patient, double> prTable[48];
// normalized prTable call it weightTable
double weightTable[48];
// table to record the different types sum in the trainData
double sumTable[48];
// the estimate table
int estiTable[48];
// the CPTs
CPT Dunetts, Sloepnea, Foriennditis, Degar, TRIMONO;
// array of all the CPTs
CPT *CPTs;
// array of the results for the trials
double results[20][20];

int readData(string filename) {
    string line;
    ifstream trainData;
    
    trainData.open(filename);
    if (trainData.is_open()) {
        while (trainData.good()) {
            Patient patients_entry;
            
            getline(trainData, line);
            if (line == "") {
                break;
            }
            
            istringstream iss(line);
            int n;
            
            iss >> n;
            patients_entry.s = n;
            
            iss >> n;
            patients_entry.f = n;
            
            iss >> n;
            patients_entry.deg = n;
            
            iss >> n;
            patients_entry.t = n;
            
            iss >> n;
            patients_entry.dun = n;
            
            if (filename == "trainData.txt") {
                patients.push_back(patients_entry);
            } else {
                tests.push_back(patients_entry);
            }
        }
        return 0;
    } else {
        return -1;
    }
}

// set the guess probability for the CPTs
void initCPTs() {
    Dunetts.condition1 = new int[3];
    Dunetts.name = "Dunetts Syndrome";
    Foriennditis.condition1 = new int[3];
    Foriennditis.name = "Foriennditis";
    Degar.condition1 = new int[3];
    Degar.name = "Degar spots";
    
    for (int i = 0; i < 3; i++) {
        Dunetts.condition1[i] = i;
        Foriennditis.condition1[i] = i;
        Degar.condition1[i] = i;
    }
    
    Dunetts.probability = new double[3];
    Dunetts.probability[0] = 0.5;
    Dunetts.probability[1] = 0.25;
    Dunetts.probability[2] = 0.25;
    
    
    Sloepnea.condition1 = new int[6];
    Sloepnea.name = "Sloepnea";
    Sloepnea.condition1[0] = 1;
    Sloepnea.condition1[1] = 1;
    Sloepnea.condition1[2] = 1;
    Sloepnea.condition1[3] = 0;
    Sloepnea.condition1[4] = 0;
    Sloepnea.condition1[5] = 0;
    
    Sloepnea.condition2 = new int[6];
    Sloepnea.condition2[0] = 0;
    Sloepnea.condition2[1] = 1;
    Sloepnea.condition2[2] = 2;
    Sloepnea.condition2[3] = 0;
    Sloepnea.condition2[4] = 1;
    Sloepnea.condition2[5] = 2;
    
    Sloepnea.probability = new double[6];
    Sloepnea.probability[0] = 0.05;
    Sloepnea.probability[1] = 0.1;
    Sloepnea.probability[2] = 0.1;
    Sloepnea.probability[3] = 0.1;
    Sloepnea.probability[4] = 0.5;
    Sloepnea.probability[5] = 0.5;
    
    
    Foriennditis.probability = new double[3];
    Foriennditis.probability[0] = 0.1;
    Foriennditis.probability[1] = 0.8;
    Foriennditis.probability[2] = 0.2;
    
    
    Degar.probability = new double[3];
    Degar.probability[0] = 0.1;
    Degar.probability[1] = 0.2;
    Degar.probability[2] = 0.8;
    
    
    TRIMONO.name = "TRIMONO-HT/S";
    TRIMONO.condition1 = new int[2];
    TRIMONO.condition1[0] = 1;
    TRIMONO.condition1[1] = 0;
    TRIMONO.probability = new double[2];
    TRIMONO.probability[0] = 0.1;
    TRIMONO.probability[1] = 0.9;
    
    CPTs = new CPT[5];
    CPTs[0] = Dunetts;
    CPTs[1] = Sloepnea;
    CPTs[2] = Foriennditis;
    CPTs[3] = Degar;
    CPTs[4] = TRIMONO;
    
}

// print the probability in the 5 CPTs
void printCPTs() {
    for (int i = 0; i < 5; i++) {
        cout << CPTs[i].name << endl;
        
        if (i == 4) {
            cout << CPTs[i].probability[0] << endl;
            cout << CPTs[i].probability[1] << endl;
            cout << "--------------------" << endl;
            continue;
        }
        
        for (int j = 0; j < 3; j++) {
            cout << CPTs[i].probability[j] << endl;
        }
        if (i == 1) {
            for (int k = 3; k < 6; k++) {
                cout << CPTs[i].probability[k] << endl;
            }
        }
        cout << "--------------------" << endl;
    }
}


// update the probability table
void calculatePr() {
    double probability;
    double pr1, pr2, pr3, pr4;
    int index;
    
    for (int a = 0; a < 2; a++) {
        for (int b = 0; b < 2; b++) {
            for (int c = 0; c < 2; c++) {
                for (int d = 0; d < 2; d++) {
                    index = a * 24 + b * 12 + c * 6 + d * 3;
                    
                    for (int e = 0; e < 3; e++) {
                        Patient p;
                        pr1 = (a == 1) ? TRIMONO.probability[0]:TRIMONO.probability[1];
                        if (a == 1) {
                            pr2 = (b == 1) ? Sloepnea.probability[e]:(1 - Sloepnea.probability[e]);
                        } else {
                            pr2 = (b == 1) ? Sloepnea.probability[e + 3]:(1 - Sloepnea.probability[e + 3]);
                        }
                        pr3 = (c == 1) ? Foriennditis.probability[e]:(1 - Foriennditis.probability[e]);
                        pr4 = (d == 1) ? Degar.probability[e]:(1 - Degar.probability[e]);
                        
                        probability = pr1 * pr2 * pr3 * pr4 * Dunetts.probability[e];
                        
                        p.s = b;
                        p.f = c;
                        p.deg = d;
                        p.t = a;
                        p.dun = e;
                        
                        pair<Patient, double> prTable_entry = make_pair(p, probability);
                        prTable[index] = prTable_entry;
                        
                        index ++;
                    }
                    
                    // normalize
                    double sum = prTable[index - 3].second + prTable[index - 2].second + prTable[index - 1].second;
                    if (sum == 0) {
                        weightTable[index - 3] = 0;
                        weightTable[index - 2] = 0;
                        weightTable[index - 1] = 0;
                    } else {
                        weightTable[index - 3] = prTable[index - 3].second/sum;
                        weightTable[index - 2] = prTable[index - 2].second/sum;
                        weightTable[index - 1] = prTable[index - 1].second/sum;
                    }
                }
            }
        }
    }
}

void printPrTable() {
    cout << "Probability Table: " << endl;
    for (int i = 0; i < 48; i++) {
        cout << prTable[i].first.s << " " << prTable[i].first.f << " " << prTable[i].first.deg << " " << prTable[i].first.t << " " << prTable[i].first.dun << " " << prTable[i].second << endl;
    }
    cout << "--------------------" << endl;
    cout << "Weight Table: " << endl;
    for (int i = 0; i < 48; i++) {
        cout << prTable[i].first.s << " " << prTable[i].first.f << " " << prTable[i].first.deg << " " << prTable[i].first.t << " " << prTable[i].first.dun << " " << weightTable[i] << endl;
    }
    cout << "--------------------" << endl;
}

// get the weight of that data entry
double getWeight(Patient p) {
    int counter = 0;
    double weight = 0.0;
    while (counter < 48) {
        if (p.t != prTable[counter].first.t) {
            counter = counter + 24;
            continue;
        }
        
        if (p.s != prTable[counter].first.s) {
            counter = counter + 12;
            continue;
        }
        
        if (p.f != prTable[counter].first.f) {
            counter = counter + 6;
            continue;
        }
        
        if (p.deg != prTable[counter].first.deg) {
            counter = counter + 3;
            continue;
        }
        
        if (p.dun == prTable[counter].first.dun) {
            sumTable[counter] = sumTable[counter] + 1;
            weight = prTable[counter].second;
            break;
        } else {
            if (p.dun == -1) {
                sumTable[counter] = sumTable[counter] + weightTable[counter];
                sumTable[counter + 1] = sumTable[counter + 1] + weightTable[counter + 1];
                sumTable[counter + 2] = sumTable[counter + 2] + weightTable[counter + 2];
                
                weight = prTable[counter].second + prTable[counter + 1].second + prTable[counter + 2].second;
                break;
            } else {
                counter++;
            }
        }
    }
    //    cout << "Weight: " << weight << endl;
    return weight;
}

void printSumTable() {
    cout << "Sum Table: " << endl;
    for (int i = 0; i < 48 ; i++) {
        cout << sumTable[i] << endl;
    }
    cout << "--------------------" << endl;
}

// update the CPTs using the new parameters
void updateCPTs(double totalSum) {
    double sumTri = 0.0;
    double sumSlo[6];
    double sumFor[3];
    double sumDeg[3];
    double sumDun[3];
    double sumTriDun[6];
    
    for (int j = 0; j < 3; j++) {
        sumDun[j] = 0.0;
        sumSlo[j] = 0.0;
        sumFor[j] = 0.0;
        sumDeg[j] = 0.0;
        sumTriDun[j] = 0.0;
    }
    
    for (int k = 0; k < 3; k++) {
        sumSlo[k + 3] = 0.0;
        sumTriDun[k + 3] = 0.0;
    }
    
    for (int i = 0; i < 48; i++) {
        Patient p = prTable[i].first;
        
        // TRIMONO
        if (p.t == 1) sumTri = sumTri + sumTable[i];
        
        // Foriennditis
        if (p.f == 1){
            if (p.dun == 0) {
                sumFor[0] = sumFor[0] + sumTable[i];
            } else if (p.dun == 1) {
                sumFor[1] = sumFor[1] + sumTable[i];
            } else {
                sumFor[2] = sumFor[2] + sumTable[i];
            }
        }
        
        // Degar
        if (p.deg == 1){
            if (p.dun == 0) {
                sumDeg[0] = sumDeg[0] + sumTable[i];
            } else if (p.dun == 1) {
                sumDeg[1] = sumDeg[1] + sumTable[i];
            } else {
                sumDeg[2] = sumDeg[2] + sumTable[i];
            }
        }
        
        // Dunetts
        if (p.dun == 0) {
            sumDun[0] = sumDun[0] + sumTable[i];
        } else if (p.dun == 1) {
            sumDun[1] = sumDun[1] + sumTable[i];
        } else {
            sumDun[2] = sumDun[2] + sumTable[i];
        }
        
        // Sloepnea
        if (p.s == 1) {
            if (p.t == 1) {
                if (p.dun == 0) {
                    sumSlo[0] = sumSlo[0] + sumTable[i];
                } else if (p.dun == 1) {
                    sumSlo[1] = sumSlo[1] + sumTable[i];
                } else {
                    sumSlo[2] = sumSlo[2] + sumTable[i];
                }
            } else {
                if (p.dun == 0) {
                    sumSlo[3] = sumSlo[3] + sumTable[i];
                } else if (p.dun == 1) {
                    sumSlo[4] = sumSlo[4] + sumTable[i];
                } else {
                    sumSlo[5] = sumSlo[5] + sumTable[i];
                }
            }
        }
        
        // Sloepnea 2
        if (p.t == 1) {
            if (p.dun == 0) {
                sumTriDun[0] = sumTriDun[0] + sumTable[i];
            } else if (p.dun == 1) {
                sumTriDun[1] = sumTriDun[1] + sumTable[i];
            } else {
                sumTriDun[2] = sumTriDun[2] + sumTable[i];
            }
            
        } else {
            if (p.dun == 0) {
                sumTriDun[3] = sumTriDun[3] + sumTable[i];
            } else if (p.dun == 1) {
                sumTriDun[4] = sumTriDun[4] + sumTable[i];
            } else {
                sumTriDun[5] = sumTriDun[5] + sumTable[i];
            }
        }
    }
    
    TRIMONO.probability[0] = sumTri/totalSum;
    TRIMONO.probability[1] = 1 - TRIMONO.probability[0];
    
    Dunetts.probability[0] = sumDun[0]/totalSum;
    Dunetts.probability[1] = sumDun[1]/totalSum;
    Dunetts.probability[2] = sumDun[2]/totalSum;
    
    Foriennditis.probability[0] = sumFor[0]/sumDun[0];
    Foriennditis.probability[1] = sumFor[1]/sumDun[1];
    Foriennditis.probability[2] = sumFor[2]/sumDun[2];
    
    Degar.probability[0] = sumDeg[0]/sumDun[0];
    Degar.probability[1] = sumDeg[1]/sumDun[1];
    Degar.probability[2] = sumDeg[2]/sumDun[2];
    
    Sloepnea.probability[0] = sumSlo[0]/sumTriDun[0];
    Sloepnea.probability[1] = sumSlo[1]/sumTriDun[1];
    Sloepnea.probability[2] = sumSlo[2]/sumTriDun[2];
    Sloepnea.probability[3] = sumSlo[3]/sumTriDun[3];
    Sloepnea.probability[4] = sumSlo[4]/sumTriDun[4];
    Sloepnea.probability[5] = sumSlo[5]/sumTriDun[5];
    
}

// estimate the accuracy of the CPTs
double estimate() {
    vector<Patient>::iterator test_it;
    double miss = 0.0;
    double hit = 0.0;
    int counter = 0;
    
    for (int i = 0; i < 48; i++) {
        estiTable[i] = 0;
    }
    
    while (counter < 48) {
        int max = 0;
        if (weightTable[counter] < weightTable[counter + 1]) max = 1;
        if (weightTable[counter + max] < weightTable[counter + 2]) max = 2;
        
        estiTable[counter + max] = 1;
        
        counter = counter + 3;
    }
    
    for (test_it = tests.begin(); test_it != tests.end(); test_it++) {
        int counter = 0;
        while (counter < 48) {
            if (test_it->t != prTable[counter].first.t) {
                counter = counter + 24;
                continue;
            }
            
            if (test_it->s != prTable[counter].first.s) {
                counter = counter + 12;
                continue;
            }
            
            if (test_it->f != prTable[counter].first.f) {
                counter = counter + 6;
                continue;
            }
            
            if (test_it->deg != prTable[counter].first.deg) {
                counter = counter + 3;
                continue;
            }
            
            if (test_it->dun == prTable[counter].first.dun) {
                if (estiTable[counter] == 1) {
                    hit++;
                } else {
                    miss++;
                }
                break;
            } else {
                counter++;
            }
        }
    }
    return hit/(hit + miss);
}

// generate the double random number
double getRandom() {
    return WHITE_NOISE * ((double)rand() / (double)RAND_MAX);
}

// add white noise to the parameters
void whiteNoise() {
    double dunettsSum = 0.0;
    double temp;
    
    TRIMONO.probability[0] = TRIMONO.probability[0] + getRandom();
    TRIMONO.probability[1] = TRIMONO.probability[1] + getRandom();
    TRIMONO.probability[0] = TRIMONO.probability[0]/(TRIMONO.probability[0] + TRIMONO.probability[1]);
    TRIMONO.probability[1] = 1 - TRIMONO.probability[0];
    
    for (int i = 0; i < 3; i++) {
        Dunetts.probability[i] = Dunetts.probability[i] + getRandom();
        dunettsSum = dunettsSum + Dunetts.probability[i];
        
        temp = Sloepnea.probability[i] + getRandom();
        Sloepnea.probability[i] = temp/(temp + ((1 - Sloepnea.probability[i]) + getRandom()));
        
        temp = Foriennditis.probability[i] + getRandom();
        Foriennditis.probability[i] = temp/(temp + ((1 - Foriennditis.probability[i]) + getRandom()));
        
        temp = Degar.probability[i] + getRandom();
        Degar.probability[i] = temp/(temp + ((1 - Degar.probability[i]) + getRandom()));
    }
    
    for (int j = 3; j < 6; j++) {
        temp = Sloepnea.probability[j] + getRandom();
        Sloepnea.probability[j] = temp/(temp + ((1 - Sloepnea.probability[j]) + getRandom()));
    }
    
    Dunetts.probability[0] = Dunetts.probability[0]/dunettsSum;
    Dunetts.probability[1] = Dunetts.probability[0]/dunettsSum;
    Dunetts.probability[2] = Dunetts.probability[0]/dunettsSum;
    
    
}


int main(int argc, const char * argv[])
{
    double oldTotalWeight = 0.0;
    double totalWeight = 0.0;
    double accuracy, old_accuracy;
    
    srand((unsigned)time(NULL));
    
    int result = readData("trainData.txt");
    
    if (result < 0) {
        cerr << "Read trainData error !!!" << endl;
        exit(1);
    }
    
    result = readData("testData.txt");
    
    if (result < 0) {
        cerr << "Read testData error !!!" << endl;
        exit(1);
    }
    
    // 20 trials
    for (int k = 0; k < 20; k++) {
        WHITE_NOISE = 0.2 * (k + 1);
        
        // run 20 times for mean value and standard error
        for (int j = 0; j < 20; j++) {
            
            int counter = 0;
            
            initCPTs();
            whiteNoise();
            
            
            while (true) {
                oldTotalWeight = totalWeight;
                totalWeight = 0.0;
                
                printCPTs();
                
                calculatePr();
                
                printPrTable();
                
                // initialize the sumTable
                for (int i = 0; i < 48; i++) {
                    sumTable[i] = 0.0;
                }
                
                // update the sum table and get the total sum
                vector<Patient>::iterator it;
                for (it = patients.begin(); it != patients.end(); it++) {
                    double weight = getWeight(*it);
                    totalWeight = totalWeight + weight;
                }
                
                printSumTable();
                
                cout << "TotalWeight: " << totalWeight << endl;
                
                if (counter == 0) old_accuracy = estimate();
                
                updateCPTs(2000);
                
                if (((totalWeight - oldTotalWeight) < 0.01 && (totalWeight - oldTotalWeight) > 0) || ((totalWeight - oldTotalWeight) > -0.01 && (totalWeight - oldTotalWeight) < 0)) break;
                
                counter++;
            }
            
            cout << "Iteration: " << counter << endl;
            
            accuracy = estimate();
            cout << "First Guess: " << old_accuracy << endl;
            cout << "Final Accuracy: " << accuracy << endl;
            results[k][j] = old_accuracy;
        }
    }
    
    // print the final result
    for (int m = 0; m < 20; m++) {
        for (int n = 0; n < 20; n++) {
            if (n == 19) {
                cout << results[m][n] << endl;
            } else {
                cout << results[m][n] << "\t";
            }
        }
    }
    
    return 0;
}

