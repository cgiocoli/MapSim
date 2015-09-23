#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

using namespace std;

// conversion: double or int -> string                                   

static const char fINT[] = "%i";
static const char fLONG[] = "%lli";
static const char fDP0[] = "%1.0f";
static const char fDP1[] = "%2.1f";
static const char fDP2[] = "%3.2f";
static const char fDP3[] = "%4.3f";
static const char fDP4[] = "%5.4f";
static const char fDP5[] = "%6.5f";
static const char ee3[] = "%4.3e";

template <class T> string conv (T &val, const char *fact)
{
  char VAL[20]; sprintf (VAL, fact, val);
  return string(VAL);
}

int IntSwap (int i);

int32_t LongSwap (int32_t i);

long long LongLongSwap (long long i);

float FloatSwap(float f);

double DoubleSwap(double d);

struct DATA_subhalo {
  int32_t Ngroups;
  int32_t TotNgroups;
  int32_t Nids;
  int32_t TotNids; int32_t temp;
  int32_t NFiles;
  int32_t Nsubhalos;
  int32_t TotNsubhalos;
};

//istream & operator>>(istream &input, DATA_subhalo &Data) { 
//  input.read((char *)&Data, sizeof(Data));
//  Data.Ngroups = LongSwap(Data.Ngroups);
//  Data.TotNgroups = LongSwap(Data.TotNgroups);
//  Data.Nids = LongSwap(Data.Nids);
//  Data.TotNids = LongLongSwap(Data.TotNids);
//  Data.NFiles = LongSwap(Data.NFiles);
//  Data.Nsubhalos = LongSwap(Data.Nsubhalos);
//  Data.TotNsubhalos = LongSwap(Data.TotNsubhalos);
 
//  return input;
//}

void read_SUBFIND (string &file,
		   vector<double> &mfof,
		   vector<double> &xfof, vector<double> &yfof, vector<double> &zfof,
		   vector<long> &firstsubinfof,
		   vector<long> &nsubinfof,
		   vector<double> &m200, vector<double> &r200,
		   vector<double> &msub,
		   vector<double> &xsub, vector<double> &ysub, vector<double> &zsub,
		   vector<long> &idPH, vector<double> &velDisp, vector<double> &vmax,
		   vector<double> &halfmassradius, vector<double> &rmax, vector<long> &groupN);
