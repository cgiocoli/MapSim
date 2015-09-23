#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>     
#include <stdlib.h>     
#include <ctime>
#include <CCfits/CCfits>
#include <readSUBFIND.h>
#ifdef COSMOLIB
#include <cosmo.h>
#endif

/*****************************************************************************/
/*                                                                           */
/*    This code has been developed to create maps from files of              */
/*      numerical simuations. Up to now it has been optimized to run on      */
/*      "CoDECS-like" simulations and to read gadget1 format files           */
/*                          - it runs on multiple snapshots                  */
/*                          - it reads the SUBFIND and FOF catalogues        */
/*                                                                           */
/*     it returns a list of .fits file of the 2D mass map in each plane      */
/*                          subfindinfield                                   */
/*                          fofinfield                                       */
/*       giving the positions et al. of the fof and subs present in the cone */
/*                                                                           */
/*                                                                           */
/*                              dev. by Carlo Giocoli - cgiocoli@gmail.com   */
/*****************************************************************************/

using namespace std;
using namespace CCfits;

const int bleft = 24;
const double speedcunit = 2.99792458e+3;

template <class T>
int locate (const std::vector<T> &v, const T x){
  size_t n = v.size ();
  int jl = -1;
  int ju = n;
  bool as = (v[n-1] >= v[0]);
  while (ju-jl > 1){
    int jm = (ju+jl)/2;
    if ((x >= v[jm]) == as)
      jl=jm;
    else
      ju=jm;
  }
  if (x == v[0])
    return 0;
  else if (x == v[n-1])
    return n-2;
  else
    return jl;
}

double getY(std:: vector<double> x, std:: vector<double> y,double xi){
  int nn = x.size();
  if(x[0]<x[nn-1]){         
    if(xi>x[nn-1]) return y[nn-1];
    if(xi<x[0]) return y[0];
  }  
  else{
    if(xi<x[nn-1]) return y[nn-1];
    if(xi>x[0]) return y[0];  
  }
  int i = locate (x,xi);
  i = std::min (std::max (i,0), int (nn)-2);
  double f=(xi-x[i])/(x[i+1]-x[i]);
  if(i>1 && i<nn-2){
    double a0,a1,a2,a3,f2;                                     
    f2 = f*f;
    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];
    a1 = y[i-1] - y[i] - a0;
    a2 = y[i+1] - y[i-1];
    a3 = y[i];
    return a0*f*f2+a1*f2+a2*f+a3;
  }                                                                      
  else return f*y[i+1]+(1-f)*y[i];           
} 

extern float weight(float, float, double);
// weight function for each pixel
float weight (float ixx, float ixh, double dx) {
  float DD = ixx-ixh;
  float x=fabs(DD)/dx;
  float w;    
  
   if(fabs(DD)<=0.5*dx)
     w=3./4.-x*x;
   else if(fabs(DD)>0.5*dx && fabs(DD)<=0.5*3.0*dx)
     w=0.5*((3./2.-x)*(3./2.-x));
   else w=0.;
   return w;
}

// grid points distribution function
valarray<float> gridist(vector<float> x, vector<float> y, int nn){
  std::valarray<float> grxy( nn*nn );  
  int n0 = x.size();
  // cout << "     gridist: check number of particles: = " << n0 << endl;
  // cout << " " << endl;
  // grid size
  double dl = 1/double(nn);    
  // --- - - - - - - - - - - - - - - - - - - - - - - ---
  //                               _ _ _ _ _ _
  //  The order of the points is: |_7_|_8_|_9_|
  //                              |_4_|_5_|_6_|
  //                              |_1_|_2_|_3_|
  // coordinate between 0 and 1 and mass particle = 1
  //
  // --- - - - - - - - - - - - - - - - - - - - - - - ---
  for (int i=0; i<n0; i++){
    int grx=int(x[i]/dl)+1;
    int gry=int(y[i]/dl)+1;

    //float xir = x[i]/dl;
    //float yir = y[i]/dl;
    
    //float xh = fabs(xir - float(grx));
    //float yh = fabs(yir - float(gry));
    
    int   gridpointx[9], gridpointy[9];
    float posgridx[9], posgridy[9];
    float wfx[9], wfy[9];
    
    gridpointx[0] = grx;
    gridpointy[0] = gry;
    
    for (int j=0; j<9; j++){
      gridpointx[j]=gridpointx[0]+(j%3)-1;
      gridpointy[j]=gridpointy[0]+(j/3)-1;
            
      posgridx[j]=(gridpointx[j]+0.5)*dl;
      posgridy[j]=(gridpointy[j]+0.5)*dl;

      wfx[j] = weight(x[i],posgridx[j],dl);
      wfy[j] = weight(y[i],posgridy[j],dl);

      int grxc = gridpointx[j];
      int gryc = gridpointy[j];
      
      if(grxc>=0 && grxc<nn && gryc>=0 && gryc<nn) 
	grxy[grxc+nn*gryc] = grxy[grxc+nn*gryc] + float(wfx[j])*float(wfy[j]); 
    }
  } 
  return grxy;
}

struct DATA
{
  // long substituted with int32_t
  int32_t npart[6];	
  double massarr[6];
  double time;
  double redshift;
  int32_t flag_sfr;
  int32_t flag_feedback;
  int32_t npartTotal[6];	
  int32_t flag_cooling;
  double boxsize;
  double om0;
  double oml;
  double h;
  int32_t la[bleft];
};

istream & operator>>(istream &input, DATA &Data) 
{ 
  input.read((char *)&Data, sizeof(Data));
#ifdef SWAP
  for (int i=0; i<6; i++) Data.npart[i] = LongSwap(Data.npart[i]); 
  for (int i=0; i<6; i++) Data.massarr[i] = DoubleSwap(Data.massarr[i]); 
  Data.time = DoubleSwap(Data.time);
  Data.redshift = DoubleSwap(Data.redshift);
  Data.flag_sfr = LongSwap(Data.flag_sfr);
  Data.flag_feedback = LongSwap(Data.flag_feedback);
  for (int i=0; i<6; i++) Data.npartTotal[i] = LongSwap(Data.npartTotal[i]); 
  Data.flag_cooling = LongSwap(Data.flag_cooling);
  Data.boxsize = DoubleSwap(Data.boxsize);
  Data.om0 = DoubleSwap(Data.om0);
  Data.oml = DoubleSwap(Data.oml);
  Data.h = DoubleSwap(Data.h);
  for (int i=0; i<bleft; i++) Data.la[i] = IntSwap(Data.la[i]); 
#endif
  return input;
};

void readParameters(int *npix, double *boxl,
		    double *zs, double *fov, string *filredshiftlist,
		    string *filsnaplist, string *pathsnap,
		    string *idc, 
		    int *seedcenter, int *seedface, int *seedsign,
		    string *subfiles, string *simulation, int *nfiles, string *partinplanes){ 

  string butstr;
  ifstream inputf;
  inputf.open("INPUT");
  if(inputf.is_open()){
    inputf >> butstr; // number of pixels
    inputf >> *npix;
    inputf >> butstr; // boxl
    inputf >> *boxl;
    inputf >> butstr; // source redshift
    inputf >> *zs;
    inputf >> butstr; // field of view in degrees
    inputf >> *fov;
    inputf >> butstr; // file with the redshift list it may contain three columns: snap 1/(1+z) z
    inputf >> *filredshiftlist;
    inputf >> butstr; // file with the snapshot list available 
    inputf >> *filsnaplist;
    inputf >> butstr; // path where the snaphosts are located
    inputf >> *pathsnap;
    inputf >> butstr; // simulation name (prefix infront at the snap file)
    inputf >> *simulation;
    inputf >> butstr;  // number of files per snapshot
    inputf >> *nfiles;
    inputf >> butstr; // path and file name of the comoving distance file (if not available you may use CosmoLib)
    inputf >> *idc;
    inputf >> butstr; // path of SUBFIND if NO I will not look for them
    inputf >> *subfiles;
    inputf >> butstr; // seed for the random location of the center
    inputf >> *seedcenter;
    inputf >> butstr; // seed for the random selection of the dice face
    inputf >> *seedface;
    inputf >> butstr; // seed for the selection of the sign of the coordinates
    inputf >> *seedsign;  
    inputf >> butstr; // which particles in the planes (ALL: one file for all, or NO: one for each)
    inputf >> *partinplanes;  
    inputf.close();
  }else{
    cout << " INPUT file does not exsit ... I will stop here!!! " << endl;
    exit(1);
  }
};

void getPolar(double x, double y, double z, double *ra, double *dec, double *d){
  *d = sqrt(x*x+y*y+z*z);
  *dec = asin(x/(*d));
  *ra = atan2(y,z);
}

// if could read a cone file created by Mauro Roncarelli for the randomization 
void readCone(string fname, double *boxl, double *buta, double *butb, double *zs,double *Ds,double *butc,
	      int *ncubes, int *nplanes,vector<int> *ibox,vector<int> *isnap,vector<double> *iDfinal,
	      vector<string> *Reflection1,vector<string> *Reflection2,vector<string> *Reflection3, 
	      vector<string> *Axes1, vector<string> *Axes2, vector<string> *Axes3,
	      vector<double> *xM0, vector<double> *yM0, vector<double> *zM0){

  ifstream ifile;
  ifile.open(fname.c_str());
  string buts;
  ifile >> buts; ifile >> buts;
  ifile >> *boxl;
  ifile >> buts; ifile >> buts;
  ifile >> *buta;
  ifile >> buts; ifile >> buts;
  ifile >> *butb;
  ifile >> buts; ifile >> buts;
  ifile >> *zs;
  ifile >> buts; ifile >> buts;
  ifile >> *Ds;
  ifile >> buts; ifile >> buts;
  ifile >> *butc;
  ifile >> buts; ifile >> buts;
  ifile >> *ncubes;
  ifile >> buts; ifile >> buts;
  ifile >> *nplanes;
  ifile >> buts;
  for(int i=0;i < *nplanes;i++){
    int a,b;
    double c;
    ifile >> a >> b >> c;
    ibox->push_back(a);
    isnap->push_back(b);
    iDfinal->push_back(c);
  }
  ifile >> buts;
  for(int i=0;i<*ncubes;i++){
    string a,b,c,d,e,f;
    double x0,y0,z0;
    double butr;
    ifile >> a >> b >> c >> d >> e >> f >> x0 >> y0 >> z0 >> butr;
    Reflection1->push_back(a);
    Reflection2->push_back(b);
    Reflection3->push_back(c);
    Axes1->push_back(d);
    Axes2->push_back(e);
    Axes3->push_back(f);
    xM0->push_back(x0);
    yM0->push_back(y0);
    zM0->push_back(z0);
  }
  ifile.close();
}

int main(){ 
  // to set in the parameter file
  int noSNAP=1;
  cout << "   ------------------------------------------------------ " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -           2D Mapping Simulation Snapshot           - " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -               collapsing one dimension             - " << endl;
  cout << "   ------------------------------------------------------ " << endl;
  // ... masses of the different type of particles
  double m0,m1,m2,m3,m4,m5;
  // ******************** to be read in the INPUT file ********************
  // ... project - number of pixels
  int npix;
  // ... set by hand the redshift of the source and than loop only up to 
  // ... the needed snaphost when creating the light cone!
  double zs,Ds;
  string filredshiftlist,filsnaplist;
  // ... loop on different snapshots
  string pathsnap; // = "/dati1/cgiocoli/CoDECS/";
  double boxl; // Mpc/h
  string idc; // comoving distance file
  int seedcenter, seedface, seedsign,nfiles;
  double fov;
  string subfiles;
  string simulation; // Simulation Name  
  string partinplanes; // ALL all part in one plane, NO each part type in different planes

  readParameters(&npix,&boxl,&zs,&fov,
		 &filredshiftlist,&filsnaplist,
		 &pathsnap,&idc,&seedcenter,
		 &seedface,&seedsign,
		 &subfiles,&simulation,&nfiles,&partinplanes);

  // check the existence of a file .cone (from Mauro Roncarelli for the randmization)
  int yes = system("ls -1 *.cone > cone.file");
  ifstream ifcone;
  ifcone.open("cone.file");
  string fname="NO";
  ifcone >> fname;
  ifcone.close();
  int nplanes,ncubes;
  vector<int> isnap,ibox;
  vector<double> iDfinal;
  vector<string> Reflection1,Axes1;
  vector<string> Reflection2,Axes2;
  vector<string> Reflection3,Axes3;
  vector<double> xM0,yM0,zM0;
  if(fname!="NO"){
    double buta,butb,butc,butd;
    readCone(fname,&boxl,&buta,&butb,&zs,&Ds,&butc,&ncubes,&nplanes,&ibox,&isnap,&iDfinal,
    	     &Reflection1,&Reflection2,&Reflection3,&Axes1,&Axes2,&Axes3,
	     &xM0,&yM0,&zM0);
    /*
      for(int i=0;i < nplanes;i++){
      cout << ibox[i] << "  " << isnap[i] << "  " << iDfinal[i] << endl;
      }
    */
    for(int i=0;i < ncubes;i++){
      cout << xM0[i] << "  " << yM0[i] << "  " << zM0[i] << endl;
      cout << Reflection1[i] << "  " << Reflection2[i] << "  " << Reflection3[i] << endl;
      cout << Axes1[i] << "  " << Axes2[i] << "  " << Axes3[i] << endl;
    }
  }
  double Omega_matter,Omega_lambda,Omega_baryon,hubble;
#ifdef COSMOLIB // if run with COSMOLIB on, you need a file cormo.par!
  /// read a file with the cosmological parameters
  /// and update the comoving distance table
  string fcosmo;
  ifstream ifcosmo;
  ifcosmo.open("cosmo.par");
  if(ifcosmo.is_open()){
    string bstring;
    ifcosmo >> bstring;
    ifcosmo >> Omega_matter;
    ifcosmo >> bstring;
    ifcosmo >> Omega_lambda;
    ifcosmo >> bstring;
    ifcosmo >> Omega_baryon;
    ifcosmo >> bstring;
    ifcosmo >> hubble;
  }else{
    cout << " you activate COSMOLIB however to run I need " << endl;
    cout << " an input file called cosmo.par containing " <<  endl;
    cout << " the cosmological parameters " << endl;
    cout << " check it out ... I will STOP here!!! " << endl;
    exit(1);
  }
  cout << " Assuming the following cosmology:" << Omega_matter << "  " << Omega_lambda << "  " << Omega_baryon << "  " << hubble << endl;
  COSMOLOGY cosmo;
  cosmo.setOmega_matter(Omega_matter, true);
  cosmo.setOmega_lambda(Omega_lambda);
  cosmo.setOmega_baryon(Omega_baryon);
  cosmo.sethubble(hubble);
#endif
  string snpix = conv(npix,fINT);
  // ... read the redshift list and the snap_available
  // ... to build up the light-cone
  ifstream redlist;
  redlist.open(filredshiftlist.c_str());
  vector <int> tsnaplist;
  vector <double> dtsnaplist;
  vector <double> tredlist;
  int nmax = 1024;
  vector<double> snapToRedshift(nmax); // no way to have more than 1024 snaphosts
  for(int i=0;i<nmax;i++){
    snapToRedshift[i] = -1;
  }
  if(redlist.is_open()){
    int buta;
    double butb,butc;
    while(redlist >> buta >> butb >> butc){
      tsnaplist.push_back(buta);
      tredlist.push_back(butc);
      dtsnaplist.push_back(buta);
      if(buta>nmax){
	cout << " check nmax variable and increase it! " << endl;
	cout << " for now I will STOP here!!! " << endl;
	exit(1);
      }
      snapToRedshift[buta] = butc;
    }
  }else{
    cout << " redshift list file redshift_list.txt does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }
  
  ifstream snaplist;
  snaplist.open(filsnaplist.c_str());
  vector <int> lsnap;
  vector<double> lred;
  vector<double> lD;
  if(snaplist.is_open()){
    int s;
    double zn;
    int check=0;
    while(snaplist >> s){
      zn = getY(dtsnaplist,tredlist,double(s));
      if(zn>zs) check=1;
      if(check==0){
	lsnap.push_back(s);
	lred.push_back(zn);
      }
    }
  }else{
    cout << " snaplist.txt does not exist in the code dir " << endl;
    cout << "   I will STOP here!!! " << endl;
    exit(1);
  }
  snaplist.close();
  int nsnaps = lsnap.size();

  cout << "  " << endl;
  cout << " opening path for snapshots: " << endl;
  cout << pathsnap << endl;
  cout << " " << endl;
  cout << " I will look for comoving distance " << endl;
  cout << "      file = " << idc << endl;
  cout << " " << endl;
  
  ifstream infiledc;
  vector<double> zl, dl;
#ifdef COSMOLIB
  double zi=0;
  double dli;
  while(zi<zs+1){
    dli = cosmo.angDist(0.,zi)*(1+zi)*hubble;  // we want it in comoving Mpc/h
    dl.push_back(dli);
    zl.push_back(zi);
    cout << zi << "  " << dli << endl;
    zi+=0.001;
  }
#else
  infiledc.open(idc.c_str());
  if(infiledc.is_open()){
    double zi,dli;
    while(infiledc >> zi >> dli){
      zl.push_back(zi);
      dl.push_back(dli*speedcunit);
      cout << zi << "  " << dli*speedcunit << endl;
    }
    infiledc.close();
  }
  else{
    cout << "  " << endl;
    cout << " the comoving distance file: " << idc << endl;
    cout << " does not exists " << endl;
    cout << " I will STOP here!!! " << endl;
    exit(1);
  }
  
  if(zs>zl[zl.size()-1]){
    cout << " source redshift larger than the highest available redshift in the comoving distance file " << endl;
    cout << "  that is = " << zl[zl.size()-1] << endl;
    cout << " I will STOP here !!! " << endl;
    exit(1);
  }
#endif
  int nreplications;
  if(fname=="NO"){
    // cout << " Ds = " << Ds << endl;
    Ds = getY(zl,dl,zs);  // comoving distance of the last plane
    // cout << " Ds = " << Ds << endl;
    // exit(1);
    nreplications = int(Ds/boxl)+1;
  }else{
    nreplications=ncubes;
    // cout << " Ds = " << Ds << endl;
    Ds = getY(zl,dl,zs);  // comoving distance of the last plane
    // cout << " Ds = " << Ds << endl;
    // exit(1);
  }
  cout << "  nreplications " << nreplications << "  " << Ds << "  " << std:: endl;
  vector<int> replication;
  vector<int> fromsnap;
  vector<double> lD2;
  if(fname=="NO"){
    for(int j=0;j<nreplications;j++){
      for(int i=0;i<nsnaps;i++){
	double ldbut = getY(zl,dl,lred[i]);
	if(ldbut>=j*boxl && ldbut<=(j+1.)*boxl){
	  std:: cout << " simulation snapshots = " << ldbut << "  " << lred[i] << "  " << j << " from snap " 
		     << lsnap[i] << "  " << boxl*j << "  " << boxl*(j+1) 
		     << std:: endl;
	  replication.push_back(j);
	  fromsnap.push_back(lsnap[i]);
	  lD.push_back(ldbut);
	}
      }
    }
    cout << "  " << endl;
    cout << "  reorganazing the planes " << std:: endl;
    cout << " " << endl;

    for(int i=0;i<lD.size();i++){
      if(i<(lD.size()-1)) lD2.push_back(lD[i+1]);
      else lD2.push_back(Ds);
    }
    for(int i=0;i<lD.size();i++){
      for(int k=1;k<=512;k++){
	if(lD[i]<double(k)*boxl && lD2[i]>double(k)*boxl){
	  lD[i+1] = lD2[i]; 
	  lD2[i]=double(k)*boxl;
	}
      }
      if(lD[i]<513*boxl && lD2[i]>513*boxl){
	cout << " exiting ... increase the number of replications by hand in the file it is now 512 !!!! " << endl;
	exit(1);
      }
    } 
  }
  std:: cout << " Comoving Distance of the last plane " << Ds << std:: endl;
  std:: cout << "  " << endl;
  vector<double> zsimlens(nsnaps);
  cout << " nsnaps = " << nsnaps << endl;
  if(fname=="NO"){
    for(int i=0;i<nsnaps;i++){
      if(i<nsnaps-1){
	if(lD[i+1]-lD2[i]>boxl*1e-9){
	  fromsnap[i] = -fromsnap[i];
	}
      }
      double dlbut = (lD[i] + lD2[i])*0.5;    
      // half distance between the two!
      zsimlens[i] = getY(dl,zl,dlbut);
      std:: cout << zsimlens[i] << " planes = " << lD[i] << "  " << lD2[i] << "  " << replication[i] 
		 << " from snap " << fromsnap[i] << std:: endl;
    }
  }
  vector<double> bfromsnap,blD,blD2,bzsimlens,blred;
  vector<int> breplication, blsnap;
  int pl=0;
  vector<int> pll;
  if(fname=="NO"){
    for(int i=0;i<nsnaps;i++){
      bfromsnap.push_back(fabs(fromsnap[i]));
      blD.push_back(lD[i]);
      blD2.push_back(lD2[i]);
      bzsimlens.push_back(zsimlens[i]);
      breplication.push_back(replication[i]);
      blsnap.push_back(lsnap[i]);
      blred.push_back(lred[i]);
      if(fromsnap[i]<=0){
	bfromsnap.push_back(-fromsnap[i]);
	blD.push_back(lD2[i]);
	blD2.push_back(lD[i+1]);
	double dlbut = (lD[i+1] + lD2[i])*0.5;    
	// half distance between the two!
	bzsimlens.push_back(getY(dl,zl,dlbut));
	breplication.push_back(replication[i+1]);
	blsnap.push_back(lsnap[i]);
	blred.push_back(lred[i]);
      }   
    }
    cout << "  " << endl;
    cout << "  re-reorganazing the planes " << std:: endl;
    cout << " " << endl; 
    // here I have to re-set what has been read for the .cone file!
  }
  ofstream planelist;
  planelist.open("planes_list.txt");
  if(fname=="NO"){
    nsnaps = bfromsnap.size();  
    for(int i=0;i<nsnaps;i++){
      std:: cout << bzsimlens[i] << " bplanes = " << blD[i] << "  " << blD2[i] << "  " << breplication[i] << " from snap " << bfromsnap[i] << std:: endl;
      pl++;
      planelist <<  pl << "   " << bzsimlens[i] << "   " << blD[i] << "   " << blD2[i] << "   " << breplication[i] << "   " << bfromsnap[i] << "   " << blred[i] << std:: endl;
      pll.push_back(pl);
    }
    
    if(blD2[nsnaps-1]<Ds){
      // we need to add one more snaphost
      snaplist.open(filsnaplist.c_str());
      double s;
      while(snaplist >> s){
	if(s<bfromsnap[nsnaps-1]){
	  bfromsnap.push_back(s);
	  double dlbut = (Ds+blD2[nsnaps-1])*0.5;
	  double zbut = getY(dl,zl,dlbut);
	  blD.push_back(blD2[nsnaps-1]);
	  blD2.push_back(Ds);
	  bzsimlens.push_back(zbut);
	  breplication.push_back(breplication[nsnaps-1]+1);
	  double zn = getY(dtsnaplist,tredlist,double(s));
	  blsnap.push_back(s);
	  lred.push_back(zn);
	  nsnaps++;
	  snaplist.close();	
	}
      }
      pl++;
      std:: cout << bzsimlens[nsnaps-1] << " bplanes = " << blD[nsnaps-1] << "  " << blD2[nsnaps-1] << "  " << breplication[nsnaps-1] << " from snap " << bfromsnap[nsnaps-1] << std:: endl;	
      planelist << pl << "   " << bzsimlens[nsnaps-1] << "   " << blD[nsnaps-1] << "   " << blD2[nsnaps-1] << "   " << breplication[nsnaps-1] << "   " << bfromsnap[nsnaps-1] << "   " << blred[nsnaps-1] << std:: endl;
      pll.push_back(pl);	 
    } 
  }else{
    // need the redefine the variables considering Mauro's file
    nsnaps = nplanes;
    bzsimlens.resize(nsnaps);
    blD.resize(nsnaps);
    blD2.resize(nsnaps);
    breplication.resize(nsnaps);
    bfromsnap.resize(nsnaps);
    blred.resize(nsnaps);
    blD[0] = 0;
    blsnap.resize(nsnaps);
    for(int i=0;i<nsnaps-1;i++){
      blD[i+1] = iDfinal[i];
    }
    for(int i=0;i<nsnaps;i++){
      blD2[i] = iDfinal[i];
      breplication[i]=ibox[i]-1;
      bfromsnap[i] = isnap[i];
      blred[i] = snapToRedshift[isnap[i]];
      blsnap[i] = isnap[i];
      pll.push_back(i+1);
      if(blred[i]<-0.5){
	cout << " error in the snapToRedshift mapping ... check the variable " << endl;
	cout << " I will STOP here!!! " << endl;
	exit(1);
      }
      double dlbut = (blD[i]+blD2[i])*0.5;
      double zbut = getY(dl,zl,dlbut);
      bzsimlens[i] = zbut;
      planelist << i+1 << "   " << bzsimlens[i] << "   " << blD[i] << "   " << blD2[i] << "   " << breplication[i] << "   " << bfromsnap[i] << "   " << blred[i] << std:: endl;
    }
  }
  planelist.close();

  std:: cout << "  " << endl; 
  // randomizzation of the box realizations :
  int nrandom = breplication[nsnaps-1]+1;
  vector<double> x0(nrandom), y0(nrandom), z0(nrandom); // ramdomizing the center of the simulation [0,1]
  vector<int> face(nrandom); // face of the dice
  vector<int> sgnX(nrandom), sgnY(nrandom),sgnZ(nrandom);
  if(fname=="NO"){
    for(int i=0;i<nrandom;i++){
      srand(seedcenter+i*13);
      x0[i] = rand() / float(RAND_MAX);
      y0[i] = rand() / float(RAND_MAX);
      z0[i] = rand() / float(RAND_MAX);
      cout << "  " << endl;
      cout << " random centers  for the box " << i << " = " << x0[i] << "  " << y0[i] << "  " << z0[i] << endl;
      face[i] = 7;
      srand(seedface+i*5);
      while(face[i]>6 || face[i]<1) face[i] = int(1+rand() / float(RAND_MAX)*5.+0.5);
      std:: cout << " face of the dice " << face[i] << std:: endl;
      sgnX[i] = 2;
      srand(seedsign+i*8);
      while(sgnX[i] > 1 || sgnX[i] < 0) sgnX[i] = int(rand() / float(RAND_MAX)+0.5);
      sgnY[i] = 2;
      while(sgnY[i] > 1 || sgnY[i] < 0) sgnY[i] = int(rand() / float(RAND_MAX)+0.5);
      sgnZ[i] = 2;
      while(sgnZ[i] > 1 || sgnZ[i] < 0) sgnZ[i] = int(rand() / float(RAND_MAX)+0.5);
      if(sgnX[i]==0) sgnX[i]=-1;
      if(sgnY[i]==0) sgnY[i]=-1;
      if(sgnZ[i]==0) sgnZ[i]=-1;
      std:: cout << " signs of the coordinates = " << sgnX[i] << "  " << sgnY[i] << " " << sgnZ[i] << endl;
    }
  }else{
    for(int i=0;i<nrandom;i++){
      // reflection or not reflection
      sgnX[i] = 1;
      sgnY[i] = 1;
      sgnZ[i] = 1;
      if(Reflection1[i]=="T") sgnX[i] = -1;
      if(Reflection2[i]=="T") sgnY[i] = -1;
      if(Reflection3[i]=="T") sgnZ[i] = -1;
      std:: cout << " signs of the coordinates = " << sgnX[i] << "  " << sgnY[i] << " " << sgnZ[i] << endl;
      // choose the face!
      if(Axes1[i]=="x" && Axes2[i]=="y" && Axes3[i]=="z") face[i] = 1;
      if(Axes1[i]=="x" && Axes2[i]=="z" && Axes3[i]=="y") face[i] = 2;
      if(Axes1[i]=="y" && Axes2[i]=="z" && Axes3[i]=="x") face[i] = 3;
      if(Axes1[i]=="y" && Axes2[i]=="x" && Axes3[i]=="z") face[i] = 4;
      if(Axes1[i]=="z" && Axes2[i]=="x" && Axes3[i]=="y") face[i] = 5;
      if(Axes1[i]=="z" && Axes2[i]=="y" && Axes3[i]=="x") face[i] = 6;
      std:: cout << " face of the dice " << face[i] << std:: endl;
      // last thing to do!
      x0[i] = xM0[i]/boxl;
      y0[i] = yM0[i]/boxl;
      z0[i] = zM0[i]/boxl;
      cout << "  " << endl;
      cout << " random centers  for the box " << i << " = " << x0[i] << "  " << y0[i] << "  " << z0[i] << endl;
    }
  }
  std:: cout << "  " << endl;
  cout << "  " << endl;
  cout << " set the field of view to be square in degrees " << endl;
  double h0,fovradiants;
  double om0, omL0;
  fovradiants = fov/180.*M_PI;
  if(fovradiants*Ds>boxl){
    std:: cout << " field view too large ... I will STOP here!!! " << std:: endl;
    std:: cout << " maximum value allowed " << boxl/Ds*180./M_PI << " in degrees " << std:: endl;
    exit(1);
  }
  cout << "  " << endl;

  // loop on snapshots ****************************************
  cout << " now loop on " << nsnaps << " snapshots " << endl;
  cout << "  " << endl;
  for(int nsnap=0;nsnap<nsnaps;nsnap++){
    if(blD2[nsnap]-blD[nsnap] < 0){
      cout << " comoving distanze of the starting point " << blD[nsnap] << endl;
      cout << " comoving distanze of the final    point " << blD2[nsnap] << endl;
      cout << " please check this out! I will STOP here!!! " << endl;
      exit(1);
    }
    int rcase = breplication[nsnap];
    std::valarray<float> mapxytot( npix*npix );
    // type 0
    int ntotxy0;
    ntotxy0 = 0;
    std::valarray<float> mapxytot0( npix*npix );
    // type 1
    int ntotxy1;
    ntotxy1 = 0;
    std::valarray<float> mapxytot1( npix*npix );
    // type 2
    int ntotxy2;
    ntotxy2 = 0;
    std::valarray<float> mapxytot2( npix*npix );
    // type 3
    int ntotxy3;
    ntotxy3 = 0;
    std::valarray<float> mapxytot3( npix*npix );
    // type 4
    int ntotxy4;
    ntotxy4 = 0;
    std::valarray<float> mapxytot4( npix*npix );
    // type 5
    int ntotxy5;
    ntotxy5 = 0;
    std::valarray<float> mapxytot5( npix*npix );
  
    string snapnum;
    if( blsnap[nsnap]<10) snapnum = "00"+conv(blsnap[nsnap],fINT);
    else if(blsnap[nsnap]>=10 && blsnap[nsnap]<100 ) snapnum = "0"+conv(blsnap[nsnap],fINT);
    else snapnum = conv(blsnap[nsnap],fINT); 
    string File = pathsnap+"/snapdir_"+snapnum+"/"
      +simulation+"_snap_"+snapnum;
    string snappl;
    if( pll[nsnap]<10) snappl = "00"+conv(pll[nsnap],fINT);
    else if(pll[nsnap]>=10 && pll[nsnap]<100 ) snappl = "0"+conv(pll[nsnap],fINT);
    else snappl = conv(pll[nsnap],fINT);     
    string Filesub;

    if(subfiles!="NO"){
      Filesub = subfiles +"/groups_"+snapnum+"/subhalo_tab_"+snapnum;
    }else{
      cout << "         " << endl;
      cout << " I won't read SUBFIND catalogues " << endl;
      cout << "         " << endl;
    }
   
    float num_float1, num_float2, num_float3; 
        
    vector<double> mfof,xfof,yfof,zfof,m200,r200;
    vector<double> xfofc,yfofc,zfofc;
    vector<long> firstsubinfof,nsubinfof;
    vector<double> msub,xsub,ysub,zsub;
    vector<double> xsubc,ysubc,zsubc;
    vector<long> idPH,groupN;
    vector<double> velDisp,vmax,halfmassradius,rmax;
    if(subfiles!="NO"){
      read_SUBFIND(Filesub,mfof,xfofc,yfofc,zfofc,
		   firstsubinfof,nsubinfof,
		   m200,r200,
		   msub,xsubc,ysubc,zsubc,
		   idPH,velDisp,vmax,
		   halfmassradius,rmax,groupN);
    }
    int nfof = mfof.size();
    int nsub = msub.size();
    xfof.resize(nfof);
    yfof.resize(nfof);
    zfof.resize(nfof);
    xsub.resize(nsub);
    ysub.resize(nsub);
    zsub.resize(nsub);
    // to check that the SUBFIND files have been read correctly
    /*
      cout << mfof[0] << "  " << m200[0] << "  " << msub[0] << "  " << halfmassradius[0] << endl;
      cout << idPH[3255] << "  " << firstsubinfof[3255] << "  " << nsubinfof[3255] << endl;
      cout << xfof[0] << "  " << yfof[0] << "  " << zfof[0] << endl;
      cout << xsub[0] << "  " << ysub[0] << "  " << zsub[0] << endl;
      
      cout << mfof[1] << "  " << m200[1] << "  " << msub[1] << "  " << halfmassradius[1] << endl;
      cout << idPH[10255] << "  " << firstsubinfof[10255] << "  " << nsubinfof[10255] << endl;    
      cout << xfof[1] << "  " << yfof[1] << "  " << zfof[1] << endl;
      cout << xsub[1] << "  " << ysub[1] << "  " << zsub[1] << endl;
      
      cout << mfof[2] << "  " << m200[2] << "  " << msub[2] << "  " << halfmassradius[2] << endl;
      cout << idPH[nfof-1] << "  " << firstsubinfof[nfof-1] << "  " << nsubinfof[nfof-1] << endl;
      cout << xfof[2] << "  " << yfof[2] << "  " << zfof[2] << endl;
      cout << xsub[2] << "  " << ysub[2] << "  " << zsub[2] << endl;
      
      cout << "                  " << endl;
      cout << " number of fof  = " << nfof << endl;
      cout << " number of sub  = " << nsub << endl;
      cout << "                  " << endl;
      exit(1);
    */
    // fof in the three projections
    vector<int> idfof;
    vector<double> fofx, fofy, fofz;
    vector<double> zifof;
    vector<double> mfofsel;
    vector<double> m200sel,r200sel;
    // sub in the three projections
    vector<int> idsub;
    vector<int> idfirstsubsel;
    vector<int> idphfof;
    vector<double> subx, suby, subz;
    vector<double> zisub;
    vector<double> msubsel;
    vector<double> velDispsub;
    vector<double> vmaxsub,rmaxsub;
    vector<double> halfmassradiussub;

    // redshift and dl of the simulation
    double zsim, dlsim;
    if(noSNAP==0){
    }else{
      nfiles=1;
    }
    for (unsigned int ff=0; ff<nfiles; ff++){ 
      // map for each mass type
      std::valarray<float> mapxy0( npix*npix );
      std::valarray<float> mapxy1( npix*npix );
      std::valarray<float> mapxy2( npix*npix );
      std::valarray<float> mapxy3( npix*npix );
      std::valarray<float> mapxy4( npix*npix );
      std::valarray<float> mapxy5( npix*npix );

      // GADGET has 6 different particle type
      vector<float> xx0(0), yy0(0), zz0(0);
      vector<float> xx1(0), yy1(0), zz1(0);
      vector<float> xx2(0), yy2(0), zz2(0);
      vector<float> xx3(0), yy3(0), zz3(0);
      vector<float> xx4(0), yy4(0), zz4(0);
      vector<float> xx5(0), yy5(0), zz5(0);
      
      string file_in = File+"."+conv(ff,fINT);
      ifstream fin(file_in.c_str());
      if (!fin) {cerr <<"Error in opening the file: "<<file_in<<"!\n\a"; exit(1);}
      	
      cout <<"    reading the input file: "<<file_in<<endl;
      
      int32_t blockheader;
      
      fin.read((char *)&blockheader, sizeof(blockheader));
      
      cout << sizeof(blockheader) << endl;
      
      DATA data; fin >>data; 
      cout << sizeof(data) << endl;
      
      fin.read((char *)&blockheader, sizeof(blockheader));    
      cout << sizeof(blockheader) << endl;
      
      fin.read((char *)&blockheader, sizeof(blockheader));
      cout << sizeof(blockheader) << endl;
      
      // total number of particles as the sum of all of them
      int dim = data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5];
      cout << " .......................................................... " << endl;
      cout << "   number of particles in this snapshot: " << endl;
      cout << data.npart[0] << " " << data.npart[1] << " " << data.npart[2] 
	   << " " << data.npart[3] << " " << data.npart[4] << " " << data.npart[5] << endl;
      
      if(ff==0){
	// compute the comoving angular diamater distance at simulation redshift
	zsim = data.redshift;
	dlsim = getY(zl,dl,zsim);
	h0 = data.h;
	cout << "  " << endl;
	cout << "      __________________ COSMOLOGY __________________  " << endl;
	cout << " " << endl;
	om0 = data.om0;
	omL0 = data.oml;
	cout << "      Omegam = " << data.om0 << " " << "Omegal = " << data.oml << endl;
	cout << "           h = " << data.h   << " " << "BoxSize = " << data.boxsize << endl;
	cout << "      redshift = " << zsim <<   " " << "Dl (comoving) = " << dlsim << endl;
	if(abs(boxl - data.boxsize/1.e+3)>1.e-2 ){
	  cout << " set boxl and data.size differ ... check it! " << std:: endl;
	  cout << "  boxl = " << boxl << "  " << " data.boxsize = " << data.boxsize/1.e+3 << endl;
	  exit(1);
	}
	
	cout << "      _______________________________________________  " << endl;
	cout << " " << endl;
	cout << "   total number of partiles in the simulation: " << endl; 
	cout << data.npartTotal[0] << " " << data.npartTotal[1] << " " << data.npartTotal[2] 
	     << " " << data.npartTotal[3] << " " << data.npartTotal[4] << " " << data.npartTotal[5] << endl;
	cout << " " << endl;
	cout << "   particle type mass array: " << endl; 
	cout << data.massarr[0] << " " << data.massarr[1] << " " << data.massarr[2] 
	     << " " << data.massarr[3] << " " <<  data.massarr[4] << " " <<  data.massarr[5] << endl; 
	m0 = data.massarr[0];
	m1 = data.massarr[1];
	m2 = data.massarr[2];
	m3 = data.massarr[3];
	m4 = data.massarr[4];
	m5 = data.massarr[5];
      
	
	if(subfiles!="NO"){	
	  // rescale fof and sub unit with respect to the new center and the boxsize
	  for(int l=0;l<nfof;l++){
	    // mirror with respect to the center
	    xfofc[l] = sgnX[rcase]*(xfofc[l]/data.boxsize);
	    yfofc[l] = sgnY[rcase]*(yfofc[l]/data.boxsize);
	    zfofc[l] = sgnZ[rcase]*(zfofc[l]/data.boxsize);
	    // wrapping periodic condition
	    if(xfofc[l]>1.) xfofc[l] = xfofc[l] - 1.;
	    if(yfofc[l]>1.) yfofc[l] = yfofc[l] - 1.;
	    if(zfofc[l]>1.) zfofc[l] = zfofc[l] - 1.;
	    if(xfofc[l]<0.) xfofc[l] = 1. + xfofc[l];
	    if(yfofc[l]<0.) yfofc[l] = 1. + yfofc[l];
	    if(zfofc[l]<0.) zfofc[l] = 1. + zfofc[l];
	    // face of the dice
	    switch (face[rcase]){
	    case(1):
	      xfof[l] = xfofc[l];
	      yfof[l] = yfofc[l];
	      zfof[l] = zfofc[l];
	      break;
	    case(2):
	      xfof[l] = xfofc[l];
	      yfof[l] = zfofc[l];
	      zfof[l] = yfofc[l];
	      break;
	    case(3):
	      xfof[l] = yfofc[l];
	      yfof[l] = zfofc[l];
	      zfof[l] = xfofc[l];
	      break;
	    case(4):
	      xfof[l] = yfofc[l];
	      yfof[l] = xfofc[l];
	      zfof[l] = zfofc[l];
	      break;
	    case(5):
	      xfof[l] = zfofc[l];
	      yfof[l] = xfofc[l];
	      zfof[l] = yfofc[l];
	      break;
	    case(6):
	      xfof[l] = zfofc[l];
	      yfof[l] = yfofc[l];
	      zfof[l] = xfofc[l];
	      break;
	    }
	    // recenter
	    xfof[l] = xfof[l] - x0[rcase];
	    yfof[l] = yfof[l] - y0[rcase];
	    zfof[l] = zfof[l] - z0[rcase];
	    /*
	      switch (face[rcase]){
	      case(1):
	      xfof[l] = sgnX[rcase]*(xfofc[l]/data.boxsize - x0[rcase]);
	      yfof[l] = sgnY[rcase]*(yfofc[l]/data.boxsize - y0[rcase]);
	      zfof[l] = sgnZ[rcase]*(zfofc[l]/data.boxsize - z0[rcase]);
	      break;
	      case(2):
	      xfof[l] = sgnX[rcase]*(xfofc[l]/data.boxsize - x0[rcase]);
	      yfof[l] = sgnZ[rcase]*(zfofc[l]/data.boxsize - z0[rcase]);
	      zfof[l] = sgnY[rcase]*(yfofc[l]/data.boxsize - y0[rcase]);	    
	      break;
	      case(3):
	      xfof[l] = sgnY[rcase]*(yfofc[l]/data.boxsize - y0[rcase]);
	      yfof[l] = sgnZ[rcase]*(zfofc[l]/data.boxsize - z0[rcase]);
	      zfof[l] = sgnX[rcase]*(xfofc[l]/data.boxsize - x0[rcase]);	    
	      break;
	      case(4):
	      xfof[l] = sgnY[rcase]*(yfofc[l]/data.boxsize - y0[rcase]);
	      yfof[l] = sgnX[rcase]*(xfofc[l]/data.boxsize - x0[rcase]);
	      zfof[l] = sgnZ[rcase]*(zfofc[l]/data.boxsize - z0[rcase]);	    
	      break;
	      case(5):
	      xfof[l] = sgnZ[rcase]*(zfofc[l]/data.boxsize - z0[rcase]);
	      yfof[l] = sgnX[rcase]*(xfofc[l]/data.boxsize - x0[rcase]);
	      zfof[l] = sgnY[rcase]*(yfofc[l]/data.boxsize - y0[rcase]);	    
	      break;
	      case(6):
	      xfof[l] = sgnZ[rcase]*(zfofc[l]/data.boxsize - z0[rcase]);
	      yfof[l] = sgnY[rcase]*(yfofc[l]/data.boxsize - y0[rcase]);
	      zfof[l] = sgnX[rcase]*(xfofc[l]/data.boxsize - x0[rcase]);	    
	      break;
	      }
	    */
	    // wrapping periodic condition again!
	    if(xfof[l]>1.) xfof[l] = xfof[l] - 1.;
	    if(yfof[l]>1.) yfof[l] = yfof[l] - 1.;
	    if(zfof[l]>1.) zfof[l] = zfof[l] - 1.;
	    if(xfof[l]<0.) xfof[l] = 1. + xfof[l];
	    if(yfof[l]<0.) yfof[l] = 1. + yfof[l];
	    if(zfof[l]<0.) zfof[l] = 1. + zfof[l];
	    zfof[l]+=double(rcase); // pile the cones
	  }
	  
	  for(int l=0;l<nsub;l++){
	    xsubc[l] = sgnX[rcase]*(xsubc[l]/data.boxsize);
	    ysubc[l] = sgnY[rcase]*(ysubc[l]/data.boxsize);
	    zsubc[l] = sgnZ[rcase]*(zsubc[l]/data.boxsize);
	    // wrapping periodic condition 
	    if(xsubc[l]>1.) xsubc[l] = xsubc[l] - 1.;
	    if(ysubc[l]>1.) ysubc[l] = ysubc[l] - 1.;
	    if(zsubc[l]>1.) zsubc[l] = zsubc[l] - 1.;
	    if(xsubc[l]<0.) xsubc[l] = 1. + xsubc[l];
	    if(ysubc[l]<0.) ysubc[l] = 1. + ysubc[l];
	    if(zsubc[l]<0.) zsubc[l] = 1. + zsubc[l];
	    
	    switch (face[rcase]){
	    case(1):
	      xsub[l] = xsubc[l];
	      ysub[l] = ysubc[l];
	      zsub[l] = zsubc[l];
	      break;
	    case(2):
	      xsub[l] = xsubc[l];
	      ysub[l] = zsubc[l];
	      zsub[l] = ysubc[l];
	      break;
	    case(3):
	      xsub[l] = ysubc[l];
	      ysub[l] = zsubc[l];
	      zsub[l] = xsubc[l];
	      break;
	    case(4):
	      xsub[l] = ysubc[l];
	      ysub[l] = xsubc[l];
	      zsub[l] = zsubc[l];
	      break;
	    case(5):
	      xsub[l] = zsubc[l];
	      ysub[l] = xsubc[l];
	      zsub[l] = ysubc[l];
	      break;
	    case(6):
	      xsub[l] = zsubc[l];
	      ysub[l] = ysubc[l];
	      zsub[l] = xsubc[l];
	      break;	    
	    }
	    // recenter
	    xsub[l] = xsub[l] - x0[rcase];
	    ysub[l] = ysub[l] - y0[rcase];
	    zsub[l] = zsub[l] - z0[rcase];
	    /*
	      switch (face[rcase]){
	      case(1):
	      xsub[l] = sgnX[rcase]*(xsubc[l]/data.boxsize - x0[rcase]);
	      ysub[l] = sgnY[rcase]*(ysubc[l]/data.boxsize - y0[rcase]);
	      zsub[l] = sgnZ[rcase]*(zsubc[l]/data.boxsize - z0[rcase]);
	      break;
	      case(2):
	      xsub[l] = sgnX[rcase]*(xsubc[l]/data.boxsize - x0[rcase]);
	      ysub[l] = sgnZ[rcase]*(zsubc[l]/data.boxsize - z0[rcase]);
	      zsub[l] = sgnY[rcase]*(ysubc[l]/data.boxsize - y0[rcase]);
	      break;
	      case(3):
	      xsub[l] = sgnY[rcase]*(ysubc[l]/data.boxsize - y0[rcase]);
	      ysub[l] = sgnZ[rcase]*(zsubc[l]/data.boxsize - z0[rcase]);
	      zsub[l] = sgnX[rcase]*(xsubc[l]/data.boxsize - x0[rcase]);
	      break;
	      case(4):
	      xsub[l] = sgnY[rcase]*(ysubc[l]/data.boxsize - y0[rcase]);
	      ysub[l] = sgnX[rcase]*(xsubc[l]/data.boxsize - x0[rcase]);
	      zsub[l] = sgnZ[rcase]*(zsubc[l]/data.boxsize - z0[rcase]);
	      break;
	      case(5):
	      xsub[l] = sgnZ[rcase]*(zsubc[l]/data.boxsize - z0[rcase]);
	      ysub[l] = sgnX[rcase]*(xsubc[l]/data.boxsize - x0[rcase]);
	      zsub[l] = sgnY[rcase]*(ysubc[l]/data.boxsize - y0[rcase]);
	      break;
	      case(6):
	      xsub[l] = sgnZ[rcase]*(zsubc[l]/data.boxsize - z0[rcase]);
	      ysub[l] = sgnY[rcase]*(ysubc[l]/data.boxsize - y0[rcase]);
	      zsub[l] = sgnX[rcase]*(xsubc[l]/data.boxsize - x0[rcase]);
	      break;	    
	      }
	    */
	    // wrapping periodic condition again!
	    if(xsub[l]>1.) xsub[l] = xsub[l] - 1.;
	    if(ysub[l]>1.) ysub[l] = ysub[l] - 1.;
	    if(zsub[l]>1.) zsub[l] = zsub[l] - 1.;
	    if(xsub[l]<0.) xsub[l] = 1. + xsub[l];
	    if(ysub[l]<0.) ysub[l] = 1. + ysub[l];
	    if(zsub[l]<0.) zsub[l] = 1. + zsub[l];
	    zsub[l]+=double(rcase); // pile the cones
	  }
	  
	  double xmin=double(*min_element(xfof.begin(), xfof.end()));
	  double xmax=double(*max_element(xfof.begin(), xfof.end()));  
	  double ymin=double(*min_element(yfof.begin(), yfof.end()));
	  double ymax=double(*max_element(yfof.begin(), yfof.end()));  
	  double zmin=double(*min_element(zfof.begin(), zfof.end()));
	  double zmax=double(*max_element(zfof.begin(), zfof.end()));  
	  
	  cout << " min max fof " << endl;
	  cout << xmin << "   " << xmax << endl;
	  cout << ymin << "   " << ymax << endl;
	  cout << zmin << "   " << zmax << endl;
	  
	  xmin=double(*min_element(xsub.begin(), xsub.end()));
	  xmax=double(*max_element(xsub.begin(), xsub.end()));  
	  ymin=double(*min_element(ysub.begin(), ysub.end()));
	  ymax=double(*max_element(ysub.begin(), ysub.end()));  
	  zmin=double(*min_element(zsub.begin(), zsub.end()));
	  zmax=double(*max_element(zsub.begin(), zsub.end()));  
	  
	  cout << " min max sub " << endl;
	  cout << xmin << "   " << xmax << endl;
	  cout << ymin << "   " << ymax << endl;
	  cout << zmin << "   " << zmax << endl;
	  cout << "  " << endl;
	  
	  // now the unit are in [0,1]
	  // select for each case only fof in the field of view
	  for(int l=0;l<nfof;l++){
	    // 0.5 is because I set the eye in the center of the box!
	    double di = sqrt(pow(xfof[l]-0.5,2)+pow(yfof[l]-0.5,2)+pow(zfof[l],2))*data.boxsize/1.e+3;
	    if(di>=blD[nsnap] && di<blD2[nsnap]){
	      double rai,deci,dd;
	      getPolar(xfof[l]-0.5,yfof[l]-0.5,zfof[l],&rai,&deci,&dd);
	      if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){
		idfof.push_back(l);
		mfofsel.push_back(mfof[l]);
		m200sel.push_back(m200[l]);
		r200sel.push_back(r200[l]);
		fofx.push_back(rai);  // need to be in radiants
		fofy.push_back(deci); // need to be in radiants
		fofz.push_back(di); 
		zifof.push_back(getY(dl,zl,di));
	      }
	    }
	  }
	  cout << " number of fof in the field " << fofz.size() << endl;
	  for(int l=0;l<nsub;l++){
	    double di = sqrt(pow(xsub[l]-0.5,2)+pow(ysub[l]-0.5,2)+pow(zsub[l],2))*data.boxsize/1.e+3;
	    if(di>=blD[nsnap] && di<blD2[nsnap]){
	      double rai,deci,dd;
	      getPolar(xsub[l]-0.5,ysub[l]-0.5,zsub[l],&rai,&deci,&dd);
	      if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
		idsub.push_back(l);
		idphfof.push_back(groupN[l]);
		// idfirstsubsel.push_back(groupN[l]);
		velDispsub.push_back(velDisp[l]);
		vmaxsub.push_back(vmax[l]);
		rmaxsub.push_back(rmax[l]);
		halfmassradiussub.push_back(halfmassradius[l]);
		msubsel.push_back(msub[l]);
		subx.push_back(rai);  // need to be in radiants
		suby.push_back(deci); // need to be in radiants
		subz.push_back(di);
		zisub.push_back(getY(dl,zl,di));
	      }
	    }	  
	  }
	  cout << " number of subs in the field " << subz.size() << endl;      
	}
      }
      if(noSNAP==0){
	// type0
	for (int pp=0; pp<data.npart[0]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;
#ifdef SWAP	  
	  xb = sgnX[rcase]*((FloatSwap(num_float1)/data.boxsize));
	  yb = sgnY[rcase]*((FloatSwap(num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*((FloatSwap(num_float3)/data.boxsize));
#else
	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));
#endif
	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx0.push_back(x);
	  yy0.push_back(y);
	  zz0.push_back(z);
	}
	
	// type1
	for (int pp=data.npart[0]; pp<data.npart[0]+data.npart[1]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;
#ifdef SWAP	  
	  xb = sgnX[rcase]*((FloatSwap(num_float1)/data.boxsize));
	  yb = sgnY[rcase]*((FloatSwap(num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*((FloatSwap(num_float3)/data.boxsize));
#else
	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));
#endif
	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx1.push_back(x);
	  yy1.push_back(y);
	  zz1.push_back(z);
	}
	
	// type2
	for (int pp=data.npart[0]+data.npart[1]; pp<data.npart[0]+data.npart[1]+data.npart[2]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;
#ifdef SWAP	  
	  xb = sgnX[rcase]*((FloatSwap(num_float1)/data.boxsize));
	  yb = sgnY[rcase]*((FloatSwap(num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*((FloatSwap(num_float3)/data.boxsize));
#else
	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));
#endif
	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx2.push_back(x);
	  yy2.push_back(y);
	  zz2.push_back(z);
	}
	
	// type3
	for (int pp=data.npart[0]+data.npart[1]+data.npart[2]; 
	     pp<data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;
#ifdef SWAP	  
	  xb = sgnX[rcase]*((FloatSwap(num_float1)/data.boxsize));
	  yb = sgnY[rcase]*((FloatSwap(num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*((FloatSwap(num_float3)/data.boxsize));
#else
	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));
#endif
	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx3.push_back(x);
	  yy3.push_back(y);
	  zz3.push_back(z);
	}
	
	// type4
	for (int pp=data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]; 
	     pp<data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;
#ifdef SWAP	  
	  xb = sgnX[rcase]*((FloatSwap(num_float1)/data.boxsize));
	  yb = sgnY[rcase]*((FloatSwap(num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*((FloatSwap(num_float3)/data.boxsize));
#else
	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));
#endif
	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx4.push_back(x);
	  yy4.push_back(y);
	  zz4.push_back(z);
	}
	
	// type5
	for (int pp=data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]; 
	     pp<data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;
#ifdef SWAP	  
	  xb = sgnX[rcase]*((FloatSwap(num_float1)/data.boxsize));
	  yb = sgnY[rcase]*((FloatSwap(num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*((FloatSwap(num_float3)/data.boxsize));
#else
	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));
#endif
	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx5.push_back(x);
	  yy5.push_back(y);
	  zz5.push_back(z);
	}
	fin.clear(); fin.close();     
      }
    
      int n0 = xx0.size();
      int n1 = xx1.size();
      int n2 = xx2.size();
      int n3 = xx3.size();
      int n4 = xx4.size();
      int n5 = xx5.size();
      
      cout << "  " << endl;
      cout << n0 <<"   type (0) particles selected until now"<<endl;
      cout << n1 <<"   type (1) particles selected until now"<<endl;
      cout << n2 <<"   type (2) particles selected until now"<<endl;
      cout << n3 <<"   type (3) particles selected until now"<<endl;
      cout << n4 <<"   type (4) particles selected until now"<<endl;
      cout << n5 <<"   type (5) particles selected until now"<<endl;
      cout << "  " << endl;
      
      int totPartxy0;
      int totPartxy1;
      int totPartxy2;
      int totPartxy3;
      int totPartxy4;
      int totPartxy5;
      
      if(n0>0){
	// quadrate box
	double xmin=double(*min_element(xx0.begin(), xx0.end()));
	double xmax=double(*max_element(xx0.begin(), xx0.end()));  
	double ymin=double(*min_element(yy0.begin(), yy0.end()));
	double ymax=double(*max_element(yy0.begin(), yy0.end()));  
	double zmin=double(*min_element(zz0.begin(), zz0.end()));
	double zmax=double(*max_element(zz0.begin(), zz0.end()));  
	cout << " " << endl;
	cout << " n0 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  0 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 0 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n0;l++){
	  double di = sqrt(pow(xx0[l]-0.5,2) + pow(yy0[l]-0.5,2) + pow(zz0[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx0[l]-0.5,yy0[l]-0.5,zz0[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx0[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy0[l]-0.5)/fovinunitbox+0.5);
	    }	  
	  }
	}
	totPartxy0=xs.size();
	// cout << " n0: totPartxy0 " << totPartxy0 << endl;
	ntotxy0+=totPartxy0;
	if(totPartxy0>0){
	  mapxy0 = gridist(xs,ys,npix);
	}
	// re-normalize to the total mass!
	double mtot0=0;
	if(totPartxy0>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot0 += mapxy0[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy0[l]=mapxy0[l]/mtot0*totPartxy0*m0;
	  }
	  //std:: cout << " total mass in the map " << mtot0*m0 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy0*m0 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n1>0){
	// quadrate box
	double xmin=double(*min_element(xx1.begin(), xx1.end()));
	double xmax=double(*max_element(xx1.begin(), xx1.end()));  
	double ymin=double(*min_element(yy1.begin(), yy1.end()));
	double ymax=double(*max_element(yy1.begin(), yy1.end()));  
	double zmin=double(*min_element(zz1.begin(), zz1.end()));
	double zmax=double(*max_element(zz1.begin(), zz1.end()));  
	cout << " " << endl;
	cout << " n1 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  1 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 1 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n1;l++){
	  double di = sqrt(pow(xx1[l]-0.5,2) + pow(yy1[l]-0.5,2) + pow(zz1[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx1[l]-0.5,yy1[l]-0.5,zz1[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx1[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy1[l]-0.5)/fovinunitbox+0.5);
	    }	  
	  }
	}
	totPartxy1=xs.size();
	// cout << " n1: totPartxy1 " << totPartxy1 << endl;
	ntotxy1+=totPartxy1;
	if(totPartxy1>0){
	  mapxy1 = gridist(xs,ys,npix);
	}
	// re-normalize to the total mass!
	double mtot1=0;
	if(totPartxy1>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot1 += mapxy1[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy1[l]=mapxy1[l]/mtot1*totPartxy1*m1;
	  }
	  //std:: cout << " total mass in the map " << mtot1*m1 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy1*m1 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n2>0){
	// quadrate box
	double xmin=double(*min_element(xx2.begin(), xx2.end()));
	double xmax=double(*max_element(xx2.begin(), xx2.end()));  
	double ymin=double(*min_element(yy2.begin(), yy2.end()));
	double ymax=double(*max_element(yy2.begin(), yy2.end()));  
	double zmin=double(*min_element(zz2.begin(), zz2.end()));
	double zmax=double(*max_element(zz2.begin(), zz2.end()));  
	cout << " " << endl;
	cout << " n2 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  2 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 2 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n2;l++){
	  double di = sqrt(pow(xx2[l]-0.5,2) + pow(yy2[l]-0.5,2) + pow(zz2[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx2[l]-0.5,yy2[l]-0.5,zz2[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx2[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy2[l]-0.5)/fovinunitbox+0.5);
	    }	  
	  }
	}
	totPartxy2=xs.size();
	// cout << " n2: totPartxy2 " << totPartxy2 << endl;
	ntotxy2+=totPartxy2;
	if(totPartxy2>0){
	  mapxy2 = gridist(xs,ys,npix);
	}
	// re-normalize to the total mass!
	double mtot2=0;
	if(totPartxy2>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot2 += mapxy2[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy2[l]=mapxy2[l]/mtot2*totPartxy2*m2;
	  }
	  //std:: cout << " total mass in the map " << mtot2*m2 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy2*m2 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n3>0){
	// quadrate box
	double xmin=double(*min_element(xx3.begin(), xx3.end()));
	double xmax=double(*max_element(xx3.begin(), xx3.end()));  
	double ymin=double(*min_element(yy3.begin(), yy3.end()));
	double ymax=double(*max_element(yy3.begin(), yy3.end()));  
	double zmin=double(*min_element(zz3.begin(), zz3.end()));
	double zmax=double(*max_element(zz3.begin(), zz3.end()));  
	cout << " " << endl;
	cout << " n3 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  3 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 3 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n3;l++){
	  double di = sqrt(pow(xx3[l]-0.5,2) + pow(yy3[l]-0.5,2) + pow(zz3[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx3[l]-0.5,yy3[l]-0.5,zz3[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx3[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy3[l]-0.5)/fovinunitbox+0.5);
	    }	  
	  }
	}
	totPartxy3=xs.size();
	// cout << " n3: totPartxy3 " << totPartxy3 << endl;
	ntotxy3+=totPartxy3;
	if(totPartxy3>0){
	  mapxy3 = gridist(xs,ys,npix);
	}
	// re-normalize to the total mass!
	double mtot3=0;
	if(totPartxy3>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot3 += mapxy3[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy3[l]=mapxy3[l]/mtot3*totPartxy3*m3;
	  }
	  //std:: cout << " total mass in the map " << mtot3*m3 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy3*m3 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n4>0){
	// quadrate box
	double xmin=double(*min_element(xx4.begin(), xx4.end()));
	double xmax=double(*max_element(xx4.begin(), xx4.end()));  
	double ymin=double(*min_element(yy4.begin(), yy4.end()));
	double ymax=double(*max_element(yy4.begin(), yy4.end()));  
	double zmin=double(*min_element(zz4.begin(), zz4.end()));
	double zmax=double(*max_element(zz4.begin(), zz4.end()));  
	cout << " " << endl;
	cout << " n4 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  4 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 4 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n4;l++){
	  double di = sqrt(pow(xx4[l]-0.5,2) + pow(yy4[l]-0.5,2) + pow(zz4[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx4[l]-0.5,yy4[l]-0.5,zz4[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx4[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy4[l]-0.5)/fovinunitbox+0.5);
	    }	  
	  }
	}
	totPartxy4=xs.size();
	// cout << " n4: totPartxy4 " << totPartxy4 << endl;
	ntotxy4+=totPartxy4;
	if(totPartxy4>0){
	  mapxy4 = gridist(xs,ys,npix);
	}
	// re-normalize to the total mass!
	double mtot4=0;
	if(totPartxy4>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot4 += mapxy4[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy4[l]=mapxy4[l]/mtot4*totPartxy4*m4;
	  }
	  //std:: cout << " total mass in the map " << mtot4*m4 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy4*m4 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n5>0){
	// quadrate box
	double xmin=double(*min_element(xx5.begin(), xx5.end()));
	double xmax=double(*max_element(xx5.begin(), xx5.end()));  
	double ymin=double(*min_element(yy5.begin(), yy5.end()));
	double ymax=double(*max_element(yy5.begin(), yy5.end()));  
	double zmin=double(*min_element(zz5.begin(), zz5.end()));
	double zmax=double(*max_element(zz5.begin(), zz5.end()));  
	cout << " " << endl;
	cout << " n5 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  5 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 5 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n5;l++){
	  double di = sqrt(pow(xx5[l]-0.5,2) + pow(yy5[l]-0.5,2) + pow(zz5[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx5[l]-0.5,yy5[l]-0.5,zz5[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx5[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy5[l]-0.5)/fovinunitbox+0.5);
	    }	  
	  }
	}
	totPartxy5=xs.size();
	// cout << " n5: totPartxy5 " << totPartxy5 << endl;
	ntotxy5+=totPartxy5;
	if(totPartxy5>0){
	  mapxy5 = gridist(xs,ys,npix);
	}
	// re-normalize to the total mass!
	double mtot5=0;
	if(totPartxy5>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot5 += mapxy5[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy5[l]=mapxy5[l]/mtot5*totPartxy5*m5;
	  }
	  //std:: cout << " total mass in the map " << mtot5*m5 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy5*m5 << std:: endl; 
	  //exit(1);
	}
      }
      std:: cout << " maps done! " << std:: endl;
      // sum up all the maps
      for (int i=0; i<npix;i++)for(int j=0;j<npix;j++){
	  mapxytot[i+npix*j] +=( 
				mapxy0[i+npix*j]+
				mapxy1[i+npix*j]+
				mapxy2[i+npix*j]+
				mapxy3[i+npix*j]+
				mapxy4[i+npix*j]+
				mapxy5[i+npix*j]);
	  
	  mapxytot0[i+npix*j] += mapxy0[i+npix*j];
	  mapxytot1[i+npix*j] += mapxy1[i+npix*j];
	  mapxytot2[i+npix*j] += mapxy2[i+npix*j];
	  mapxytot3[i+npix*j] += mapxy3[i+npix*j];
	  mapxytot4[i+npix*j] += mapxy4[i+npix*j];
	  mapxytot5[i+npix*j] += mapxy5[i+npix*j];
	}
      std:: cout << " done map*tot " << std:: endl;
    }
    if(noSNAP==0){
      if(partinplanes=="ALL"){
	/**
	 * write image array(s) to FITS files all particles in a FITS file!
	 */
	long naxis = 2;
	long naxes[2]={ npix,npix };
	string fileoutput;
	// fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	fileoutput = simulation+"."+snappl+".plane_"+snpix+".fits";
	std::auto_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	std::vector<long> naxex( 2 );
	naxex[0]=npix;
	naxex[1]=npix;
	PHDU *phxy=&ffxy->pHDU();
	phxy->write( 1, npix*npix, mapxytot );
	// phxy->addKey ("x0",x0," unit of the boxsize");
	// phxy->addKey ("y0",y0," unit of the boxsize");  
	phxy->addKey ("REDSHIFT",zsim," "); 
	phxy->addKey ("PHYSICALSIZE",fov," "); 
	phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	phxy->addKey ("nparttype0",ntotxy0," "); 
	phxy->addKey ("nparttype1",ntotxy1," "); 
	phxy->addKey ("nparttype2",ntotxy2," "); 
	phxy->addKey ("nparttype3",ntotxy3," "); 
	phxy->addKey ("nparttype4",ntotxy4," "); 
	phxy->addKey ("nparttype5",ntotxy5," "); 
	phxy->addKey ("HUBBLE",h0," "); 
	phxy->addKey ("OMEGAMATTER",om0," "); 
	phxy->addKey ("OMEGALAMBDA",omL0," "); 
	phxy->addKey ("m0",m0," "); 
	phxy->addKey ("m1",m1," "); 
	phxy->addKey ("m2",m2," "); 
	phxy->addKey ("m3",m3," "); 
	phxy->addKey ("m4",m4," "); 
	phxy->addKey ("m5",m5," "); 
      }else{
	/**
	 * write image array(s) to FITS files each particle type in different planes
	 */
	// type0
	if(ntotxy0>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	  fileoutput = simulation+"."+snappl+".ptype0_plane_"+snpix+".fits";
	  std::auto_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot0 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype0",ntotxy0," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m0",m0," "); 
	}
	// type1
	if(ntotxy1>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	  fileoutput = simulation+"."+snappl+".ptype1_plane_"+snpix+".fits";
	  std::auto_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot1 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype1",ntotxy1," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m1",m1," "); 
	}
	// type2
	if(ntotxy2>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	  fileoutput = simulation+"."+snappl+".ptype2_plane_"+snpix+".fits";
	  std::auto_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot2 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype2",ntotxy2," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m2",m2," "); 
	}
	// type3
	if(ntotxy3>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	  fileoutput = simulation+"."+snappl+".ptype3_plane_"+snpix+".fits";
	  std::auto_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot3 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype3",ntotxy3," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m3",m3," "); 
	}
	// type4
	if(ntotxy4>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	  fileoutput = simulation+"."+snappl+".ptype4_plane_"+snpix+".fits";
	  std::auto_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot4 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype4",ntotxy4," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m4",m4," "); 
	}
	// type5
	if(ntotxy5>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	  fileoutput = simulation+"."+snappl+".ptype5_plane_"+snpix+".fits";
	  std::auto_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot5 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype5",ntotxy5," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m5",m5," "); 
	}     
      }
    }
    if(subfiles!="NO"){
      // write a file with the fof and sub in the filed of view
      string filefof,filesub;
      ofstream fofinfield,subinfield;
      
      // filefof = "fofinfield_"+simulation+"."+snapnum+".dat";
      filefof = "fofinfield_"+simulation+"."+snappl+".dat";
      cout << "  " << endl;
      cout << "writing " << mfofsel.size() << " fof in " << filefof << endl;      
      fofinfield.open(filefof.c_str());    
      for(int l=0;l<mfofsel.size();l++){
	fofinfield << idfof[l] << "  " 
		   << mfofsel[l] << "  "
		   << fofx[l] << "  " 
		   << fofy[l] << "  " 
		   << zifof[l] << "  " 
		   << fofz[l] << "  " 
		   << m200sel[l] << "  " << r200sel[l] << endl;
      }
      fofinfield.close();
      
      //filesub = "subinfield_"+simulation+"."+snapnum+".dat";
      filesub = "subinfield_"+simulation+"."+snappl+".dat";
      cout << "  " << endl;
      cout << "writing " << msubsel.size() << " subhaloes in " << filesub << endl;      
      subinfield.open(filesub.c_str());
      for(int l=0;l<msubsel.size();l++){
	subinfield << idsub[l] << "  "
		   << idphfof[l] << "  "
		   << subx[l] << "  " 
		   << suby[l] << "  " 
		   << subz[l] << "  " 
		   << zisub[l] << "  " 
		   << msubsel[l] << "  " 
		   << velDispsub[l] << "  "
		   << vmaxsub[l] << "  "  
		   << halfmassradiussub[l] << "  " 
		   << rmax[l] << endl;
      }
      subinfield.close();
    }
  }
  cout << " end of work ... ;-)  " << endl;
}
