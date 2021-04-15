#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<vector>
#include<iterator>
#include<algorithm>
#include<fstream>
#include<TH1F.h>
#include<TH2.h>
#include<TF1.h>
#include<TH1D.h>
#include<TH2D.h>
#include"TGraphErrors.h"
#include"TGraph2D.h"
#include"TStyle.h"
#include"TTree.h"
#include"TVector3.h"
#include"TFile.h"
#include"TMath.h"
#include"TROOT.h"
#include"TSystem.h"
#include"TFitResult.h"
#include"TFitResultPtr.h"
#include"TRandom3.h"
#include"Math/Minimizer.h"
#include"Math/Factory.h"
#include"Math/Functor.h"
#include"libReadWCSim/libReadWCSimProjectHeaders.h"

//***********************Constant Values*******************
const TVector3 detextents(370.096, 0.0, 0.0);
const TVector3 photonstartpos(150.0, 0.0, 0.0); //cm
const TVector3 sourcepos(150.0, 0.0, 0.0); //cm//used for the calculation of chi-square
const double energy = 3.09999;//eV
const double energy1 = 3.02438;//3.17948;//3.02438;//eV
const double wavelength  = 400.07807;//nm
const double wavelength1  = 409.94912;//389.95116;//409.94912;//nm
const double pr_toffs = 0.0;//ns offset added to the propagation(to test).
const double digi_toffs = 950.0;//-11.39//950 default embedded in WCSim;
const double digitoffs = 950.0;//950.014;//ns//used for chisquare test
const double ref_index = 1.34419;
const double ref_index1 = 1.34331;//1.34515;//1.34331;//for calculation of group velocity
const double refIndex  = 1.34419; // for chisquare test
const double oneoverc = 0.03335646; // ns/cm
const double solc = 29.9792; //ns/cm //speed of light - c
const double groupveloc = 21.6273;//fitted one//21.70855; //cm/ns:average of one above and 
                                       //one below energy of actual energy  
                                       //from table of energies in source files in WCSim.

const int numofit = 2;//set the number of iterations you want to perform

//***********************Classes*********************************
/// a singleton class to hold the pmt time information needed
/// to do the fit for the timing calibration
class PMTTimeInfo {

  //protected:
  // fInstance=NULL; 

private:

  static PMTTimeInfo * fInstance;
  PMTTimeInfo();
  PMTTimeInfo(const PMTTimeInfo&);
  PMTTimeInfo& operator = (const PMTTimeInfo&);
  std::vector< int > pmts;// vector of pmt numbers (keys to the map) that are being used (Data Member)
  std::vector<int> ranpmts;//for randomly chosen PMTs
  std::vector< TVector3 > pmtspos; // vector of pmt positions (Data Member)
  std::map< int , std::vector< double > > fpmt_times; // map between pmt number and a vector of digitizer times (one per flash) (Data Member)
  std::vector< double >randToffs;

public:

  static PMTTimeInfo* GetInstance();

  //Some member functions for PMTs and corresponding positions & time
  void add_pmts( const std::vector<int> & addpmts ){ pmts = addpmts ; }
  void add_ranpmts( const std::vector<int> & addranpmts){ranpmts.clear(); ranpmts = addranpmts;}
  void add_pmts_pos( const std::vector< TVector3 > & addpmtpos ){ pmtspos = addpmtpos ; }
  void add_rantime(const std::vector<double> & addrantoffs){randToffs.clear(); randToffs = addrantoffs;} 
  void add_time( int ipmt, double atime ){
    // check if pmt is in map already:
    if ( fpmt_times.count( ipmt ) > 0 ){
      fpmt_times[ipmt].push_back( atime );
    } else {
      fpmt_times[ipmt] = std::vector<double>{atime};
    }
  }
  
  const std::vector<int>& ranpmts2() {return ranpmts; }
  const std::vector<int>& pmts2() { return pmts; }
  const std::vector<TVector3 >& pmtpos2() { return pmtspos; }
  const std::vector<double>& times( int ipmt ) { return fpmt_times[ipmt]; }
  const std::vector<double>& randToffs2() { return randToffs; }
  
};

PMTTimeInfo * PMTTimeInfo::fInstance = NULL;

PMTTimeInfo *PMTTimeInfo::GetInstance()
{
  if(fInstance == NULL)
    fInstance = new PMTTimeInfo();
  return fInstance;
}

PMTTimeInfo::PMTTimeInfo(){}

void printchi2(double sx1,
	       double sy1,
	       double sz1,
	       double groupvelocity1,
	       double prtoffs1,
	       double chi2tot1)
{
  std::cout<<"sx: "<<sx1
	   <<"\tsy: "<<sy1
	   <<"\tsz: "<<sz1
	   <<"\tgroupvelocity: "<<groupvelocity1
	   <<"\tprtoffs: "<<prtoffs1
	   <<"\tchi2tot: "<<chi2tot1
	   <<std::endl;
}

void print2chi2(double pmt1,
		double pmt2,
		double pmt3,
		double pmt4,
		double pmt5,
		double pmt6,
		double pmt7,
		double pmt8,
		double pmt9,
		double pmt10)
{
  std::cout<<"pmt1: "<<pmt1
	   <<" | pmt2: "<<pmt2
	   <<" | pmt3: "<<pmt3
	   <<" | pmt4: "<<pmt4
	   <<" | pmt5: "<<pmt5
	   <<" | pmt6: "<<pmt6
	   <<" | pmt7: "<<pmt7
	   <<" | pmt8: "<<pmt8
	   <<" | pmt9: "<<pmt9
	   <<" | pmt10: "<<pmt10
	   <<std::endl;
}

//********************Calculation for chi square********************
std::vector<int> chosenPMTs(std::vector<int> &allrandomPMTs){
  std::vector<int>fillinvector = allrandomPMTs;
  int fillinvectorsize=fillinvector.size();
  for(int i=0; i<=fillinvectorsize-1; ++i){
    // std::cout<<"-;@ "<<i<<":"<<fillinvector[i]<<std::endl;
  }
  return fillinvector;
}


double TCalibChi2(const double * pars1){
  
  PMTTimeInfo * fti = PMTTimeInfo::GetInstance();
  const std::vector<int> & pmts =  fti->pmts2();
  const std::vector<TVector3> & pmtspos = fti->pmtpos2();
  const std::vector<double>& randToffs = fti->randToffs2();

  // //std::cout<<"numofPMTs:"<<numofpmts<<std::endl;
  // for (int ipmt = 0 ; ipmt<3; ++ipmt )
  //   int pmt = pmts[ ipmt ];
  // const std::vector<double>& times = fti->times( pmt );
  // for(int i=0; i<5; ++i)
  //   std::cout<<i<<": "<<times[i]<<std::endl;  
  // 
  // 

  double sx = pars1[0];
  double sy = pars1[1];
  double sz = pars1[2];
  //double n  = pars[ 3 ];
  double groupVel = pars1[3];
  double prtoffs = pars1[4];

  /* std::cout<<"pars1[0]:"<<pars1[0]
	   <<" | pars1[1]:"<<pars1[1]
	   <<" | pars1[2]:"<<pars1[2]
	   <<" | pars1[3]:"<<pars1[3]
	   <<" | pars1[4]:"<<pars1[4]
	   <<std::endl;
  */

  //double chi2val = 0.0;
  long double chi2tot = 0.0;
  int numofpmts = pmts.size();
  //std::cout<<"numofPMTs:"<<numofpmts<<std::endl;
  for (int ipmt = 0 ; ipmt<numofpmts; ++ipmt ){
    int pmt = pmts[ ipmt ];
    const TVector3 & rpmt = pmtspos[ipmt];
    double dist2 = 
      (rpmt.X() - sx)*(rpmt.X() - sx) +
      (rpmt.Y() - sy)*(rpmt.Y() - sy) +
      (rpmt.Z() - sz)*(rpmt.Z() - sz) ;
    double dist = std::sqrt( dist2 );
    double propagtime = dist/groupVel;
    //double propagtime = dist*n*oneoverc;
    //double tsum = 0;
    //double parsum = 0;
    //std::cout<<"propagtime: "<<propagtime<<std::endl;
    const std::vector<double>& times = fti->times( pmt );
    int numoftimes = times.size();
    for ( int itime = 0; itime < numoftimes; ++itime ){
      if (times[itime] < 0.01) continue;
      double dt = ( times[itime] - propagtime - prtoffs - randToffs[ipmt]);
      //double dt = ( times[itime] - propagtime );
      //tsum +=  times[itime]*times[itime];
      //parsum += pars[ipmt+11];
      chi2tot += dt * dt;
      //chi2val = dt*dt;	  
    }
    //std::cout << "parsum " << parsum << ", t0 " << t0 << ", n " << n << ", tsum " << tsum  << ", Dist " << dist << ", chi2tot " << chi2tot << std::endl;
  }
  //printchi2(sx,sy,sz,groupVel,prtoffs,chi2tot );	
  // may need to divide this by sigma^2 and to get reduced chi2
  // would need to find degrees of freedom as total number of flashes * npmts
  return chi2tot;
}


double TCalibChi2toffs(const double * pars2){
  
  PMTTimeInfo * fti = PMTTimeInfo::GetInstance();
  const std::vector<int> & ranpmts = fti->ranpmts2();
  const std::vector<int> & pmts =  fti->pmts2();
  const std::vector<TVector3> & pmtspos = fti->pmtpos2();
  int numofranpmts = ranpmts.size();

  double randtoffs [numofranpmts];
  for (int i=0; i<numofranpmts; ++i){
    randtoffs[i]=pars2[i];
  }
  double sx2 = pars2[numofranpmts];//149.388;
  double sy2 = pars2[numofranpmts+1];//0.00462268;
  double sz2 = pars2[numofranpmts+2];//-0.0023173;
  double groupVel2 = pars2[numofranpmts+3];//21.6273; 
  double prtoffs2 = pars2[numofranpmts+4];//-0.203427;

  long double chi2tot2 = 0.0;
  int pmtindex = 0;

  //loop over PMTs
  for (int ipmt = 0 ; ipmt<numofranpmts; ++ipmt ){
    pmtindex = ranpmts[ipmt];
    const TVector3 & rpmt2 = pmtspos[pmtindex];
    double dist22 = 
      (rpmt2.X() - sx2)*(rpmt2.X() - sx2) +
      (rpmt2.Y() - sy2)*(rpmt2.Y() - sy2) +
      (rpmt2.Z() - sz2)*(rpmt2.Z() - sz2) ;
    double dist2 = std::sqrt( dist22 );
    double propagtime2 = dist2/groupVel2;
    const std::vector<double>& times2 = fti->times( pmtindex );
    int numoftimes2 = times2.size();

    for ( int itime = 0; itime < numoftimes2; ++itime ){
      if (times2[itime] < 0.01) continue;
      double dt2 = ( times2[itime] - propagtime2 - prtoffs2 - randtoffs[ipmt]);
      chi2tot2 += dt2 * dt2;	  
    }
  }
  //std::cout<<"In minimization function"<<std::endl;
  //print2chi2(rantoffs[0], rantoffs[1],rantoffs[2],rantoffs[3],rantoffs[4],rantoffs[5],rantoffs[6],rantoffs[7],rantoffs[8],rantoffs[9]);	
  // may need to divide this by sigma^2 and to get reduced chi2
  // would need to find degrees of freedom as total number of flashes * npmts
  return chi2tot2;
}

//******************************************************************************
//-----------------------------------------------------------------------------*
//******************************************************************************



void mpmt_analysis_truetime(char * ininfile = "/home/jashan/software/triumfLab/runwcsim/infile_notoffs_noraysc.txt",
			    char * infile1 = "blanktoffs.txt",
			    bool verbose = false)
{
  //std::cout<<"check1"<<std::endl;

  std::ifstream fin;
  std::cout<<"Opening File "<<infile1<<std::endl;
  fin.open(infile1);
  if(!fin.is_open()){
    std::cerr<<"Failed to open file"<<infile1<<std::endl;
    return;
  }

  std::string junk;
  Long64_t Npts = 0;
  const int MAX = 20000;
  int PMT_toffs_tmp[MAX];
  double time_toffs_tmp[MAX];
  int PMT_toffs[15808];
  double time_toffs[15808];
  while(fin>>PMT_toffs_tmp[Npts]>>time_toffs_tmp[Npts]){
    Npts++;
  }
  std::cout<<"Read in "<<Npts<<" pmt time offsets"<<std::endl;
  fin.close();

  //***NOTE***
  //The index in WCSim start from 1
  //The index from the simulation output and the text file included starts from 0
  //for instance: tubeId starts indexing from 1 whereas we start looping over PMTs from 0 
  for (int i = 0; i<15808; i++){
    time_toffs[i]=time_toffs_tmp[i+1];
    PMT_toffs[i]=PMT_toffs_tmp[i+1];
  } //just to set up the index


  char * wcsimdirenv;
  wcsimdirenv  = getenv ("WCSIMDIR");
  if(wcsimdirenv!=NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
  }else{
    gSystem->Load("./libWCSimRoot.so");
  }


  TFile *file = 0;
  //Get the pointer to the tree from the file
  TTree * tree = NULL;
  std::ifstream in(ininfile);
  std::string infilename;
  while(in>>infilename){
    if(file ==0){
      file = new TFile(infilename.c_str(),"read");
      if(!file->IsOpen()){
	std::cout<<"Error, could not open file: "<<infilename.c_str()<<std::endl;
	return;
      }
      tree = (TTree*)file->Get("wcsimT");
      break;
    }
  }
  in.close();

  
  if ( tree == NULL ){
    std::cout<<"Could not find tree in file "<<infilename<<endl;
    return;
  } else {
    std::cout<<"Okay found tree in file "<<infilename<<endl;
  }


  //Get the number of events
  std::cout<<"About to get tree entries"<<std::endl;
  int nevent = tree->GetEntries();
  std::cout<<"Got tree entries."<<std::endl;
  printf("nevent %d\n",nevent);


  //create WCSimRootEvent to put stuff from the tree in
  WCSimRootEvent * wcsimrootsuperevent = new WCSimRootEvent();

  //Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  //Force deletion to prevent memory leak
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);


  //Geometry tree - only need 1 'event'
  TTree* geotree = (TTree*)file->Get("wcsimGeoT");
  WCSimRootGeom *geo = 0;
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose)std::cout<<"Geotree has"<<geotree->GetEntries()<<"entries"<<std::endl;
  if(geotree->GetEntries() == 0){
    exit(9);
  }
  geotree->GetEntry(0);
  double det_r = geo->GetWCCylRadius();
  double det_l = geo->GetWCCylLength();
  double det_x0 = geo->GetWCOffset(0);
  double det_y0 = geo->GetWCOffset(1);
  double det_z0 = geo->GetWCOffset(2);
  double det_xmin = 9e99;
  double det_xmax = -9e99;
  double det_ymin = 9e99;
  double det_ymax = -9e99;
  double det_zmin = 9e99;
  double det_zmax = -9e99;
  double PMTradius = geo->GetWCPMTRadius();
  int npmts = geo->GetWCNumPMT();
  for(int ipmt = 0; ipmt<npmts; ++ipmt){
    WCSimRootPMT pmt = geo->GetPMT(ipmt);
    if (pmt.GetPosition(0)>det_xmax) det_xmax = pmt.GetPosition(0);
    if (pmt.GetPosition(0)<det_xmin) det_xmin = pmt.GetPosition(0);
    if (pmt.GetPosition(1)>det_ymax) det_ymax = pmt.GetPosition(1);
    if (pmt.GetPosition(1)<det_ymin) det_ymin = pmt.GetPosition(1);
    if (pmt.GetPosition(2)>det_zmax) det_zmax = pmt.GetPosition(2);
    if (pmt.GetPosition(2)<det_zmin) det_zmin = pmt.GetPosition(2);
  }

 std::cout<<"PMT Radius = "<< PMTradius <<"cms"<<std::endl;
  std::cout<<"Detector center = ("<<det_x0<<","<<det_y0<<","<<det_z0<<")"<<std::endl;
  std::cout<<"Detector xextent = ("<<det_xmin<<","<<det_xmax<<")"
	   <<"yextent = ("<<det_ymin<<", "<<det_ymax<<")"
	   <<"zextent = ("<<det_zmin<<", "<<det_zmax<<")"<<std::endl;
  std::cout<<"Refractive Index is: "<<ref_index<<std::endl;

  //Options tree  - only need 1 event
  if(0){
    TTree*opttree = (TTree*)file->Get("wcsimRootOptionsT");
    WCSimRootOptions *opt = 0;
    opttree->SetBranchAddress("wcsimrootoptions", &opt);
    if(verbose)std::cout<<"Opttree has "<<opttree->GetEntries()<<"entries"<<std::endl;
    if(opttree->GetEntries()==0){
      exit(9);
    }
    opttree->GetEntry(0);
    opt->Print();
  }

  int perpmtplots = 0;
  int fittingplots = 0;
  std::cout<<"Enter 1 for plotting individual per PMT plots"<<std::endl;
  std::cin>>perpmtplots;
  if (perpmtplots == 1){
    std::cout<<"Enter 1 for fitting individual per PMT plots"<<std::endl;
    std::cin>>fittingplots;
  }
    
  //start with the main subevent, as it contains most of the info
  WCSimRootTrigger* wcsimrootevent;
  TFile* Plots = new TFile("raysc_100_notoffs_50evts.root","RECREATE");
  if (Plots->IsOpen())printf("File is opened successfully \n");

  TH1D *hits = new TH1D("Hits","PMT Hits; Hits; counts/bin", 250, 0., 10000.);
  TH1D *time = new TH1D("Time","PMT Times; Time; counts/bin", 200, 900., 1100.);
  TH1D * digioff_cor = new TH1D("digioff_cor", "Digital offset correction to time; T-digitoffs; counts/bin", 200, -50., 150.);
  TH1D *cortime1D = new TH1D("cortime1D", "Corrected time; T-digitoffs-propag_time; counts/bin", 300, -50., 150.);
  TH1D *cortime_offset = new TH1D("cortime_offset", "Corrected time using random time offset ; T-digitoffs-propag_time-time_toffs; counts/bin", 600, -50., 100.);
  TH1D *cortime_rtoffs = new TH1D("cortime_rtoffs", "Corrected time using random time offset ; T-digitoffs-time_toffs; counts/bin", 600, -50., 100.);
  TH1D *randtoffs = new TH1D("randtoffs", "random time offset of the PMTs; random timeoffset; counts/bin", 200.,-1.0, 1.0);
TH1D *minrandtoffs = new TH1D("minrandtoffs", "minimized random time offset of the PMTs; random timeoffset; counts/bin", 200.,-1.0, 1.0);
 // TH1D *minimizedXpos = new TH1D("minimizedXpos","minimized value of X-position over series of iterations; minimized X-pos; counts/bin", 200, 140.0, 160.0);
 // TH1D *minimizedYpos = new TH1D("minimizedYpos","minimized value of Y-position over series of iterations; minimized Y-pos; counts/bin", 200, -10.0, 10.0);
 // TH1D *minimizedZpos = new TH1D("minimizedZpos","minimized value of Y-position over series of iterations; minimized Y-pos; counts/bin", 200, -10.0, 10.0);
 // TH1D *minimizedGroupV = new TH1D("minimizedGroupV","minimized group velocity over series of iterations; minimized group velocity; counts/bin", 100., 15.0, 25.0); 
 // TH1D *minimizedproToffs = new TH1D("minimizedproToffs","minimized group velocity over series of iterations; minimized propagation time offset; counts/bin", 100., -10.0, 10.0); 
  TH1D *pro_time = new TH1D("PropagationTime","Propagation Time; Propagation time; counts/bin", 300, -50., 150.);
  TH1D *truetimes = new TH1D("truetime","True time when a photon hit and a photo-electron was produced; TrueTime; counts/bin",200, -50., 150.);
  TH1D *timedist = new TH1D("timedist", "Time distribution; True time - Propag_time; counts/bin", 1000, -1., 1.);
  TH1D *timedist1 = new TH1D("timedist1", "Time distribution using group velocity; True time - Propag_time; counts/bin", 1000, -1., 1.);
  TH1D *timedist_ene = new TH1D("timedist_ene", "Time distribution using group velocity using energy; True time - Propag_time_ene; counts/bin", 1000, -1., 1.);
TH1D *timedist_exacthitpos = new TH1D("timedist_exacthitpos", "Time distribution using exact hit position of the photon; True time - PhHitpos; counts/bin", 1000, -1., 1.);
  TH1D *htcalibfit = new TH1D("htcalibfit","Time Calibration peak using gaussian fit; fit peak for corrected time - random time-offset; counts/bin",600, -50., 100.);
  TH1D *charge = new TH1D("Charge","Charge; Charge; counts/bin",110, -10, 100);
  TH2D *TimeCharge = new TH2D("TimeCharge","Time-Charge; Time; Charge",160, 900., 1100., 70, -10., 60.);
  //Total charge for all the events
  TH1D *hqtot = new TH1D("hqtot", "Total Charge over 200 events; Total Charge(pC); counts/bin", 375, 18000., 31000.);
  TH1D *hpmtnum = new TH1D("hpmtnum", "Times PMT was hit; PMT number; count/PMT", npmts, -0.5, float(npmts)-0.5);
  TH1D *hnhitpmtperev = new TH1D("hnhitpmtperev","number of times each PMT is hit in each event; Nhits per PMT", 10, -5., 5.);
  TH2D * AvChargeVsProbability = new TH2D("AvChargeVsProbability", "Average Charge vs Event Probability; Average Charge; Event Probability", 240, 0., 12., 100, 0., 1.);
  TH1D *chisquarefit = new TH1D("chisquarefit","Chi2 likelihood using fitter; ((corDigiTime)^2); counts/bin", 600, -50., 100.);
  TH2D *chi2perpos = new TH2D("chi2perposition","Chi2 likelihood using fitter per position; (Position(cm); (chi2)", 74., -370., 370., 10000., 1.0, 1000.0 );
  TH2D *cortimeperpos = new TH2D("cortimeperposition","Corrected Time  per position; Position; (corTruetime(ns))", 74., -370., 370., 5500., -50., 500.  );
  TH2D *cortime2perpos = new TH2D("cortime2perposition","Corrected Time  per position; Position; (corTruetime)^2(ns)^2", 74., -370., 370., 15000., -50., 1000.  );
  TH1D *eventsperpmt_truecut = new TH1D("eventsperpmt_truecut", "Number of events over cut on true time; number of events; counts/bin", 310, -10., 300.);
  TH1D*chi2truetime = new TH1D("chi2truetime","Digitizer time in the chi2 calculation; digitizerTime(ns); counts/bin", 200., 900., 1100.);
  TH2D*cortrtimevsprtime = new TH2D("cortrtimevsprtime","Corrected True time vs propagation Time; Propagation Time(ns); Corrected True Time(ns)", 350, 0., 35., 100, 0., 1.);
  TH2D *RIperpmt = new TH2D ("RIperpmt","Distribution of Refractive Index over PMTs; PMT Number; RefractiveIndex", 15808., 0.0, 15808.0, 800., 1., 5.);
  TH1D*RIforidealcortr = new TH1D("RIforidealcortr", "Refractive Index Distribution for Ideal Corrected True Time(~0); Rrefractive Index; counts/bin", 3800., 1.2, 5.);
  

  //plots to see how PMTs are arranged -- one for side, top and bottom separately
  TH2D *pmtnum_side = new TH2D("pmtnum_side","Arrangement of PMTs on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *pmtnum_top = new TH2D("pmtnum_top","Arrangement of PMTs on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *pmtnum_bot = new TH2D("pmtnum_bot","Arrangements of PMTs on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);


  //plots to see eventsin PMTs -- one for side, top and bottom separately
  TH2D *ev_side = new TH2D("Events_side","Number of events observed in each PMT on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *ev_top = new TH2D("Events_top","Number of events observed in each PMT on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *ev_bot = new TH2D("Events_bot","Number of events observed in each PMT on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);


  //plots to see eventsin PMTs -- one for side, top and bottom separately
  TH2D *evsidecut_trueT = new TH2D("evsidecut_trueT","Number of events observed in each PMT on sideafter applying cut over true time; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *evtopcut_trueT = new TH2D("evtopcut_trueT","Number of events observed in each PMT on top after applying cut over true time; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *evbotcut_trueT = new TH2D("evbotcut_trueT","Number of events observed in each PMT on bottom after applying cut over true time; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  //plots of peak time for hits on each PMT -- one for side, top and bottom separately
  TH2D *htside_t = new TH2D("htside_t","Final time in each PMT as observed on digitizer on side; #phi(radians); y(cm)",1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *httop_t = new TH2D("httop_t","Final time in each PMT as observed on digitizer on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *htbot_t = new TH2D("htbot_t","Final time in each PMT as observed on digitizer on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  //plots of corrected time(corrected for digital time-offset, propagation time, random time-offset) for hits on each PMT -- one for side, top and bottom separately
  TH2D *htside_cor = new TH2D("htside_cor","Corrcted time(T-prop_time-Digi_Toffs-randomTime-offset) in each PMT on side; #phi(radians); y(cm)",1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *httop_cor = new TH2D("httop_cor","Corrected time (T-prop_time-Digi_Toffs-randomTime-offset) in each PMT on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *htbot_cor = new TH2D("htbot_cor","Corrected time (T-prop_time-Digi_Toffs-randomTime-offset) in each PMT on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  //plots of corrected time(corrected for digital time-offset, propagation,time) for hits on each PMT -- one for side, top and bottom separately
  TH2D *htside_cor1 = new TH2D("htside_cor1","Corrected time (T-prop_time-Digi_Toffs) in each PMT on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *httop_cor1 = new TH2D("httop_cor1","Corrected time (T-prop_time-Digi_Toffs) in each PMT on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *htbot_cor1 = new TH2D("htbot_cor1","Corrected time (T-prop_time-Digi_Toffs) in each PMT on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  //plots of corrected true time(corrected for propagation time) for hits on each PMT -- one for side, top and bottom separately
  TH2D *truetside = new TH2D("truetside","Corrected True time in each PMT on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *truettop = new TH2D("truettop","Corrected True time in each PMT on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *truetbot = new TH2D("truetbot","Corrected True time in each PMT on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

 //plots of corrected true time(corrected for propagation time of photon took to reach from source to exact hit position) for hits on each PMT -- one for side, top and bottom separately
  TH2D *truexactside = new TH2D("truexactside","Time a photon took to reach from source to exact hit position in each PMT on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *truexacttop = new TH2D("truexacttop","Time a photon took to reach from source to exact hit position in each PMT on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *truexactbot = new TH2D("truexactbot","Time a photon took to reach from source to exact hit position in each PMT on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  //plots for average charge in PMTs -- one for side, top and bottom separately
  TH2D *hqav_side = new TH2D("hqav_side","Average Charge in PMTs on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *hqav_top = new TH2D("hqav_top","Average Charge in PMTs on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *hqav_bot = new TH2D("hqav_bot","Average Charge in PMTs on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  TH2D *trueventside1 = new TH2D("trueventside1", "Total number of events experienced on side by each PMT in first peak; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax );
  TH2D *trueventtop1 = new TH2D("trueventtop1","Total number of events experienced on top by each PMT in first peak; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *trueventbot1 = new TH2D("trueventbot1","Total number of events experienced on bot by each PMT in first peak; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  TH2D *trueventside2 = new TH2D("trueventside2", "Total number of events experienced on side by each PMT in second peak; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax );
  TH2D *trueventtop2 = new TH2D("trueventtop2","Total number of events experienced on top by each PMT in second peak; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *trueventbot2 = new TH2D("trueventbot2","Total number of events experienced on bot by each PMT in second peak; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  TH2D *truetside_gv = new TH2D("truetside_gv","Corrected True time using group velocity in each PMT on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *truettop_gv = new TH2D("truettop_gv","Corrected True time using group velocity in each PMT on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *truetbot_gv = new TH2D("truetbot_gv","Corrected True time using group velocity in each PMT on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  TH2D *percentagehitsinpmt_side = new TH2D("percentagehitsinpmt_side","Percentage hits in each PMT on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  TH2D *percentagehitsinpmt_top = new TH2D("percentagehitsinpmt_top","Percentage hits in each PMT on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  TH2D *percentagehitsinpmt_bot = new TH2D("percentagehitsinpmt_bot","Percentage hits in each PMT on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

  /*//plots for average charge in PMTs -- one for side, top and bottom separately
    TH2D *meandtside = new TH2D("meandtside","True times in PMTs on side; #phi(radians); y(cm)", 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
    TH2D *meandttop = new TH2D("meandttop","True times in PMTs on top; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
    TH2D *meandtbot = new TH2D("meandtbot","True times in PMTs on bottom; x(cm); z(cm)", 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  */

    
  TH1D **timeperpmt;
  TH1D **timeperpmtcor_pr;
  TH1D **tperpmtcor;
  TH1D **correctimeperpmt;
  TH1D **truetimeperpmt;
  TDirectory *curdir = Plots->GetDirectory("/");
  TDirectory *timeplotdir = Plots->mkdir("TimePerPMT");
  timeplotdir->cd();
  timeperpmt = new TH1D*[npmts];
  timeperpmtcor_pr = new TH1D*[npmts];
  tperpmtcor = new TH1D*[npmts];
  correctimeperpmt = new TH1D*[npmts];
  truetimeperpmt = new TH1D*[npmts];
  char hname[100];
  char htitle[100];
  char pmtndir[10];
  for (int ipmt = 0; ipmt<npmts; ++ipmt){
    if(ipmt%200 ==0){
      sprintf(pmtndir, "%05d", ipmt);
      (timeplotdir->mkdir(pmtndir))->cd();
    }
    sprintf(hname,"DigiTime_pmt%05d",ipmt);
    sprintf(htitle,"Time for tubeID%05d; DigitizerTime(ns); counts/bin", ipmt);
    timeperpmt[ipmt] = new TH1D(hname, htitle, 400, 920., 1120.);

    sprintf(hname,"corTimePMT_pr%05d",ipmt);
    sprintf(htitle,"Corrected Time for tubeID(propagationCorrection)%05d; DigiTime-propagationTime(ns); counts/bin", ipmt);
    timeperpmtcor_pr[ipmt] = new TH1D(hname, htitle, 400, -10., 100.);

    sprintf(hname,"cor_prtime_pmt%05d",ipmt);
    sprintf(htitle,"Corrected Time for tubeID%05d; DigiTime-digiToffs-prop_time-randomToffs(ns); counts/bin", ipmt);
    tperpmtcor[ipmt] = new TH1D(hname, htitle, 400, -0., 100.); 

    sprintf(hname,"fullycor_time_pmt%05d",ipmt);
    sprintf(htitle,"Fully Corrected Time for tubeID%05d; DigiTime-relprop_time-digiToffs-randomToffs(ns); counts/bin", ipmt);
    correctimeperpmt[ipmt] = new TH1D(hname, htitle, 400, -60., 120.);
   
    sprintf(hname,"Truetime_pmt%05d",ipmt);
    sprintf(htitle,"True Time for tubeID%05d; True time - propagationTime(ns); counts/bin", ipmt);
    truetimeperpmt[ipmt] = new TH1D(hname, htitle, 110, -10., 100.);    
  }
  curdir->cd();
    
  TH1D **qperpmt;
  TDirectory *curdir1 = Plots->GetDirectory("/");
  TDirectory *qplotdir = Plots->mkdir("ChargePerPMT");
  qplotdir->cd();
  qperpmt = new TH1D*[npmts];
  for (int ipmt = 0; ipmt<npmts; ++ipmt){
    if(ipmt%200 ==0){
      sprintf(pmtndir, "%05d", ipmt);
      (qplotdir->mkdir(pmtndir))->cd();
    }
    sprintf(hname,"q_pmt%05d",ipmt);
    sprintf(htitle,"Charge for tubeID%05d; charge(adc); counts/bin", ipmt);
    qperpmt[ipmt] = new TH1D(hname, htitle, 100, -50., 50.);
  }
  curdir1->cd();

  TH1D **cortimeperevent;
  TDirectory *curdir2 = Plots->GetDirectory("/");
  TDirectory *tplotdir = Plots->mkdir("CorTimeperevent");
  tplotdir->cd();
  cortimeperevent = new TH1D*[nevent];
  char evdir[10];
  for (int ev = 0; ev<nevent; ++ev){
    if(ev%20 ==0){
      sprintf(evdir, "%05d", ev);
      (tplotdir->mkdir(evdir))->cd();
    }
    sprintf(hname,"corTime_event%05d", ev);
    sprintf(htitle,"CorTimeperEvent%05d; T-propag_time-time_toffs; counts/bin", ev);
    cortimeperevent[ev] = new TH1D(hname, htitle, 100, -50., 150.);
  }
  curdir2->cd();


  //int npos = (det_xmax - det_xmin);  
  TH2D **meandtside;
  TH2D **meandttop;
  TH2D **meandtbot;
  TDirectory *curdir3 = Plots->GetDirectory("/");
  TDirectory *dtplotdir = Plots->mkdir("meandt");
  dtplotdir->cd();
  meandtside = new TH2D*[30800];
  meandttop = new TH2D*[30800];
  meandtbot = new TH2D*[30800];
  char posdir[10];
  for(int pos = -400; pos <= 370; pos+=10){
    //std::cout<<"pos"<<pos<<std::endl;
    if(pos%100 == 0){
      //std::cout<<"pos1"<<pos<<std::endl;
      sprintf(posdir, "pos%05d", pos);
      (dtplotdir->mkdir(posdir))->cd();
    }
    sprintf(hname, "meandtside_pos%05d", pos);
    sprintf(htitle, "meandtside%05d; phi(radians); y(cm)", pos);
    meandtside[pos] = new TH2D(hname, htitle, 1800, -TMath::Pi(), TMath::Pi(), 500, det_ymin, det_ymax);
  
    sprintf(hname, "meandttop_pos%05d", pos);
    sprintf(htitle, "meandttop%05d; x(cm); z(cm)", pos);
    meandttop[pos] = new TH2D(hname, htitle, 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);

    sprintf(hname, "meantbot_pos%05d", pos);
    sprintf(htitle, "meandtbot%05d; x(cm); z(cm)", pos);
    meandtbot[pos] = new TH2D(hname, htitle, 500, det_xmin, det_xmax, 500, det_zmin, det_zmax);
  }
  curdir3->cd();


  TH1D **RIperPMT;
  TDirectory *curdir4 = Plots->GetDirectory("/");
  TDirectory *RIplotdir = Plots->mkdir("RIperPMT");
  RIplotdir->cd();
  RIperPMT = new TH1D*[npmts];
  char pmtdir[10];
  for (int ipmt = 0; ipmt<npmts; ++ipmt){
    if(ipmt%200==0){
      sprintf(pmtdir, "%05d", ipmt);
      (RIplotdir->mkdir(pmtdir))->cd();
    }
    sprintf(hname,"RI_pmt%05d", ipmt);
    sprintf(htitle,"Refractive Index distribution%05d; Refractive Index; counts/bin", ipmt);
    RIperPMT[ipmt] = new TH1D(hname, htitle, 200, 1., 3.0);
  }
  curdir4->cd();


  TH1D**minRanTOffs;
  TDirectory *curdir5 = Plots->GetDirectory("/");
  TDirectory *Rantoffsdir = Plots->mkdir("minRanTOffs");
  Rantoffsdir->cd();
  minRanTOffs = new TH1D*[numofit];
  for(int it = 0; it<numofit; ++it){
    sprintf(hname, "minRantoffs%05d", it);
    sprintf(htitle, "Minimized Random Time-offset of all the PMTs; minimized random timeoffset; counts/bin", it);
    minRanTOffs[it] = new TH1D(hname, htitle, 200., -1., 1.);
  }
  curdir5->cd();


  TH1D**corTruetime_it;
  TDirectory *curdir7 = Plots->GetDirectory("/");
  TDirectory *cortruetimedir = Plots->mkdir("corTruetime_it");
  cortruetimedir->cd();
  corTruetime_it = new TH1D*[numofit];
  for(int it = 0; it<numofit; ++it){
    sprintf(hname, "corTruetime%05d", it);
    sprintf(htitle, "Corrected True Time; True time - propagation time; counts/bin", it);
    corTruetime_it[it] = new TH1D(hname, htitle, 1000., -1., 1.);
  }
  curdir7->cd();

      
  int num_trig = 0;

  //Some pointer arrays to be filled and used at different points in algorithm
  int *counthitsperpmt = new int [npmts];
  double *q_av = new double[npmts];
  double *countevperpmt = new double[npmts];
  double *propag_time = new double[npmts];//propagation time calculated in digihits
  double *relpropag_time =new double[npmts];//propagation time relative to the propagation time of the nearest pmt to the source
  double *pr_time = new double[npmts];//propagation time for variable positions for chisq
  double *prop_time = new double[npmts];//propagation time calculated in truetimes.
  double *prop_time1 = new double[npmts];//propagation time for truetimes using group velocity.
  double *prop_time_ene = new double[npmts];
  double *PhHit_prop_time =new double[npmts];
  double *calcdist = new double [npmts];
  double *totqperpmt = new double[npmts];
  double *hitprob = new double[npmts];
  double *correc_time =  new double[npmts];
  double *cor_time = new double[npmts];
  double *totalq = new double[npmts];
  double *hitinPMT = new double[npmts];
  double *hitinPMT_Tprompt = new double[npmts];//for calculating hits in the prompt digitime 
  double *hitinPMT_corTprompt = new double[npmts];//for calculating hits in prompt cortime
  double *hitinPMT_T = new double[npmts];//for calculating hits in the reflected digitime
  double *hitinPMT_corT = new double[npmts];//for calculating hits in the reflected cortime
  double *hitinPMT_truecut = new double[npmts];//for making cut for the PMTs over true time
  double *tmaxbin = new double[npmts];
  double *cortmaxbin = new double[npmts];
  double *digitime = new double[npmts];
  double *meancor_time = new double[npmts];
  double *hitsperpmt1 = new double[npmts];
  double *cortruetime = new double[npmts];
  double *cortruetime1 = new double[npmts];
  double *truemaxbin = new double[npmts];
  double *truetime1 = new double[npmts];
  // double *cortimesigma;
  double *truetmedian = new double[npmts];
  double chisqtot = 0.;
  TF1 *fitcor;
  double *chisqneum = new double[npmts];
  double *chisqdenom =  new double[npmts];
  double *chisqperpmt = new double[npmts];
  double *chisqneum_truecut = new double[npmts];
  double *chisqdenom_truecut =  new double[npmts];
  double *chisqperpmt_truecut = new double[npmts];
  double *truelobound = new double[npmts];
  double *truehibound = new double[npmts];
  int *totalCount = new int[npmts];
  double *countTruehitsperpmtincut = new double[npmts];
  double *countTruehitsperpmt = new double[npmts];
  //int *randPMTs = new int[10];
  int *countsinpmt1 = new int[npmts];//for counting the number of hits in true prompt region first peak
  int *countsinpmt2 = new int[npmts];//for counting the number of hits in true prompt region second peak
  std::vector<int>includePMT; 
  std::vector<int>includeranPMT; 
  std::vector<TVector3> posPMT;
  double meanref_index = (ref_index + ref_index1)/2;//mean of ref_index for group velocity
  double energyderivative = (ref_index1-ref_index)/std::log(energy1/energy); 
  double groupvelene = solc/(meanref_index+energyderivative);
  std::cout<<"groupvelocity_ene: "<<groupvelene<<std::endl;
  double lambdaovern = wavelength/ref_index;
  double derivative = (ref_index1-ref_index)/(wavelength1-wavelength);
  double covern = solc/ref_index;
  double groupvel = (covern*(1.0+(lambdaovern*derivative)));
  std::cout<<"groupvelocity: "<<groupvel<<std::endl;
  double *percentagehitsinPMTs = new double[npmts];
  //int *PMTsnoticed = new int[arraylength]; 


  //Now loop over events
  in.open(ininfile);
  while(in>>infilename){
    file =new TFile(infilename.c_str(),"read");
    if(!file->IsOpen()){
      std::cout<<"Error, could not open input file:"<<infilename.c_str()<<std::endl;
      return;
    }
    tree = (TTree*)file->Get("wcsimT");
    //Set the branch address for reading from the tree
    TBranch *branch= tree->GetBranch("wcsimrootevent");
    branch->SetAddress(&wcsimrootsuperevent);
    //Force deletion to prevent memory leak
    tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
    int nevent = tree->GetEntries();
    printf("nevent in file %s is %d\n", infilename.c_str(), nevent);


    //*********EVENT LOOP begins here************
    //for(int ev = 0; ev<nevent; ev++)
   for(int ev = 0; ev<nevent; ev++)
      {
	for (int ipmt=0; ipmt<npmts; ++ipmt){
	  counthitsperpmt[ipmt]=0;
	}

	//std::cout<<"check5"<<std::endl;

	double totalcharge = 0.;


	printf("eventNumber %d\n",ev);
	//Read the event from the tree into the WCSinRootEvent instance
	int itree = tree->LoadTree(ev);
	//cout<<"LoadTree replies"<<itree<<endl;
	int nbytes = tree->GetEntry(ev,1);
	if(verbose)
	  std::cout<<"GetEntry replies"<<nbytes<<std::endl;
	wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
	if (verbose){
	  printf("****************************************************");
	  printf("Evt, date %f %f\n",wcsimrootevent->GetHeader()->GetEvtNum(),
		 wcsimrootevent->GetHeader()->GetDate());
	  printf("Mode %d\n", wcsimrootevent->GetMode());
	  printf("Number of subevents %d\n",wcsimrootsuperevent->GetNumberOfSubEvents());
	  printf("Vtxvol %f\n", wcsimrootevent->GetVtxvol());
	  printf("Vtx %d %d\n", wcsimrootevent->GetVtx(0), wcsimrootevent->GetVtx(1));
	}

	if (verbose){
	  printf("Jmu %d\n", wcsimrootevent->GetJmu());
	  printf("Npar %d\n", wcsimrootevent->GetNpar());
	  printf("Ntrack %d\n", wcsimrootevent->GetNtrack());
	}

	//Now read the tracks in the events
	//Get the number of tracks
	int ntrack = wcsimrootevent->GetNtrack();
	if (verbose) printf("ntracks=%d\n",ntrack);

	//Loop through elements in the TClonesArray of WCSimTracks
	for (int i=0; i<ntrack; i++)
	  {
	    TObject *element = (wcsimrootevent->GetTracks())->At(i);
	    WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
	    if(verbose){
	      printf("Track input: %d\n",wcsimroottrack->GetIpnu());
	      printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());

	      for(int j=0; j<3; j++)
		printf("Track dir: %d %f\n",j, wcsimroottrack->GetDir(j));
		  
	      printf("Track energy: %f\n", wcsimroottrack->GetE());
	      printf("Track momentum: %f\n",wcsimroottrack->GetP());
	      printf("Track mass: %f\n",wcsimroottrack->GetM());
	    }
	  }//end of loop over tracks

	//Now look at the cherenkov hits
	int ncherenkovhits = wcsimrootevent->GetNcherenkovhits();
	int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();

	hits->Fill(ncherenkovdigihits);
	if(verbose){
	  printf("node id: %i\n", ev);
	  printf("Ncherenkovhits %d\n", ncherenkovhits);
	  printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
	  std::cout<<"RAW HITs: "<<std::endl;
	}

	//std::cout<<"check5b"<<std::endl;

	hitinPMT[npmts] = 0;
	hitinPMT_Tprompt[npmts] = 0.;
	hitinPMT_corTprompt[npmts] = 0.;
	hitinPMT_truecut[npmts] = 0.;
	hitinPMT_T[npmts] = 0.;
	hitinPMT_corT[npmts] = 0;


	for (int i =0; i<ncherenkovdigihits; i++){
	  TObject *Digi = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
	  WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Digi);
	  int tubeId = wcsimrootcherenkovdigihit->GetTubeId();
	  ( counthitsperpmt[tubeId-1]) ++;
	  WCSimRootPMT pos = geo->GetPMT(tubeId-1);
	  double posX = pos.GetPosition(0);
	  double posY = pos.GetPosition(1);
	  double posZ = pos.GetPosition(2);
	  double q = wcsimrootcherenkovdigihit->GetQ();
	  double t = wcsimrootcherenkovdigihit->GetT();

	      
	  double ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
	  hitsperpmt1[tubeId-1]= ncherenkovdigihits;
	  hitinPMT[tubeId-1] += 1; 
	  charge->Fill(q);
	  time->Fill(t);
	  digitime[tubeId-1]=t;
	  TimeCharge->Fill(t,q);
	  totqperpmt[tubeId-1] = q;
	  totalq[tubeId-1] +=q;
	  hpmtnum->Fill(float(tubeId-1));
	  if(tubeId<=npmts && tubeId>0){

	    TVector3 pmtpos(pos.GetPosition(0), pos.GetPosition(1), pos.GetPosition(2));
	    TVector3 rprop = pmtpos - photonstartpos;//corresponds to real WCSim source position
	    double dist= rprop.Mag();
	    propag_time[tubeId-1] = (dist*ref_index)*oneoverc+pr_toffs;
	    TVector3 nearpro = detextents - photonstartpos;
	    double nearprodist = nearpro.Mag();
	    double nearprotime = nearprodist*ref_index*oneoverc;
	    relpropag_time[tubeId-1] = propag_time[tubeId -1] - nearprotime;
	    digioff_cor->Fill(t-digi_toffs);
	    
	    if (perpmtplots == 1){
	      qperpmt[tubeId-1]->Fill(q);
	      timeperpmt[tubeId-1]->Fill(t);
	      timeperpmtcor_pr[tubeId-1]->Fill(t-digi_toffs-propag_time[tubeId-1]);
	      tperpmtcor[tubeId-1]->Fill(t-propag_time[tubeId-1]-time_toffs[tubeId-1]);
	      TVector3 rpropag = pmtpos - sourcepos;
	      double dist1 = rpropag.Mag();
	      pr_time[tubeId-1] = ((dist1*refIndex)*oneoverc)+pr_toffs;
	      
	      tmaxbin[tubeId-1]=timeperpmt[tubeId-1]->GetXaxis()->GetBinCenter(timeperpmt[tubeId-1]->GetMaximumBin());
	      double lobound = tmaxbin[tubeId-1]-2.0;
	      double hibound = tmaxbin[tubeId-1]+2.0;
	      if (digitime[tubeId-1]>=lobound && digitime[tubeId-1]<=hibound){
		hitinPMT_Tprompt[tubeId-1]+=1;
	      }
	      else{
		hitinPMT_T[tubeId-1]+=1;
	      }
	    }//end of if(perpmtplots)
	    
	    correctimeperpmt[tubeId-1]->Fill(t-digi_toffs-relpropag_time[tubeId-1]-time_toffs[tubeId-1]);
	    cortime1D->Fill(t-digi_toffs-relpropag_time[tubeId-1]);
	    cortime_offset->Fill(t-digi_toffs-propag_time[tubeId-1]-time_toffs[tubeId-1]);
	    cortime_rtoffs->Fill(t-digi_toffs-time_toffs[tubeId-1]);
	    pro_time->Fill(propag_time[tubeId-1]);
	    correc_time[tubeId-1] = t-digi_toffs-relpropag_time[tubeId-1]-time_toffs[tubeId-1];
	   
	

	  
	    cortmaxbin[tubeId-1] = correctimeperpmt[tubeId-1]->GetXaxis()->GetBinCenter(correctimeperpmt[tubeId-1]->GetMaximumBin());
	    double corlobound = cortmaxbin[tubeId-1]-3.5;
	    double corhibound = cortmaxbin[tubeId-1]+3.5;
	    if (correc_time[tubeId-1]>=corlobound && correc_time[tubeId-1]<=corhibound){
	      hitinPMT_corTprompt[tubeId-1]+=1;
	    }
	    else{
	      hitinPMT_corT[tubeId-1]+=1;
	    }

	    cortimeperevent[ev]->Fill(correc_time[tubeId-1]);
	  }	
	}//end of loop over ncherenkovdigihits
       	  
	//Grab the big arrays of times and parent IDs
	TClonesArray * timeArray = wcsimrootevent->GetCherenkovHitTimes();
	int totalPe = 0;
	countsinpmt1[npmts] = 0;
	countsinpmt2[npmts] = 0;

	//Loop through elements in the TClonesArray of WCSimRootCherenkovHits
	for(int i=0; i<ncherenkovhits; i++ )
	  {
	    TObject* Hit = (wcsimrootevent->GetCherenkovHits())->At(i);
	    WCSimRootCherenkovHit *wcsimrootcherenkovhit = 
	      dynamic_cast<WCSimRootCherenkovHit*>(Hit);
	    int tubeNumber = wcsimrootcherenkovhit->GetTubeID();
	    int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	    int peForTube = wcsimrootcherenkovhit->GetTotalPe(1);
	    WCSimRootPMT pmt = geo->GetPMT(tubeNumber-1);
	    totalPe += peForTube;

	    double noverc = ref_index*oneoverc;//phase velocity used for calculation of propagation time
	    double covern = solc/ref_index;
	    double lambdaovern = wavelength/ref_index;
	    double meanref_index = (ref_index + ref_index1)/2;//mean of ref_index for group velocity
	    double energyderivative = (ref_index1-ref_index)/std::log(energy1/energy); 
	    double derivative = (ref_index1-ref_index)/(wavelength1-wavelength);
	    double groupvel = (covern*(1.0+(lambdaovern*derivative)));//group velocity for calculation of propagation time
	    double groupvelene = solc/(meanref_index+energyderivative);
	    groupvelene = 21.70855;
	    //std::cout<<"group velocity is:"<<groupvel<<std::endl;
	    //std::cout<<"phase velocity is:"<<covern<<std::endl;

	    WCSimRootPMT pos = geo->GetPMT(tubeNumber-1);
	    TVector3 pmtpos1(pos.GetPosition(0), pos.GetPosition(1), pos.GetPosition(2));
	    TVector3 pmtOrientation(pos.GetOrientation(0), pos.GetOrientation(1), pos.GetOrientation(2));
	    // pmt position given by WCSim is center of sphere. Shift that to center of cathode
	    TVector3 radius_unitvector = PMTradius*pmtOrientation;//  Rnhat radius of PMt times unit vector direction //Eqn2
	    pmtpos1 += radius_unitvector;
	    TVector3 distdiff = pmtpos1 - photonstartpos;//Eqn1  rN
	    TVector3 r_cvector = (distdiff - radius_unitvector);//length from source position to center of the PMT sphere.//Eqn1-Eqn2
	    double r_c = r_cvector.Mag();

	    TVector3 r_cPh = ((r_c - PMTradius)/r_c)*r_cvector;// rc distance between source and top the PMT where the actual photon is hit.
	    //r_cPh += radius_unitvector;
	    double phHitdist = r_cPh.Mag();
	    //std::cout<< "pmtradius:" << PMTradius << "|| pmtorientation_X:" << pmtOrientation.X()  << " || pmtorientation_Y:" << pmtOrientation.Y() << " || pmtorientation_Z:" << pmtOrientation.Z() << " || radius_unitvector_X:" << radius_unitvector.X() << " || radius_unitvector_Y:" << radius_unitvector.Y() << " || radius_unitvector_Z:" << radius_unitvector.Z() << " || distdiff_X:" << distdiff.X() << " || distdiff_Y:" << distdiff.Y() << " || distdiff_Z:" << distdiff.Z() << " || r_c: " << r_c << " || r_cvector_X: " << r_cvector.X() << " || r_cvector_Y: " << r_cvector.Y()  << " || r_cvector_Z: " << r_cvector.Z()<< " || phHitdist: " << phHitdist << std::endl;
	    double distMag = distdiff.Mag();
	    prop_time[tubeNumber-1] = (distMag*noverc)+pr_toffs;//for chisquare calculations
	    prop_time1[tubeNumber-1] = (distMag/groupvel);
	    prop_time_ene[tubeNumber-1] = (distMag/groupvelene);
	    PhHit_prop_time[tubeNumber-1] = (phHitdist/groupvel);
	    calcdist[tubeNumber-1]= solc/distMag;
	    //if(i<10) //Only printfirst 10 tubes
	    if(verbose) printf("Total Pe: %d times(",peForTube);
	    for(int j = timeArrayIndex; j < timeArrayIndex +peForTube; j++)
	      {
		WCSimRootCherenkovHitTime* HitTime = 
		  dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(j));
		if (verbose)printf("%6.2f",HitTime->GetTruetime());
		if(verbose)std::cout<<")"<<std::endl;
		truetimes->Fill(HitTime->GetTruetime());
		countTruehitsperpmt[tubeNumber-1]++;
		cortruetime[tubeNumber-1]=HitTime->GetTruetime()-prop_time[tubeNumber-1];
		cortruetime1[tubeNumber-1]=HitTime->GetTruetime()-prop_time1[tubeNumber-1]+10.;
		truetime1[tubeNumber-1] = HitTime->GetTruetime();
		RIperPMT[tubeNumber-1]->Fill(calcdist[tubeNumber-1]*truetime1[tubeNumber-1]);
		RIforidealcortr->Fill(calcdist[tubeNumber-1]*truetime1[tubeNumber-1]);
		truetimeperpmt[tubeNumber-1]->Fill(HitTime->GetTruetime());
		timedist->Fill(HitTime->GetTruetime()-prop_time[tubeNumber-1]);
		timedist1->Fill(HitTime->GetTruetime()-prop_time1[tubeNumber-1]);
		timedist_ene->Fill(HitTime->GetTruetime()-prop_time_ene[tubeNumber-1]);
		timedist_exacthitpos->Fill(HitTime->GetTruetime()-PhHit_prop_time[tubeNumber-1]);
		cortrtimevsprtime->Fill(prop_time[tubeNumber-1],cortruetime[tubeNumber-1]);
	      } 
	   
	  }//end of loop over cherenkov hits
	
	//Look at digitized hit information
	//Get the number of Digitized hits
	//Loop over sub events           
	if(verbose)std::cout<<"Digitized hits: "<<std::endl;
	for (int index = 0; index<wcsimrootsuperevent->GetNumberOfEvents(); index++)
	  {
	    if(index > 0) continue;
	    wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
	    if(verbose)std::cout<<"sub event number = "<<index<<"\n";
	    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
	    if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
	    if(ncherenkovdigihits>0) num_trig++;
	    for(int i = 0; i<ncherenkovdigihits; i++)
	      {
		//Loop through elelments in TCloneArray of WCSimRootCherenkovDigihits
		TObject *element  = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
		WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
		  dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
		if (verbose){
		  if(i<10)
		    printf("q, t,ttubeid: %f %f %d \n",wcsimrootcherenkovdigihit->GetQ(),
			   wcsimrootcherenkovdigihit->GetT(),wcsimrootcherenkovdigihit->GetTubeId());
		}
	      }//end of loop over cherenkov digihits
	  }//end of loop over trigger
      
	//reinitialize super event between loops
	wcsimrootsuperevent->ReInitialize();
	for(int ipmt=0; ipmt<npmts; ++ipmt){
	  WCSimRootPMT pos = geo->GetPMT(ipmt);
	  double posX = pos.GetPosition(0);
	  double posY = pos.GetPosition(1);
	  double posZ = pos.GetPosition(2);	      
	  TVector3 pmtpos(pos.GetPosition(0), pos.GetPosition(1), pos.GetPosition(2));
	  TVector3 rpropag = pmtpos - sourcepos;
	  double dist1 = rpropag.Mag();
	  //fill plot of number of times each PMT is hit in each event
	  hnhitpmtperev->Fill(float(counthitsperpmt[ipmt]));
	  if (perpmtplots == 1){
	    if (posY>=det_ymax-10.0){
	      ev_top->Fill(posX, posZ, counthitsperpmt[ipmt]);
	    }
	    else if (posY<=det_ymin+10.0){
	      ev_bot->Fill(posX, posZ, counthitsperpmt[ipmt]);
	    }
	    else{
	      double ang = atan2(posZ, posX);
	      ev_side->Fill(ang, posY, counthitsperpmt[ipmt]);
	    }
	    countevperpmt[ipmt] +=(double)counthitsperpmt[ipmt];
	  }//end of if(perpmtplots)
	  totalcharge += (double)totqperpmt[ipmt];
	      
	}//end of loop over PMTs
	hqtot->Fill(totalcharge);
      }//end of loop over events
  }//end of loop over files
  Plots->cd();


  for(int ev = 0; ev<nevent; ev++){
    //std::cout<<"event"<<ev<<std::endl;
    //printf("eventNumber %d\n",ev);
    int itree = tree->LoadTree(ev);
    int nbytes = tree->GetEntry(ev,1);
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    int ncherenkovhits = wcsimrootevent->GetNcherenkovhits();
    TClonesArray* timeArray = wcsimrootevent->GetCherenkovHitTimes();
   
    for (int i =0; i<ncherenkovhits; i++){
      TObject *Hit = (wcsimrootevent->GetCherenkovHits())->At(i);
      WCSimRootCherenkovHit *wcsimrootcherenkovhit = dynamic_cast<WCSimRootCherenkovHit*>(Hit);
      int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
      int peForTube = wcsimrootcherenkovhit->GetTotalPe(1);
      int tubeID = wcsimrootcherenkovhit->GetTubeID();
      int pmt = tubeID-1;
      for (int j = timeArrayIndex; j<timeArrayIndex+peForTube; j++)
	{
	  WCSimRootCherenkovHitTime *HitTime = dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(j));
   	  if(hitinPMT_corTprompt[pmt]>0){	
	    double t = HitTime->GetTruetime();
	    double timediff = t - prop_time[pmt];
	    double timeinpmt =0.;
	    if (timediff>=0. && timediff<=0.35)
	      {
		countsinpmt1[pmt] +=1; 
	      }
	    else if(timediff>=0.35 && timediff<= 0.75){
	      countsinpmt2[pmt] +=1;
	    }
	  } 
	}
    }
  }



  //fitting parameters to fit a gaussian peak and get the mean value
  double rawtubemean[MAX];
  double cortubemean[MAX];
  cortubemean[0] = 0.;
  double cortubeerr[MAX];
  double zeroes[MAX];
  double pmtnum[MAX];
  char fcn_name[40];
  char fcor_name[40];
  char ftrue_name[40];
  char qtubemean[MAX];
  double qtubemeane[MAX];
  double qtubedistance[MAX];
  double qctubemean[MAX];
  double truetubemean[MAX];
  double qctubemeane[MAX];
  double qctubeang[MAX];
  double numpmts_corT = 0.;
  double numpmts_trueT = 0;


  int totalBins =0;
  for (int ipmt=0; ipmt<npmts; ++ipmt){
    if (fittingplots == 1){
      WCSimRootPMT pos = geo->GetPMT(ipmt);
      double posX = pos.GetPosition(0);
      double posY = pos.GetPosition(1);
      double posZ = pos.GetPosition(2);
      TVector3 pmtpos(pos.GetPosition(0), pos.GetPosition(1), pos.GetPosition(2));
      TVector3 rprop = pmtpos - photonstartpos;
      double dist = rprop.Mag();
      qtubedistance[ipmt] = dist;
      qtubemean[ipmt] = qperpmt[ipmt]->GetMean();
      qtubemeane[ipmt] = qperpmt[ipmt]->GetMeanError();
      //correct for tube orientation
      TVector3 normal (pos.GetOrientation(0),pos.GetOrientation(1), pos.GetOrientation(2));
      double corfactor = fabs(std::cos(acos(-1.0)-normal.Angle(rprop)));
      qctubemean[ipmt]= qtubemean[ipmt]*dist*dist/100.0/100.0;
      qctubemeane[ipmt]= qtubemeane[ipmt]*dist*dist/100.0/100.0;
      qctubeang[ipmt]=corfactor;

      //Gaussian fit for corrected time per PMT
      sprintf(fcn_name,"gaussian_%05d",ipmt);
      TF1*fcn = new TF1(fcn_name,"gaus", 900., 1000.);
      int maxbin = timeperpmtcor_pr[ipmt]->GetMaximumBin();
      double xofmaxbin = timeperpmtcor_pr[ipmt]->GetXaxis()->GetBinCenter(maxbin);
      double maxval = timeperpmtcor_pr[ipmt]->GetBinContent(maxbin);
      double xmax = xofmaxbin + 3.5;
      double xmin = xofmaxbin - 3.5;
      fcn->SetParameter(0, maxval);
      fcn->SetParameter(1, xofmaxbin);
      fcn->SetParameter(2, 2.0);
      timeperpmtcor_pr[ipmt]->Fit(fcn,"Q","", xmin, xmax);
      //mean = fcn->GetParameter(1)
      //sigma = fcn->GetParameter(2)
      cortubemean[ipmt] = fcn->GetParameter(1)-digi_toffs;
      cortubeerr[ipmt] = fcn->GetParError(1);
      zeroes[ipmt] = 0;
    }//end of if(fittingplots)
    cortubemean[ipmt] = timeperpmtcor_pr[ipmt]->GetMean();
    pmtnum[ipmt] = ipmt;

    /*   //Gaussian fit for true time per PMT
	 sprintf(ftrue_name,"true_gaus_%05d",ipmt);
	 TF1*fcn_true = new TF1(ftrue_name,"gaus", 0., 20.);
	 int truemaxbin = truetimeperpmt[ipmt]->GetMaximumBin();
	 double truexofmaxbin = truetimeperpmt[ipmt]->GetXaxis()->GetBinCenter(truemaxbin);
	 double truemaxval = truetimeperpmt[ipmt]->GetBinContent(truemaxbin);
	 double truexmax = truexofmaxbin - 1.5;
	 double truexmin = truexofmaxbin + 1.5;
	 fcn_true->SetParameter(0, truemaxval);
	 fcn_true->SetParameter(1,truexofmaxbin);
	 fcn_true->SetParameter(2,2.0);
	 truetimeperpmt[ipmt]->Fit(fcn_true,"Q","", truexmin, truexmax);
	 //mean = fcn->GetParameter(1)
	 //rms = fcn->GetParameter(2)
	 truetubemean[ipmt] = fcn_true->GetParameter(1);
    */


    truemaxbin[ipmt]  = truetimeperpmt[ipmt]->GetXaxis()->GetBinCenter(truetimeperpmt[ipmt]->GetMaximumBin());
    truelobound[ipmt] = truemaxbin[ipmt]-2.0;
    truehibound[ipmt] = truemaxbin[ipmt]+2.0;
    /*std::cout <<"truemaxbin:"<<truemaxbin[ipmt]
	      <<"truelobound:"<<truelobound[ipmt]
	      <<"truehibound:"<<truehibound[ipmt]<<std::endl;
    */
    if (truetime1[ipmt]>=truelobound[ipmt] && truetime1[ipmt]<=truehibound[ipmt]){
      hitinPMT_truecut[ipmt]+=1;
    }

    totalBins = truetimeperpmt[ipmt]->GetXaxis()->GetNbins();
    //std::cout<<"totalbins:"<<totalBins<<std::endl;
    for(int i = 1; i<=totalBins; ++i){
      //int initialBin = truetimeperpmt[ipmt]->GetBin(i);
      double bincenter = truetimeperpmt[ipmt]->GetBinCenter(i);
      int bincontent = truetimeperpmt[ipmt]->GetBinContent(i);
      /*std::cout<<"bincenter: "<<bincenter
	<<"bincontent: "<<truetimeperpmt[ipmt]->GetBinContent(i)<<std::endl;*/
      if (bincenter>=truelobound[ipmt] && bincenter<=truehibound[ipmt]){
	totalCount[ipmt] += truetimeperpmt[ipmt]->GetBinContent(i);
	//std::cout<<"totalcount: "<<totalCount<<std::endl;
      }
    }
    //std::cout<<"totalCount: "<<totalCount[ipmt]
    //	     <<" || hitinPMT_truecut[ipmt]"<<hitinPMT_truecut[ipmt]<<std::endl;

    if (verbose) {
      std::cout<<"ipmt = "<<ipmt<<"PMT_toffs[ipmt] = "<<PMT_toffs[ipmt]
	       <<" Time_toffs: "<<time_toffs[ipmt]
	       <<" Cortime: "<<timeperpmtcor_pr[ipmt]->GetMean()-time_toffs[ipmt]
	       <<" corTimeperpmt: "<<timeperpmtcor_pr[ipmt]->GetMean()
	       <<" countEvents: "<<countevperpmt[ipmt]
	       <<" Average Charge:"<<totqperpmt[ipmt]<<std::endl;
    }

    //std::cout<<"PhHit_prop_time: "<<PhHit_prop_time[ipmt]<<std::endl;
    if(hitinPMT_corTprompt[ipmt]>0){
      numpmts_corT+=1;
    }   

    if (hitinPMT_truecut[ipmt]>150){
      numpmts_trueT+=1;
    }
   
    if (perpmtplots == 1){
      q_av[ipmt] = (double)totalq[ipmt]/200.;
      hitprob[ipmt] = (double)countevperpmt[ipmt]/200.;//total 200 events
      eventsperpmt_truecut->Fill(hitinPMT_truecut[ipmt]); 
    }//end of if(perpmtplots)
    
    htcalibfit->Fill(cortubemean[ipmt]-time_toffs[ipmt]); 
    AvChargeVsProbability->Fill(q_av[ipmt], hitprob[ipmt]);
    randtoffs->Fill(time_toffs[ipmt]);
    pmtnum[ipmt] = ipmt;
    RIperpmt->Fill(ipmt, RIperPMT[ipmt]->GetMean());


    percentagehitsinPMTs[ipmt] = totalCount[ipmt]/countTruehitsperpmt[ipmt];
    /*std::cout<<ipmt<<" || percentagehitsinPMTs: "<<percentagehitsinPMTs[ipmt]
	     <<" || countTruehitsperpmtincut[ipmt]: "<<countTruehitsperpmtincut[ipmt]
	     <<" || countTruehitsperpmt[ipmt]: "<<countTruehitsperpmt[ipmt]<<std::endl;
    */

    if (perpmtplots == 1){
      WCSimRootPMT pos = geo->GetPMT(ipmt);
      double posX = pos.GetPosition(0);
      double posY = pos.GetPosition(1);
      double posZ = pos.GetPosition(2);
      TVector3 pmtpos(pos.GetPosition(0), pos.GetPosition(1), pos.GetPosition(2));
      if(posY >=det_ymax-10.0){
	pmtnum_top->Fill(posX, posZ, pmtnum[ipmt]);
	httop_t->Fill(posX, posZ, timeperpmt[ipmt]->GetMean());
	httop_cor->Fill(posX, posZ, correctimeperpmt[ipmt]->GetMean());
	httop_cor1->Fill(posX, posZ, timeperpmtcor_pr[ipmt]->GetMean() - digi_toffs);
	hqav_top->Fill(posX, posZ, qperpmt[ipmt]->GetMean());
	evtopcut_trueT->Fill(posX, posZ, hitinPMT_truecut[ipmt]);
	truettop->Fill(posX, posZ, truetmedian[ipmt]);
	truexacttop->Fill(posX, posZ, PhHit_prop_time[ipmt]);
	trueventtop1->Fill(posX, posZ, countsinpmt1[ipmt]);
	trueventtop2->Fill(posX, posZ, countsinpmt2[ipmt]);
	//truettop_gv->Fill(posX, posZ, cortruetime1[ipmt]);
	percentagehitsinpmt_top->Fill(posX, posZ, percentagehitsinPMTs[ipmt]);

      }
      else if(posY <=det_ymin+10.0){
	pmtnum_bot->Fill(posX, posZ, pmtnum[ipmt]);
	htbot_t->Fill(posX, posZ, timeperpmt[ipmt]->GetMean());
	htbot_cor->Fill(posX, posZ, correctimeperpmt[ipmt]->GetMean());
	htbot_cor1->Fill(posX, posZ, timeperpmtcor_pr[ipmt]->GetMean() - digi_toffs);
	hqav_bot->Fill(posX, posZ, qperpmt[ipmt]->GetMean());
	evbotcut_trueT->Fill(posX, posZ, hitinPMT_truecut[ipmt]);
	truetbot->Fill(posX, posZ, truetmedian[ipmt]);
	truexactbot->Fill(posX, posZ, PhHit_prop_time[ipmt]);
	trueventbot1->Fill(posX, posZ, countsinpmt1[ipmt]);
	trueventbot2->Fill(posX, posZ, countsinpmt2[ipmt]);
	//truetbot_gv->Fill(posX, posZ, cortruetime1[ipmt]);
	percentagehitsinpmt_bot->Fill(posX, posZ, percentagehitsinPMTs[ipmt]);
      }
      else{
	double ang = atan2(posZ, posX);
	pmtnum_side->Fill(ang, posY, pmtnum[ipmt]);
	htside_t->Fill(ang, posY, timeperpmt[ipmt]->GetMean());
	htside_cor->Fill(ang, posY, correctimeperpmt[ipmt]->GetMean());
	htside_cor1->Fill(ang, posY, timeperpmtcor_pr[ipmt]->GetMean() - digi_toffs);
	hqav_side->Fill(ang, posY, qperpmt[ipmt]->GetMean());
	evsidecut_trueT->Fill(ang, posY, hitinPMT_truecut[ipmt]);
	truetside->Fill(ang, posY, truetmedian[ipmt]);
	truexactside->Fill(ang, posY, PhHit_prop_time[ipmt]);
	trueventside1->Fill(ang, posY, countsinpmt1[ipmt]);
	trueventside2->Fill(ang, posY, countsinpmt2[ipmt]);
	//truetside_gv->Fill(ang, posY, cortruetime1[ipmt]);
	percentagehitsinpmt_side->Fill(ang, posY, percentagehitsinPMTs[ipmt]);
      }
    }//end of if(perpmtplots)
  }

 

  //********************************************************************************//
  //**************************Minimizing the parameters*****************************//
  //********************************************************************************//
  std::cout<<"Allocating memory for fit parameters"<<std::endl;
  //2D array for plotting a graph
  //Double_t *mintoffsperpmt = new Double_t[npmts+1][numofit+1];
  double **mintoffsperpmt = new double*[npmts];  
  for(int i=0; i<npmts; ++i){
    mintoffsperpmt[i] = new double[numofit]; 
  }

  const int arrayLength = numofit;
  //Int_t *arrayLength = new Int_t[numofit];
  Double_t *minimizedXpos = new Double_t[numofit];
  Double_t *minimizedYpos = new Double_t[numofit];
  Double_t *minimizedZpos = new Double_t[numofit];
  Double_t *minimizedGroupV = new Double_t[numofit];
  Double_t *minimizedproToffs = new Double_t[numofit];
  Double_t x_arr[numofit];
  
  double lambdaovern_it = wavelength/ref_index;
  double derivative_it = (ref_index1-ref_index)/(wavelength1-wavelength);
  double covern_it = solc/ref_index;
  double groupvelocity = (covern_it*(1.0+(lambdaovern_it*derivative_it)));

  PMTTimeInfo *fti = PMTTimeInfo::GetInstance();
  std::vector<double>minpars;  
  std::vector<double>rantoffs;


  std::cout<<"Finished allocating memory for fit parameters"<<std::endl;
  bool dofit=false;
  if (dofit){

  for (int i=0; i<15807; ++i){
    rantoffs.push_back(time_toffs[i]);
  }
  fti->add_rantime(rantoffs);
  minpars.push_back( sourcepos.X() );
  minpars.push_back( sourcepos.Y() );
  minpars.push_back( sourcepos.Z() );
  minpars.push_back( groupvelocity );
  minpars.push_back( pr_toffs );


  for(int i =0; i<npmts; ++i){
    includePMT.push_back(i);
      WCSimRootPMT pos = geo->GetPMT(i);
      double posX = pos.GetPosition(0);
      double posY = pos.GetPosition(1);
      double posZ = pos.GetPosition(2);
      TVector3 pmtpos(pos.GetPosition(0), pos.GetPosition(1), pos.GetPosition(2));
      posPMT.push_back( pmtpos );
  }

  fti->add_pmts( includePMT );
  fti->add_pmts_pos( posPMT );
  // second loop over events, to fill PMTTimeInfo for the pmts we have selected
  for(int ev = 0; ev<nevent; ev++){
    int itree = tree->LoadTree(ev);
    int nbytes = tree->GetEntry(ev,1);
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    int ncherenkovhits = wcsimrootevent->GetNcherenkovhits();
    TClonesArray* timeArray = wcsimrootevent->GetCherenkovHitTimes();
    for (int i =0; i<ncherenkovhits; i++){
      TObject *Hit = (wcsimrootevent->GetCherenkovHits())->At(i);
      WCSimRootCherenkovHit *wcsimrootcherenkovhit = dynamic_cast<WCSimRootCherenkovHit*>(Hit);
      int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
      int peForTube = wcsimrootcherenkovhit->GetTotalPe(1);
      int tubeID = wcsimrootcherenkovhit->GetTubeID();
      int pmt = tubeID-1;
      for (int j = timeArrayIndex; j<timeArrayIndex+peForTube; j++)
	{
	  WCSimRootCherenkovHitTime *HitTime = dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(j));
   	  if(hitinPMT_corTprompt[pmt]>0){	
	    double t = HitTime->GetTruetime();
	    double timeinpmt =0.;
	    if (t>=truelobound[pmt] && t<=truehibound[pmt])
	      {
		timeinpmt = t; 
	      }
	    fti->add_time( pmt, timeinpmt );
	  }
	}
    }
  }

  double *prop_time_1 = new double[npmts]; 


  //---------------------------------------------------------------//
  //--------------------Iterations begin here----------------------//
  //---------------------------------------------------------------//

  for(int it=0; it<numofit; ++it){
    std::cout<<"\viteration: "<<it+1<<"\v"<<std::endl;

    std::cout<<"minpars[0]: "<<minpars[0]
	     <<" || minpars[1]: "<<minpars[1]
	     <<" || minpars[2]: "<<minpars[2]<<std::endl;
    const TVector3 photonstartposit(minpars[0], minpars[1], minpars[2]);
    for(int ev = 0; ev<nevent; ev++){
      int itree = tree->LoadTree(ev);
      int nbytes = tree->GetEntry(ev,1);
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
      int ncherenkovhits = wcsimrootevent->GetNcherenkovhits();
      TClonesArray* timeArray = wcsimrootevent->GetCherenkovHitTimes();
      for (int i =0; i<ncherenkovhits; i++){
	TObject *Hit = (wcsimrootevent->GetCherenkovHits())->At(i);   

	WCSimRootCherenkovHit *wcsimrootcherenkovhit = dynamic_cast<WCSimRootCherenkovHit*>(Hit);
	int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	int peForTube = wcsimrootcherenkovhit->GetTotalPe(1);
	int tubeID = wcsimrootcherenkovhit->GetTubeID();
	int pmt = tubeID-1;

	WCSimRootPMT pos_1 = geo->GetPMT(pmt);
	TVector3 pmtpos_1(pos_1.GetPosition(0), pos_1.GetPosition(1), pos_1.GetPosition(2));
	TVector3 distdiff_1 = pmtpos_1 - photonstartposit;
	double distMag_1 = distdiff_1.Mag();
	prop_time_1[pmt] = (distMag_1/groupvel);

	for (int j = timeArrayIndex; j<timeArrayIndex+peForTube; j++)
	  {
	    WCSimRootCherenkovHitTime *HitTime = dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(j));
	    double t = HitTime->GetTruetime();
	    double TrueTime =0.;
	    if (t>=truelobound[pmt] && t<=truehibound[pmt])
	      {
		TrueTime = t;
	      }
	    corTruetime_it[it]->Fill(TrueTime - prop_time_1[pmt]);

	  }
      }
    }

    /*double *prop_time_1 = new double[npmts]; 
    const TVector3 photonstartposit(minpars[0], minpars[1], minpars[2]);
    for(int i=0; i<npmts; ++i){
    WCSimRootPMT pos_1 = geo->GetPMT(i);
    TVector3 pmtpos_1(pos_1.GetPosition(0), pos_1.GetPosition(1), pos_1.GetPosition(2));
    TVector3 distdiff_1 = pmtpos_1 - photonstartposit;
    double distMag_1 = distdiff_1.Mag();
    prop_time_1[i] = (distMag_1/groupvel);
    corTruetime_it[it]->Fill(prop_time_1[i]);
    }*/
   
    
    ROOT::Math::Minimizer *min1 =  ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    // min->SetMaxFunctioCalls(1000000);
    //min->SetMaxIterations(100000);
    //min->SetTolerance(0.001);
 
    // calculate number of parameters:
    int position_pars = 3; // three coordinates location of source
    //int ref_index_pars = 1; // refractive index
    int groupV = 1;
    //int gtime_offset_par = 1; // global time offset parameter
    int protimeoffset = 1; //offset in the propagation time of photons from source to PMTs
    int total_num_parameters = position_pars + groupV + protimeoffset;
    ROOT::Math::Functor fcn1(&TCalibChi2, total_num_parameters);
    const std::vector< int > pmts_used = fti->pmts2();
    std::vector< std::string > parnames1;
    std::vector< double > step1;
    std::vector< double > pars1;

    parnames1.push_back( "source_x" ); 
    pars1.push_back( minpars[0] ); 
    step1.push_back(1.0);
  
    parnames1.push_back( "source_y" ); 
    pars1.push_back( minpars[1] );  
    step1.push_back(1.0);
  
    parnames1.push_back( "source_z" ); 
    pars1.push_back( minpars[2] );  
    step1.push_back(1.0);

    //parnames.push_back( "ref_index");
    //pars.push_back( refIndex ); 
    //step.push_back(0.01);

    parnames1.push_back("groupvelocity");
    pars1.push_back( minpars[3] );
    step1.push_back(0.01);

    parnames1.push_back("protimeoffset");
    pars1.push_back(minpars[4]);
    step1.push_back(0.01);

    min1->SetFunction(fcn1);
    int numofpars =pars1.size();
    for (int i=0; i<total_num_parameters; ++i ){
      min1->SetVariable(i, parnames1[i], pars1[i], step1[i] );
      std::cout << "It: " << it << " || param: " << parnames1[i] << " || val: " << pars1[i] << std::endl;
    }

    min1->SetLimitedVariable(0, parnames1[0], pars1[0], step1[0], -370.0,  370.0);
    min1->SetLimitedVariable(1, parnames1[1], pars1[1], step1[1], -50.0,  50.0);
    min1->SetLimitedVariable(2, parnames1[2], pars1[2], step1[2], -25.0, 25.0);
    min1->SetLimitedVariable(3, parnames1[3], pars1[3], step1[3], 20.0, 30.0);
    min1->SetLimitedVariable(4, parnames1[4], pars1[4], step1[4], -10.0, 10.0);
    min1->Minimize();
    //min1->PrintResults();

    //std::cout<<"check1"<<std::endl;

    const double *xarr1 = min1->X();
    minpars.clear();
    for(int i=0; i<5; ++i){
      minpars.push_back(xarr1[i]);
	//std::cout<<"Min val1:"<<xarr1[i]<<std::endl;    
    }

    int minparsize = minpars.size();
    /*for (int i=0; i<minparsize; ++i){
      std::cout<<minpars[i]<<std::endl;
      }*/

    minimizedXpos[it]=(minpars[0]);
    minimizedYpos[it]=(minpars[1]);
    minimizedZpos[it]=(minpars[2]);
    minimizedGroupV[it]=(minpars[3]);
    minimizedproToffs[it]=(minpars[4]);

    std::vector<int>allrandomPMTs;
    int allrandomPMTsize = allrandomPMTs.size();
    std::vector<int>top;
    std::vector<int>bot;
    std::vector<int>side1;
    std::vector<int>midside;
    std::vector<int>side2;
    std::vector<int>randomnumber;
    std::vector<double>mintoffs (15808,0);//a vector holding the minimized value of time offsets
    randomnumber.reserve(10);
    int randomnumbersize =randomnumber.size();
    std::vector<int>randomPMTstop;
    std::vector<int>randomPMTsside1;
    std::vector<int>randomPMTsmidside;
    std::vector<int>randomPMTsside2;
    std::vector<int>randomPMTsbot;
    std::vector<int>rantoffstemp;
    std::vector<int>nexttop;
    std::vector<int>nextbot;
    std::vector<int>nextside1;
    std::vector<int>nextmidside;
    std::vector<int>nextside2;

    //fill in first section usingtop PMTs and PMTs on the side(topmost)
    //3161PMTs in this part of detector
    for(int i=0; i<=576; ++i){
      //if (percentagehitsinPMTs[i]>=0.6 && percentagehitsinPMTs[i]<=1)
	top.push_back(i);
      
    }
    for(int i=13224; i<=15807; ++i){
      //if (percentagehitsinPMTs[i]>=0.6 && percentagehitsinPMTs[i]<=1)
	top.push_back(i);
      
    }
    //fill in second section using the PMTs from the side below the topmost
    //3162PMTs in this section
    for(int i=577; i<=3738; ++i){
      //if (percentagehitsinPMTs[i]>=0.6 && percentagehitsinPMTs[i]<=1)
	side1.push_back(i);
      
    }
    //compose the third sections of the PMTs right from the middle of the side of detector
    //3162PMTs in this section 
    for(int i=3739; i<=6900; ++i){
      //if (percentagehitsinPMTs[i]>=0.6 && percentagehitsinPMTs[i]<=1)
	midside.push_back(i);
      
    }
    //the fourth section is composed of the PMTs from the bottom of detector
    //3162PMTs in this part of detector 
    for(int i=6901; i<=10062; ++i){
      //if (percentagehitsinPMTs[i]>=0.6 && percentagehitsinPMTs[i]<=1)
	side2.push_back(i);
      
    }
    //fifth section includes the PMTs from the bottom most of side of the detector and the bottom itself
    //3161PMTs in this part of detector 
    for(int i=10063; i<=13223; ++i){
      //if (percentagehitsinPMTs[i]>=0.6 && percentagehitsinPMTs[i]<=1)
	bot.push_back(i);
      
    }

    int topsize = top.size();
    int side1size =side1.size();
    int midsidesize = midside.size();
    int side2size = side2.size();
    int botsize = bot.size();
    std::cout<<"topsize:"<<topsize
	     <<" || side1size:"<<side1size
	     <<" || midsidesize:"<<midsidesize
	     <<" || side2size:"<<side2size
	     <<" || botsize:"<<botsize<<std::endl;
    

    //keep ON the loop until it finishes randomizing and using all the PMTs in different sections on a
    //count 2 PMTs from each section.
    double rand =0;  
    int i0=1;
    TRandom3*random = new TRandom3(0);//choose seed by itself
    //loop over choosing random PMTs and minimization begins here 
    while(side2size>2){
      /*std::cout<<"\v"<<"%%%%%%%%||"
	       <<i0<<"||%%%%%%%%"
	       <<"\v"<<std::endl;
      */   
      randomnumber.clear();
      for (int i1 =0; i1<=1; ++i1){
	rand = random->Uniform(topsize);
	int rand1 = static_cast<int>(rand);
	randomnumber.push_back(rand1);
      }
    
      randomnumbersize =randomnumber.size();
      std::sort(randomnumber.begin(),randomnumber.end());
      //erase the repeated elementsds of randomnumber vector
      randomnumber.erase(std::unique(randomnumber.begin(),randomnumber.end()),randomnumber.end());
      randomnumbersize=randomnumber.size();
      //if the randomnumber vector has less than 2 elements,
      //this will be resolved in next loop.      
      //now chosing the next random number
      //check to keep 2 elements in randomnumber vector
      if(randomnumbersize<2){
	//if randomnumber vector has less than 2 elements, jump into this loop
	while(randomnumbersize<2){
	  //std::cout<<"in while loop of randomnumber"<<std::endl;
	  //to find the next random number for making randomnumber vector of 10
	  rand = random->Uniform(topsize);
	  int rand1 = static_cast<int>(rand);
	  randomnumber.insert(randomnumber.begin()+1, rand1);
	  //double check for the random number chosen in a second turn
	  std::sort(randomnumber.begin(),randomnumber.end());
	  randomnumber.erase(std::unique(randomnumber.begin(),randomnumber.end()),randomnumber.end());
	  randomnumbersize=randomnumber.size();
	}//end of while loop over randomnumbersize
      }
    
      randomPMTstop.clear();
      randomPMTsside1.clear();
      randomPMTsmidside.clear();
      randomPMTsside2.clear();
      randomPMTsbot.clear();
      allrandomPMTs.clear();
      int x0=0;
      //std::cout<<"allrandomPMTsize:"<<allrandomPMTsize<<std::endl;
      for (int i2=0; i2<=1; ++i2){
	//std::cout<<"random numbers: "<<randomnumber[i2]<<std::endl;
	x0=randomnumber[i2];
	randomPMTstop.push_back(top[x0]);
	allrandomPMTs.push_back(top[x0]);
      }
      for (int i2=0; i2<=1; ++i2){
	x0=randomnumber[i2];
	allrandomPMTs.push_back(side1[x0]);
	randomPMTsside1.push_back(side1[x0]);
      }  
      for (int i2=0; i2<=1; ++i2){
	x0=randomnumber[i2];
	allrandomPMTs.push_back(midside[x0]);
	randomPMTsmidside.push_back(midside[x0]);
      }
      for (int i2=0; i2<=1; ++i2){
	x0=randomnumber[i2];
	allrandomPMTs.push_back(side2[x0]);
	randomPMTsside2.push_back(side2[x0]);
      }
      for (int i2=0; i2<=1; ++i2){
	x0=randomnumber[i2];
	allrandomPMTs.push_back(bot[x0]);
	randomPMTsbot.push_back(bot[x0]);
      }
      includeranPMT.clear();
      int randomPMTstopsize=randomPMTstop.size();
      //all random PMTs is just to know the total PMTs in vector
      //std::cout<<"***now all in order***"<<std::endl;
      allrandomPMTsize=allrandomPMTs.size();
      for(int j1=0; j1<=allrandomPMTsize-1; ++j1){
	//std::cout<<j1<<":-"<<allrandomPMTs[j1]<<std::endl;
	includeranPMT.push_back(allrandomPMTs[j1]);
      }
      fti->add_ranpmts(allrandomPMTs);

      ROOT::Math::Minimizer *min2 =  ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
      // min->SetMaxFunctioCalls(1000000);
      //min->SetMaxIterations(100000);
      //min->SetTolerance(0.001);
     
      // calculate number of parameters:
      int position_pars2 = 3; // three coordinates location of source
      int groupV2 = 1;
      int protimeoffset2 = 1; //offset in the propagation time of photons from source to PMTs
      int PMTs = allrandomPMTs.size();
      int total_num_parameters2 = position_pars2 + groupV2 + protimeoffset2 + PMTs;
      ROOT::Math::Functor fcn2(&TCalibChi2toffs, total_num_parameters2);
      const std::vector< int > pmts_used2 = fti->ranpmts2();
      std::vector< std::string > parnames2;
      std::vector< double > step2;
      std::vector< double > pars2;

      //std::cout<<"allrandomPMTsize:"<<PMTs<<std::endl;
      for(int i=0; i<PMTs; ++i){
	int pmt = allrandomPMTs[i];
	std::ostringstream os2;
	os2<<"PMT("<<pmt<<")_toffs";
	parnames2.push_back(os2.str());
	pars2.push_back(rantoffs[pmt]);
	step2.push_back(0.001);
      }

      parnames2.push_back( "source_x" ); 
      pars2.push_back( minpars[0]); 
      step2.push_back(1.0);
  
      parnames2.push_back( "source_y" ); 
      pars2.push_back( minpars[1] );  
      step2.push_back(1.0);
  
      parnames2.push_back( "source_z" ); 
      pars2.push_back( minpars[2] );  
      step2.push_back(1.0);

      parnames2.push_back("groupvelocity");
      pars2.push_back( minpars[3] );
      step2.push_back(0.01);

      parnames2.push_back("protimeoffset");
      pars2.push_back( minpars[4] );
      step2.push_back(0.01);

      min2->SetFunction(fcn2);
      for (int i=0; i<PMTs; ++i){
	min2->SetLimitedVariable(i, parnames2[i], pars2[i], step2[i], -10.0,  10.0);
	//std::cout << "It: " << it << " i: " << i << " param: " << parnames2[i] << " val: " << pars2[i] << std::endl;
      }
      for (int i=PMTs; i<=PMTs+4; ++i ){
	min2->SetVariable(i, parnames2[i], pars2[i], step2[i] );
	min2->FixVariable(i);
	//std::cout << "It: " << it << " i: " << i << " param: " << parnames2[i] << " val: " << pars2[i] << std::endl;
      }
      min2->Minimize();
      //min2->PrintResults();
      const double *xarr2 = min2->X();
      for(int i=0; i<PMTs; ++i){
	mintoffs[allrandomPMTs[i]] = (xarr2[i]);
	//std::cout<<"Min val2:"<<mintoffs[i]<<"::"<<xarr2[i]<<std::endl;    
      }
      std::sort(randomPMTstop.begin(),randomPMTstop.end());
      //clear up the already chosen PMTs 
      //so that these are not repeated in the second selection
      nexttop.reserve(topsize-2);
      nextside1.reserve(side1size-2);
      nextmidside.reserve(midsidesize-2);
      nextside2.reserve(side2size-2);
      nextbot.reserve(botsize-2);
      nexttop.clear();
      nextside1.clear();
      nextmidside.clear();
      nextside2.clear();
      nextbot.clear();
      std::set_difference(top.begin(), top.end(), randomPMTstop.begin(), randomPMTstop.end(), std::back_inserter(nexttop));
      std::set_difference(side1.begin(), side1.end(), randomPMTsside1.begin(), randomPMTsside1.end(), std::back_inserter(nextside1));
      std::set_difference(midside.begin(), midside.end(), randomPMTsmidside.begin(), randomPMTsmidside.end(), std::back_inserter(nextmidside));
      std::set_difference(side2.begin(), side2.end(), randomPMTsside2.begin(), randomPMTsside2.end(), std::back_inserter(nextside2));
      std::set_difference(bot.begin(), bot.end(), randomPMTsbot.begin(), randomPMTsbot.end(), std::back_inserter(nextbot));

      //std::cout<<"check1"<<std::endl;
      int nextbotsize =nextbot.size();
      int nextside1size = nextside1.size();
      int nextmidsidesize = nextmidside.size();
      int nextside2size =nextside2.size();
      int nexttopsize = nexttop.size();  
  
      /*std::cout<<"nexttopsize"<<nexttopsize
	       <<" || nextside1size"<<nextside1size
	       <<" || nextmidsidesize"<<nextmidsidesize
	       <<" || nextside2size:"<<nextside2size
	       <<" || nextbotsize:"<<nextbotsize<<std::endl;
      */      
      top.clear();
      side1.clear();
      midside.clear();
      side2.clear();
      bot.clear();

      for(int i3=0; i3<=nexttopsize-1; ++i3){
	top.push_back(nexttop[i3]);
	bot.push_back(nextbot[i3]);
      }
      //std::cout<<"check3"<<std::endl;
      for(int i3=0; i3<=nextside2size-1; ++i3){
	side1.push_back(nextside1[i3]);
	midside.push_back(nextmidside[i3]);
	side2.push_back(nextside2[i3]);      
      }  

      i0+=1;
      topsize=top.size();
      side1size=side1.size();
      midsidesize=midside.size();
      side2size=side2.size();
      botsize=bot.size();
      min2->Clear();
      //std::cout<<"~!# in while loop #!~"<<std::endl;
    }//end of while(side2size==2) loop***************

    std::cout<<"~!# out of while loop #!~"<<std::endl;
    std::vector<int>remainingPMTs;
    for(int i=0; i<topsize; ++i){
      remainingPMTs.push_back(top[i]);
    }

    for(int i=0; i<topsize; ++i){
      remainingPMTs.push_back(bot[i]);
    }

    for(int i=0; i<side1size; ++i){
      remainingPMTs.push_back(side1[i]);
    }

    for(int i=0; i<midsidesize; ++i){
      remainingPMTs.push_back(midside[i]);
    }

    for(int i=0; i<side2size; ++i){
      remainingPMTs.push_back(side2[i]);
    }

    fti->add_ranpmts(remainingPMTs);

    //std::cout<<"remaining PMTs:"<<remainingPMTs.size()<<std::endl;
    ROOT::Math::Minimizer *min3 =  ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    // min->SetMaxFunctioCalls(1000000);
    //min->SetMaxIterations(100000);
      //min->SetTolerance(0.001);
     
      // calculate number of parameters:
      int position_pars3 = 3; // three coordinates location of source
      int groupV3 = 1;
      int protimeoffset3 = 1; //offset in the propagation time of photons from source to PMTs
      int remainPMTs = remainingPMTs.size();
      int total_num_parameters3 = position_pars3 + groupV3 + protimeoffset3 + remainPMTs;
      ROOT::Math::Functor fcn3(&TCalibChi2toffs, total_num_parameters3);
      const std::vector< int > pmts_used3 = fti->ranpmts2();
      std::vector< std::string > parnames3;
      std::vector< double > step3;
      std::vector< double > pars3;

      for(int i=0; i<remainPMTs; ++i){
	std::ostringstream os3;
	os3<<"PMT("<<i<<")_toffs";
	parnames3.push_back(os3.str());
	pars3.push_back(rantoffs[remainingPMTs[i]]);
	step3.push_back(0.001);
      }

      parnames3.push_back( "source_x" ); 
      pars3.push_back( minpars[0]); 
      step3.push_back(1.0);
  
      parnames3.push_back( "source_y" ); 
      pars3.push_back( minpars[1] );  
      step3.push_back(1.0);
  
      parnames3.push_back( "source_z" ); 
      pars3.push_back( minpars[2] );  
      step3.push_back(1.0);

      parnames3.push_back("groupvelocity");
      pars3.push_back( minpars[3] );
      step3.push_back(0.01);

      parnames3.push_back("protimeoffset");
      pars3.push_back( minpars[4] );
      step3.push_back(0.01);

      min3->SetFunction(fcn3);
      for (int i=0; i<remainPMTs; ++i){
	min3->SetLimitedVariable(i, parnames3[i], pars3[i], step3[i], -10.0,  10.0);
      }
      for (int i=remainPMTs; i<=remainPMTs+4; ++i ){
	//std::cout<<"checknn"<<std::endl;
	min3->SetVariable(i, parnames3[i], pars3[i], step3[i] );
	min3->FixVariable(i);
      }

      min3->Minimize();
      //min3->PrintResults();
      const double *xarr3 = min3->X();
      for(int i=0; i<remainPMTs; ++i){
	mintoffs[remainingPMTs[i]]=(xarr3[i]);
	//std::cout<<"Min val3:"<<mintoffs[i]<<"::"<<xarr3[i]<<std::endl;    
      }

      int mintoffsize = mintoffs.size();
      rantoffs.clear();
      //std::sort(mintoffs.begin(), mintoffs.end());
      std::cout<<"mintoffsize:"<<mintoffsize<<std::endl;
      x_arr[it] =it;//assuming the minimization is on 0.01 scale over each iteration
      for (int i=0; i<mintoffsize; ++i){
	minRanTOffs[it]->Fill(mintoffs[i]);
	if (it == 0){
	minrandtoffs->Fill(mintoffs[i]);
	}
	rantoffs.push_back(mintoffs[i]);
	//mintoffsperpmt[i][it] = mintoffs[i];
	//std::cout<<i<<": "<<mintoffs[i]<<std::endl;
      }  
     fti->add_rantime(mintoffs);
     mintoffs.clear();
  }// end of loop over iterations - it
  
  }// end loop that decides whether to run above fit code  
  /* const std::vector<int> & pmts =  fti->pmts2();
  const std::vector<TVector3> & pmtspos = fti->pmtpos2();
  
  for (int i = -370; i <= 370; i++){
    double value;
    if (i%10==0){
      double sx = static_cast<double>(i);
      //double sx = 150.0;
      double sy = 0.0;// pars1[ 1 ];
      double sz = 0.0;//pars1[ 2 ];
      double n  = pars1[ 3 ];
      double chi2val = 0.0;
      double chi2tot = 0.0;
      double initprotime = 0.;
      int noPMTs = pmts.size();
     
      //std::cout<<"no.ofPMTs: "<<noPMTs<<std::endl;
      for (int ipmt = 0; ipmt<noPMTs; ++ipmt){
	double meandt = 0.;
	double sumdt = 0.;
	double middt = 0.;
	double meandt1 = 0.;
	double sumdt1 = 0.;
	double middt1 = 0.;
	int numoftruetimes = 0;
	int numofcortruetimes = 0;
	int pmt = pmts[ ipmt ];
	double noverc = ref_index*oneoverc;//phase velocity used for calculation of propagation time
	double covern = solc/ref_index;
	double meanref_index = (ref_index + ref_index1)/2;//mean of ref_index for group velocity
	double energyderivative = (ref_index1-ref_index)/std::log(energy1/energy); 
	double lambdaovern = wavelength/ref_index;
	double derivative = (ref_index1-ref_index)/(wavelength1-wavelength);
	//double derivative=((ref_index1-ref_index)/(energy-energy1))*energy*energy1;
	double groupvel = (covern*(1.0+(lambdaovern*derivative)));//group velocity for calculation of propagation time
	double groupvelene = solc/(meanref_index+energyderivative);
	groupvelene = 21.70855;
	double protime =0.;//distance of propagation
	double protime1 =0.;//distance of propagation using group velocity
	double proptime_ene = 0.;
	double relprotime = 0.; //relative propagation time
	//std::cout<<"check2"<<std::endl;
	const TVector3 & rpmt = pmtspos[ipmt];
	double dist2 = 
	  (rpmt.X() - sx)*(rpmt.X() - sx) +
	  (rpmt.Y() - sy)*(rpmt.Y() - sy) +
	  (rpmt.Z() - sz)*(rpmt.Z() - sz) ;
	double dist = std::sqrt( dist2 );
	//std::cout<<"check3"<<std::endl;
	protime = dist * n * oneoverc;
	protime1 = dist/groupvel;
	proptime_ene = dist/groupvelocity;
	//std::vector<double>truetimes;
	std::vector<double>cortruetimes;
	const std::vector<double>& times = fti->times( pmt );
	int numoftimes = times.size();
	for ( int itime = 0; itime < numoftimes; ++itime ){
	  if( times[itime] < 0.01 ) continue;
	  chi2truetime->Fill(times[itime]);
	  double dt = ( times[itime] - protime );
	  double dt1 = (times[itime] - protime1);
	  double dt_ene =(times[itime] - proptime_ene);
	  //truetimes.push_back(dt);
	  if (i==150.){
	    cortruetimes.push_back(dt1);
	    sumdt1 +=dt1;
	  }
	  sumdt += dt;
	  chi2val = dt_ene*dt_ene;
	  cortimeperpos->Fill(i, dt_ene);
	  chi2tot += chi2val;
	  chisquarefit->Fill(chi2val); 
	  cortime2perpos->Fill(i, chi2val); 
	}
	//meandt1 = sumdt/numofcortruetimes;
	if(i==150.){
	  std::sort(cortruetimes.begin(), cortruetimes.end());
	  numofcortruetimes = cortruetimes.size();
	  //std::sort(truetimes.begin(), truetimes.end());
	  //	numoftruetimes = truetimes.size();
	  //std::cout<<ipmt<<" TrueTime vector"<<numofcortruetimes<<std::endl;
	  //if (numoftruetimes%2!=0)
	  //middt = truetimes[(numoftruetimes/2)+1];
	  //std::cout<<"middt"<<middt<<std::endl;
	  //
	  //else if(numoftruetimes%2==0)
	  //middt = (truetimes[(numoftruetimes/2)]+truetimes[(numoftruetimes/2)+1])/2;
	  //std::cout<<"middt"<<middt<<std::endl;
	  //
	  if (numofcortruetimes%2!=0){
	    middt1 = cortruetimes[(numofcortruetimes/2)+1];
	  }
	  else if(numofcortruetimes%2==0){
	    middt1 = (cortruetimes[(numofcortruetimes/2)]+cortruetimes[(numofcortruetimes/2)+1])/2;
	  }

	  WCSimRootPMT pos = geo->GetPMT(ipmt);
	  double posX = pos.GetPosition(0);
	  double posY = pos.GetPosition(1);
	  double posZ = pos.GetPosition(2);
	  TVector3 pmtpos(pos.GetPosition(0), pos.GetPosition(1), pos.GetPosition(2));
	  if(posY >=det_ymax-10.0){
	    //meandttop[150]->Fill(posX, posZ, middt);
	    truettop_gv->Fill(posX, posZ, middt1);

	  }
	  else if(posY <=det_ymin+10.0){
	    //meandtbot[150]->Fill(posX, posZ, middt);
	    truetbot_gv->Fill(posX, posZ, middt1);

	  }
	  else{
	    double ang = atan2(posZ, posX);
	    //meandtside[150]->Fill(ang, posY, middt);
	    truetside_gv->Fill(ang, posY, middt1);

	  }
	  //std::cout<<"ipmt:"<<pmts[ipmt]<<" chi2valperpmt:"<<chi2val<<std::endl;
	  //std::cout<<ipmt<<" Propagation time:"<<propdist<<std::endl;
	  //std::cout << "i = " << i << std::endl;
	}
      }
      value = chi2tot/1e+06;
      std::cout<<i<<"\tchi2: "<<chi2tot<<std::endl;
      chi2perpos->Fill(i, value);
    }
  }*/

  /* for (int i=0; i<npmts; ++i){
    for(int j=0; j<numofit; ++j){
      std::cout<<" "<<mintoffsperpmt[i][j];
    }
    std::cout<<"\n"<<std::endl;
    }*/

  /* TGraph **minToffsperPMT;
  TDirectory *curdir6 = Plots->GetDirectory("/");
  TDirectory *minToffsplotdir = Plots->mkdir("minToffsperPMT");
  minToffsplotdir->cd();
  minToffsperPMT = new TGraph*[npmts];
  for (int ipmt = 0; ipmt<npmts; ++ipmt){
    if(ipmt%200==0){
      sprintf(pmtdir, "%05d", ipmt);
      (minToffsplotdir->mkdir(pmtdir))->cd();
    }
    minToffsperPMT[ipmt] = new TGraph(arrayLength, x_arr, mintoffsperpmt[ipmt]);
    minToffsperPMT[ipmt]->SetTitle("minimized Random Timeoffset");
    minToffsperPMT[ipmt]->SetMinimum(-1.0);
    minToffsperPMT[ipmt]->SetMaximum(2.0);
    minToffsperPMT[ipmt]->Write();
  }
  curdir6->cd();
  */
  /*
  TGraph *minXpos = new TGraph(arrayLength, x_arr, minimizedXpos);
  TGraph *minYpos = new TGraph(arrayLength, x_arr, minimizedYpos);
  TGraph *minZpos = new TGraph(arrayLength, x_arr, minimizedZpos);
  TGraph *minGroupV = new TGraph(arrayLength, x_arr, minimizedGroupV);
  TGraph *minproToffs = new TGraph(arrayLength, x_arr, minimizedproToffs);
 
  minXpos->SetTitle("minimized X-position");
  minXpos->SetMinimum(145.0);
  minXpos->SetMaximum(155.0);
  minXpos->Write();

  minYpos->SetTitle("minimized Y-position");
  minYpos->SetMinimum(-1.0);
  minYpos->SetMaximum(2.0);
  minYpos->Write();

  minZpos->SetTitle("minimized Z-position");
  minZpos->SetMinimum(-1.0);
  minZpos->SetMaximum(2.0);
  minZpos->Write();

  minGroupV->SetTitle("minimized Group Velocity");
  minGroupV->SetMinimum(15.0);
  minGroupV->SetMaximum(25.0);
  minGroupV->Write();

  minproToffs->SetTitle("minimized propagation Time-offset");
  minproToffs->SetMinimum(-1.0);
  minproToffs->SetMaximum(2.0);
  minproToffs->Write();
  
  htside_t->SetMinimum(900.);
  httop_t->SetMinimum(900.);
  htbot_t->SetMinimum(900.);
  truetside_gv->SetMinimum(0.);
  truettop_gv->SetMinimum(0.);
  truetbot_gv->SetMinimum(0.);
  */
  Plots->Write();
  Plots->Close();
    
  std::cout<<"num_trig "<<num_trig<<std::endl;

}
  
