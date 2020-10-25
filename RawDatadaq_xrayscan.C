#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <pmonitor/pmonitor.h>
#include <TTree.h>
#include "TCanvas.h"
#include <Event/Event.h>
#include <Event/EventTypes.h> 
#include <TProfile.h>
#include "TF1.h"
#include <iostream>
#include <fstream>

using namespace std;

int nEvents = 0;

//==========================
// GENERAL CONTROL VARIABLES
const int   kControl_NumberOfChannels = 512; // 128*4 (256 + 256)
const int   kControl_NumberofSamples  = 30;  // 30 for SRS
const int   kControl_MaxADCvalue      = 4096;// 4096
const float kControl_NoiseCut = 80.0; // value used to decide if channel was hit or not (centered RMS = StdDev)
const float kControl_RejectFlag1 = 800.; // hardwire-value for flaging bad events (type1) check code below
const int   kControl_SpyChannel = 113; //113;//95; // spy channel for wave inspection
int   kControl_EventNumber = 0;
bool  kControl_Sgn[kControl_NumberOfChannels]; // Channel shows activity?
TProfile *kControl_QA_RAWMEANpc;
TH2F     *kControl_QA_RAWRMSpc;
TH1F     *kControl_QA_RAWWAVEhits;
TH1F     *kControl_QA_RAWWAVEnoise;
TH2F     *kControl_QA_RAWWAVEspy;

//==========================
// PEDESTAL CONTROL VARIABLES
const int   kPedestal_NumberOfSamplesForFastPedestal = 5;
float kPedestal_MEAN[kControl_NumberOfChannels];
const bool  kPedestal_QA = true;
  TProfile *kPedestal_QA_MEANpc;

//==========================
// HIT CONTROL VARIABLES
float kHit_PI[kControl_NumberOfChannels];
float kHit_PH[kControl_NumberOfChannels];
float kHit_PHT[kControl_NumberOfChannels];
float kHit_PHA[kControl_NumberOfChannels];
float kHit_PHAT[kControl_NumberOfChannels];
const bool kHit_QA = true;
  TH1F *kHit_QA_PI;
  TH1F *kHit_QA_PH;
  TH1F *kHit_QA_PHT;
  TH2F *kHit_QA_PHPHT;
  TH1F *kHit_QA_PHA;
  TH1F *kHit_QA_PHAT;
  TProfile *kHit_QA_CHANNELMAP;

//==========================
// CLUSTER CONTROL VARIABLES 
TH1F *kCluster_QX;
TH1F *kCluster_QY;
TH1F *kCluster_XCENTROID;
TH1F *kCluster_YCENTROID;
TH1F *kCluster_XSTD_DEV;
TH1F *kCluster_YSTD_DEV;
TH2F *kCluster_XYCENTROID;

int pinit() {
  for (int i=0; i<kControl_NumberOfChannels; i++) {
    kPedestal_MEAN[i]=0;
    kHit_PH[i]=0;
    kHit_PHT[i]=0;
    kControl_Sgn[i]=0;
    kHit_PHA[i]=0;
    kHit_PHAT[i]=0;
    kHit_PI[i]=0;
  }
  return 0;
}

int phsave(char const*) {
  return 0;
}

int pcontrol() {
  return 0;
}

int pcontrol(int) {
  return 0;
}

void CreateHistograms() {
  // creates new histogram objects
  std::cout << "CreateHistograms INIT" << std::endl;

  kControl_QA_RAWMEANpc = new TProfile( "kControl_QA_RAWMEANpc", "kControl_QA_RAWMEANpc;channel;meanRawADC",
					kControl_NumberOfChannels, -0.5, kControl_NumberOfChannels-0.5 );
  kControl_QA_RAWRMSpc  = new TH2F( "kControl_QA_RAWRMSpc",  "kControl_QA_RAWRMSpc;channel;stddevRawADC",
				    kControl_NumberOfChannels, -0.5, kControl_NumberOfChannels-0.5, 100, 0, 1000 );
  kControl_QA_RAWWAVEhits  = new TH1F( "kControl_QA_RAWWAVEhits",  "kControl_QA_RAWWAVEhits;sample",  kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );  
  kControl_QA_RAWWAVEnoise = new TH1F( "kControl_QA_RAWWAVEnoise", "kControl_QA_RAWWAVEnoise;sample", kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );  
  kControl_QA_RAWWAVEspy   = new TH2F( "kControl_QA_RAWWAVEspy", Form("kControl_QA_RAWWAVEspy chn=%d;event id;sample",kControl_SpyChannel),
				       100, -0.5, 99.5, kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );
  
  if(kPedestal_QA) { // only when QA was requested because it is expensive
    kPedestal_QA_MEANpc    = new TProfile( "kPedestal_QA_MEANpc",    "kPedestal_QA_MEANpc;channel;meanADC",
					   kControl_NumberOfChannels, -0.5, kControl_NumberOfChannels-0.5 );
  }

  if(kHit_QA) {
    kHit_QA_PI   = new TH1F( "kHit_QA_PI",   "kHit_QA_PI;pulse integral", 400, -100, 100 );
    kHit_QA_PH   = new TH1F( "kHit_QA_PH",   "kHit_QA_PH;pulse MAX height", 400, -100, 2000 );
    kHit_QA_PHT  = new TH1F( "kHit_QA_PHT",  "kHit_QA_PHT;pulse MAX height trigger", kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );
    kHit_QA_PHPHT= new TH2F( "kHit_QA_PHPHT","kHit_QA_PHPHT;pulse MAX height;pulse MAX height trigger",
			     100,-100,2000,kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5);
    kHit_QA_PHA  = new TH1F( "kHit_QA_PHA",  "kHit_QA_PHA;pulse AVG height", 400, -100, 2000 );
    kHit_QA_PHAT = new TH1F( "kHit_QA_PHAT", "kHit_QA_PHAT;pulse AVG height trigger", kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );
    kHit_QA_CHANNELMAP = new TProfile( "kHit_QA_CHANNELMAP", "kHit_QA_CHANNELMAP;channel", kControl_NumberOfChannels, -0.5, kControl_NumberOfChannels-0.5 );
  }

  kCluster_QX = new TH1F("kCluster_QX","kCluster_QX;sumQx;events",1000,0,30000);
  kCluster_QY = new TH1F("kCluster_QY","kCluster_QY;sumQy;events",1000,0,30000);
  kCluster_XCENTROID = new TH1F("kCluster_XCENTROID","kCluster_XCENTROID;mm;events",5000,0,100);
  kCluster_YCENTROID = new TH1F("kCluster_YCENTROID","kCluster_YCENTROID;mm;events",1000,0,100);
  kCluster_XSTD_DEV = new TH1F("kCluster_XSTD_DEV","kCluster_XSTD_DEV;mm;events",1000,0,10);
  kCluster_YSTD_DEV = new TH1F("kCluster_YSTD_DEV","kCluster_YSTD_DEV;mm;events",1000,0,10);
  kCluster_XYCENTROID = new TH2F("kCluster_XYCENTROID","kCluster_XYCENTROID;mm;mm",200,0,100,200,0,100);
}
//========================================================================
int check_signals(Packet *p) {
  int flag = 0; // good event
  int nsamp =  p->iValue(0,"NSAMPLES");
  for(int idx=0; idx!=kControl_NumberOfChannels; ++idx) { // loop over all channels
    kControl_Sgn[ idx ] = false; // clear
    int c = idx%128; // unpacking particulars of SRS: channel
    int h = idx/128; // unpacking particulars of SRS: hybrid
    // computing MEAN and STD-DEV
    double sx = 0;
    double sxx = 0;
    for(int i = 0; i < nsamp; i++) { // loop over samples
      float raw = p->iValue ( c, i, 2*h ); // unpack data
      sx += raw;
      sxx += raw*raw;
    }
    double ped_mean = sx/nsamp;
    double ped_rms = sqrt( sxx/nsamp - ped_mean*ped_mean );
    if(ped_rms>kControl_RejectFlag1) flag = 1; // flag event that contains a big spread in any channel
    kControl_QA_RAWMEANpc->Fill( idx, ped_mean );
    kControl_QA_RAWRMSpc->Fill(  idx, ped_rms );
    // if STD-DEV is higher than kControl_NoiseCut then declare activity in this channel
    kControl_Sgn[ idx ] = ped_rms > kControl_NoiseCut;
    for(int i = 0; i < nsamp; i++) { // loop over samples
      float raw = p->iValue ( c, i, 2*h ); // unpack data
      if( kControl_Sgn[ idx ] ) kControl_QA_RAWWAVEhits->Fill( i, raw );
      else kControl_QA_RAWWAVEnoise->Fill( i, raw );
      if(idx==kControl_SpyChannel)
	kControl_QA_RAWWAVEspy->Fill( kControl_EventNumber, i, raw );
    }
  }
  return flag;
}
//========================================================================
void pedestal_calculation(Packet *p) {
  // The following strategy comes from Bob:
  // For each channel we compute the MEAN and RMS of the adc values for all samples
  // based on the StandardDeviation (RMS for TH1F in old code), we decide if the channel was hit or not.
  // the pedestal is the MEAN when there is no hit
  // when there is a hit, the approach is different
  int nsamp =  p->iValue(0,"NSAMPLES");
  for(int idx=0; idx!=kControl_NumberOfChannels; ++idx) { // loop over all channels
    int c = idx%128; // unpacking particulars of SRS: channel
    int h = idx/128; // unpacking particulars of SRS: hybrid
    // computing MEAN and RMS
    double sx = 0;
    // if there is no activity in this channel compute mean over all samples
    // otherwise just use the mean over few samples around the highest fluctuation
    int highest = 0;
    int imax = 0;
    for(int i = 0; i < nsamp; i++) { // loop over samples
      float raw = p->iValue ( c, i, 2*h ); // unpack data
      sx += raw;
      if(highest < raw) {
	highest = raw;
	imax = i; //saves the positions of the highest fluctuation
      }
    }
    double ped_mean = sx/nsamp;
    if(kControl_Sgn[idx]) { // there is activity in this cell so we take the mean over few samples around the highest fluctuation
      sx = 0;
      int ilow = imax;
      if(ilow>nsamp-kPedestal_NumberOfSamplesForFastPedestal) ilow = nsamp-kPedestal_NumberOfSamplesForFastPedestal;
      for(int i = ilow; i < ilow+kPedestal_NumberOfSamplesForFastPedestal; i++) {
	float raw = p->iValue ( c, i, 2*h ); // unpack data
	sx += raw;
      }
      double ped_mean = sx/kPedestal_NumberOfSamplesForFastPedestal;
    }
    kPedestal_MEAN[ h*128 + c ] = ped_mean;
    if(kPedestal_QA) {
      kPedestal_QA_MEANpc->Fill( idx, ped_mean );
    }
  }
}
//========================================================================
void hit_calculation(Packet *p) {
  int nsamp =  p->iValue(0,"NSAMPLES");

  for(int idx=0; idx!=kControl_NumberOfChannels; ++idx) { // loop over all channels
    kHit_PH[idx] = 0;  // pulse amplitude
    kHit_PHT[idx] = 0; // pulse trigger time
    kHit_PHA[idx] = 0; // pulse amplitude wrt three samples
    kHit_PHAT[idx] = 0;// pulse trigger time wrt three samples
    kHit_PI[idx] = 0;  // pulse integral over full window
    if(!kControl_Sgn[idx]) continue; // discard the channels with no activity above noise
    int c = idx%128; // unpacking particulars of SRS: channel
    int h = idx/128; // unpacking particulars of SRS: hybrid
    double ped = kPedestal_MEAN[idx];
    float minRAW = 1E+37; // initialize to a riduculus high value
    int whenHappened = 0;
    for(int i = 0; i < nsamp; i++) { // loop over samples
      float raw = p->iValue ( c, i, 2*h ); // unpack data
      kHit_PI[idx] += (ped - raw);
      if(raw<minRAW) {
	minRAW = raw;
	whenHappened = i;
      }
    }
    kHit_PH[idx] = ped - minRAW;
    kHit_PHT[idx] = whenHappened;
    int iSample = whenHappened;
    if(iSample < 1 || iSample > nsamp - 1 )
      continue; // discard wave where maximum happens to close to edges
    // computing height as average of three samples around peak
    float val_nm1 = ped - p->iValue ( c, iSample-1, 2*h ); // unpack data
    float val_n   = ped - p->iValue ( c, iSample,   2*h ); // unpack data
    float val_np1 = ped - p->iValue ( c, iSample+1, 2*h ); // unpack data
    kHit_PHA[idx] = val_nm1 + val_n + val_np1;
    kHit_PHAT[idx] = ((iSample-1)*val_nm1 + iSample*val_n + (iSample+1)*val_np1) / kHit_PHA[idx];
    kHit_PHA[idx] /= 3;
    if(kHit_QA) {
      kHit_QA_PI->Fill( kHit_PI[idx] );
      kHit_QA_PH->Fill( kHit_PH[idx] );
      kHit_QA_PHA->Fill( kHit_PHA[idx] );
      kHit_QA_PHT->Fill( kHit_PHT[idx] );
      kHit_QA_PHPHT->Fill( kHit_PH[idx], kHit_PHT[idx] );
      kHit_QA_PHAT->Fill( kHit_PHAT[idx] );
      kHit_QA_CHANNELMAP->Fill( idx, kHit_PH[idx] );
    }
  }
}
//========================================================================
void cluster_calculation(Packet *p) {
  if(kControl_NumberOfChannels!=512) return;
  int map_x[256], map_y[256];
  for(int i=0; i!=256; ++i) {
    map_x[i] = i*0.4 - 0.2;
    map_y[i] = i*0.4 - 0.2;
  }
  double sqx_xx=0;
  double sqx_x=0;
  double sqx=0;
  double sqy_yy=0;
  double sqy_y=0;
  double sqy=0;
  for(int idx=0; idx!=256; ++idx) { // for COMPAS half the channels is one detector
    if(kControl_Sgn[idx] ) {
      sqx_xx+= kHit_PH[idx] * map_x[idx] * map_x[idx]; // chargeX*posX*posX
      sqx_x += kHit_PH[idx] * map_x[idx]; // chargeX*posX
      sqx   += kHit_PH[idx]; // chargeX
    }
    if(kControl_Sgn[idx+256]) {    
      sqy_yy+= kHit_PH[idx+256] * map_y[idx] * map_y[idx]; // chargeY*posY*posY
      sqy_y += kHit_PH[idx+256] * map_y[idx]; // chargeY*posY
      sqy   += kHit_PH[idx+256]; // chargeY
    }
  }
  double mX = sqx_x / sqx;
  double mY = sqy_y / sqy;
  double mXX = sqx_xx / sqx;
  double mYY = sqy_yy / sqy;
  double sX = sqrt( mXX-mX*mX );
  double sY = sqrt( mYY-mY*mY );
  //if((mX>39 &&mX<49)&& sX>0.5 &&sX<3.0){ //Vero
  if(sX>0.5 &&sX<3.0){ //Vero
    kCluster_QX->Fill( sqx );
    kCluster_QY->Fill( sqy );
    kCluster_XCENTROID->Fill( mX );
    kCluster_YCENTROID->Fill( mY );
    kCluster_XSTD_DEV->Fill( sX );
    kCluster_YSTD_DEV->Fill( sY );
    kCluster_XYCENTROID->Fill( mX, mY);
  }
}
/*
//=======================================================================
void resolution_calculation(){
  int peakbin = kCluster_XCENTROID->GetMaximumBin();
  float peakpos=kCluster_XCENTROID->GetBinCenter(peakbin);

  TF1 *f1 = new TF1("f1", "gaus");
  f1->SetParLimits(1, peakpos-1.5,peakpos+1.5 ); // search within 1.5 mm
  f1->SetParLimits(2, 0.1 - 0.05, 0.1 + 0.05);
  float minnn = peakpos - 0.5; //qx->GetMean()-1.5;
  float maxxx = peakpos + 0.5; //qx->GetMean()+1.5;
  kCluster_XCENTROID->GetXaxis()->SetRangeUser(minnn, maxxx);
  kCluster_XCENTROID->Fit("f1", "0");
  kResidual_MeanX = f1->GetParameter(1);
  kResidual_MeanX_er=f1->GetParError(1);
  kResidual_SigmaX = f1->GetParameter(2);
  kResidual_SigmaX_er=f1->GetParError(2);
  delete f1;
  
}
//========================================================================
void commom_mode_calculation(Packet *p){
  //calculating pulse height
  pulse_height(p);
  //calculating CM correction curve
  int Ns=0;
  float ph_phtm1=0;
  float ph_pht=0;
  
  for (int h = 0; h< 4; h++)
    {
      Ns=0;
      for (int i = 0; i < p->iValue(0,"NSAMPLES"); i++)
	{
	  if (Ntb[h][i] != 0) ph_tb_av[h][i]=ph_tb_av[h][i]/Ntb[h][i];

	  Ns++;
	  ph_tb_av_os[h]+= ph_tb_av[h][i]; //After all the averaging above, maybe the final pedestal is not centered on zero, so
	}                                  //this offset it meant to re-center ph_tb_av on zero, just like all other pedestal corrected pulses,
      if (Ns !=0) ph_tb_av_os[h]/=Ns;      //which are, by definition, centered on zero.
      //cout<<"ped_cm_os_"<<h<<": "<<ph_tb_av_os[h]<<endl;
    }


  //Finding Scale factor for CM correction...

  float scale[50]; //50 Bob
  for (int i=0; i<50; i++)
    {
      scale[i]=i*0.1;
    }


  float ph_cut_for_calc = ph_cut_min;
  //if (debug) ph_cut_for_calc = -300.0; //incl. all pads, except exclude extreme outliers for debug mode

  for (int h = 0; h< 4; h++)
    {
      for (int c = 0; c< 128; c++)
	{
	  if (kHit_PHA[h*128 + c] > ph_cut_for_calc)  //Include all pads for purpose of diagnostics, for debug mode only
	    {
	      residual_sum_min = 100000000000000;

	      for (int s = 0; s<50; s++)
		{
		  residual_sum[s]=0.0;
		  for (int i = 0; i < p->iValue(0,"NSAMPLES"); i++)
		    {
		      if ((kPedestal_MEAN[ h*128 + c ] - p->iValue ( c, i, 2*h)) < ph_cut_min) //only baseline samples are used to decide scaling factor
			{
			  //do the fitting with both curves centered on zero...
			  residual_sum[s] += abs( (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (scale[s]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
			}

		    }

		  if (residual_sum_min > residual_sum[s])
		    {
		      residual_sum_min = residual_sum[s];
		      residual_sum_min_scale[h][c] = scale[s];
		    }
		}
	    }
	}
    }
  //CM correcting pulse shape...
  int Nps=0;
  NNhits_cm=0;
  vector<float> pedvect;
  pedvect.clear();

  for ( int h = 0; h< 4; h++)
    {
      for ( int c = 0; c< 128; c++)
	{
	  peds=0.0;

	  if (kHit_PHA[h*128 + c] == 0.0)  //Non-Hits (we make all non-hits contribute to ped. distr., though this may not be the best way to define the ped.-------------
	    {
	      Nps=0;
	      for (int i = 0; i < p->iValue(0,"NSAMPLES"); i++)
		{
		  Nps++;
		  peds += ( (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );

		  if (h==0) h_ped_apv0_cm->Fill( (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		  if (h==1) h_ped_apv1_cm->Fill( (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		  if (h==2) h_ped_apv2_cm->Fill( (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		  if (h==3) h_ped_apv3_cm->Fill( (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );

		  if (h==0) h_pulse_ped_cm_apv0->Fill(i, (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		  if (h==1) h_pulse_ped_cm_apv1->Fill(i, (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		  if (h==2) h_pulse_ped_cm_apv2->Fill(i, (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		  if (h==3) h_pulse_ped_cm_apv3->Fill(i, (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		}
	      if (Nps !=0) kPedestal_MEAN[ h*128 + c ] = peds/Nps;//kPedestal_MEAN[ h*128 + c ] = peds/Nps;  //New CM corrected pedestal (centered on zero)
	    }
	  else  //Hits -------------           
	    {
	      NNhits_cm++;

	      for (int i = 0; i < p->iValue(0,"NSAMPLES"); i++)
		{
		  h_pulse_cm[NNhits_cm]->Fill(i, ((kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h]))) );

		  h_pulse_cmref[NNhits_cm]->Fill(i, (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		  hit_apvch_cm[NNhits_cm]=h*128 + c;

		  if (ph_cm[h*128 + c]< ((kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h]))) )
		    {
		      ph_cm[h*128 + c]= ((kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])));
		      pht_cm[h*128 + c]= i;
		    }

		  pedvect.push_back( (kPedestal_MEAN[ h*128 + c ] - p->iValue(c, i, 2*h)) - (residual_sum_min_scale[h][c]* (ph_tb_av[h][i]-ph_tb_av_os[h])) );
		}


	      //sort(pedvect.begin(), pedvect.end());
	      for (int i=0; i<kPedestal_NumberOfSamplesForFastPedestal; i++)  //kNumberOfSamplesForFastPedestal = 9, ie the 9 lowest amp. samples of a hit pad determine ped level.
		{
		  peds+=pedvect[i];

		  //cout<<"pedvect_"<<i<<": "<<pedvect[i]<<endl;

		  if (h==0) h_ped_apv0_cm->Fill(pedvect[i]);
		  if (h==1) h_ped_apv1_cm->Fill(pedvect[i]);
		  if (h==2) h_ped_apv2_cm->Fill(pedvect[i]);
		  if (h==3) h_ped_apv3_cm->Fill(pedvect[i]);
		}
	      if (kPedestal_NumberOfSamplesForFastPedestal !=0) kPedestal_MEAN[ h*128 + c ] = peds/kPedestal_NumberOfSamplesForFastPedestal;//kPedestal_MEAN[ h*128 + c ] = peds/kNumberOfSamplesForFastPedestal;  //New CM corrected pedestal (centered on zero)
	      //***how about ped mean Vs apv ch and ped rms vs apv ch
	    }
	}
    }
  //Pulse Heigth amp. based on new CM corrected pulse...
  int Nt2=0;
  for ( int j=1; j< NNhits_cm+1; j++)  //order of pads hit NNhits starts at 1
    {
      //cout<<"hit_apvch: "<<hit_apvch_cm[j]<<endl;

      if (kHit_PHA[hit_apvch_cm[j]]>ph_cut_min)  ///---Hits
	{
	  Nt2=0;
	  //cout<<"Nt2: "<<Nt2<<endl;
	  
	  //h_pulse_test_3->Reset();
	  //sprintf(title, "Pulse-Test-3(CM Correction), APV: %d, x: %d, y: %d", hit_apvch_cm[j], xmap[hit_apvch_cm[j]], ymap[hit_apvch_cm[j]]);
	  //h_pulse_test_3->SetTitle(title);
	  
	  for (int i = 0; i < p->iValue(0,"NSAMPLES"); i++)
	    {
	      float v3 = h_pulse_cm[j]->GetBinContent(i) - kPedestal_MEAN[hit_apvch_cm[j]]; //amp. corrected for new ped.
	      
	      if (i==pht_cm[hit_apvch_cm[j]]-1) ph_phtm1= v3;
	      if (i==pht_cm[hit_apvch_cm[j]]) ph_pht= v3;
	      
	      //h_pulse_test_3->Fill(i, v3);
	      
	      if ((ph_pht - ph_phtm1)/ph_pht > 0.4)
		{
		  if (i >= pht_cm[hit_apvch_cm[j]] && i <= pht_cm[hit_apvch_cm[j]]+2) //pulseheight amp. averageed over three (or less) peak samples
		    {
		      Nt2++;
		      pha_cm[hit_apvch_cm[j]]+=v3;
		      phat_cm[hit_apvch_cm[j]]+=(v3*i);
		    }
		}
	      else
		{
		  if (i >= pht_cm[hit_apvch_cm[j]]-1 && i <= pht_cm[hit_apvch_cm[j]]+1) //pulseheight amp. averageed over three (or less) peak samples
		    {
		      Nt2++;
		      pha_cm[hit_apvch_cm[j]]+=v3;
		      phat_cm[hit_apvch_cm[j]]+=(v3*i);
		    }
		}
	    }
	  if (pha_cm[hit_apvch_cm[j]] != 0) phat_cm[hit_apvch_cm[j]]=phat_cm[hit_apvch_cm[j]]/pha_cm[hit_apvch_cm[j]];  //peaktime weighted by amp. of sample
	  if (Nt2 !=0) pha_cm[hit_apvch_cm[j]]=pha_cm[hit_apvch_cm[j]]/Nt2; //pulseheight amp. averageed over three peak samples
	  
	  if (pha_cm[hit_apvch_cm[j]]>ph_cut_min)
	    {
	      
	      if (hit_apvch_cm[j]<128) h_hit_peaktime_apv0->Fill(phat_cm[hit_apvch_cm[j]]);
	      if (hit_apvch_cm[j]>=128 && hit_apvch_cm[j]<256) h_hit_peaktime_apv1->Fill(phat_cm[hit_apvch_cm[j]]);
	      if (hit_apvch_cm[j]>=256 && hit_apvch_cm[j]<384) h_hit_peaktime_apv2->Fill(phat_cm[hit_apvch_cm[j]]);
	      if (hit_apvch_cm[j]>=384) h_hit_peaktime_apv3->Fill(phat_cm[hit_apvch_cm[j]]);
	      
	    }
	}
	
      else //Non-Hits---------
	{
	  //Here we effectively require all thre samples to be above ph_cut_min
	  //Stricter cut on pulse height (reduce the possibility of identifying fluctuations as hits)
	  pha_cm[hit_apvch_cm[j]]=0.0;
	  phat_cm[hit_apvch_cm[j]]=0.0;
	}
	  
    }

}
*/
//=======================================================================
int process_event (Event *e){
  if ( e->getEvtType() == 12) {
    cout<<"N_EvtProc: "<<nEvents<<", Event: "<<e->getEvtSequence()<<endl;
  }
  //-------Reset all N counters at beginning of run--------
  if ( e->getEvtType() == 9) { // This if statement must come before cutting on event type below---  9=begin event, 1=all middle events, 12=last event
    nEvents = 0;
  }
  //************************************************************************************
  if ( e->getEvtType() != 1) return 0;  //9=begin event, 1=all middle events, 12=last event //MOVE DOWN: No other event types will be accepted below this line!!!???
  //************************************************************************************
  if ( e->getEvtSequence() <0) return 0;  //skip first n events: e->getEvtSequence() < n  //Bob
  //************************************************************************************
  kControl_EventNumber = e->getEvtSequence();
  //************************************************************************************
  nEvents++;
  if(nEvents%1000 == 0) {
    cout << "-----Events_processed: " << nEvents << ", Event_Number: " << kControl_EventNumber <<endl;
  }
  //here we fill histograms and do things
  Packet *p = e->getPacket(1010);  //this is the package for SRS
  if (p) {
    // fill the "trace" histograms for the 4 hybrids
    //    cout << "ENTRE AL PAQUETE" << endl;

    //1. checks overall health of the event and creates activity map which is used everywhere
    int ret = check_signals(p);
    if(ret!=0) {
      delete p;
      return 0; // skip event
    }
    //2. pedestal calculation
    pedestal_calculation(p);

    //3. common mode calculation
    if(0) {
      // commom_mode_calculation(p); // CALLS pulse_height(p) already
    }

    //3. compute hits
    hit_calculation(p);

    //4. clusterizer
    cluster_calculation(p);

    //===================before going out of the packet we need to delete it============================
    
    delete p;
  }
  return 0;
}
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
int main(int nn, char** arg){
  //TString filename = "/data/gem/xrayscan/COMPASS_4GEM_trial1/xrayscan_0000000205-0000.evt"; 
  //TString filename = "/data/gem/xrayscan/COMPASS_4GEM_trial2/COMAPSS_4GEM__0000000003-0000.evt"; 
  // TString filename = "/data/gem/xrayscan/COMPASS_4GEM_trial3/COMAPSS_4GEM__0000000108-0000.evt";//use this for slides
  TString file_in=arg[1]; 
  //  TString filename =Form("/data/gem/xrayscan/COMPASS_4GEM_trial3/%s",file_in.Data()); 
  TString filename =Form("/data/gem/xrayscan/COMPASS_4GEM_trial1/COMPASS/DetLayer3/%s",file_in.Data()); 

  std::cout << filename.Data() << std::endl;
  pfileopen( filename.Data() );

  //  CreateTree();
  CreateHistograms();
  //Fill_xymap();
  //pstart(); // DO NOT DO THIS
  //prun(500);
  prun();

  //calculation resolution and writing a file with real_motor_position,mean,mean_error,sigma,sigma_error 
  //resolution();
  
  //writing histograms 
  TFile *outhistos = new TFile( Form("out_resul%s.txt",file_in.Data()), "RECREATE" );

  kControl_QA_RAWMEANpc->Write();
  kControl_QA_RAWRMSpc->Write();
  kControl_QA_RAWWAVEhits->Write();
  kControl_QA_RAWWAVEnoise->Write();
  kControl_QA_RAWWAVEspy->Write();
  
  if(kPedestal_QA) {
    kPedestal_QA_MEANpc->Write();
  }
  if(kHit_QA) {
    kHit_QA_PI->Write();
    kHit_QA_PH->Write();
    kHit_QA_PHT->Write();
    kHit_QA_PHPHT->Write();
    kHit_QA_PHA->Write();
    kHit_QA_PHAT->Write();
    kHit_QA_CHANNELMAP->Write();
  }
  kCluster_QX->Write();
  kCluster_QY->Write();
  kCluster_XCENTROID->Write();
  kCluster_YCENTROID->Write();
  kCluster_XSTD_DEV->Write();
  kCluster_YSTD_DEV->Write();
  kCluster_XYCENTROID->Write();

  // file in usualmente es: compass_det3_0000000035-0000.evt
  // adding to the beginning of the file
  ofstream file_resul;
  file_resul.open(Form("out_resul%s.txt",file_in.Data()));
  //file_resul << file_in.Data() << " ";
  //file_resul << kResidual_MeanX<<" "<<kResidual_MeanX_er<<" "<<kResidual_SigmaX<<" "<<kResidual_SigmaX_er<<" \n";
  file_resul.close();
  
  outhistos->Close();

  
 return 0;

}
