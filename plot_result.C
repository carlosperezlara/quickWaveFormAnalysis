#include <TFile.h>
#include <iostream>
#include "TCanvas.h"
void plot_result(){
  //histograms to plot
  TH2F *h_centroid_motor = new TH2F("centroid_motor", "Centroid vs Motor Pos.", 1000, 0, 100, 1000, 0, 100);
  //  TH2F *h_centroid_motor = new TH2F("centroid_motor", "Centroid vs Motor Pos.", 2000, 0+0.025, 100+0.025, 2000, 0, 100);
  h_centroid_motor->GetXaxis()->SetTitle("Motor Position (mm)");
  h_centroid_motor->GetYaxis()->SetTitle("Centroid (mm)"); 

  TH1F *h_centroid = new TH1F("centroid","Centroid(mm)", 1000, 0, 100);
  TH1F *h_res = new TH1F("h_res","resolution(mm)", 100, 0, 1);

  TH2F *hresid_motor = new TH2F("hresid_motor", "Residual vs Motor Pos.", 1000, 0, 100, 1000, 0, 100);
  hresid_motor->GetXaxis()->SetTitle("Motor Position (mm)");
  hresid_motor->GetYaxis()->SetTitle("Residual (mm)");

  TH1F *h_residual = new TH1F("h_residual","residual(mm)", 1000, 0, 100);
  h_residual->GetXaxis()->SetTitle("Residual (mm)");

  TH2F *hsigma_pos = new TH2F("hsigma_pos", "Resolution vs Motor Pos.", 1000, 0, 100, 100, 0, 1.);
  hsigma_pos->GetXaxis()->SetTitle("Motor Position (mm)");
  hsigma_pos->GetYaxis()->SetTitle("Resolution (mm)");
  
  
  float x_mean;
  float x_mean_er;
  float x_res;
  float x_res_er;
  ifstream infile;
  ifstream infile2;
  float motor_x_pos;
  //  infile.open("resul_trial3.txt");
  //  infile2.open("position_trial2.txt");
  infile2.open("motor_position_trial1.txt");
  infile.open("resul_trial1.txt");
  while(infile && infile2){
  infile>>x_mean>>x_mean_er>>x_res>>x_res_er;
  infile2>>motor_x_pos;
  //std:cout<<"value "<<x_mean<<"\n";
  h_centroid->Fill(x_mean);
  h_res->Fill(x_res);
  h_centroid_motor->Fill(motor_x_pos,x_mean);
  hresid_motor->Fill(motor_x_pos,(x_mean-motor_x_pos));
  h_residual->Fill(x_mean-motor_x_pos);
  hsigma_pos->Fill(motor_x_pos,x_res);
  //cout<<"motor_position "<<motor_x_pos<<" residual "<<(x_mean-motor_x_pos)<<endl;
  //motor_x_pos=motor_x_pos+0.1;
  }

  TCanvas *n1=new TCanvas("centroid","centroid",900,700);
  n1->Divide(2,3);
  n1->cd(1);
  h_centroid->Draw();
  n1->cd(2);
  h_res->Draw();
  n1->cd(3);
  h_centroid_motor->Draw("colz");
  n1->cd(4);
  hresid_motor->Draw("colz");
  n1->cd(5);
  h_residual->Draw();
  n1->cd(6);
  hsigma_pos->Draw("colz");
  }
