/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/tree_makerAlgorithm.cc
 * 
 *  @brief  Implementation of the tree_maker algorithm class.
 * 
 *  $Log: $
 */

#include <iostream>
#include <array>
#include <vector>
#include <iomanip>
#include <cmath>
#include <numeric>
#include "math.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TArrow.h"
#include "TMarker.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "tree_makerAlgorithm.h"
#include "TROOT.h"
#define PI 3.14159265


using namespace pandora;
using namespace lar_content;

namespace lar_reco
{



StatusCode tree_makerAlgorithm::Run()
{

    	const CaloHitList *pCaloHitList(nullptr);
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
	

	float wire_pitch = 0.48;
	float howbigbinx_mm = 2.5;//previously 1.5
	float howbigbinx = (float (howbigbinx_mm))/10;
	float howbigbinz_mm = 4.8;
	float howbigbinz = (float (howbigbinz_mm))/10;
	int number_x_bin = 21;
	int number_y_bin = 11;
	int total_bins = number_x_bin*number_y_bin;
	int centre=0;
	
        std::vector<float> drifttime_vector; // drift-time position
	//CaloHitList CaloHitVector;
	std::vector<float> zplane_vector;    // z-plane position
	std::vector<float> ex_vector;        // half-size of the wave width
	std::vector<float> ey_vector;        // half wire pitch
	std::vector<float> charge_vector;    // deposited integrated charge 
	std::vector<int> id_vector;          // pdg number
	std::vector<float> distance_from_axis_vector;
	std::vector<float> energy_bin;
	std::vector<float> percentage_vector;
        std::vector<float> binx_vector;
        std::vector<float> biny_vector;
	std::vector<float> track_weight_vector;
	std::vector<float> full_weight_vector;
	std::vector<float> energy_bin_edge;

	//int idx(0);
 
	for(const CaloHit *const pCaloHit : *pCaloHitList)
	{
		if(pCaloHit->GetHitType()!=TPC_VIEW_W)
		{
			continue;
		}
               // idx++;
               // CaloHitVector.push_back(pCaloHit);
                drifttime_vector.push_back(pCaloHit->GetPositionVector().GetX());
	 	zplane_vector.push_back(pCaloHit->GetPositionVector().GetZ());
		ex_vector.push_back((pCaloHit->GetCellSize1())/2);
		ey_vector.push_back((wire_pitch)/2);
		charge_vector.push_back(pCaloHit->GetInputEnergy());
			//finding which MCParticle has given the biggest contribution
		const MCParticleWeightMap &map = pCaloHit->GetMCParticleWeightMap();
		float biggest = 0;
		int important_code = 0;
		float weight_sum = 0;
		float weight_track = 0;
//this->weight(const MCParticleWeightMap &map)
		for (const auto iter : map)
		{
    			const MCParticle *pMCParticle = iter.first;
    			int code = pMCParticle->GetParticleId(); // Code pdg (Particle Data Group)
    			float weight = iter.second;
			weight_sum += weight;
			if(std::fabs(code)!= 11 && code!= 22)
			{
				weight_track += weight;
			}			
			if(weight > biggest)
			{
				biggest = weight;
				important_code = code;
			}
		}
		float percentage=-1;
		
		if(!map.empty())
		{
			percentage=weight_track/weight_sum;
			full_weight_vector.push_back(weight_sum);
			track_weight_vector.push_back(weight_track);
				
		}
		else
		{
			full_weight_vector.push_back(-9999);
			track_weight_vector.push_back(-9999);
		}
		percentage_vector.push_back(percentage);
		id_vector.push_back(important_code);

	}
	float* ey = new float [drifttime_vector.size()]();
	
	std::copy(ey_vector.begin(), ey_vector.end(), ey);
	

	int offset_x = (number_x_bin -1)/2;
	int offset_y = (number_y_bin -1)/2;
	
	float min_x = float(offset_x)*howbigbinx+(howbigbinx/2);
	float min_y = float(offset_y)*howbigbinz+(howbigbinz/2);
	float max_x = float(offset_x)*howbigbinx+(howbigbinx/2);
	float max_y = float(offset_y)*howbigbinz+(howbigbinz/2);


	int answer=1;
	while(answer==1)
	{

	std::cout << "Do you want to see a histogram related to this event or do you want to skip to next event? Enter 1 to see a histogram or 0 to skip to to the next event" << std::endl;
	std::cin >> answer;
	if(answer==0)
	{
		break;
	}
	std::cout << "Which event do you want to see? Please insert a number between 1 and " << drifttime_vector.size() << std::endl;
	std::cin >> centre;

	//for(int l=0; l<drifttime_vector.size(); l++)
	//{
		//centre=l;	

		std::vector<float> drifttime_bis_vector;
		std::vector<float> zplane_bis_vector;
		
		int id_centre=id_vector.at(centre);//UNUSED
                int is_track = this->IsTrack(id_centre);//UNUSED
//Filling the small histogram with gaussian distribution
		std::string title;
		title = "position vs time with integrated charge, Calohit number " + std::to_string(centre);
		TH2* hgrid = new TH2F(this->safe_name("hgrid"), title.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid->GetXaxis()->SetDecimals();//messo per le inesattezze
		hgrid->GetYaxis()->SetDecimals();
		hgrid->GetYaxis()->SetTitleOffset(1.2);
		//hgrid->GetXaxis()->SetNdivisions(-11);//messo dopo per tentativi
		this->vector_shifter(drifttime_vector,zplane_vector,drifttime_bis_vector,zplane_bis_vector,centre);
		this->histogram_filler_gaussian(ex_vector,hgrid,charge_vector,drifttime_bis_vector,zplane_bis_vector);
//Filling other 2 small histogram with total weight and track-weight
		TH2* hgridtotal = new TH2F(this->safe_name("hgridtotal"), title.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		TH2* hgridtrack = new TH2F(this->safe_name("hgridtrack"), title.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		this->histogram_filler(hgridtotal,drifttime_bis_vector,zplane_bis_vector,full_weight_vector);
		this->histogram_filler(hgridtrack,drifttime_bis_vector,zplane_bis_vector,track_weight_vector);
///create new histogram filtered
		TH2* hgrid3 = new TH2F(this->safe_name("hgrid3"), title.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid3->GetXaxis()->SetDecimals();//messo per le inesattezze
		hgrid3->GetYaxis()->SetDecimals();
		hgrid3->GetYaxis()->SetTitleOffset(1.2);
		float bin_content_edge=0;

		for(int j=1; j<=number_x_bin; j++)
		{
			for(int k=1;k<=number_y_bin; k++)
			{

				////edge detection filter//////////////////////
				if(j>=2 && k >=2)
				{
					bin_content_edge=-(hgrid->GetBinContent(j-1,k+1))-(hgrid->GetBinContent(j,k+1))-(hgrid->GetBinContent(j+1,k+1))-(hgrid->GetBinContent(j-1,k))+60*(hgrid->GetBinContent(j,k))-(hgrid->GetBinContent(j+1,k))-(hgrid->GetBinContent(j-1,k-1))-(hgrid->GetBinContent(j,k-1))-(hgrid->GetBinContent(j+1,k-1));	
				}
				if(j==1 && k>=2)
				{
					bin_content_edge=-(hgrid->GetBinContent(j,k+1))-(hgrid->GetBinContent(j+1,k+1))+60*(hgrid->GetBinContent(j,k))-(hgrid->GetBinContent(j+1,k))-(hgrid->GetBinContent(j,k-1))-(hgrid->GetBinContent(j+1,k-1));
				}
				if(j>=2 && k==1)
				{
					bin_content_edge=-(hgrid->GetBinContent(j-1,k+1))-(hgrid->GetBinContent(j,k+1))-(hgrid->GetBinContent(j+1,k+1))-(hgrid->GetBinContent(j-1,k))+60*(hgrid->GetBinContent(j,k))-(hgrid->GetBinContent(j+1,k));	
				}
				if(j==1 && k==1)
				{
					bin_content_edge=-(hgrid->GetBinContent(j,k+1))-(hgrid->GetBinContent(j+1,k+1))+60*(hgrid->GetBinContent(j,k))-(hgrid->GetBinContent(j+1,k));	
				}
				energy_bin_edge.push_back(bin_content_edge);
				//std::cout << "energy bin edge " << bin_content_edge << std::endl;
				///end of the filter
				///filling new filter-histo
				hgrid3->Fill(hgrid->GetXaxis()->GetBinCenter(j),hgrid->GetYaxis()->GetBinCenter(k),bin_content_edge);
	
			}
		
		}
//Studying some properties of the small histogram
		float bins_occupied = 0;
		float av_distance=0;
		float av_energy_bin=0;
		float tot_energy=0;
		this->histogram_study(number_x_bin,number_y_bin,bins_occupied,av_distance,av_energy_bin,binx_vector,biny_vector,energy_bin,hgrid, total_bins,tot_energy);
		float rms_energy = this->rms_calculator(energy_bin, av_energy_bin);
		float total_weight_histo = this->histogram_reader(number_x_bin,number_y_bin,hgridtotal);
		float track_weight_histo = this->histogram_reader(number_x_bin,number_y_bin,hgridtrack);
		float ratio_weight_histo = track_weight_histo/total_weight_histo;
//PCA calculations
		int isGoodCaloHit = 1;
		float energy_sum=0;
		energy_sum=std::accumulate(energy_bin.begin(), energy_bin.end(), 0.0);
		float x_average_weighted=0;
    		float z_average_weighted=0;
		matrix covariance(0,0,0);
		this->average_weighted(binx_vector,biny_vector,energy_bin,energy_sum, x_average_weighted,z_average_weighted);
		this->covariance_matrix(binx_vector,biny_vector,energy_bin,x_average_weighted,z_average_weighted,energy_sum,covariance);
	
		if (std::fabs(covariance.GetElement2()) < std::numeric_limits<float>::epsilon()&&std::fabs(covariance.GetElement3()) > std::numeric_limits<float>::epsilon())
		{
		    isGoodCaloHit = 0;
		}
		/*/ this is a tool to search and plot stuff with specific features, to be used in the for loop
		if (isGoodCaloHit!=0)
		{
			delete hgrid;
			drifttime_bis_vector.clear();
			zplane_bis_vector.clear();
			energy_bin.clear();
			binx_vector.clear();
			biny_vector.clear();
			continue;
		}
		/*/

		float major_axis_weighted = 0;
		float minor_axis_weighted = 0;
		this->axes_calculator(covariance,major_axis_weighted,minor_axis_weighted);
		float axis_ratio_weighted=minor_axis_weighted/major_axis_weighted; //UNUSED	
		float centroid_distance_weighted = sqrt((x_average_weighted*x_average_weighted)+(z_average_weighted*z_average_weighted)); //UNUSED
		float eigenvector_major_weighted = 0;
		if (std::fabs(covariance.GetElement2()) > std::numeric_limits<float>::epsilon()&&isGoodCaloHit>0)
		{
			this->eigenvector_calculator(covariance,major_axis_weighted,eigenvector_major_weighted);
		}
		this->distance_from_axis_calculator(eigenvector_major_weighted,x_average_weighted,z_average_weighted,binx_vector,biny_vector,distance_from_axis_vector);
		float average_distance_from_axis = accumulate( distance_from_axis_vector.begin(),distance_from_axis_vector.end(), 0.0)/ distance_from_axis_vector.size();
		float sd_distance_from_axis = this->standard_deviation(distance_from_axis_vector,average_distance_from_axis);
		

		float half_dimension_x=(float(number_x_bin)*(howbigbinx))/2;
//Crossing major axis
		float crossed_energy=0;
			
		if (isGoodCaloHit)
		{
			crossed_energy=this->crossed_energy_calculator(binx_vector,biny_vector,energy_bin,eigenvector_major_weighted,x_average_weighted,z_average_weighted,howbigbinx,howbigbinz);
		}
		float ratio_crossed_energy = float (crossed_energy)/tot_energy;

//Rotation method
		std::string title2;
		title2 = "position vs time with integrated charge 2, CaloHit " + std::to_string(centre);
		TH2* hgrid4 = new TH2F(this->safe_name("hgrid4"), title2.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		float angle = PI/2;
		if (std::fabs(covariance.GetElement2()) > std::numeric_limits<float>::epsilon()&&isGoodCaloHit>0)
		{
			angle=atan(1/eigenvector_major_weighted);
		}
		if (isGoodCaloHit==0&&covariance.GetElement3()>covariance.GetElement1())
		{
			angle=0;
		}		
                float rot_charge=0;
		float unrotated_charge=0;
		float deviation=this->rotation_method(number_x_bin, number_y_bin,hgrid,hgrid4,rot_charge,angle,unrotated_charge);
		float diff_charge=this->difference_charge(number_x_bin, number_y_bin, rot_charge, energy_sum);


		std::cout <<"just to use them " << rms_energy << ratio_weight_histo << axis_ratio_weighted << centroid_distance_weighted <<sd_distance_from_axis <<ratio_crossed_energy << deviation << diff_charge << is_track <<std::endl;

		std::cout << "bins_occupied " << bins_occupied << " ,deviation " << deviation << " ,diff_charge " << diff_charge << " ,major axis " << major_axis_weighted << " ,minor axis " << minor_axis_weighted <<  " ,isGoodCaloHit? " << isGoodCaloHit << " , angle " << angle << " ,tot energy " << tot_energy <<  std::endl;
		std::cout << "centroid x " << x_average_weighted << " , z " << z_average_weighted << std::endl;

		std::cout << "MATRIX " << std::endl;
		std::cout << covariance.GetElement1() << " " << covariance.GetElement2() << std::endl;
		std::cout << covariance.GetElement2() << " " << covariance.GetElement3() << std::endl;

//////MONTECARLO PARTICLE//////////////////////////
		std::string title_montecarlo;
		title_montecarlo = "Cheated Histogram, Calohit number " + std::to_string(centre);
		TH2* hgrid5 = new TH2F(this->safe_name("hgrid5"), title_montecarlo.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid5->GetXaxis()->SetDecimals();//messo per le inesattezze
		hgrid5->GetYaxis()->SetDecimals();
		hgrid5->GetYaxis()->SetTitleOffset(1.2);
		for(int i=0; i<drifttime_vector.size(); i++)
		{
			//int binx = hgrid5->GetXaxis()->FindBin(drifttime_bis_vector.at(i));
			//int biny = hgrid5->GetYaxis()->FindBin(zplane_bis_vector.at(i));	
			if(std::fabs(id_vector.at(i))==11 || id_vector.at(i)==22)
			{
				hgrid5->Fill(drifttime_bis_vector.at(i), zplane_bis_vector.at(i),0.1);
			}
			else
			{
				hgrid5->Fill(drifttime_bis_vector.at(i), zplane_bis_vector.at(i),1);
			}		
		}
//////////////////////////////////////////////////////////////////

///////////////////////PERCENTAGE OF TRACKNESS//////////////////////////////////
		std::string title_percentage;
		title_percentage = "Percentage of Trackness Histogram, Calohit number " + std::to_string(centre);
		TH2* hgrid7 = new TH2F(this->safe_name("hgrid7"), title_percentage.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid7->GetXaxis()->SetDecimals();//messo per le inesattezze
		hgrid7->GetYaxis()->SetDecimals();
		hgrid7->GetYaxis()->SetTitleOffset(1.2);

		for(int i=0; i<drifttime_vector.size(); i++)
		{
			if(percentage_vector.at(i)==0)
			{
				percentage_vector.at(i)=0.00001;
			}		
			hgrid7->Fill(drifttime_bis_vector.at(i), zplane_bis_vector.at(i),percentage_vector.at(i));

		}

/////////PDG PARTICLE NUMBER PRINTED//////////////////////////////////
		std::string title_pdg;
		title_pdg = "PDG Histogram, Calohit number " + std::to_string(centre);
		TH2* hgrid6 = new TH2F(this->safe_name("hgrid6"), title_pdg.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid6->GetXaxis()->SetDecimals();//messo per le inesattezze
		hgrid6->GetYaxis()->SetDecimals();
		hgrid6->GetYaxis()->SetTitleOffset(1.2);
		for(int i=0; i<drifttime_vector.size(); i++)
		{
			hgrid6->Fill(drifttime_bis_vector.at(i), zplane_bis_vector.at(i),id_vector.at(i));
		}
////////////////////////////////////////////////////////////////////////		
		PandoraMonitoringApi::Pause(this->GetPandora());
		//TGraphErrors *gr2 = new TGraphErrors(drifttime_vector.size(),drifttime_bis,zplane_bis,ex,ey);
		////normal one
		TCanvas *c1 = new TCanvas("c1","canvas", 900,900);
		c1->SetRightMargin(0.15);
		c1->SetLeftMargin(0.15);
		c1->SetTopMargin(0.15);
		c1->SetBottomMargin(0.15);
		hgrid->SetStats(0);
		hgrid->SetXTitle("Time [cm]");
		hgrid->SetYTitle("W Plane Position [cm]");
		hgrid->GetXaxis()->SetNdivisions(number_x_bin);
		hgrid->GetXaxis()->SetTickLength(0);
   		hgrid->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid->GetYaxis()->SetTitleOffset(2);
		hgrid->GetXaxis()->SetTitleOffset(2);		
		hgrid->Draw("colz");

		TPad *grid = new TPad("grid","",0,0,1,1);
		grid->SetRightMargin(0.15);
		grid->SetLeftMargin(0.15);
		grid->SetTopMargin(0.15);
		grid->SetBottomMargin(0.15);
		grid->Draw();
		grid->cd();
		grid->SetGrid();
		grid->SetFillStyle(4000);
		grid->SetFrameFillStyle(0);
		/////////filter
/*
		TCanvas *c2 = new TCanvas("c2","canvas", 500,500);
		c2->SetRightMargin(0.15);
		c2->SetLeftMargin(0.15);
		c2->SetTopMargin(0.15);
		c2->SetBottomMargin(0.15);
		
		hgrid3->SetStats(0);
		hgrid3->SetXTitle("Time [cm]");
		hgrid3->SetYTitle("W Plane Position [cm]");
		hgrid3->GetXaxis()->SetNdivisions(number_x_bin);
		hgrid3->GetXaxis()->SetTickLength(0);
   		hgrid3->GetYaxis()->SetNdivisions(-number_y_bin);
		//hgrid3->Draw("colz text");

		TPad *grid = new TPad("grid","",0,0,1,1);
		grid->SetRightMargin(0.15);
		grid->SetLeftMargin(0.15);
		grid->SetTopMargin(0.15);
		grid->SetBottomMargin(0.15);
		grid->Draw();
		grid->cd();
		grid->SetGrid();
		grid->SetFillStyle(4000);
		grid->SetFrameFillStyle(0);*/


		TH2 *hgrid2 = new TH2C("hgrid2", "", number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid2->Draw();
		hgrid2->SetStats(0);
		hgrid2->GetXaxis()->SetNdivisions(-number_x_bin);
		hgrid2->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid2->GetYaxis()->SetLabelOffset(999.);
		hgrid2->GetXaxis()->SetLabelOffset(999.);
	

		TMarker *centroid_weighted = new TMarker (double (x_average_weighted),double(z_average_weighted),20);//full circle
		centroid_weighted->SetMarkerSize(2.0);
		centroid_weighted->SetMarkerColor(1);//black
		centroid_weighted->Draw("same");
		TF1 *fa2 = new TF1("fa1","[0]*(x)",-half_dimension_x+(howbigbinx/2),half_dimension_x+(howbigbinx/2));
		fa2->SetParameter(0,eigenvector_major_weighted);
		fa2->SetLineColor(2); //red
		fa2->SetLineWidth(5);
   		//fa2->Draw("same");
		//fa3 is passing through centroid
		TF1 *fa3 = new TF1("fa3","[0]*(x-[1])+[2]",-half_dimension_x+(howbigbinx/2),half_dimension_x+(howbigbinx/2));
		if (isGoodCaloHit==0&&covariance.GetElement3()>covariance.GetElement1())
		{
		TLine *line = new TLine(0,-min_y,0,min_y);
		line->SetLineColor(1); //black
		line->SetLineWidth(5);  		
		line->Draw("same");
		}
		else
		{
		fa3->SetParameter(0,eigenvector_major_weighted);
		fa3->SetParameter(1,x_average_weighted);
		fa3->SetParameter(2,z_average_weighted);
		fa3->SetLineColor(1); //black
		fa3->SetLineWidth(5);  		
		fa3->Draw("same");
		}
/////////rotated
		TCanvas *c3 = new TCanvas("c3","canvas", 900,900);
		c3->SetRightMargin(0.15);
		c3->SetLeftMargin(0.15);
		c3->SetTopMargin(0.15);
		c3->SetBottomMargin(0.15);
		
		hgrid4->SetStats(0);
		hgrid4->SetXTitle("Time [cm]");
		hgrid4->SetYTitle("W Plane Position [cm]");
		hgrid4->GetXaxis()->SetNdivisions(number_x_bin);
		hgrid4->GetXaxis()->SetTickLength(0);
   		hgrid4->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid4->GetYaxis()->SetTitleOffset(2);
		hgrid4->GetXaxis()->SetTitleOffset(2);		
		hgrid4->Draw("colz");

		TPad *grid3 = new TPad("grid3","",0,0,1,1);
		grid3->SetRightMargin(0.15);
		grid3->SetLeftMargin(0.15);
		grid3->SetTopMargin(0.15);
		grid3->SetBottomMargin(0.15);
		grid3->Draw();
		grid3->cd();
		grid3->SetGrid();
		grid3->SetFillStyle(4000);
		grid3->SetFrameFillStyle(0);
		hgrid2->Draw();
/////////////////////////////////////
		TCanvas *c4 = new TCanvas("c4","canvas", 900,900);
		c4->SetRightMargin(0.15);
		c4->SetLeftMargin(0.15);
		c4->SetTopMargin(0.15);
		c4->SetBottomMargin(0.15);
		
		hgrid5->SetStats(0);
		hgrid5->SetXTitle("Time [cm]");
		hgrid5->SetYTitle("W Plane Position [cm]");
		hgrid5->GetXaxis()->SetNdivisions(number_x_bin);
		hgrid5->GetXaxis()->SetTickLength(0);
   		hgrid5->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid5->GetYaxis()->SetTitleOffset(2);
		hgrid5->GetXaxis()->SetTitleOffset(2);		
		hgrid5->Draw("colz");

		TPad *grid4 = new TPad("grid4","",0,0,1,1);
		grid4->SetRightMargin(0.15);
		grid4->SetLeftMargin(0.15);
		grid4->SetTopMargin(0.15);
		grid4->SetBottomMargin(0.15);
		grid4->Draw();
		grid4->cd();
		grid4->SetGrid();
		grid4->SetFillStyle(4000);
		grid4->SetFrameFillStyle(0);
		hgrid2->Draw();
///////////////////////////////////////////////////////////////////
		TCanvas *c5 = new TCanvas("c5","canvas", 900,900);
		c5->SetRightMargin(0.15);
		c5->SetLeftMargin(0.15);
		c5->SetTopMargin(0.15);
		c5->SetBottomMargin(0.15);

		hgrid6->SetStats(0);
		hgrid6->SetXTitle("Time [cm]");
		hgrid6->SetYTitle("W Plane Position [cm]");
		hgrid6->GetXaxis()->SetNdivisions(number_x_bin);
		hgrid6->GetXaxis()->SetTickLength(0);
   		hgrid6->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid6->GetYaxis()->SetTitleOffset(2);
		hgrid6->GetXaxis()->SetTitleOffset(2);
		hgrid6->Draw("text");
		
		TPad *grid5 = new TPad("grid5","",0,0,1,1);
		grid5->SetRightMargin(0.15);
		grid5->SetLeftMargin(0.15);
		grid5->SetTopMargin(0.15);
		grid5->SetBottomMargin(0.15);
		grid5->Draw();
		grid5->cd();
		grid5->SetGrid();
		grid5->SetFillStyle(4000);
		grid5->SetFrameFillStyle(0);
		hgrid2->Draw();
		
//////////////////////////////////////////////////////////
		TCanvas *c6 = new TCanvas("c6","canvas", 900,900);
		c6->SetRightMargin(0.15);
		c6->SetLeftMargin(0.15);
		c6->SetTopMargin(0.15);
		c6->SetBottomMargin(0.15);

		hgrid7->SetStats(0);
		hgrid7->SetXTitle("Time [cm]");
		hgrid7->SetYTitle("W Plane Position [cm]");
		hgrid7->GetXaxis()->SetNdivisions(number_x_bin);
		hgrid7->GetXaxis()->SetTickLength(0);
   		hgrid7->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid7->GetYaxis()->SetTitleOffset(2);
		hgrid7->GetXaxis()->SetTitleOffset(2);
		hgrid7->Draw("colz");
		hgrid7->SetMaximum(1);
		hgrid7->SetMinimum(0);
		
		TPad *grid6 = new TPad("grid6","",0,0,1,1);
		grid6->SetRightMargin(0.15);
		grid6->SetLeftMargin(0.15);
		grid6->SetTopMargin(0.15);
		grid6->SetBottomMargin(0.15);
		grid6->Draw();
		grid6->cd();
		grid6->SetGrid();
		grid6->SetFillStyle(4000);
		grid6->SetFrameFillStyle(0);
		hgrid2->Draw();

/////////////////////////////////////////////////////////

		TCanvas *c7 = new TCanvas("c7","canvas", 900,900);
		c7->SetRightMargin(0.15);
		c7->SetLeftMargin(0.15);
		c7->SetTopMargin(0.15);
		c7->SetBottomMargin(0.15);
		
		hgrid->Draw("text");
		TPad *grid7 = new TPad("grid7","",0,0,1,1);
		grid7->SetRightMargin(0.15);
		grid7->SetLeftMargin(0.15);
		grid7->SetTopMargin(0.15);
		grid7->SetBottomMargin(0.15);
		grid7->Draw();
		grid7->cd();
		grid7->SetGrid();
		grid7->SetFillStyle(4000);
		grid7->SetFrameFillStyle(0);
		hgrid2->Draw();


///////////////////////////////////////////////////////////////
	
	PandoraMonitoringApi::Pause(this->GetPandora());
	delete hgrid;//normal
	delete hgrid2;//weird
	delete hgrid3;//filtered
	delete hgrid4;//rotated
	delete hgrid5;//cheated
	delete hgrid6;//pdg number
	delete hgrid7;//percentage
	delete hgridtotal;
	delete hgridtrack;
	delete c1;
	//delete c2;
	delete c3;
	delete c4;
	delete c5;
	delete c6;
	delete c7;
	binx_vector.clear();
	biny_vector.clear();
	energy_bin.clear();
	drifttime_bis_vector.clear();
	zplane_bis_vector.clear();
	//} //MEGALOOP FOR o questo o while
}//while loop
	delete [] ey;
	percentage_vector.clear();
    return STATUS_CODE_SUCCESS;

}
//------------------------------------------------------------------------------------------------------------------------------------------

int tree_makerAlgorithm::IsTrack(int id_centre)
{
    if(std::fabs(id_centre)==11 || id_centre==22)
    { 
        return 0;
    }
    else
    {
        return 1;
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::vector_shifter(const FloatVector &drifttime_vector,const FloatVector &zplane_vector,FloatVector &drifttime_bis_vector,FloatVector &zplane_bis_vector,int centre)
{
    for(int i=0; i<drifttime_vector.size(); i++)
    {
	drifttime_bis_vector.push_back(drifttime_vector.at(i)-drifttime_vector.at(centre));//+(howbigbinx/2));
        zplane_bis_vector.push_back(zplane_vector.at(i)-zplane_vector.at(centre));//+(howbigbinx/2));
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::histogram_filler_gaussian(const FloatVector &ex_vector,TH2* hgrid,const FloatVector &charge_vector,const FloatVector &drifttime_bis_vector,const FloatVector &zplane_bis_vector)
{
    for(int i=0; i<drifttime_bis_vector.size(); i++)
    {
        int binx = hgrid->GetXaxis()->FindBin(drifttime_bis_vector.at(i));
	int biny = hgrid->GetYaxis()->FindBin(zplane_bis_vector.at(i));
	float binYcentre = hgrid->GetYaxis()->GetBinCenter(biny);
	int binx_left = hgrid->GetXaxis()->FindBin(drifttime_bis_vector.at(i)-ex_vector.at(i));
	int binx_right = hgrid->GetXaxis()->FindBin(drifttime_bis_vector.at(i)+ex_vector.at(i));
	int bindiff_left = binx-binx_left;
	int bindiff_right = binx_right-binx; 

 	for (int a = 0; a<=(bindiff_left+bindiff_right);a++)
	{  
			
	    float bin_low_edge = float (hgrid->GetXaxis()->GetBinLowEdge(binx-bindiff_left+a));
	    float bin_up_edge = float (hgrid->GetXaxis()->GetBinUpEdge(binx-bindiff_left+a));
	    float distance_low = (drifttime_bis_vector.at(i))-bin_low_edge;
	    float distance_up = (drifttime_bis_vector.at(i))-bin_up_edge;
	    float integral_low = -0.5*(erff(distance_low/((sqrt(2))*ex_vector.at(i))));
	    float integral_up = -0.5*(erff(distance_up/((sqrt(2))*ex_vector.at(i))));
	    float integral = integral_up-integral_low;		
	    hgrid->Fill(hgrid->GetXaxis()->GetBinCenter(binx-bindiff_left+a), binYcentre,charge_vector.at(i)*integral);
	}			
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::histogram_filler(TH2* hgrid, const FloatVector &drifttime_bis_vector,const FloatVector &zplane_bis_vector,const FloatVector &track_weight_vector)
{
    for(int i=0; i<drifttime_bis_vector.size(); i++)
    {
	hgrid->Fill(drifttime_bis_vector.at(i), zplane_bis_vector.at(i),track_weight_vector.at(i));
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

TString tree_makerAlgorithm::safe_name(const TString &initial_name)
{
 if(gROOT && gROOT->FindObjectAny(initial_name.Data()))
	delete gROOT->FindObjectAny(initial_name.Data());
return TString(initial_name);
} 
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::histogram_study(int number_x_bin, int number_y_bin,float &bins_occupied,float &av_distance,float &av_energy_bin,FloatVector &binx_vector,FloatVector &biny_vector,FloatVector &energy_bin,TH2* hgrid,float total_bins,float &tot_energy)
{
    float tot_distance=0;
    double bin_content = 0;
    int filled_bins = 0;
    float distance = 0;
    float distanceX=0;
    float distanceY=0;
    for(int j=1; j<=number_x_bin; j++)
    {
	for(int k=1;k<=number_y_bin; k++)
        {
	    bin_content=hgrid->GetBinContent(j,k);
	    if(bin_content!=0)
	    {
		filled_bins++;
		distanceX=hgrid->GetXaxis()->GetBinCenter(j);
		distanceY=hgrid->GetYaxis()->GetBinCenter(k);
		distance=sqrt((distanceX*distanceX)+(distanceY*distanceY));
		tot_distance=tot_distance+distance;
		energy_bin.push_back(bin_content);
		binx_vector.push_back(distanceX);
		biny_vector.push_back(distanceY);						
	    }	
	}	
    }
    ///RATIO
    bins_occupied = float(filled_bins)/float(total_bins);
    ////DISTANCE		
    av_distance=tot_distance/float(filled_bins);
    ////ENERGY per BIN
    tot_energy = accumulate(energy_bin.begin(),energy_bin.end(), 0.0)/energy_bin.size();		
    av_energy_bin=tot_energy/float(filled_bins);
}
//------------------------------------------------------------------------------------------------------------------------------------------
float tree_makerAlgorithm::histogram_reader(int number_x_bin, int number_y_bin,TH2* hgrid)
{
    double bin_content = 0;
    for(int j=1; j<=number_x_bin; j++)
    {
	for(int k=1;k<=number_y_bin; k++)
        {
	    if(bin_content<0)
	    {
		continue;
	    }
	    bin_content=bin_content+hgrid->GetBinContent(j,k);
        }
    }
    return float(bin_content);
}
//------------------------------------------------------------------------------------------------------------------------------------------
float tree_makerAlgorithm::rms_calculator(const FloatVector &energy_bin, float av_energy_bin)
{
    std::vector<float> intermediate_energy_vector;

    for (unsigned int element = 0; element < energy_bin.size(); element++)
    {
        intermediate_energy_vector.push_back(energy_bin.at(element)-av_energy_bin);
    }		
    float rms_energy = sqrt( ( std::inner_product( intermediate_energy_vector.begin(), intermediate_energy_vector.end(), intermediate_energy_vector.begin(), 0.0 ) )/ (intermediate_energy_vector.size()-1));
    intermediate_energy_vector.clear();
    return rms_energy =0;
}
//------------------------------------------------------------------------------------------------------------------------------------------
tree_makerAlgorithm::matrix::matrix(float element1, float element2, float element3) : 
    m_element1(element1),
    m_element2(element2),
    m_element3(element3)
{
}
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::average_weighted(const FloatVector &binx_vector, const FloatVector &biny_vector, const FloatVector &energy_bin, float energy_sum, float &x_average_weighted,float &z_average_weighted)
{
    for(unsigned int element = 0; element<binx_vector.size(); element++)
    {
        x_average_weighted=x_average_weighted+((binx_vector.at(element)*energy_bin.at(element))/energy_sum);
        z_average_weighted=z_average_weighted+((biny_vector.at(element)*energy_bin.at(element))/energy_sum);	
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::covariance_matrix(const FloatVector &binx_vector, const FloatVector &biny_vector, const FloatVector &energy_bin, float x_average_weighted,float z_average_weighted,float energy_sum, matrix &covariance)
{
    std::vector<float> b_x_weighted_vector;
    std::vector<float> b_z_weighted_vector;
    float element_1_weighted=0;
    float element_2_weighted=0;
    float element_3_weighted=0;

    for(unsigned int element = 0; element<binx_vector.size(); element++)
    {	
        b_x_weighted_vector.push_back((binx_vector.at(element)-x_average_weighted));
	b_z_weighted_vector.push_back((biny_vector.at(element)-z_average_weighted));	
	//std::cout << "element 1 " << energy_bin.at(element)*b_x_weighted_vector.at(element)*b_x_weighted_vector.at(element) << std::endl; //togli
	std::cout << "element 2 " << energy_bin.at(element)*b_x_weighted_vector.at(element)*b_z_weighted_vector.at(element) << std::endl; //togli
	//std::cout << "element 3 " << energy_bin.at(element)*b_z_weighted_vector.at(element)*b_z_weighted_vector.at(element) << std::endl; //togli
	element_1_weighted=element_1_weighted+(energy_bin.at(element)*b_x_weighted_vector.at(element)*b_x_weighted_vector.at(element));
        element_3_weighted=element_3_weighted+(energy_bin.at(element)*b_z_weighted_vector.at(element)*b_z_weighted_vector.at(element));
	element_2_weighted=element_2_weighted+(energy_bin.at(element)*b_x_weighted_vector.at(element)*b_z_weighted_vector.at(element));	
                        
    }
    b_x_weighted_vector.clear();
    b_z_weighted_vector.clear();

    float denominator_weighted=1/(energy_sum-1);
    element_1_weighted=element_1_weighted*denominator_weighted;
    element_2_weighted=element_2_weighted*denominator_weighted;
    element_3_weighted=element_3_weighted*denominator_weighted;

    matrix newMatrix(element_1_weighted, element_2_weighted, element_3_weighted);

    covariance = newMatrix;

}
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::axes_calculator(const matrix &covariance, float &major_axis_weighted, float &minor_axis_weighted)
{
    float element_b_weighted=covariance.GetElement1()+covariance.GetElement3();
    float delta_weighted=(element_b_weighted*element_b_weighted)-4*((covariance.GetElement1()*covariance.GetElement3())-(covariance.GetElement2()*covariance.GetElement2()));
    major_axis_weighted=(element_b_weighted+sqrt(delta_weighted))/2;
    minor_axis_weighted=(element_b_weighted-sqrt(delta_weighted))/2;
}
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::eigenvector_calculator(const matrix &covariance, float major_axis_weighted, float &eigenvector_major_weighted)
{
    float A_weighted = covariance.GetElement1()-major_axis_weighted;
    eigenvector_major_weighted = -A_weighted/covariance.GetElement2();
}
//------------------------------------------------------------------------------------------------------------------------------------------
void tree_makerAlgorithm::distance_from_axis_calculator(float eigenvector_major_weighted,float x_average_weighted,float z_average_weighted,const FloatVector &binx_vector,const FloatVector &biny_vector, FloatVector &distance_from_axis_vector)
{
    float denominator_distance = sqrt((eigenvector_major_weighted*eigenvector_major_weighted)+1);
    float q = -(eigenvector_major_weighted*x_average_weighted)+z_average_weighted;
    float numerator_distance=0;		
    float distance_from_axis=0;
    for (unsigned int element = 0; element<binx_vector.size(); element++)
    {
	numerator_distance = std::abs(biny_vector.at(element)-(eigenvector_major_weighted*binx_vector.at(element))-q);			
	distance_from_axis= numerator_distance/denominator_distance;
        distance_from_axis_vector.push_back(distance_from_axis);
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
float tree_makerAlgorithm::standard_deviation(const FloatVector &distance_from_axis_vector,float average_distance_from_axis)
{
    std::vector<float> intermediate_sd_vector;
    for (unsigned int element = 0; element < distance_from_axis_vector.size(); element++)
        {
            intermediate_sd_vector.push_back(distance_from_axis_vector.at(element)-average_distance_from_axis);
	}
	float sd_distance_from_axis = sqrt( ( std::inner_product( intermediate_sd_vector.begin(), intermediate_sd_vector.end(), intermediate_sd_vector.begin(), 0.0 ) )/ (intermediate_sd_vector.size()-1));
	intermediate_sd_vector.clear();
	return sd_distance_from_axis;
}
//------------------------------------------------------------------------------------------------------------------------------------------
float tree_makerAlgorithm::crossed_energy_calculator(const FloatVector &binx_vector,const FloatVector &biny_vector,const FloatVector &energy_bin,float eigenvector_major_weighted,float x_average_weighted,float z_average_weighted,float howbigbinx, float howbigbinz)
{	
	float crossed_energy=0;
        for(unsigned int element=0; element<binx_vector.size(); element++)
	{
	    int crossing_counter=0;
	    float binx_left=binx_vector.at(element)-(howbigbinx/2);
            float binz_down=biny_vector.at(element)-(howbigbinz/2);
	    float binx_right=binx_vector.at(element)+(howbigbinx/2);
	    float binz_up=biny_vector.at(element)+(howbigbinz/2);
	    if(std::fabs(eigenvector_major_weighted) > std::numeric_limits<float>::epsilon())
            {
	        float inverse_eigen=1/eigenvector_major_weighted;	
	        float z_point_left=(eigenvector_major_weighted*(binx_left-x_average_weighted))+z_average_weighted;
	        float z_point_right=(eigenvector_major_weighted*(binx_right-x_average_weighted))+z_average_weighted;
	        float x_point_down=(inverse_eigen*(binz_down-z_average_weighted))+x_average_weighted;
	        float x_point_up=(inverse_eigen*(binz_up-z_average_weighted))+x_average_weighted;

	        if(z_point_left <= binz_up && z_point_left >= binz_down)
	            {
		        crossing_counter++;
		    }
	        if(z_point_right <= binz_up && z_point_right >= binz_down)
		    {
		        crossing_counter++;
		    }
	        if(x_point_down <= binx_right && x_point_down >= binx_left)
		    {
	                crossing_counter++;
		    }
	        if(x_point_up <= binx_right && x_point_up >= binx_left)
		    {
	                crossing_counter++;
		    }
	    }
	    if(std::fabs(eigenvector_major_weighted) < std::numeric_limits<float>::epsilon())
	    {
		crossing_counter=2;
	    }
	    if(crossing_counter == 2)
		{
		    crossed_energy=crossed_energy+energy_bin.at(element);
		}
        }
    return crossed_energy;
}
//------------------------------------------------------------------------------------------------------------------------------------------
float tree_makerAlgorithm::rotation_method(int number_x_bin, int number_y_bin,TH2 * hgrid,TH2* hgrid4,float &rot_charge,float angle,float &unrotated_charge)
{		
    float deviation=0;
    float contribute=0;
    for(int j=1; j<=number_x_bin; j++)
    {
        for(int k=1;k<=number_y_bin; k++)
	{
	    float old_x= hgrid->GetXaxis()->GetBinCenter(j);				
	    float old_y= hgrid->GetYaxis()->GetBinCenter(k);
	    float new_x=std::cos(angle)*old_x-std::sin(angle)*old_y;
	    float new_y=std::sin(angle)*old_x+std::cos(angle)*old_y;
	    hgrid4->Fill(new_x,new_y,hgrid->GetBinContent(j,k));
	    unrotated_charge=unrotated_charge+hgrid->GetBinContent(j,k);
	}
    }
	
    for(int j=1; j<=number_x_bin; j++)
    {
        for(int k=1;k<=number_y_bin; k++)
	{ 		
	    if((j==((number_x_bin-1)/2)) || (j==((number_x_bin-1)/2)+1) || (j==((number_x_bin-1)/2)+2))
	    {
		contribute=(hgrid4->GetBinContent(j,k));
            }
	    else
	    {
		contribute=(-hgrid4->GetBinContent(j,k));
	    }
	    deviation=deviation+contribute;
            rot_charge=rot_charge+hgrid4->GetBinContent(j,k);
	}
    }

    deviation=deviation/rot_charge;
    return deviation;
}
//------------------------------------------------------------------------------------------------------------------------------------------
float tree_makerAlgorithm::difference_charge(int number_x_bin, int number_y_bin, float rot_charge, float energy_sum)
{
    float diff_charge=0;
    for(int j=1; j<=number_x_bin; j++)
    {
        for(int k=1;k<=number_y_bin; k++)
	{
	    diff_charge=(energy_sum-rot_charge)/energy_sum;
	}
    }
    return diff_charge;
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode tree_makerAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_reco
