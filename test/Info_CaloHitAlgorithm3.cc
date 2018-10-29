/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/Info_CaloHitAlgorithm3.cc
 * 
 *  @brief  Implementation of the calohit_properties algorithm class.
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
#include "Info_CaloHitAlgorithm3.h"
#define PI 3.14159265


using namespace pandora;
using namespace lar_content;

namespace lar_reco
{



StatusCode Info_CaloHitAlgorithm3::Run()
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
	int centre=0;
	//int total_bins = number_x_bin*number_y_bin;
	
        std::vector<float> drifttime_vector;
	std::vector<float> zplane_vector;
	std::vector<float> ex_vector;
	std::vector<float> ey_vector;
	std::vector<float> charge_vector;
	std::vector<int> id_vector;
	std::vector<float> energy_bin;
	std::vector<float> distance_from_axis_vector;
	std::vector<float> intermediate_sd_vector;
	std::vector<float> energy_bin_edge;
	
 
	for(const CaloHit *const pCaloHit : *pCaloHitList)
	{
		if(pCaloHit->GetHitType()!=TPC_VIEW_W)
		{
			continue;
		}
                drifttime_vector.push_back(pCaloHit->GetPositionVector().GetX());
	 	zplane_vector.push_back(pCaloHit->GetPositionVector().GetZ());
		ex_vector.push_back((pCaloHit->GetCellSize1())/2);
		ey_vector.push_back((wire_pitch)/2);
		charge_vector.push_back(pCaloHit->GetInputEnergy());
			//finding which MCParticle has given the biggest contribution
		const MCParticleWeightMap &map = pCaloHit->GetMCParticleWeightMap();
		float biggest = 0;
		int important_code = 0;
		for (const auto iter : map)
		{
    			const MCParticle *pMCParticle = iter.first;
    			int code = pMCParticle->GetParticleId(); // Code pdg (Particle Data Group)
    			float weight = iter.second;
			if(weight > biggest)
			{
				biggest = weight;
				important_code = code;
			}
		}
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

	//float min_x = float(offset_x)*howbigbinx;
	//float min_y = float(offset_y)*howbigbinz;
	//float max_x = float(offset_x+1)*howbigbinx;
	//float max_y = float(offset_y+1)*howbigbinz;

	//std::cout << "Which event do you want to see? Please insert a number between 1 and " << drifttime_vector.size() << std::endl;
	//std::cin >> centre;

	for(int l=0; l<drifttime_vector.size(); l++)
	{
		centre=l;	
		int id_centre=id_vector.at(centre);

		std::vector<float> drifttime_bis_vector;
		std::vector<float> zplane_bis_vector;
	
	        ////TRACK OR SHOWER IDENTIFICATION////////////////////////////////////////////////

		//int is_track=0; //if 0 is shower, if 1 is track
	
		std::cout << "Code PDG for the selected CaloHit is: " << id_centre << std::endl;

		if(std::fabs(id_centre)==11 || id_centre==22)
		{
			std::cout << " therefore your event is shower-like" << std::endl;
			//is_track = 0;
		}
		else
		{
			std::cout << " therefore your event is track-like" << std::endl;
			//is_track = 1;
		}
	/////////////////////////////////////////////////////////////////////
		TH2* hgrid = new TH2F("hgrid", "position vs time with integrated charge", number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid->GetXaxis()->SetDecimals();//messo per le inesattezze
		hgrid->GetYaxis()->SetDecimals();
		hgrid->GetYaxis()->SetTitleOffset(1.2);
		//hgrid->GetXaxis()->SetNdivisions(-11);//messo dopo per tentativi


		for(int i=0; i<drifttime_vector.size(); i++)
		{
			drifttime_bis_vector.push_back(drifttime_vector.at(i)-drifttime_vector.at(centre));//+(howbigbinx/2));
			zplane_bis_vector.push_back(zplane_vector.at(i)-zplane_vector.at(centre));//+(howbigbinx/2));
		///GAUSSIAN DISTRIBUTION
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
	///create new histogram filtered
		TH2* hgrid3 = new TH2F("hgrid3", "position vs time with integrated charge", number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid3->GetXaxis()->SetDecimals();//messo per le inesattezze
		hgrid3->GetYaxis()->SetDecimals();
		hgrid3->GetYaxis()->SetTitleOffset(1.2);

/////////
					//HOW MUCH OF THE HISTOGRAM IS OCCUPIED////////////////////////////////////////////////////
		float bin_content = 0;
		float bin_content_edge=0;
		float distanceX=0;
		float distanceY=0;
		//int filled_bins=0;
		int total_bins=0;
		float tot_energy=0;
		float tot_energy2=0;
		std::vector<float> binx_vector;
		std::vector<float> biny_vector;
		for(int j=1; j<=number_x_bin; j++)
		{
			for(int k=1;k<=number_y_bin; k++)
			{
				bin_content=hgrid->GetBinContent(j,k);
				total_bins++;
				////edge detection filter
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
				//end of that
				if(bin_content!=0)
				{
					//filled_bins++;
					distanceX=hgrid->GetXaxis()->GetBinCenter(j);
					distanceY=hgrid->GetYaxis()->GetBinCenter(k);
					energy_bin.push_back(bin_content);
					binx_vector.push_back(distanceX);
					biny_vector.push_back(distanceY);
					tot_energy2=tot_energy2+bin_content;						
				}	
			}
		
		}
tot_energy = accumulate(energy_bin.begin(),energy_bin.end(), 0.0);

		
	//PCA ANALYSIS//////////////////////////////////////////////////////////////////////////
		float half_dimension_x=(float(number_x_bin)*(howbigbinx))/2;

/* 
////OLD VERSION
		std::vector<float> pca_x_vector;
		std::vector<float> pca_z_vector;
		std::vector<float> pca_energy_vector;
		float x_sum=0;
		float z_sum=0;
		float energy_sum=0;
		float half_dimension_x=(float(number_x_bin)*(howbigbinx))/2;
		float half_dimension_z=(float(number_y_bin)*(howbigbinz))/2;

		//I only select CaloHit inside my window
		for(int m=0; m<drifttime_vector.size(); m++)
		{
			if(drifttime_bis_vector.at(m)>=(-half_dimension_x+(howbigbinx/2)) && drifttime_bis_vector.at(m)<=(half_dimension_x+(howbigbinx/2)) && zplane_bis_vector.at(m)>=(-half_dimension_z+(howbigbinz/2)) && zplane_bis_vector.at(m)<=(half_dimension_z+(howbigbinz/2)))
			{

				pca_x_vector.push_back(drifttime_bis_vector.at(m));
				pca_z_vector.push_back(zplane_bis_vector.at(m));
				pca_energy_vector.push_back(charge_vector.at(m));
				
			}
		
		}

		///RMS energy

                float average = accumulate( pca_energy_vector.begin(),pca_energy_vector.end(), 0.0/ pca_energy_vector.size());
		float energy_rms = sqrt( ( std::inner_product( pca_energy_vector.begin(), pca_energy_vector.end(), pca_energy_vector.begin(), 0.0 ) )/ pca_energy_vector.size());
		
		std::cout <<"the energy RMS is " << energy_rms << " whilst the energy average is " << average << std::endl;

		/////END RMS ENERGY
		
		//I calculate the average of the two dimensions
		x_sum=std::accumulate(pca_x_vector.begin(), pca_x_vector.end(), 0.0);
		z_sum=std::accumulate(pca_z_vector.begin(), pca_z_vector.end(), 0.0);
		energy_sum=std::accumulate(pca_energy_vector.begin(), pca_energy_vector.end(), 0.0);
		
		float x_average= x_sum/float(pca_x_vector.size()); //the x and z coordinates of the centroid
		float z_average= z_sum/float(pca_z_vector.size()); 
		float x_average_weighted=0;
		float y_average_weighted=0;

		float* pca_x = new float [pca_x_vector.size()]();
		float* pca_z = new float [pca_z_vector.size()]();
		float* pca_energy = new float [pca_energy_vector.size()]();
		std::copy(pca_x_vector.begin(), pca_x_vector.end(), pca_x);
		std::copy(pca_z_vector.begin(), pca_z_vector.end(), pca_z);
		std::copy(pca_energy_vector.begin(), pca_energy_vector.end(), pca_energy);
		float* b_x = new float [pca_x_vector.size()]();
		float* b_z = new float [pca_z_vector.size()]();
		float* b_x_weighted = new float [pca_x_vector.size()]();
		float* b_z_weighted = new float [pca_z_vector.size()]();
		float element_1=0;
		float element_2=0;
		float element_3=0;
		float element_1_weighted=0;
		float element_2_weighted=0;
		float element_3_weighted=0;
		//I do the unweighted verson plus I prepare for the weighted
		for(int n=0; n<pca_x_vector.size(); n++)
		{
			b_x[n]=pca_x_vector[n]-x_average;
			b_z[n]=pca_z_vector[n]-z_average;
			element_1=element_1+(b_x[n]*b_x[n]);
			element_3=element_3+(b_z[n]*b_z[n]);
			element_2=element_2+(b_x[n]*b_z[n]);
			x_average_weighted=x_average_weighted+((pca_x[n]*pca_energy[n])/energy_sum);
			y_average_weighted=y_average_weighted+((pca_z[n]*pca_energy[n])/energy_sum);	
		}

		//Weighted

		for(int p=0; p<pca_x_vector.size(); p++)
		{
			b_x_weighted[p]=pca_energy[p]*(pca_x_vector[p]-x_average_weighted);
			b_z_weighted[p]=pca_energy[p]*(pca_z_vector[p]-y_average_weighted);
			element_1_weighted=element_1_weighted+(b_x_weighted[p]*b_x_weighted[p]);
			element_3_weighted=element_3_weighted+(b_z_weighted[p]*b_z_weighted[p]);
			element_2_weighted=element_2_weighted+(b_x_weighted[p]*b_z_weighted[p]);	
		}
		float vector_pca_size=float(pca_x_vector.size());
		float denominator=1/(vector_pca_size-1);
		float denominator_weighted=1/(energy_sum-1);
		element_1=element_1*denominator;
		element_2=element_2*denominator;
		element_3=element_3*denominator;
		element_1_weighted=element_1_weighted*denominator_weighted;
		element_2_weighted=element_2_weighted*denominator_weighted;
		element_3_weighted=element_3_weighted*denominator_weighted;

		
		float element_b=element_1+element_3;
		float delta=(element_b*element_b)-4*((element_1*element_3)-(element_2*element_2));
		float major_axis=(element_b+sqrt(delta))/2;
		float minor_axis=(element_b-sqrt(delta))/2;

		float element_b_weighted=element_1_weighted+element_3_weighted;
		float delta_weighted=(element_b_weighted*element_b_weighted)-4*((element_1_weighted*element_3_weighted)-(element_2_weighted*element_2_weighted));
		float major_axis_weighted=(element_b_weighted+sqrt(delta_weighted))/2;
		float minor_axis_weighted=(element_b_weighted-sqrt(delta_weighted))/2;*/

		float energy_sum=0;

           //I calculate the average of the two dimensions
		energy_sum=std::accumulate(energy_bin.begin(), energy_bin.end(), 0.0);
		float x_average_weighted=0;
		float z_average_weighted=0;
		std::vector<float> b_x_weighted_vector;
		std::vector<float> b_z_weighted_vector;
		float element_1_weighted=0;
		float element_2_weighted=0;
		float element_3_weighted=0;
		//I do the unweighted verson plus I prepare for the weighted
		for(unsigned int element = 0; element<binx_vector.size(); element++)
		{
			x_average_weighted=x_average_weighted+((binx_vector.at(element)*energy_bin.at(element))/energy_sum);
			z_average_weighted=z_average_weighted+((biny_vector.at(element)*energy_bin.at(element))/energy_sum);	
		}
		//Weighted
		for(unsigned int element = 0; element<binx_vector.size(); element++)
		{
			b_x_weighted_vector.push_back((binx_vector.at(element)-x_average_weighted));
			b_z_weighted_vector.push_back((biny_vector.at(element)-z_average_weighted));			
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

		
		float element_b_weighted=element_1_weighted+element_3_weighted;
		float delta_weighted=(element_b_weighted*element_b_weighted)-4*((element_1_weighted*element_3_weighted)-(element_2_weighted*element_2_weighted));
		float major_axis_weighted=(element_b_weighted+sqrt(delta_weighted))/2;
		//float minor_axis_weighted=(element_b_weighted-sqrt(delta_weighted))/2;

		//eigenvector

		//float A = element_1-major_axis;
		//float E = element_1-minor_axis;

		//float eigenvector_major = -A/element_2;
		//float eigenvector_minor = -E/element_2;

		float A_weighted = element_1_weighted-major_axis_weighted;
		//float E_weighted = element_1_weighted-minor_axis_weighted;
		/*
		std::cout << "MATRIX" << std::endl;
		std::cout << element_1 << " " << element_2 << std::endl;
		std::cout << element_2 << " " << element_3 << std::endl;
		std::cout << "eigenvalue 1 " << major_axis << " and eigenvalue 2 " << minor_axis << std:: endl;
		std::cout << "eigenspace" << std::endl;
		std::cout << 1 << " " << 1 << std::endl;
		std::cout << eigenvector_major << " " << eigenvector_minor << std::endl;
		std::cout << "check orthogonality " << 1 + (eigenvector_major)*(eigenvector_minor) << std::endl;

		std::cout << "MATRIX_weighted" << std::endl;
		std::cout << element_1_weighted << " " << element_2_weighted << std::endl;
		std::cout << element_2_weighted << " " << element_3_weighted << std::endl;
		std::cout << "eigenvalue 1_weighted " << major_axis_weighted << " and eigenvalue 2_weighted " << minor_axis_weighted << std:: endl;
		std::cout << "eigenspace_weighted" << std::endl;
		std::cout << 1 << " " << 1 << std::endl;
		std::cout << eigenvector_major_weighted << " " << eigenvector_minor_weighted << std::endl;
		std::cout << "check orthogonality " << 1 + (eigenvector_major_weighted)*(eigenvector_minor_weighted) << std::endl;*/

		//DISTANCE CALOHIT FROM WEIGHTED MAJOR AXIS
		/* OLD VERSION
		float denominator_distance = sqrt((eigenvector_major_weighted*eigenvector_major_weighted)+1);
		float numerator_distance=0;		
		float distance_from_axis=0;
		std::vector<float> distance_from_axis_vector;
		for (int kl=0; kl<pca_x_vector.size(); kl++)
		{
			numerator_distance = std::abs(pca_z[kl]-(eigenvector_major_weighted*pca_x[kl]));			
			distance_from_axis= numerator_distance/denominator_distance;
			distance_from_axis_vector.push_back(distance_from_axis);

			std::cout << "distance " << distance_from_axis << std::endl;//" and x_distance " << x_distance << std::endl;
		} 
		
		float* distance_from_axis_array = new float [distance_from_axis_vector.size()]();
		std::copy(distance_from_axis_vector.begin(), distance_from_axis_vector.end(), distance_from_axis_array);

		float average_distance_from_axis = accumulate( distance_from_axis_vector.begin(),distance_from_axis_vector.end(), 0.0)/ distance_from_axis_vector.size();
		std::vector<float> intermediate_sd_vector;
		float* intermediate_sd = new float [distance_from_axis_vector.size()]();
		for (int yu=0; yu<distance_from_axis_vector.size(); yu++)
		{
			intermediate_sd[yu]=distance_from_axis_array[yu]-average_distance_from_axis;
			intermediate_sd_vector.push_back(intermediate_sd[yu]);
		}
		float sd_distance_from_axis = sqrt( ( std::inner_product( intermediate_sd_vector.begin(), intermediate_sd_vector.end(), intermediate_sd_vector.begin(), 0.0 ) )/ (intermediate_sd_vector.size()-1));

		std::cout << "the average distance from axis is " << average_distance_from_axis << " cm and the rms is " << sd_distance_from_axis << std::endl;*/

				//DISTANCE CALOHIT FROM WEIGHTED MAJOR AXIS

		float eigenvector_major_weighted = 0;
		if(element_2_weighted!=0)
		{
		
			eigenvector_major_weighted = -A_weighted/element_2_weighted;
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
			float average_distance_from_axis = accumulate( distance_from_axis_vector.begin(),distance_from_axis_vector.end(), 0.0)/ distance_from_axis_vector.size();
		
			for (unsigned int element = 0; element < distance_from_axis_vector.size(); element++)
			{
				intermediate_sd_vector.push_back(distance_from_axis_vector.at(element)-average_distance_from_axis);
			}

			float sd_distance_from_axis = sqrt( ( std::inner_product( intermediate_sd_vector.begin(), intermediate_sd_vector.end(), intermediate_sd_vector.begin(), 0.0 ) )/ (intermediate_sd_vector.size()-1));
			std::cout << "the average distance from axis is " << average_distance_from_axis << " cm and the rms is " << sd_distance_from_axis << std::endl;

		}
		//END OF PCA ANALYSIS/////////////////////////////////////////////////////////////////
		////ANDYS' SUGGESTION


	//float angle_degree=atan (1/eigenvector_major_weighted) * 180 / PI;
			float deviation=0;
			//float normalization_1=1/(float(number_y_bin)*3);
			//float normalization_2=1/float(number_x_bin*(number_x_bin-3));
			
	float angle=atan(1/eigenvector_major_weighted);
	std::cout << "the angle between pca major axis and y axis is " << angle << std::endl;
		TH2* hgrid4 = new TH2F("hgrid4", "position vs time with integrated charge", number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid4->GetXaxis()->SetDecimals();//messo per le inesattezze
		hgrid4->GetYaxis()->SetDecimals();
		hgrid4->GetYaxis()->SetTitleOffset(1.2);
		float unrotated_charge =0;
		float rot_charge=0;
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
float contribute=0;
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
std::cout <<"deviation is " << deviation << std::endl;	
std::cout <<"total rotated charge is " << rot_charge << " whilst unratated is " << unrotated_charge << std::endl;
	//ANDY ANALYIS/////////////////////////////////////////////
	int crossed_bins=0;
	float crossed_energy=0;
	for(unsigned int element=0; element<binx_vector.size(); element++)
	{
		int crossing_counter=0;
		float binx_left=binx_vector.at(element)-(howbigbinx/2);
		float binz_down=biny_vector.at(element)-(howbigbinz/2);
		float binx_right=binx_vector.at(element)+(howbigbinx/2);
		float binz_up=biny_vector.at(element)+(howbigbinz/2);
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
		if(crossing_counter == 2)
		{
			crossed_bins++;
			crossed_energy=crossed_energy+energy_bin.at(element);
		}
	}
float crossed_ratio=crossed_energy/tot_energy;

std::cout << "crossed bins " << crossed_bins << " , total energy " << tot_energy << " tot energy 2 " << tot_energy2 << " ,crossed energy " << crossed_energy << " , and ratio " << crossed_ratio << std::endl;

///////////END OF ANDY ANALYSIS/////////////////////////


		PandoraMonitoringApi::Pause(this->GetPandora());
		//TGraphErrors *gr2 = new TGraphErrors(drifttime_vector.size(),drifttime_bis,zplane_bis,ex,ey);
		////normal one
		TCanvas *c1 = new TCanvas("c1","canvas", 500,500);
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

	
		TH2 *hgrid2 = new TH2C("hgrid2", "position vs time with integrated charge", number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid2->Draw();
		hgrid2->SetStats(0);
		hgrid2->GetXaxis()->SetNdivisions(-number_x_bin);
		hgrid2->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid2->GetYaxis()->SetLabelOffset(999.);
		hgrid2->GetXaxis()->SetLabelOffset(999.);

		//gr2->Draw("P,same");
	


		/*
		TArrow *major = new TArrow (double (x_average),double (z_average),double (x_average+0.3),double (z_average+(eigenvector_major*0.3)), 0.05, "|>");
		TArrow *minor = new TArrow (double (x_average),double (z_average),double (x_average+0.3),double (z_average+(eigenvector_minor*0.3)), 0.05, "|>");
		TArrow *major_weighted = new TArrow (double (x_average_weighted),double (y_average_weighted),double (x_average_weighted+0.3),double (y_average_weighted+(eigenvector_major_weighted*0.3)), 0.05, "|>");
		TArrow *minor_weighted = new TArrow (double (x_average_weighted),double (y_average_weighted),double (x_average_weighted+0.3),double (y_average_weighted+(eigenvector_minor_weighted*0.3)), 0.05, "|>");
		//major->SetAngle(60);		
		major->SetLineWidth(2);
                major->SetFillColor(2);
		//minor->SetAngle(60);
		minor->SetLineWidth(2);
                minor->SetFillColor(2);
		major->Draw("same");
		minor->Draw("same");
		major_weighted->Draw("same");
		minor_weighted->Draw("same");*/


		TMarker *centroid_weighted = new TMarker (double (x_average_weighted),double(z_average_weighted),20);//full circle
		centroid_weighted->SetMarkerSize(2.0);
		centroid_weighted->SetMarkerColor(1);//black
		centroid_weighted->Draw("same");
		/*TF1 *fa1 = new TF1("fa1","[0]*(x-[1])+[2]",-half_dimension_x+(howbigbinx/2),half_dimension_x+(howbigbinx/2));
		fa1->SetParameter(0,eigenvector_major);
		fa1->SetParameter(1,x_average);
		fa1->SetParameter(2,z_average);
		fa1->SetLineStyle(2);//dash
		fa1->SetLineColor(6); //pourple
		fa1->SetLineWidth(5);   		
		fa1->Draw("same");*/
		//fa2 is passing through zero
		TF1 *fa2 = new TF1("fa1","[0]*(x)",-half_dimension_x+(howbigbinx/2),half_dimension_x+(howbigbinx/2));
		fa2->SetParameter(0,eigenvector_major_weighted);
		fa2->SetLineColor(2); //red
		fa2->SetLineWidth(5);
   		//fa2->Draw("same");
		//fa3 is passing through centroid
		TF1 *fa3 = new TF1("fa3","[0]*(x-[1])+[2]",-half_dimension_x+(howbigbinx/2),half_dimension_x+(howbigbinx/2));
		fa3->SetParameter(0,eigenvector_major_weighted);
		fa3->SetParameter(1,x_average_weighted);
		fa3->SetParameter(2,z_average_weighted);
		fa3->SetLineColor(1); //black
		fa3->SetLineWidth(5);  		
		fa3->Draw("same");
/////////rotated
		TCanvas *c3 = new TCanvas("c3","canvas", 500,500);
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
	
	
	PandoraMonitoringApi::Pause(this->GetPandora());

	delete c1;//it is a root bug: it creates automatically a c1 canvas but then you have to close it
	//delete c2;
	delete c3;
	binx_vector.clear();
	biny_vector.clear();
	energy_bin.clear();
	}
    return STATUS_CODE_SUCCESS;

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode Info_CaloHitAlgorithm3::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_reco
