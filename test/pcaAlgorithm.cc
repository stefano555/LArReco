/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/pcaAlgorithm.cc
 * 
 *  @brief  Implementation of the calohit_properties algorithm class.
 * 
 *  $Log: $
 */
#include <iostream>
#include <array>
#include <vector>
#include <iomanip>
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
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "pcaAlgorithm.h"


using namespace pandora;
using namespace lar_content;

namespace lar_reco
{

pcaAlgorithm::pcaAlgorithm() :
    m_treeName(""),
    m_fileName("")
{
    // Called when algorithm is created
}

pcaAlgorithm::~pcaAlgorithm()
{
    // Called when algorithm is destroyed
 PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "RECREATE"));
}

StatusCode pcaAlgorithm::Run()
{

    	const CaloHitList *pCaloHitList(nullptr);
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
	

	float wire_pitch = 0.3;
	float howbigbinx_mm = 1.5;
	float howbigbinx = (float (howbigbinx_mm))/10;
	int howbigbinz_mm = 3;
	float howbigbinz = (float (howbigbinz_mm))/10;
	int number_x_bin = 21;
	int number_y_bin = 21;
	int centre=0;
	int total_bins = number_x_bin*number_y_bin;
	
        std::vector<float> drifttime_vector; // drift-time position
	std::vector<float> zplane_vector;    // z-plane position
	std::vector<float> ex_vector;        // half-size of the wave width
	std::vector<float> ey_vector;        // half wire pitch
	std::vector<float> charge_vector;    // deposited integrated charge 
	std::vector<int> id_vector;          // pdg number
	std::vector<float> bin_ratio;        // ration occupied vs total bins in the small histogram
	std::vector<int> is_track_vector;    // 0 for shower 1 for track
	std::vector<float> av_distance_vector;
	std::vector<float> av_energy_eV_bin_vector;
	
 
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
		const MCParticleWeightMap map = pCaloHit->GetMCParticleWeightMap();
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
	float* drifttime = new float [drifttime_vector.size()]();
	float* zplane = new float [drifttime_vector.size()]();
	float* ex = new float [drifttime_vector.size()]();
	float* ey = new float [drifttime_vector.size()]();
	float* charge = new float [drifttime_vector.size()]();
	int* id = new int [drifttime_vector.size()]();
	

	std::copy(drifttime_vector.begin(), drifttime_vector.end(), drifttime);
	std::copy(zplane_vector.begin(), zplane_vector.end(), zplane);
	std::copy(ex_vector.begin(), ex_vector.end(), ex);
	std::copy(ey_vector.begin(), ey_vector.end(), ey);
	std::copy(charge_vector.begin(), charge_vector.end(), charge);
	std::copy(id_vector.begin(), id_vector.end(), id);
	float* drifttime_bis = new float [drifttime_vector.size()]();
	float* zplane_bis = new float [drifttime_vector.size()]();

	int offset_x = (number_x_bin -1)/2;
	int offset_y = (number_y_bin -1)/2;

	float min_x = float(offset_x)*howbigbinx;
	float min_y = float(offset_y)*howbigbinz;
	float max_x = float(offset_x+1)*howbigbinx;
	float max_y = float(offset_y+1)*howbigbinz;

	//std::cout << "Which event do you want to see? Please insert a number between 1 and " << drifttime_vector.size() << std::endl;
	//std::cin >> centre;
	//int track_counter = 0;
	//int shower_counter = 0;
	for(int l=0; l<drifttime_vector.size(); l++)
	{
		centre=l;	
		int id_centre=id[centre];
	
	        ////TRACK OR SHOWER IDENTIFICATION////////////////////////////////////////////////

		int is_track=0; //if 0 is shower, if 1 is track
		if(id_centre==11 || id_centre==22)
		{
			//shower_counter++;
			is_track = 0;
		}
		else
		{
			//track_counter++;
			is_track = 1;
		}
		is_track_vector.push_back(is_track);
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "is_track_vector", &is_track_vector));//domande: questo e' un vettore, se fai solo is_track non va...perche'? cosa cambia tra caricare
//vettore nel for cicle e fuori?
	/////////////////////////////////////////////////////////////////////
		TH2* hgrid = new TH2F("hgrid", "position vs time with integrated charge", number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);



		for(int i=0; i<drifttime_vector.size(); i++)
		{
			drifttime_bis[i] = drifttime[i] - drifttime[centre-1] + (howbigbinx/2);
			zplane_bis[i] = zplane[i] - zplane[centre-1] + (howbigbinz/2);
		///GAUSSIAN DISTRIBUTION
			int binx = hgrid->GetXaxis()->FindBin(drifttime_bis[i]);
			int biny = hgrid->GetYaxis()->FindBin(zplane_bis[i]);
			float binYcentre = hgrid->GetYaxis()->GetBinCenter(biny);
			int binx_left = hgrid->GetXaxis()->FindBin(drifttime_bis[i]-ex[i]);
			int binx_right = hgrid->GetXaxis()->FindBin(drifttime_bis[i]+ex[i]); 
			int bindiff_left = binx-binx_left;
			int bindiff_right = binx_right-binx; 

 			for (int a = 0; a<=(bindiff_left+bindiff_right);a++)
			{  
			
				float bin_low_edge = float (hgrid->GetXaxis()->GetBinLowEdge(binx-bindiff_left+a));
				float bin_up_edge = float (hgrid->GetXaxis()->GetBinUpEdge(binx-bindiff_left+a));
				float distance_low = (drifttime_bis[i])-bin_low_edge;
				float distance_up = (drifttime_bis[i])-bin_up_edge;
				float integral_low = -0.5*(erff(distance_low/((sqrt(2))*ex[i])));
				float integral_up = -0.5*(erff(distance_up/((sqrt(2))*ex[i])));
				float integral = integral_up-integral_low;
				
				hgrid->Fill(hgrid->GetXaxis()->GetBinCenter(binx-bindiff_left+a), binYcentre,charge[i]*integral);
			}
		
				
		}

	//PCA ANALYSIS//////////////////////////////////////////////////////////////////////////


		std::vector<float> pca_x_vector;
		std::vector<float> pca_z_vector;
		float x_sum=0;
		float z_sum=0;
		float half_dimension_x=(float(number_x_bin)*(howbigbinx))/2;
		float half_dimension_z=(float(number_y_bin)*(howbigbinz))/2;

		//I only select CaloHit inside my window
		for(int m=0; m<drifttime_vector.size(); m++)
		{
			if(drifttime_bis[m]>=-half_dimension_x && drifttime_bis[m]<=half_dimension_x && zplane_bis[m]>=-half_dimension_z && zplane_bis[m]<=half_dimension_z)
			{
				pca_x_vector.push_back(drifttime_bis[m]);
				pca_z_vector.push_back(zplane_bis[m]);
			}
			
		}
		
		//I calculate the average of the two dimensions
		x_sum=std::accumulate(pca_x_vector.begin(), pca_x_vector.end(), 0.0);
		z_sum=std::accumulate(pca_z_vector.begin(), pca_z_vector.end(), 0.0);
		
		float x_average= x_sum/float(pca_x_vector.size()); 
		float z_average= z_sum/float(pca_z_vector.size()); 

		float* pca_x = new float [pca_x_vector.size()]();
		float* pca_z = new float [pca_z_vector.size()]();
		std::copy(pca_x_vector.begin(), pca_x_vector.end(), pca_x);
		std::copy(pca_z_vector.begin(), pca_z_vector.end(), pca_z);
		float* b_x = new float [pca_x_vector.size()]();
		float* b_z = new float [pca_z_vector.size()]();
		float element_1=0;
		float element_2=0;
		float element_3=0;

		for(int n=0; n<drifttime_vector.size(); n++)
		{
			b_x[n]=pca_x_vector[n]-x_average;
			b_z[n]=pca_z_vector[n]-z_average;
			element_1=element_1+(b_x[n]*b_x[n]);
			element_3=element_3+(b_z[n]*b_z[n]);
			element_2=element_2+(b_x[n]*b_z[n]);	
		}

		element_1=(element_1)/(1/(drifttime_vector.size()-1));
		element_2=(element_2)/(1/(drifttime_vector.size()-1));
		element_3=(element_3)/(1/(drifttime_vector.size()-1));
		
		float element_b=element_1+element_3;
		float delta=(element_b*element_b)-4*((element_1*element_3)-(element_2*element_2));
		float major_axis=(element_b+sqrt(delta))/2;
		float minor_axis=(element_b-sqrt(delta))/2;
		std::cout << "major axis is " << major_axis << "minor axis is " << minor_axis << std::endl;

		//END OF PCA ANALYSIS/////////////////////////////////////////////////////////////////

	//HOW MUCH OF THE HISTOGRAM IS OCCUPIED////////////////////////////////////////////////////
		double bin_content = 0;
		int filled_bins = 0;
		float bins_occupied = 0;
		float distance = 0;
		float distanceX=0;
		float distanceY=0;
		float tot_distance=0;
		float av_distance=0;
		float energy_eV_bin=0;
		float tot_energy=0;
		float av_energy_eV_bin=0;
		
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
					energy_eV_bin=bin_content*170*23.6;
					tot_energy=tot_energy+energy_eV_bin;
				}	
			}
		
		}
		///RATIO
		bins_occupied = float(filled_bins)/float(total_bins);
		bin_ratio.push_back(bins_occupied);
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bin_ratio", &bin_ratio));
		////DISTANCE		
		av_distance=tot_distance/float(filled_bins);
		av_distance_vector.push_back(av_distance);
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "av_distance", &av_distance_vector));
		////ENERGY IN eV
		av_energy_eV_bin=tot_energy/float(filled_bins);
		av_energy_eV_bin_vector.push_back(av_energy_eV_bin);
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "av_energy_eV_bin", &av_energy_eV_bin_vector));
		//std::cout << "av_energy_eV_bin is " << av_energy_eV_bin << " and av_distance is " << av_distance << std::endl;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		delete hgrid;

		//std::cout << "the filled/total bins ration is " << bins_occupied << std::endl;
	////////////////////////////////////////////////////////////////////////////////////

		//PandoraMonitoringApi::Pause(this->GetPandora());
		/*TGraphErrors *gr2 = new TGraphErrors(drifttime_vector.size(),drifttime_bis,zplane_bis,ex,ey);
		TCanvas *c1 = new TCanvas();
		hgrid->SetStats(0);
		hgrid->SetXTitle("Time [cm]");
		hgrid->SetYTitle("W Plane Position [cm]");
		hgrid->GetXaxis()->SetNdivisions(-number_x_bin);
   		hgrid->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid->Draw("colz");

		TPad *grid = new TPad("grid","",0,0,1,1);
		grid->Draw();
		grid->cd();
		grid->SetGrid();
		grid->SetFillStyle(4000);
		grid->SetFrameFillStyle(0);

		TH2 *hgrid2 = new TH2C("hgrid", "position vs time with integrated charge", number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		hgrid2->Draw();
		hgrid2->SetStats(0);
		hgrid2->GetXaxis()->SetNdivisions(-number_x_bin);
		hgrid2->GetYaxis()->SetNdivisions(-number_y_bin);
		hgrid2->GetYaxis()->SetLabelOffset(999.);
		hgrid2->GetXaxis()->SetLabelOffset(999.);

		gr2->Draw("P,same");*/
	
	//PandoraMonitoringApi::Pause(this->GetPandora());

	//delete c1;//it is a root bug: it creates automatically a c1 canvas but then you have to close it
	} //END OF THE MEGA FOR LOOP
	PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
	//std::cout << "number of track-like calohits was " << track_counter << " number of shower-like events was " << shower_counter << std::endl;

    return STATUS_CODE_SUCCESS;

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode pcaAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace workshop_content
