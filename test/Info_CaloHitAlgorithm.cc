/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/Info_CaloHitAlgorithm.cc
 * 
 *  @brief  Implementation of the info_calohit algorithm class.
 * 
 *  $Log: $
 */
#include <iostream>
#include <array>
#include <vector>
#include <iomanip>
#include "math.h"
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
#include "Info_CaloHitAlgorithm.h"

////se ci sono problemi sono alla riga 98 e 99 sti cazzo di int float non si capisce niente
using namespace pandora;
using namespace lar_content;

namespace lar_reco
{
StatusCode Info_CaloHitAlgorithm::Run()
{
	
    	const CaloHitList *pCaloHitList(nullptr);
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

	float wire_pitch = 0.3;
	float howbigbinx_mm = 1.5;
	float howbigbinx = howbigbinx_mm/10;
	int howbigbinx_submm = howbigbinx_mm*10;
	int howbigbinz_mm = 3;
	float howbigbinz = (float (howbigbinz_mm))/10;

        std::vector<float> drifttime_vector;
	std::vector<float> zplane_vector;
	std::vector<float> ex_vector;
	std::vector<float> ey_vector;
	std::vector<int> id_particle_vector;
	
	
 
	for(const CaloHit *const pCaloHit : *pCaloHitList)
	{

		if(pCaloHit->GetHitType()!=TPC_VIEW_W)
		{
			continue;
		}
		
		
		//std::cout << "InputHit - HitType: " << pCaloHit->GetHitType() << ", the position in time (cm) is: " << pCaloHit->GetPositionVector().GetX() << std::endl;
                drifttime_vector.push_back(pCaloHit->GetPositionVector().GetX());
	 	zplane_vector.push_back(pCaloHit->GetPositionVector().GetZ());
		ex_vector.push_back((pCaloHit->GetCellSize1())/2);
		ey_vector.push_back((wire_pitch)/2);
		
	}
	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//FINDING THE PERFECT RANGE AND BINING

	//I look for the minimum and maximum values in both drift time (cm) and z plane (cm) and the maximum width (cm) as well
	float min_time = *(std::min_element(drifttime_vector.begin(), drifttime_vector.end()));
	float max_time = *(std::max_element(drifttime_vector.begin(), drifttime_vector.end()));
	float min_pos = *(std::min_element(zplane_vector.begin(), zplane_vector.end()));	
	float max_pos = *(std::max_element(zplane_vector.begin(), zplane_vector.end()));
	float max_width = *(std::max_element(ex_vector.begin(), ex_vector.end()));
	
	//I define an initial range. Since the width at each side of the central point is width/2, I stay at width-distance from the border
	float min_time_border = min_time - max_width - howbigbinx;
	float max_time_border = max_time + max_width + howbigbinx;
	float init_mintime_rounded_float = nearbyintf(10*(min_time_border)); //rounded and in mm!
	float init_maxtime_rounded_float = nearbyintf(10*(max_time_border)); //rounded and in mm!
	int init_mintime_rounded = int (init_mintime_rounded_float);
	int init_maxtime_rounded = int (init_maxtime_rounded_float);


	float min_pos_border = min_pos - (5*wire_pitch);
	float max_pos_border = max_pos + (5*wire_pitch);
	float init_minpos_rounded_float = nearbyintf(10*(min_pos_border));
	float init_maxpos_rounded_float = nearbyintf(10*(max_pos_border));
	int init_minpos_rounded = int (init_minpos_rounded_float);
	int init_maxpos_rounded = int (init_maxpos_rounded_float);


	int initial_range = init_maxtime_rounded - init_mintime_rounded; //rounded and in mm!
	int initial_range_pos = init_maxpos_rounded - init_minpos_rounded;

	//in order to have a perfect bining, I have to round to mm

	float modulus_submm = float((initial_range)*10 % howbigbinx_submm);
	float toadd = howbigbinx_mm-modulus_submm/10;//rounded and in mm!
	float final_mintime_border_float = float (init_mintime_rounded) - (toadd/2);
	float final_maxtime_border_float = float (init_maxtime_rounded) + (toadd/2);
	float final_mintime_border = final_mintime_border_float/10; //rounded and in cm!
	float final_maxtime_border = final_maxtime_border_float/10; //rounded and in cm!
	float bin_time = (final_maxtime_border - final_mintime_border) / howbigbinx;

	float toaddpos = float (howbigbinz_mm-(initial_range_pos % howbigbinz_mm));
	float final_minpos_border_float = float (init_minpos_rounded) - (toaddpos/2);
	float final_maxpos_border_float = float (init_maxpos_rounded) + (toaddpos/2);
	float final_minpos_border = final_minpos_border_float/10; //rounded and in cm!
	float final_maxpos_border = final_maxpos_border_float/10; //rounded and in cm!
	float bin_pos = (final_maxpos_border - final_minpos_border) / howbigbinz;
	//////////////////////////////////////////////////////////////////////////////////////////////////

	float* drifttime = new float [drifttime_vector.size()]();
	float* zplane = new float [drifttime_vector.size()]();
	float* ex = new float [drifttime_vector.size()]();
	float* ey = new float [drifttime_vector.size()]();
	
	

	std::copy(drifttime_vector.begin(), drifttime_vector.end(), drifttime);
	std::copy(zplane_vector.begin(), zplane_vector.end(), zplane);
	std::copy(ex_vector.begin(), ex_vector.end(), ex);
	std::copy(ey_vector.begin(), ey_vector.end(), ey);
	

	TCanvas *c1 = new TCanvas();//new
	TH2* h1 = new TH2F("h1", "position vs time with integrated charge", bin_time, final_mintime_border, final_maxtime_border, bin_pos,final_minpos_border, final_maxpos_border);

	for(const CaloHit *const pCaloHit2 : *pCaloHitList)
	{


/////////parte nuova!!!
		if(pCaloHit2->GetHitType()!=TPC_VIEW_W)
		{
			continue;
		}
//////////////////

		int binx = h1->GetXaxis()->FindBin(pCaloHit2->GetPositionVector().GetX());//-> function related to a pointer . related to object
		int biny = h1->GetYaxis()->FindBin(pCaloHit2->GetPositionVector().GetZ());      
		float binYcentre = h1->GetYaxis()->GetBinCenter(biny);
		int binx_left = h1->GetXaxis()->FindBin(pCaloHit2->GetPositionVector().GetX()-(pCaloHit2->GetCellSize1())/2);
		int binx_right = h1->GetXaxis()->FindBin(pCaloHit2->GetPositionVector().GetX()+(pCaloHit2->GetCellSize1())/2);

		int bindiff_left = binx-binx_left;
		int bindiff_right = binx_right-binx;

		for (int a = 0; a<=(bindiff_left+bindiff_right);a++)
		{  
			
			float bin_low_edge = float (h1->GetXaxis()->GetBinLowEdge(binx-bindiff_left+a));
			float bin_up_edge = float (h1->GetXaxis()->GetBinUpEdge(binx-bindiff_left+a));
			float distance_low = (pCaloHit2->GetPositionVector().GetX())-bin_low_edge;
			float distance_up = (pCaloHit2->GetPositionVector().GetX())-bin_up_edge;
			float integral_low = -0.5*(erff(distance_low/(sqrt(2)*((pCaloHit2->GetCellSize1())/2))));
			float integral_up = -0.5*(erff(distance_up/(sqrt(2)*((pCaloHit2->GetCellSize1())/2))));
			float integral = integral_up-integral_low;
			//std::cout << "element " << pCaloHit2 << " charge " << pCaloHit2->GetInputEnergy() << " bin " << a << " smeared charge " << (pCaloHit2->GetInputEnergy())*(integral) << " integral " << integral << std::endl; 
			//when you ask for pCaloHit2 you get something like 0x2c60060 because it is a pointer and you get the address of the memory
			
			h1->Fill(h1->GetXaxis()->GetBinCenter(binx-bindiff_left+a), binYcentre,(pCaloHit2->GetInputEnergy())*(integral));
		}
	}

	
	TGraphErrors *gr = new TGraphErrors(drifttime_vector.size(),drifttime,zplane,ex,ey);
	h1->SetStats(0);
	h1->SetXTitle("Time [cm]");
	h1->SetYTitle("W Plane Position [cm]");
	
	h1->Draw("colz");
	
	gr->Draw("P,same");


	delete[]drifttime;
	delete[]zplane;
	delete[]ex;
	delete[]ey;
	PandoraMonitoringApi::Pause(this->GetPandora());
	//delete h1;//new
	delete c1;//new
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode Info_CaloHitAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here

    return STATUS_CODE_SUCCESS;
}

} // namespace workshop_content

