/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/calohit_propertiesAlgorithm.cc
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
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "calohit_propertiesAlgorithm.h"
#include "TROOT.h"

#define PI 3.14159265

using namespace pandora;
using namespace lar_content;

namespace lar_reco
{

calohit_propertiesAlgorithm::calohit_propertiesAlgorithm() :
    m_treeName(""),
    m_fileName(""),
    m_eventNumber(-1)
{
    // Called when algorithm is created
}

calohit_propertiesAlgorithm::~calohit_propertiesAlgorithm()
{
    // Called when algorithm is destroyed
 PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "RECREATE"));
}

StatusCode calohit_propertiesAlgorithm::Run()
{
        m_eventNumber++;

    	const CaloHitList *pCaloHitList(nullptr);
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
	

	float wire_pitch = 0.48;
	float howbigbinx_mm = 2.5;//previously 1.5
	float howbigbinx = (float (howbigbinx_mm))/float (10);
	float howbigbinz_mm = 4.8;
	float howbigbinz = (float (howbigbinz_mm))/float (10);
	int number_x_bin = 21;
	int number_y_bin = 11;
	int centre=0;
	int total_bins = number_x_bin*number_y_bin;

	//typedef std::vector<float> FloatVector;

        std::vector<float> drifttime_vector; // drift-time position
        CaloHitList CaloHitVector;
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

       int idx(0);
	for(const CaloHit *const pCaloHit : *pCaloHitList)
	{
		if(pCaloHit->GetHitType()!=TPC_VIEW_W)
		{
			continue;
		}
                idx++;
                CaloHitVector.push_back(pCaloHit);

// Energy of calo hit 5 in vector = caloHitVector.at(4)->GetEnergy();
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

	int offset_x = (number_x_bin -1)/2;
	int offset_y = (number_y_bin -1)/2;

	float min_x = float(offset_x)*howbigbinx+(howbigbinx/2);
	float min_y = float(offset_y)*howbigbinz+(howbigbinz/2);
	float max_x = float(offset_x)*howbigbinx+(howbigbinx/2);
	float max_y = float(offset_y)*howbigbinz+(howbigbinz/2);

	

	for(int l=0; l<drifttime_vector.size(); l++)
	{
		std::vector<float> drifttime_bis_vector;
		std::vector<float> zplane_bis_vector;
		centre=l;	
		
		int id_centre=id_vector.at(centre);

                int is_track = this->IsTrack(id_centre);

		float position_x=drifttime_vector.at(l);
		float position_y=zplane_vector.at(l);		
		float single_charge=charge_vector.at(l);
		
//Filling the small histogram with gaussian distribution
		std::string title;
		title = "position vs time with integrated charge, CaloHit number " + std::to_string(centre);
		TH2* hgrid = new TH2F(this->safe_name("hgrid"), title.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		this->vector_shifter(drifttime_vector,zplane_vector,drifttime_bis_vector,zplane_bis_vector,centre);
		this->histogram_filler_gaussian(ex_vector,hgrid,charge_vector,drifttime_bis_vector,zplane_bis_vector);
//Filling other 2 small histogram with total weight and track-weight
		TH2* hgrid2 = new TH2F(this->safe_name("hgrid2"), title.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		TH2* hgrid3 = new TH2F(this->safe_name("hgrid3"), title.c_str(), number_x_bin, -min_x, max_x, number_y_bin, -min_y, max_y);
		this->histogram_filler(hgrid2,drifttime_bis_vector,zplane_bis_vector,full_weight_vector);
		this->histogram_filler(hgrid3,drifttime_bis_vector,zplane_bis_vector,track_weight_vector);
//Studying some properties of the small histogram
		float bins_occupied = 0;
		float av_distance=0;
		float av_energy_bin=0;
		float tot_energy=0;
		this->histogram_study(number_x_bin,number_y_bin,bins_occupied,av_distance,av_energy_bin,binx_vector,biny_vector,energy_bin,hgrid, total_bins,tot_energy);
		float rms_energy = this->rms_calculator(energy_bin, av_energy_bin);
		float total_weight_histo = this->histogram_reader(number_x_bin,number_y_bin,hgrid2);
		float track_weight_histo = this->histogram_reader(number_x_bin,number_y_bin,hgrid3);
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
		
		float element1=covariance.GetElement1();//test
		float element2=covariance.GetElement2();//test
		float element3=covariance.GetElement3();//test
		float major_axis_weighted = 0;
		float minor_axis_weighted = 0;
		this->axes_calculator(covariance,major_axis_weighted,minor_axis_weighted);
		float axis_ratio_weighted=minor_axis_weighted/major_axis_weighted;	
		float centroid_distance_weighted = sqrt((x_average_weighted*x_average_weighted)+(z_average_weighted*z_average_weighted));
		float eigenvector_major_weighted = 0;
		if (std::fabs(covariance.GetElement2()) > std::numeric_limits<float>::epsilon()&&isGoodCaloHit>0)
		{
			this->eigenvector_calculator(covariance,major_axis_weighted,eigenvector_major_weighted);
		}

		this->distance_from_axis_calculator(eigenvector_major_weighted,x_average_weighted,z_average_weighted,binx_vector,biny_vector,distance_from_axis_vector);
		float average_distance_from_axis = accumulate( distance_from_axis_vector.begin(),distance_from_axis_vector.end(), 0.0)/ distance_from_axis_vector.size();
		float sd_distance_from_axis = this->standard_deviation(distance_from_axis_vector,average_distance_from_axis);
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
		float unrotated_charge=0;
                float rot_charge=0;
		float deviation=this->rotation_method(number_x_bin, number_y_bin,hgrid,hgrid4,rot_charge,angle,unrotated_charge);
		float diff_charge=this->difference_charge(number_x_bin, number_y_bin, rot_charge, energy_sum);
//Filling a Root Tree
		
		float percentage_fill=percentage_vector.at(l);

                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventNumber", m_eventNumber));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "CaloHitNumber", l));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bin_ratio", bins_occupied));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "av_distance", av_distance));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rms_energy", rms_energy));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "is_track", is_track));
		//PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "index", centre));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "av_distance", av_distance));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "av_energy_bin", av_energy_bin));			
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "position_x", position_x));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "position_y", position_y));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "single_charge", single_charge));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "axis_ratio_weighted", axis_ratio_weighted));
	        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroid_distance_weighted", centroid_distance_weighted));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "average_distance_from_axis", average_distance_from_axis));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "sd_distance_from_axis", sd_distance_from_axis));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ratio_crossed_energy", ratio_crossed_energy));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isGoodCaloHit", isGoodCaloHit));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "deviation", deviation));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "diff_charge", diff_charge));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "percentage_trackness", percentage_fill));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ratio_weight_histo", ratio_weight_histo));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "unrotated_charge", unrotated_charge));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rot_charge", rot_charge));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "angle", angle));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "element1", element1));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "element2", element2));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "element3", element3));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eigenvector_major_weighted", eigenvector_major_weighted));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "major_eigenvalue", major_axis_weighted));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minor_eigenvalue", minor_axis_weighted));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroid_x", x_average_weighted));//test
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroid_y", z_average_weighted));//test
		PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
		
//Deleting objects before exiting for loop		
		delete hgrid;
		delete hgrid4;
		drifttime_bis_vector.clear();
		zplane_bis_vector.clear();
		energy_bin.clear();
		binx_vector.clear();
		biny_vector.clear();
		
	} //END OF THE MEGA FOR LOOP
	percentage_vector.clear();

	return STATUS_CODE_SUCCESS;

}

//------------------------------------------------------------------------------------------------------------------------------------------

int calohit_propertiesAlgorithm::IsTrack(int id_centre)
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
void calohit_propertiesAlgorithm::vector_shifter(const FloatVector &drifttime_vector,const FloatVector &zplane_vector,FloatVector &drifttime_bis_vector,FloatVector &zplane_bis_vector,int centre)
{
    for(int i=0; i<drifttime_vector.size(); i++)
    {
	drifttime_bis_vector.push_back(drifttime_vector.at(i)-drifttime_vector.at(centre));//+(howbigbinx/2));
        zplane_bis_vector.push_back(zplane_vector.at(i)-zplane_vector.at(centre));//+(howbigbinx/2));
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
void calohit_propertiesAlgorithm::histogram_filler_gaussian(const FloatVector &ex_vector,TH2* hgrid,const FloatVector &charge_vector,const FloatVector &drifttime_bis_vector,const FloatVector &zplane_bis_vector)
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
void calohit_propertiesAlgorithm::histogram_filler(TH2* hgrid, const FloatVector &drifttime_bis_vector,const FloatVector &zplane_bis_vector,const FloatVector &track_weight_vector)
{
    for(int i=0; i<drifttime_bis_vector.size(); i++)
    {
	hgrid->Fill(drifttime_bis_vector.at(i), zplane_bis_vector.at(i),track_weight_vector.at(i));
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

TString calohit_propertiesAlgorithm::safe_name(const TString &initial_name)
{
 if(gROOT && gROOT->FindObjectAny(initial_name.Data()))
	delete gROOT->FindObjectAny(initial_name.Data());
return TString(initial_name);
} 
//------------------------------------------------------------------------------------------------------------------------------------------
void calohit_propertiesAlgorithm::histogram_study(int number_x_bin, int number_y_bin,float &bins_occupied,float &av_distance,float &av_energy_bin,FloatVector &binx_vector,FloatVector &biny_vector,FloatVector &energy_bin,TH2* hgrid,float total_bins, float &tot_energy)
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
float calohit_propertiesAlgorithm::histogram_reader(int number_x_bin, int number_y_bin,TH2* hgrid)
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
float calohit_propertiesAlgorithm::rms_calculator(const FloatVector &energy_bin, float av_energy_bin)
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
calohit_propertiesAlgorithm::matrix::matrix(float element1, float element2, float element3) : 
    m_element1(element1),
    m_element2(element2),
    m_element3(element3)
{
}
//------------------------------------------------------------------------------------------------------------------------------------------
void calohit_propertiesAlgorithm::average_weighted(const FloatVector &binx_vector, const FloatVector &biny_vector, const FloatVector &energy_bin, float energy_sum, float &x_average_weighted,float &z_average_weighted)
{
    for(unsigned int element = 0; element<binx_vector.size(); element++)
    {
        x_average_weighted=x_average_weighted+((binx_vector.at(element)*energy_bin.at(element))/energy_sum);
        z_average_weighted=z_average_weighted+((biny_vector.at(element)*energy_bin.at(element))/energy_sum);	
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
void calohit_propertiesAlgorithm::covariance_matrix(const FloatVector &binx_vector, const FloatVector &biny_vector, const FloatVector &energy_bin, float x_average_weighted,float z_average_weighted,float energy_sum, matrix &covariance)
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
void calohit_propertiesAlgorithm::axes_calculator(const matrix &covariance, float &major_axis_weighted, float &minor_axis_weighted)
{
    float element_b_weighted=covariance.GetElement1()+covariance.GetElement3();
    float delta_weighted=(element_b_weighted*element_b_weighted)-4*((covariance.GetElement1()*covariance.GetElement3())-(covariance.GetElement2()*covariance.GetElement2()));
    major_axis_weighted=(element_b_weighted+sqrt(delta_weighted))/2;
    minor_axis_weighted=(element_b_weighted-sqrt(delta_weighted))/2;
}
//------------------------------------------------------------------------------------------------------------------------------------------
void calohit_propertiesAlgorithm::eigenvector_calculator(const matrix &covariance, float major_axis_weighted, float &eigenvector_major_weighted)
{
    float A_weighted = covariance.GetElement1()-major_axis_weighted;
    eigenvector_major_weighted = -A_weighted/covariance.GetElement2();
}
//------------------------------------------------------------------------------------------------------------------------------------------
void calohit_propertiesAlgorithm::distance_from_axis_calculator(float eigenvector_major_weighted,float x_average_weighted,float z_average_weighted,const FloatVector &binx_vector,const FloatVector &biny_vector, FloatVector &distance_from_axis_vector)
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
float calohit_propertiesAlgorithm::standard_deviation(const FloatVector &distance_from_axis_vector,float average_distance_from_axis)
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
float calohit_propertiesAlgorithm::crossed_energy_calculator(const FloatVector &binx_vector,const FloatVector &biny_vector,const FloatVector &energy_bin,float eigenvector_major_weighted,float x_average_weighted,float z_average_weighted,float howbigbinx, float howbigbinz)
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
float calohit_propertiesAlgorithm::rotation_method(int number_x_bin, int number_y_bin,TH2 * hgrid,TH2* hgrid4,float &rot_charge,float angle,float &unrotated_charge)
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
float calohit_propertiesAlgorithm::difference_charge(int number_x_bin, int number_y_bin, float rot_charge, float energy_sum)
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
StatusCode calohit_propertiesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace workshop_content


/* notes
        //PCA ANALYSIS//////////////////////////////////////////////////////////////////////////

		//float half_dimension_x=(float(number_x_bin)*(howbigbinx))/2;
		//float half_dimension_z=(float(number_y_bin)*(howbigbinz))/2;

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
*/
/*

for (const float element : floatVector)
{
    minusMeanFloatVecotr.push_back(element - average);
}

*/
