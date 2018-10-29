/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/calohit_propertiesAlgorithm.h
 * 
 *  @brief  Header file for the calohit_properties algorithm class.
 * 
 *  $Log: $
 */
#ifndef WORKSHOP_calohit_properties_ALGORITHM_H
#define WORKSHOP_calohit_properties_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "TString.h"
#include "TH2.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

namespace lar_reco
{

/**
 *  @brief  calohit_propertiesAlgorithm class
 */
class calohit_propertiesAlgorithm : public pandora::Algorithm
{
   
    /**
     *  @brief  Constructor
     */
    calohit_propertiesAlgorithm();
    /**
     *  @brief  Destructor
     */
    ~calohit_propertiesAlgorithm();

public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Run();
    class matrix
    {
    // Access specifier
    public:
        matrix(float element1, float element2, float element3);

        float GetElement1() const;
        void SetElement1(float element);
        float GetElement2() const;
        void SetElement2(float element);
        float GetElement3() const;
        void SetElement3(float element);

    // Data Members
    private:
        float m_element1;
        float m_element2;
        float m_element3;
    };	
    class calohitclass
    {
    // Access specifier
    public:
      	calohitclass(pandora::CaloHitList mylist);
	
	pandora::CaloHitList GetList() const;
    // Data Members
    private:
        pandora::CaloHitList m_mylist;
    };
    typedef std::map<const pandora::CaloHit *,calohitclass> MyMap;
    typedef std::map<const int,const pandora::CaloHit *> CaloMap;
    class simulated
    {
    // Access specifier
    public:
        simulated();

        pandora::IntVector GetPdg() const;
        void SetPdg(pandora::IntVector element);
        pandora::FloatVector GetWeight() const;
        void SetWeight(pandora::FloatVector element);

    // Data Members
    private:
	pandora::IntVector m_pdg;
	pandora::FloatVector m_weight;
    };

    /**
     *  @brief  Is track
     */
    int IsTrack(int id_centre);
    /**
     *  @brief  Is track
     */
    void vector_shifter(const pandora::FloatVector &drifttime_vector,const pandora::FloatVector &zplane_vector,pandora::FloatVector &drifttime_bis_vector,pandora::FloatVector &zplane_bis_vector,int centre);
    /**
     *  @brief  Is track
     */
    void histogram_filler_gaussian(const pandora::FloatVector &ex_vector,TH2* hgrid,const pandora::FloatVector &charge_vector,const pandora::FloatVector &drifttime_bis_vector,const pandora::FloatVector &zplane_bis_vector);
    /**
     *  @brief  Is track
     */
    void histogram_filler(TH2* hgrid, const pandora::FloatVector &drifttime_bis_vector,const pandora::FloatVector &zplane_bis_vector,const pandora::FloatVector &track_weight_vector);
    /**
     *  @brief  Is track
     */
    TString safe_name(const TString &initial_name);
    /**
     *  @brief  Is track
     */
    void histogram_study(int number_x_bin,int number_y_bin, float &bins_occupied,float &av_distance,float &av_energy_bin,pandora::FloatVector &binx_vector,pandora::FloatVector &biny_vector,pandora::FloatVector &energy_bin,TH2* hgrid,float total_bins, float &tot_energy);
    /**
     *  @brief  Is track
     */
    float histogram_reader(int number_x_bin, int number_y_bin,TH2* hgrid);
    /**
     *  @brief  Is track
     */
    float rms_calculator(const pandora::FloatVector &energy_bin, float av_energy_bin);
    /**
     *  @brief  it finds the weigthed average of the x and y positions which will be used to calculate the covariance matrix
     */
    void average_weighted(const pandora::FloatVector &binx_vector, const pandora::FloatVector &biny_vector, const pandora::FloatVector &energy_bin, float energy_sum, float &x_average_weighted,float &z_average_weighted);
    /**
     *  @brief  it calculate the covariance matrix: element 1 equals (1,1), element 2 equals (1,2) and (2,1), element 3 equals (2,2)
     */
    void covariance_matrix(const pandora::FloatVector &binx_vector, const pandora::FloatVector &biny_vector, const pandora::FloatVector &energy_bin,float x_average_weighted,float z_average_weighted, float energy_sum, matrix &covariance);
    /**
     *  @brief  it calculates the eigenvalues of the covariance matrix, called major axis and minor axis
     */
    void axes_calculator(const matrix &covariance, float &major_axis_weighted, float &minor_axis_weighted);
    /**
     *  @brief  it takes as input the covariance matrix and one eigenvalue and it calculates the eigenvector
     */
    void eigenvector_calculator(const matrix &covariance, float major_axis_weighted, float &eigenvector_major_weighted);
    /**
     *  @brief  it fills a vector with the distance of each single bin from the major pca axis passing through the centroid
     */
    void distance_from_axis_calculator(float eigenvector_major_weighted,float x_average_weighted,float z_average_weighted,const pandora::FloatVector &binx_vector,const pandora::FloatVector &biny_vector, pandora::FloatVector &distance_from_axis_vector);
    /**
     *  @brief  given a vector of elements and the average, it calcultas the standard deviation
     */
    float standard_deviation(const pandora::FloatVector &distance_from_axis_vector,float average_distance_from_axis);
    /**
     *  @brief  it returns the amount of energy of all the bin crossed by the pca major axis
     */
    float crossed_energy_calculator(const pandora::FloatVector &binx_vector,const pandora::FloatVector &biny_vector,const pandora::FloatVector &energy_bin, float eigenvector_major_weighted,float x_average_weighted,float z_average_weighted,float howbigbinx, float howbigbinz);
    /**
     *  @brief  it applies the rotation method
     */
    float rotation_method(int number_x_bin, int number_y_bin,TH2* hgrid,TH2* hgrid4,float &rot_charge,float angle,float &unrotated_charge);
    /**
     *  @brief  it returns the amount of energy lost during rotation
     */
    float difference_charge(int number_x_bin, int number_y_bin, float rot_charge, float energy_sum);
    /**
     *  @brief  it returns the amount of energy lost during rotation
     */
    void IsParentAMuon(const pandora::MCParticle *pMCParticle, bool &hasParentMuon);
    /**
     *  @brief  
     */
    void chisquare_calculator(const float slope,const float calo_centroid_x,const float calo_centroid_z,const float wire_pitch, float &chisquare,const MyMap &calohitmap,const int centre, const CaloMap &inttocalo);
    /**
     *  @brief  it calculates the chi square value of the pca major axis
     */
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
	std::string m_treeName;
	int m_file_counter;
	std::string m_fileName;
	
        int m_eventNumber;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float calohit_propertiesAlgorithm::matrix::GetElement1() const 
{
    return m_element1;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void calohit_propertiesAlgorithm::matrix::SetElement1(float element) 
{
    m_element1 = element;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float calohit_propertiesAlgorithm::matrix::GetElement2() const 
{
    return m_element2;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void calohit_propertiesAlgorithm::matrix::SetElement2(float element) 
{
    m_element2 = element;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float calohit_propertiesAlgorithm::matrix::GetElement3() const 
{
    return m_element3;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void calohit_propertiesAlgorithm::matrix::SetElement3(float element) 
{
    m_element3 = element;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline pandora::IntVector calohit_propertiesAlgorithm::simulated::GetPdg() const
{
    return m_pdg;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void calohit_propertiesAlgorithm::simulated::SetPdg(pandora::IntVector element)
{
    m_pdg = element;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline pandora::FloatVector calohit_propertiesAlgorithm::simulated::GetWeight() const
{
    return m_weight;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void calohit_propertiesAlgorithm::simulated::SetWeight(pandora::FloatVector element)
{
    m_weight = element;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline pandora::CaloHitList calohit_propertiesAlgorithm::calohitclass::GetList() const
{
    return m_mylist;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline pandora::Algorithm *calohit_propertiesAlgorithm::Factory::CreateAlgorithm() const
{
    return new calohit_propertiesAlgorithm();
}

} // namespace workshop_content

#endif // #ifndef WORKSHOP_calohit_properties_ALGORITHM_H
