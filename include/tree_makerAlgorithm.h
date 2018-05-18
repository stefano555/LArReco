/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/Info_CaloHitAlgorithm.h
 * 
 *  @brief  Header file for the info_calohit algorithm class.
 * 
 *  $Log: $
 */
#ifndef WORKSHOP_TREE_MAKER_ALGORITHM_H
#define WORKSHOP_TREE_MAKER_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_reco
{

/**
 *  @brief  Info_CaloHitAlgorithm class
 */
class tree_makerAlgorithm : public pandora::Algorithm
{
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
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
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

    // Member variables here
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

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float tree_makerAlgorithm::matrix::GetElement1() const 
{
    return m_element1;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void tree_makerAlgorithm::matrix::SetElement1(float element) 
{
    m_element1 = element;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float tree_makerAlgorithm::matrix::GetElement2() const 
{
    return m_element2;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void tree_makerAlgorithm::matrix::SetElement2(float element) 
{
    m_element2 = element;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float tree_makerAlgorithm::matrix::GetElement3() const 
{
    return m_element3;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void tree_makerAlgorithm::matrix::SetElement3(float element) 
{
    m_element3 = element;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline pandora::Algorithm *tree_makerAlgorithm::Factory::CreateAlgorithm() const
{
    return new tree_makerAlgorithm();
}

} // namespace workshop_content

#endif // #ifndef WORKSHOP_TREE_MAKER_ALGORITHM_H
