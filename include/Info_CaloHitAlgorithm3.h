/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/Info_CaloHitAlgorithm.h
 * 
 *  @brief  Header file for the info_calohit algorithm class.
 * 
 *  $Log: $
 */
#ifndef WORKSHOP_INFO_CALOHIT_ALGORITHM3_H
#define WORKSHOP_INFO_CALOHIT_ALGORITHM3_H 1

#include "Pandora/Algorithm.h"

namespace lar_reco
{

/**
 *  @brief  Info_CaloHitAlgorithm class
 */
class Info_CaloHitAlgorithm3 : public pandora::Algorithm
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

    // Member variables here
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *Info_CaloHitAlgorithm3::Factory::CreateAlgorithm() const
{
    return new Info_CaloHitAlgorithm3();
}

} // namespace workshop_content

#endif // #ifndef WORKSHOP_INFO_CALOHIT_ALGORITHM2_H
