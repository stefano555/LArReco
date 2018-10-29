/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/pcaAlgorithm.h
 * 
 *  @brief  Header file for the calohit_properties algorithm class.
 * 
 *  $Log: $
 */
#ifndef WORKSHOP_pca_ALGORITHM_H
#define WORKSHOP_pca_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

namespace lar_reco
{

/**
 *  @brief  pcaAlgorithm class
 */
class pcaAlgorithm : public pandora::Algorithm
{
    /**
     *  @brief  Constructor
     */
    pcaAlgorithm();
    /**
     *  @brief  Destructor
     */
    ~pcaAlgorithm();

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
	std::string m_treeName;
	std::string m_fileName;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *pcaAlgorithm::Factory::CreateAlgorithm() const
{
    return new pcaAlgorithm();
}

} // namespace workshop_content

#endif // #ifndef WORKSHOP_calohit_properties_ALGORITHM_H
