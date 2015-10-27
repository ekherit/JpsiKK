// =====================================================================================
//
//       Filename:  RootMCTopo.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 17:03:52
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootMCTopo.h"
StatusCode RootMCTopo::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem("indexmc", m_idxmc, 0, 100);
  status = tuple->addIndexedItem("pdgid", m_idxmc, m_pdgid);
  status = tuple->addIndexedItem("motheridx", m_idxmc, m_motheridx);
  status = tuple->addIndexedItem("idx", m_idxmc, m_idx);
  status = tuple->addItem("hash", m_hash);
  return status;
}

void RootMCTopo::init(void)
{
  m_idxmc=0;
}


virtual void fill(Event::McParticleCol * mcParticleCol)
{
	//check the MC information
	//if(!mcParticleCol)
	//{
	//	log << MSG::ERROR << "Could not retrieve McParticelCol" << endreq;
	//	return StatusCode::FAILURE;
	//}
	//Fill MC TOPO INFORMATION
	bool psipDecay = false;
	int rootIndex = -1;
	Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
	int m_numParticle = 0;
	for (; iter_mc != mcParticleCol->end(); iter_mc++)
	{
		if ((*iter_mc)->primaryParticle()) continue;
		if (!(*iter_mc)->decayFromGenerator()) continue;
		if ((*iter_mc)->particleProperty()==100443)
		{
			psipDecay = true;
			rootIndex = (*iter_mc)->trackIndex();
		}
		if (!psipDecay) continue;
		int pdgid = (*iter_mc)->particleProperty();
		int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
		fMCTopo.m_pdgid[m_numParticle] = pdgid;
		fMCTopo.m_motheridx[m_numParticle] = mcidx;
		fMCTopo.m_idx[m_numParticle] = (*iter_mc)->trackIndex()-rootIndex;
		fMCTopo.m_hash=0; //no hash calculation now
		m_numParticle += 1;
	}
	fMCTopo.m_idxmc = m_numParticle;
}
