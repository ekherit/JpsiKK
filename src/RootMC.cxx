// =====================================================================================
//
//       Filename:  RootMC.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 17:03:06
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootMC.h"

StatusCode RootMC::init_tuple(void)
{
  StatusCode status;
  status = tuple->addItem ("psip_decay", psip_decay); //flag for psip decay
  status = tuple->addItem ("jpsi_decay", jpsi_decay); //flag for jpsi decay 
  status = tuple->addItem ("KK", KK);               //KK event
  status = tuple->addItem ("uu", uu);               //mu mu event
  status = tuple->addItem ("oo", oo);               //other event
  status = tuple->addItem ("ntrack", ntrack,0,4); 
  status = tuple->addIndexedItem ("id",    ntrack, pid);
  status = tuple->addIndexedItem ("q",     ntrack, q);
  status = tuple->addIndexedItem ("E",     ntrack, E);
  status = tuple->addIndexedItem ("p",     ntrack, p);
  status = tuple->addIndexedItem ("px",    ntrack, px);
  status = tuple->addIndexedItem ("py",    ntrack, py);
  status = tuple->addIndexedItem ("pz",    ntrack, pz);
  status = tuple->addIndexedItem ("pt",    ntrack, pt);
  status = tuple->addIndexedItem ("theta", ntrack, theta);
  status = tuple->addIndexedItem ("phi",   ntrack, phi);
  return status;
}

void RootMC::init(void)
{
  ntrack=4;
}

void RootMC::fill(const std::vector<CLHEP::HepLorentzVector> & Pkf,  Event::McParticleCol * mcParticleCol)
{
	//Fill my mc truth information
	ntrack=4;
	HepLorentzVector MCPpion[2];
	HepLorentzVector MCPkaon_or_muon[2];
	bool pi_minus(false);
	bool pi_plus(false);
	bool K_minus(false);
	bool K_plus(false);
	bool mu_minus(false);
	bool mu_plus(false);
	int mytrack=0;
	psip_decay = 0;
	jpsi_decay = 0;
	KK = 0;
	uu = 0;
	oo = 0;
	bool psipDecay = false;
	int rootIndex = -1;
	Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
	for (iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
	{
		if ((*iter_mc)->primaryParticle()) continue;
		if (!(*iter_mc)->decayFromGenerator()) continue;
		//if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
		if ((*iter_mc)->particleProperty()==100443)
		{
			psipDecay = true;
			psip_decay=1;
			rootIndex = (*iter_mc)->trackIndex();
		}
		if (!psipDecay) continue;
		if ((*iter_mc)->particleProperty()==443)
		{
			jpsi_decay=1;
		}
		if (jpsi_decay!=1) continue;
		if((*iter_mc)->particleProperty() == +211) 
		{
			pi_plus=true;
			MCPpion[1] = (*iter_mc)->initialFourMomentum();
			pid[1]=211;
			mytrack++;
		}
		if((*iter_mc)->particleProperty() == -211) 
		{
			MCPpion[0] = (*iter_mc)->initialFourMomentum();
			pi_minus=true;
			pid[0]=-211;
			mytrack++;
		}
		if( ! pi_plus && !pi_minus) continue; //keep only psip to Jpsi pi pi decay
		switch((*iter_mc)->particleProperty())
		{
			case -13:
				MCPkaon_or_muon[0] = (*iter_mc)->initialFourMomentum();
				mu_minus = true;
				mytrack++;
				pid[2]=-13;
				break;
			case +13:
				MCPkaon_or_muon[1] = (*iter_mc)->initialFourMomentum();
				mu_plus = true;
				mytrack++;
				pid[3]=13;
				break;
			case -321:
				MCPkaon_or_muon[0] = (*iter_mc)->initialFourMomentum();
				K_minus=true;
				mytrack++;
				pid[2]=-321;
				break;
			case +321:
				MCPkaon_or_muon[1] = (*iter_mc)->initialFourMomentum();
				mytrack++;
				K_plus = true;
				pid[3]=321;
				break;
		};
	}
	if(K_plus && K_minus) 
	{
		KK=1;
		oo=0;
	}
	if(mu_plus && mu_minus) 
	{
		uu=1;
		oo=0;
	}
	if(KK==1 && uu==1) oo=1;
	if(mytrack!=4) oo=1;
	if(KK==1 || uu==1)
	{
		vector<HepLorentzVector> P(4);
		P[0]=MCPpion[0];
		P[1]=MCPpion[1];
		P[2]=MCPkaon_or_muon[0];
		P[3]=MCPkaon_or_muon[1];
		for(int i=0;i<4;i++)
		{
			q[i] = 0; 
			E[i] = P[i].e();
			p[i] = P[i].rho();
			px[i]= P[i].px();
			py[i]= P[i].py();
			pz[i]= P[i].pz();
			pt[i]= P[i].perp();
			theta[i]= P[i].theta();
			phi[i] = P[i].phi();
		}
	}
	else
	{
		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			if ((*iter_mc)->primaryParticle()) continue;
			if (!(*iter_mc)->decayFromGenerator()) continue;
			HepLorentzVector P = (*iter_mc)->initialFourMomentum();
			Hep3Vector p_mc = P.vect();
			for(int i=0;i<4;i++)
			{
				Hep3Vector p_rec = Pkf[i].vect();
				Hep3Vector dp = p_rec - p_mc;
				if(dp.mag()/std::min(p_rec.mag(), p_mc.mag()) < 0.05)
				{
					pid[i] = (*iter_mc)->particleProperty(); 
					//q[i] = 0; 
					E[i] = P.e();
					p[i] = P.rho();
					px[i]= P.px();
					py[i]= P.py();
					pz[i]= P.pz();
					pt[i]= P.perp();
					theta[i]= P.theta();
					phi[i] = P.phi();
				}
			}
		}
	}
}
