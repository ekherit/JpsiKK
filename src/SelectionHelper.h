// =====================================================================================
//
//       Filename:  SelectionHelper.h
//
//    Description:  Structure helps to select data by collection kinematic fit information
//                  and  other cuts
//
//        Version:  1.0
//        Created:  19.10.2015 20:56:53
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once

#include <vector>


#include "CLHEP/Vector/LorentzVector.h"
using CLHEP::HepLorentzVector;


#include "EvtRecEvent/EvtRecTrack.h"


#include "PhysConst.h"
#include "utils.h"
#include "Utils.h"
#include "SelectionConfig.h"
#include "Defs.h"


struct SelectionHelper_t
{
	int  channel;           //channel of the fit K, mu store here result of the selection
	bool pass_kinematic;    //pass kinematic cut
	bool pass_pid;          //pass pid cut
	bool pass_electron;     //pass electron cut
	bool pass;           		//total pass

	double W;               //center of mass energy,  GeV

	bool good_kinematic_fit; //result  of the kinematic fit 
	double kin_chi2;         //chi2
	std::vector<HepLorentzVector>  P;
	std::vector<EvtRecTrackIterator> tracks;
	EvtRecTrackIterator end;

	//std::vector<double>    kin_chi2;  //the kinematic fit result chi2 for all hypo
	std::vector<double>  mypid_chi2;
	std::vector<double>        prob; //the probability 

	void init(void)
	{
		channel = -1;
		good_kinematic_fit = false;
		kin_chi2 = 1e100;
		//mypid_chi2=1e100;
		W = PSIP_MASS;
		pass_kinematic = false;
		pass_pid = false;
		pass_electron = false;
		pass = false;
	}

	SelectionHelper_t(double cme= PSIP_MASS)
	{
		W = cme;
	}


	SelectionHelper_t(double cme, EvtRecTrackIterator END)
	{
		init();
		W = cme;
		end = END;
	}

	inline operator bool() const { return pass; }


	bool passElectrons(SelectionConfig & cfg)
	{
		double MIN_MOMENTUM[5] = { cfg.MIN_KAON_MOMENTUM,  cfg.MIN_MUON_MOMENTUM,  0, 0, 0}; 
		double MAX_MOMENTUM[5] = { cfg.MAX_KAON_MOMENTUM,  cfg.MAX_MUON_MOMENTUM,  0, 0, 0}; 
		double MIN_EP_RATIO[5] = { cfg.MIN_KAON_EP_RATIO,  cfg.MIN_MUON_EP_RATIO,  0, 0, 0}; 
		double MAX_EP_RATIO[5] = { cfg.MAX_KAON_EP_RATIO,  cfg.MAX_MUON_EP_RATIO,  0, 0, 0}; 
		pass_electron = false;
		for(int i=2;i<4;i++)
		{
			if(tracks[i]==end) continue;
			RecMdcTrack  * mdcTrk = (*tracks[i])->mdcTrack();
			RecEmcShower * emcTrk = (*tracks[i])->emcShower();
			double EpRatio = emcTrk->energy()/mdcTrk->p();
			if( 
					in(EpRatio, MIN_EP_RATIO[channel],  MAX_EP_RATIO[channel] ) 
					&&
					in(mdcTrk->p(), MIN_MOMENTUM[channel],  MAX_MOMENTUM[channel] )
				) pass_electron = true;
		}
		return pass_electron;
	}

	void setMyPid(void)
	{
		vector<double> & chi2 = mypid_chi2;
		chi2.resize(5);
		std::fill(chi2.begin(), chi2.end(), 0);
		for(int i=2;i<4;i++)
		{
			if(tracks[i]==end) continue;
			vector<double> chi2_tmp = get_chi2(tracks[i]);
			for(int pid =0;pid<5;pid++)
			{
				chi2[pid]+=chi2_tmp[pid];
			}
		}
	}


	bool passPid(SelectionConfig & cfg)
	{
		setMyPid();
		bool & result = pass_pid;
		result = false;
		double & chi2 = mypid_chi2[channel]; //current chi2

		//global cut
		if( chi2 > cfg.MAX_PID_CHI2)
		{
			result = false;
			return result;
		}

		switch(channel)
		{
			case ID_KAON:
				if( 
						chi2 <  mypid_chi2[ID_PION]   &&
						chi2 <  mypid_chi2[ID_MUON]
					)
				{
					result = true;
				}
					break;
			case ID_MUON:
				if( 
						chi2 <  mypid_chi2[ID_KAON]
					)
				{
					result = true;
				}
					break;
				break;
			default:
				break;
		}
		return result;
	}


	bool passKinematic(SelectionConfig & cfg)
	{
		pass_kinematic = kin_chi2 < cfg.MAX_KIN_CHI2 && good_kinematic_fit;
		return pass_kinematic;
	}


	bool totalPass(SelectionConfig & cfg)
	{
		pass = false;
		passKinematic(cfg)
		passElectrons(cfg);
		passPid(cfg);
		pass = passKinematic(cfg) && passElectrons(cfg) && passPid(cfg);
		return pass;
	}
};
