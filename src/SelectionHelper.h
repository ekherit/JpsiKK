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

	std::vector<double>    kin_chi2;  //the kinematic fit result chi2 for all hypo
	std::vector<double>  mypid_chi2;
	std::vector<double>        prob; //the probability 

	void init(void)
	{
		channel = -1;
		success = false;
		chi2 = 1e100;
		mypid_chi2=1e100;
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


	SelectionHelper_t(double cme= PSIP_MASS, EvtRecTrackIterator END)
	{
		init();
		W = cme;
		end = END;
	}

	operator bool() const { return pass; }

	void setMyPid(SelectionHelper_t & kfp)
	{
		vector<double> & chi2 = kfp.mypid_chi2;
		chi2.resize(5);
		std::fill(chi2.begin(), chi2.end(), 0);
		for(int i=2;i<4;i++)
		{
			if(kfp.tracks[i]==kfp.end) continue;
			vector<double> chi2_tmp = get_chi2(kfp.tracks[i]);
			for(int pid =0;pid<5;pid++)
			{
				chi2[pid]+=chi2_tmp[pid];
			}
		}
	}


	bool passPid(void)
	{
		double & result = pass_pid;
		result = false;
		double & chi2 = mypid_chi2[channel]; //current chi2

		//global cut
		if( chi2 < 100)
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
		}
		return result;
	}

};
