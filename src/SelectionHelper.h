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
#include <iostream>
#include <limits>


using namespace std;


#include "CLHEP/Vector/LorentzVector.h"
using CLHEP::HepLorentzVector;


#include "EvtRecEvent/EvtRecTrack.h"
#include "ParticleID/ParticleID.h"


#include "PhysConst.h"
#include "utils.h"
#include "Utils.h"
#include "SelectionConfig.h"
#include "Defs.h"
#include "KinematicFit.h"


extern  std::vector<KinematicFit_t> kinfit(const std::vector<EvtRecTrackIterator> & Tracks,  const double CENTER_MASS_ENERGY);


struct SelectionHelper_t
{
	const SelectionConfig * cfg;
	int  channel;           //channel of the fit K, mu store here result of the selection
	bool pass_kinematic;    //pass kinematic cut
	bool pass_pid;          //pass pid cut
	bool pass_electron;     //pass electron cut
	bool pass;           		//total pass



	//bool good_kinematic_fit; //result  of the kinematic fit 
	//double kin_chi2;         //chi2
	//std::vector<HepLorentzVector>  P;
	std::vector<EvtRecTrackIterator> tracks;
	//EvtRecTrackIterator end;

	std::vector <KinematicFit_t> KF;  //kinematic fit for different hypo
	//std::vector<double>    kin_chi2;
	std::vector<double>  mypid_chi2;
	std::vector<double>    pid_chi2; //ParticleId chi2
	std::vector<double>        prob; //the probability 


	void init(void)
	{
		channel = -1;
		//good_kinematic_fit = false;
		//kin_chi2 = 1e100;
		//mypid_chi2=1e100;
		pass_kinematic = false;
		pass_pid = false;
		pass_electron = false;
		pass = false;
	}

	SelectionHelper_t(const SelectionConfig & c) : cfg(&c)
	{
		init();
	}

	inline operator bool() const { return pass; }


	bool passElectrons(void)
	{
		double MIN_MOMENTUM[5] = { cfg->MIN_KAON_MOMENTUM,  cfg->MIN_MUON_MOMENTUM,  0, 0, 0}; 
		double MAX_MOMENTUM[5] = { cfg->MAX_KAON_MOMENTUM,  cfg->MAX_MUON_MOMENTUM,  0, 0, 0}; 
		double MIN_EP_RATIO[5] = { cfg->MIN_KAON_EP_RATIO,  cfg->MIN_MUON_EP_RATIO,  0, 0, 0}; 
		double MAX_EP_RATIO[5] = { cfg->MAX_KAON_EP_RATIO,  cfg->MAX_MUON_EP_RATIO,  0, 0, 0}; 
		pass_electron = false;
		for(int i=2;i<tracks.size();i++)
		{
			//if(tracks[i]==end) continue;
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
		for(int i=2;i<tracks.size();i++)
		{
			//if(tracks[i]==end) continue;
			vector<double> chi2_tmp = get_chi2(tracks[i]);
			for(int pid =0;pid<5;pid++)
			{
				chi2[pid]+=chi2_tmp[pid];
			}
		}
	}

	void setPid(void)
	{
		ParticleID * PID = ParticleID::instance();
		PID->init();
		vector<double> & chi2 = pid_chi2;
		chi2.resize(5);
		std::fill(chi2.begin(), chi2.end(), 0);
		for(int i=2;i<tracks.size();i++)
		{
			//if(tracks[i]==end) continue;
			PID->setRecTrack((*tracks[i]));
			PID->setMethod(PID->methodProbability());
			PID->setChiMinCut(4);
			PID->usePidSys(PID->useDedx() | PID->useTof1() | PID->useTof2());
			PID->identify(PID->all()); 
			PID->calculate();
			chi2[ID_KAON]     += PID->chi(3);
			chi2[ID_MUON]     += PID->chi(1);
			chi2[ID_PION]     += PID->chi(2);
			chi2[ID_ELECTRON] += PID->chi(0);
			chi2[ID_PROTON]   += PID->chi(4);
		}
	}

	bool passPid(const vector<double> & pchi2)
	{
		const double & chi2 = pchi2[channel]; //current chi2

		if( chi2 > cfg->MAX_PID_CHI2)
		{
			return false;
		}

		switch(channel)
		{
			case ID_KAON:
				if( 
						chi2 <  pchi2[ID_MUON]
//           &&chi2 <  pchi2[ID_PION]
					)
				{
					pass_pid = true;
				}
				break;
			case ID_MUON:
				if( 
						chi2 <  pchi2[ID_KAON]
					)
				{
					pass_pid = true;
				}
				break;
				break;
			default:
				break;
		}
		return pass_pid;
	}

	bool passMyPid(void)
	{
		setMyPid();
		return passPid(mypid_chi2);
	}

	bool passPid(void)
	{
		setPid();
		return passPid(pid_chi2);
	}




	bool passKinematic(void)
	{
		pass_kinematic = KF[channel].chi2 < cfg->MAX_KIN_CHI2 && KF[channel].success;
		return pass_kinematic;
	}



	void kinfit(
			TrackPair_t & pion_pair,
			EvtRecTrackIterator & track
			)
	{
		tracks.resize(3);
		tracks[0] = pion_pair.first;
		tracks[1] = pion_pair.second;
		tracks[2] = track;
		KF = ::kinfit(tracks,  cfg->CENTER_MASS_ENERGY);
	}

	void select_channel_by_kinematic_fit(void)
	{
		std::list<int> pid_list;
		pid_list.push_back(ID_KAON);
		pid_list.push_back(ID_MUON);
		//pid_list.push_back(ID_ELECTRON);
		//pid_list.push_back(ID_PROTON);
		channel = ID_PROTON; //preselect some hypotesa
		for(std::list<int>::iterator pid = pid_list.begin(); pid!=pid_list.end() ; pid++)
		{
			if(KF[*pid].chi2 < KF[channel].chi2)
			{
				channel = *pid;
			}
		}
		
		//channel = 4; //preselect some hypotesa
		//for(int i=0;i<2;i++)
		//{
		//	if(KF[i].chi2 < KF[channel].chi2)
		//	{
		//		channel = i;
		//	}
		//}
	}

	void select_channel_by_kinematic_fit_and_pid(void)
	{
		setPid();
		pass_pid = false;
		pass_kinematic = false;
		vector<double> chi2(5);
		for(int i=0;i<chi2.size();i++)
		{
			chi2[i] = pid_chi2[i] + KF[i].chi2;
		}


		//select best channel
		channel = 4;
		for(int i=0;i<5;i++)
		{
			if(chi2[i] < chi2[channel])
			{
				channel = i;
			}
		}
	}

	bool passKinPid(void)
	{
		if( KF[channel].chi2 + pid_chi2[channel] > cfg->MAX_KIN_CHI2 ||  !KF[channel].success)
		{
			pass_pid = false;
			pass_kinematic =false;
			return false;
		}
		pass_kinematic = true;
		pass_pid = true;
		//P = KF[channel].P;
		return true;
	}

	bool totalPass(void)
	{
		pass = false;
		select_channel_by_kinematic_fit(); //after this we allways has channel
		passKinematic();
		passElectrons();
		passPid();
		//clog << "channel = " << channel << endl;
		//if(!passKinematic()) return false;
		//clog << "pass_kinematic = " << pass_kinematic << endl;
		//if(!passElectrons()) return false;
		//clog << "pass_electrons = " << pass_electron << endl;
		//if(!passPid()) return false;
		//clog << "pass_pid = " << pass_pid << endl;
		//clog << "kin_chi2 = " << KF[channel].chi2 << " " ;
		//clog << "pass_pid: " << pass_pid << "  pid_chi2 = " <<  pid_chi2[channel] << " " ;
		//clog << "pass_electron: " << pass_electron << endl;
		pass = pass_kinematic && pass_electron && pass_pid;
		return pass;
	}

	bool dummyPass(void)
	{
		select_channel_by_kinematic_fit();
		setPid();
		pass = true;
		return pass;
	}

	double getKinChi2(int pid) const
	{
		return KF[pid].chi2;
	}

	double getPidChi2(int pid) const
	{
		return pid_chi2[pid];
	}

	std::vector<HepLorentzVector> getMomentum(int pid) const
	{
		return KF[pid].P;
	}

};
