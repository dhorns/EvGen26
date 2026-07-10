/*
 * EvGen26
 *
 * Object-oriented version of EvGen.cxx.
 *
 * DLH	May 5, 2026			First Version
 * DLH	July 8, 2026		Adding Li-6 and Li-7
 *
 */

#include "TF1.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "Riostream.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "physics.h"

class EvGen26
{

	public:

		EvGen26(const char* paramFile = "par/EvGen26.in");
		int Run();

	private:

		struct ABCdata
		{
			Double_t k = 0.0;
			Double_t a = 0.0;
			Double_t b = 0.0;
			Double_t c = 0.0;
		};

		struct Params
		{
			ULong_t counts = static_cast<ULong_t>(1e5);
			Double_t tgt_len = 5.00;
			Double_t vtx_rad = 0.50;
			Double_t sig_bm = 0.5;
			Double_t p_gamma = 0.0;
			Double_t Sigma = 0.0;
			Double_t phi_offset = 0.0;
			Double_t e_lo = 200.0;
			Double_t e_bite = 1.0;
			Double_t e_step = 10.0;
			UInt_t e_bins = 11;
			TString process = "compton";
			TString tgt = "p";
			Int_t chan_lo = 0;
			Int_t chan_hi = 351;
			TString te_file = "tageng855_new.dat";
		};

	private:

		TString fParamFile;
		Params fParam;
		std::vector<ABCdata> fABC;
		Double_t fEnergy[328] = {};
		Double_t fDEnergy[328] = {};

	private:

		Bool_t ReadParams();
		void ReadABC(TString process);
		Double_t ScatCTH(Double_t* x, Double_t* par) const;
		Double_t PhiDist(Double_t* x, Double_t* par) const;
		Double_t RecoilP(Double_t ke, Double_t mtgt, Double_t th, Double_t mpi) const;
		TString GenNames(UInt_t npart, UInt_t* ptag) const;

		static Double_t Sqr(Double_t x) { return x*x; }
};

EvGen26::EvGen26( const char* paramFile) : fParamFile( paramFile)
{
}

void EvGen26::ReadABC( TString process)
{

	TString infile;

	if (process == "compton") infile = "par/abc_compton.dat";
	else if (process == "pi0") infile = "par/abc_pi0.dat";
	else if (process == "eta") infile = "par/abc_eta.dat";

	std::ifstream inFile(infile);
	if (!inFile.is_open())
	{
		std::cout << "Error opening file " << infile << std::endl;
		std::exit(-1);
	}

	fABC.clear();
	ABCdata entry;
	while (inFile >> entry.k >> entry.a >> entry.b >> entry.c)
	{
		fABC.push_back( entry);
	}

	if ( fABC.empty())
	{
		std::cout << "No ABC parameters found in " << infile << std::endl;
		std::exit(-1);
	}

	std::cout << " ABC parameters read in..." << std::endl;
}

Bool_t EvGen26::ReadParams()
{
	Bool_t eflag, cflag;
	UInt_t i;
	Double_t eff, deff, jnk;
	TString string, name;

	eflag = cflag = kFALSE;

	// Defaults are already set in Params, but reset them here in case ReadParams()
	// is ever called more than once on the same object.
	fParam = Params();

	name = fParamFile;
	std::ifstream inFile(name);
	if (!inFile.is_open())
	{
		std::cout << "Error opening file " << name << std::endl;
		std::exit(1);
	}

	while (!inFile.eof())
	{
		name.ReadLine(inFile);
		if (name[0] != '#')
		{
			string = "Process: ";
			if (name.Contains(string)) fParam.process = name.Remove(0, string.Length());

			string = "Target: ";
			if (name.Contains(string)) fParam.tgt = name.Remove(0, string.Length());

			string = "TargetLength: ";
			if (name.Contains(string))
			{
				name.Remove(0, string.Length());
				fParam.tgt_len = name.Atof();
			}

			string = "BeamSpotRadius: ";
			if (name.Contains(string))
			{
				name.Remove(0, string.Length());
				fParam.vtx_rad = name.Atof();
			}

			string = "BeamSpotSigma: ";
			if (name.Contains(string))
			{
				name.Remove(0, string.Length());
				fParam.sig_bm = name.Atof();
			}

			string = "BeamPol: ";
			if (name.Contains(string))
			{
				name.Remove(0, string.Length());
				fParam.p_gamma = name.Atof();
			}

			string = "PhotAsym: ";
			if (name.Contains(string))
			{
				name.Remove(0, string.Length());
				fParam.Sigma = name.Atof();
			}

			string = "PhiOffset: ";
			if (name.Contains(string))
			{
				name.Remove(0, string.Length());
				fParam.phi_offset = name.Atof()*kD2R;
			}

			string = "Throws: ";
			if (name.Contains(string))
			{
				name.Remove(0, string.Length());
				fParam.counts = name.Atoi();
			}

			string = "Energy: ";
			if (name.Contains(string))
			{
				eflag = kTRUE;
				name.Remove(0, string.Length());
				string = name;
				string.Remove(0, string.Last(' '));
				fParam.e_bins = string.Atoi();

				string = name;
				string.Remove(string.First(' '));
				fParam.e_bite = string.Atof();

				string = name;
				string.Remove(string.Last(' '));
				string.Remove(0, string.Last(' '));
				fParam.e_step = string.Atof();

				string = name;
				string.Remove(string.Last(' '));
				string.Remove(string.Last(' '));
				string.Remove(0, string.First(' '));
				fParam.e_lo = string.Atof();
			}
			string = "TagEngFile: ";
			if (name.Contains(string))
			{
				cflag = kTRUE;
				fParam.te_file = name.Remove(0, string.Length());
				fParam.te_file.Prepend("par/");
			}
			string = "Channels: ";
			if ( name.Contains( string))
			{
				name.Remove(0, string.Length());
				string = name;
				string.Remove(0, string.Last(' '));
				fParam.chan_hi = string.Atoi();

				string = name;
				string.Remove(string.First(' '));
				fParam.chan_lo = string.Atoi();
			}
		}
	}
	inFile.close();

	if ((eflag == kTRUE) && (cflag == kTRUE))
	{
		std::cout << "You have chosen BOTH tagger channels and energies.\n";
		std::cout << "Pick one or the other.\n";
		std::exit(-1);
	}
	else if ((eflag == kFALSE) && (cflag == kFALSE)) {
		std::cout << "You have chosen neither tagger channels nor energies.\n";
		std::cout << "Pick one or the other.\n";
		std::exit(-1);
	}

	if ((fParam.process != "compton") && (fParam.process != "pi0") && (fParam.process != "eta")) {
		std::cout << "Invalid Process String.\n";
		std::exit(-1);
	}

	if ((fParam.tgt != "p") && (fParam.tgt != "c") && (fParam.tgt != "w")
				&& (fParam.tgt != "he3") && (fParam.tgt != "he4")
				&& (fParam.tgt != "li6") && (fParam.tgt != "li7")) {
		std::cout << "Invalid Target String.\n";
		std::exit(-1);
	}

	name = Form("Process = %s", (const char*) fParam.process);
	std::cout << name << std::endl;
	name = Form("Target = %s", (const char*) fParam.tgt);
	std::cout << name << std::endl;
	name = Form("Target Length = %4.2f cm", fParam.tgt_len);
	std::cout << name << std::endl;
	name = Form("Beam Spot Radius = %4.2f cm", fParam.vtx_rad);
	std::cout << name << std::endl;
	name = Form("Beam Spot Sigma = %4.2f cm", fParam.sig_bm);
	std::cout << name << std::endl;
	name = Form("Beam Polarization = %4.2f", fParam.p_gamma);
	std::cout << name << std::endl;
	name = Form("Photon Asymmetry = %4.2f", fParam.Sigma);
	std::cout << name << std::endl;
	name = Form("Phi Offset = %5.1f", fParam.phi_offset/kD2R);
	std::cout << name << std::endl;
	name = Form("Throws = %ld", fParam.counts);
	std::cout << name << std::endl;

	if (eflag == kTRUE)
	{
		name = Form(" Energy Range = %5.1f - %5.1f MeV in %4.1f MeV steps",
					fParam.e_lo, fParam.e_lo + fParam.e_step*(fParam.e_bins - 1),
					fParam.e_step);
		std::cout << name << std::endl;
		name = Form(" Energy bite %4.1f MeV", fParam.e_bite);
		std::cout << name << std::endl;
	}
	else if (cflag == kTRUE)
	{
		if ((fParam.chan_lo < 0) || (fParam.chan_hi > 328) || (fParam.chan_lo > fParam.chan_hi))
		{
			std::cout << "Invalid Tagger Channel Range.";
			std::cout << "Must be 0-351";
			std::cout << std::endl;
			std::exit(-1);
		}
		name = Form("Using tagger energy file %s", (const char*) fParam.te_file);
		std::cout << name << std::endl;
		name = Form("Tagger channels %d - %d", fParam.chan_lo, fParam.chan_hi);
		std::cout << name << std::endl;

		// Read in tagger energies from specified file
		std::ifstream teFile( fParam.te_file);
		if (!teFile.is_open())
		{
			std::cout << "Error in tagger energy file " << fParam.te_file << std::endl;
			std::exit(-1);
		}
		while (teFile >> i >> jnk >> eff >> deff >> jnk)
		{
			if (i < 328)
			{
				fEnergy[i] = eff;
				fDEnergy[i] = deff/2;
			}
		}
		teFile.close();
	}

	return eflag;
}

Double_t EvGen26::ScatCTH(Double_t* x, Double_t* par) const
{

	UInt_t i;
	Double_t scat_cth;
	Double_t aa, bb, cc;
	Double_t mpi, th, sth;
	Double_t q, ff;
	Double_t cth, ke;
	Double_t ap, alpha, mtgt;
	Double_t a, b, c;

	ap = 0;
	alpha = 0;
	mtgt = 0;

	cth = x[0];
	ke = par[0];

	scat_cth = 1;		// Default is phase space

	// For the RecoilP function
	if (fParam.process == "compton") mpi = 0;
	else if (fParam.process == "pi0") mpi = kMPI0_MEV;
	else if (fParam.process == "eta") mpi = kMETA_MEV;
	else mpi = 0;

	if (fParam.tgt == "p")
	{
		aa = 1;
		bb = 0;
		cc = 0;

		if (!fABC.empty())
		{
			if (ke < fABC.front().k)
			{
				aa = fABC.front().a;
				bb = fABC.front().b;
				cc = fABC.front().c;
			}
			if (ke > fABC.back().k)
			{
				aa = fABC.back().a;
				bb = fABC.back().b;
				cc = fABC.back().c;
			}
			for (i = 0; i + 1 < fABC.size(); i++)
			{
				if ((ke >= fABC[i].k) && (ke <= fABC[i+1].k))
				{
					aa = fABC[i].a + (fABC[i+1].a-fABC[i].a)/(fABC[i+1].k-fABC[i].k)
									*(ke - fABC[i].k);
					bb = fABC[i].b + (fABC[i+1].b-fABC[i].b)/(fABC[i+1].k-fABC[i].k)
									*(ke - fABC[i].k);
					cc = fABC[i].c + (fABC[i+1].c-fABC[i].c)/(fABC[i+1].k-fABC[i].k)
									*(ke - fABC[i].k);
					break;
				}
			}
		}
		scat_cth = aa + bb*cth + cc*Sqr(cth);
	}
	else if ((fParam.tgt == "c") || (fParam.tgt == "w") || (fParam.tgt == "o"))
	{
		sth = std::sqrt(1 - Sqr(cth));
		th = std::acos(cth);
		if ((fParam.tgt == "c") || (fParam.tgt == "w"))
		{
			ap = 1.66;
			alpha = 1.3333;
			mtgt = kM_C12_MEV;
		}
		else if (fParam.tgt == "o")
		{
			ap = 1.76;
			alpha = 2;
			mtgt = kM_O16_MEV;
		}
		q = RecoilP(ke, mtgt, th, mpi);
		ff = (1-alpha/2/(2+3*alpha)*Sqr(q*ap))*std::exp(-0.25*Sqr(q*ap));

		scat_cth = Sqr(sth*ff);
	}
	else if (fParam.tgt == "he3")
	{
		sth = std::sqrt(1 - Sqr(cth));
		th = std::acos(cth);
		b = 0.528;
		mtgt = kM_HE3_MEV;
		q = RecoilP(ke, mtgt, th, mpi);
		ff = std::exp(-b*q*q);

		scat_cth = Sqr(sth*ff);
	}
	else if (fParam.tgt == "he4")
	{
		sth = std::sqrt(1 - Sqr(cth));
		th = std::acos(cth);
		b = 0.352;
		mtgt = kM_HE4_MEV;
		q = RecoilP(ke, mtgt, th, mpi);
		ff = std::exp(-b*q*q);

		scat_cth = Sqr(sth*ff);
	}
	else if (fParam.tgt == "li6")
	{
		sth = std::sqrt(1 - Sqr(cth));
		th = std::acos(cth);
		a = 0.391;
		b = 1.063;
		c = 1.328;
		mtgt = kM_LI6_MEV;
		q = RecoilP(ke, mtgt, th, mpi);
		ff = (2*std::exp(-a*q*q)+(1-b*q*q)*std::exp(-c*q*q))/3;

		scat_cth = Sqr(sth*ff);
	}
	else if (fParam.tgt == "li7")
	{
		sth = std::sqrt(1 - Sqr(cth));
		th = std::acos(cth);
		a = 0.180;
		b = 0.694;
		mtgt = kM_LI7_MEV;
		q = RecoilP(ke, mtgt, th, mpi);
		ff = (1 - a*q*q)*std::exp(-b*q*q);

		scat_cth = Sqr(sth*ff);
	}
	return scat_cth;
}

Double_t EvGen26::PhiDist(Double_t* x, Double_t* par) const
{

	Double_t ph, p_g, sig, phi0, xs;

	ph = x[0];
	p_g = par[0];
	sig = par[1];
	phi0 = par[2];

	xs = 1 + p_g*sig*std::cos(2*(ph + phi0));

	return xs;
}

Double_t EvGen26::RecoilP(Double_t ke, Double_t mtgt, Double_t th, Double_t mpi) const
{

	TVector3 k, q, p;
	Double_t qp;
	Double_t W, A, B, C, aa, bb, cc;

	W = ke + mtgt;

	A = 2*ke*mtgt + mpi*mpi;
	B = 2*ke*std::cos(th);
	C = 2*W;

	aa = C*C - B*B;
	bb = -2*A*B;
	cc = Sqr(C*mpi) - A*A;

	qp = (-bb + std::sqrt(bb*bb - 4*aa*cc))/2/aa;
	if (qp != qp) std::cout << "<E> invalid pion momentum" << std::endl;

	k.SetXYZ(0, 0, ke);
	q.SetXYZ(qp*std::sin(th), 0, qp*std::cos(th));
	p = k - q;

	return p.Mag()/kHBARC;
}

TString EvGen26::GenNames(UInt_t npart, UInt_t* ptag) const
{
	UInt_t i, j;

	TString pstr[] = {"Px", "Py", "Pz", "Pt", "En"};
	TString beam = "X_vtx:Y_vtx:Z_vtx:Px_bm:Py_bm:Pz_bm:Pt_bm:En_bm";
	TString particles;
	TString names;

	for (i = 0; i < npart; i++)
	{
		for (j = 0; j < 5; j++)
		{
			particles.Append(pstr[j]);
			if ((i == (npart-1)) && (j == 4)) particles.Append(Form("_l%02d%02d", i+1, ptag[i]));
			else particles.Append(Form("_l%02d%02d:", i+1, ptag[i]));
		}
	}

	names = beam + ":" + particles;

	return names;
}

int EvGen26::Run()
{
	Bool_t flag;
	UInt_t i, j, stop;
	UInt_t npart, update;
	UInt_t ptag[10];
	
	Double_t ke, kthr;
	Double_t e_mid, e_low, e_high, e_bite;
	Double_t vtx_x, vtx_y, vtx_z;
	Double_t pm, pe_cm;
	Double_t qm, qth_cm, qph_cm, mom_cm;
	Double_t qx_cm, qy_cm, qz_cm, qe_cm;
	Double_t gammae, gammax, gammay, gammaz, gammath, gammaph;

	TVector3 vtx, dircos, q3mom;
	TVector3 cmBoost, labBoost, pi0Boost, etaBoost;

	TLorentzVector k, p, p1, q, p4In, p4Out, forpion, foreta;
	TLorentzVector gamma1, gamma2, photon1, photon2;
	TLorentzVector k_cm, p_cm, p1_cm, q_cm;

	TLorentzVector electron, positron;

	TString name, string, gnames;

	// Adding a comment here for no reason

	std::cout << "--------" << std::endl;
	std::cout << "EvGen26" << std::endl;
	std::cout << "--------" << std::endl;

	// Read in target and beam parameters
	flag = ReadParams();

	//	Read in fit values for Compton angular distributions
	ReadABC( fParam.process);

	update = fParam.counts/20;

	// For the case of Compton scattering
	// 	Particle 1:		Target
	// 	Particle 2:		Scattered Photon
	// For the case of neutral pion production
	// 	Particle 1:		Target
	// 	Particle 2:		Pion (not tracked in G4)
	// 	Particle 3:		Decay photon 1
	// 	Particle 4:		Decay photon 2
	// For the case of eta production
	// 	Particle 1:		Target
	// 	Particle 2:		Eta (not tracked in G4)
	// 	Particle 3:		Decay photon 1
	// 	Particle 4:		Decay photon 2

	// Particle ID #s etc.
	// Default target is proton
	ptag[0] = 14;
	pm = kMP_MEV/1000;								// Default is proton
	// 12C nucleus as target or for cell windows for LH2
	if ( ( fParam.tgt == "c") || ( fParam.tgt == "w")) {
		ptag[0] = 67;
		pm = kM_C12_MEV/1000;
	}
	else if (fParam.tgt == "he3") {
		ptag[0] = 49;
		pm = kM_HE3_MEV/1000;
	}
	else if (fParam.tgt == "he4") {
		ptag[0] = 47;
		pm = kM_HE4_MEV/1000;
	}
	else if (fParam.tgt == "li6") {
		ptag[0] = 61;
		pm = kM_LI6_MEV/1000;
	}
	else if (fParam.tgt == "li7") {
		ptag[0] = 62;
		pm = kM_LI7_MEV/1000;
	}

	// Default scattered particle is one photon
	npart = 2;
	ptag[1] = 1;
	qm = 0;

	// One pion decaying into two photons
	if ( fParam.process == "pi0") {
		npart = 4;
		ptag[1] = 7;					// pi0
		qm = kMPI0_MEV/1000;
		ptag[2] = 1;					// gamma
		ptag[3] = 1;					// gamma
	}

	// One eta decaying into two photons
	if ( fParam.process == "eta") {
		npart = 4;
		ptag[1] = 17;					// eta
		qm = kMETA_MEV/1000;
		ptag[2] = 1;					// gamma
		ptag[3] = 1;					// gamma
	}

	// Array for filling ntuple
	Float_t var[100];

	//	Minimum photon energy required for reaction (in GeV)
	//	kthr for eta production is 0.70723 GeV
	kthr = qm + qm*qm/(2*pm);

	// Set the seed for the random number generator
	gRandom->SetSeed();

	//	Generate GEANT name string for ntuple
	gnames = GenNames( npart, ptag);

	// Spit out update message 20 times
	update = fParam.counts/20;

	if ( flag == kTRUE) {
		j = 0;
		stop = fParam.e_bins;
	}
	else {
		j = fParam.chan_lo;
		stop = fParam.chan_hi+1;
	}

	// Energy range for incident photons
	while ( j < stop) {

		if ( flag == kTRUE) {
			e_mid = fParam.e_lo + j*fParam.e_step;
			e_bite = fParam.e_bite/2;
			name = Form( "out/%dcm/%s_%s_%d_in.root", (int) fParam.tgt_len,
				(const char*) fParam.process, (const char*) fParam.tgt, (int) e_mid);
		}
		else {
			name = Form( "Channel %d", j);
			std::cout << name << std::endl;
			name = Form( "out/%dcm/%s_%s_chan%d_in.root", (int) fParam.tgt_len,
				(const char*) fParam.process, (const char*) fParam.tgt, j);
			e_mid = fEnergy[j];
			e_bite = fDEnergy[j];
		}

		TFile hfile( name, "RECREATE", "MC_Ntuple_File");

		//check open hfile
		if(hfile.IsZombie()){
			std:: cout << "Error: Couldn't find output file: " << name << std::endl;
			exit(1);
		}

		e_low = e_mid - e_bite;
		e_high = e_mid + e_bite;
		name = Form( "Energy = %5.1f +/- %4.1f MeV", e_mid, e_bite);
		std::cout << name << std::endl;

		if ( e_high < (kthr*1000)) { // converts kthr to MeV
			std::cout << "Photon Energy below threshold.  Skipping..."
				<< std::endl;
			continue;
		}

		e_low /= 1000;
		e_high /= 1000;

		//	Create ntuple for kinematic variables
		//	It is absolutely necessary for Geant4!
		TNtuple *h1 = new TNtuple( "h1", "TMCUserGenerator", gnames);

		// These histograms are only for debugging.  You can comment them out if
		// you want.  However, if you are trying to calculate angular
		// efficiencies, then you will need to keep "h4", the scattered lab
		// angle.
				
		TH1F *h2;
		TH1F *h3;
		TH1F *h4;
		TH1F *h5;
		TH1F *h6;
		TH1F *h7;
		TH1F *h8;
		TH1F *h9;
		TH1F *h10;
		TH1F *h11;
		TH1F *h12;
		TH1F *h13;
		TH1F *h14;
		TH1F *h15;
		TH1F *h16;
		TH2F *h17;
		
		if ( fParam.process != "eta") {
			h2 = new TH1F( "h2", "Photon Beam Energy (MeV)", 450, 0, 450);
			h3 = new TH1F( "h3", "Scattered KE (MeV)", 300, 0, 300);
			h4 = new TH1F( "h4", "Scattered #theta (deg)", 360, 0, 180);
			h5 = new TH1F( "h5", "Scattered CM #theta", 360, 0, 180);
			h6 = new TH1F( "h6", "Scattered #phi (deg)", 720, -180, 180);
			h7 = new TH1F( "h7", "Recoil KE (MeV)", 300, 0, 300);
			h8 = new TH1F( "h8", "Recoil #theta (deg)", 360, 0, 180);
			h9 = new TH1F( "h9", "Recoil #phi (deg)", 720, -180, 180);

			// These are only for pi0 production but must be defined regardless... now useful for eta production as well!
			h10 = new TH1F( "h10", "Decay Particle 1 KE (MeV)", 300, 0, 300);
			h11 = new TH1F( "h11", "Decay Particle 1 #theta (deg)", 36, 0, 180);
			h12 = new TH1F( "h12", "Decay Particle 1 #phi (deg)", 72, -180,
				180);
			h13 = new TH1F( "h13", "Decay Particle 2 KE (MeV)", 300, 0, 300);
			h14 = new TH1F( "h14", "Decay Particle 2 #theta (deg)", 36, 0, 180);
			h15 = new TH1F( "h15", "Decay Particle 2 #phi (deg)", 72, -180,
				180);

			h16 = new TH1F( "h16", "Scattered CM std::cos(#theta)", 20, -1, 1);
		
			h17 = new TH2F( "h17", "Recoil KE vs Theta", 36, 0, 180, 300, 0, 300);		
		}
		
		if ( fParam.process == "eta") {
			h2 = new TH1F( "h2", "Photon Beam Energy (MeV)", 480, 700, 3100);
			h3 = new TH1F( "h3", "Scattered KE (MeV)", 500, 0, 2500);
			h4 = new TH1F( "h4", "Scattered #theta (deg)", 360, 0, 180);
			h5 = new TH1F( "h5", "Scattered CM #theta", 360, 0, 180);
			h6 = new TH1F( "h6", "Scattered #phi (deg)", 720, -180, 180);
			h7 = new TH1F( "h7", "Recoil KE (MeV)", 500, 0, 2500);
			h8 = new TH1F( "h8", "Recoil #theta (deg)", 360, 0, 180);
			h9 = new TH1F( "h9", "Recoil #phi (deg)", 720, -180, 180);

			// These are only for pi0 production but must be defined regardless... now useful for eta production as well!
			h10 = new TH1F( "h10", "Decay Particle 1 KE (MeV)", 500, 0, 3000);
			h11 = new TH1F( "h11", "Decay Particle 1 #theta (deg)", 36, 0, 180);
			h12 = new TH1F( "h12", "Decay Particle 1 #phi (deg)", 72, -180,
				180);
			h13 = new TH1F( "h13", "Decay Particle 2 KE (MeV)", 500, 0, 3000);
			h14 = new TH1F( "h14", "Decay Particle 2 #theta (deg)", 36, 0, 180);
			h15 = new TH1F( "h15", "Decay Particle 2 #phi (deg)", 72, -180,
				180);

			h16 = new TH1F( "h16", "Scattered CM std::cos(#theta)", 20, -1, 1);
		
			h17 = new TH2F( "h17", "Recoil KE vs Theta", 36, 0, 180, 500, 0, 2500);		
		}

		//
		// New TTree stuff
		// 
		Double_t beamE;
		Double_t scatKE, scatThCM, scatTh, scatPhi;
		Double_t recoilKE, recoilThCM, recoilTh, recoilPhi;
		TTree *egTree = new TTree( "EvGen26Tree", "EvGen26 Compton, Pi0, and Eta Kinematics");
		egTree->Branch( "BeamE", &beamE);
		egTree->Branch( "ScatKE", &scatKE);
		egTree->Branch( "ScatTheta", &scatTh);
		egTree->Branch( "ScatThetaCM", &scatThCM);
		egTree->Branch( "ScatPhi", &scatPhi);
		egTree->Branch( "RecoilKE", &recoilKE);
		egTree->Branch( "RecoilTheta", &recoilTh);
		egTree->Branch( "RecoilThetaCM", &recoilThCM);
		egTree->Branch( "RecoilPhi", &recoilPhi);
		//
		//

		// Bremsstrahlung distribution for beam energy (in GeV)
		TF1 *f1 = new TF1( "f1", "1/x", e_low, e_high);

		// Angular distributions for proton, 3/4-He, 12-C, 16-O
		// ScatCTH distribution for Compton scattering and pi0 production
		TF1 *f2;
		f2 = new TF1( "ScatCTH",
			[this](Double_t* x, Double_t* par) { return this->ScatCTH(x, par); },
			-1, 1, 1);
		// It depends on incident photon energy (in MeV!)
		f2->SetParameter( 0, e_mid);		

		// Phi distributions using beam pol and asymmetry
		TF1 *f3 = new TF1( "PhiDist",
			[this](Double_t* x, Double_t* par) { return this->PhiDist(x, par); },
			-kPI, kPI, 3);
		// It depends on beam polarization and photon asymmetry
		f3->SetParameter( 0, fParam.p_gamma);
		f3->SetParameter( 1, fParam.Sigma);
		f3->SetParameter( 2, fParam.phi_offset);

		for ( i = 1; i <= fParam.counts; i++)
		{

			if ( i && (i%update) == 0)
				std::cout << "     events analysed: " << i << std::endl;

			// Incident photon energy
			// It obviously must be greater than pion threshold
			ke = f1->GetRandom();
			if ( ke < kthr) continue;

			// Reaction vertex

			// Vertex position for windows of LH2 target
			// This is valid for both 5cm and 10cm target lengths
			if ( fParam.tgt == "w") {
				Double_t w = gRandom->Rndm();
				Double_t thick;

				if ( ( 0 <= w) && ( w <= 0.2427))
				{ 
					// 125 um Kapton window
					thick = 0.0125;
					vtx_z = thick*(-0.5 + gRandom->Rndm()) - fParam.tgt_len/2;
				}
				else if ( ( 0.2427 < w) && ( w <= 0.4854))
				{ 
					// 125 um Kapton window
					thick = 0.0125;
					vtx_z = thick*(-0.5 + gRandom->Rndm()) + fParam.tgt_len/2;
				}
				else if ( ( 0.4854 < w) && (w <= 0.5242))
				{ 
					// 20 um Kapton window
					thick = 0.0020;
					vtx_z = thick*(-0.5 + gRandom->Rndm()) + fParam.tgt_len/2 + 5;
				}
				else if ( ( 0.5242 < w) && ( w <= 0.7669))
				{ 
					// 125 um Kapton window
					thick = 0.0125;
					vtx_z = thick*(-0.5 + gRandom->Rndm()) + fParam.tgt_len/2 + 7;
				
				}
				else
				{ 
					// 120 um PEN foil
					thick = 0.0120;
					vtx_z = thick*(-0.5 + gRandom->Rndm()) + fParam.tgt_len/2 + 7.1;
				}
			}
			// Non-target window vertex (everthing else)
			else vtx_z = fParam.tgt_len*(-0.5 + gRandom->Rndm());

			// The while statement cuts off the gaussian xy values of vertex
			// position so that they are inside the target.
			while ( std::sqrt( Sqr( vtx_x = gRandom->Gaus( 0, fParam.sig_bm)) + Sqr(
							vtx_y = gRandom->Gaus( 0, fParam.sig_bm))) <
					fParam.vtx_rad);

			vtx.SetXYZ( vtx_x, vtx_y, vtx_z);

			// Target is at rest
			p.SetPxPyPzE( 0, 0, 0, pm);

			// Incident beam direction using direction cosines.
			dircos.SetXYZ( 0, 0, 1);
			k.SetPxPyPzE( dircos.X()*ke, dircos.Y()*ke, dircos.Z()*ke, ke);

			// Total incoming 4-momentum
			p4In = k + p;

			// These 3-vectors boosts between the lab and CM frames
			labBoost = p4In.BoostVector();
			cmBoost = -p4In.BoostVector();

			// Boost the initial-state particles' 4-momentum to the CM frame
			k_cm = k;
			k_cm.Boost( cmBoost);
			p_cm = p;
			p_cm.Boost( cmBoost);

			// Scattered CM energy
			qe_cm = (2*pm*ke + qm*qm)/(2*std::sqrt( pm*pm+2*pm*ke)); 

			// CM momentum is the magnitude of the 3-momentum of either
			// CM particle
			mom_cm = std::sqrt( qe_cm*qe_cm - qm*qm);

			// Target CM energy
			pe_cm = std::sqrt( mom_cm*mom_cm + pm*pm);

			// Pick CM angular distributions for scattered particle:
			// 	Theta is from ABC fit function
			
			qth_cm = std::acos( f2->GetRandom());

			// 	Phi is from a function with p_gamma and Sigma
			qph_cm = f3->GetRandom();

			// Momentum components
			qx_cm = mom_cm*std::sin( qth_cm)*std::cos( qph_cm);
			qy_cm = mom_cm*std::sin( qth_cm)*std::sin( qph_cm);
			qz_cm = mom_cm*std::cos( qth_cm);

			q_cm.SetPxPyPzE( qx_cm, qy_cm, qz_cm, qe_cm);

			// CM momentum conservation
			p1_cm = k_cm + p_cm - q_cm;

			// Boost CM values of final-state particles back to lab frame
			q = q_cm;
			q.Boost( labBoost);
			p1 = p1_cm;
			p1.Boost( labBoost);

			// This next section is for the two decay particles from the pi0
			
			if ( fParam.process == "pi0") {

				// Find the boost to pion rest frame
				q3mom = q.Vect();
				forpion.TLorentzVector::SetPxPyPzE( 0, 0, q3mom.Mag(), q.E());
				pi0Boost = forpion.BoostVector();

				// Pion decay into 2 photons
				gammae = qm/2;

				// Pick random theta and phi from phase space for photon 1
				gammath = std::acos( -1 + 2*gRandom->Rndm());
				gammaph = kPI*( -1 + 2*gRandom->Rndm());

				// Calculate momentum components of photon 1
				gammax = gammae*std::sin( gammath)*std::cos( gammaph);
				gammay = gammae*std::sin( gammath)*std::sin( gammaph);
				gammaz = gammae*std::cos( gammath);

				// 4-momenta for the two gammas (#2 is just in the opposite
				// direction but with same energy in the pion rest frame)
				gamma1.SetPxPyPzE( gammax, gammay, gammaz, gammae);
				gamma2.SetPxPyPzE( -gammax, -gammay, -gammaz, gammae);

				// Boost photons to frame of pion
				photon1 = gamma1;
				photon1.Boost( pi0Boost);
				photon2 = gamma2;
				photon2.Boost( pi0Boost);

				// Rotate the angles to go from the direction of pion back to lab
				// frame
				photon1.RotateY( q.Theta());
				photon1.RotateZ( q.Phi());
				photon2.RotateY( q.Theta());
				photon2.RotateZ( q.Phi());
			}
			
			// This section is for the two decay particles from the eta
			
			if ( fParam.process == "eta") {
				
				// Find the boost to eta rest frame
				q3mom = q.Vect();
				foreta.TLorentzVector::SetPxPyPzE( 0, 0, q3mom.Mag(), q.E());
				etaBoost = foreta.BoostVector();
				
				// Eta decay into 2 photons
				gammae = qm/2;
				
				// Pick random theta and phi from phase space for photon 1
				gammath = std::cos( -1 + 2*gRandom->Rndm());
				gammaph = kPI*( -1 + 2*gRandom->Rndm());
				
				// Calculate momentum components of photon 1
				gammax = gammae*std::sin( gammath)*std::cos( gammaph);
				gammay = gammae*std::sin( gammath)*std::sin( gammaph);
				gammaz = gammae*std::cos( gammath);
				
				// 4-momenta for the two photons (#2 is just in the opposite
				// direction but with same energy in the eta rest frame)
				gamma1.SetPxPyPzE( gammax, gammay, gammaz, gammae);
				gamma2.SetPxPyPzE( -gammax, -gammay, -gammaz, gammae);
				
				// Boost photons to frame of eta
				photon1 = gamma1;
				photon1.Boost( etaBoost);
				photon2 = gamma2;
				photon2.Boost( etaBoost);
				
				// Rotate the angles to go from the direction of eta back to lab
				// frame
				photon1.RotateY( q.Theta());
				photon1.RotateZ( q.Phi());
				photon2.RotateY( q.Theta());
				photon2.RotateZ( q.Phi());
			}

			// These next lines store the particle properties in the ntuple
			// variable.

			// Interaction vertex position
			var[0] = vtx.X();
			var[1] = vtx.Y();
			var[2] = vtx.Z();

			// Incident photon beam
			var[3] = dircos.X();
			var[4] = dircos.Y();
			var[5] = dircos.Z();
			var[6] = ke;
			var[7] = ke;

			// Recoil nucleus
			var[8] = p1.Px()/p1.P();
			var[9] = p1.Py()/p1.P();
			var[10] = p1.Pz()/p1.P();
			var[11] = p1.P();
			var[12] = p1.E();

			// Scattered Particle
			var[13] = q.Px()/q.P();
			var[14] = q.Py()/q.P();
			var[15] = q.Pz()/q.P();
			var[16] = q.P();
			var[17] = q.E();

			if ( fParam.process == "pi0") {

				// Photon 1
				var[18] = photon1.Px()/photon1.P();
				var[19] = photon1.Py()/photon1.P();
				var[20] = photon1.Pz()/photon1.P();
				var[21] = photon1.P();
				var[22] = photon1.E();

				// Photon 2
				var[23] = photon2.Px()/photon2.P();
				var[24] = photon2.Py()/photon2.P();
				var[25] = photon2.Pz()/photon2.P();
				var[26] = photon2.P();
				var[27] = photon2.E();
			}
			
			// Does the same thing as above but for eta process
			if ( fParam.process == "eta") {

				// Photon 1
				var[18] = photon1.Px()/photon1.P();
				var[19] = photon1.Py()/photon1.P();
				var[20] = photon1.Pz()/photon1.P();
				var[21] = photon1.P();
				var[22] = photon1.E();

				// Photon 2
				var[23] = photon2.Px()/photon2.P();
				var[24] = photon2.Py()/photon2.P();
				var[25] = photon2.Pz()/photon2.P();
				var[26] = photon2.P();
				var[27] = photon2.E();
			}

			// Fill ntuple
			h1->Fill( var);

			// Fill histograms with quantities (in MeV)
			h2->Fill( ke*1000);
			h3->Fill( 1000*(q.E() - q.M()));
			h4->Fill( q.Theta()/kD2R);
			h5->Fill( q_cm.Theta()/kD2R);
			h16->Fill( std::cos( q_cm.Theta()));
			h6->Fill( q.Phi()/kD2R);
			h7->Fill( 1000*(p1.E() - p1.M()));
			h8->Fill( p1.Theta()/kD2R);
			h9->Fill( p1.Phi()/kD2R);
			h17->Fill( p1.Theta()/kD2R, 1000*(p1.E()-p1.M()));

			if ( fParam.process == "pi0") {

				h10->Fill( 1000*photon1.E());
				h11->Fill( photon1.Theta()/kD2R);
				h12->Fill( photon1.Phi()/kD2R);
				h13->Fill( 1000*photon2.E());
				h14->Fill( photon2.Theta()/kD2R);
				h15->Fill( photon2.Phi()/kD2R);

			}
			
			// Does the same thing as above but for eta process
			if ( fParam.process == "eta") {

				h10->Fill( 1000*photon1.E());
				h11->Fill( photon1.Theta()/kD2R);
				h12->Fill( photon1.Phi()/kD2R);
				h13->Fill( 1000*photon2.E());
				h14->Fill( photon2.Theta()/kD2R);
				h15->Fill( photon2.Phi()/kD2R);

			}

			// Tree variables
			beamE = ke*1000;
			scatKE = 1000*(q.E() - q.M());
			scatTh = q.Theta()/kD2R;
			scatThCM = q_cm.Theta()/kD2R;
			scatPhi = q.Phi()/kD2R;
			recoilKE = 1000*(p1.E() - p1.M());
			recoilTh = p1.Theta()/kD2R;
			recoilThCM = 180 - scatThCM/kD2R;
			recoilPhi = p1.Phi()/kD2R;

			// Fill Tree
			egTree->Fill();
		}

		// Write tree to file
		hfile.cd();
		egTree->Write();

		// Write histograms to file
		hfile.Write();

		// This isn't really necessary, but can be used for debugging.
//		h1->Print();

		// Close file
		hfile.Close();

		j++;
	}

	return 0;

}

int main()
{
	EvGen26 evgen;
	return evgen.Run();
}
