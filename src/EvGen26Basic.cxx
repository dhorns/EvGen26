//
// EvGen26Basic.cxx
//
// Object-oriented version of EvGenBasic.cxx
// DLH original with help from ChatGPT
//

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
#include "Riostream.h"

#include "physics.h"

class EvGen26Basic
{

	public:

		EvGen26Basic( const char* parFile = "par/EvGen26Basic.in");
		~EvGen26Basic() = default;

		Int_t Run();

	private:

		void SetDefaults();
		void ReadParameters( const char* parFile);
		void ValidateParameters() const;
		void PrintParameters() const;
		void SetParticleInfo();
		void GenerateEnergyBin( UInt_t bin);

		TString GenNames() const;
		TString OutputFileName( Double_t eMid) const;

		static Double_t Sqr(Double_t x) { return x*x; }
		static Double_t Momentum(Double_t e, Double_t m) { return std::sqrt(e*e - m*m); }

	private:

		// Parameters read from file, with defaults
		TString  fParticle;
		Double_t fTargetLength;
		Double_t fVertexRadius;
		Int_t    fCounts;

		Double_t fEnergyBite;
		Double_t fEnergyLow;
		Double_t fEnergyStep;
		UInt_t   fEnergyBins;

		Double_t fTheta;
		Double_t fThetaBite;

		Double_t fPhi;
		Double_t fPhiBite;

		// Fixed/internal generation settings
		Double_t fBeamSigma;
		UInt_t   fNPart;
		UInt_t   fPTag[10];
		Double_t fParticleMass;

};

int main()
{

	EvGen26Basic evgen;
	return evgen.Run();

}

EvGen26Basic::EvGen26Basic( const char* parFile)
{

	SetDefaults();
	ReadParameters( parFile);
	ValidateParameters();
	SetParticleInfo();
	PrintParameters();

}

void EvGen26Basic::SetDefaults()
{

	fParticle = "g";
	fTargetLength = 2.00;
	fVertexRadius = 0.50;
	fCounts = static_cast<Int_t>(1e5);
	
	fEnergyBite = 10.0;
	fEnergyLow = 150.0;
	fEnergyStep = 10.0;
	fEnergyBins = 1;
	
	fThetaBite = 180.0;
	fTheta = 90.0;
	
	fPhiBite = 360.0;
	fPhi = 180.0;
	
	fBeamSigma = 0.5;
	
	fNPart = 1;
	for ( UInt_t i = 0; i < 10; ++i) fPTag[i] = 0;
	fPTag[0] = 1;
	fParticleMass = 0.0;

}

void EvGen26Basic::ReadParameters( const char* parFile)
{

	TString name;
	TString string;

	std::ifstream inFile( parFile);
	if ( !inFile.is_open())
	{
		std::cout << "Error opening file " << parFile << std::endl;
		exit(-1);
	}

	while ( !inFile.eof())
	{
		name.ReadLine( inFile);

		if ( name.Length() == 0) continue;
		if ( name[0] == '#') continue;

		string = "Particle: ";
		if ( name.Contains(string)) fParticle = name.Remove( 0, string.Length());

		string = "TargetLength: ";
		if ( name.Contains(string))
		{
			name.Remove( 0, string.Length());
			fTargetLength = name.Atof();
		}

		string = "BeamSpotRadius: ";
		if ( name.Contains( string))
		{
			name.Remove( 0, string.Length());
			fVertexRadius = name.Atof();
		}

		string = "Throws: ";
		if ( name.Contains( string)) {
			name.Remove( 0, string.Length());
			fCounts = name.Atoi();
		}

		string = "ParticleEnergy: ";
		if ( name.Contains( string))
		{
			name.Remove( 0, string.Length());

			string = name;
			string.Remove( 0, string.Last(' '));
			fEnergyBins = string.Atoi();

			string = name;
			string.Remove( string.First(' '));
			fEnergyBite = string.Atof();

			string = name;
			string.Remove( string.Last(' '));
			string.Remove( 0, string.Last(' '));
			fEnergyStep = string.Atof();

			string = name;
			string.Remove( string.Last(' '));
			string.Remove( string.Last(' '));
			string.Remove( 0, string.First(' '));
			fEnergyLow = string.Atof();
		}

		string = "LabTheta: ";
		if ( name.Contains( string))
		{
			name.Remove( 0, string.Length());

			string = name;
			string.Remove( string.Last(' '));
			fTheta = string.Atof();

			string = name;
			string.Remove( 0, string.First(' '));
			fThetaBite = string.Atof();
		}

		string = "LabPhi: ";
		if ( name.Contains( string))
		{
			name.Remove( 0, string.Length());

			string = name;
			string.Remove( string.Last(' '));
			fPhi = string.Atof();

			string = name;
			string.Remove( 0, string.First(' '));
			fPhiBite = string.Atof();
		}
	}

	inFile.close();
}

void EvGen26Basic::ValidateParameters() const
{
	if ((fParticle != "g") && (fParticle != "p") && (fParticle != "d") &&
		(fParticle != "he3") && (fParticle != "he4"))
	{
		std::cout << "Particle \"" << fParticle << "\" not supported" << std::endl;
		exit( -1);
	}

	if ( fCounts <= 0)
	{
		std::cout << "Number of throws must be positive" << std::endl;
		exit( -1);
	}

	if ( fEnergyBins == 0)
	{
		std::cout << "Number of energy bins must be positive" << std::endl;
		exit( -1);
	}
}

void EvGen26Basic::SetParticleInfo()
{
	// Particle ID #s are standard GEANT numbers.
	// Default particle is a photon.
	fNPart = 1;
	fPTag[0] = 1;
	fParticleMass = 0.0;

	if ( fParticle == "p")
	{
		fPTag[0] = 14;
		fParticleMass = kMP_MEV;
	}
	else if ( fParticle == "d")
	{
		fPTag[0] = 45;
		fParticleMass = kM_DEUT_MEV;
	}
	else if ( fParticle == "he3")
	{
		fPTag[0] = 49;
		fParticleMass = kM_HE3_MEV;
	}
	else if ( fParticle == "he4")
	{
		fPTag[0] = 47;
		fParticleMass = kM_HE4_MEV;
	}

	// To GeV
	fParticleMass /= 1000.0;

}

void EvGen26Basic::PrintParameters() const
{

	TString name;

	std::cout << "--------" << std::endl;
	std::cout << "EvGen26Basic" << std::endl;
	std::cout << "--------" << std::endl;
	
	name = "Particle = " + fParticle;
	std::cout << name << std::endl;
	
	name = Form( "Target Length = %4.2f cm", fTargetLength);
	std::cout << name << std::endl;
	
	name = Form( "Beam Spot Radius = %4.2f cm", fVertexRadius);
	std::cout << name << std::endl;
	
	name = Form( "Throws = %d", fCounts);
	std::cout << name << std::endl;
	
	name = Form( "Particle Energy: %5.1f +/- %5.1f MeV in steps of %5.1f MeV for %2d bins",
			fEnergyLow, fEnergyBite, fEnergyStep, fEnergyBins);
	std::cout << name << std::endl;
	
	name = Form( "Particle Theta = %5.1f +/- %5.1f deg", fTheta, fThetaBite/2.0);
	std::cout << name << std::endl;
	
	name = Form( "Particle Phi = %6.1f +/- %5.1f deg", fPhi, fPhiBite/2.0);
	std::cout << name << std::endl;

}

Int_t EvGen26Basic::Run()
{

	gRandom->SetSeed();

	for ( UInt_t j = 0; j < fEnergyBins; ++j) GenerateEnergyBin( j);

	return 0;
}

void EvGen26Basic::GenerateEnergyBin( UInt_t bin)
{

	Double_t eMid = fEnergyLow + bin*fEnergyStep;
	TString name = Form( "Energy = %5.1f +/- %4.1f MeV", eMid, fEnergyBite/2.0);
	std::cout << name << std::endl;
	
	name = OutputFileName( eMid);
	TFile hfile( name, "RECREATE", "MC_Ntuple_File");
	
	TString gnames = GenNames();
	
	// Create ntuple for kinematic variables.
	// This is absolutely necessary for Geant4.
	TNtuple* h1 = new TNtuple( "h1", "TMCUserGenerator", gnames);
	
	// These histograms are only for debugging.
	TH1F* h2 = new TH1F( "h2", "Particle KE (MeV)", 1000, 0, 1000);
	TH1F* h3 = new TH1F( "h3", "Particle Momentum (MeV/c)", 2000, 0, 2000);
	TH1F* h4 = new TH1F( "h4", "Particle #theta (deg)", 180, 0, 180);
	TH1F* h5 = new TH1F( "h5", "Particle #phi (deg)", 400, -200, 200);
	
	Float_t var[100];
	TVector3 vtx;
	TVector3 dircos;
	
	Int_t update = fCounts/20;
	if (update <= 0) update = 1;
	
	Double_t cth_hi = std::cos( ( fTheta - fThetaBite/2.0)*kD2R);
	Double_t cth_lo = std::cos( ( fTheta + fThetaBite/2.0)*kD2R);

	for ( Int_t i = 1; i <= fCounts; ++i)
	{
		if ( ( i % update) == 0)
 			std::cout << "     events analysed: " << i << std::endl;

		// Particle energy.
		Double_t ke = gRandom->Uniform( eMid - fEnergyBite/2.0, eMid + fEnergyBite/2.0);
		ke /= 1000.0;

		// Particle theta and phi.
		Double_t th = std::acos( gRandom->Uniform(cth_lo, cth_hi));
//		Double_t th = gRandom->Uniform( fTheta - fThetaBite/2.0, fTheta + fThetaBite/2.0)*kD2R;
		Double_t ph = gRandom->Uniform( fPhi - fPhiBite/2.0, fPhi + fPhiBite/2.0)*kD2R;

		Double_t vtx_x = 0.0;
		Double_t vtx_y = 0.0;
		Double_t vtx_z = fTargetLength*( -0.5 + gRandom->Rndm());

		while (std::sqrt(Sqr(vtx_x = gRandom->Gaus(0.0, fBeamSigma)) +
			Sqr(vtx_y = gRandom->Gaus(0.0, fBeamSigma))) < fVertexRadius);

		vtx.SetXYZ( vtx_x, vtx_y, vtx_z);

		// Interaction vertex position.
		var[0] = vtx.X();
		var[1] = vtx.Y();
		var[2] = vtx.Z();

		// Incident photon beam, in this case turned off.
		var[3] = 0;
		var[4] = 0;
		var[5] = 0;
		var[6] = 0;
		var[7] = 0;

		// Direction cosines.
		dircos.SetXYZ( std::sin( th)*std::cos( ph), std::sin( th)*std::sin( ph), std::cos( th));

		// Energy and momentum.
		Double_t energy = ke + fParticleMass;
		Double_t momentum = Momentum( energy, fParticleMass);

		// Particle variables.
		var[8] = dircos.X();
		var[9] = dircos.Y();
		var[10] = dircos.Z();
		var[11] = momentum;
		var[12] = energy;

		h1->Fill(var);

		// Fill histograms with quantities in MeV.
		h2->Fill( ke*1000.0);
		h3->Fill( momentum*1000.0);
		h4->Fill( th/kD2R);
		h5->Fill( ph/kD2R);

	}

	hfile.Write();
//	h1->Print();
	hfile.Close();
}

TString EvGen26Basic::GenNames() const
{

	TString pstr[] = {"Px", "Py", "Pz", "Pt", "En"};
	TString beam = "X_vtx:Y_vtx:Z_vtx:Px_bm:Py_bm:Pz_bm:Pt_bm:En_bm";
	TString particles;
	TString namesgit remote add origin git@github.com:dhorns/EvGen26.git;

	for ( UInt_t i = 0; i < fNPart; ++i)
	{
		for ( UInt_t j = 0; j < 5; ++j)
		{
			particles.Append( pstr[j]);
			if ( ( i == (fNPart - 1)) && (j == 4))
				particles.Append( Form( "_l%02d%02d", i+1, fPTag[i]));
			else
				particles.Append( Form( "_l%02d%02d:", i+1, fPTag[i]));
		}
	}

	names = beam + ":" + particles;
	return names;

}

TString EvGen26Basic::OutputFileName( Double_t eMid) const
{

	return Form( "out/basic/%s_%d_%d_%d_in.root",
		static_cast<const char*>( fParticle),
		static_cast<Int_t>( eMid),
		static_cast<Int_t>( fTheta),
		static_cast<Int_t>( fPhi));

}
