#include <TGeoManager.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TGeoNode.h>
#include <TMaterial.h>
#include <TMixture.h>
#include <TShape.h>
#include <TString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <stdlib.h>
#include <TROOT.h>
#include <TGenPhaseSpace.h>
#include <TStyle.h>
#include "TGeoPatternFinder.h"
#include "TSystem.h"
#include <TGeoNavigator.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include "TApplication.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TColor.h"
#include <TRandom3.h>
#include <time.h>
#include "TTree.h"
#include <random>
#include <chrono>

using namespace std;

double funcl(double x1, double x2, double x3, double x4, double x5){
	//           zp         zt         Elab       thetalab_rad   mp/mt
	// units mb/sr
	double ruthlab=1.296*pow(x1*x2/x3,2)*(pow(1/sin(x4/2),4)-2.*x5*x5);
	return(ruthlab);
}

TVector3 ToBeamFrame(const TVector3& v,
																					const TVector3& ex,
																					const TVector3& ey,
																					const TVector3& ez)
{
	return TVector3(v.Dot(ex), v.Dot(ey), v.Dot(ez));
}

// Struttura di supporto per restituire tutti i dati utili del colpo calcolato
struct HitData {
	double Ta;
	double Pa;
	double cos_TA;
	int iii;        // Indice della strip frontale
	int jjj;        // Indice della strip posteriore
	TVector3 vbina; // Direzione apparente della particella
};

HitData ProcessDetectorHit(
																											TVector3 vhp, TVector3 vCenter, TGeoRotation* rot,
																											int numStrips, double step, TRandom3* ran,
																											TVector3 beamOrigin, TVector3 ex_beam, TVector3 ey_beam, TVector3 ez_beam,
																											double normalY, double directionSignX)
{
	HitData result;
	double deg = 180. / TMath::Pi();
	
	// 1. Trasformazione nel sistema locale
	TVector3 vd = vhp - vCenter;
	double vvd[3] = {vd(0), vd(1), vd(2)};
	double vvdr[3] = {0., 0., 0.};
	rot->MasterToLocal(vvd, vvdr);
	
	// 2. Ricerca del pixel (front e back)
	int ii = 0, jj = 0;
	for (ii = 0; ii < numStrips; ii++) {
		double limit1 = directionSignX < 0 ? (-2.48 + ii * step) : (2.48 - ii * step);
		double limit2 = directionSignX < 0 ? (-2.48 + (ii + 1) * step) : (2.48 - (ii + 1) * step);
		
		if (directionSignX < 0) {
			if (vvdr[0] >= limit1 && vvdr[0] < limit2) break;
		} else {
			if (vvdr[0] < limit1 && vvdr[0] >= limit2) break;
		}
	}
	for (jj = 0; jj < numStrips; jj++) {
		if (vvdr[2] >= -2.48 + jj * step && vvdr[2] < -2.48 + (jj + 1) * step) break;
	}
	
	// 3. Randomizzazione all'interno del pixel
	double cf = directionSignX < 0 ? (-2.48 + ii * step + step * ran->Rndm()) : (2.48 - ii * step - step * ran->Rndm());
	double cb = -2.48 + jj * step + step * ran->Rndm();
	
	// Salviamo gli indici delle strip per gli istogrammi
	result.iii = ii;
	result.jjj = jj;
	
	// 4. Ritorno alle coordinate globali (apparenti)
	double avdr[3] = {cf, vvdr[1], cb};
	double avd[3] = {0., 0., 0.};
	double avdN[3] = {0., normalY, 0.};
	double avN[3] = {0., 0., 0.};
	
	rot->LocalToMaster(avdr, avd);
	rot->LocalToMaster(avdN, avN);
	
	// 5. Calcolo cinematica finale
	TVector3 vhit = vCenter + TVector3(avd[0], avd[1], avd[2]);
	TVector3 vbin_lab = vhit - beamOrigin;
	TVector3 vbin = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
	
	result.Ta = vbin.Theta() * deg;
	result.Pa = vbin.Phi() * deg;
	result.vbina = vbin; // Salviamo la direzione per lo smearing
	
	TVector3 ubin = vbin.Unit();
	TVector3 ubinN = ToBeamFrame(TVector3(avN[0], avN[1], avN[2]), ex_beam, ey_beam, ez_beam).Unit();
	result.cos_TA = ubinN.Dot(ubin);
	
	return result;
}

void geometry_and_physics_aau_15() {
	
	// Inizio cronometro per il Setup Generale
	auto start_total = std::chrono::high_resolution_clock::now();
	
	// TApplication theApp("App", NULL, NULL);
	
	gStyle->SetPalette(1);
	
	gSystem->Load("libPhysics.so");
	
	gSystem->Load("libGeom.so");
	
	TGeoManager *sc = new TGeoManager("scattcham","Scattering Chamber");
	
	TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum",0,0,0);
	TGeoMaterial *matSi = new TGeoMaterial("Si",28.086,14,2.321);
	
	TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1,matVacuum);
	TGeoMedium *Si = new TGeoMedium("Root Material",2,matSi);
	
	TGeoVolume *top = sc->MakeSphere("TOP",Vacuum,0.,100.,0.,180.,0.,360.);
	sc->SetTopVolume(top);
	
	TGeoVolume *psdA = sc->MakeBox("PSDA",Si,2.48,0.05,2.48);
	TGeoVolume *psdB = sc->MakeBox("PSDB",Si,2.48,0.05,2.48);
	
	// from deg to rad
	double dtr=TMath::Pi()/180.;
	
	// close geometry
	
	// this is detector A1
	double theta1=57.20951105*dtr;
	double   phi1=(45+180)*dtr;
	double  dist1=5.917993748;
	double    px1=dist1*sin(theta1)*cos(phi1);
	double    py1=dist1*sin(theta1)*sin(phi1);
	double    pz1=dist1*cos(theta1);
	TVector3 vA1(px1,py1,pz1); //position of det. A1 wrt chamber center
	TGeoRotation *rot1 = new TGeoRotation("rot1",-45,-34,0);
	TGeoCombiTrans *pos1 = new TGeoCombiTrans(px1,py1,pz1,rot1);
	
	// this is detector A2
	double theta2=34.06306845*dtr;
	double   phi2=(45+270-4.122298)*dtr;
	double  dist2=5.960733486;
	double    px2=dist2*sin(theta2)*cos(phi2);
	double    py2=dist2*sin(theta2)*sin(phi2);
	double    pz2=dist2*cos(theta2);
	TVector3 vA2(px2,py2,pz2); //position of det. A2 wrt chamber center
	TGeoRotation *rot2 = new TGeoRotation("rot2",45,-56,0);
	TGeoCombiTrans *pos2 = new TGeoCombiTrans(px2,py2,pz2,rot2);
	
	//this is detector A3
	double theta3=54.19850734*dtr;
	double   phi3=(45-7.546523)*dtr;
	double  dist3=5.820723714;
	double    px3=dist3*sin(theta3)*cos(phi3);
	double    py3=dist3*sin(theta3)*sin(phi3);
	double    pz3=dist3*cos(theta3);
	TVector3 vA3(px3,py3,pz3); //position of det. A3 wrt chamber center
	TGeoRotation *rot3 = new TGeoRotation("rot3",-45,34,0);
	TGeoCombiTrans *pos3 = new TGeoCombiTrans(px3,py3,pz3,rot3);
	
	// this is detector A4
	double theta4=33.48644304*dtr;
	double   phi4=(45+90+5.797854)*dtr;
	double  dist4=5.920746921;
	double    px4=dist4*sin(theta4)*cos(phi4);
	double    py4=dist4*sin(theta4)*sin(phi4);
	double    pz4=dist4*cos(theta4);
	TVector3 vA4(px4,py4,pz4); //position of det. A4 wrt chamber center
	TGeoRotation *rot4 = new TGeoRotation("rot4",45,56,0);
	TGeoCombiTrans *pos4 = new TGeoCombiTrans(px4,py4,pz4,rot4);
	
	// far geometry
	
	// this is detector B1
	double theta5=5.178370225*dtr;
	double   phi5=(45+180-1.140064)*dtr;
	double  dist5=66.82273415;
	double    px5=dist5*sin(theta5)*cos(phi5);
	double    py5=dist5*sin(theta5)*sin(phi5);
	double    pz5=dist5*cos(theta5);
	TVector3 vB1(px5,py5,pz5); //position of det. B1 wrt chamber center
	TGeoRotation *rot5 = new TGeoRotation("rot5",-45,90,0);
	TGeoCombiTrans *pos5 = new TGeoCombiTrans(px5,py5,pz5,rot5);
	
	// this is detector B2
	double theta6=4.36381047*dtr;
	double   phi6=(45+270+4.001077)*dtr;
	double  dist6=68.74930347;
	double    px6=dist6*sin(theta6)*cos(phi6);
	double    py6=dist6*sin(theta6)*sin(phi6);
	double    pz6=dist6*cos(theta6);
	TVector3 vB2(px6,py6,pz6); //position of det. B2 wrt chamber center
	TGeoRotation *rot6 = new TGeoRotation("rot6",45,90,0);
	TGeoCombiTrans *pos6 = new TGeoCombiTrans(px6,py6,pz6,rot6);
	
	// this is detector B3
	double theta7=5.306408794*dtr;
	double   phi7=(45-1.1124)*dtr;
	double  dist7=66.83643692;
	double    px7=dist7*sin(theta7)*cos(phi7);
	double    py7=dist7*sin(theta7)*sin(phi7);
	double    pz7=dist7*cos(theta7);
	TVector3 vB3(px7,py7,pz7); //position of det. B3 wrt chamber center
	TGeoRotation *rot7 = new TGeoRotation("rot7",-45,-90,0);
	TGeoCombiTrans *pos7 = new TGeoCombiTrans(px7,py7,pz7,rot7);
	
	// this is detector B4
	double theta8=3.955043244*dtr;
	double   phi8=(45+90-2.479051)*dtr;
	double  dist8=68.71364311;
	double    px8=dist8*sin(theta8)*cos(phi8);
	double    py8=dist8*sin(theta8)*sin(phi8);
	double    pz8=dist8*cos(theta8);
	TVector3 vB4(px8,py8,pz8); //position of det. B4 wrt chamber center
	TGeoRotation *rot8 = new TGeoRotation("rot8",45,-90,0);
	TGeoCombiTrans *pos8 = new TGeoCombiTrans(px8,py8,pz8,rot8);
	
	// detector placing in space
	
	top->AddNode(psdA,1,pos1);
	top->AddNode(psdA,2,pos2);
	top->AddNode(psdA,3,pos3);
	top->AddNode(psdA,4,pos4);
	
	top->AddNode(psdB,5,pos5);
	top->AddNode(psdB,6,pos6);
	top->AddNode(psdB,7,pos7);
	top->AddNode(psdB,8,pos8);
	
	gGeoManager->CloseGeometry();
	top->SetLineColor(kRed);
	gGeoManager->SetTopVisible();
	
	
	// get the nodeid of each detector
	
	gGeoManager->SetCurrentPoint(px1,py1,pz1);
	gGeoManager->FindNode();
	TGeoNode *nodeA1 = gGeoManager->GetCurrentNode();
	Int_t idA1 = gGeoManager->GetCurrentNodeId();
	cout << "A1 " << idA1 << endl;
	
	gGeoManager->SetCurrentPoint(px2,py2,pz2);
	gGeoManager->FindNode();
	TGeoNode *nodeA2 = gGeoManager->GetCurrentNode();
	Int_t idA2 = gGeoManager->GetCurrentNodeId();
	cout << "A2 " << idA2 << endl;
	
	gGeoManager->SetCurrentPoint(px3,py3,pz3);
	gGeoManager->FindNode();
	TGeoNode *nodeA3 = gGeoManager->GetCurrentNode();
	Int_t idA3 = gGeoManager->GetCurrentNodeId();
	cout << "A3 " << idA3 << endl;
	
	gGeoManager->SetCurrentPoint(px4,py4,pz4);
	gGeoManager->FindNode();
	TGeoNode *nodeA4 = gGeoManager->GetCurrentNode();
	Int_t idA4 = gGeoManager->GetCurrentNodeId();
	cout << "A4 " << idA4 << endl;
	
	gGeoManager->SetCurrentPoint(px5,py5,pz5);
	gGeoManager->FindNode();
	TGeoNode *nodeB1 = gGeoManager->GetCurrentNode();
	Int_t idB1 = gGeoManager->GetCurrentNodeId();
	cout << "B1 " << idB1 << endl;
	
	gGeoManager->SetCurrentPoint(px6,py6,pz6);
	gGeoManager->FindNode();
	TGeoNode *nodeB2 = gGeoManager->GetCurrentNode();
	Int_t idB2 = gGeoManager->GetCurrentNodeId();
	cout << "B2 " << idB2 << endl;
	
	gGeoManager->SetCurrentPoint(px7,py7,pz7);
	gGeoManager->FindNode();
	TGeoNode *nodeB3 = gGeoManager->GetCurrentNode();
	Int_t idB3 = gGeoManager->GetCurrentNodeId();
	cout << "B3 " << idB3 << endl;
	
	gGeoManager->SetCurrentPoint(px8,py8,pz8);
	gGeoManager->FindNode();
	TGeoNode *nodeB4 = gGeoManager->GetCurrentNode();
	Int_t idB4 = gGeoManager->GetCurrentNodeId();
	cout << "B4 " << idB4 << endl;
	
	
	
	
	long sec;
	time( &sec );
	TRandom3* ran = new TRandom3((unsigned)sec);
	
	Double_t     u=931.49410242/1000.;
	// unita' di massa (rif. 12C) in GeV/c
	
	// 197Au(4He,4He)197Au
	
	Double_t     mp=4.001506179127; Double_t     zp=2;
	// massa proiettile
	
	Double_t     mt=196.96656879; Double_t     zt=79;
	// massa target
	
	Double_t     m1=4.001506179127;
	// massa prima part. uscente 4He
	
	Double_t     m2=196.96656879;
	// massa seconda part. uscente 197Au
	
	Double_t     rm=mp/mt;
	
	// Q della reazione = 0 perchè scatt el
	Double_t     q2=0.000000000;
	
	//(Momentum, Energy units are Gev/C, GeV)
	Double_t masses[2] = { u*m1, u*m2 } ;
	
	
	Double_t deg=180./TMath::Pi();
	Double_t radius=100.;
	// Int_t geve=0; // including modulation for Rutherford
	
	Double_t step=0.31; //strip width in cm
	Double_t stepA=0.155; //strip width in cm for 32 strip detectors (A group)
	
	Int_t ninc=1e06; //pps
	
	Double_t ntar=1; //atoms per mbarn. 100=target thickness in ug/cm2
	
	Double_t testvar=1;
	Double_t Ta_true=0.;
	
	// theta and phi in degrees
	Double_t Ta=0.;
	Double_t Pa=0.;
	Double_t Tau=0.;
	Double_t Pau=0.;
	
	// theta 4He in the cm
	Double_t Tacm=0.;
	
	// detected energy including det. res. 0.5%
	Double_t Ea=0.;
	Double_t Eau=0.;
	
	// energy loss
	Double_t E1=0.;
	Double_t E2=0.;
	Double_t E3=0.;
	Double_t Eapp=0.;
	Double_t R1=0.;
	Double_t R2=0.;
	Double_t R3=0.;
	Double_t aAu=0.069423;
	Double_t bAu=1.73209;
	Double_t cAu=-2.51568;
	Double_t aAl=0.29978;
	Double_t bAl=3.83143;
	Double_t cAl=-6.64501;
	Double_t aSi=0.339179;
	Double_t bSi=4.46653;
	Double_t cSi=-7.99929;
	
	
	//statistical weight
	Double_t WW;
	
	//weight to get counts per second
	Double_t Wcps=1.;
	
	//flagging particles in detectors (=1 particle detected, =0 particle undetected)
	Int_t AUflag=0;
	Int_t ADflag=0;
	Int_t ALflag=0;
	Int_t ARflag=0;
	Int_t BUflag=0;
	Int_t BDflag=0;
	Int_t BLflag=0;
	Int_t BRflag=0;
	
	Double_t					cos_TA=0.;
	
	
	TFile *fout = new TFile("sim_aau_15.root", "RECREATE");
	
	TTree *sim = new TTree("sim", "sim");
	
	sim->Branch("WW",  &WW,  "WW/D");
	
	sim->Branch("Wcps",  &Wcps,  "Wcps/D");
	
	sim->Branch("Ea", &Ea, "Ea/D"); //26Al energy in MeV
	sim->Branch("Ta", &Ta, "Ta/D"); //26Al theta  in deg
	sim->Branch("Pa", &Pa, "Pa/D"); //26Al phi    in deg
	
	sim->Branch("Eau", &Eau, "Eau/D"); //2H energy in MeV
	sim->Branch("Tau", &Tau, "Tau/D"); //2H theta  in deg
	sim->Branch("Pau", &Pau, "Pau/D"); //2H phi    in deg
	
	sim->Branch("Tacm", &Tacm, "Tacm/D"); //26Al theta in cm  in deg
	
	sim->Branch("BUflag", &BUflag, "BUflag/I");
	sim->Branch("BDflag", &BDflag, "BDflag/I");
	sim->Branch("BRflag", &BRflag, "BRflag/I");
	sim->Branch("BLflag", &BLflag, "BLflag/I");
	sim->Branch("AUflag", &AUflag, "AUflag/I");
	sim->Branch("ADflag", &ADflag, "ADflag/I");
	sim->Branch("ARflag", &ARflag, "ARflag/I");
	sim->Branch("ALflag", &ALflag, "ALflag/I");
	
	sim->Branch("cos_TA", &cos_TA, "cos_TA/D");
	
	
	
	TFile *fhisto=new TFile("h_sim_aau_15.root","RECREATE");
	
	
	// creating arrays of histograms
	
	char pos_names[] = "UDRL";
	
	TH1F *AF[4][32][32];
	for (int a=0;a<4;a++) {
		for (int b=0;b<32;b++){
			for(int c=0;c<32;c++){
				AF[a][b][c] = new TH1F (TString::Format("h_A%c_F%d_B%d", pos_names[a], b, c), TString::Format("h_A%d_F%d_B%d", pos_names[a], b, c), 20000,0,20);
			}
		}
	}
	TH1F *BF[4][16][16];
	for (int a=0;a<4;a++) {
		for (int b=0;b<16;b++){
			for(int c=0;c<16;c++){
				BF[a][b][c] = new TH1F (TString::Format("h_B%c_F%d_B%d", pos_names[a], b, c), TString::Format("h_B%c_F%d_B%d", pos_names[a], b, c), 20000,0,20);
			}
		}
	}
	
	
	Int_t iii=0,jjj=0;
	
	//fired events
	Int_t nnfired=0;
	
	//good events --> if an 26Al has reached detectors
	Int_t ngood=0;
	
	//multiple hit counter (two particles in the same detector)
	Int_t nmulth=0;
	
	//lost alphas counter: there cannot be lost alphas with this trigger
	Int_t lost4Hes=0;
	
	
	TGenPhaseSpace event;
	
	fout->cd();
	
	// Fine setup e inizio Loop di simulazione
	auto end_setup = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> setup_duration = end_setup - start_total;
	cout << "--- Setup completato in: " << setup_duration.count() << " secondi ---" << endl;
	
	auto start_loop = std::chrono::high_resolution_clock::now();
	
	for (int nfired=1; nfired<=ninc; nfired++) {
		
		nnfired++;
		
		Double_t rp = ran->Rndm();
		
		// energia e impulso del proiettile
		Double_t ep = (15 - 0.00597849) / 1000.;        // energy loss
		Double_t pp = sqrt(2. * u * mp * ep);       // momentum of 4He at that energy
		
		// inclinazione del fascio (in gradi)
		Double_t theta_x = 0.0 * TMath::DegToRad();;   // tilt attorno all'asse Y (bending in X)
		Double_t theta_y = 0.0 * TMath::DegToRad();;   // tilt attorno all'asse X (bending in Y)
		
		// shift globale del fascio (in cm)
		Double_t beamShiftX = 0.0;
		Double_t beamShiftY = 0.0;
		
		// beam spot (diametro 0.2 cm → raggio 0.1 cm)
		const Double_t R = 0.1;
		
		// thickness half-target (in um)
		Double_t     thickness=0.023043861;
		
		//-----------------------------------------
		// 4-VETTORI PER CINEMATICA
		//-----------------------------------------
		TLorentzVector target(0.0, 0.0, 0.0, u * mt);
		TLorentzVector beam(0.0, 0.0, pp, u*mp + ep);
		
		// applica inclinazione al fascio
		beam.RotateY(theta_x);
		beam.RotateX(theta_y);
		
		// SDR fascio
		TVector3 ez_beam = beam.Vect().Unit();
		
		TVector3 ref(0,1,0);
		if (fabs(ez_beam.Dot(ref)) > 0.99)
			ref = TVector3(1,0,0);
		
		TVector3 ex_beam = (ref.Cross(ez_beam)).Unit();
		TVector3 ey_beam = ez_beam.Cross(ex_beam).Unit();
		
		
		Double_t rx = ran->Rndm();     // uniformi in [0,1]
		Double_t ry = ran->Rndm();
		
		Double_t dx = -R + rx * (2*R);
		Double_t dy = -R + ry * (2*R);
		
		if (dx*dx + dy*dy >= R*R) continue;   // accept–reject circolare
		
		// posizione finale del beam nel target
		Double_t postx = beamShiftX + dx;
		Double_t posty = beamShiftY + dy;
		
		TVector3 beamOrigin(postx, posty, 0);
		
		TLorentzVector W = beam + target;
		
		event.SetDecay(W, 2, masses);
		
		Double_t weight = event.Generate();
		if (weight == 0) continue;
		
		WW=weight;
		
		TLorentzVector *p1 = event.GetDecay(0); //4He
		TLorentzVector *p2 = event.GetDecay(1); //197Au
		
		
		TVector3 pp1=p1->Vect();
		TVector3 pp2=p2->Vect();
		
		
		TVector3 dir1=pp1.Unit();
		TVector3 dir2=pp2.Unit();
		
		TVector3 dir1_beam = ToBeamFrame(dir1, ex_beam, ey_beam, ez_beam);
		Ta_true = dir1_beam.Theta();   // radianti
		
		// particles tracking
		
		Double_t xx1=dir1(0);
		Double_t yy1=dir1(1);
		Double_t zz1=dir1(2);
		gGeoManager->SetCurrentPoint(postx,posty,0);
		gGeoManager->SetCurrentDirection(xx1,yy1,zz1);
		TGeoNode *current1 = gGeoManager->GetCurrentNode();
		gGeoManager->FindNode();
		gGeoManager->FindNextBoundary(radius);
		Double_t snext1 = gGeoManager->GetStep();
		TGeoNode *newNode1 = gGeoManager->Step();
		Bool_t istate1 = gGeoManager->IsStepEntering();
		Int_t idnode1 = gGeoManager->GetCurrentNodeId();
		
		TGeoNode *current11 = gGeoManager->GetCurrentNode();
		const Double_t *cpoint1 = gGeoManager->GetCurrentPoint();
		TVector3 vhp1(cpoint1[0],cpoint1[1],cpoint1[2]); //hit position of particle 1
		
		
		Double_t xx2=dir2(0);
		Double_t yy2=dir2(1);
		Double_t zz2=dir2(2);
		gGeoManager->SetCurrentPoint(postx,posty,0);
		gGeoManager->SetCurrentDirection(xx2,yy2,zz2);
		TGeoNode *current2 = gGeoManager->GetCurrentNode();
		gGeoManager->FindNode();
		gGeoManager->FindNextBoundary(radius);
		Double_t snext2 = gGeoManager->GetStep();
		TGeoNode *newNode2 = gGeoManager->Step();
		Bool_t istate2 = gGeoManager->IsStepEntering();
		Int_t idnode2 = gGeoManager->GetCurrentNodeId();
		
		TGeoNode *current22 = gGeoManager->GetCurrentNode();
		const Double_t *cpoint2 = gGeoManager->GetCurrentPoint();
		TVector3 vhp2(cpoint2[0],cpoint2[1],cpoint2[2]); //hit position of particle 2
		
		
		
		
		//flagging undetected particles
		Int_t lost1flag=0;
		Int_t lost2flag=0;
		
		//flagging particles in detectors (=1 particle detected, =0 particle undetected)
		AUflag=0;
		ADflag=0;
		ALflag=0;
		ARflag=0;
		BUflag=0;
		BDflag=0;
		BLflag=0;
		BRflag=0;
		
		
		//apparent direction of the particles from detection pixel
		TVector3 vbina;
		TVector3 vbinAu;
		
		
		//4-momenta
		TLorentzVector pb1; //4He 4-momentum in GeV/c
		TLorentzVector pb2; //197Au 4-momentum in GeV/c
		
		
		//vector momenta (parte spaziale del 4-momentum)
		TVector3 ppb1;
		TVector3 ppb2;
		
		// I detect the 4He particle only
		if(istate1 == 1) {
			
			ngood++;
			if(idnode1 == idnode2) {
				nmulth++;
			}
			
			int lost1flag = 0;
			HitData res; // Creiamo una singola variabile per accogliere i risultati
			
			// ==================================================
			// IDENTIFICAZIONE DEL RIVELATORE
			// ==================================================
			if      (idnode1 == idA1) { res = ProcessDetectorHit(vhp1, vA1, rot1, 32, stepA, ran, beamOrigin, ex_beam, ey_beam, ez_beam, -1.0, -1.0); ALflag = 1; }
			else if (idnode1 == idA2) { res = ProcessDetectorHit(vhp1, vA2, rot2, 32, stepA, ran, beamOrigin, ex_beam, ey_beam, ez_beam, -1.0, -1.0); ADflag = 1; }
			else if (idnode1 == idA3) { res = ProcessDetectorHit(vhp1, vA3, rot3, 32, stepA, ran, beamOrigin, ex_beam, ey_beam, ez_beam,  1.0,  1.0); ARflag = 1; }
			else if (idnode1 == idA4) { res = ProcessDetectorHit(vhp1, vA4, rot4, 32, stepA, ran, beamOrigin, ex_beam, ey_beam, ez_beam,  1.0,  1.0); AUflag = 1; }
			else if (idnode1 == idB1) { res = ProcessDetectorHit(vhp1, vB1, rot5, 16, step,  ran, beamOrigin, ex_beam, ey_beam, ez_beam,  1.0, -1.0); BLflag = 1; }
			else if (idnode1 == idB2) { res = ProcessDetectorHit(vhp1, vB2, rot6, 16, step,  ran, beamOrigin, ex_beam, ey_beam, ez_beam,  1.0, -1.0); BDflag = 1; }
			else if (idnode1 == idB3) { res = ProcessDetectorHit(vhp1, vB3, rot7, 16, step,  ran, beamOrigin, ex_beam, ey_beam, ez_beam, -1.0,  1.0); BRflag = 1; }
			else if (idnode1 == idB4) { res = ProcessDetectorHit(vhp1, vB4, rot8, 16, step,  ran, beamOrigin, ex_beam, ey_beam, ez_beam, -1.0,  1.0); BUflag = 1; }
			else {
				// Se la particella non colpisce nessun rivelatore noto
				lost1flag = 1;
				lost4Hes++;
			}
			
			// ==================================================
			// ENERGIA, SMEARING E RIEMPIMENTO ISTOGRAMMI
			// ==================================================
			if (lost1flag == 0) {
				// Estraiamo i dati dalla nostra struttura per comodità
				Ta = res.Ta;
				Pa = res.Pa;
				cos_TA = res.cos_TA;
				iii = res.iii;
				jjj = res.jjj;
				vbina = res.vbina;
				
				Double_t e1 = (p1->Energy() - u * m1); // 4He energy in GeV
				Double_t mp1 = pp1.Mag(); // 4He momentum in GeV/c
				
				// Smearing della risoluzione (usando random C++11)
				std::random_device rd1;
				std::mt19937 generator(rd1());
				std::normal_distribution<double> distribution1(mp1, 0.0001 * mp1);
				Double_t mp1s = distribution1(generator);
				Double_t dmp1 = (mp1s - mp1) / mp1; // variazione percentuale dell'impulso
				
				TVector3 dirb1 = vbina.Unit(); // apparent direction
				TVector3 vpb1 = mp1s * dirb1; // 4He momentum (vector) in GeV/c
				
				pb1.SetPx(vpb1(0));
				pb1.SetPy(vpb1(1));
				pb1.SetPz(vpb1(2));
				pb1.SetE(u * m1 + e1 * (1 + 2 * dmp1));
				ppb1 = pb1.Vect();
				
				Eapp = e1 * (1 + 2 * dmp1) * 1000.; // detected energy in MeV
				Double_t epm = 1000. * ep;
				
				testvar = funcl(zp, zt, epm, Ta_true, rm);
				Wcps = weight * ntar * funcl(zp, zt, epm, Ta_true, rm);
				
				// Calcolo perdite di energia nei vari strati
				R1 = (aAu * Eapp * Eapp + bAu * Eapp + cAu) - (thickness / cos(Ta_true));
				E1 = (-bAu + sqrt(bAu * bAu - 4 * aAu * (cAu - R1))) / (2 * aAu);
				
				R2 = (aAl * E1 * E1 + bAl * E1 + cAl) - (0.3 / cos_TA);
				E2 = (-bAl + sqrt(bAl * bAl - 4 * aAl * (cAl - R2))) / (2 * aAl);
				
				R3 = (aSi * E2 * E2 + bSi * E2 + cSi) - (0.5 / cos_TA);
				E3 = (-bSi + sqrt(bSi * bSi - 4 * aSi * (cSi - R3))) / (2 * aSi);
				
				Ea = E3;
				
				// Riempimento degli array di istogrammi in base al rivelatore colpito
				if      (idnode1 == idA1) AF[3][iii][jjj]->Fill(Ea, Wcps);
				else if (idnode1 == idA2) AF[1][iii][jjj]->Fill(Ea, Wcps);
				else if (idnode1 == idA3) AF[2][iii][jjj]->Fill(Ea, Wcps);
				else if (idnode1 == idA4) AF[0][iii][jjj]->Fill(Ea, Wcps);
				else if (idnode1 == idB1) BF[3][iii][jjj]->Fill(Ea, Wcps);
				else if (idnode1 == idB2) BF[1][iii][jjj]->Fill(Ea, Wcps);
				else if (idnode1 == idB3) BF[2][iii][jjj]->Fill(Ea, Wcps);
				else if (idnode1 == idB4) BF[0][iii][jjj]->Fill(Ea, Wcps);
				else {
					cout << "panic: rivelatore sconosciuto nel riempimento istogrammi" << endl;
				}
			} // fine if lost1flag == 0
			
			// Salvataggio sul file TTree
			fout->cd();
			sim->Fill();
			
			// Stampa di controllo ogni 50.000 eventi
			if (nnfired % 50000 == 0) {
				cout << nnfired << " " << Wcps << " " << testvar << " " << (180 / 3.141592654) * Ta_true << " " << weight << endl;
			}
		}
	} // FINE CICLO FOR PRINCIPALE
	
	// Fine Loop e Calcolo Tempi
	auto end_loop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> loop_duration = end_loop - start_loop;
	std::chrono::duration<double> total_duration = end_loop - start_total;
	
	cout << "--- Simulazione conclusa ---" << endl;
	cout << "Tempo Loop: " << loop_duration.count() << " secondi" << endl;
	cout << "Tempo Totale: " << total_duration.count() << " secondi" << endl;
	cout << "Media per evento: " << (loop_duration.count() / ninc) * 1e6 << " microsecondi/evento" << endl;
	
	// ==================================================
	// SALVATAGGIO FINALE E CHIUSURA FILE
	// ==================================================
	cout << "good events = " << ngood << "     multihit = " << nmulth << "     # proj. = " << nnfired << endl;
	
	fout->cd();
	sim->Write();
	fout->Write();
	
	fhisto->cd();
	for (int a = 0; a < 4; a++) {
		for (int b = 0; b < 32; b++) {
			for(int c = 0; c < 32; c++) {
				AF[a][b][c]->Write();
			}
		}
	}
	for (int a = 0; a < 4; a++) {
		for (int b = 0; b < 16; b++) {
			for(int c = 0; c < 16; c++) {
				BF[a][b][c]->Write();
			}
		}
	}
	
	fout->Close();
	fhisto->Write();
	fhisto->Close();
	
	cout << "Finished" << endl;
}
