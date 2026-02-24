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

void geometry_and_physics_aau_15() {

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
	
	Int_t ninc=1e07; //pps
	
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
	
	
	
	Int_t iii,jjj;
	
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
			
			
			if( idnode1==idnode2 )
				{
				nmulth++;
				}
			
			//searching for 4He apparent trajectory
	
			if      (idnode1==idA1) {
				//cout << "4He in A1" << endl;
				TVector3 vdA1=vhp1-vA1;
				//vector from det. A1 center to hit point
				
				Double_t vvdA1[3]={vdA1(0),vdA1(1),vdA1(2)};
				Double_t vvdA1r[3]={0.,0.,0.};
				
				//cout << vdA1(0) << " " << vdA1(1) << " " << vdA1(2) << " " << vdA1.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot1->MasterToLocal(vvdA1,vvdA1r);
				//for A1, B1 +largeZ, B16 -largeZ
				//        F1 -largeX, F16 +largeX
				
				TVector3 vdA1r(vvdA1r[0],vvdA1r[1],vvdA1r[2]);
				
				//cout << vdA1r(0) << " " << vdA1r(1) << " " << vdA1r(2) << " " << vdA1r.Mag() << endl;
				
				Double_t detA1f=vvdA1r[0];
				Double_t detA1b=vvdA1r[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<32; ii++){
					if(detA1f>=-2.48+ii*stepA && detA1f<-2.48+(ii+1)*stepA) break;
				}
				for (jj=0; jj<32; jj++){
					if(detA1b>=-2.48+jj*stepA && detA1b<-2.48+(jj+1)*stepA) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfA1=-2.48+ii*stepA+stepA*rxxs;
				Double_t cbA1=-2.48+jj*stepA+stepA*ryys;
				
				iii=ii;
				jjj=jj;
				
				//binned (apparent->a) position on detector A1
				Double_t avdA1r[3]={cfA1,vvdA1r[1],cbA1};
				Double_t avdA1[3]={0.,0.,0.};
				Double_t avdN[3]={0.,-1.,0.};
				Double_t avN[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot1->LocalToMaster(avdA1r,avdA1);
				rot1->LocalToMaster(avdN,avN);
				
				//this is the vector from det.A1 center to apparent hit point
				TVector3 vdA1a(avdA1[0],avdA1[1],avdA1[2]);
				TVector3 vdN(avN[0],avN[1],avN[2]);
				
				//cout << vdA1a(0) << " " << vdA1a(1) << " " << vdA1a(2) << " " << vdA1a.Mag() << endl;
				
				//apparent trajectory of 4He
				TVector3 vhit=vA1+vdA1a;		// punto di hit globale
				TVector3 vbin_lab = vhit - beamOrigin;   // direzione reale dalla reazione
				TVector3 vbinA1 = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
				Ta=(vbinA1.Theta())*deg;
				Pa=(vbinA1.Phi())*deg;
				vbina=vbinA1;
				TVector3 ubinA1=vbinA1.Unit();
				TVector3 ubinN=ToBeamFrame(vdN, ex_beam, ey_beam, ez_beam).Unit();
				cos_TA=ubinN.Dot(ubinA1);
				ALflag=1;
			}
			
			else if (idnode1==idA2) {
				//cout << "4He in A2" << endl;
				TVector3 vdA2=vhp1-vA2;
				//vector from det. A2 center to hit point
				
				Double_t vvdA2[3]={vdA2(0),vdA2(1),vdA2(2)};
				Double_t vvdA2r[3]={0.,0.,0.};
				
				//cout << vdA2(0) << " " << vdA2(1) << " " << vdA2(2) << " " << vdA2.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot2->MasterToLocal(vvdA2,vvdA2r);
				//for A2, B1 +largeZ, B16 -largeZ
				//        F1 -largeX, F16 +largeX
				
				TVector3 vdA2r(vvdA2r[0],vvdA2r[1],vvdA2r[2]);
				
				//cout << vdA2r(0) << " " << vdA2r(1) << " " << vdA2r(2) << " " << vdA2r.Mag() << endl;
				
				Double_t detA2f=vvdA2r[0];
				Double_t detA2b=vvdA2r[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<32; ii++){
					if(detA2f>=-2.48+ii*stepA && detA2f<-2.48+(ii+1)*stepA) break;
				}
				for (jj=0; jj<32; jj++){
					if(detA2b>=-2.48+jj*stepA && detA2b<-2.48+(jj+1)*stepA) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfA2=-2.48+ii*stepA+stepA*rxxs;
				Double_t cbA2=-2.48+jj*stepA+stepA*ryys;
				
				iii=ii;
				jjj=jj;
				
				//binned (apparent->a) position on detector A2
				Double_t avdA2r[3]={cfA2,vvdA2r[1],cbA2};
				Double_t avdA2[3]={0.,0.,0.};
				Double_t avdN[3]={0.,-1.,0.};
				Double_t avN[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot2->LocalToMaster(avdA2r,avdA2);
				rot2->LocalToMaster(avdN,avN);
				
				//this is the vector from det.A2 center to apparent hit point
				TVector3 vdA2a(avdA2[0],avdA2[1],avdA2[2]);
				TVector3 vdN(avN[0],avN[1],avN[2]);
				
				//cout << vdA2a(0) << " " << vdA2a(1) << " " << vdA2a(2) << " " << vdA2a.Mag() << endl;
				
				//apparent trajectory of 4He
				TVector3 vhit=vA2+vdA2a;		// punto di hit globale
				TVector3 vbin_lab = vhit - beamOrigin;   // direzione reale dalla reazione
				TVector3 vbinA2 = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
				Ta=(vbinA2.Theta())*deg;
				Pa=(vbinA2.Phi())*deg;
				vbina=vbinA2;
				TVector3 ubinA2=vbinA2.Unit();
				TVector3 ubinN=ToBeamFrame(vdN, ex_beam, ey_beam, ez_beam).Unit();
				cos_TA=ubinN.Dot(ubinA2);
				ADflag=1;
			}
			
			else if (idnode1==idA3) {
				//cout << "4He in A3" << endl;
				TVector3 vdA3=vhp1-vA3;
				//vector from det. A3 center to hit point
				
				Double_t vvdA3[3]={vdA3(0),vdA3(1),vdA3(2)};
				Double_t vvdA3r[3]={0.,0.,0.};
				
				//cout << vdA3(0) << " " << vdA3(1) << " " << vdA3(2) << " " << vdA3.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot3->MasterToLocal(vvdA3,vvdA3r);
				//for A3, B1 +largeZ, B16 -largeZ
				//        F1 +largeX, F16 -largeX
				
				TVector3 vdA3r(vvdA3r[0],vvdA3r[1],vvdA3r[2]);
				
				//cout << vdA3r(0) << " " << vdA3r(1) << " " << vdA3r(2) << " " << vdA3r.Mag() << endl;
				
				Double_t detA3f=vvdA3r[0];
				Double_t detA3b=vvdA3r[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<32; ii++){
					if(detA3f<2.48-ii*stepA && detA3f>=2.48-(ii+1)*stepA) break;
				}
				for (jj=0; jj<32; jj++){
					if(detA3b>=-2.48+jj*stepA && detA3b<-2.48+(jj+1)*stepA) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfA3=2.48-ii*stepA-stepA*rxxs;
				Double_t cbA3=-2.48+jj*stepA+stepA*ryys;
				
				iii=ii;
				jjj=jj;
				
				//binned (apparent->a) position on detector A3
				Double_t avdA3r[3]={cfA3,vvdA3r[1],cbA3};
				Double_t avdA3[3]={0.,0.,0.};
				Double_t avdN[3]={0.,1.,0.};
				Double_t avN[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot3->LocalToMaster(avdA3r,avdA3);
				rot3->LocalToMaster(avdN,avN);
				
				//this is the vector from det.A3 center to apparent hit point
				TVector3 vdA3a(avdA3[0],avdA3[1],avdA3[2]);
				TVector3 vdN(avN[0],avN[1],avN[2]);
				
				//cout << vdA3a(0) << " " << vdA3a(1) << " " << vdA3a(2) << " " << vdA3a.Mag() << endl;
				
				//apparent trajectory of 4He
				TVector3 vhit=vA3+vdA3a;		// punto di hit globale
				TVector3 vbin_lab = vhit - beamOrigin;   // direzione reale dalla reazione
				TVector3 vbinA3 = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
				Ta=(vbinA3.Theta())*deg;
				Pa=(vbinA3.Phi())*deg;
				vbina=vbinA3;
				TVector3 ubinA3=vbinA3.Unit();
				TVector3 ubinN=ToBeamFrame(vdN, ex_beam, ey_beam, ez_beam).Unit();
				cos_TA=ubinN.Dot(ubinA3);
				ARflag=1;
			}
			
			else if (idnode1==idA4) {
				//cout << "4He in A4" << endl;
				TVector3 vdA4=vhp1-vA4;
				//vector from det. A4 center to hit point
				
				Double_t vvdA4[3]={vdA4(0),vdA4(1),vdA4(2)};
				Double_t vvdA4r[3]={0.,0.,0.};
				
				//cout << vdA4(0) << " " << vdA4(1) << " " << vdA4(2) << " " << vdA4.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot4->MasterToLocal(vvdA4,vvdA4r);
				//for A4, B1 +largeZ, B16 -largeZ
				//        F1 +largeX, F16 -largeX
				
				TVector3 vdA4r(vvdA4r[0],vvdA4r[1],vvdA4r[2]);
				
				//cout << vdA4r(0) << " " << vdA4r(1) << " " << vdA4r(2) << " " << vdA4r.Mag() << endl;
				
				Double_t detA4f=vvdA4r[0];
				Double_t detA4b=vvdA4r[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<32; ii++){
					if(detA4f<2.48-ii*stepA && detA4f>=2.48-(ii+1)*stepA) break;
				}
				for (jj=0; jj<32; jj++){
					if(detA4b>=-2.48+jj*stepA && detA4b<-2.48+(jj+1)*stepA) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfA4=2.48-ii*stepA-stepA*rxxs;
				Double_t cbA4=-2.48+jj*stepA+stepA*ryys;
				
				iii=ii;
				jjj=jj;
				
				//binned (apparent->a) position on detector A4
				Double_t avdA4r[3]={cfA4,vvdA4r[1],cbA4};
				Double_t avdA4[3]={0.,0.,0.};
				Double_t avdN[3]={0.,1.,0.};
				Double_t avN[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot4->LocalToMaster(avdA4r,avdA4);
				rot4->LocalToMaster(avdN,avN);
				
				//this is the vector from det.A4 center to apparent hit point
				TVector3 vdA4a(avdA4[0],avdA4[1],avdA4[2]);
				TVector3 vdN(avN[0],avN[1],avN[2]);
				
				//cout << vdA4a(0) << " " << vdA4a(1) << " " << vdA4a(2) << " " << vdA4a.Mag() << endl;
				
				//apparent trajectory of 4He
				TVector3 vhit=vA4+vdA4a;		// punto di hit globale
				TVector3 vbin_lab = vhit - beamOrigin;   // direzione reale dalla reazione
				TVector3 vbinA4 = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
				Ta=(vbinA4.Theta())*deg;
				Pa=(vbinA4.Phi())*deg;
				vbina=vbinA4;
				TVector3 ubinA4=vbinA4.Unit();
				TVector3 ubinN=ToBeamFrame(vdN, ex_beam, ey_beam, ez_beam).Unit();
				cos_TA=ubinN.Dot(ubinA4);
				AUflag=1;
				
			}
			
			else if (idnode1==idB1) {
				//cout << "4He in B1" << endl;
				TVector3 vdB1=vhp1-vB1;
				//vector from det. B1 center to hit point
				
				Double_t vvdB1[3]={vdB1(0),vdB1(1),vdB1(2)};
				Double_t vvdB1r[3]={0.,0.,0.};
				
				//cout << vdB1(0) << " " << vdB1(1) << " " << vdB1(2) << " " << vdB1.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot5->MasterToLocal(vvdB1,vvdB1r);
				//for B1, B1 +largeZ, B16 -largeZ
				//        F1 -largeX, F16 +largeX
				
				TVector3 vdB1r(vvdB1r[0],vvdB1r[1],vvdB1r[2]);
				
				//cout << vdB1r(0) << " " << vdB1r(1) << " " << vdB1r(2) << " " << vdB1r.Mag() << endl;
				
				Double_t detB1f=vvdB1r[0];
				Double_t detB1b=vvdB1r[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<16; ii++){
					if(detB1f>=-2.48+ii*step && detB1f<-2.48+(ii+1)*step) break;
				}
				for (jj=0; jj<16; jj++){
					if(detB1b>=-2.48+jj*step && detB1b<-2.48+(jj+1)*step) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfB1=-2.48+ii*step+step*rxxs;
				Double_t cbB1=-2.48+jj*step+step*ryys;
				
				iii=ii;
				jjj=jj;
				
				//binned (apparent->a) position on detector B1
				Double_t avdB1r[3]={cfB1,vvdB1r[1],cbB1};
				Double_t avdB1[3]={0.,0.,0.};
				Double_t avdN[3]={0.,1.,0.};
				Double_t avN[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot5->LocalToMaster(avdB1r,avdB1);
				rot5->LocalToMaster(avdN,avN);
				
				//this is the vector from det.B1 center to apparent hit point
				TVector3 vdB1a(avdB1[0],avdB1[1],avdB1[2]);
				TVector3 vdN(avN[0],avN[1],avN[2]);
				
				//cout << vdB1a(0) << " " << vdB1a(1) << " " << vdB1a(2) << " " << vdB1a.Mag() << endl;
				
				//apparent trajectory of 4He
				TVector3 vhit=vB1+vdB1a;		// punto di hit globale
				TVector3 vbin_lab = vhit - beamOrigin;   // direzione reale dalla reazione
				TVector3 vbinB1 = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
				Ta=(vbinB1.Theta())*deg;
				Pa=(vbinB1.Phi())*deg;
				vbina=vbinB1;
				TVector3 ubinB1=vbinB1.Unit();
				TVector3 ubinN=ToBeamFrame(vdN, ex_beam, ey_beam, ez_beam).Unit();
				cos_TA=ubinN.Dot(ubinB1);
				BLflag=1;
			}
			
			else if (idnode1==idB2) {
				//cout << "4He in B2" << endl;
				TVector3 vdB2=vhp1-vB2;
				//vector from det. B2 center to hit point
				
				Double_t vvdB2[3]={vdB2(0),vdB2(1),vdB2(2)};
				Double_t vvdB2r[3]={0.,0.,0.};
				
				//cout << vdB2(0) << " " << vdB2(1) << " " << vdB2(2) << " " << vdB2.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot6->MasterToLocal(vvdB2,vvdB2r);
				//for B2, B1 +largeZ, B16 -largeZ
				//        F1 -largeX, F16 +largeX
				
				TVector3 vdB2r(vvdB2r[0],vvdB2r[1],vvdB2r[2]);
				
				//cout << vdB2r(0) << " " << vdB2r(1) << " " << vdB2r(2) << " " << vdB2r.Mag() << endl;
				
				Double_t detB2f=vvdB2r[0];
				Double_t detB2b=vvdB2r[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<16; ii++){
					if(detB2f>=-2.48+ii*step && detB2f<-2.48+(ii+1)*step) break;
				}
				for (jj=0; jj<16; jj++){
					if(detB2b>=-2.48+jj*step && detB2b<-2.48+(jj+1)*step) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfB2=-2.48+ii*step+step*rxxs;
				Double_t cbB2=-2.48+jj*step+step*ryys;
				
				iii=ii;
				jjj=jj;
				
				//binned (apparent->a) position on detector B2
				Double_t avdB2r[3]={cfB2,vvdB2r[1],cbB2};
				Double_t avdB2[3]={0.,0.,0.};
				Double_t avdN[3]={0.,1.,0.};
				Double_t avN[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot6->LocalToMaster(avdB2r,avdB2);
				rot6->LocalToMaster(avdN,avN);
				
				//this is the vector from det.B2 center to apparent hit point
				TVector3 vdB2a(avdB2[0],avdB2[1],avdB2[2]);
				TVector3 vdN(avN[0],avN[1],avN[2]);
				
				//cout << vdB2a(0) << " " << vdB2a(1) << " " << vdB2a(2) << " " << vdB2a.Mag() << endl;
				
				//apparent trajectory of 4He
				TVector3 vhit=vB2+vdB2a;		// punto di hit globale
				TVector3 vbin_lab = vhit - beamOrigin;   // direzione reale dalla reazione
				TVector3 vbinB2 = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
				Ta=(vbinB2.Theta())*deg;
				Pa=(vbinB2.Phi())*deg;
				vbina=vbinB2;
				TVector3 ubinB2=vbinB2.Unit();
				TVector3 ubinN=ToBeamFrame(vdN, ex_beam, ey_beam, ez_beam).Unit();
				cos_TA=ubinN.Dot(ubinB2);
				BDflag=1;
			}
			
			else if (idnode1==idB3) {
				//cout << "4He in B3" << endl;
				TVector3 vdB3=vhp1-vB3;
				//vector from det. B3 center to hit point
				
				Double_t vvdB3[3]={vdB3(0),vdB3(1),vdB3(2)};
				Double_t vvdB3r[3]={0.,0.,0.};
				
				//cout << vdB3(0) << " " << vdB3(1) << " " << vdB3(2) << " " << vdB3.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot7->MasterToLocal(vvdB3,vvdB3r);
				//for B3, B1 +largeZ, B16 -largeZ
				//        F1 +largeX, F16 -largeX
				
				TVector3 vdB3r(vvdB3r[0],vvdB3r[1],vvdB3r[2]);
				
				//cout << vdB3r(0) << " " << vdB3r(1) << " " << vdB3r(2) << " " << vdB3r.Mag() << endl;
				
				Double_t detB3f=vvdB3r[0];
				Double_t detB3b=vvdB3r[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<16; ii++){
					if(detB3f<2.48-ii*step && detB3f>=2.48-(ii+1)*step) break;
				}
				for (jj=0; jj<16; jj++){
					if(detB3b>=-2.48+jj*step && detB3b<-2.48+(jj+1)*step) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfB3=2.48-ii*step-step*rxxs;
				Double_t cbB3=-2.48+jj*step+step*ryys;
				
				iii=ii;
				jjj=jj;
				
				//binned (apparent->a) position on detector B3
				Double_t avdB3r[3]={cfB3,vvdB3r[1],cbB3};
				Double_t avdB3[3]={0.,0.,0.};
				Double_t avdN[3]={0.,-1.,0.};
				Double_t avN[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot7->LocalToMaster(avdB3r,avdB3);
				rot7->LocalToMaster(avdN,avN);
				
				//this is the vector from det.B3 center to apparent hit point
				TVector3 vdB3a(avdB3[0],avdB3[1],avdB3[2]);
				TVector3 vdN(avN[0],avN[1],avN[2]);
				
				//cout << vdB3a(0) << " " << vdB3a(1) << " " << vdB3a(2) << " " << vdB3a.Mag() << endl;
				
				//apparent trajectory of 4He
				TVector3 vhit=vB3+vdB3a;		// punto di hit globale
				TVector3 vbin_lab = vhit - beamOrigin;   // direzione reale dalla reazione
				TVector3 vbinB3 = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
				Ta=(vbinB3.Theta())*deg;
				Pa=(vbinB3.Phi())*deg;
				vbina=vbinB3;
				TVector3 ubinB3=vbinB3.Unit();
				TVector3 ubinN=ToBeamFrame(vdN, ex_beam, ey_beam, ez_beam).Unit();
				cos_TA=ubinN.Dot(ubinB3);
				BRflag=1;
			}
			
			else if (idnode1==idB4) {
				//cout << "4He in B4" << endl;
				TVector3 vdB4=vhp1-vB4;
				//vector from det. B4 center to hit point
				
				Double_t vvdB4[3]={vdB4(0),vdB4(1),vdB4(2)};
				Double_t vvdB4r[3]={0.,0.,0.};
				
				//cout << vdB4(0) << " " << vdB4(1) << " " << vdB4(2) << " " << vdB4.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot8->MasterToLocal(vvdB4,vvdB4r);
				//for B4, B1 +largeZ, B16 -largeZ
				//        F1 +largeX, F16 -largeX
				
				TVector3 vdB4r(vvdB4r[0],vvdB4r[1],vvdB4r[2]);
				
				//cout << vdB4r(0) << " " << vdB4r(1) << " " << vdB4r(2) << " " << vdB4r.Mag() << endl;
				
				Double_t detB4f=vvdB4r[0];
				Double_t detB4b=vvdB4r[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<16; ii++){
					if(detB4f<2.48-ii*step && detB4f>=2.48-(ii+1)*step) break;
				}
				for (jj=0; jj<16; jj++){
					if(detB4b>=-2.48+jj*step && detB4b<-2.48+(jj+1)*step) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfB4=2.48-ii*step-step*rxxs;
				Double_t cbB4=-2.48+jj*step+step*ryys;
				
				iii=ii;
				jjj=jj;
				
				//binned (apparent->a) position on detector B4
				Double_t avdB4r[3]={cfB4,vvdB4r[1],cbB4};
				Double_t avdB4[3]={0.,0.,0.};
				Double_t avdN[3]={0.,-1.,0.};
				Double_t avN[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot8->LocalToMaster(avdB4r,avdB4);
				rot8->LocalToMaster(avdN,avN);
				
				//this is the vector from det.B4 center to apparent hit point
				TVector3 vdB4a(avdB4[0],avdB4[1],avdB4[2]);
				TVector3 vdN(avN[0],avN[1],avN[2]);
				
				//cout << vdB4a(0) << " " << vdB4a(1) << " " << vdB4a(2) << " " << vdB4a.Mag() << endl;
				
				//apparent trajectory of 4He
				TVector3 vhit=vB4+vdB4a;		// punto di hit globale
				TVector3 vbin_lab = vhit - beamOrigin;   // direzione reale dalla reazione
				TVector3 vbinB4 = ToBeamFrame(vbin_lab, ex_beam, ey_beam, ez_beam);
				Ta=(vbinB4.Theta())*deg;
				Pa=(vbinB4.Phi())*deg;
				vbina=vbinB4;
				TVector3 ubinB4=vbinB4.Unit();
				TVector3 ubinN=ToBeamFrame(vdN, ex_beam, ey_beam, ez_beam).Unit();
				cos_TA=ubinN.Dot(ubinB4);
				BUflag=1;
			}
			
			else{
				// cout << "4He lost in space" << endl;
				lost1flag=1;
				lost4Hes++;
			}
			
			
			
			
			//searching for 197Au apparent trajectory
	/*
			if      (idnode2==idA1){
				//cout << "197Au in A1" << endl;
				TVector3 vdA1x=vhp2-vA1;
				//vector from det. A1 center to hit point
				
				Double_t vvdA1x[3]={vdA1x(0),vdA1x(1),vdA1x(2)};
				Double_t vvdA1xr[3]={0.,0.,0.};
				
				//cout << vdA1x(0) << " " << vdA1x(1) << " " << vdA1x(2) << " " << vdA1x.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot1->MasterToLocal(vvdA1x,vvdA1xr);
				//for A1, B1 +largeZ, B16 -largeZ
				//        F1 -largeX, F16 +largeX
				
				TVector3 vdA1xr(vvdA1xr[0],vvdA1xr[1],vvdA1xr[2]);
				
				//cout << vdA1xr(0) << " " << vdA1xr(1) << " " << vdA1xr(2) << " " << vdA1xr.Mag() << endl;
				
				Double_t detA1xf=vvdA1xr[0];
				Double_t detA1xb=vvdA1xr[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<32; ii++){
					if(detA1xf>=-2.48+ii*stepA && detA1xf<-2.48+(ii+1)*stepA) break;
				}
				for (jj=0; jj<32; jj++){
					if(detA1xb<2.48-jj*stepA && detA1xb>=2.48-(jj+1)*stepA) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfA1x=-2.48+ii*stepA+stepA*rxxs;
				Double_t cbA1x=2.48-jj*stepA-stepA*ryys;
				
				//binned (apparent->a) position on detector A1
				Double_t avdA1xr[3]={cfA1x,vvdA1xr[1],cbA1x};
				Double_t avdA1x[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot1->LocalToMaster(avdA1xr,avdA1x);
				
				//this is the vector from det.A1 center to apparent hit point
				TVector3 vdA1xa(avdA1x[0],avdA1x[1],avdA1x[2]);
				
				//cout << vdA1xa(0) << " " << vdA1xa(1) << " " << vdA1xa(2) << " " << vdA1xa.Mag() << endl;
				
				//apparent trajectory of 197Au
				TVector3 vbinA1x=vA1+vdA1xa;
				Tau=(vbinA1x.Theta())*deg;
				Pau=(vbinA1x.Phi())*deg;
				vbinAu=vbinA1x;
				
			}
			else if (idnode2==idA2) {
				//cout << "197Au in A2" << endl;
				TVector3 vdA2x=vhp2-vA2;
				//vector from det. A2 center to hit point
				
				Double_t vvdA2x[3]={vdA2x(0),vdA2x(1),vdA2x(2)};
				Double_t vvdA2xr[3]={0.,0.,0.};
				
				//cout << vdA2x(0) << " " << vdA2x(1) << " " << vdA2x(2) << " " << vdA2x.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot2->MasterToLocal(vvdA2x,vvdA2xr);
				//for A2, B1 +largeZ, B16 -largeZ
				//        F1 -largeX, F16 +largeX
				
				TVector3 vdA2xr(vvdA2xr[0],vvdA2xr[1],vvdA2xr[2]);
				
				//cout << vdA2xr(0) << " " << vdA2xr(1) << " " << vdA2xr(2) << " " << vdA2xr.Mag() << endl;
				
				Double_t detA2xf=vvdA2xr[0];
				Double_t detA2xb=vvdA2xr[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<32; ii++){
					if(detA2xf>=-2.48+ii*stepA && detA2xf<-2.48+(ii+1)*stepA) break;
				}
				for (jj=0; jj<32; jj++){
					if(detA2xb<2.48-jj*stepA && detA2xb>=2.48-(jj+1)*stepA) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfA2x=-2.48+ii*stepA+stepA*rxxs;
				Double_t cbA2x=2.48-jj*stepA-stepA*ryys;
				
				//binned (apparent->a) position on detector A2
				Double_t avdA2xr[3]={cfA2x,vvdA2xr[1],cbA2x};
				Double_t avdA2x[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot2->LocalToMaster(avdA2xr,avdA2x);
				
				//this is the vector from det.A2 center to apparent hit point
				TVector3 vdA2xa(avdA2x[0],avdA2x[1],avdA2x[2]);
				
				//cout << vdA2xa(0) << " " << vdA2xa(1) << " " << vdA2xa(2) << " " << vdA2xa.Mag() << endl;
				
				//apparent trajectory of 197Au
				TVector3 vbinA2x=vA2+vdA2xa;
				Tau=(vbinA2x.Theta())*deg;
				Pau=(vbinA2x.Phi())*deg;
				vbinAu=vbinA2x;
				
			}
			else if (idnode2==idA3) {
				//cout << "197Au in A3" << endl;
				TVector3 vdA3x=vhp2-vA3;
				//vector from det. A3 center to hit point
				
				Double_t vvdA3x[3]={vdA3x(0),vdA3x(1),vdA3x(2)};
				Double_t vvdA3xr[3]={0.,0.,0.};
				
				//cout << vdA3x(0) << " " << vdA3x(1) << " " << vdA3x(2) << " " << vdA3x.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot3->MasterToLocal(vvdA3x,vvdA3xr);
				//for A3, B1 +largeZ, B16 -largeZ
				//        F1 +largeX, F16 -largeX
				
				TVector3 vdA3xr(vvdA3xr[0],vvdA3xr[1],vvdA3xr[2]);
				
				//cout << vdA3xr(0) << " " << vdA3xr(1) << " " << vdA3xr(2) << " " << vdA3xr.Mag() << endl;
				
				Double_t detA3xf=vvdA3xr[0];
				Double_t detA3xb=vvdA3xr[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<32; ii++){
					if(detA3xf<2.48-ii*stepA && detA3xf>=2.48-(ii+1)*stepA) break;
				}
				for (jj=0; jj<32; jj++){
					if(detA3xb<2.48-jj*stepA && detA3xb>=2.48-(jj+1)*stepA) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfA3x=2.48-ii*stepA-stepA*rxxs;
				Double_t cbA3x=2.48-jj*stepA-stepA*ryys;
				
				//binned (apparent->a) position on detector A3
				Double_t avdA3xr[3]={cfA3x,vvdA3xr[1],cbA3x};
				Double_t avdA3x[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot3->LocalToMaster(avdA3xr,avdA3x);
				
				//this is the vector from det.A3 center to apparent hit point
				TVector3 vdA3xa(avdA3x[0],avdA3x[1],avdA3x[2]);
				
				//cout << vdA3xa(0) << " " << vdA3xa(1) << " " << vdA3xa(2) << " " << vdA3xa.Mag() << endl;
				
				//apparent trajectory of 197Au
				TVector3 vbinA3x=vA3+vdA3xa;
				Tau=(vbinA3x.Theta())*deg;
				Pau=(vbinA3x.Phi())*deg;
				vbinAu=vbinA3x;
				
			}
			else if (idnode2==idA4) {
				//cout << "197Au in A4" << endl;
				TVector3 vdA4x=vhp2-vA4;
				//vector from det. A4 center to hit point
				
				Double_t vvdA4x[3]={vdA4x(0),vdA4x(1),vdA4x(2)};
				Double_t vvdA4xr[3]={0.,0.,0.};
				
				//cout << vdA4x(0) << " " << vdA4x(1) << " " << vdA4x(2) << " " << vdA4x.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot4->MasterToLocal(vvdA4x,vvdA4xr);
				//for A4, B1 +largeZ, B16 -largeZ
				//        F1 +largeX, F16 -largeX
				
				TVector3 vdA4xr(vvdA4xr[0],vvdA4xr[1],vvdA4xr[2]);
				
				//cout << vdA4xr(0) << " " << vdA4xr(1) << " " << vdA4xr(2) << " " << vdA4xr.Mag() << endl;
				
				Double_t detA4xf=vvdA4xr[0];
				Double_t detA4xb=vvdA4xr[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<32; ii++){
					if(detA4xf<2.48-ii*stepA && detA4xf>=2.48-(ii+1)*stepA) break;
				}
				for (jj=0; jj<32; jj++){
					if(detA4xb<2.48-jj*stepA && detA4xb>=2.48-(jj+1)*stepA) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfA4x=2.48-ii*stepA-stepA*rxxs;
				Double_t cbA4x=2.48-jj*stepA-stepA*ryys;
				
				//binned (apparent->a) position on detector A4
				Double_t avdA4xr[3]={cfA4x,vvdA4xr[1],cbA4x};
				Double_t avdA4x[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot4->LocalToMaster(avdA4xr,avdA4x);
				
				//this is the vector from det.A4 center to apparent hit point
				TVector3 vdA4xa(avdA4x[0],avdA4x[1],avdA4x[2]);
				
				//cout << vdA4xa(0) << " " << vdA4xa(1) << " " << vdA4xa(2) << " " << vdA4xa.Mag() << endl;
				
				//apparent trajectory of 197Au
				TVector3 vbinA4x=vA4+vdA4xa;
				Tau=(vbinA4x.Theta())*deg;
				Pau=(vbinA4x.Phi())*deg;
				vbinAu=vbinA4x;
				
			}
			else if (idnode2==idB1) {
				//cout << "197Au in B1" << endl;
				TVector3 vdB1x=vhp2-vB1;
				//vector from det. B1 center to hit point
				
				Double_t vvdB1x[3]={vdB1x(0),vdB1x(1),vdB1x(2)};
				Double_t vvdB1xr[3]={0.,0.,0.};
				
				//cout << vdB1x(0) << " " << vdB1x(1) << " " << vdB1x(2) << " " << vdB1x.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot5->MasterToLocal(vvdB1x,vvdB1xr);
				//for B1, B1 +largeZ, B16 -largeZ
				//        F1 -largeX, F16 +largeX
				
				TVector3 vdB1xr(vvdB1xr[0],vvdB1xr[1],vvdB1xr[2]);
				
				//cout << vdB1xr(0) << " " << vdB1xr(1) << " " << vdB1xr(2) << " " << vdB1xr.Mag() << endl;
				
				Double_t detB1xf=vvdB1xr[0];
				Double_t detB1xb=vvdB1xr[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<16; ii++){
					if(detB1xf>=-2.48+ii*step && detB1xf<-2.48+(ii+1)*step) break;
				}
				for (jj=0; jj<16; jj++){
					if(detB1xb<2.48-jj*step && detB1xb>=2.48-(jj+1)*step) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfB1x=-2.48+ii*step+step*rxxs;
				Double_t cbB1x=2.48-jj*step-step*ryys;
				
				//binned (apparent->a) position on detector B1
				Double_t avdB1xr[3]={cfB1x,vvdB1xr[1],cbB1x};
				Double_t avdB1x[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot5->LocalToMaster(avdB1xr,avdB1x);
				
				//this is the vector from det.B1 center to apparent hit point
				TVector3 vdB1xa(avdB1x[0],avdB1x[1],avdB1x[2]);
				
				//cout << vdB1xa(0) << " " << vdB1xa(1) << " " << vdB1xa(2) << " " << vdB1xa.Mag() << endl;
				
				//apparent trajectory of 197Au
				TVector3 vbinB1x=vB1+vdB1xa;
				Tau=(vbinB1x.Theta())*deg;
				Pau=(vbinB1x.Phi())*deg;
				vbinAu=vbinB1x;
				
			}
			else if (idnode2==idB2) {
				//cout << "197Au in B2" << endl;
				TVector3 vdB2x=vhp2-vB2;
				//vector from det. B2 center to hit point
				
				Double_t vvdB2x[3]={vdB2x(0),vdB2x(1),vdB2x(2)};
				Double_t vvdB2xr[3]={0.,0.,0.};
				
				//cout << vdB2x(0) << " " << vdB2x(1) << " " << vdB2x(2) << " " << vdB2x.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot6->MasterToLocal(vvdB2x,vvdB2xr);
				//for B2, B1 +largeZ, B16 -largeZ
				//        F1 -largeX, F16 +largeX
				
				TVector3 vdB2xr(vvdB2xr[0],vvdB2xr[1],vvdB2xr[2]);
				
				//cout << vdB2xr(0) << " " << vdB2xr(1) << " " << vdB2xr(2) << " " << vdB2xr.Mag() << endl;
				
				Double_t detB2xf=vvdB2xr[0];
				Double_t detB2xb=vvdB2xr[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<16; ii++){
					if(detB2xf>=-2.48+ii*step && detB2xf<-2.48+(ii+1)*step) break;
				}
				for (jj=0; jj<16; jj++){
					if(detB2xb<2.48-jj*step && detB2xb>=2.48-(jj+1)*step) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfB2x=-2.48+ii*step+step*rxxs;
				Double_t cbB2x=2.48-jj*step-step*ryys;
				
				//binned (apparent->a) position on detector B2
				Double_t avdB2xr[3]={cfB2x,vvdB2xr[1],cbB2x};
				Double_t avdB2x[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot6->LocalToMaster(avdB2xr,avdB2x);
				
				//this is the vector from det.B2 center to apparent hit point
				TVector3 vdB2xa(avdB2x[0],avdB2x[1],avdB2x[2]);
				
				//cout << vdB2xa(0) << " " << vdB2xa(1) << " " << vdB2xa(2) << " " << vdB2xa.Mag() << endl;
				
				//apparent trajectory of 197Au
				TVector3 vbinB2x=vB2+vdB2xa;
				Tau=(vbinB2x.Theta())*deg;
				Pau=(vbinB2x.Phi())*deg;
				vbinAu=vbinB2x;
				
			}
			else if (idnode2==idB3) {
				//cout << "197Au in B3" << endl;
				TVector3 vdB3x=vhp2-vB3;
				//vector from det. B3 center to hit point
				
				Double_t vvdB3x[3]={vdB3x(0),vdB3x(1),vdB3x(2)};
				Double_t vvdB3xr[3]={0.,0.,0.};
				
				//cout << vdB3x(0) << " " << vdB3x(1) << " " << vdB3x(2) << " " << vdB3x.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot7->MasterToLocal(vvdB3x,vvdB3xr);
				//for B3, B1 +largeZ, B16 -largeZ
				//        F1 +largeX, F16 -largeX
				
				TVector3 vdB3xr(vvdB3xr[0],vvdB3xr[1],vvdB3xr[2]);
				
				//cout << vdB3xr(0) << " " << vdB3xr(1) << " " << vdB3xr(2) << " " << vdB3xr.Mag() << endl;
				
				Double_t detB3xf=vvdB3xr[0];
				Double_t detB3xb=vvdB3xr[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<16; ii++){
					if(detB3xf<2.48-ii*step && detB3xf>=2.48-(ii+1)*step) break;
				}
				for (jj=0; jj<16; jj++){
					if(detB3xb<2.48-jj*step && detB3xb>=2.48-(jj+1)*step) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfB3x=2.48-ii*step-step*rxxs;
				Double_t cbB3x=2.48-jj*step-step*ryys;
				
				//binned (apparent->a) position on detector B3
				Double_t avdB3xr[3]={cfB3x,vvdB3xr[1],cbB3x};
				Double_t avdB3x[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot7->LocalToMaster(avdB3xr,avdB3x);
				
				//this is the vector from det.B3 center to apparent hit point
				TVector3 vdB3xa(avdB3x[0],avdB3x[1],avdB3x[2]);
				
				//cout << vdB3xa(0) << " " << vdB3xa(1) << " " << vdB3xa(2) << " " << vdB3xa.Mag() << endl;
				
				//apparent trajectory of 197Au
				TVector3 vbinB3x=vB3+vdB3xa;
				Tau=(vbinB3x.Theta())*deg;
				Pau=(vbinB3x.Phi())*deg;
				vbinAu=vbinB3x;
				
			}
			else if (idnode2==idB4) {
				//cout << "197Au in B4" << endl;
				TVector3 vdB4x=vhp2-vB4;
				//vector from det. B4 center to hit point
				
				Double_t vvdB4x[3]={vdB4x(0),vdB4x(1),vdB4x(2)};
				Double_t vvdB4xr[3]={0.,0.,0.};
				
				//cout << vdB4x(0) << " " << vdB4x(1) << " " << vdB4x(2) << " " << vdB4x.Mag() << endl;
				
				//this is the transformation to planar configuration
				rot8->MasterToLocal(vvdB4x,vvdB4xr);
				//for B4, B1 +largeZ, B16 -largeZ
				//        F1 +largeX, F16 -largeX
				
				TVector3 vdB4xr(vvdB4xr[0],vvdB4xr[1],vvdB4xr[2]);
				
				//cout << vdB4xr(0) << " " << vdB4xr(1) << " " << vdB4xr(2) << " " << vdB4xr.Mag() << endl;
				
				Double_t detB4xf=vvdB4xr[0];
				Double_t detB4xb=vvdB4xr[2];
				
				Int_t ii,jj;
				
				//search for the pixel f->front, b->back
				for (ii=0; ii<16; ii++){
					if(detB4xf<2.48-ii*step && detB4xf>=2.48-(ii+1)*step) break;
				}
				for (jj=0; jj<16; jj++){
					if(detB4xb<2.48-jj*step && detB4xb>=2.48-(jj+1)*step) break;
				}
				Double_t rxxs = ran->Rndm();
				Double_t ryys = ran->Rndm();
				Double_t cfB4x=2.48-ii*step-step*rxxs;
				Double_t cbB4x=2.48-jj*step-step*ryys;
				
				//binned (apparent->a) position on detector B4
				Double_t avdB4xr[3]={cfB4x,vvdB4xr[1],cbB4x};
				Double_t avdB4x[3]={0.,0.,0.};
				
				//we go back to the real position in the space
				rot8->LocalToMaster(avdB4xr,avdB4x);
				
				//this is the vector from det.B4 center to apparent hit point
				TVector3 vdB4xa(avdB4x[0],avdB4x[1],avdB4x[2]);
				
				//cout << vdB4xa(0) << " " << vdB4xa(1) << " " << vdB4xa(2) << " " << vdB4xa.Mag() << endl;
				
				//apparent trajectory of 197Au
				TVector3 vbinB4x=vB4+vdB4xa;
				Tau=(vbinB4x.Theta())*deg;
				Pau=(vbinB4x.Phi())*deg;
				vbinAu=vbinB4x;
				
			}
			else{
				//cout << "197Au lost in space" << endl;
				
				lost2flag=1;
			}
			
			*/
			
			//associating the energy accounting for det. res.
			
			std::default_random_engine generator;
			
			if(lost1flag==0)
				{
				
				Double_t e1=(p1->Energy()-u*m1); //4He energy in GeV
				Double_t mp1=pp1.Mag(); //4He momentum in GeV/c
				std::random_device rd1;
				std::mt19937 generator(rd1());
				std::normal_distribution<double> distribution1(/*mean=*/mp1, /*stddev=*/0.0001*mp1);
				Double_t mp1s = distribution1(generator);
				Double_t dmp1=(mp1s-mp1)/mp1; //variazione percentuale dell'impulso
				TVector3 dirb1=vbina.Unit(); //apparent direction
				TVector3 vpb1=mp1s*dirb1; //4He momentum (vector) in GeV/c
				
				pb1.SetPx(vpb1(0));
				pb1.SetPy(vpb1(1));
				pb1.SetPz(vpb1(2));
				pb1.SetE(u*m1+e1*(1+2*dmp1));
				ppb1=pb1.Vect();
				
				// Ea=e1*(1+2*dmp1)*1000.;
				Eapp=e1*(1+2*dmp1)*1000.; //detected energy of 4He if we negliect energy losses
				
				Double_t epm=1000.*ep;
				
				testvar=funcl(zp,zt,epm,Ta_true,rm);
				Wcps=weight*ntar*funcl(zp,zt,epm,Ta_true,rm);
				// Wcps=1.;
				
					// energy loss in half target of 197Au
					R1=(aAu*Eapp*Eapp+bAu*Eapp+cAu)-(thickness/cos(Ta_true));
					E1=(-bAu+sqrt(bAu*bAu-4*aAu*(cAu-R1)))/(2*aAu);
					// cout << "R1: " << R1 << endl;
					// cout << "E1: " << E1 << endl;
					
					// energy loss in 0.3um Al in dead layer
					R2=(aAl*E1*E1+bAl*E1+cAl)-(0.3/cos_TA);
					E2=(-bAl+sqrt(bAl*bAl-4*aAl*(cAl-R2)))/(2*aAl);
					// cout << "R2: " << R2 << endl;
					// cout << "E2: " << E2 << endl;
					
					// energy loss in 0.5um Si in dead layer
					R3=(aSi*E2*E2+bSi*E2+cSi)-(0.5/cos_TA);
					E3=(-bSi+sqrt(bSi*bSi-4*aSi*(cSi-R3)))/(2*aSi);
					// cout << "R2: " << R2 << endl;
					// cout << "E2: " << E2 << endl;
					
					Ea=E3;
				
				
				if      (idnode1==idA1) {
					AF[3][iii][jjj]->Fill(Ea, Wcps);
				}
				else if (idnode1==idA2) {
					AF[1][iii][jjj]->Fill(Ea, Wcps);
				}
				else if (idnode1==idA3) {
					AF[2][iii][jjj]->Fill(Ea, Wcps);
				}
				else if (idnode1==idA4) {
					AF[0][iii][jjj]->Fill(Ea, Wcps);
				}
				else if (idnode1==idB1) {
					BF[3][iii][jjj]->Fill(Ea, Wcps);
				}
				else if (idnode1==idB2) {
					BF[1][iii][jjj]->Fill(Ea, Wcps);
				}
				else if (idnode1==idB3) {
					BF[2][iii][jjj]->Fill(Ea, Wcps);
				}
				else if (idnode1==idB4) {
					BF[0][iii][jjj]->Fill(Ea, Wcps);
				}
				else{
					cout << "panic" << endl;
				}
				}
			
			/*
			
			if(lost2flag==0)
				{
				
				Double_t e2=(p2->Energy()-u*m2); //197Au energy in GeV
				Double_t mp2=pp2.Mag(); //197Au momentum in GeV/c
				std::random_device rd2;
				std::mt19937 generator(rd2());
				std::normal_distribution<double> distribution2(mp2, 0.0025*mp2);
				Double_t mp2s = distribution2(generator);
				Double_t dmp2=(mp2s-mp2)/mp2; //variazione percentuale dell'impulso
				TVector3 dirb2=vbinAu.Unit(); //apparent direction
				TVector3 vpb2=mp2s*dirb2; //197Au momentum (vector) in GeV/c
				
				pb2.SetPx(vpb2(0));
				pb2.SetPy(vpb2(1));
				pb2.SetPz(vpb2(2));
				pb2.SetE(u*m2+e2*(1+2*dmp2));
				ppb2=pb2.Vect();
				Eau=e2*(1+2*dmp2)*1000.; //detected energy of 197Au
				
				fhisto->cd();
				h2_t_p_Au->Fill(Pau,Tau,weight);
				h2_e_t_Au->Fill(Eau,Tau,weight);
				
				h2_t_p_Au_w->Fill(Pau,Tau,Wcps);
				h2_e_t_Au_w->Fill(Eau,Tau,Wcps);
				}
			
			if(lost2flag==1)
				{
				Eau=0;
				Tau=0;
				Pau=0;
				}
			*/
			
			fout->cd();
			sim->Fill();
			
			
			if(nnfired%50000==0) cout << nnfired << " " << Wcps << " " << testvar << " " << (180/3.141592654)*Ta_true << " " << weight << endl;
			
		}
		
	}
	
	cout << "good events = " << ngood << "     multihit = " << nmulth << "     # proj. = " << nnfired << endl;
	
	fout->cd();
	
	sim->Write();
	
	fout->Write();
	
	fhisto->cd();
	
	for (int a=0;a<4;a++) {
		for (int b=0;b<32;b++){
			for(int c=0;c<32;c++){
				AF[a][b][c]->Write();
			}
		}
	}
	
	for (int a=0;a<4;a++) {
		for (int b=0;b<16;b++){
			for(int c=0;c<16;c++){
				BF[a][b][c]->Write();
			}
		}
	}
	
	
	fout->Close();
	fhisto->Write();
	fhisto->Close();
	
	cout << "Finished" << endl;
	
}
