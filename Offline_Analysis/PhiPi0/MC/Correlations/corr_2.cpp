#include<cmath>
void corr_2(){

    TCanvas *c1 = new TCanvas();



  TChain* chain=new TChain("h1");
  chain->Add("PhiPi0_4s_GMC_0.root");

/*
	c1->cd();
	h1->SetTitle("Plot");
	h1->SetMarkerColor(kPink);
	h1->Draw("Dzeromas : Deltam","Mygenfla != 1");
	h1->SetMarkerColor(kBlue);
	h1->Draw("Dzeromas : Deltam","Mygenfla == 1","same");*/

/*
	c1->cd();
	h1->SetTitle("Plot");
	h1->SetMarkerColor(kPink);
	h1->Draw("Deltam : Dzeromas","Mygenfla != 1");
	h1->SetMarkerColor(kBlue);
	h1->Draw("Deltam : Dzeromas","Mygenfla == 1","same");
*/

/*
	c1->cd();	h1->SetTitle("Plot");
	h1->SetMarkerColor(kPink);
	h1->Draw("Pizmomen : Dstps","Mygenfla != 1");
	h1->SetMarkerColor(kBlue);
	h1->Draw("Pizmomen : Dstps","Mygenfla == 1");//,"same");
	h1->SetTitle("Plot");
*/
	h1->Draw("Dzeromas","Mygenfla == 1");


    gPad->Update();

}

