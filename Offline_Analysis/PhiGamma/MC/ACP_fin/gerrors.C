void gerrors() {
   //Draw a graph with error bars
   // To see the output of this macro, click begin_html <a href="gif/gerrors.gif">here</a>. end_html
   //Author: Rene Brun
   
   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);

 //  c1->SetFillColor(42);
//   c1->SetGrid();
//   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);


//2Da
   const Int_t n = 10;
   Float_t y[n]  = {5.0,-5.1,2.4,-11.0,13.1,10.0};
   Float_t x[n]  = {-0.1,0.9,1.9,2.9,3.9,4.9};
   Float_t ey[n] = {7.3,7.6,7.3,7.2,7.4,7.3};
   Float_t ex[n] = {0,0,0,0,0,0};
   TGraphErrors *gr = new TGraphErrors(n,x,y,ex,ey);

//2Db
   const Int_t n3 = 10;
   Float_t y3[n3]  = {5.8,-2.2,-3.3,-16.8,8.6,4.4};
   Float_t x3[n3]  = {0.0,1.0,2.0,3.0,4.0,5.0};
   Float_t ey3[n3] = {7.0,7.1,6.7,6.8,6.7,7.0};
   Float_t ex3[n3] = {0,0,0,0,0,0};
   TGraphErrors *gr3 = new TGraphErrors(n3,x3,y3,ex3,ey3);

//3D
   const Int_t n4 = 10;
   Float_t y4[n4]  = {7.9,-1.2,1.5,-14.5,13.5,3.8};
   Float_t x4[n4]  = {0.1,1.1,2.1,3.1,4.1,5.1};
   Float_t ey4[n4] = {6.2,6.4,6.2,6.2,6.2,6.2};
   Float_t ex4[n4] = {0,0,0,0,0,0};
   TGraphErrors *gr4 = new TGraphErrors(n4,x4,y4,ex4,ey4);



  const Int_t n2 = 10;
   Float_t y2[n2]  = {0};
   Float_t x2[n2]  = {0};
   Float_t ey2[n2] = {0};
   Float_t ex2[n2] = {0};
   TGraphErrors *gr2 = new TGraphErrors(n2,x2,y2,ex2,ey2);




   gr->SetTitle("Measured A_{CP}^{uncorr}");
   gr->SetMarkerColor(9);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");


  gr2->SetTitle("Measured A_{CP}^{uncorr}");
   gr2->SetMarkerColor(0);
   gr2->SetMarkerStyle(21);


  gr4->SetTitle("Measured A_{CP}^{uncorr}");
   gr4->SetMarkerColor(2);
   gr4->SetMarkerStyle(21);

   gr3->SetMarkerColor(8);
   gr3->SetMarkerStyle(21);


   gr->SetMarkerSize(1.4);
   gr2->SetMarkerSize(1.6);
   gr3->SetMarkerSize(1.4);
   gr4->SetMarkerSize(1.4);

//   gr2->Draw("AP","SAME");

     TMultiGraph *mg = new TMultiGraph();
     mg->SetTitle("Measured A_{CP}^{uncorr}");
     mg->Add(gr,"p");

  //   mg->Add(gr2,"p");
     mg->Add(gr3,"p");
     mg->Add(gr4,"p");
     mg->Add(gr2,"p");
     mg->Draw("A");

  leg = new TLegend(0.6,0.7,0.89,0.89);

  leg->AddEntry(gr,"#DeltaM-cos #theta_{hel} Fit (MC)  ","lep");
//  leg->AddEntry(h1,"D^{#pm}_{s}","f");
//  leg->AddEntry(h2,"D^{#pm}","f");
  leg->AddEntry(gr3,"M_{D}-cos #theta_{hel} Fit (MC)  ","lep");
  leg->AddEntry(gr4,"3-D Fit (MC)","lep");
  leg->Draw();



   c1->Update();
}









