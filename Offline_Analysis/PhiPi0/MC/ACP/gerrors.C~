void gerrors() {
   //Draw a graph with error bars
   // To see the output of this macro, click begin_html <a href="gif/gerrors.gif">here</a>. end_html
   //Author: Rene Brun
   
   TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);

 //  c1->SetFillColor(42);
//   c1->SetGrid();
//   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);

   const Int_t n = 10;
   Float_t y[n]  = {-0.09,-0.82,0.83,1.63,-0.29,-0.28};
   Float_t x[n]  = {-0.05,0.95,1.95,2.95,3.95,4.95};
   Float_t ey[n] = {0.84,0.85,0.84,0.84,0.83,0.84};
   Float_t ex[n] = {0,0,0,0,0,0};
   TGraphErrors *gr = new TGraphErrors(n,x,y,ex,ey);

   const Int_t n3 = 10;
   Float_t y3[n3]  = {0.16,-0.96,1.04,1.50,-0.25,-0.22};
   Float_t x3[n3]  = {0.05,1.05,2.05,3.05,4.05,5.05};
   Float_t ey3[n3] = {0.83,0.84,0.83,0.83,0.82,0.83};
   Float_t ex3[n3] = {0,0,0,0,0,0};
   TGraphErrors *gr3 = new TGraphErrors(n3,x3,y3,ex3,ey3);


  const Int_t n1 = 10;
   Float_t y1[n1]  = {-0.47};
   Float_t x1[n1]  = {5.95};
   Float_t ey1[n1] = {0.86};
   Float_t ex1[n1] = {0};
   TGraphErrors *gr1 = new TGraphErrors(n1,x1,y1,ex1,ey1);

  const Int_t n4 = 10;
   Float_t y4[n4]  = {-0.57};
   Float_t x4[n4]  = {6.05};
   Float_t ey4[n4] = {0.85};
   Float_t ex4[n4] = {0};
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

   gr1->SetTitle("Measured A_{CP}^{uncorr}");
   gr1->SetMarkerColor(2);
   gr1->SetMarkerStyle(21);

  gr2->SetTitle("Measured A_{CP}^{uncorr}");
   gr2->SetMarkerColor(0);
   gr2->SetMarkerStyle(21);


   gr3->SetMarkerColor(8);
   gr3->SetMarkerStyle(21);


   gr4->SetMarkerColor(6);
   gr4->SetMarkerStyle(21);



   gr->SetMarkerSize(1.4);
   gr1->SetMarkerSize(1.4);
   gr2->SetMarkerSize(1.6);
   gr3->SetMarkerSize(1.4);
   gr4->SetMarkerSize(1.4);

//   gr1->Draw("AP","SAME");
//   gr2->Draw("AP","SAME");

     TMultiGraph *mg = new TMultiGraph();
     mg->SetTitle("Measured A_{CP}^{uncorr}");
     mg->Add(gr,"p");
     mg->Add(gr1,"p");
  //   mg->Add(gr2,"p");
     mg->Add(gr3,"p");
     mg->Add(gr4,"p");
     mg->Add(gr2,"p");
     mg->Draw("A");

  leg = new TLegend(0.6,0.7,0.89,0.89);

  leg->AddEntry(gr,"MC (1-D Fit)","lep");
//  leg->AddEntry(h1,"D^{#pm}_{s}","f");
//  leg->AddEntry(h2,"D^{#pm}","f");
  leg->AddEntry(gr1,"Data (1-D Fit)","lep");
  leg->AddEntry(gr3,"MC (2-D Fit)","lep");
  leg->AddEntry(gr4,"Data (2-D Fit)","lep");
  leg->Draw();



   c1->Update();
}
