void testdraw()
{
   Int_t n=20;
   Double_t x[n],y[n];
   for (Int_t i=0; i < n; i++) {
      x[i]=i*0.1;
      y[i]=10*sin(x[i]+0.2);
   }
   TGraph *gr1 = new TGraph (n,x,y);
   TAxis *axis = gr1->GetXaxis();

   axis->SetLimits(0.,5.);                 // along X
   gr1->SetMaximum(20.);   // along          
   gr1->SetMinimum(-20.);  //   Y     

   gr1->Draw("AC*");
}
