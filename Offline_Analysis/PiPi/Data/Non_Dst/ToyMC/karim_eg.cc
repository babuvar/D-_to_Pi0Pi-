void karim_eg(Int_t nsig_toy_in=1000, Int_t nbkg_toy_in=1000) {

  Int_t ntoys(0), save_it(1);
  Int_t ntoys_in=1;  
  ntoys=ntoys_in;

  char *s = new char[1];

  RooRandom::randomGenerator()->SetSeed(1);

  // VARIABLES
  RooRealVar de("de","de",1.68,2.06,"GeV");
  de.setBins(40);

  // PDF (Signal):
  RooRealVar coreFractor("coreFrac", "fraction of core Gaussian", 0.83763);
  RooRealVar SigMean("SigMean", "mean Signal E", 1.862, "GeV"),
    SignalSig1("SignalSig1", "delta E resolution (1)", 0.020, "GeV"),
    SignalSig2("SignalSig2", "delta E resolution (2)", 0.020, "GeV");
  RooGaussian Sig1("signal1", "signal Gaussian (1)", 
		   de, SigMean, SignalSig1);
  RooGaussian Sig2("signal2", "signal Gaussian (2)", 
		   de, SigMean, SignalSig2);
  RooAddPdf signal("signal", "signal double Gaussian", 
		   Sig1, Sig2, coreFractor);

  // PDF (BCKG):
  RooRealVar dELin("dELin", "delta E slope", -3, -10., 10.0);
  RooPolynomial bkgDeltaE("bkgDeltaE","background Delta E", 
			  de, RooArgList(dELin));

  RooRealVar nSig("nSig", "Number of Signal Events", 200000, 10000.0, 1000000.);
  RooRealVar nBkg("nBkg", "Number of Background Events",2000000., 0.0, 10000000);

  RooAddPdf model("model", "Model", 
		  RooArgList(bkgDeltaE,signal), 
		  RooArgList(nBkg,nSig)); 

  // from here, just testing toy generating...
  Double_t nsig_toy(0), nbkg_toy(0), delin_toy(0);
  nbkg_toy = nbkg_toy_in;
  delin_toy = -0.30;
  nsig_toy =nsig_toy_in;  

  dELin.setVal(delin_toy);
  dELin.setConstant(kTRUE);

  char t_dat[5]=".dat";
  char str_toy[1];
  sprintf(str_toy, "%d", nsig_toy_in);
  for (Int_t i = 0; i<ntoys; i++) {

    cout << " TOY # " << i << endl;

    RooDataSet *ToyMC = new RooDataSet("ToyMC", "ToyMC",RooArgSet(de));

    Int_t nsig_toy_ex = RooRandom::randomGenerator()->Poisson(nsig_toy);
    if (nsig_toy_ex>0) {
      RooDataSet *ToyMC4 = signal.generate(RooArgSet(de), nsig_toy_ex);
      ToyMC->append(*ToyMC4);
    }

    Int_t nbkg_toy_ex = RooRandom::randomGenerator()->Poisson(nbkg_toy);
    if (nbkg_toy_ex>0) {
      RooDataSet *ToyMC3 = bkgDeltaE.generate(RooArgSet(de), nbkg_toy_ex);
      ToyMC->append(*ToyMC3);
    }

    std::cout << " PRODUCE: nsig= " << nsig_toy_ex 
	      << ", nbkg= " << nbkg_toy_ex
	      << " -> total= " << ToyMC->numEntries()
	      << std::endl;

    if (save_it){
      char t_file[50];
      strcpy(t_file,"./samples.dat");
      //strcat(t_file,t_dat);
      cout << " t_file = " << t_file << endl;
      ToyMC->write(t_file);
    }
  }

  //dELin.setConstant(kFALSE);
  RooDataSet *data = RooDataSet::read("./samples.dat",RooArgSet(de));
  nSig.setVal(3);
  nSig.setError(0.1);  
  nBkg.setVal(11.0);

  RooFitResult* fitRes_toy = model.fitTo(*data,"e");
 
  RooPlot *dEPlot = data->plotOn(de.frame(30));
  model.plotOn(dEPlot);
  /*
  model.plotOn(dEPlot,Components(RooArgList(bkgDeltaE)),
	       LineStyle(kDashed));  
  model.plotOn(dEPlot,Components(RooArgList(signal)),
	       LineColor(kRed),LineStyle(4));
  */
  model.paramOn(dEPlot,data);
	
  TCanvas *c = new TCanvas("c","Delta E fit", 500, 400);
  c->Divide(1,1);
  c->cd(1); 
  dEPlot->SetTitle("");
  dEPlot->Draw();
  c->Print("de_data_show.eps");

}
