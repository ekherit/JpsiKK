{
  TFile f("sample.root");
  ((TTree*)f.Get("event"))->MakeClass("RootEvent");
  ((TTree*)f.Get("mc"))->MakeClass("RootMC");
  ((TTree*)f.Get("mdc"))->MakeClass("RootMdc");
  ((TTree*)f.Get("emc"))->MakeClass("RootEmc");
  ((TTree*)f.Get("tof"))->MakeClass("RootTof");
  ((TTree*)f.Get("dedx"))->MakeClass("RootDedx");
  ((TTree*)f.Get("muc"))->MakeClass("RootMuc");
  ((TTree*)f.Get("neutral"))->MakeClass("RootNeutral");
  //((TTree*)f.Get("mctopo"))->MakeClass("RootMCTopo");
}
