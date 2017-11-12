void compileMacro(TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  
  //
  gSystem->CompileMacro("./Lattice.cxx",opt.Data());  // load class MHasSHA256
  //
}
