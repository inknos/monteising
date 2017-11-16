void compileMacro(TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  
  //
  //gSystem->Load("libgomp");
  TString cmd( gSystem->GetMakeSharedLib() );
  cmd.ReplaceAll("g++","g++ -fopenmp");
  gSystem->SetMakeSharedLib(cmd); 
  gSystem->CompileMacro("./Lattice.cxx",opt.Data());  // load class MHasSHA256
  //
}
