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
  //TString cmd( gSystem->GetMakeSharedLib() );
  //cmd.ReplaceAll("g++","g++ -fopenmp");
  //gSystem->SetMakeSharedLib(cmd);

  gSystem->CompileMacro("./Lattice.cxx", opt.Data());  // load class Lattice
  //gInterpreter->GenerateDictionary("vector<Track&gt","Track.h;vector");
  //gInterpreter->GenerateDictionary("Lattice","Lattice.h");
  gSystem->CompileMacro("./DrawLattice.cxx", opt.Data());
  //
}
