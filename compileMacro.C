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
  /*  //Parallel
  TString cmd( gSystem->GetMakeSharedLib() );
  cmd.ReplaceAll("g++","g++ -fopenmp");
  gSystem->SetMakeSharedLib(cmd);
  */

  TString workingDir(gSystem->WorkingDirectory());

  gSystem->SetBuildDir(workingDir + TString("/build"), true);

  TString includePath("-I$ROOTSYS/include -I");
  includePath += workingDir + TString("/include");
  gSystem->SetIncludePath(includePath);

  TString srcPath(workingDir);
  srcPath += TString("/src/");
  //
  gSystem->CompileMacro(srcPath + TString("Lattice.cxx"), opt.Data());  // load class Lattice
  gSystem->CompileMacro(srcPath + TString("DrawLattice.cxx"), opt.Data());
  gSystem->CompileMacro(srcPath + TString("Block.cxx"), opt.Data());
  gSystem->CompileMacro(srcPath + TString("SimulationLattice.cxx"), opt.Data());
  gSystem->CompileMacro(srcPath + TString("AnalysisLattice.cxx"), opt.Data());
  //eja
}
