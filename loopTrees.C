int loopTrees( char *name1 = "Tracks.root", char *name2 = "output.root")
 {
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("StUtilities");
  gSystem->Load("StIOMaker");
  gSystem->Load("StarClassLibrary");
  gSystem->Load("StEvent");
  gSystem->Load("MyAnalysisMaker");

  printf("Input file:  %s\n", name1);
  printf("Output file: %s\n", name2);

  MyAnalysisMaker  LoopMaker("Loop");
  return LoopMaker.doLoop(name1, name2);

  gROOT->Reset();
}
