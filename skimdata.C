#include "SkimNew.C"

void skimdata(){
    SkimNew *s1 = new SkimNew();
    s1->ToggleBranches();
    s1->Loop();
    cout << "All Finished!" << endl;
}