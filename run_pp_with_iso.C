#include "pp_with_iso.C"

void run_pp_with_iso(){
    pp_with_iso *s1 = new pp_with_iso();
    s1->Loop(true);
    cout << "All Finished!" << endl;
}