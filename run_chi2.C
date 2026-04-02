#include "chisquaretest.h"
#include "newchisquaretest.C"
#include "newchisquaretestpp.C"
void run_chi2(bool firsttime = 0)
{

    // Don't forget to change everything in the header file
    // Don't forget to change params in header!!! Like binning ...
    // shift_-0.5_0.2_smear_0.0_0.015_modified_signal_21_21_100.root
    // shift_-0.2_0.0_smear_0.0_0.006_modified_signal_21_21_100.root
    // jobID_3_shift_-0.5_0.2_smear_0.0_0.015_modified_signal_21_21_1000.root

    newchisquaretest(1, 2, firsttime);
    // newchisquaretest(2, 2, firsttime);
    // newchisquaretest(3, 2, firsttime);
    // newchisquaretest(4, 2, firsttime);
    // newchisquaretest(5, 2, firsttime);
    // newchisquaretest(6, 2, firsttime);
    // newchisquaretest(7, 2, firsttime);
    // newchisquaretest(8, 2, firsttime);
    // newchisquaretest(9, 2, firsttime);

    // don't run pp until I regenerate template (since now cent == 4 means 0-100)
    // need to add 1D pT as syst FA_nominal_inclusive_1D_pT

    //newchisquaretestpp(1, 2);
    //newchisquaretestpp(2, 2);
    //newchisquaretestpp(3, 2);
    //newchisquaretestpp(4, 2);
    //newchisquaretestpp(5, 2);
    //newchisquaretestpp(6, 2);
    //newchisquaretestpp(7, 2);
    //newchisquaretestpp(8, 2);

    cout << "Finished" << endl;
}