#include "chisquaretest.h"
#include "newchisquaretest.C"
#include "newchisquaretestpp.C"
void run_chi2()
{

    // Don't forget to change everything in the header file
    // Don't forget to change params in header!!! Like binning ...
    // shift_-0.5_0.2_smear_0.0_0.015_modified_signal_21_21_100.root
    // shift_-0.2_0.0_smear_0.0_0.006_modified_signal_21_21_100.root
    // jobID_3_shift_-0.5_0.2_smear_0.0_0.015_modified_signal_21_21_1000.root

    /*newchisquaretest(1, 0);
    newchisquaretest(2, 0);
    newchisquaretest(3, 0);
    newchisquaretest(4, 0);
    newchisquaretest(5, 0);
    newchisquaretest(6, 0);
    newchisquaretest(7, 0);


    newchisquaretest(1, 1);
    newchisquaretest(2, 1);
    newchisquaretest(3, 1);
    newchisquaretest(4, 1);
    newchisquaretest(5, 1);
    newchisquaretest(6, 1);
    newchisquaretest(7, 1);*/


    newchisquaretestpp(1, 0);
    newchisquaretestpp(2, 0);
    newchisquaretestpp(3, 0);
    newchisquaretestpp(4, 0);
    newchisquaretestpp(5, 0);
    newchisquaretestpp(6, 0);
    newchisquaretestpp(7, 0);


    newchisquaretestpp(1, 1);
    newchisquaretestpp(2, 1);
    newchisquaretestpp(3, 1);
    newchisquaretestpp(4, 1);
    newchisquaretestpp(5, 1);
    newchisquaretestpp(6, 1);
    newchisquaretestpp(7, 1);


    cout << "Finished" << endl;

}