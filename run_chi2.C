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
    cout << "we are here " << endl;
    newchisquaretest(1, 1);*/

    newchisquaretest(1, 2);
    newchisquaretest(2, 2);
    newchisquaretest(3, 2);
    newchisquaretest(4, 2);
    newchisquaretest(5, 2);

    newchisquaretest(7, 2);
    newchisquaretest(8, 2);
    newchisquaretest(9, 2);

    /*newchisquaretestpp(1, 2);
    newchisquaretestpp(2, 2);
    newchisquaretestpp(3, 2);
    newchisquaretestpp(4, 2);
    newchisquaretestpp(5, 2);

    newchisquaretestpp(7, 2);*/


    cout << "Finished" << endl;

}