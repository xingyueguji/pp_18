#include "chisquaretest.h"
#include "newchisquaretest.C"
#include "newchisquaretestpp.C"
void run_chi2(bool isNew = 1)
{

    // Don't forget to change everything in the header file
    // Don't forget to change params in header!!! Like binning ...
    // shift_-0.5_0.2_smear_0.0_0.015_modified_signal_21_21_100.root
    // shift_-0.2_0.0_smear_0.0_0.006_modified_signal_21_21_100.root
    // jobID_3_shift_-0.5_0.2_smear_0.0_0.015_modified_signal_21_21_1000.root

    newchisquaretest(0, 0, isNew);
    newchisquaretest(1, 1, isNew);
    newchisquaretest(0, 1, isNew);
    newchisquaretest(1, 0, isNew);

    newchisquaretestpp(0, 1, isNew);
    newchisquaretestpp(0, 0, isNew);
    newchisquaretestpp(1, 1, isNew);
    newchisquaretestpp(1, 0, isNew);
}