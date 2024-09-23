void get_run_groups(){
	auto df = ROOT::RDF::FromCSV("filtered_output.csv");

    // Extract the column "col1" into a std::vector
    auto run_number = df.Take<Long64_t>("run_number");
    auto lumi = df.Take<Double_t>("recorded_lumi");

    // Now you have std::vector<int> and std::vector<float>
    // You can use them as you like

	Double_t sum_lumi = std::accumulate(lumi.begin(), lumi.end(),0.0f);
	cout << "Sum Lumi is " << sum_lumi << endl;

	Double_t expected_lumi_per_group = sum_lumi / 25;
	cout << "Expected lumi per group is " << expected_lumi_per_group << endl;

	Double_t num_runs = run_number->size();
	cout << "Total num of runs is " << num_runs << endl;

	Double_t current_accu_lumi = 0;

	Int_t counter = 0;

	Int_t RunrangeUpper[22] = {};
	Int_t RunrangeLower[22] = {};

	for (int i = 0; i < num_runs; ++i){
		if (current_accu_lumi == 0) RunrangeUpper[counter] = run_number->at(i);
		current_accu_lumi += lumi->at(i);

		if (current_accu_lumi >= expected_lumi_per_group){
			cout << "lumi for this group is " << current_accu_lumi << endl;
			current_accu_lumi = 0;
			RunrangeLower[counter] = run_number->at(i);
			counter++;
		}

	}
	
	for (int i = 0; i < 22; i++){
		cout << RunrangeUpper[21-i] << ",";
	}
	cout << "    " << endl;

	for (int i = 0; i < 22; i++){
		cout << RunrangeLower[21-i] << ",";
	}
	cout << "    " << endl;

}