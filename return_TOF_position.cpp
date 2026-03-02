#include "./return_TOF_position.h"

using namespace std;

TOF_reconstructor::TOF_reconstructor() : _v_eff(181.974),
	scintillator_bias({-0.0253795,
		0.307779,
		-0.0568409,
		0.437488,
		-0.137346,
		-0.350588,
		-0.153047,
		-0.374798
	}), 
	y_positions({60.025, 42.875, 25.725, 8.575, -8.575, -25.725, -42.875, -60.025}),
	x_dimensions({40.92, 94.0, 112.0, 123.0, 123.0, 112.0, 94.0, 40.92}),
	v_eff_uncertainty(33.7644),
	sigma_sipm_i({0.333076, 0.360282, 0.356124, 0.304019, 0.226525, 0.243068, 0.306572, 0.273243}),
	sigma_sipm_i_uncertainties({0.0101467, 0.0490062, 0.0703601, 0.0993561, 0.133235, 0.103015, 0.0575791, 0.0123247}),
	corresponding_pairs({{0, 8}, {1, 10}, {2, 11}, {3, 12}, {4, 13}, {5, 14}, {6, 15}, {7, 16}, {8, 0}, {10, 1}, {11, 2}, {12, 3}, {13, 4}, {14, 5}, {15, 6}, {16, 7}}),
	TOF_IDs({0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16})

{}

TOF_reconstructor::TOF_reconstructor(double v_eff)
	: 	TOF_reconstructor()
		{
			_v_eff = v_eff;
		} 

void TOF_reconstructor::setVeff(double v){ _v_eff = v; }

double TOF_reconstructor::GetVeff() const { return _v_eff; }

double TOF_reconstructor::GetScintDimensionX(int i){
	return x_dimensions[i];
}
double TOF_reconstructor::GetScintPositionY(int i){
	return y_positions[i];
}
double TOF_reconstructor::Get_scint_xmax(int i){
	return x_dimensions[i]/2;
}
double TOF_reconstructor::Get_scint_xmin(int i){
	return -x_dimensions[i]/2;
}
double TOF_reconstructor::Get_ymax(){
	double scintillator_block_halfheight = 16.25 / 2.0;
	return y_positions[0] + scintillator_block_halfheight;
}
double TOF_reconstructor::Get_ymin(){
	double scintillator_block_halfheight = 16.25 / 2.0;
	return y_positions[7] - scintillator_block_halfheight;
}
vector<T5_hit> TOF_reconstructor::Return_position(const vector<int>& hit_mpmt_ids,
						  const vector<int>& hit_pmt_ids, 
						  const vector<double>& hit_pmt_times) {
	bool verbose = false;

	vector<T5_hit> all_hits;

	pair<double,double> position;
	set<int> hit_IDs;
	vector<vector<double>> TOF_times(16);	// Create 16 vectors where to store all the times detected by T5 SiPMs, some will be empty, some will have multiple hits, ideally paired -- later will check for hits in paired detectors
	bool pair_hit = false;
	bool valid_hit = false;
	double trigger_time = 0;
	for (int i = 0; i < hit_pmt_ids.size(); i++){
		if (hit_mpmt_ids.at(i) == 132 && hit_pmt_ids.at(i) == 19){
			trigger_time = hit_pmt_times.at(i);	// Take the first trigger event (should be only one, hopefully), and set it as trigger time, then break the for loop
			break;
		}
	}
	for (int i = 0; i < hit_mpmt_ids.size(); i++){
		auto mPMT_id = hit_mpmt_ids.at(i);	
		if (mPMT_id != 132) continue;
		auto PMT_id = hit_pmt_ids.at(i);
		if (! TOF_IDs.count(PMT_id)) continue;	// skip pmt IDs which are not T5 IDs
		const auto SiPM = distance(TOF_IDs.begin(), TOF_IDs.find(PMT_id));	// find the SiPM number -- different to channel number due to faulty channels on mpmt card 132
		TOF_times[SiPM].push_back(hit_pmt_times.at(i) - trigger_time);	// correct the measured T5 sipm time by the trigger time
		hit_IDs.insert(PMT_id);	// add the saved time ID to a set, to check for pair hits later
	}
	for (const auto& pair : corresponding_pairs){
		if (hit_IDs.count(pair.first) && hit_IDs.count(pair.second)){
			pair_hit = true;	// Loop over the saved pmt IDs, check if they are in a set of corresponding pairs, if not set bool to false
			break;
		} 
	}
	if(pair_hit == false){
		if (verbose)cout << "Pair of SiPMs was not hit, no reconstruction possible" << endl; // if no pair hit found -- return error value
		return {{false, {-999, -999}, -999, -999}};
	}
	else if(pair_hit == true){
		for (int i = 0; i < 8; i++){	// Loop over all saved times in the corresponding vectors, compare all the times
			for (const auto& sipm_time_a : TOF_times[i]){
				if (sipm_time_a > -140 || sipm_time_a < -170) {
					if (verbose) cout << "WARNING: SiPM time is out of expected hit times" << endl; 
					// warning that the detected time is outside the expected event peak, as detected in low energy events -- not tested in high energy events
				}
				for (const auto& sipm_time_b : TOF_times[i+8]){
					if (sipm_time_b > -140 || sipm_time_b < -170) {
						if (verbose) cout << "WARNING: SiPM time is out of expected hit times" << endl;
					}
					auto time_diff = sipm_time_a - sipm_time_b;
					double avg_time = (sipm_time_a + sipm_time_b) / 2;
					valid_hit = true;
					position.first = (time_diff - scintillator_bias.at(i)) * _v_eff / 2.0;
					position.second = y_positions[i];
					
					if (abs(position.first) > x_dimensions[i]/2){
						if (verbose)cout << "Reconstructed position is out of bounds" << endl;
						// If the reconstruction is outside of the corresponding scintillator bounds, set hit to false, return error value
						position.first = -999;
						valid_hit = false;
					}
					// verbose checks
					if (verbose) cout << "X coordinate is:\t" << position.first;
					
					//calculate the position uncertainty
					double uncertainty = sqrt(pow(v_eff_uncertainty, 2) * pow(time_diff/2, 2) + pow(sigma_sipm_i[i], 2) * pow(_v_eff/2, 2));
					if (verbose) cout << "\tX uncertainty is:\t" << uncertainty << endl;

					// save the hit to the vector of all hits in the event
					all_hits.push_back({valid_hit, position, uncertainty, avg_time});
				}
			}
		}

	}
	if (all_hits.size() != 1 && verbose)cout << "WARNING: More than 1 scintillator was hit" << endl;


	return all_hits;
}
bool TOF_reconstructor::HasMultiHits(const std::vector<T5_hit>& hits) {
    // If the vector is too small, we can return false immediately
    if (hits.size() < 2) return false;

    int valid_count = 0;
    for (const auto& hit : hits) {
        if (hit.valid_hit) {
            valid_count++;
        }
        // "Smart" Optimization: Stop as soon as we confirm > 1
        if (valid_count > 1) {
            return true;
        }
    }
    return false;
}
bool TOF_reconstructor::IsEventValid(const std::vector<T5_hit>& hits) {
    // If the vector is too small, we can return false immediately
    if (hits.size() < 1) return false;

    for (const auto& hit : hits) {
        if (hit.valid_hit) {
		return true;
        }
        // "Smart" Optimization: Stop as soon as we confirm > 1
    }
    return false;
}

