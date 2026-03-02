#ifndef RETURN_TOF_POSITION_H
#define RETURN_TOF_POSITION_H

#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <math.h>


struct T5_hit{
	bool valid_hit = false; // Tells if the hit was out of bounds or was reconstructed inside
	std::pair<double, double> position; //X and Y coordinates of the T5 hit
	double uncertainty; //Position uncertainty, assumed that v_eff and sigma_sipm are uncorrelated
	double hit_time; // Time of the detection, calculated as average time of the two SiPMs
	
};

// The class TOF reconstructor has several useful functions. The main one is the Return_position, which takes the data, analyzes it and returns a vcector of the T5 hit structure that contains the validity of the hit, the position, uncertainty of the hit and the hit time of all hits. It also has the method HasMultiHits, which checks if the event has munltiple valid hits (i.e. several hits in the expected timeframe of time - trigger between -140 and -170 ns). This method looks at the vector of T5 hits returned by Return_position. Then there is IsEventValid method, which checks, if there are any valid hits inside of the event. 

class TOF_reconstructor{
	public:
		explicit TOF_reconstructor(double v_eff);
		TOF_reconstructor();
		void setVeff(double v);
		double GetVeff() const;
		double GetScintDimensionX(int i);
		double GetScintPositionY(int i);
		double Get_scint_xmin(int i);
		double Get_scint_xmax(int i);
		double Get_ymax();
		double Get_ymin();
		std::vector<T5_hit> Return_position(const std::vector<int>& hit_mpmt_ids,
						    const std::vector<int>& hit_pmt_ids, 
						    const std::vector<double>& hit_pmt_times);
		static bool HasMultiHits(const std::vector<T5_hit>& hits);
		static bool IsEventValid(const std::vector<T5_hit>& hits);
	private:
		double _v_eff;
		const std::vector<double> scintillator_bias;
		const std::vector<double> y_positions;
		const std::vector<double> x_dimensions;
		const double v_eff_uncertainty;
		const std::vector<double> sigma_sipm_i;
		const std::vector<double> sigma_sipm_i_uncertainties;
		const std::set<std::pair<int, int>> corresponding_pairs;
    		const std::set<int> TOF_IDs;


};



#endif 
