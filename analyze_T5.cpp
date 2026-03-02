#include <vector>
#include <iostream>
#include <string>
#include <iomanip>

#include <ROOT/RVec.hxx>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TF2.h>

#include <nlohmann/json.hpp>

#include "./utils.h"
#include "./return_TOF_position.h"
#include "buffer.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;

using namespace ROOT;

int main(int argc, char** argv){
	
	if (argc < 2){
		cerr << "One argument expected, but was not provided. Usage: ./analyze_T5 <run_number>" << endl;
		return -1;
	}
//	TApplication app("app", &argc, argv);

	string arg = argv[1];
	auto run_number = std::stoi(arg);
	TString filename = "WCTE_data/charged_particle/WCTE_offline_R" + arg + "S0_VME_matched.root";
	auto file = TFile::Open(filename, "READ");
	
	if (!file || file-> IsZombie()){
		cerr << "ERROR, file did not open" << endl;
		return -1;
	}

	auto tree = file->Get<TTree>("WCTEReadoutWindows");

	nlohmann::json

	vector<float>* arr_bm_times = nullptr;
	vector<float>* arr_bm_charges = nullptr;
	vector<int>* arr_bm_time_ids = nullptr;
	vector<int>* arr_bm_charge_ids = nullptr;

	vector<double>* arr_pmt_times = nullptr;
	vector<int>* arr_pmt_ids = nullptr;
	vector<int>* arr_mpmt_ids = nullptr;

	tree-> SetBranchStatus("*", 0);
	tree-> SetBranchStatus("beamline_pmt_qdc_ids", 1);
	tree-> SetBranchStatus("beamline_pmt_tdc_ids", 1);
	tree-> SetBranchStatus("beamline_pmt_tdc_times", 1);
	tree-> SetBranchStatus("beamline_pmt_qdc_charges", 1);
	tree-> SetBranchStatus("hit_pmt_times", 1);
	tree-> SetBranchStatus("hit_mpmt_card_ids", 1);
	tree-> SetBranchStatus("hit_pmt_channel_ids", 1);

	tree-> SetBranchAddress("beamline_pmt_qdc_ids", &arr_bm_charge_ids);
	tree-> SetBranchAddress("beamline_pmt_tdc_ids", &arr_bm_time_ids);
	tree-> SetBranchAddress("beamline_pmt_tdc_times", &arr_bm_times);
	tree-> SetBranchAddress("beamline_pmt_qdc_charges", &arr_bm_charges);
	tree-> SetBranchAddress("hit_pmt_times", &arr_pmt_times);
	tree-> SetBranchAddress("hit_mpmt_card_ids", &arr_mpmt_ids);
	tree-> SetBranchAddress("hit_pmt_channel_ids", &arr_pmt_ids);


	Cuts cut;
	TOF_reconstructor recon;
	Histograms hists;
	setup_histograms(hists, recon);
	int n_pass_cut = 0; 
	int n_events = tree->GetEntries();
	int verb = 1000;
	for(int i = 0; i < n_events; i++){
		tree->GetEntry(i);
		//Print progress

		if (i % verb == 0) cout << "\rAnalyzed " << i << " of " << n_events << std::setprecision(2) << std::fixed <<
			" events (" << static_cast<float>(i)/n_events * 100 << " %)" << std::flush << endl; 

		RVecI bm_time_ids(arr_bm_time_ids->data(), arr_bm_time_ids->size());
		RVecI bm_charge_ids(arr_bm_charge_ids->data(), arr_bm_charge_ids->size());
		RVecF bm_times(arr_bm_times->data(), arr_bm_times->size());
		RVecF bm_charges(arr_bm_charges->data(), arr_bm_charges->size());

		RVecD pmt_times(arr_pmt_times->data(), arr_pmt_times->size());
		RVecI pmt_ids(arr_pmt_ids->data(), arr_pmt_ids->size());
		RVecI mpmt_ids(arr_mpmt_ids->data(), arr_mpmt_ids->size());

		
		if (!cut.hit_T0_T1(bm_time_ids, bm_charge_ids) ||
		    !cut.hit_T4(bm_charge_ids, bm_charges) ||
		    !cut.did_not_hit_HC(bm_charge_ids, bm_charges) ||
		    !cut.hit_T5(mpmt_ids, pmt_ids))
			continue;

		n_pass_cut++;

		auto mask_T5_board = (mpmt_ids == cut.get_T5_board());
		auto T5_board_ids = pmt_ids[mask_T5_board];
		auto T5_board_times = pmt_times[mask_T5_board];

		auto positions = recon.Return_position(*arr_mpmt_ids, *arr_pmt_ids, *arr_pmt_times);

		for (const auto& hit : positions){
			if (!hit.valid_hit) continue;
			hists.fill("positions", hit.position.first, hit.position.second);
		}

	}
	auto hist = hists.get_histogram_2D("positions");
	TF2* gaus_2D = new TF2("gaus_2D", "bigaus", recon.Get_scint_xmin(3), recon.Get_scint_xmax(3), recon.Get_ymin(), recon.Get_ymax());
	gaus_2D->SetParameters(130, 0, 40, 0, 40, 0);
	hist->Fit(gaus_2D, "R");

	gaus_2D = (TF2*)hist->GetFunction("gaus_2D");

	double volume = gaus_2D->GetParameter(0);
	double sig_x  = gaus_2D->GetParameter(2);
	double sig_y  = gaus_2D->GetParameter(4);
	double rho    = gaus_2D->GetParameter(5); // Correlation factor

	// 3. Calculate the TRUE mathematical peak height of the bigaus function
	double denominator = 2.0 * TMath::Pi() * sig_x * sig_y * std::sqrt(1.0 - rho*rho);
	double peak_amplitude = volume / denominator;

	// 4. Define your contour levels! 
	// 1-sigma drops to e^(-0.5)
	// 2-sigma drops to e^(-2.0)
	// 3-sigma drops to e^(-4.5)

	// Let's draw all 3 levels to make it look incredibly professional:
	double contours[3];
	contours[0] = peak_amplitude * std::exp(-4.5); // 3-sigma (widest, lowest)
	contours[1] = peak_amplitude * std::exp(-2.0); // 2-sigma
	contours[2] = peak_amplitude * std::exp(-0.5); // 1-sigma (tightest, highest)

	// 5. Apply the contours to your TF2
	// The arguments are: (number_of_levels, array_of_levels)
	gaus_2D->SetContour(3, contours);

	// 6. Make the contour lines stand out against the color map
	gaus_2D->SetLineColor(kRed);
	gaus_2D->SetLineWidth(2);
	gaus_2D->SetLineStyle(1); // Solid lines

	for (int i = 0; i < 8; i++){
		TString h_name = "positions_" + std::to_string(i); 
		hists. hist_projectX("positions", h_name.Data(), i+1, i+1);
		hists. get_histogram(h_name.Data())->Fit("gaus", "QR", "", recon.Get_scint_xmin(i), recon. Get_scint_xmax(i));
	}
	TString plots_directory = "plots/Run_" + std::to_string(run_number);
	gSystem->Exec("mkdir -p " + plots_directory);
	gSystem->cd(plots_directory);
	hists.print_all();

	cout << endl << n_pass_cut << " events out of " << n_events << " passed cuts" << endl;


	//	app.Run();

	return 0;
}
