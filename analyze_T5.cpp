#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

#include <ROOT/RVec.hxx>
#include <TApplication.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// #include <nlohmann/json.hpp>

#include "return_TOF_position.h"
#include "utils.h"
// #include "buffer.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace ROOT;

namespace {

void write_mask_summary_csv(const std::string &summary_path, int run_number,
                            const std::string &input_file,
                            long long total_events, long long good_events,
                            long long suspicious_events) {
    const bool file_exists = std::ifstream(summary_path).good();
    std::ofstream summary_file(summary_path,
                               std::ios::out | std::ios::app);
    if (!summary_file.is_open()) {
        std::cerr << "ERROR: Could not open summary CSV: " << summary_path
                  << std::endl;
        return;
    }

    if (!file_exists) {
        summary_file
            << "run_number,input_file,total_events,good_events,good_fraction_pct,"
               "events_with_multiple_main_window_hits_on_same_scintillator\n";
    }

    const double good_fraction_pct =
        total_events > 0 ? (100.0 * good_events) / total_events : 0.0;

    summary_file << run_number << ',' << input_file << ',' << total_events
                 << ',' << good_events << ',' << std::fixed
                 << std::setprecision(3) << good_fraction_pct << ','
                 << suspicious_events << '\n';
    summary_file.close();
}

} // namespace

int main(int argc, char **argv) {

    int run_number = -1;
    vector<string> input_paths;
    TString output_path = "output.root";
    bool debug = false;

    int opt;
    while ((opt = getopt(argc, argv, "r:i:o:d")) != -1) {
        switch (opt) {
        case 'r':
            run_number = std::stoi(optarg);
            break;
        case 'i':
            input_paths.push_back(optarg);
            break;
        case 'o':
            output_path = optarg;
            break;
        case 'd':
            debug = true;
            break;
        default:
            cerr << "Usage: " << argv[0]
                 << " -r <run_number> [-i <input_file>] [-o "
                    "<output_file>] [-d]"
                 << endl;
            return -1;
        }
    }
    // debug print the input paths
    cout << "Input paths: " << endl;
    for (size_t input_idx = 0; input_idx < input_paths.size(); ++input_idx) {
        cout << input_paths[input_idx] << endl;
    }

    RUN_NUMBER = run_number;

    if (input_paths.empty()) {
        throw std::runtime_error("No input files provided");
    }

    for (size_t input_idx = 0; input_idx < input_paths.size(); ++input_idx) {
        TString filename = input_paths[input_idx];
        TString current_output_path = output_path;

        TString base_name = gSystem->BaseName(filename);
        base_name.ReplaceAll(".root", "_T5.root");

        if (!current_output_path.EndsWith("/")) {
            current_output_path += "/";
        }
        current_output_path = output_path + base_name;

        auto file = TFile::Open(filename, "READ");

        if (!file || file->IsZombie()) {
            cerr << "ERROR, file did not open: " << filename << endl;
            continue;
        }

        // std::ifstream
        // config_file("/home/frantisek/scripts/config.json");
        // nlohmann::json config = nlohmann::json:: parse(config_file);
        // auto& run_config = config[std::to_string(run_number)];
        // BEAM_MOMENTUM = run_config.value("Beam momentum (MeV/c)", 0);
        // cout << "Config file loaded, beam momentum is " <<
        // BEAM_MOMENTUM << "MeV/c" << endl ;

        bool hardware_processed_data = false;
        string tree_name = "";
        if (file->GetListOfKeys()->Contains("ProcessedWaveforms")) {
            tree_name = "ProcessedWaveforms";
            hardware_processed_data = true;
        } else if (file->GetListOfKeys()->Contains("WCTEReadoutWindows")) {
            tree_name = "WCTEReadoutWindows";
        } else {
            cerr << "ERROR: Input file does not contain recognized "
                    "tree "
                    "(WCTEReadoutWindows or ProcessedWaveforms)"
                 << endl;
            return -1;
        }
        TTreeReader tree(tree_name.c_str(), file);

        // vector<float>* arr_bm_times = nullptr;
        // vector<float>* arr_bm_charges = nullptr;
        // vector<int>* arr_bm_time_ids = nullptr;
        // vector<int>* arr_bm_charge_ids = nullptr;

        std::unique_ptr<TTreeReaderValue<vector<double>>> arr_pmt_times;
        std::unique_ptr<TTreeReaderValue<vector<float>>> arr_pmt_charges;
        std::unique_ptr<TTreeReaderValue<vector<int>>> arr_pmt_ids;
        std::unique_ptr<TTreeReaderValue<vector<int>>> arr_mpmt_ids;

        // needed for array input of the processed waveforms file

        const int MAX_HITS = 4000;

        std::unique_ptr<TTreeReaderValue<int>> n_pmt_times;
        std::unique_ptr<TTreeReaderValue<int>> n_pmt_charges;
        std::unique_ptr<TTreeReaderValue<int>> n_pmt_chans;
        std::unique_ptr<TTreeReaderValue<int>> n_pmt_cards;

        std::unique_ptr<TTreeReaderArray<double>> processed_hit_times;
        std::unique_ptr<TTreeReaderArray<float>> processed_hit_charges;
        std::unique_ptr<TTreeReaderArray<int>> processed_hit_cards;
        std::unique_ptr<TTreeReaderArray<int>> processed_hit_chans;

        if (hardware_processed_data) {
            n_pmt_times =
                std::make_unique<TTreeReaderValue<int>>(tree, "nhit_time");
            n_pmt_charges =
                std::make_unique<TTreeReaderValue<int>>(tree, "nhit_charge");
            n_pmt_chans =
                std::make_unique<TTreeReaderValue<int>>(tree, "nhit_chan");
            n_pmt_cards =
                std::make_unique<TTreeReaderValue<int>>(tree, "nhit_card");

            processed_hit_times =
                std::make_unique<TTreeReaderArray<double>>(tree, "hit_time");
            processed_hit_charges =
                std::make_unique<TTreeReaderArray<float>>(tree, "hit_charge");
            processed_hit_cards =
                std::make_unique<TTreeReaderArray<int>>(tree, "hit_card");
            processed_hit_chans =
                std::make_unique<TTreeReaderArray<int>>(tree, "hit_chan");

        } else {
            arr_pmt_times = std::make_unique<TTreeReaderValue<vector<double>>>(
                tree, "hit_pmt_times");
            arr_pmt_charges =
                std::make_unique<TTreeReaderValue<vector<float>>>(
                    tree, "hit_pmt_charges");
            arr_pmt_ids = std::make_unique<TTreeReaderValue<vector<int>>>(
                tree, "hit_pmt_channel_ids");
            arr_mpmt_ids = std::make_unique<TTreeReaderValue<vector<int>>>(
                tree, "hit_mpmt_card_ids");
        }

        // tree-> SetBranchAddress("beamline_pmt_qdc_ids",
        // &arr_bm_charge_ids); tree->
        // SetBranchAddress("beamline_pmt_tdc_ids", &arr_bm_time_ids);
        // tree-> SetBranchAddress("beamline_pmt_tdc_times",
        // &arr_bm_times); tree->
        // SetBranchAddress("beamline_pmt_qdc_charges",
        // &arr_bm_charges);

        Cuts cut;
        TOF_reconstructor recon;
        // Histograms hists;
        // setup_histograms(hists, recon);
        int n_pass_cut = 0;
        int n_T5_valid_events = 0;
        long long total_events = tree.GetEntries();
        long long good_events = 0;
        long long suspicious_main_window_events = 0;
        auto n_events = total_events;
        if (debug) {
            cout << "Debug mode enabled: limiting to 5000 events" << endl;
            n_events = std::min(n_events, 5000LL);
        }
        int verb = 1000;
        int n_events_with_multiple_valid_hits = 0;
        int n_events_with_valid_hits_in_expected_window = 0;
        int n_events_with_multiple_scint_hits = 0;
        int n_events_with_multiple_valid_hits_had_one_in_expected_window = 0;
        int n_invalid_hits = 0;
        int n_events_out_of_bounds = 0;

        vector<event_T5_detection> all_T5_hits;

        while (tree.Next()) {
            auto i = tree.GetCurrentEntry();

            if (i >= n_events)
                break;

            vector<double> vec_hit_time;
            vector<double> vec_hit_charge;
            vector<int> vec_hit_chan;
            vector<int> vec_hit_card;

            event_T5_detection detections;

            if (hardware_processed_data) {
                if (**n_pmt_times >= MAX_HITS) {
                    // continue and push back dummy
                    // detections
                    std::cerr << "ERROR: Too many hits in event "
                              << tree.GetCurrentEntry() << " nhits "
                              << **n_pmt_times << " (max " << MAX_HITS << ")"
                              << std::endl;
                    detections.event_nr = tree.GetCurrentEntry();
                    all_T5_hits.push_back(detections);
                    continue;
                }
                // change array output into std vectors
                vec_hit_time.assign((*processed_hit_times).begin(),
                                    (*processed_hit_times).end());
                vec_hit_charge.assign((*processed_hit_charges).begin(),
                                      (*processed_hit_charges).end());
                vec_hit_card.assign((*processed_hit_cards).begin(),
                                    (*processed_hit_cards).end());
                vec_hit_chan.assign((*processed_hit_chans).begin(),
                                    (*processed_hit_chans).end());
            } else {
                vec_hit_time.assign((**arr_pmt_times).begin(),
                                    (**arr_pmt_times).end());
                vec_hit_charge.assign((**arr_pmt_charges).begin(),
                                      (**arr_pmt_charges).end());
                vec_hit_card.assign((**arr_mpmt_ids).begin(),
                                    (**arr_mpmt_ids).end());
                vec_hit_chan.assign((**arr_pmt_ids).begin(),
                                    (**arr_pmt_ids).end());
            }

            if (i % verb == 0)
                cout << "\rAnalyzed " << i << " of " << n_events
                     << std::setprecision(2) << std::fixed << " events ("
                     << static_cast<float>(i) / n_events * 100 << " %)"
                     << std::flush << endl;

            // RVecI bm_time_ids(arr_bm_time_ids->data(),
            // arr_bm_time_ids->size()); RVecI
            // bm_charge_ids(arr_bm_charge_ids->data(),
            // arr_bm_charge_ids->size()); RVecF
            // bm_times(arr_bm_times->data(), arr_bm_times->size());
            // RVecF bm_charges(arr_bm_charges->data(),
            // arr_bm_charges->size());

            RVecD pmt_times(vec_hit_time.data(), vec_hit_time.size());
            RVecD pmt_charges(vec_hit_charge.data(), vec_hit_charge.size());
            RVecI pmt_ids(vec_hit_chan.data(), vec_hit_chan.size());
            RVecI mpmt_ids(vec_hit_card.data(), vec_hit_card.size());

            if (!cut.hit_T5(mpmt_ids, pmt_ids)) {
                detections.event_nr = i;
                all_T5_hits.push_back(detections);
                continue;
            }

            n_pass_cut++;

            // auto mask_T5_board = (mpmt_ids == cut.get_T5_board());
            // auto T5_board_ids = pmt_ids[mask_T5_board];
            // auto T5_board_times = pmt_times[mask_T5_board];

            detections = recon.Return_position(i, vec_hit_card, vec_hit_chan,
                                               vec_hit_time, vec_hit_charge);

            std::array<int, 8> main_window_hits_per_scintillator = {};
            for (const auto &hit : detections.T5_hits) {
                if (hit.is_in_time_window && hit.is_valid_hit &&
                    hit.scintillator_id >= 0 &&
                    hit.scintillator_id < static_cast<int>(main_window_hits_per_scintillator.size())) {
                    ++main_window_hits_per_scintillator[hit.scintillator_id];
                }
            }
            const bool has_multiple_hits_on_same_scintillator =
                std::any_of(main_window_hits_per_scintillator.begin(),
                            main_window_hits_per_scintillator.end(),
                            [](int count) { return count > 1; });
            if (has_multiple_hits_on_same_scintillator) {
                ++suspicious_main_window_events;
                std::cerr << "WARNING: event " << i
                          << " has multiple main-window hits on the same "
                             "scintillator; review the reconstruction output."
                          << std::endl;
            }

            if (detections.IsClean && detections.HasInTimeWindow) {
                ++good_events;
            }

            if (detections.HasValidHit)
                n_T5_valid_events++;
            if (detections.HasMultipleValidHits) {
                n_events_with_multiple_valid_hits++;
                if (detections.HasInTimeWindow)
                    n_events_with_multiple_valid_hits_had_one_in_expected_window++;
                if (detections.HasMultipleScintillatorsHit)
                    n_events_with_multiple_scint_hits++;
            }
            if (detections.HasInTimeWindow)
                n_events_with_valid_hits_in_expected_window++;
            if (detections.HasHit && !detections.HasValidHit) {
                n_invalid_hits++;
            }
            if (detections.HasOutOfBounds)
                n_events_out_of_bounds++;

            all_T5_hits.push_back(detections);
        }

        file->Close();

        // ROOT output-file writing is disabled for file-by-file Python-driven
        // analysis runs that only need the summary statistics.
        // TFile *output_file = TFile::Open(current_output_path, "RECREATE");
        // if (!output_file || output_file->IsZombie()) {
        //     cerr << "ERROR: Did not open output file: " << current_output_path
        //          << endl;
        //     continue;
        // }

        // output_file->cd();
        TH1D *h_main_window_hits_per_scintillator =
            new TH1D("T5_main_window_hits_per_scintillator",
                     "Main-window hits per scintillator;Scintillator ID;Count",
                     8, -0.5, 7.5);
        TH1D *h_suspicious_events =
            new TH1D("T5_suspicious_multi_hits_same_scintillator",
                     "Events with >1 main-window hit on the same scintillator;"
                     "Event count;Entries",
                     2, -0.5, 1.5);
        h_suspicious_events->SetBinContent(2, suspicious_main_window_events);
        h_suspicious_events->SetBinContent(1, total_events - suspicious_main_window_events);

        for (const auto &event : all_T5_hits) {
            for (const auto &hit : event.T5_hits) {
                if (hit.is_in_time_window && hit.is_valid_hit &&
                    hit.scintillator_id >= 0 &&
                    hit.scintillator_id < 8) {
                    h_main_window_hits_per_scintillator->Fill(hit.scintillator_id);
                }
            }
        }

        // h_main_window_hits_per_scintillator->Write();
        // h_suspicious_events->Write();

        // TTree *output_tree = new TTree("T5_Events", "Reconstructed T5 events");
        // --- Single-value branches (per event) ---
        int b_n_main_particles = 0;
        int b_event_nr = 0;
        double b_main_hit_time = 0;
        double b_main_hit_charge = 0;
        double b_main_position_x = 0;
        double b_main_position_y = 0;
        bool b_T5HitMask;

        // output_tree->Branch("event_nr", &b_event_nr, "event_nr/I");
        // output_tree->Branch("T5_hit_mask", &b_T5HitMask, "T5_hit_mask/O");
        // output_tree->Branch("T5_n_main_bunch_particles", &b_n_main_particles,
        //                     "T5_n_main_bunch_particles/I");
        // output_tree->Branch("T5_hit_time", &b_main_hit_time, "T5_hit_time/D");
        // output_tree->Branch("T5_hit_charge", &b_main_hit_charge,
        //                     "T5_hit_charge/D");
        // output_tree->Branch("T5_hit_pos_x", &b_main_position_x,
        //                     "T5_hit_pos_x/D");
        // output_tree->Branch("T5_hit_pos_y", &b_main_position_y,
        //                     "T5_hit_pos_y/D");

        // --- Vector branches (multiple hits per event) ---
        // Additional hits -- Any hit other than the main hit
        std::vector<double> *b_additional_hit_pos_x = new std::vector<double>();
        std::vector<double> *b_additional_hit_pos_y = new std::vector<double>();
        std::vector<double> *b_additional_hit_time = new std::vector<double>();
        std::vector<double> *b_additional_hit_charge =
            new std::vector<double>();
        std::vector<int> *b_main_hit_scintillator_id = new std::vector<int>();
        std::vector<int> *b_additional_hit_scintillator_id =
            new std::vector<int>();

        // output_tree->Branch("T5_additional_hit_pos_x", &b_additional_hit_pos_x);
        // output_tree->Branch("T5_additional_hit_pos_y", &b_additional_hit_pos_y);
        // output_tree->Branch("T5_additional_hit_time", &b_additional_hit_time);
        // output_tree->Branch("T5_additional_hit_charge",
        //                     &b_additional_hit_charge);
        // output_tree->Branch("T5_hit_scintillator_id",
        //                     &b_main_hit_scintillator_id);
        // output_tree->Branch("T5_additional_hit_scintillator_id",
        //                     &b_additional_hit_scintillator_id);

        for (const auto &event : all_T5_hits) {
            b_event_nr = event.event_nr;
            b_T5HitMask = false;

            if (!event.IsClean || !event.HasInTimeWindow)
                b_T5HitMask = true;

            b_n_main_particles = 0;
            b_main_hit_charge = 0;
            
            b_main_hit_time = 0;
            b_main_position_x = 0;
            b_main_position_y = 0;

            b_additional_hit_time->clear();
            b_additional_hit_charge->clear();
            b_additional_hit_pos_x->clear();
            b_additional_hit_pos_y->clear();
            b_main_hit_scintillator_id->clear();
            b_additional_hit_scintillator_id->clear();

            if (!event.HasValidHit) {
                // output_tree->Fill();
                continue;
            }
            bool main_hit_happened = false;
            for (const auto &hit : event.T5_hits) {
                if (hit.quality == HitQuality::AccidentalCoincidence)
                    continue;
                if (hit.is_in_time_window && !main_hit_happened) {
                    b_main_hit_time = hit.raw_time;
                    b_main_hit_charge = hit.hit_charge;
                    b_main_position_x = hit.position_x;
                    b_main_position_y = hit.position_y;
                    b_main_hit_scintillator_id->push_back(hit.scintillator_id);
                    b_n_main_particles++;
                    main_hit_happened = true;
                } else if (hit.is_in_time_window) {
                    b_n_main_particles++;
                    b_additional_hit_time->push_back(hit.raw_time);
                    b_additional_hit_charge->push_back(hit.hit_charge);
                    b_additional_hit_pos_x->push_back(hit.position_x);
                    b_additional_hit_pos_y->push_back(hit.position_y);
                    b_additional_hit_scintillator_id->push_back(hit.scintillator_id);
                }

                else {
                    b_additional_hit_time->push_back(hit.raw_time);
                    b_additional_hit_charge->push_back(hit.hit_charge);
                    b_additional_hit_pos_x->push_back(hit.position_x);
                    b_additional_hit_pos_y->push_back(hit.position_y);
                    b_additional_hit_scintillator_id->push_back(hit.scintillator_id);
                }
            }

            // output_tree->Fill();
        }
        // output_tree->Write();

        delete b_additional_hit_charge;
        delete b_additional_hit_time;
        delete b_additional_hit_pos_x;
        delete b_additional_hit_pos_y;
        delete b_main_hit_scintillator_id;
        delete b_additional_hit_scintillator_id;

        write_mask_summary_csv((std::string(current_output_path.Data()) +
                                    ".summary.csv"),
                               run_number, filename.Data(), total_events,
                               good_events, suspicious_main_window_events);

        // output_file->Close();
    } // end loop over input files

    return 0;
}
