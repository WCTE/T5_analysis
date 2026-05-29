# T5 detector analysis

This project reconstructs T5 detector hits from WCTE input ROOT files and writes a simplified per-event ROOT tree for downstream analysis.

## Dependencies

The code is built against ROOT 6 (the Makefile uses `root-config`).

Typical setup:

- ROOT 6.x
- C++ compiler with C++17 support
- GNU Make

## Build

From the project root:

```sh
make
```

This produces the executable `analyze_T5`.

## Usage

The current executable accepts the following command-line options:

```sh
./analyze_T5 -r <run_number> -i <input_file.root> [-o <output_dir>] [-d]
```

Options:

- `-r <run_number>`: required run number used for the analysis metadata.
- `-i <input_file.root>`: input ROOT file to analyze. You can pass this option multiple times.
- `-o <output_dir>`: output directory prefix used to write one reconstructed file per input file.
- `-d`: enable debug mode (currently limits the event loop to the first 5000 entries).

Example:

```sh
./analyze_T5 -r 1361 -i /path/to/run_1361.root -o output/
```

The program recognizes the following input trees:

- `ProcessedWaveforms`
- `WCTEReadoutWindows`

If neither tree is found, the program stops with an error.

## What the program does

For each input file, the analysis:

1. Opens the input ROOT file.
2. Detects the appropriate TTree.
3. Reconstructs T5 hits using the T5 SiPM timing information.
4. Writes a new ROOT file containing a simplified event tree.

## Output tree structure

The output ROOT file contains one tree named `T5_Events` with the following branches:

| Branch name | Type | Description |
| :--- | :---: | :--- |
| `event_nr` | `Int_t` | Event number from the input tree |
| `T5_hit_mask` | `Bool_t` | Flag indicating whether the event is considered problematic / not clean |
| `T5_n_main_bunch_particles` | `Int_t` | Number of reconstructed main-bunch hits written for this event |
| `T5_hit_time` | `Double_t` | Time of the primary reconstructed hit |
| `T5_hit_charge` | `Double_t` | Charge of the primary reconstructed hit |
| `T5_hit_pos_x` | `Double_t` | Reconstructed x position of the primary hit [mm] |
| `T5_hit_pos_y` | `Double_t` | Reconstructed y position of the primary hit [mm] |
| `T5_additional_hit_pos_x` | `vector<double>` | x positions of additional reconstructed hits |
| `T5_additional_hit_pos_y` | `vector<double>` | y positions of additional reconstructed hits |
| `T5_additional_hit_time` | `vector<double>` | Times of additional reconstructed hits |
| `T5_additional_hit_charge` | `vector<double>` | Charges of additional reconstructed hits |

## Notes

- The program currently writes one output file per input file, using the input file basename with `_T5.root` appended.
- The output directory path passed to `-o` is used as the base location for those generated files.
- The reconstruction logic is implemented mainly in `return_TOF_position.cpp` and uses the detector geometry constants defined in `return_TOF_position.h`.
