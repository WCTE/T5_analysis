#include "buffer.h"
#include "return_TOF_position.h"

void setup_histograms(Histograms &hists, TOF_reconstructor& recon){
	int n_scints = 8;
	hists.book2D("positions", "Reconstructed T5 positions;x[mm];y[mm];count",
			50,
			-recon. GetScintDimensionX(3)/2,
			recon. GetScintDimensionX(3)/2,
			8,
			recon. GetScintPositionY(7) - 16.25/2.0,
			recon. GetScintPositionY(0) + 16.25/2.0);
	for(int i = 0; i < n_scints; i++){
		hists. book1D(Form("positions_%i", i), Form("Reconstructed positions in scintillator %i;x[mm];count", i),
			      50, -recon. GetScintDimensionX(i)/2, recon. GetScintDimensionX(i)/2);
	}
}
