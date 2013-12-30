




#ifndef _LABELING_
#define _LABELING_



const int MAXLABELCONNECTED = 2000;
const int BG = 0;

void SetImage(int dim_x, int dim_y, double* space_density_1);
/// void SetSpaceDensity2D(Int_t i, Int_t j, Double_t value);
/// Double_t GetSpaceDensity2D(Int_t i, Int_t j);
/// Int_t SecondStep(unsigned int last_i, unsigned int NewLabel, int* C );
int LabelingDuePassi();
int LabelingDuePassi8();	/// Unused for now

#endif
