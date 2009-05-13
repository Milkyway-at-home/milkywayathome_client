#ifndef COORDS_H
#define COORDS_H

#define rmat00 -0.054875539726
#define rmat01 -0.873437108010
#define rmat02 -0.483834985808
#define rmat10  0.494109453312
#define rmat11 -0.444829589425
#define rmat12  0.746982251810
#define rmat20 -0.867666135858
#define rmat21 -0.198076386122
#define rmat22  0.455983795705


#define stripe_separation 2.5
#define survey_center_dec 32.5
#define survey_center_ra 185.0

#define F_A_NODE 	95.0	//survey_center_ra - 90.0
#define F_A_NODE_RAD	(float)1.65806281566619873046875

#define D_A_NODE	95.0
#define D_A_NODE_RAD	D_A_NODE * D_DEG2RAD

double d_get_incl(int wedge);
double d_get_incl_rad(int wedge);

float f_get_incl(int wedge);
float f_get_incl_rad(int wedge);

void d_lbr2xyz(const double *lbr, double *xyz);

#endif
