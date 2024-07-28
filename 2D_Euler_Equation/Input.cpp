/* Explicit Two-Derivative Two-Step Runge-Kutta Methods
 * Copyright (C) 2024
 * Author QinXueyu
 */
// #pragma once
#include "xBurgers.hpp"
using namespace std;

void CADSolver::InputControl_xy(void)
{
	unsigned short nSkip = 8, iSkip = 0;
	mydouble AuxVal_a;
	NumberThreads = 1; // OpenMP

	nCell = 150;
	nCellYo = nCell;
	TimeScheme = 231; // case 223: 	TDTSRK223();	// case 224:    TDTSRK224();	// case 235:    TDTSRK235();
	CFL = 1.0;
	StopTime = 10.;	  // 总时间
	StepMax = 500000; // 总步数
	Case_name = 83;	  // sod:1; 2D_Remain:11;  Blast_wave: 2;  order_x: 80; order_xy: 82; vortex: 83;
					// DoubleMach 22; High_jet 23;
					// Couette_flow 26;  Viscous shock tube: 24;  Kelvin_Helmholtz 25; Laminar 27;
	BoundaryScheme = 0; // periodic 0; Expand 1; reflection_x 2;  Viscous shock tube: 24;
	SpaceScheme = 72;	// up1:1; WENO5-JS 50;  WENO5-Z 52; TENO5 53; TENO5-A 54; WENO7-JS 70;  WENO7-Z 72;
	SpaceScheme_cd = 6; // up1:1; WENO5-JS负50; WENO5-JS正 51;  center 4-th,  6-th
	char_scheme = 0;
	Local_Lax_Friedrichs = 1; // Local 1; Global 0;
	Case_1D_or_2D = 2;		  // 1D 1; 2D 2;
	Hy_flag = 0;			  // 混合格式 1 开启

	//___________________positivity_____________________
	positivity_flag = 0; // 保正性
	positivity_flag_LW = 0;
	//___________________Viscous_____________________
	Viscous_flag = 0; // Non_viscous: 0 ;   Viscous:1
	Viscous_flag_dt = 0;
	Viscous_cd = 6;
	LW_Viscous_cd = 4;
	//________________________________________
	nb = 10; // 扩展单元数
	nCell_2L = nCell + 2 * nb, nCellYo_2L = nCellYo + 2 * nb;
	nCell_L = nCell + nb, nCellYo_L = nCellYo + nb;
	nVar = 4;
	nDegree = 3;
	ska = 1;
	show_n = 50;
	T_max_u = -10.;
	T_max_v = -10.;
	cout << "---------------- Input -----------------" << endl;
	Intial_zone();
	Intial_date();
}

void CADSolver::Intial_zone(void)
{
	switch (Case_name)
	{

	case 82: // x,y  order
		StopTime = 2.;
		nx_begin = 0., nx_end = 2.;
		ny_begin = 0., ny_end = 2.;
		BoundaryScheme = 0;
		T_max_u = 1.0;
		T_max_v = 1.0;
		break;

	case 83: // 2D vortex
		nx_begin = -5., nx_end = 5.;
		ny_begin = -5., ny_end = 5.;
		BoundaryScheme = 0;
		StopTime = 10.;
		Re = 100000.;
		Pr = 0.25;
		T_max_u = 3.0;
		T_max_v = 3.0;
		break;

	case 1:
		nx_begin = 0., nx_end = 1.;
		ny_begin = 0., ny_end = 1.;
		BoundaryScheme = 1;
		T_max_u = 2.0;
		T_max_v = 0.0;
		StopTime = 0.2;
		break;

	case 11:
		nx_begin = 0., nx_end = 1.;
		ny_begin = 0., ny_end = 1.;
		BoundaryScheme = 1;
		StopTime = 0.3;
		T_max_u = 2.0;
		T_max_v = 2.0;
		break;

	case 7: // case_Shock_vortex_2D
		nx_begin = 0., nx_end = 2.;
		ny_begin = 0., ny_end = 1.;
		BoundaryScheme = 1;
		StopTime = 0.75;
		T_max_u = 2.5;
		T_max_v = 1.5;
		break;

	case 2:
		nx_begin = 0., nx_end = 10.;
		ny_begin = 0., ny_end = 1.;
		BoundaryScheme = 2;
		StopTime = 0.38;
		break;

	case 80: // x order
		StopTime = 2.;
		nx_begin = 0., nx_end = 2.;
		ny_begin = 0., ny_end = 2.;
		BoundaryScheme = 0;
		T_max_u = 2.0;
		break;

	case 81: // y order
		StopTime = 2.;
		nx_begin = 0., nx_end = 2.;
		ny_begin = 0., ny_end = 2.;
		BoundaryScheme = 0;
		break;

	default:
		cout << "Intial_zone is not avaliable " << endl;
		break;
	}
}

void CADSolver::Intial_date(void)
{
	Gamma = 1.4, S_ct = 1.0e-10;

	if (Case_name == 23)
		Gamma = 5. / 3., S_ct = 1.0e-2; // High_jet case   {D_log}=Log({D})

	if (Case_name == 2)
		S_ct = 1.0e-5; // Blast_wave

	if (Case_name == 11)
		S_ct = 1.0e-14;

	if (Case_name == 1)
		S_ct = 1.0e-10;

	// if (Case_name == 24)
	// 	S_ct = 1.0e-7; // Viscous shock tube: 24;  S_ct = 1.0e-1

	if (Case_name == 27)
		S_ct = 1.0e-1; // Laminar: 27;  S_ct = 1.0e-1

	Gamma_1 = Gamma - 1.0;
	C_v = 1. / (Gamma * Gamma_1);
	C_p = Gamma * C_v;
	PI_Number = 4.0 * atan(1.0);
	Ts = 288.15;
	Tc = 110.4;

	D_utx0 = new mydouble **[nCell_2L];
	D_uty0 = new mydouble **[nCell_2L];

	D_ut = new mydouble *[nCell_2L];
	D_vt = new mydouble *[nCell_2L];
	D_Tt = new mydouble *[nCell_2L];
	D_utt = new mydouble *[nCell_2L];
	D_vtt = new mydouble *[nCell_2L];
	D_Ttt = new mydouble *[nCell_2L];

	D_pt = new mydouble *[nCell_2L];
	D_mu = new mydouble *[nCell_2L];
	D_ka = new mydouble *[nCell_2L];

	D_duvpc = new mydouble **[nCell_2L];
	D_du = new mydouble **[nCell_2L];
	D_du0 = new mydouble **[nCell_2L];
	D_du1 = new mydouble **[nCell_2L];
	D_ddu1 = new mydouble **[nCell_2L];
	D_ddu = new mydouble **[nCell_2L];
	D_ddu0 = new mydouble **[nCell_2L];
	D_du_up1 = new mydouble **[nCell_2L];
	D_un = new mydouble **[nCell_2L];
	D_un0 = new mydouble **[nCell_2L];
	D_un1 = new mydouble **[nCell_2L];
	D_ddu_up1 = new mydouble **[nCell_2L];
	D_dddu_up1 = new mydouble **[nCell_2L];
	D_dddu = new mydouble **[nCell_2L];
	flagh_hy = new mydouble **[nCell_2L];

	D_Coord = new mydouble **[nCell_2L];
	D_Exact = new mydouble **[nCell_2L];
	D_out = new mydouble **[nCell_2L];

	D_fxz = new mydouble **[nCell_2L];
	D_fxf = new mydouble **[nCell_2L];
	D_gx = new mydouble **[nCell_2L];
	D_gy = new mydouble **[nCell_2L];
	D_gx_t = new mydouble **[nCell_2L];
	D_gy_t = new mydouble **[nCell_2L];

	D_fxz_up1 = new mydouble **[nCell_2L];
	D_fyz_up1 = new mydouble **[nCell_2L];

	D_fx_char = new mydouble **[nCell_2L];
	D_fy_char = new mydouble **[nCell_2L];
	D_fx_up1_char = new mydouble **[nCell_2L];
	D_fy_up1_char = new mydouble **[nCell_2L];

	D_fyz = new mydouble **[nCell_2L];
	D_fyf = new mydouble **[nCell_2L];
	D_fxat = new mydouble ***[nCell_2L];
	D_fyat = new mydouble ***[nCell_2L];

	D_dun1 = new mydouble **[nCell_2L];
	D_ddun1 = new mydouble **[nCell_2L];
	D_duns1 = new mydouble **[nCell_2L];
	D_dduns1 = new mydouble **[nCell_2L];

	D_un2 = new mydouble **[nCell_2L];
	D_un3 = new mydouble **[nCell_2L];
	D_dun2 = new mydouble **[nCell_2L];
	D_ddun2 = new mydouble **[nCell_2L];

	// ***D_du2, ***D_du3,   ***D_du4;
	// 定义上面的变量
	D_du2 = new mydouble **[nCell_2L];
	D_du3 = new mydouble **[nCell_2L];
	D_du4 = new mydouble **[nCell_2L];

	for (long iDegree = 0; iDegree < nCell_2L; iDegree++)
	{

		D_ut[iDegree] = new mydouble[nCellYo_2L];
		D_vt[iDegree] = new mydouble[nCellYo_2L];
		D_Tt[iDegree] = new mydouble[nCellYo_2L];

		D_utt[iDegree] = new mydouble[nCellYo_2L];
		D_vtt[iDegree] = new mydouble[nCellYo_2L];
		D_Ttt[iDegree] = new mydouble[nCellYo_2L];

		D_pt[iDegree] = new mydouble[nCellYo_2L];
		D_mu[iDegree] = new mydouble[nCellYo_2L];
		D_ka[iDegree] = new mydouble[nCellYo_2L];

		D_dun1[iDegree] = new mydouble *[nCellYo_2L];
		D_ddun1[iDegree] = new mydouble *[nCellYo_2L];
		D_duns1[iDegree] = new mydouble *[nCellYo_2L];
		D_dduns1[iDegree] = new mydouble *[nCellYo_2L];

		D_un2[iDegree] = new mydouble *[nCellYo_2L];
		D_un3[iDegree] = new mydouble *[nCellYo_2L];
		D_dun2[iDegree] = new mydouble *[nCellYo_2L];
		D_ddun2[iDegree] = new mydouble *[nCellYo_2L];

		D_du2[iDegree] = new mydouble *[nCellYo_2L];
		D_du3[iDegree] = new mydouble *[nCellYo_2L];
		D_du4[iDegree] = new mydouble *[nCellYo_2L];

		D_utx0[iDegree] = new mydouble *[nCellYo_2L];
		D_uty0[iDegree] = new mydouble *[nCellYo_2L];
		D_duvpc[iDegree] = new mydouble *[nCellYo_2L];
		D_du[iDegree] = new mydouble *[nCellYo_2L];
		D_du0[iDegree] = new mydouble *[nCellYo_2L];
		D_du1[iDegree] = new mydouble *[nCellYo_2L];
		D_ddu1[iDegree] = new mydouble *[nCellYo_2L];
		D_ddu[iDegree] = new mydouble *[nCellYo_2L];
		D_ddu0[iDegree] = new mydouble *[nCellYo_2L];
		D_du_up1[iDegree] = new mydouble *[nCellYo_2L];
		D_ddu_up1[iDegree] = new mydouble *[nCellYo_2L];
		D_dddu_up1[iDegree] = new mydouble *[nCellYo_2L];
		D_dddu[iDegree] = new mydouble *[nCellYo_2L];

		flagh_hy[iDegree] = new mydouble *[nCellYo_2L];
		D_fx_char[iDegree] = new mydouble *[nCellYo_2L];
		D_fy_char[iDegree] = new mydouble *[nCellYo_2L];

		D_fx_up1_char[iDegree] = new mydouble *[nCellYo_2L];
		D_fy_up1_char[iDegree] = new mydouble *[nCellYo_2L];
		D_gx[iDegree] = new mydouble *[nCellYo_2L];
		D_gy[iDegree] = new mydouble *[nCellYo_2L];
		D_gx_t[iDegree] = new mydouble *[nCellYo_2L];
		D_gy_t[iDegree] = new mydouble *[nCellYo_2L];

		D_Coord[iDegree] = new mydouble *[nCellYo_2L];
		D_Exact[iDegree] = new mydouble *[nCellYo_2L];
		D_out[iDegree] = new mydouble *[nCellYo_2L];
		D_un[iDegree] = new mydouble *[nCellYo_2L];
		D_un0[iDegree] = new mydouble *[nCellYo_2L];
		D_un1[iDegree] = new mydouble *[nCellYo_2L];
		D_fxz[iDegree] = new mydouble *[nCellYo_2L];
		D_fxf[iDegree] = new mydouble *[nCellYo_2L];
		D_fyz[iDegree] = new mydouble *[nCellYo_2L];
		D_fyf[iDegree] = new mydouble *[nCellYo_2L];
		D_fxz_up1[iDegree] = new mydouble *[nCellYo_2L];
		D_fyz_up1[iDegree] = new mydouble *[nCellYo_2L];
		D_fxat[iDegree] = new mydouble **[nCellYo_2L];
		D_fyat[iDegree] = new mydouble **[nCellYo_2L];

		for (long jDegree = 0; jDegree < nCellYo_2L; jDegree++)
		{

			D_dun1[iDegree][jDegree] = new mydouble[nVar];
			D_ddun1[iDegree][jDegree] = new mydouble[nVar];
			D_duns1[iDegree][jDegree] = new mydouble[nVar];
			D_dduns1[iDegree][jDegree] = new mydouble[nVar];

			D_un2[iDegree][jDegree] = new mydouble[nVar];
			D_un3[iDegree][jDegree] = new mydouble[nVar];
			D_dun2[iDegree][jDegree] = new mydouble[nVar];
			D_ddun2[iDegree][jDegree] = new mydouble[nVar];

			D_utx0[iDegree][jDegree] = new mydouble[nVar];
			D_uty0[iDegree][jDegree] = new mydouble[nVar];
			D_du[iDegree][jDegree] = new mydouble[nVar];
			D_du0[iDegree][jDegree] = new mydouble[nVar];
			D_du1[iDegree][jDegree] = new mydouble[nVar];
			D_ddu1[iDegree][jDegree] = new mydouble[nVar];
			D_ddu_up1[iDegree][jDegree] = new mydouble[nVar];
			D_dddu_up1[iDegree][jDegree] = new mydouble[nVar];
			D_dddu[iDegree][jDegree] = new mydouble[nVar];

			D_du2[iDegree][jDegree] = new mydouble[nVar];
			D_du3[iDegree][jDegree] = new mydouble[nVar];
			D_du4[iDegree][jDegree] = new mydouble[nVar];

			flagh_hy[iDegree][jDegree] = new mydouble[nVar + 3];

			D_ddu[iDegree][jDegree] = new mydouble[nVar];
			D_ddu0[iDegree][jDegree] = new mydouble[nVar];
			D_un[iDegree][jDegree] = new mydouble[nVar];
			D_un0[iDegree][jDegree] = new mydouble[nVar];
			D_un1[iDegree][jDegree] = new mydouble[nVar];
			D_du_up1[iDegree][jDegree] = new mydouble[nVar];
			D_fxz[iDegree][jDegree] = new mydouble[nVar];
			D_fxf[iDegree][jDegree] = new mydouble[nVar];
			D_fyz[iDegree][jDegree] = new mydouble[nVar];
			D_fyf[iDegree][jDegree] = new mydouble[nVar];
			D_fxz_up1[iDegree][jDegree] = new mydouble[nVar];
			D_fyz_up1[iDegree][jDegree] = new mydouble[nVar];
			D_duvpc[iDegree][jDegree] = new mydouble[6];
			D_gx[iDegree][jDegree] = new mydouble[nVar];
			D_gy[iDegree][jDegree] = new mydouble[nVar];
			D_gx_t[iDegree][jDegree] = new mydouble[nVar];
			D_gy_t[iDegree][jDegree] = new mydouble[nVar];

			D_fx_char[iDegree][jDegree] = new mydouble[nVar];
			D_fy_char[iDegree][jDegree] = new mydouble[nVar];
			D_fx_up1_char[iDegree][jDegree] = new mydouble[nVar];
			D_fy_up1_char[iDegree][jDegree] = new mydouble[nVar];

			D_Coord[iDegree][jDegree] = new mydouble[2];
			D_Exact[iDegree][jDegree] = new mydouble[5];
			D_out[iDegree][jDegree] = new mydouble[5];

			D_fxat[iDegree][jDegree] = new mydouble *[nVar];
			D_fyat[iDegree][jDegree] = new mydouble *[nVar];

			for (long i = 0; i < nVar; i++)
			{
				D_fxat[iDegree][jDegree][i] = new mydouble[8];
				D_fyat[iDegree][jDegree][i] = new mydouble[8];
			}
		}
	}
}

void CADSolver::Delete_date(void)
{

	delete[] D_utx0;
	delete[] D_uty0;
	delete[] D_duvpc, D_un, D_un0, D_un1;
	delete[] D_du, D_du0, D_du1, D_ddu, D_ddu0, D_fxat, D_fyat;
}

CSolver::CSolver(void) {}

CSolver::~CSolver(void) {}

mydouble CADSolver::AddTest(mydouble a, mydouble b)
{
	return a + b;
}