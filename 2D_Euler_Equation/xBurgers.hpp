
/* Explicit Two-Derivative Two-Step Runge-Kutta Methods
 * Copyright (C) 2024
 * Author QinXueyu
 */

// #pragma once

#include <cmath>
#include <math.h>
#include <cstdlib>
#include <limits>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <iomanip>
#include <time.h>
#include <omp.h>
#include <chrono>

using namespace std;

typedef double mydouble;

class CSolver
{

public:
	CSolver(void);

	~CSolver(void);

	virtual void Driver(void);
};

class CADSolver : public CSolver
{

protected:
	mydouble nx_end, nx_begin, ny_end, ny_begin;				////////
	mydouble dx, dy;											////////
	unsigned long nCell_L, nCellYo_L, nCell_2L, nCellYo_2L, nb; // 扩展后
	unsigned long scheme_name1, show_n;							////////
	mydouble ska, D_max_u, D_max_v, T_max_u, T_max_v;			////////
	unsigned long SpaceScheme, SpaceScheme_cd, char_scheme;		////////
	unsigned long Case_name, BoundaryScheme, positivity_flag, positivity_flag_LW;
	unsigned long Local_Lax_Friedrichs, Case_1D_or_2D;
	mydouble C_v, C_p, Re, Pr, Ts, Tc, Ma;
	mydouble **D_mu, **D_ka;

	unsigned long Hy_flag, Viscous_flag, Viscous_flag_dt, LW_Viscous_cd, Viscous_cd;
	// ----------------------------------

	mydouble ***D_utx0, ***D_uty0, ***D_duvpc, ***D_du1, ***D_ddu1, ***D_du_up1, ***D_ddu_up1, ***D_dddu_up1;
	mydouble ***D_un0, ***D_du0, ***D_ddu0, ***D_un, ***D_du, ***D_ddu, ***D_dddu;
	mydouble ***D_du2, ***D_du3, ***D_du4;

	mydouble ***D_un1, ***D_dun1, ***D_ddun1, ***D_un2, ***D_dun2, ***D_ddun2, ***D_un3;
	mydouble ***D_uns1, ***D_duns1, ***D_dduns1;

	mydouble ***D_fxz, ***D_fxf, ***D_fyz, ***D_fyf, ***D_fxz_up1, ***D_fyz_up1;
	mydouble ***D_fx_char, ***D_fy_char, ***D_fx_up1_char, ***D_fy_up1_char;
	mydouble ***D_gx, ***D_gy, ***D_gx_t, ***D_gy_t, **D_ut, **D_vt, **D_pt, **D_Tt, **D_utt, **D_vtt, **D_Ttt;

	mydouble ****D_fxat, ****D_fyat, ***flagh_hy;
	mydouble ***D_Coord, ***D_Exact, ***D_out;
	mydouble S_d, S_u, S_v, S_p, S_ct;
	// ----------------------------------

	unsigned short NumberThreads;
	unsigned short nCellBegin, nCellEnd;
	unsigned short RungeKuttaStep, iRK_Step;
	unsigned long nCell, nPoint, nCellYo, nPointYo;
	unsigned long StepMax, IntStepMax;
	unsigned long ExtIter, IntIter, TotalIter;
	unsigned short iStep, iVar;
	unsigned long iCell, iPoint, jCell, jPoint;

	unsigned short iDegree, jDegree, kDegree, pDegree, qDegree;
	unsigned long iIterate, SaveStep, iSave;

	unsigned short ProblemName;

	unsigned short SolverName;
	string SolverNameChar;
	string SolverLinear;
	string SolverBurgers;
	string SolverEuler;
	string SolverNonLinear;

	unsigned short BCStyle;
	string BCStyleChar;
	string BCPeriodic;
	string BCExtraplated;
	string BCReflective;

	unsigned short TimeScheme;

	mydouble T_flag;

	mydouble Alpha;
	mydouble Aspeed;
	mydouble StopTime;
	mydouble Visc;
	mydouble CFL, amax_uc, amax_vc;

	mydouble Reynolds;

	unsigned short nDegree;
	unsigned short Degree_minus_One;
	unsigned short nDegreeSource;
	unsigned long nPointSP;

	mydouble TotalVol;

	mydouble DT, DT_w1, DTn1, DT_w2;
	mydouble DTPhysics;
	mydouble TimeNow;

	unsigned short nDim, iDim;
	unsigned short nVar;

	mydouble PI_Number;
	mydouble Gamma;
	mydouble Gamma_1;
	mydouble OffsetSol;

	mydouble CPUTime_One;
	mydouble CPUTime_Two, CPUTime_Three, CPUTime_Four, CPUTime_Five;

public:
	CADSolver(void);

	~CADSolver(void);

	void Driver(void);

	void Initial_date(void);

	void MesherAndAllocation_xy(void);
	void Mesher_out(void);

	mydouble scheme_df_dx(mydouble *left, mydouble *right, long scheme_name, long scheme_i);
	mydouble scheme_df_dx_char(mydouble *left, mydouble *right, long scheme_name, long scheme_i);
	mydouble scheme_df_dx_char_weno5(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_char_Teno5(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_char_Teno5A(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_char_Teno6(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_char_Teno7(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_char_Teno8(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_char_Teno88(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_char_weno5Z(mydouble *left, mydouble *right, long scheme_i);

	mydouble scheme_df_dx_char_weno7Z(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_char_weno7Z_Hy(mydouble *left, mydouble *right, long scheme_i, long iCell, long JCell, long nn);

	mydouble scheme_df_dx_char_weno70(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_weno70(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_weno7Z(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_weno5(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_TENO5(mydouble *left, mydouble *right, long scheme_i);
	mydouble scheme_df_dx_TENO7(mydouble *left, mydouble *right, long scheme_i);

	mydouble scheme_df_dx_char_weno7_z(mydouble *left, long scheme_i);
	mydouble scheme_df_dx_char_weno7_f(mydouble *left, long scheme_i);

	mydouble scheme_ALW_df_dx(mydouble *left, long scheme_name, long scheme_i);

	void Non_viscous_Splitting(void);
	void Splitting(void);
	void Splitting_char(void);
	void Splitting_char_n(void);

	void Splitting_char_x(mydouble nx, mydouble ny);
	void Splitting_char_y(mydouble nx, mydouble ny);
	void Splitting_char_x_n(mydouble nx, mydouble ny);
	void Splitting_char_y_n(mydouble nx, mydouble ny);

	void df_dx_dy(void);
	void df_dx_dy0(void);
	void df_dx_dy1(void);
	void comput_du_1up_positivity(void);
	void comput_du_1up_positivity0(void);

	void comput_ddu_1up_positivity(void);
	void comput_dddu_1up_positivity(void);

	void comput_err1(void);

	void uu_to_cc(void);

	void Boundary_du(void);
	void Boundary_du_periodic(void);
	void Boundary_du_Expand(void);
	void Boundary_du_reflection_x(void);

	void Boundary_ddu(void);
	void Boundary_ddu_periodic(void);
	void Boundary_ddu_Expand(void);
	void Boundary_ddu_reflection_x(void);

	void Boundary_u_t(void);
	void Boundary_u_t_periodic(void);
	void Boundary_u_tt(void);
	void Boundary_u_tt_periodic(void);

	void Boundary_uu(void);
	void Boundary_uu_Expand(void);
	void Boundary_uu_Expand_Shock_vortex(void);
	void Boundary_uu_periodic(void);
	void Boundary_uu_reflection_x(void);
	void Boundary_uu_Viscous_shock(void);
	void Boundary_uu_Couette_flow(void);
	void Boundary_uu_DoubleMach(void);
	void Boundary_uu_High_jet(void);
	void Boundary_uu_Laminar(void);
	void Boundary_uu_Laminar0(void);

	void Boundary_dut_periodic(void);
	void Boundary_dutxy_all(void);
	void Boundary_dutxy_periodic(void);
	void Boundary_dutxy_expand(void);
	void Boundary_dutxy_reflection_x(void);

	void Boundary_dutxy_LWA_periodic(void);

	void ComputeDT(void);
	void shock_sensor_Hy(void);
	void Viscous_Splitting(void);
	void Viscous_mu(void);
	void Viscous_Splitting_LW(void);
	void Viscous_Splitting_LW_cd(void);

	void Viscous_Schemes_LW(void);
	void Viscous_Schemes_LW_cd(void);
	void Viscous_Splitting_LW3(void);
	void Viscous_Schemes_LW3(void);
	void Viscous_Splitting_LW3_cd(void);
	void Viscous_Schemes_LW3_cd(void);

	void Viscous_Schemes(void);
	void Write_hy(long ii);
	void comput_du_ALW(void);
	void scheme_du_dx(void);
	void scheme_du_dx_New(void);

	void store_du_duu(void);

	void comput_du_LW_New(void);
	void comput_du_LW_Old(void);
	void comput_du_LW_New0(void);

	void comput_ddu_LW_New(void);
	void scheme_ddu_dx_New(void);

	void TimeIntegration(unsigned short iRK_Step);
	void TimeIntegration_RK11(unsigned short iRK_Step);

	void TimeIntegration_TD23_ALW(unsigned short iRK_Step);
	void TimeIntegration_TDTDRK_initial(unsigned short iRK_Step, double a_int);

	void TimeIntegration_LW3_ALW(void);
	void TimeIntegration_TDMSRK223_ALW(unsigned short iRK_Step);
	void TimeIntegration_TDMSRK224_ALW(unsigned short iRK_Step);
	void TimeIntegration_TDMSRK235_ALW(unsigned short iRK_Step);
	void TimeIntegration_TDMSRK326_ALW(unsigned short iRK_Step);

	void TimeIntegration_TDMSRK224_ALW_Hy(unsigned short iRK_Step);

	void TimeIntegration_TD23_ALW_Hy(unsigned short iRK_Step);
	void TimeIntegration_TD24_ALW(unsigned short iRK_Step);
	void TimeIntegration_TD34_ALW(unsigned short iRK_Step);
	void TimeIntegration_TD35_ALW(unsigned short iRK_Step);
	void TimeIntegration_TD46_ALW(unsigned short iRK_Step);

	void TimeIntegration_TDMS23_ALW(unsigned short iRK_Step);
	void TimeIntegration_TDMS34_ALW(unsigned short iRK_Step);
	void TimeIntegration_TDMS34_ALW_1(unsigned short iRK_Step);

	void RK33(void);
	void RK54(void);
	void RK65(void);
	void RK76(void);

	void RK22(void);
	void RK11(void);

	void TD23(void);
	void TD23_old(void);

	void TDMSRK223(void);
	void TDMSRK224(void);
	void TDMSRK235(void);
	void TDMSRK326(void);

	void TDMSRK224_Hy(void);

	void TD23_Hy(void);
	void TD24(void);
	void TD24_Old(void);
	void LW3_New(void);
	void LW3_New_0(void);
	void TD23_New(void);
	void TDTDRK_New_initial(double a_int);
	
	void TD24_New(void);
	void TD34_New(void);
	void TD35_New(void);
	void TD46_New(void);

	void TD34_Old(void);
	void TDMS23(void);
	void TDMS23_Old(void);
	void TDMS34(void);
	void TDMS34_Old(void);

	void TDMS23_New(void);
	void TDMS34_New(void);
	void TimeIntegration_RK22(unsigned short iRK_Step);
	void TimeIntegration_RK54(unsigned short iRK_Step);
	void TimeIntegration_RK65(unsigned short iRK_Step);
	void TimeIntegration_RK76(unsigned short iRK_Step);

	void InputControl_xy(void);
	void Intial_date(void);
	void Intial_zone(void);
	void Delete_date(void);

	void DeleteStore(void);
	void case_sod(mydouble a, mydouble b);
	void case_sod_2D(mydouble a, mydouble b);
	void case_Shock_vortex_2D(mydouble a, mydouble b);
	void case_Viscous_shock(mydouble a, mydouble b);
	void case_Blast_wave(mydouble a, mydouble b);
	void case_Kelvin_Helmholtz(mydouble a, mydouble b);
	void case_Couette_flow(mydouble a, mydouble b);
	void case_Laminar(mydouble a, mydouble b);
	void case_High_jet(mydouble a, mydouble b);
	void case_DoubleMach(mydouble a, mydouble b);
	void out_date(void);
	void out_date_expand(void);
	void out_date_1D(void);
	mydouble AddTest(mydouble a, mydouble b);
};
