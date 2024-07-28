// #pragma once
#include "xBurgers.hpp"
using namespace std;

// *(double(*)[5])(*(D_fxz+1)): [5]

void CADSolver::ComputeDT(void)
{
	// #pragma omp barrier
	mydouble dtlocal;
	mydouble chavel_x, chavel_y, uu, vv, cc, max_uc, max_vc, dt_vis, mindt_vis, dt_k, mindt_k;
	mydouble local_x, local_y;
	mydouble artv, dd2;
	max_uc = 0., max_vc = 0., mindt_vis = 1E10, mindt_k = 1E10;

	if (dx < dy)
		dd2 = dx * dx;
	else
		dd2 = dy * dy;

#pragma omp parallel for shared(max_uc, max_vc, ska, dt_vis, mindt_vis, dt_k, mindt_k, dd2)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb - 1; iCell < nCell_L; iCell++)
		{
#pragma omp critical
			{
				if (D_duvpc[iCell][jCell][4] > 1.E-18)
				{
					if (max_uc < fabs(D_duvpc[iCell][jCell][1]) + D_duvpc[iCell][jCell][4])
						max_uc = fabs(D_duvpc[iCell][jCell][1]) + D_duvpc[iCell][jCell][4];
					if (max_vc < fabs(D_duvpc[iCell][jCell][2]) + D_duvpc[iCell][jCell][4])
						max_vc = fabs(D_duvpc[iCell][jCell][2]) + D_duvpc[iCell][jCell][4];
				}
				else
					ska = 10;

				if (Viscous_flag_dt == 1)
				{
					// dt_vis = CFL * dd / ( D_mu[iCell][jCell] / D_duvpc[iCell][jCell][0]);
					dt_vis = CFL * dd2 / D_mu[iCell][jCell];
					if (dt_vis < mindt_vis)
						mindt_vis = dt_vis;
					dt_k = CFL * dd2 / D_ka[iCell][jCell];
					if (dt_k < mindt_k)
						mindt_k = dt_k;
				}

				if (D_duvpc[iCell][jCell][0] > -1.E5)
					1;
				else
					ska = 10;
			}
		}
	}
	// #pragma omp barrier
	// amax_uc = max_uc, amax_vc = max_vc;
	// max_uc = 3., max_vc = 3.;
	if (T_max_u > 0)
	{
		max_uc = T_max_u; // 固定时间步长
		max_vc = T_max_v;
	}

	if (Case_1D_or_2D == 1)
		max_vc = 0.;
	DTn1 = DT;

	DT = CFL / (max_uc / dx + max_vc / dy);

	// cout << "Step = " << ExtIter << " mindt_k = " << mindt_k << endl;
	// cout << "Step = " << ExtIter << " mindt_vis = " << mindt_vis << endl;

	// cout << "Step = " << ExtIter << " max_uc = " << max_uc << endl;
	// cout << "Step = " << ExtIter << " max_vc = " << max_vc << endl;

	if ((TimeNow + DT) > StopTime - 1.E-10)
	{
		DT = StopTime - TimeNow;
		T_flag = 10.;
		cout << "Step = " << ExtIter << "  c DT = " << DT  << endl;
	}
	// DT_w2 = DT_w1;
	// DT_w1 = DTn1 / DT;
	// if ((DT_w1 < 0.8 || DT_w1 > 1.3) && ExtIter > 1)
	// {
	// 	cout << "Step = " << ExtIter << " DT_w1 = " << DT_w1 << endl;
	// }
	if (ExtIter < 5)
		cout << "Step = " << ExtIter << "  c DT = " << DT <<  "dx=" << dx << "    dy=" << dy << endl;
}

void CADSolver::uu_to_cc(void) // uu_to_cc
{
	// #pragma omp barrier
	mydouble pressure, density, velocity_x, velocity_y, cc;
	long jCell, iCell, nn;
// #pragma omp parallel for num_threads(NumberThreads)
#pragma omp parallel for private(pressure, density, velocity_x, velocity_y, cc, jCell, iCell, nn)
	for (jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu[5], uv2;
		for (iCell = 0; iCell < nCell_2L; iCell++)
		{
			for (nn = 0; nn < 4; nn++)
				uuu[nn] = D_un[iCell][jCell][nn];

			density = uuu[0], velocity_x = uuu[1] / density;
			velocity_y = uuu[2] / density;
			uv2 = velocity_x * velocity_x + velocity_y * velocity_y;
			pressure = Gamma_1 * (uuu[3] - 0.5 * density * uv2);
			cc = sqrt(Gamma * pressure / density);
			if (cc > 1.E-18)
				1;
			else
			{
				// cout << "cc." << cc;
				cc = 1E-10;
			}
			D_duvpc[iCell][jCell][0] = density;
			D_duvpc[iCell][jCell][1] = velocity_x;
			D_duvpc[iCell][jCell][2] = velocity_y;
			D_duvpc[iCell][jCell][3] = pressure;
			D_duvpc[iCell][jCell][4] = cc;
			D_duvpc[iCell][jCell][5] = (uuu[3] / density - 0.5 * uv2) / C_v;
		}
	}
	// #pragma omp barrier
}

//___________________du_char__________________
void CADSolver::Non_viscous_Splitting(void)
{

	switch (char_scheme)
	{
	case 0:
		Splitting();
		df_dx_dy();
		break;
	case 1:
		Splitting();
		Splitting_char_n();
		break;

	default:
		cout << "Char_scheme is not avaliable for this Solver." << endl;
		break;
	}
}

void CADSolver::Splitting_char_n(void)
{
	// #pragma omp barrier

	mydouble ss1 = 0., ss2 = 0., pp1 = 0., pp2 = 0.;

	Splitting_char_x_n(1.0, 0.);
	if (Case_1D_or_2D == 2)
		Splitting_char_y_n(0., 1.0);

		// #pragma omp barrier

#pragma omp parallel for private(ss1, ss2, pp1, pp2)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				// mydouble ***D_fx_char, ***D_fy_char, ***D_fx_up1_char, ***D_fy_up1_char;
				ss1 = (D_fx_char[iCell][jCell][nn] - D_fx_char[iCell - 1][jCell][nn]) / dx;
				if (Case_1D_or_2D == 2)
					ss2 = (D_fy_char[iCell][jCell][nn] - D_fy_char[iCell][jCell - 1][nn]) / dy;
				else
					ss2 = 0.;

				D_du[iCell][jCell][nn] = ss1 + ss2;

				if (positivity_flag == 1)
				{
					pp1 = (D_fx_up1_char[iCell][jCell][nn] - D_fx_up1_char[iCell - 1][jCell][nn]) / dx;
					if (Case_1D_or_2D == 2)
						pp2 = (D_fy_up1_char[iCell][jCell][nn] - D_fy_up1_char[iCell][jCell - 1][nn]) / dy;
					else
						pp2 = 0.;
					D_du_up1[iCell][jCell][nn] = pp1 + pp2;
				}
			}
		}
	}
	// #pragma omp barrier
	if (positivity_flag == 1)
		comput_du_1up_positivity();

	if (Hy_flag == 1)
		shock_sensor_Hy();
}

void CADSolver::Splitting_char_x_n(mydouble nx, mydouble ny)
{
	// Book: Riemann Solversand Numerical Methods for Fluid Dynamics, Third Edition ||  Toro PP.108 275 355
	// Book: I Do Like CFD ||  pp.69
	// Report : Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes
	// for Hyperbolic Conservation Laws. || Chi-Wang Shu pp.31  NASA/CR-97-206253 ICASE Report No. 97-65

#pragma omp parallel for
	for (long jCell = nb - 4; jCell < nCellYo_L + 5; jCell++)
	{
		mydouble fxz[5], fxz_up1[5], fyz[5], RM[4][4], LM[4][4];
		mydouble a1, a2, u1, u2, v1, v2, c1, c2, H1, H2, r1, r2, r0, u0, v0, H0, c0, uv, uvxy, lmax, Zc, c0_2, c0c0_1, c0u0;
		// mydouble lamda[4][4] = {0}, lamda_z[4][4] = {0}, lamda_f[4][4] = {0};
		mydouble qz[15] = {0}, qf[15] = {0}, wz[4][4] = {0}, wf[4][4] = {0};
		long n_w = 5, mm;
		for (long iCell = nb - 4; iCell < nCell_L + 5; iCell++)
		{
			a1 = D_duvpc[iCell][jCell][0];
			u1 = D_duvpc[iCell][jCell][1];
			v1 = D_duvpc[iCell][jCell][2];
			c1 = D_duvpc[iCell][jCell][4];
			H1 = c1 * c1 / Gamma_1 + 0.5 * (u1 * u1 + v1 * v1); // ! 总焓

			a2 = D_duvpc[iCell + 1][jCell][0];
			u2 = D_duvpc[iCell + 1][jCell][1];
			v2 = D_duvpc[iCell + 1][jCell][2];
			c2 = D_duvpc[iCell + 1][jCell][4];
			H2 = c2 * c2 / Gamma_1 + 0.5 * (u2 * u2 + v2 * v2); // ! 总焓
			// !        Roe 平均
			r1 = sqrt(a1), r2 = sqrt(a2), r0 = r1 + r2;
			u0 = (r1 * u1 + r2 * u2) / r0;
			v0 = (r1 * v1 + r2 * v2) / r0;
			H0 = (r1 * H1 + r2 * H2) / r0;
			uv = u0 * u0 + v0 * v0;
			c0 = sqrt(Gamma_1 * (H0 - uv * 0.5)); // ! j+1/2点处的声速

			// uvxy = nx * u0 + ny * v0;

			Zc = Gamma_1 / (c0 * c0), c0_2 = 1. / (2. * c0);
			c0u0 = c0 * u0;
			RM[0][0] = 1.0, RM[0][1] = 1.0, RM[0][2] = 1.0, RM[0][3] = 0.0;
			RM[1][0] = -c0 + u0, RM[1][1] = u0, RM[1][2] = c0 + u0, RM[1][3] = 0.0;
			RM[2][0] = v0, RM[2][1] = v0, RM[2][2] = v0, RM[2][3] = 1;
			RM[3][0] = H0 - c0u0, RM[3][1] = uv * 0.5, RM[3][2] = H0 + c0u0, RM[3][3] = v0;

			LM[0][0] = 0.25 * uv * Zc + c0_2 * u0;
			LM[0][1] = -(1. + c0u0 * Zc) * c0_2;
			LM[0][2] = -(c0 * v0 * Zc) * c0_2;
			LM[0][3] = Zc * 0.5;

			LM[1][0] = 1. - uv * Zc * 0.5, LM[1][1] = u0 * Zc, LM[1][2] = v0 * Zc, LM[1][3] = -Zc;
			LM[2][0] = 0.25 * uv * Zc - c0_2 * u0;
			LM[2][1] = (1. - c0u0 * Zc) * c0_2;
			LM[2][2] = -c0 * v0 * Zc * c0_2;
			LM[2][3] = Zc * 0.5;
			LM[3][0] = -v0, LM[3][1] = 0.0, LM[3][2] = 1, LM[3][3] = 0.;

			// RM[0][0] = 1.0, RM[0][1] = 1.0, RM[0][2] = 1.0, RM[0][3] = 0.0;
			// RM[1][0] = -c0 * nx + u0, RM[1][1] = u0, RM[1][2] = c0 * nx + u0, RM[1][3] = -ny;
			// RM[2][0] = -c0 * ny + v0, RM[2][1] = v0, RM[2][2] = c0 * ny + v0, RM[2][3] = nx;
			// RM[3][0] = H0 - c0 * uvxy, RM[3][1] = uv / 2., RM[3][2] = H0 + c0 * uvxy, RM[3][3] = -ny * u0 + nx * v0;
			// LM[0][0] = 1. / 4. * ((2. * uvxy) / c0 + uv * Zc);
			// LM[0][1] = -((nx + c0 * u0 * Zc) / (2. * c0));
			// LM[0][2] = -((ny + c0 * v0 * Zc) / (2. * c0));
			// LM[0][3] = Zc / 2.;
			// LM[1][0] = 1. - (uv * Zc) / 2., LM[1][1] = u0 * Zc, LM[1][2] = v0 * Zc, LM[1][3] = -Zc;
			// LM[2][0] = 1. / 4. * (-((2. * uvxy) / c0) + uv * Zc);
			// LM[2][1] = (nx - c0 * u0 * Zc) / (2. * c0);
			// LM[2][2] = (ny - c0 * v0 * Zc) / (2. * c0);
			// LM[2][3] = Zc / 2.;
			// LM[3][0] = ny * u0 - nx * v0, LM[3][1] = -ny, LM[3][2] = nx, LM[3][3] = 0.;

			for (long i = 0; i < 4; i++)
			{
				for (long m = -n_w; m < n_w + 1; m++)
				{
					mm = m + n_w;
					qz[mm] = 0, qf[mm] = 0;
					for (long k = 0; k < 4; k++)
					{
						qz[mm] += LM[i][k] * D_fxz[iCell + m][jCell][k];
						qf[mm] += LM[i][k] * D_fxf[iCell + m][jCell][k];
					}
				}
				// if (Hy_flag == 1)
				// 	fxz[i] = scheme_df_dx_char_weno7Z_Hy(qz, qf, n_w, iCell, jCell, i);
				// else
				fxz[i] = scheme_df_dx_char(qz, qf, SpaceScheme, n_w);
				if (positivity_flag == 1)
					fxz_up1[i] = qz[n_w] + qf[n_w + 1];
			}

			for (long i = 0; i < 4; i++)
			{
				D_fx_char[iCell][jCell][i] = 0.;
				for (long k = 0; k < 4; k++)
					D_fx_char[iCell][jCell][i] += RM[i][k] * fxz[k];

				if (positivity_flag == 1)
				{
					D_fx_up1_char[iCell][jCell][i] = 0.;
					for (long k = 0; k < 4; k++)
						D_fx_up1_char[iCell][jCell][i] += RM[i][k] * fxz_up1[k];
				}
			}
		}
	}
}

void CADSolver::Splitting_char_y_n(mydouble nx, mydouble ny)
{
	// Book: Riemann Solversand Numerical Methods for Fluid Dynamics, Third Edition  Toro PP.108 275 355
	// Book: I Do Like CFD   pp.69

#pragma omp parallel for
	for (long jCell = nb - 4; jCell < nCellYo_L + 5; jCell++)
	{
		mydouble fxz[5], fyz[5], fyz_up1[5], RM[4][4], LM[4][4];
		mydouble a1, a2, u1, u2, v1, v2, c1, c2, H1, H2, r1, r2, r0, u0, v0, H0, c0, uv, uvxy, lmax, Zc, c0_2, c0c0_1, c0v0;
		// mydouble lamda[4][4] = {0}, lamda_z[4][4] = {0}, lamda_f[4][4] = {0};
		mydouble qz[15] = {0}, qf[15] = {0}, wz[4][4] = {0}, wf[4][4] = {0};
		long n_w = 5, mm;
		for (long iCell = nb - 4; iCell < nCell_L + 5; iCell++)
		{

			a1 = D_duvpc[iCell][jCell][0];
			u1 = D_duvpc[iCell][jCell][1];
			v1 = D_duvpc[iCell][jCell][2];
			c1 = D_duvpc[iCell][jCell][4];
			H1 = c1 * c1 / Gamma_1 + 0.5 * (u1 * u1 + v1 * v1); // ! 总焓

			a2 = D_duvpc[iCell][jCell + 1][0];
			u2 = D_duvpc[iCell][jCell + 1][1];
			v2 = D_duvpc[iCell][jCell + 1][2];
			c2 = D_duvpc[iCell][jCell + 1][4];
			H2 = c2 * c2 / Gamma_1 + 0.5 * (u2 * u2 + v2 * v2); // ! 总焓
			// !        Roe 平均
			r1 = sqrt(a1), r2 = sqrt(a2), r0 = r1 + r2;
			u0 = (r1 * u1 + r2 * u2) / r0;
			v0 = (r1 * v1 + r2 * v2) / r0;
			H0 = (r1 * H1 + r2 * H2) / r0;
			uv = u0 * u0 + v0 * v0;
			c0 = sqrt(Gamma_1 * (H0 - uv * 0.5)); // ! j+1/2点处的声速

			// uvxy = nx * u0 + ny * v0;
			Zc = Gamma_1 / (c0 * c0), c0_2 = 1. / (2. * c0);
			c0v0 = c0 * v0;
			RM[0][0] = 1.0, RM[0][1] = 1.0, RM[0][2] = 1.0, RM[0][3] = 0.0;
			RM[1][0] = u0, RM[1][1] = u0, RM[1][2] = u0, RM[1][3] = -1.0;
			RM[2][0] = -c0 + v0, RM[2][1] = v0, RM[2][2] = c0 + v0, RM[2][3] = 0.0;
			RM[3][0] = H0 - c0v0, RM[3][1] = uv * 0.5, RM[3][2] = H0 + c0v0, RM[3][3] = -u0;

			LM[0][0] = 0.25 * uv * Zc + c0_2 * v0;
			LM[0][1] = -(c0 * u0 * Zc) * c0_2;
			LM[0][2] = -(1.0 + c0v0 * Zc) * c0_2;
			LM[0][3] = Zc * 0.5;

			LM[1][0] = 1. - uv * Zc * 0.5, LM[1][1] = u0 * Zc, LM[1][2] = v0 * Zc, LM[1][3] = -Zc;
			LM[2][0] = 0.25 * uv * Zc - c0_2 * v0;
			LM[2][1] = -c0 * u0 * Zc * c0_2;
			LM[2][2] = (1.0 - c0v0 * Zc) * c0_2;
			LM[2][3] = Zc * 0.5;
			LM[3][0] = u0, LM[3][1] = -1.0, LM[3][2] = 0.0, LM[3][3] = 0.;

			// RM[0][0] = 1.0, RM[0][1] = 1.0, RM[0][2] = 1.0, RM[0][3] = 0.0;
			// RM[1][0] = -c0 * nx + u0, RM[1][1] = u0, RM[1][2] = c0 * nx + u0, RM[1][3] = -ny;
			// RM[2][0] = -c0 * ny + v0, RM[2][1] = v0, RM[2][2] = c0 * ny + v0, RM[2][3] = nx;
			// RM[3][0] = H0 - c0 * uvxy, RM[3][1] = uv / 2., RM[3][2] = H0 + c0 * uvxy, RM[3][3] = -ny * u0 + nx * v0;

			// LM[0][0] = 1. / 4. * ((2. * uvxy) / c0 + uv * Zc);
			// LM[0][1] = -((nx + c0 * u0 * Zc) / (2. * c0));
			// LM[0][2] = -((ny + c0 * v0 * Zc) / (2. * c0));
			// LM[0][3] = Zc / 2.;

			// LM[1][0] = 1. - (uv * Zc) / 2., LM[1][1] = u0 * Zc, LM[1][2] = v0 * Zc, LM[1][3] = -Zc;

			// LM[2][0] = 1. / 4. * (-((2. * uvxy) / c0) + uv * Zc);
			// LM[2][1] = (nx - c0 * u0 * Zc) / (2. * c0);
			// LM[2][2] = (ny - c0 * v0 * Zc) / (2. * c0);
			// LM[2][3] = Zc / 2.;
			// LM[3][0] = ny * u0 - nx * v0, LM[3][1] = -ny, LM[3][2] = nx, LM[3][3] = 0.;

			for (long i = 0; i < 4; i++)
			{
				for (long m = -n_w; m < n_w + 1; m++)
				{
					mm = m + n_w;
					qz[mm] = 0, qf[mm] = 0;
					for (long k = 0; k < 4; k++)
					{
						qz[mm] += LM[i][k] * D_fyz[iCell][jCell + m][k];
						qf[mm] += LM[i][k] * D_fyf[iCell][jCell + m][k];
					}
				}

				fyz[i] = scheme_df_dx_char(qz, qf, SpaceScheme, n_w);
				if (positivity_flag == 1)
					fyz_up1[i] = qz[n_w] + qf[n_w + 1];
			}

			for (long i = 0; i < 4; i++)
			{
				D_fy_char[iCell][jCell][i] = 0.;
				for (long k = 0; k < 4; k++)
					D_fy_char[iCell][jCell][i] += RM[i][k] * fyz[k];

				if (positivity_flag == 1)
				{
					D_fy_up1_char[iCell][jCell][i] = 0.;
					for (long k = 0; k < 4; k++)
						D_fy_up1_char[iCell][jCell][i] += RM[i][k] * fyz_up1[k];
				}
			}
		}
	}
}

mydouble CADSolver::scheme_df_dx_char(mydouble *fz, mydouble *ff, long scheme_name, long i)
{
	mydouble sz1, sf1, ss1 = 0., zz = 0.;
	// long ii = i;

	switch (scheme_name)
	{
	case 1:
		sz1 = (fz[i]);	   // 正通量
		sf1 = (ff[i + 1]); // 负通量
		zz = sz1 + sf1;
		break;

	case 50:
		// WENO5-JS
		zz = scheme_df_dx_char_weno5(fz, ff, i);
		break;

	case 52:
		// WENO5-Z
		zz = scheme_df_dx_char_weno5Z(fz, ff, i);
		break;
	case 53:
		// TENO5
		zz = scheme_df_dx_char_Teno5(fz, ff, i);
		break;

	case 54:
		// TENO5A
		zz = scheme_df_dx_char_Teno5A(fz, ff, i);
		break;
	case 72:
		// WENO7-Z
		zz = scheme_df_dx_char_weno7Z(fz, ff, i);
		break;
	case 70:
		// WENO7-JS
		// sz1 = scheme_df_dx_char_weno7_z(fz, i);
		// sf1 = scheme_df_dx_char_weno7_f(ff, i);
		// zz = sz1 + sf1;

		zz = scheme_df_dx_char_weno70(fz, ff, i);
		break;

	case 6:
		// TENO5
		zz = scheme_df_dx_char_Teno6(fz, ff, i);
		break;

	case 73:
		// TENO5
		zz = scheme_df_dx_char_Teno7(fz, ff, i);
		break;
	case 8:
		// TENO5
		zz = scheme_df_dx_char_Teno8(fz, ff, i);
		break;

	default:
		cout << "1up Space_Scheme  char is not avaliable for this Solver." << endl;
		break;
	}

	return zz;
}

mydouble CADSolver::scheme_df_dx_char_weno5Z(mydouble *ffz, mydouble *fff, long i)
{
	mydouble yy = 0.0, ss1, ss2;
	// /*  //weno_5-js  正通量//
	mydouble ep = 1.E-15, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	mydouble sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, tau;
	sn0 = 13. / 12. * pow((*(ffz + i - 2) - 2. * (*(ffz + i - 1)) + *(ffz + i)), 2) + 0.25 * pow((*(ffz + i - 2) - 4. * (*(ffz + i - 1)) + 3. * (*(ffz + i))), 2);
	sn1 = 13. / 12. * pow((*(ffz + i - 1) - 2. * (*(ffz + i)) + *(ffz + i + 1)), 2) + 0.25 * pow((*(ffz + i - 1) - *(ffz + i + 1)), 2);
	sn2 = 13. / 12. * pow((*(ffz + i) - 2. * (*(ffz + i + 1)) + *(ffz + i + 2)), 2) + 0.25 * pow((3. * (*(ffz + i)) - 4. * (*(ffz + i + 1)) + *(ffz + i + 2)), 2);

	tau = abs(sn2 - sn0);
	az0 = C03 * (1.0 + pow(tau / (ep + sn0), 2));
	az1 = C13 * (1.0 + pow(tau / (ep + sn1), 2));
	az2 = C23 * (1.0 + pow(tau / (ep + sn2), 2));
	// az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
	W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
	q03 = 2. / 6. * (*(ffz + i - 2)) - 7. / 6. * (*(ffz + i - 1)) + 11. / 6. * (*(ffz + i));
	q13 = -1. / 6. * (*(ffz + i - 1)) + 5. / 6. * (*(ffz + i)) + 2. / 6. * (*(ffz + i + 1));
	q23 = 2. / 6. * (*(ffz + i)) + 5. / 6. * (*(ffz + i + 1)) - 1. / 6. * (*(ffz + i + 2));
	ss1 = W0 * q03 + W1 * q13 + W2 * q23;
	// weno_5-js 负通量
	//  double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	//  double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
	// C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	sn0 = 13. / 12. * pow((*(fff + i + 1) - 2. * (*(fff + i + 2)) + *(fff + i + 3)), 2) + 0.25 * pow((3. * (*(fff + i + 1)) - 4. * (*(fff + i + 2)) + *(fff + i + 3)), 2);
	sn1 = 13. / 12. * pow((*(fff + i) - 2. * (*(fff + i + 1)) + *(fff + i + 2)), 2) + 0.25 * pow((*(fff + i) - *(fff + i + 2)), 2);
	sn2 = 13. / 12. * pow((*(fff + i - 1) - 2. * (*(fff + i)) + *(fff + i + 1)), 2) + 0.25 * pow((*(fff + i - 1) - 4. * (*(fff + i)) + 3. * (*(fff + i + 1))), 2);

	tau = abs(sn2 - sn0);
	az0 = C03 * (1.0 + pow(tau / (ep + sn0), 2));
	az1 = C13 * (1.0 + pow(tau / (ep + sn1), 2));
	az2 = C23 * (1.0 + pow(tau / (ep + sn2), 2));
	// az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
	W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
	q03 = 11. / 6. * (*(fff + i + 1)) - 7. / 6. * (*(fff + i + 2)) + 2. / 6. * (*(fff + i + 3));
	q13 = 2. / 6. * (*(fff + i)) + 5. / 6. * (*(fff + i + 1)) - 1. / 6. * (*(fff + i + 2));
	q23 = -1. / 6. * (*(fff + i - 1)) + 5. / 6. * (*(fff + i)) + 2. / 6. * (*(fff + i + 1));
	ss2 = W0 * q03 + W1 * q13 + W2 * q23;
	yy = ss1 + ss2;
	// */

	return yy;
}

mydouble CADSolver::scheme_df_dx_char_weno5(mydouble *ffz, mydouble *fff, long i)
{
	mydouble yy = 0.0, ss1, ss2;
	// /*  //weno_5-js  正通量//
	mydouble ep = 1.E-15, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	mydouble sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23;
	sn0 = 13. / 12. * pow((*(ffz + i - 2) - 2. * (*(ffz + i - 1)) + *(ffz + i)), 2) + 0.25 * pow((*(ffz + i - 2) - 4. * (*(ffz + i - 1)) + 3. * (*(ffz + i))), 2);
	sn1 = 13. / 12. * pow((*(ffz + i - 1) - 2. * (*(ffz + i)) + *(ffz + i + 1)), 2) + 0.25 * pow((*(ffz + i - 1) - *(ffz + i + 1)), 2);
	sn2 = 13. / 12. * pow((*(ffz + i) - 2. * (*(ffz + i + 1)) + *(ffz + i + 2)), 2) + 0.25 * pow((3. * (*(ffz + i)) - 4. * (*(ffz + i + 1)) + *(ffz + i + 2)), 2);
	az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
	W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
	q03 = 2. / 6. * (*(ffz + i - 2)) - 7. / 6. * (*(ffz + i - 1)) + 11. / 6. * (*(ffz + i));
	q13 = -1. / 6. * (*(ffz + i - 1)) + 5. / 6. * (*(ffz + i)) + 2. / 6. * (*(ffz + i + 1));
	q23 = 2. / 6. * (*(ffz + i)) + 5. / 6. * (*(ffz + i + 1)) - 1. / 6. * (*(ffz + i + 2));
	ss1 = W0 * q03 + W1 * q13 + W2 * q23;
	// weno_5-js 负通量
	//  double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	//  double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
	// C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	sn0 = 13. / 12. * pow((*(fff + i + 1) - 2. * (*(fff + i + 2)) + *(fff + i + 3)), 2) + 0.25 * pow((3. * (*(fff + i + 1)) - 4. * (*(fff + i + 2)) + *(fff + i + 3)), 2);
	sn1 = 13. / 12. * pow((*(fff + i) - 2. * (*(fff + i + 1)) + *(fff + i + 2)), 2) + 0.25 * pow((*(fff + i) - *(fff + i + 2)), 2);
	sn2 = 13. / 12. * pow((*(fff + i - 1) - 2. * (*(fff + i)) + *(fff + i + 1)), 2) + 0.25 * pow((*(fff + i - 1) - 4. * (*(fff + i)) + 3. * (*(fff + i + 1))), 2);
	az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
	W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
	q03 = 11. / 6. * (*(fff + i + 1)) - 7. / 6. * (*(fff + i + 2)) + 2. / 6. * (*(fff + i + 3));
	q13 = 2. / 6. * (*(fff + i)) + 5. / 6. * (*(fff + i + 1)) - 1. / 6. * (*(fff + i + 2));
	q23 = -1. / 6. * (*(fff + i - 1)) + 5. / 6. * (*(fff + i)) + 2. / 6. * (*(fff + i + 1));
	ss2 = W0 * q03 + W1 * q13 + W2 * q23;
	yy = ss1 + ss2;
	// */

	return yy;
}

mydouble CADSolver::scheme_df_dx_char_Teno5A(mydouble *ffz, mydouble *fff, long i)
{
	mydouble yy = 0.0, ss1, ss2;

	mydouble varepsilon = 1.0e-12;
	mydouble a1 = 10.5, a2 = 3.5, cr = 0.25, csi = 1.0e-3;
	mydouble ct = 1.0e-7;

	mydouble q0, q1, q2, q3, IS0, IS1, IS2, Is3, eps, eta0, eta1, eta2, rm;
	mydouble alpha0, alpha1, alpha2, alpha3;
	mydouble omega0, omega1, omega2, omega3;
	mydouble isumAlpha, absTerm0, absTerm1, d0, d1, d2, fhatp, fhatm, dsum;

	mydouble rp = 1e-12;
	mydouble is1par = 13.0 / 12.0;
	mydouble is2par = 0.25;

	mydouble a[3][3] = {{2.0 / 6.0, -7.0 / 6.0, 11.0 / 6.0},
						{-1.0 / 6.0, 5.0 / 6.0, 2.0 / 6.0},
						{2.0 / 6.0, 5.0 / 6.0, -1.0 / 6.0}};

	mydouble cw[4] = {0.1, 0.6, 0.3};
	// polinomia
	q0 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
	q1 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
	q2 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];

	// smoothness index
	IS0 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
	IS1 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
	IS2 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);

	// alpha
	absTerm0 = abs(IS0 - IS2);
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;

	// adaptive ct
	eps = 0.9 * cr / (1.0 - 0.9 * cr) * pow(csi, 2);
	eta0 = (2.0 * abs((ffz[i + 0] - ffz[i - 1]) * (ffz[i - 1] - ffz[i - 2])) + eps) /
		   ((ffz[i + 0] - ffz[i - 1]) * (ffz[i + 0] - ffz[i - 1]) + (ffz[i - 1] - ffz[i - 2]) * (ffz[i - 1] - ffz[i - 2]) + eps);
	eta1 = (2.0 * abs((ffz[i + 1] - ffz[i + 0]) * (ffz[i + 0] - ffz[i - 1])) + eps) /
		   ((ffz[i + 1] - ffz[i + 0]) * (ffz[i + 1] - ffz[i + 0]) + (ffz[i + 0] - ffz[i - 1]) * (ffz[i + 0] - ffz[i - 1]) + eps);
	eta2 = (2.0 * abs((ffz[i + 2] - ffz[i + 1]) * (ffz[i + 1] - ffz[i + 0])) + eps) /
		   ((ffz[i + 2] - ffz[i + 1]) * (ffz[i + 2] - ffz[i + 1]) + (ffz[i + 1] - ffz[i + 0]) * (ffz[i + 1] - ffz[i + 0]) + eps);

	rm = 1.0 - min(1.0, min(eta0, min(eta1, eta2)) / cr);

	ct = 1.0 * pow(10.0, -floor(a1 - a2 * (1.0 - pow((1.0 - rm), 4) * (1.0 + 4.0 * rm))));

	// delta
	d0 = omega0 < ct ? 0. : 1.;
	d1 = omega1 < ct ? 0. : 1.;
	d2 = omega2 < ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;

	// WENO plus reconstruction
	fhatp = omega0 * q0 + omega1 * q1 + omega2 * q2;

	// --- WENO minus part --- !

	// polinomia
	q0 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
	q1 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
	q2 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];

	// smoothness index
	IS0 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
	IS1 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
	IS2 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);

	// alpha
	absTerm0 = abs(IS0 - IS2);
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;

	// adaptive ct
	eta0 = (2.0 * abs((fff[i + 1] - fff[i + 2]) * (fff[i + 2] - fff[i + 3])) + eps) /
		   ((fff[i + 1] - fff[i + 2]) * (fff[i + 1] - fff[i + 2]) + (fff[i + 2] - fff[i + 3]) * (fff[i + 2] - fff[i + 3]) + eps);
	eta1 = (2.0 * abs((fff[i + 0] - fff[i + 1]) * (fff[i + 1] - fff[i + 2])) + eps) /
		   ((fff[i + 0] - fff[i + 1]) * (fff[i + 0] - fff[i + 1]) + (fff[i + 1] - fff[i + 2]) * (fff[i + 1] - fff[i + 2]) + eps);
	eta2 = (2.0 * abs((fff[i - 1] - fff[i + 0]) * (fff[i + 0] - fff[i + 1])) + eps) /
		   ((fff[i - 1] - fff[i + 0]) * (fff[i - 1] - fff[i + 0]) + (fff[i + 0] - fff[i + 1]) * (fff[i + 0] - fff[i + 1]) + eps);

	rm = 1.0 - min(1.0, min(eta0, min(eta1, eta2)) / cr);

	ct = 1.0 * pow(10.0, -floor(a1 - a2 * (1.0 - pow((1.0 - rm), 4) * (1.0 + 4.0 * rm))));

	// delta

	d0 = omega0 < ct ? 0. : 1.;
	d1 = omega1 < ct ? 0. : 1.;
	d2 = omega2 < ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;

	// WENO minus reconstruction
	fhatm = omega0 * q0 + omega1 * q1 + omega2 * q2;

	yy = (fhatp + fhatm);
	return yy;
}

mydouble CADSolver::scheme_df_dx_char_Teno5(mydouble *ffz, mydouble *fff, long i)
{
	mydouble yy = 0.0;

	mydouble varepsilon = 1.0e-16;

	// mydouble S_ct = 1.0e-7;
	// mydouble S_ct = 1.0e-3; // High_jet case

	mydouble q0, q1, q2, q3, IS0, IS1, IS2, Is3, eps, eta0, eta1, eta2, rm;
	mydouble alpha0, alpha1, alpha2, alpha3;
	mydouble omega0, omega1, omega2, omega3;
	mydouble isumAlpha, absTerm0, absTerm1, d0, d1, d2, fhatp, fhatm, dsum;

	mydouble is1par = 13.0 / 12.0;
	mydouble is2par = 0.25;

	mydouble a[3][3] = {{2.0 / 6.0, -7.0 / 6.0, 11.0 / 6.0},
						{-1.0 / 6.0, 5.0 / 6.0, 2.0 / 6.0},
						{2.0 / 6.0, 5.0 / 6.0, -1.0 / 6.0}};

	mydouble cw[4] = {0.1, 0.6, 0.3};
	// polinomia
	q0 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
	q1 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
	q2 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];

	// smoothness index
	IS0 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
	IS1 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
	IS2 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);

	// alpha
	absTerm0 = abs(IS0 - IS2);
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;

	// delta
	d0 = omega0 < S_ct ? 0. : 1.;
	d1 = omega1 < S_ct ? 0. : 1.;
	d2 = omega2 < S_ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;

	// WENO plus reconstruction
	fhatp = omega0 * q0 + omega1 * q1 + omega2 * q2;

	// --- WENO minus part --- !

	// polinomia
	q0 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
	q1 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
	q2 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];

	// smoothness index
	IS0 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
	IS1 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
	IS2 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);

	// alpha
	absTerm0 = abs(IS0 - IS2);
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;

	// delta

	d0 = omega0 < S_ct ? 0. : 1.;
	d1 = omega1 < S_ct ? 0. : 1.;
	d2 = omega2 < S_ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;

	// WENO minus reconstruction
	fhatm = omega0 * q0 + omega1 * q1 + omega2 * q2;

	yy = (fhatp + fhatm);
	return yy;
}

mydouble CADSolver::scheme_df_dx_char_Teno6(mydouble *ffz, mydouble *fff, long i)
{
	mydouble yy = 0.0;

	mydouble varepsilon = 1.0e-16;

	// mydouble S_ct = 1.0e-7;
	// mydouble S_ct = 1.0e-3; // High_jet case

	mydouble q0, q1, q2, q3, IS0, IS1, IS2, IS3, eps, eta0, eta1, eta2, rm;
	mydouble alpha0, alpha1, alpha2, alpha3;
	mydouble omega0, omega1, omega2, omega3;
	mydouble isumAlpha, absTerm0, absTerm1, d0, d1, d2, d3, fhatp, fhatm, dsum;

	mydouble is1par = 13.0 / 12.0;
	mydouble is2par = 0.25;

	mydouble a[4][4] = {{2.0 / 6.0, -7.0 / 6.0, 11.0 / 6.0, 0},
						{-1.0 / 6.0, 5.0 / 6.0, 2.0 / 6.0, 0},
						{2.0 / 6.0, 5.0 / 6.0, -1.0 / 6.0, 0},
						{3.0 / 12.0, 13.0 / 12.0, -5.0 / 12.0, 1.0 / 12.0}};

	// mydouble cw[4] = {0.1, 0.6, 0.3};
	mydouble cw[4] = {1. / 20., 9. / 20., 6. / 20., 4. / 20.};
	// polinomia

	q0 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
	q1 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
	q2 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];
	q3 = a[3][0] * ffz[i + 0] + a[3][1] * ffz[i + 1] + a[3][2] * ffz[i + 2] + a[3][3] * ffz[i + 3];
	// smoothness index
	IS0 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
	IS1 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
	IS2 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);
	IS3 = 1. / 36. * pow(-11. * ffz[i + 0] + 18.0 * ffz[i + 1] - 9. * ffz[i + 2] + 2. * ffz[i + 3], 2) +
		  13. / 12. * pow(2. * ffz[i + 0] - 5.0 * ffz[i + 1] + 4. * ffz[i + 2] - ffz[i + 3], 2) +
		  781. / 720. * pow(-ffz[i + 0] + 3.0 * ffz[i + 1] - 3. * ffz[i + 2] + ffz[i + 3], 2);

	// alpha
	// absTerm0 = abs(IS0 - IS2);
	// = β3 − 1/6(β0 + β2 + 4β1),
	absTerm0 = IS3 - (IS1 + IS0 + 4. * IS2) / 6.;
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);
	alpha3 = 1.0 + pow(absTerm0 / (IS3 + varepsilon), 6);

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;
	omega3 = alpha3 * isumAlpha;

	// delta
	d0 = omega0 < S_ct ? 0. : 1.;
	d1 = omega1 < S_ct ? 0. : 1.;
	d2 = omega2 < S_ct ? 0. : 1.;
	d3 = omega3 < S_ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;
	omega3 = d3 * cw[3] * dsum;

	// WENO plus reconstruction
	fhatp = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3;

	// --- WENO minus part --- !

	// polinomia
	q0 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
	q1 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
	q2 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];
	q3 = a[3][0] * fff[i + 1] + a[3][1] * fff[i + 0] + a[3][2] * fff[i - 1] + a[3][3] * fff[i - 2];

	// smoothness index
	IS0 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
	IS1 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
	IS2 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);
	IS3 = 1. / 36. * pow(-11. * fff[i + 1] + 18.0 * fff[i + 0] - 9. * fff[i - 1] + 2. * fff[i - 2], 2) +
		  13. / 12. * pow(2. * fff[i + 1] - 5.0 * fff[i + 0] + 4. * fff[i - 1] - fff[i - 2], 2) +
		  781. / 720. * pow(-fff[i + 1] + 3.0 * fff[i + 0] - 3. * fff[i - 1] + fff[i - 2], 2);
	// alpha

	absTerm0 = IS3 - (IS1 + IS0 + 4. * IS2) / 6.;
	// absTerm0 = abs(IS0 - IS2);
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);
	alpha3 = 1.0 + pow(absTerm0 / (IS3 + varepsilon), 6);

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;
	omega3 = alpha3 * isumAlpha;

	// delta

	d0 = omega0 < S_ct ? 0. : 1.;
	d1 = omega1 < S_ct ? 0. : 1.;
	d2 = omega2 < S_ct ? 0. : 1.;
	d3 = omega3 < S_ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;
	omega3 = d3 * cw[3] * dsum;

	// WENO minus reconstruction
	fhatm = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3;

	yy = (fhatp + fhatm);
	return yy;
}

mydouble CADSolver::scheme_df_dx_char_Teno7(mydouble *ffz, mydouble *fff, long i)
{
	mydouble yy = 0.0;

	mydouble varepsilon = 1.0e-16;

	// mydouble S_ct = 1.0e-7;
	// mydouble S_ct = 1.0e-3; // High_jet case

	mydouble q0, q1, q2, q3, q4, IS0, IS1, IS2, IS3, IS4, eps, eta0, eta1, eta2, rm;
	mydouble alpha0, alpha1, alpha2, alpha3, alpha4;
	mydouble omega0, omega1, omega2, omega3, omega4;
	mydouble isumAlpha, absTerm0, absTerm1, d0, d1, d2, d3, d4, fhatp, fhatm, dsum, tao5;

	mydouble is1par = 13.0 / 12.0;
	mydouble is2par = 0.25;

	mydouble a[5][4] = {{2.0 / 6.0, -7.0 / 6.0, 11.0 / 6.0, 0},
						{-1.0 / 6.0, 5.0 / 6.0, 2.0 / 6.0, 0},
						{2.0 / 6.0, 5.0 / 6.0, -1.0 / 6.0, 0},
						{3.0 / 12.0, 13.0 / 12.0, -5.0 / 12.0, 1.0 / 12.0},
						{-3.0 / 12.0, 13.0 / 12.0, -23.0 / 12.0, 25.0 / 12.0}};

	// mydouble cw[4] = {0.1, 0.6, 0.3};
	// mydouble cw[4] = {1. / 20., 9. / 20., 6. / 20., 4. / 20.};
	// mydouble cw[5] = {3. / 35., 18. / 35., 9. / 35., 4. / 35., 1. / 35.};
	mydouble cw[5] = {18. / 35., 9. / 35., 3. / 35., 4. / 35., 1. / 35.};
	// polinomia

	// q0 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
	// q1 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
	// q2 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];

	q2 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
	q0 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
	q1 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];

	q3 = a[3][0] * ffz[i + 0] + a[3][1] * ffz[i + 1] + a[3][2] * ffz[i + 2] + a[3][3] * ffz[i + 3];
	q4 = a[4][0] * ffz[i - 3] + a[4][1] * ffz[i - 2] + a[4][2] * ffz[i - 1] + a[4][3] * ffz[i];

	// smoothness index
	// IS0 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
	// IS1 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
	// IS2 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);

	IS2 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
	IS0 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
	IS1 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);

	IS3 = 1. / 36. * pow(-11. * ffz[i + 0] + 18.0 * ffz[i + 1] - 9. * ffz[i + 2] + 2. * ffz[i + 3], 2) +
		  13. / 12. * pow(2. * ffz[i + 0] - 5.0 * ffz[i + 1] + 4. * ffz[i + 2] - ffz[i + 3], 2) +
		  781. / 720. * pow(-ffz[i + 0] + 3.0 * ffz[i + 1] - 3. * ffz[i + 2] + ffz[i + 3], 2);

	IS4 = 1. / 36. * pow(-2. * ffz[i - 3] + 9.0 * ffz[i - 2] - 18. * ffz[i - 1] + 11. * ffz[i + 0], 2) +
		  13. / 12. * pow(-ffz[i - 3] + 4.0 * ffz[i - 2] - 5. * ffz[i - 1] + 2. * ffz[i + 0], 2) +
		  781. / 720. * pow(-ffz[i - 3] + 3.0 * ffz[i - 2] - 3. * ffz[i - 1] + ffz[i + 0], 2);

	// tao=(32154783380.0_8*v(4)**2 + 18133963560.0_8*v(5)**2 + 2927992563.0_8*v(6)**2 -&
	// &969999969.0_8*v(6)*v(7) + 84070496.0_8*v(7)**2 - 12546315963.0_8*v(6)*v(3) +&
	// &1902531828.0_8*v(7)*v(3) + 18133963560.0_8*v(3)**2 + 4550242446.0_8*v(6)*v(2) -&
	// &676871859.0_8*v(7)*v(2) - 14296379553.0_8*v(3)*v(2) + 2927992563.0_8*v(2)**2 -&
	// &676871859.0_8*v(6)*v(1) + 99022657.0_8*v(7)*v(1) + 2283428883.0_8*v(3)*v(1) -&
	// &969999969.0_8*v(2)*v(1) + 84070496.0_8*v(1)**2 +&
	// &3.0_8*v(5)*(-4765459851.0_8*v(6) + 761142961.0_8*v(7) + 11273559435.0_8*v(3) - &
	//    &4182105321.0_8*v(2) + 634177276.0_8*v(1)) - &
	// &4.0_8*v(4)*(11857967655.0_8*v(5) - 4520834943.0_8*v(6) + 701563133.0_8*v(7) + &
	//    &11857967655.0_8*v(3) - 4520834943.0_8*v(2) + 701563133.0_8*v(1)))/59875200.0_8

	// 以上代码是Fortran， v(4)对应ffz[i + 0]，v(5)对应ffz[i + 1]，v(6)对应ffz[i + 2]，v(7)对应ffz[i + 3]，v(3)对应ffz[i -1]， v(2)对应ffz[i - 2]，v(1)对应ffz[i - 3]
	// 把tao 转化为 tao5， v(4)对应ffz[i + 0]，v(5)对应ffz[i + 1]，v(6)对应ffz[i + 2]，v(7)对应ffz[i + 3]，v(3)对应ffz[i -1]， v(2)对应ffz[i - 2]，v(1)对应ffz[i - 3]
	tao5 = (32154783380.0 * pow(ffz[i + 0], 2) + 18133963560.0 * pow(ffz[i + 1], 2) + 2927992563.0 * pow(ffz[i + 2], 2) -
			969999969.0 * ffz[i + 2] * ffz[i + 3] + 84070496.0 * pow(ffz[i + 3], 2) - 12546315963.0 * ffz[i + 2] * ffz[i - 1] +
			1902531828.0 * ffz[i + 3] * ffz[i - 1] + 18133963560.0 * pow(ffz[i - 1], 2) + 4550242446.0 * ffz[i + 2] * ffz[i - 2] -
			676871859.0 * ffz[i + 3] * ffz[i - 2] - 14296379553.0 * ffz[i - 1] * ffz[i - 2] + 2927992563.0 * pow(ffz[i - 2], 2) -
			676871859.0 * ffz[i + 2] * ffz[i - 3] + 99022657.0 * ffz[i + 3] * ffz[i - 3] + 2283428883.0 * ffz[i - 1] * ffz[i - 3] -
			969999969.0 * ffz[i - 2] * ffz[i - 3] + 84070496.0 * pow(ffz[i - 3], 2) +
			3.0 * ffz[i + 1] * (-4765459851.0 * ffz[i + 2] + 761142961.0 * ffz[i + 3] + 11273559435.0 * ffz[i - 1] - 4182105321.0 * ffz[i - 2] + 634177276.0 * ffz[i - 3]) -
			4.0 * ffz[i + 0] * (11857967655.0 * ffz[i + 1] - 4520834943.0 * ffz[i + 2] + 701563133.0 * ffz[i + 3] + 11857967655.0 * ffz[i - 1] - 4520834943.0 * ffz[i - 2] + 701563133.0 * ffz[i - 3])) /
		   59875200.0;

	// alpha
	// absTerm0 = abs(IS0 - IS2);
	// = β3 − 1/6(β0 + β2 + 4β1),
	// absTerm0 = abs(IS3 - (IS1 + IS0 + 4. * IS2) / 6.);
	absTerm0 = abs(tao5 - (IS1 + IS2 + 4. * IS0) / 6.);

	// absTerm1 = abs(IS4 - (IS1 + IS0 + 4. * IS2) / 6.);
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);
	alpha3 = 1.0 + pow(absTerm0 / (IS3 + varepsilon), 6);
	alpha4 = 1.0 + pow(absTerm0 / (IS4 + varepsilon), 6);
	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3 + alpha4);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;
	omega3 = alpha3 * isumAlpha;
	omega4 = alpha4 * isumAlpha;

	// delta
	d0 = omega0 < S_ct ? 0. : 1.;
	d1 = omega1 < S_ct ? 0. : 1.;
	d2 = omega2 < S_ct ? 0. : 1.;
	d3 = omega3 < S_ct ? 0. : 1.;
	d4 = omega4 < S_ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3] + d4 * cw[4]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;
	omega3 = d3 * cw[3] * dsum;
	omega4 = d4 * cw[4] * dsum;

	// WENO plus reconstruction
	fhatp = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3 + omega4 * q4;

	// --- WENO minus part --- !

	// polinomia
	// q0 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
	// q1 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
	// q2 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];

	q2 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
	q0 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
	q1 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];
	q3 = a[3][0] * fff[i + 1] + a[3][1] * fff[i + 0] + a[3][2] * fff[i - 1] + a[3][3] * fff[i - 2];
	q4 = a[4][0] * fff[i + 4] + a[4][1] * fff[i + 3] + a[4][2] * fff[i + 2] + a[4][3] * fff[i + 1];

	// smoothness index
	// IS0 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
	// IS1 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
	// IS2 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);

	IS2 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
	IS0 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
	IS1 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);

	IS3 = 1. / 36. * pow(-11. * fff[i + 1] + 18.0 * fff[i + 0] - 9. * fff[i - 1] + 2. * fff[i - 2], 2) +
		  13. / 12. * pow(2. * fff[i + 1] - 5.0 * fff[i + 0] + 4. * fff[i - 1] - fff[i - 2], 2) +
		  781. / 720. * pow(-fff[i + 1] + 3.0 * fff[i + 0] - 3. * fff[i - 1] + fff[i - 2], 2);
	IS4 = 1. / 36. * pow(-2. * fff[i + 4] + 9.0 * fff[i + 3] - 18. * fff[i + 2] + 11. * fff[i + 1], 2) +
		  13. / 12. * pow(-fff[i + 4] + 4.0 * fff[i + 3] - 5. * fff[i + 2] + 2. * fff[i + 1], 2) +
		  781. / 720. * pow(-fff[i + 4] + 3.0 * fff[i + 3] - 3. * fff[i + 2] + fff[i + 1], 2);
	// alpha
	// tao5= (32154783380.0 * pow(ffz[i + 0],2) + 18133963560.0 * pow(ffz[i + 1],2) + 2927992563.0 * pow(ffz[i + 2],2) -
	// 969999969.0 * ffz[i + 2] * ffz[i + 3] + 84070496.0 * pow(ffz[i + 3],2) - 12546315963.0 * ffz[i + 2] * ffz[i - 1] +
	// 1902531828.0 * ffz[i + 3] * ffz[i - 1] + 18133963560.0 * pow(ffz[i - 1],2) + 4550242446.0 * ffz[i + 2] * ffz[i - 2] -
	// 676871859.0 * ffz[i + 3] * ffz[i - 2] - 14296379553.0 * ffz[i - 1] * ffz[i - 2] + 2927992563.0 * pow(ffz[i - 2],2) -
	// 676871859.0 * ffz[i + 2] * ffz[i - 3] + 99022657.0 * ffz[i + 3] * ffz[i - 3] + 2283428883.0 * ffz[i - 1] * ffz[i - 3] -
	// 969999969.0 * ffz[i - 2] * ffz[i - 3] + 84070496.0 * pow(ffz[i - 3],2) +
	// 3.0 * ffz[i + 1] * (-4765459851.0 * ffz[i + 2] + 761142961.0 * ffz[i + 3] + 11273559435.0 * ffz[i - 1] -
	//     4182105321.0 * ffz[i - 2] + 634177276.0 * ffz[i - 3]) -
	//  4.0 * ffz[i + 0] * (11857967655.0 * ffz[i + 1] - 4520834943.0 * ffz[i + 2] + 701563133.0 * ffz[i + 3] +
	//     11857967655.0 * ffz[i - 1] - 4520834943.0 * ffz[i - 2] + 701563133.0 * ffz[i - 3]))/59875200.0;

	// fff和ffz关于i+0.5对称  ffz[i + 4] 对应 fff[i - 3]  ffz[i + 3] 对应 fff[i - 2]  ffz[i + 2] 对应 fff[i - 1]  ffz[i + 1] 对应 fff[i + 0]  ffz[i + 0] 对应 fff[i + 1]  ffz[i - 1] 对应 fff[i + 2]  ffz[i - 2] 对应 fff[i + 3]  ffz[i - 3] 对应 fff[i + 4]
	tao5 = (32154783380.0 * pow(fff[i + 1], 2) + 18133963560.0 * pow(fff[i + 0], 2) + 2927992563.0 * pow(fff[i - 1], 2) -
			969999969.0 * fff[i - 1] * fff[i - 2] + 84070496.0 * pow(fff[i - 2], 2) - 12546315963.0 * fff[i - 1] * fff[i + 2] +
			1902531828.0 * fff[i - 2] * fff[i + 2] + 18133963560.0 * pow(fff[i + 2], 2) + 4550242446.0 * fff[i - 1] * fff[i + 3] -
			676871859.0 * fff[i - 2] * fff[i + 3] - 14296379553.0 * fff[i + 2] * fff[i + 3] + 2927992563.0 * pow(fff[i + 3], 2) -
			676871859.0 * fff[i - 1] * fff[i + 4] + 99022657.0 * fff[i - 2] * fff[i + 4] + 2283428883.0 * fff[i + 2] * fff[i + 4] -
			969999969.0 * fff[i + 3] * fff[i + 4] + 84070496.0 * pow(fff[i + 4], 2) +
			3.0 * fff[i + 0] * (-4765459851.0 * fff[i - 1] + 761142961.0 * fff[i - 2] + 11273559435.0 * fff[i + 2] - 4182105321.0 * fff[i + 3] + 634177276.0 * fff[i + 4]) -
			4.0 * fff[i + 1] * (11857967655.0 * fff[i + 0] - 4520834943.0 * fff[i - 1] + 701563133.0 * fff[i - 2] + 11857967655.0 * fff[i + 2] - 4520834943.0 * fff[i + 3] + 701563133.0 * fff[i + 4])) /
		   59875200.0;

	// absTerm0 = abs(IS3 - (IS1 + IS0 + 4. * IS2) / 6.);
	// absTerm1 = abs(IS4 - (IS1 + IS0 + 4. * IS2) / 6.);

	absTerm0 = abs(tao5 - (IS1 + IS2 + 4. * IS0) / 6.);
	// absTerm0 = abs(IS0 - IS2);
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);
	alpha3 = 1.0 + pow(absTerm0 / (IS3 + varepsilon), 6);
	alpha4 = 1.0 + pow(absTerm0 / (IS4 + varepsilon), 6);

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3 + alpha4);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;
	omega3 = alpha3 * isumAlpha;
	omega4 = alpha4 * isumAlpha;

	// delta

	d0 = omega0 < S_ct ? 0. : 1.;
	d1 = omega1 < S_ct ? 0. : 1.;
	d2 = omega2 < S_ct ? 0. : 1.;
	d3 = omega3 < S_ct ? 0. : 1.;
	d4 = omega4 < S_ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3] + d4 * cw[4]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;
	omega3 = d3 * cw[3] * dsum;
	omega4 = d4 * cw[4] * dsum;

	// WENO minus reconstruction
	fhatm = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3 + omega4 * q4;

	yy = (fhatp + fhatm);
	return yy;
}

mydouble CADSolver::scheme_df_dx_char_Teno8(mydouble *ffz, mydouble *fff, long i)
{
	mydouble yy = 0.0;

	mydouble varepsilon = 1.0e-40;

	// mydouble S_ct = 1.0e-7;
	// mydouble S_ct = 1.0e-3; // High_jet case

	mydouble q0, q1, q2, q3, q4, q5, IS0, IS1, IS2, IS3, IS4, IS5, eps, eta0, eta1, eta2, rm;
	mydouble alpha0, alpha1, alpha2, alpha3, alpha4, alpha5;
	mydouble omega0, omega1, omega2, omega3, omega4, omega5;
	mydouble isumAlpha, absTerm0, absTerm1, d0, d1, d2, d3, d4, d5, fhatp, fhatm, dsum, tau8;

	mydouble is1par = 13.0 / 12.0;
	mydouble is2par = 0.25;

	mydouble a[6][5] = {{2.0 / 6.0, -7.0 / 6.0, 11.0 / 6.0, 0, 0},
						{-1.0 / 6.0, 5.0 / 6.0, 2.0 / 6.0, 0, 0},
						{2.0 / 6.0, 5.0 / 6.0, -1.0 / 6.0, 0, 0},
						{3.0 / 12.0, 13.0 / 12.0, -5.0 / 12.0, 1.0 / 12.0, 0},
						{-3.0 / 12.0, 13.0 / 12.0, -23.0 / 12.0, 25.0 / 12.0, 0},
						{12.0 / 60.0, 77.0 / 60.0, -43.0 / 60.0, 17.0 / 60.0, -3.0 / 60.0}};

	// mydouble cw[4] = {0.1, 0.6, 0.3};
	// mydouble cw[4] = {1. / 20., 9. / 20., 6. / 20., 4. / 20.};
	// mydouble cw[5] = {3. / 35., 18. / 35., 9. / 35., 4. / 35., 1. / 35.};
	mydouble cw[6] = {4. / 70., 30. / 70., 18. / 70., 12. / 70., 1. / 70., 5. / 70.};
	// polinomia

	q0 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
	q1 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
	q2 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];
	q3 = a[3][0] * ffz[i + 0] + a[3][1] * ffz[i + 1] + a[3][2] * ffz[i + 2] + a[3][3] * ffz[i + 3];
	q4 = a[4][0] * ffz[i - 3] + a[4][1] * ffz[i - 2] + a[4][2] * ffz[i - 1] + a[4][3] * ffz[i];
	q5 = a[5][0] * ffz[i + 0] + a[5][1] * ffz[i + 1] + a[5][2] * ffz[i + 2] + a[5][3] * ffz[i + 3] + a[5][4] * ffz[i + 4];
	// smoothness index
	IS0 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
	IS1 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
	IS2 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);
	IS3 = 1. / 36. * pow(-11. * ffz[i + 0] + 18.0 * ffz[i + 1] - 9. * ffz[i + 2] + 2. * ffz[i + 3], 2) +
		  13. / 12. * pow(2. * ffz[i + 0] - 5.0 * ffz[i + 1] + 4. * ffz[i + 2] - ffz[i + 3], 2) +
		  781. / 720. * pow(-ffz[i + 0] + 3.0 * ffz[i + 1] - 3. * ffz[i + 2] + ffz[i + 3], 2);

	IS4 = 1. / 36. * pow(-2. * ffz[i - 3] + 9.0 * ffz[i - 2] - 18. * ffz[i - 1] + 11. * ffz[i + 0], 2) +
		  13. / 12. * pow(-ffz[i - 3] + 4.0 * ffz[i - 2] - 5. * ffz[i - 1] + 2. * ffz[i + 0], 2) +
		  781. / 720. * pow(-ffz[i - 3] + 3.0 * ffz[i - 2] - 3. * ffz[i - 1] + ffz[i + 0], 2);

	IS5 = 1. / 144. * pow(-25. * ffz[i + 0] + 48.0 * ffz[i + 1] - 36. * ffz[i + 2] + 16. * ffz[i + 3] - 3. * ffz[i + 4], 2) +
		  13. / 1728. * pow(35. * ffz[i + 0] - 104.0 * ffz[i + 1] + 114. * ffz[i + 2] - 56. * ffz[i + 3] + 11. * ffz[i + 4], 2) +
		  781. / 2880. * pow(-5. * ffz[i + 0] + 18.0 * ffz[i + 1] - 24. * ffz[i + 2] + 14. * ffz[i + 3] - 3. * ffz[i + 4], 2) -
		  1. / 4320. * (35. * ffz[i + 0] - 104.0 * ffz[i + 1] + 114. * ffz[i + 2] - 56. * ffz[i + 3] + 11. * ffz[i + 4]) *
			  (ffz[i + 0] - 4. * ffz[i + 1] + 6. * ffz[i + 2] - 4. * ffz[i + 3] + ffz[i + 4]) +
		  32803. / 30240. * pow(ffz[i + 0] - 4. * ffz[i + 1] + 6. * ffz[i + 2] - 4. * ffz[i + 3] + ffz[i + 4], 2);

	tau8 = 1. / 62270208000. * 75349098471. * ffz[i + 4] - 1078504915264. * ffz[i + 3] + 3263178215782. * ffz[i + 2] - 5401061230160. * ffz[i + 1] + 5274436892970. * ffz[i] - 3038037798592. * ffz[i - 1] + 956371298594. * ffz[i - 2] - 127080660272. * ffz[i - 3] + ffz[i + 3] * (3944861897609. * ffz[i + 3] - 24347015748304. * ffz[i + 2] + 41008808432890. * ffz[i + 1] - 40666174667520. * ffz[i] + 23740865961334. * ffz[i - 1] - 7563868580208. * ffz[i - 2] + 1016165721854. * ffz[i - 3]) + ffz[i + 2] * (38329064547231. * ffz[i + 2] - 131672853704480. * ffz[i + 1] + 132979856899250. * ffz[i] - 78915800051952. * ffz[i - 1] + 25505661974314. * ffz[i - 2] - 3471156679072. * ffz[i - 3]) + ffz[i + 1] * (115451981835025. * ffz[i + 1] - 238079153652400. * ffz[i] + 144094750348910. * ffz[i - 1] - 47407534412640. * ffz[i - 2] + 6553080547830. * ffz[i - 3]) + ffz[i] * (125494539510175. * ffz[i] - 155373333547520. * ffz[i - 1] + 52241614797670. * ffz[i - 2] - 7366325742800. * ffz[i - 3]) + ffz[i - 1] * (49287325751121. * ffz[i - 1] - 33999931981264. * ffz[i - 2] + 4916835566842. * ffz[i - 3]) + ffz[i - 2] * (6033767706599. * ffz[i - 2] - 1799848509664. * ffz[i - 3]) + 139164877641. * ffz[i - 3] * ffz[i - 3];
	// tau8 = 1. / 62270208000. *   75349098471. * ffz[i + 4] - 1078504915264. * ffz[i + 3] + 3263178215782. * ffz[i + 2]
	// - 5401061230160. * ffz[i + 1] + 5274436892970. * ffz[i] - 3038037798592. * ffz[i - 1]
	// + 956371298594. * ffz[i - 2] - 127080660272. * ffz[i - 3]  + ffz[i + 3]* (3944861897609. * ffz[i + 3]
	// - 24347015748304. * ffz[i + 2] + 41008808432890. * ffz[i + 1] - 40666174667520. * ffz[i]
	// + 23740865961334. * ffz[i - 1] - 7563868580208. * ffz[i - 2] + 1016165721854. * ffz[i - 3])
	// + ffz[i + 2]*(38329064547231. * ffz[i + 2] - 131672853704480. * ffz[i + 1] + 132979856899250. * ffz[i]
	// - 78915800051952. * ffz[i - 1] + 25505661974314. * ffz[i - 2] - 3471156679072. * ffz[i - 3])
	// + ffz[i + 1] *(115451981835025. * ffz[i + 1] - 238079153652400. * ffz[i] + 144094750348910. * ffz[i - 1]
	// - 47407534412640. * ffz[i - 2] + 6553080547830. * ffz[i - 3]) + ffz[i]*(125494539510175. * ffz[i]
	// - 155373333547520. * ffz[i - 1] + 52241614797670. * ffz[i - 2] - 7366325742800. * ffz[i - 3])
	// + ffz[i - 1] *(49287325751121. * ffz[i - 1] - 33999931981264. * ffz[i - 2] + 4916835566842. * ffz[i - 3])
	// + ffz[i - 2]*(6033767706599. * ffz[i - 2] - 1799848509664. * ffz[i - 3]) + 139164877641. * ffz[i - 3] * ffz[i - 3];
	// // alpha
	// absTerm0 = abs(IS0 - IS2);
	// = β3 − 1/6(β0 + β2 + 4β1),
	tau8 = abs(tau8);
	absTerm0 = abs(tau8 - (IS1 + IS0 + 4. * IS2) / 6.);
	// absTerm0 = tau8;
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);
	alpha3 = 1.0 + pow(absTerm0 / (IS3 + varepsilon), 6);
	alpha4 = 1.0 + pow(absTerm0 / (IS4 + varepsilon), 6);
	alpha5 = 1.0 + pow(absTerm0 / (IS5 + varepsilon), 6);
	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;
	omega3 = alpha3 * isumAlpha;
	omega4 = alpha4 * isumAlpha;
	omega5 = alpha5 * isumAlpha;

	// delta
	d0 = omega0 < S_ct ? 0. : 1.;
	d1 = omega1 < S_ct ? 0. : 1.;
	d2 = omega2 < S_ct ? 0. : 1.;
	d3 = omega3 < S_ct ? 0. : 1.;
	d4 = omega4 < S_ct ? 0. : 1.;
	d5 = omega5 < S_ct ? 0. : 1.;

	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3] + d4 * cw[4] + d5 * cw[5]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;
	omega3 = d3 * cw[3] * dsum;
	omega4 = d4 * cw[4] * dsum;
	omega5 = d5 * cw[5] * dsum;

	// WENO plus reconstruction
	// fhatp = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3 + omega4 * q4;
	fhatp = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3 + omega4 * q4 + omega5 * q5;

	// --- WENO minus part --- !

	// polinomia
	q0 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
	q1 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
	q2 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];
	q3 = a[3][0] * fff[i + 1] + a[3][1] * fff[i + 0] + a[3][2] * fff[i - 1] + a[3][3] * fff[i - 2];
	q4 = a[4][0] * fff[i + 4] + a[4][1] * fff[i + 3] + a[4][2] * fff[i + 2] + a[4][3] * fff[i + 1];
	q5 = a[5][0] * fff[i + 1] + a[5][1] * fff[i + 0] + a[5][2] * fff[i - 1] + a[5][3] * fff[i - 2] + a[5][4] * fff[i - 3];

	// smoothness index
	IS0 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
	IS1 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
	IS2 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);
	IS3 = 1. / 36. * pow(-11. * fff[i + 1] + 18.0 * fff[i + 0] - 9. * fff[i - 1] + 2. * fff[i - 2], 2) +
		  13. / 12. * pow(2. * fff[i + 1] - 5.0 * fff[i + 0] + 4. * fff[i - 1] - fff[i - 2], 2) +
		  781. / 720. * pow(-fff[i + 1] + 3.0 * fff[i + 0] - 3. * fff[i - 1] + fff[i - 2], 2);
	IS4 = 1. / 36. * pow(-2. * fff[i + 4] + 9.0 * fff[i + 3] - 18. * fff[i + 2] + 11. * fff[i + 1], 2) +
		  13. / 12. * pow(-fff[i + 4] + 4.0 * fff[i + 3] - 5. * fff[i + 2] + 2. * fff[i + 1], 2) +
		  781. / 720. * pow(-fff[i + 4] + 3.0 * fff[i + 3] - 3. * fff[i + 2] + fff[i + 1], 2);
	IS5 = 1. / 144. * pow(-25. * fff[i + 1] + 48.0 * fff[i + 0] - 36. * fff[i - 1] + 16. * fff[i - 2] - 3. * fff[i - 3], 2) +
		  13. / 1728. * pow(35. * fff[i + 1] - 104.0 * fff[i + 0] + 114. * fff[i - 1] - 56. * fff[i - 2] + 11. * fff[i - 3], 2) +
		  781. / 2880. * pow(-5. * fff[i + 1] + 18.0 * fff[i + 0] - 24. * fff[i - 1] + 14. * fff[i - 2] - 3. * fff[i - 3], 2) -
		  1. / 4320. * (35. * fff[i + 1] - 104.0 * fff[i + 0] + 114. * fff[i - 1] - 56. * fff[i - 2] + 11. * fff[i - 3]) *
			  (fff[i + 1] - 4. * fff[i + 0] + 6. * fff[i - 1] - 4. * fff[i - 2] + fff[i - 3]) +
		  32803. / 30240. * pow(fff[i + 1] - 4. * fff[i + 0] + 6. * fff[i - 1] - 4. * fff[i - 2] + fff[i - 3], 2);
	// alpha
	// fff和ffz关于i+0.5对称  ffz[i + 4] 对应 fff[i - 3]  ffz[i + 3] 对应 fff[i - 2]  ffz[i + 2] 对应 fff[i - 1]  ffz[i + 1] 对应 fff[i + 0]  ffz[i + 0] 对应 fff[i + 1]  ffz[i - 1] 对应 fff[i + 2]  ffz[i - 2] 对应 fff[i + 3]  ffz[i - 3] 对应 fff[i + 4]
	tau8 = 1. / 62270208000. * 75349098471. * fff[i - 3] - 1078504915264. * fff[i - 2] + 3263178215782. * fff[i - 1] - 5401061230160. * fff[i + 0] + 5274436892970. * fff[i + 1] - 3038037798592. * fff[i + 2] + 956371298594. * fff[i + 3] - 127080660272. * fff[i + 4] + fff[i - 2] * (3944861897609. * fff[i - 2] - 24347015748304. * fff[i - 1] + 41008808432890. * fff[i + 0] - 40666174667520. * fff[i + 1] + 23740865961334. * fff[i + 2] - 7563868580208. * fff[i + 3] + 1016165721854. * fff[i + 4]) + fff[i - 1] * (38329064547231. * fff[i - 1] - 131672853704480. * fff[i + 0] + 132979856899250. * fff[i + 1] - 78915800051952. * fff[i + 2] + 25505661974314. * fff[i + 3] - 3471156679072. * fff[i + 4]) + fff[i + 0] * (115451981835025. * fff[i + 0] - 238079153652400. * fff[i + 1] + 144094750348910. * fff[i + 2] - 47407534412640. * fff[i + 3] + 6553080547830. * fff[i + 4]) + fff[i + 1] * (125494539510175. * fff[i + 1] - 155373333547520. * fff[i + 2] + 52241614797670. * fff[i + 3] - 7366325742800. * fff[i + 4]) + fff[i + 2] * (49287325751121. * fff[i + 2] - 33999931981264. * fff[i + 3] + 4916835566842. * fff[i + 4]) + fff[i + 3] * (6033767706599. * fff[i + 3] - 1799848509664. * fff[i + 4]) + 139164877641. * fff[i + 4] * fff[i + 4];

	tau8 = abs(tau8);
	absTerm0 = abs(tau8 - (IS1 + IS0 + 4. * IS2) / 6.);
	// absTerm1 = abs(IS4 - (IS1 + IS0 + 4. * IS2) / 6.);
	// absTerm0 = tau8;
	alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
	alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
	alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);
	alpha3 = 1.0 + pow(absTerm0 / (IS3 + varepsilon), 6);
	alpha4 = 1.0 + pow(absTerm0 / (IS4 + varepsilon), 6);
	alpha5 = 1.0 + pow(absTerm0 / (IS5 + varepsilon), 6);

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;
	omega3 = alpha3 * isumAlpha;
	omega4 = alpha4 * isumAlpha;
	omega5 = alpha5 * isumAlpha;

	// delta

	d0 = omega0 < S_ct ? 0. : 1.;
	d1 = omega1 < S_ct ? 0. : 1.;
	d2 = omega2 < S_ct ? 0. : 1.;
	d3 = omega3 < S_ct ? 0. : 1.;
	d4 = omega4 < S_ct ? 0. : 1.;
	d5 = omega5 < S_ct ? 0. : 1.;

	// dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3] + d4 * cw[4]);
	dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3] + d4 * cw[4] + d5 * cw[5]);

	omega0 = d0 * cw[0] * dsum;
	omega1 = d1 * cw[1] * dsum;
	omega2 = d2 * cw[2] * dsum;
	omega3 = d3 * cw[3] * dsum;
	omega4 = d4 * cw[4] * dsum;
	omega5 = d5 * cw[5] * dsum;

	// WENO minus reconstruction
	// fhatm = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3 + omega4 * q4;
	fhatm = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3 + omega4 * q4 + omega5 * q5;

	yy = (fhatp + fhatm);
	return yy;
}

mydouble CADSolver::scheme_df_dx_char_weno7Z(mydouble *ffz, mydouble *fff, long i)
{
	double yy = 0.0, ss1, ss2;

	//         /*   weno7-Z================================
	double S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33,
		ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, tau;
	double CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
		   ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
		   ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
		   ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
		   ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
		   b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
		   b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
		   d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
	double e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
		   e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
		   e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
		   e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-30; //   !! WENO-JS

	// ! 7th order WENO scheme
	// ! 1  阶导数      *(ffz + i - 3)
	S10 = ax11 * *(ffz + i - 3) + ax12 * *(ffz + i - 2) + ax13 * *(ffz + i - 1) + ax14 * *(ffz + i);
	S11 = ax21 * *(ffz + i - 2) - *(ffz + i - 1) + ax23 * *(ffz + i) + ax24 * *(ffz + i + 1);
	S12 = ax31 * *(ffz + i - 1) + ax32 * *(ffz + i) + *(ffz + i + 1) + ax34 * *(ffz + i + 2);
	S13 = ax41 * *(ffz + i) + ax42 * *(ffz + i + 1) + ax43 * *(ffz + i + 2) + ax44 * *(ffz + i + 3);
	//  ! 2 阶导数
	S20 = -*(ffz + i - 3) + b12 * *(ffz + i - 2) + b13 * *(ffz + i - 1) + b14 * *(ffz + i);
	S21 = *(ffz + i - 1) + b22 * *(ffz + i) + *(ffz + i + 1);
	S22 = *(ffz + i) + b22 * *(ffz + i + 1) + *(ffz + i + 2);
	S23 = b41 * *(ffz + i) + b42 * *(ffz + i + 1) + b43 * *(ffz + i + 2) - *(ffz + i + 3);
	// ! 3 阶导数
	S30 = -*(ffz + i - 3) + c12 * (*(ffz + i - 2) - *(ffz + i - 1)) + *(ffz + i);
	S31 = -*(ffz + i - 2) + c12 * (*(ffz + i - 1) - *(ffz + i)) + *(ffz + i + 1);
	S32 = -*(ffz + i - 1) + c12 * (*(ffz + i) - *(ffz + i + 1)) + *(ffz + i + 2);
	S33 = -*(ffz + i) + c12 * (*(ffz + i + 1) - *(ffz + i + 2)) + *(ffz + i + 3);

	S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
	S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
	S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
	S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
	// tau = abs(S0 - S1 - S2 + S3);
	tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
	// !-------WENO Z----------------------
	ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
	ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
	ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
	ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
	// !-------WENO J-S----------------------
	// ax0 = CC0 / ((ep + S0) * (ep + S0));
	// ax1 = CC1 / ((ep + S1) * (ep + S1));
	// ax2 = CC2 / ((ep + S2) * (ep + S2));
	// ax3 = CC3 / ((ep + S3) * (ep + S3));
	// !-----------------------------------------------
	am = ax0 + ax1 + ax2 + ax3;
	// *(flagh + i) = abs(ax0 / am - CC0) + abs(ax1 / am - CC1) + abs(ax2 / am - CC2) + abs(ax3 / am - CC3);
	// !  4阶差分格式的通量
	q0 = e11 * *(ffz + i - 3) + e12 * *(ffz + i - 2) + e13 * *(ffz + i - 1) + e14 * *(ffz + i);
	q1 = e21 * *(ffz + i - 2) + e22 * *(ffz + i - 1) + e23 * *(ffz + i) + e24 * *(ffz + i + 1);
	q2 = e31 * *(ffz + i - 1) + e32 * *(ffz + i) + e33 * *(ffz + i + 1) + e34 * *(ffz + i + 2);
	q3 = e41 * *(ffz + i) + e42 * *(ffz + i + 1) + e43 * *(ffz + i + 2) + e44 * *(ffz + i + 3);
	// !  由4个4阶差分格式组合成1个7阶差分格式
	// !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
	ss1 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

	// -------------------------------
	// !      7th order WENO scheme
	// ! 1  阶导数   *(fff + i - 2)
	S10 = ax11 * *(fff + i + 4) + ax12 * *(fff + i + 3) + ax13 * *(fff + i + 2) + ax14 * *(fff + i + 1);
	S11 = ax21 * *(fff + i + 3) - *(fff + i + 2) + ax23 * *(fff + i + 1) + ax24 * *(fff + i);
	S12 = ax31 * *(fff + i + 2) + ax32 * *(fff + i + 1) + *(fff + i) + ax34 * *(fff + i - 1);
	S13 = ax41 * *(fff + i + 1) + ax42 * *(fff + i) + ax43 * *(fff + i - 1) + ax44 * *(fff + i - 2);
	// ! 2 阶导数
	S20 = -*(fff + i + 4) + b12 * *(fff + i + 3) + b13 * *(fff + i + 2) + b14 * *(fff + i + 1);
	S21 = *(fff + i + 2) + b22 * *(fff + i + 1) + *(fff + i);
	S22 = *(fff + i + 1) + b22 * *(fff + i) + *(fff + i - 1);
	S23 = b41 * *(fff + i + 1) + b42 * *(fff + i) + b43 * *(fff + i - 1) - *(fff + i - 2);
	// ! 3 阶导数
	S30 = -*(fff + i + 4) + c12 * (*(fff + i + 3) - *(fff + i + 2)) + *(fff + i + 1);
	S31 = -*(fff + i + 3) + c12 * (*(fff + i + 2) - *(fff + i + 1)) + *(fff + i);
	S32 = -*(fff + i + 2) + c12 * (*(fff + i + 1) - *(fff + i)) + *(fff + i - 1);
	S33 = -*(fff + i + 1) + c12 * (*(fff + i) - *(fff + i - 1)) + *(fff + i - 2);

	S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
	S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
	S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
	S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
	// tau = abs(S0 - S1 - S2 + S3);
	tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
	// !-------WENO Z----------------------
	ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
	ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
	ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
	ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
	// !-------WENO J-S----------------------
	// ax0 = CC0 / ((ep + S0) * (ep + S0));
	// ax1 = CC1 / ((ep + S1) * (ep + S1));
	// ax2 = CC2 / ((ep + S2) * (ep + S2));
	// ax3 = CC3 / ((ep + S3) * (ep + S3));
	// !-----------------------------------------------
	am = ax0 + ax1 + ax2 + ax3;
	// !  4阶差分格式的通量
	q0 = e11 * *(fff + i + 4) + e12 * *(fff + i + 3) + e13 * *(fff + i + 2) + e14 * *(fff + i + 1);
	q1 = e21 * *(fff + i + 3) + e22 * *(fff + i + 2) + e23 * *(fff + i + 1) + e24 * *(fff + i);
	q2 = e31 * *(fff + i + 2) + e32 * *(fff + i + 1) + e33 * *(fff + i) + e34 * *(fff + i - 1);
	q3 = e41 * *(fff + i + 1) + e42 * *(fff + i) + e43 * *(fff + i - 1) + e44 * *(fff + i - 2);
	// !  由4个4阶差分格式组合成1个7阶差分格式
	ss2 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;
	yy = ss1 + ss2;
	//   */

	return yy;
}

mydouble CADSolver::scheme_df_dx_char_weno7Z_Hy(mydouble *ffz, mydouble *fff, long i, long iCell, long jCell, long nn)
{
	double yy = 0.0, ss1, ss2;
	//         /*   weno7-Z================================
	double S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33,
		ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, tau;
	double CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
		   ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
		   ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
		   ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
		   ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
		   b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
		   b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
		   d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
	double e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
		   e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
		   e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
		   e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-6; //   !! WENO-JS

	// ! 7th order WENO scheme
	// ! 1  阶导数      *(ffz + i - 3)
	S10 = ax11 * *(ffz + i - 3) + ax12 * *(ffz + i - 2) + ax13 * *(ffz + i - 1) + ax14 * *(ffz + i);
	S11 = ax21 * *(ffz + i - 2) - *(ffz + i - 1) + ax23 * *(ffz + i) + ax24 * *(ffz + i + 1);
	S12 = ax31 * *(ffz + i - 1) + ax32 * *(ffz + i) + *(ffz + i + 1) + ax34 * *(ffz + i + 2);
	S13 = ax41 * *(ffz + i) + ax42 * *(ffz + i + 1) + ax43 * *(ffz + i + 2) + ax44 * *(ffz + i + 3);
	//  ! 2 阶导数
	S20 = -*(ffz + i - 3) + b12 * *(ffz + i - 2) + b13 * *(ffz + i - 1) + b14 * *(ffz + i);
	S21 = *(ffz + i - 1) + b22 * *(ffz + i) + *(ffz + i + 1);
	S22 = *(ffz + i) + b22 * *(ffz + i + 1) + *(ffz + i + 2);
	S23 = b41 * *(ffz + i) + b42 * *(ffz + i + 1) + b43 * *(ffz + i + 2) - *(ffz + i + 3);
	// ! 3 阶导数
	S30 = -*(ffz + i - 3) + c12 * (*(ffz + i - 2) - *(ffz + i - 1)) + *(ffz + i);
	S31 = -*(ffz + i - 2) + c12 * (*(ffz + i - 1) - *(ffz + i)) + *(ffz + i + 1);
	S32 = -*(ffz + i - 1) + c12 * (*(ffz + i) - *(ffz + i + 1)) + *(ffz + i + 2);
	S33 = -*(ffz + i) + c12 * (*(ffz + i + 1) - *(ffz + i + 2)) + *(ffz + i + 3);

	S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
	S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
	S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
	S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
	// tau = abs(S0 - S1 - S2 + S3);
	tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
	// !-------WENO Z----------------------
	ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
	ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
	ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
	ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
	// !-------WENO J-S----------------------
	// ax0 = CC0 / ((ep + S0) * (ep + S0));
	// ax1 = CC1 / ((ep + S1) * (ep + S1));
	// ax2 = CC2 / ((ep + S2) * (ep + S2));
	// ax3 = CC3 / ((ep + S3) * (ep + S3));
	// !-----------------------------------------------
	am = ax0 + ax1 + ax2 + ax3;
	// flagh_hy[iCell][jCell][5] = 0.;
	flagh_hy[iCell][jCell][nn] = 0.;
	flagh_hy[iCell][jCell][nn] = abs(ax0 / am - CC0) + abs(ax1 / am - CC1) + abs(ax2 / am - CC2) + abs(ax3 / am - CC3);
	// !  4阶差分格式的通量
	q0 = e11 * *(ffz + i - 3) + e12 * *(ffz + i - 2) + e13 * *(ffz + i - 1) + e14 * *(ffz + i);
	q1 = e21 * *(ffz + i - 2) + e22 * *(ffz + i - 1) + e23 * *(ffz + i) + e24 * *(ffz + i + 1);
	q2 = e31 * *(ffz + i - 1) + e32 * *(ffz + i) + e33 * *(ffz + i + 1) + e34 * *(ffz + i + 2);
	q3 = e41 * *(ffz + i) + e42 * *(ffz + i + 1) + e43 * *(ffz + i + 2) + e44 * *(ffz + i + 3);
	// !  由4个4阶差分格式组合成1个7阶差分格式
	// !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
	ss1 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

	// -------------------------------
	// !      7th order WENO scheme
	// ! 1  阶导数   *(fff + i - 2)
	S10 = ax11 * *(fff + i + 4) + ax12 * *(fff + i + 3) + ax13 * *(fff + i + 2) + ax14 * *(fff + i + 1);
	S11 = ax21 * *(fff + i + 3) - *(fff + i + 2) + ax23 * *(fff + i + 1) + ax24 * *(fff + i);
	S12 = ax31 * *(fff + i + 2) + ax32 * *(fff + i + 1) + *(fff + i) + ax34 * *(fff + i - 1);
	S13 = ax41 * *(fff + i + 1) + ax42 * *(fff + i) + ax43 * *(fff + i - 1) + ax44 * *(fff + i - 2);
	// ! 2 阶导数
	S20 = -*(fff + i + 4) + b12 * *(fff + i + 3) + b13 * *(fff + i + 2) + b14 * *(fff + i + 1);
	S21 = *(fff + i + 2) + b22 * *(fff + i + 1) + *(fff + i);
	S22 = *(fff + i + 1) + b22 * *(fff + i) + *(fff + i - 1);
	S23 = b41 * *(fff + i + 1) + b42 * *(fff + i) + b43 * *(fff + i - 1) - *(fff + i - 2);
	// ! 3 阶导数
	S30 = -*(fff + i + 4) + c12 * (*(fff + i + 3) - *(fff + i + 2)) + *(fff + i + 1);
	S31 = -*(fff + i + 3) + c12 * (*(fff + i + 2) - *(fff + i + 1)) + *(fff + i);
	S32 = -*(fff + i + 2) + c12 * (*(fff + i + 1) - *(fff + i)) + *(fff + i - 1);
	S33 = -*(fff + i + 1) + c12 * (*(fff + i) - *(fff + i - 1)) + *(fff + i - 2);

	S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
	S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
	S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
	S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
	// tau = abs(S0 - S1 - S2 + S3);
	tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
	// !-------WENO Z----------------------
	ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
	ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
	ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
	ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
	// !-------WENO J-S----------------------
	// ax0 = CC0 / ((ep + S0) * (ep + S0));
	// ax1 = CC1 / ((ep + S1) * (ep + S1));
	// ax2 = CC2 / ((ep + S2) * (ep + S2));
	// ax3 = CC3 / ((ep + S3) * (ep + S3));
	// !-----------------------------------------------
	am = ax0 + ax1 + ax2 + ax3;
	flagh_hy[iCell][jCell][nn] += abs(ax0 / am - CC0) + abs(ax1 / am - CC1) + abs(ax2 / am - CC2) + abs(ax3 / am - CC3);
	// !  4阶差分格式的通量
	q0 = e11 * *(fff + i + 4) + e12 * *(fff + i + 3) + e13 * *(fff + i + 2) + e14 * *(fff + i + 1);
	q1 = e21 * *(fff + i + 3) + e22 * *(fff + i + 2) + e23 * *(fff + i + 1) + e24 * *(fff + i);
	q2 = e31 * *(fff + i + 2) + e32 * *(fff + i + 1) + e33 * *(fff + i) + e34 * *(fff + i - 1);
	q3 = e41 * *(fff + i + 1) + e42 * *(fff + i) + e43 * *(fff + i - 1) + e44 * *(fff + i - 2);
	// !  由4个4阶差分格式组合成1个7阶差分格式
	ss2 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;
	yy = ss1 + ss2;
	//   */

	return yy;
}

mydouble CADSolver::scheme_df_dx_char_weno70(mydouble *ffz, mydouble *fff, long i)
{

	mydouble a[4][4] = {{-3.0 / 12.0, 13.0 / 12.0, -23.0 / 12.0, 25.0 / 12.0},
						{1.0 / 12.0, -5.0 / 12.0, 13.0 / 12.0, 3.0 / 12.0},
						{-1.0 / 12.0, 7.0 / 12.0, 7.0 / 12.0, -1.0 / 12.0},
						{3.0 / 12.0, 13.0 / 12.0, -5.0 / 12.0, 1.0 / 12.0}};
	mydouble Cw[4] = {1.0 / 35.0, 12.0 / 35.0, 18.0 / 35.0, 4.0 / 35.0};

	mydouble q0, q1, q2, q3, Is0, Is1, Is2, Is3, alpha0, alpha1, alpha2, alpha3, varepsilon = 1E-15;
	mydouble omega0, omega1, omega2, omega3, isumAlpha, absTerm0, absTerm1, fhatp, fhatm, yy;

	// Plus part
	// polynomials
	q0 = a[0][0] * ffz[i - 3] + a[0][1] * ffz[i - 2] + a[0][2] * ffz[i - 1] + a[0][3] * ffz[i + 0];
	q1 = a[1][0] * ffz[i - 2] + a[1][1] * ffz[i - 1] + a[1][2] * ffz[i + 0] + a[1][3] * ffz[i + 1];
	q2 = a[2][0] * ffz[i - 1] + a[2][1] * ffz[i + 0] + a[2][2] * ffz[i + 1] + a[2][3] * ffz[i + 2];
	q3 = a[3][0] * ffz[i + 0] + a[3][1] * ffz[i + 1] + a[3][2] * ffz[i + 2] + a[3][3] * ffz[i + 3];

	// smoothness index
	Is0 = ffz[i - 3] * (547.0 * ffz[i - 3] - 3882.0 * ffz[i - 2] + 4642.0 * ffz[i - 1] - 1854.0 * ffz[i + 0]) + ffz[i - 2] * (7043.0 * ffz[i - 2] - 17246.0 * ffz[i - 1] + 7042.0 * ffz[i + 0]) + ffz[i - 1] * (11003.0 * ffz[i - 1] - 9402.0 * ffz[i + 0]) + ffz[i + 0] * (2107.0 * ffz[i + 0]);
	Is1 = ffz[i - 2] * (267.0 * ffz[i - 2] - 1642.0 * ffz[i - 1] + 1602.0 * ffz[i + 0] - 494.0 * ffz[i + 1]) + ffz[i - 1] * (2843.0 * ffz[i - 1] - 5966.0 * ffz[i + 0] + 1922.0 * ffz[i + 1]) + ffz[i + 0] * (3443.0 * ffz[i + 0] - 2522.0 * ffz[i + 1]) + ffz[i + 1] * (547.0 * ffz[i + 1]);
	Is2 = ffz[i - 1] * (547.0 * ffz[i - 1] - 2522.0 * ffz[i + 0] + 1922.0 * ffz[i + 1] - 494.0 * ffz[i + 2]) + ffz[i + 0] * (3443.0 * ffz[i + 0] - 5966.0 * ffz[i + 1] + 1602.0 * ffz[i + 2]) + ffz[i + 1] * (2843.0 * ffz[i + 1] - 1642.0 * ffz[i + 2]) + ffz[i + 2] * (267.0 * ffz[i + 2]);
	Is3 = ffz[i + 0] * (2107.0 * ffz[i + 0] - 9402.0 * ffz[i + 1] + 7042.0 * ffz[i + 2] - 1854.0 * ffz[i + 3]) + ffz[i + 1] * (11003.0 * ffz[i + 1] - 17246.0 * ffz[i + 2] + 4642.0 * ffz[i + 3]) + ffz[i + 2] * (7043.0 * ffz[i + 2] - 3882.0 * ffz[i + 3]) + ffz[i + 3] * (547.0 * ffz[i + 3]);

	// alpha
	absTerm0 = abs(Is0 - Is3);
	absTerm1 = abs(Is0 - Is1 - Is2 + Is3);

	// WENO Z
	// alpha0 = Cw[0] * (1.0 + pow(absTerm0 / (Is0 + varepsilon), 2));
	// alpha1 = Cw[1] * (1.0 + pow(absTerm1 / (Is1 + varepsilon), 2));
	// alpha2 = Cw[2] * (1.0 + pow(absTerm0 / (Is2 + varepsilon), 2));
	// alpha3 = Cw[3] * (1.0 + pow(absTerm1 / (Is3 + varepsilon), 2));
	// WENO jS
	alpha0 = Cw[0] / ((Is0 + varepsilon) * (Is0 + varepsilon));
	alpha1 = Cw[1] / ((Is1 + varepsilon) * (Is1 + varepsilon));
	alpha2 = Cw[2] / ((Is2 + varepsilon) * (Is2 + varepsilon));
	alpha3 = Cw[3] / ((Is3 + varepsilon) * (Is3 + varepsilon));

	// alpha(r) = Cw(r) / (IS(r) + varepsilon) *! < WENO JS

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;
	omega3 = alpha3 * isumAlpha;

	// WENO plus reconstruction
	fhatp = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3;

	// Minus part
	// polynomials
	q0 = a[0][0] * fff[i + 4] + a[0][1] * fff[i + 3] + a[0][2] * fff[i + 2] + a[0][3] * fff[i + 1];
	q1 = a[1][0] * fff[i + 3] + a[1][1] * fff[i + 2] + a[1][2] * fff[i + 1] + a[1][3] * fff[i + 0];
	q2 = a[2][0] * fff[i + 2] + a[2][1] * fff[i + 1] + a[2][2] * fff[i + 0] + a[2][3] * fff[i - 1];
	q3 = a[3][0] * fff[i + 1] + a[3][1] * fff[i + 0] + a[3][2] * fff[i - 1] + a[3][3] * fff[i - 2];

	// smoothness index
	Is0 = fff[i + 4] * (547.0 * fff[i + 4] - 3882.0 * fff[i + 3] + 4642.0 * fff[i + 2] - 1854.0 * fff[i + 1]) + fff[i + 3] * (7043.0 * fff[i + 3] - 17246.0 * fff[i + 2] + 7042.0 * fff[i + 1]) + fff[i + 2] * (11003.0 * fff[i + 2] - 9402.0 * fff[i + 1]) + fff[i + 1] * (2107.0 * fff[i + 1]);
	Is1 = fff[i + 3] * (267.0 * fff[i + 3] - 1642.0 * fff[i + 2] + 1602.0 * fff[i + 1] - 494.0 * fff[i + 0]) + fff[i + 2] * (2843.0 * fff[i + 2] - 5966.0 * fff[i + 1] + 1922.0 * fff[i + 0]) + fff[i + 1] * (3443.0 * fff[i + 1] - 2522.0 * fff[i + 0]) + fff[i + 0] * (547.0 * fff[i + 0]);
	Is2 = fff[i + 2] * (547.0 * fff[i + 2] - 2522.0 * fff[i + 1] + 1922.0 * fff[i + 0] - 494.0 * fff[i - 1]) + fff[i + 1] * (3443.0 * fff[i + 1] - 5966.0 * fff[i + 0] + 1602.0 * fff[i - 1]) + fff[i + 0] * (2843.0 * fff[i + 0] - 1642.0 * fff[i - 1]) + fff[i - 1] * (267.0 * fff[i - 1]);
	Is3 = fff[i + 1] * (2107.0 * fff[i + 1] - 9402.0 * fff[i + 0] + 7042.0 * fff[i - 1] - 1854.0 * fff[i - 2]) + fff[i + 0] * (11003.0 * fff[i + 0] - 17246.0 * fff[i - 1] + 4642.0 * fff[i - 2]) + fff[i - 1] * (7043.0 * fff[i - 1] - 3882.0 * fff[i - 2]) + fff[i - 2] * (547.0 * fff[i - 2]);

	// alpha
	absTerm0 = abs(Is0 - Is3);
	absTerm1 = abs(Is0 - Is1 - Is2 + Is3);

	// WENO Z
	// alpha0 = Cw[0] * (1.0 + pow(absTerm0 / (Is0 + varepsilon), 2));
	// alpha1 = Cw[1] * (1.0 + pow(absTerm1 / (Is1 + varepsilon), 2));
	// alpha2 = Cw[2] * (1.0 + pow(absTerm0 / (Is2 + varepsilon), 2));
	// alpha3 = Cw[3] * (1.0 + pow(absTerm1 / (Is3 + varepsilon), 2));

	// WENO jS
	alpha0 = Cw[0] / ((Is0 + varepsilon) * (Is0 + varepsilon));
	alpha1 = Cw[1] / ((Is1 + varepsilon) * (Is1 + varepsilon));
	alpha2 = Cw[2] / ((Is2 + varepsilon) * (Is2 + varepsilon));
	alpha3 = Cw[3] / ((Is3 + varepsilon) * (Is3 + varepsilon));

	// omega
	isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3);
	omega0 = alpha0 * isumAlpha;
	omega1 = alpha1 * isumAlpha;
	omega2 = alpha2 * isumAlpha;
	omega3 = alpha3 * isumAlpha;

	// WENO minus reconstruction
	fhatm = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3;

	// Combine plus and minus parts
	yy = (fhatp + fhatm);

	return yy;
}

mydouble CADSolver::scheme_df_dx_char_weno7_z(mydouble *ffz, long i)
{

	mydouble C0 = 1.0 / 35.0, C1 = 12.0 / 35.0, C2 = 18.0 / 35.0, C3 = 4.0 / 35.0;
	mydouble a11 = -2.0 / 6.0, a12 = 9.0 / 6.0, a13 = -18.0 / 6.0, a14 = 11.0 / 6.0;
	mydouble a21 = 1.0 / 6.0, a23 = 3.0 / 6.0, a24 = 2.0 / 6.0;
	mydouble a31 = -2.0 / 6.0, a32 = -3.0 / 6.0, a34 = -1.0 / 6.0;
	mydouble a41 = -11.0 / 6.0, a42 = 18.0 / 6.0, a43 = -9.0 / 6.0, a44 = 2.0 / 6.0;
	mydouble b12 = 4.0, b13 = -5.0, b14 = 2.0, b22 = -2.0;
	mydouble b41 = 2.0, b42 = -5.0, b43 = 4.0, c12 = 3.0;
	mydouble d12 = 13.0 / 12.0, d13 = 1043.0 / 960.0, d14 = 1.0 / 12.0;
	mydouble e11 = -3.0 / 12.0, e12 = 13.0 / 12.0, e13 = -23.0 / 12.0, e14 = 25.0 / 12.0;
	mydouble e21 = 1.0 / 12.0, e22 = -5.0 / 12.0, e23 = 13.0 / 12.0, e24 = 3.0 / 12.0;
	mydouble e31 = -1.0 / 12.0, e32 = 7.0 / 12.0, e33 = 7.0 / 12.0, e34 = -1.0 / 12.0;
	mydouble e41 = 3.0 / 12.0, e42 = 13.0 / 12.0, e43 = -5.0 / 12.0, e44 = 1.0 / 12.0;
	mydouble epsilon = 1e-15; // WENO-JS

	mydouble S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33;
	mydouble S0, S1, S2, S3;
	mydouble a0, a1, a2, a3, am, q0, q1, q2, q3, hh;

	// 1st derivatives
	S10 = a11 * ffz[i - 3] + a12 * ffz[i - 2] + a13 * ffz[i - 1] + a14 * ffz[i + 0];
	S11 = a21 * ffz[i - 2] - ffz[i - 1] + a23 * ffz[i + 0] + a24 * ffz[i + 1];
	S12 = a31 * ffz[i - 1] + a32 * ffz[i + 0] + ffz[i + 1] + a34 * ffz[i + 2];
	S13 = a41 * ffz[i + 0] + a42 * ffz[i + 1] + a43 * ffz[i + 2] + a44 * ffz[i + 3];

	// 2nd derivatives
	S20 = -ffz[i - 3] + b12 * ffz[i - 2] + b13 * ffz[i - 1] + b14 * ffz[i + 0];
	S21 = ffz[i - 1] + b22 * ffz[i + 0] + ffz[i + 1];
	S22 = ffz[i + 0] + b22 * ffz[i + 1] + ffz[i + 2];
	S23 = b41 * ffz[i + 0] + b42 * ffz[i + 1] + b43 * ffz[i + 2] - ffz[i + 3];

	// 3rd derivatives
	S30 = -ffz[i - 3] + c12 * (ffz[i - 2] - ffz[i - 1]) + ffz[i + 0];
	S31 = -ffz[i - 2] + c12 * (ffz[i - 1] - ffz[i + 0]) + ffz[i + 1];
	S32 = -ffz[i - 1] + c12 * (ffz[i + 0] - ffz[i + 1]) + ffz[i + 2];
	S33 = -ffz[i + 0] + c12 * (ffz[i + 1] - ffz[i + 2]) + ffz[i + 3];

	S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
	S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
	S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
	S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

	// WENO-JS
	a0 = C0 / ((epsilon + S0) * (epsilon + S0));
	a1 = C1 / ((epsilon + S1) * (epsilon + S1));
	a2 = C2 / ((epsilon + S2) * (epsilon + S2));
	a3 = C3 / ((epsilon + S3) * (epsilon + S3));

	am = a0 + a1 + a2 + a3;

	// Flux calculation using 4th order differences
	q0 = e11 * ffz[i - 3] + e12 * ffz[i - 2] + e13 * ffz[i - 1] + e14 * ffz[i + 0];
	q1 = e21 * ffz[i - 2] + e22 * ffz[i - 1] + e23 * ffz[i + 0] + e24 * ffz[i + 1];
	q2 = e31 * ffz[i - 1] + e32 * ffz[i + 0] + e33 * ffz[i + 1] + e34 * ffz[i + 2];
	q3 = e41 * ffz[i + 0] + e42 * ffz[i + 1] + e43 * ffz[i + 2] + e44 * ffz[i + 3];

	// Combine four 4th order differences into one 7th order difference
	hh = (a0 * q0 + a1 * q1 + a2 * q2 + a3 * q3) / am;

	return hh;
}

mydouble CADSolver::scheme_df_dx_char_weno7_f(mydouble *fff, long i)
{

	const double C0 = 1.0 / 35.0, C1 = 12.0 / 35.0, C2 = 18.0 / 35.0, C3 = 4.0 / 35.0;
	const double a11 = -2.0 / 6.0, a12 = 9.0 / 6.0, a13 = -18.0 / 6.0, a14 = 11.0 / 6.0;
	const double a21 = 1.0 / 6.0, a23 = 3.0 / 6.0, a24 = 2.0 / 6.0;
	const double a31 = -2.0 / 6.0, a32 = -3.0 / 6.0, a34 = -1.0 / 6.0;
	const double a41 = -11.0 / 6.0, a42 = 18.0 / 6.0, a43 = -9.0 / 6.0, a44 = 2.0 / 6.0;
	const double b12 = 4.0, b13 = -5.0, b14 = 2.0, b22 = -2.0;
	const double b41 = 2.0, b42 = -5.0, b43 = 4.0, c12 = 3.0;
	const double d12 = 13.0 / 12.0, d13 = 1043.0 / 960.0, d14 = 1.0 / 12.0;
	const double e11 = -3.0 / 12.0, e12 = 13.0 / 12.0, e13 = -23.0 / 12.0, e14 = 25.0 / 12.0;
	const double e21 = 1.0 / 12.0, e22 = -5.0 / 12.0, e23 = 13.0 / 12.0, e24 = 3.0 / 12.0;
	const double e31 = -1.0 / 12.0, e32 = 7.0 / 12.0, e33 = 7.0 / 12.0, e34 = -1.0 / 12.0;
	const double e41 = 3.0 / 12.0, e42 = 13.0 / 12.0, e43 = -5.0 / 12.0, e44 = 1.0 / 12.0;
	const double ep = 1e-8;

	double S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33;
	double S0, S1, S2, S3;
	double a0, a1, a2, a3, am, q0, q1, q2, q3, hh;

	// 一阶导数
	S10 = a11 * fff[i + 4] + a12 * fff[i + 3] + a13 * fff[i + 2] + a14 * fff[i + 1];
	S11 = a21 * fff[i + 3] - fff[i + 2] + a23 * fff[i + 1] + a24 * fff[i + 0];
	S12 = a31 * fff[i + 2] + a32 * fff[i + 1] + fff[i + 0] + a34 * fff[i - 1];
	S13 = a41 * fff[i + 1] + a42 * fff[i + 0] + a43 * fff[i - 1] + a44 * fff[i - 2];

	// 二阶导数
	S20 = -fff[i + 4] + b12 * fff[i + 3] + b13 * fff[i + 2] + b14 * fff[i + 1];
	S21 = fff[i + 2] + b22 * fff[i + 1] + fff[i + 0];
	S22 = fff[i + 1] + b22 * fff[i + 0] + fff[i - 1];
	S23 = b41 * fff[i + 1] + b42 * fff[i + 0] + b43 * fff[i - 1] - fff[i - 2];

	// 三阶导数
	S30 = -fff[i + 4] + c12 * (fff[i + 3] - fff[i + 2]) + fff[i + 1];
	S31 = -fff[i + 3] + c12 * (fff[i + 2] - fff[i + 1]) + fff[i + 0];
	S32 = -fff[i + 2] + c12 * (fff[i + 1] - fff[i + 0]) + fff[i - 1];
	S33 = -fff[i + 1] + c12 * (fff[i + 0] - fff[i - 1]) + fff[i - 2];

	// 一阶和二阶导数平方
	S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
	S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
	S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
	S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

	// 权重系数
	a0 = C0 / pow((ep + S0), 2);
	a1 = C1 / pow((ep + S1), 2);
	a2 = C2 / pow((ep + S2), 2);
	a3 = C3 / pow((ep + S3), 2);

	// 总权重
	am = a0 + a1 + a2 + a3;

	// 4阶差分格式的通量
	q0 = e11 * fff[i + 4] + e12 * fff[i + 3] + e13 * fff[i + 2] + e14 * fff[i + 1];
	q1 = e21 * fff[i + 3] + e22 * fff[i + 2] + e23 * fff[i + 1] + e24 * fff[i + 0];
	q2 = e31 * fff[i + 2] + e32 * fff[i + 1] + e33 * fff[i + 0] + e34 * fff[i - 1];
	q3 = e41 * fff[i + 1] + e42 * fff[i + 0] + e43 * fff[i - 1] + e44 * fff[i - 2];

	// 由4个4阶差分格式组合成1个7阶差分格式
	hh = (a0 * q0 + a1 * q1 + a2 * q2 + a3 * q3) / am;

	return hh;
}

//___________________du_no_char__________________
void CADSolver::Splitting(void)
{
	// #pragma omp barrier

	D_max_u = 0., D_max_v = 0.;
#pragma omp parallel for shared(D_max_u, D_max_v)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
#pragma omp critical
			{
				if (D_max_u < fabs(D_duvpc[iCell][jCell][1]) + D_duvpc[iCell][jCell][4])
					D_max_u = fabs(D_duvpc[iCell][jCell][1]) + D_duvpc[iCell][jCell][4];
				if (D_max_v < fabs(D_duvpc[iCell][jCell][2]) + D_duvpc[iCell][jCell][4])
					D_max_v = fabs(D_duvpc[iCell][jCell][2]) + D_duvpc[iCell][jCell][4];
			}
		}
	}
	// #pragma omp barrier
	mydouble pressure, density, velocity_x, velocity_y, cc, lmax;
	long jCell, iCell, nn;
#pragma omp parallel for private(pressure, density, velocity_x, velocity_y, cc, lmax, jCell, iCell, nn)
	for (jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu[5], fxz[5], fxf[5], fyz[5], fyf[5];
		for (iCell = 0; iCell < nCell_2L; iCell++)
		{

			density = D_duvpc[iCell][jCell][0];
			velocity_x = D_duvpc[iCell][jCell][1];
			velocity_y = D_duvpc[iCell][jCell][2];
			pressure = D_duvpc[iCell][jCell][3];
			cc = D_duvpc[iCell][jCell][4];

			for (nn = 0; nn < 4; nn++)
				uuu[nn] = D_un[iCell][jCell][nn];

			if (Local_Lax_Friedrichs == 1)
				lmax = fabs(velocity_x) + cc;
			else
				lmax = 1.0 * D_max_u;

			fxz[0] = 0.5 * (uuu[1] + lmax * uuu[0]); // 正通量
			fxz[1] = 0.5 * (uuu[1] * velocity_x + pressure + lmax * uuu[1]);
			fxz[2] = 0.5 * (uuu[1] * velocity_y + lmax * uuu[2]);
			fxz[3] = 0.5 * ((uuu[3] + pressure) * velocity_x + lmax * uuu[3]);

			fxf[0] = 0.5 * (uuu[1] - lmax * uuu[0]); // 负通量
			fxf[1] = 0.5 * (uuu[1] * velocity_x + pressure - lmax * uuu[1]);
			fxf[2] = 0.5 * (uuu[1] * velocity_y - lmax * uuu[2]);
			fxf[3] = 0.5 * ((uuu[3] + pressure) * velocity_x - lmax * uuu[3]);

			if (Local_Lax_Friedrichs == 1)
				lmax = fabs(velocity_y) + cc;
			else
				lmax = 1.0 * D_max_v;

			fyz[0] = 0.5 * (uuu[2] + lmax * uuu[0]); // 正通量
			fyz[1] = 0.5 * (uuu[2] * velocity_x + lmax * uuu[1]);
			fyz[2] = 0.5 * (uuu[2] * velocity_y + pressure + lmax * uuu[2]);
			fyz[3] = 0.5 * ((uuu[3] + pressure) * velocity_y + lmax * uuu[3]);

			fyf[0] = 0.5 * (uuu[2] - lmax * uuu[0]); // 负通量
			fyf[1] = 0.5 * (uuu[2] * velocity_x - lmax * uuu[1]);
			fyf[2] = 0.5 * (uuu[2] * velocity_y + pressure - lmax * uuu[2]);
			fyf[3] = 0.5 * ((uuu[3] + pressure) * velocity_y - lmax * uuu[3]);

			for (nn = 0; nn < 4; nn++)
			{

				D_fxz[iCell][jCell][nn] = fxz[nn];
				D_fxf[iCell][jCell][nn] = fxf[nn];
				D_fyz[iCell][jCell][nn] = fyz[nn];
				D_fyf[iCell][jCell][nn] = fyf[nn];
			}
		}
	}
	// #pragma omp barrier
}

void CADSolver::df_dx_dy(void) //   OpenMP 并行
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb - 4; jCell < nCellYo_L + 5; jCell++)
	{
		for (long iCell = nb - 4; iCell < nCell_L + 5; iCell++)
		{

			mydouble fz[15], ff[15];
			mydouble ss1 = 0., ss2 = 0., sss = 0., sa;
			long i = 5, j = 0;

			for (long nn = 0; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < 6; ww++)
				{
					fz[j] = D_fxz[iCell + ww][jCell][nn];
					ff[j] = D_fxf[iCell + ww][jCell][nn];

					j = j + 1;
				}

				ss1 = scheme_df_dx(fz, ff, SpaceScheme, i) / dx;
				j = 0;
				for (long ww = -i; ww < 6; ww++)
				{
					fz[j] = D_fyz[iCell][jCell + ww][nn];
					ff[j] = D_fyf[iCell][jCell + ww][nn];
					j = j + 1;
				}
				ss2 = scheme_df_dx(fz, ff, SpaceScheme, i) / dy;

				sss = ss1 + ss2;

				D_du[iCell][jCell][nn] = sss;
			}
		}
	}
}

mydouble CADSolver::scheme_df_dx(mydouble *fz, mydouble *ff, long scheme_name, long i)
{
	mydouble sz1, sf1, ss1 = 0.;
	// long ii = i;

	switch (scheme_name)
	{
	case 1:
		sz1 = (fz[i] - fz[i - 1]); // 正通量
		sf1 = (ff[i + 1] - ff[i]); // 负通量
		ss1 = sz1 + sf1;
		break;
	case 5:
		// 5阶 迎风
		sz1 = (-*(fz + i - 3) * 2. + *(fz + i - 2) * 15. - *(fz + i - 1) * 60. + *(fz + i) * 20. + *(fz + i + 1) * 30. - *(fz + i + 2) * 3.) / 60;
		sf1 = (*(ff + i - 2) * 3. - *(ff + i - 1) * 30. - *(ff + i) * 20. + *(ff + i + 1) * 60. - *(ff + i + 2) * 15. + *(ff + i + 3) * 2.) / 60.;
		ss1 = sz1 + sf1;
		break;

	case 7:
		sz1 = (fz[i - 4] * 3. - fz[i - 3] * 28. + fz[i - 2] * 126. - fz[i - 1] * 420. + fz[i] * 105. +
			   fz[i + 1] * 252. - fz[i + 2] * 42. + fz[i + 3] * 4.) /
			  420.; // 正通量
		sf1 = (-ff[i - 3] * 4. + ff[i - 2] * 42. - ff[i - 1] * 252. - ff[i] * 105. +
			   ff[i + 1] * 420. - ff[i + 2] * 126. + ff[i + 3] * 28. - ff[i + 4] * 3.) /
			  420.; // 负通量
		ss1 = sz1 + sf1;
		break;

	case 6:
		sz1 = (-fz[i - 3] + fz[i - 2] * 9. - fz[i - 1] * 45. +
			   fz[i + 1] * 45. - fz[i + 2] * 9. + fz[i + 3]) /
			  60.; // 正通量
		sf1 = (-ff[i - 3] + ff[i - 2] * 9. - ff[i - 1] * 45. +
			   ff[i + 1] * 45. - ff[i + 2] * 9. + ff[i + 3]) /
			  60.; // 负通量
		ss1 = sz1 + sf1;
		break;
	case 50:
		ss1 = scheme_df_dx_weno5(fz, ff, i);
		break;
	case 53:
		ss1 = scheme_df_dx_TENO5(fz, ff, i);
		break;

	case 70:
		ss1 = scheme_df_dx_weno70(fz, ff, i);
		break;

	case 72:
		ss1 = scheme_df_dx_weno7Z(fz, ff, i);
		break;

	case 73:
		ss1 = scheme_df_dx_TENO7(fz, ff, i);
		break;

	// case 72:

	// 	break;
	case 8:
		sf1 = (-*(ff + i - 3) * 5. + *(ff + i - 2) * 60. - *(ff + i - 1) * 420. - *(ff + i) * 378 + *(ff + i + 1) * 1050. - *(ff + i + 2) * 420. + *(ff + i + 3) * 140. - *(ff + i + 4) * 30. + *(ff + i + 5) * 3.) / 840.;
		sz1 = -(-*(fz + i + 3) * 5. + *(fz + i + 2) * 60. - *(fz + i + 1) * 420. - *(fz + i) * 378 + *(fz + i - 1) * 1050. - *(fz + i - 2) * 420. + *(fz + i - 3) * 140. - *(fz + i - 4) * 30. + *(fz + i - 5) * 3.) / 840.;
		ss1 = sz1 + sf1;
		break;
	default:
		cout << "1up Space_Scheme  nochar is not avaliable for this Solver." << endl;
		break;
	}

	return ss1;
}

mydouble CADSolver::scheme_df_dx_weno70(mydouble *fz, mydouble *ff, long i)
{
	mydouble S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33,
		ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, tau;
	mydouble CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
			 ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
			 ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
			 ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
			 ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
			 b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
			 b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
			 d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
	mydouble e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
			 e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
			 e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
			 e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-6; //   !! WENO-JS
	mydouble yy, sz1, sf1;
	// /*   weno7-Z================================
	// ! 7th order WENO scheme
	// ! 1  阶导数      *(fz + i - 3)
	S10 = ax11 * *(fz + i - 3) + ax12 * *(fz + i - 2) + ax13 * *(fz + i - 1) + ax14 * *(fz + i);
	S11 = ax21 * *(fz + i - 2) - *(fz + i - 1) + ax23 * *(fz + i) + ax24 * *(fz + i + 1);
	S12 = ax31 * *(fz + i - 1) + ax32 * *(fz + i) + *(fz + i + 1) + ax34 * *(fz + i + 2);
	S13 = ax41 * *(fz + i) + ax42 * *(fz + i + 1) + ax43 * *(fz + i + 2) + ax44 * *(fz + i + 3);
	//  ! 2 阶导数
	S20 = -*(fz + i - 3) + b12 * *(fz + i - 2) + b13 * *(fz + i - 1) + b14 * *(fz + i);
	S21 = *(fz + i - 1) + b22 * *(fz + i) + *(fz + i + 1);
	S22 = *(fz + i) + b22 * *(fz + i + 1) + *(fz + i + 2);
	S23 = b41 * *(fz + i) + b42 * *(fz + i + 1) + b43 * *(fz + i + 2) - *(fz + i + 3);
	// ! 3 阶导数
	S30 = -*(fz + i - 3) + c12 * (*(fz + i - 2) - *(fz + i - 1)) + *(fz + i);
	S31 = -*(fz + i - 2) + c12 * (*(fz + i - 1) - *(fz + i)) + *(fz + i + 1);
	S32 = -*(fz + i - 1) + c12 * (*(fz + i) - *(fz + i + 1)) + *(fz + i + 2);
	S33 = -*(fz + i) + c12 * (*(fz + i + 1) - *(fz + i + 2)) + *(fz + i + 3);

	S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
	S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
	S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
	S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
	// tau = abs(S0 - S1 - S2 + S3);
	tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
	// !-------WENO Z----------------------
	// ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
	// ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
	// ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
	// ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
	// !-------WENO J-S----------------------
	ax0 = CC0 / ((ep + S0) * (ep + S0));
	ax1 = CC1 / ((ep + S1) * (ep + S1));
	ax2 = CC2 / ((ep + S2) * (ep + S2));
	ax3 = CC3 / ((ep + S3) * (ep + S3));
	// !-----------------------------------------------
	am = ax0 + ax1 + ax2 + ax3;
	// !  4阶差分格式的通量
	q0 = e11 * *(fz + i - 3) + e12 * *(fz + i - 2) + e13 * *(fz + i - 1) + e14 * *(fz + i);
	q1 = e21 * *(fz + i - 2) + e22 * *(fz + i - 1) + e23 * *(fz + i) + e24 * *(fz + i + 1);
	q2 = e31 * *(fz + i - 1) + e32 * *(fz + i) + e33 * *(fz + i + 1) + e34 * *(fz + i + 2);
	q3 = e41 * *(fz + i) + e42 * *(fz + i + 1) + e43 * *(fz + i + 2) + e44 * *(fz + i + 3);
	// !  由4个4阶差分格式组合成1个7阶差分格式
	// !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
	sz1 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

	// -------------------------------
	// i = ii; ////////////////////////
	// !      7th order WENO scheme
	// ! 1  阶导数   *(ff + i - 2)
	S10 = ax11 * *(ff + i + 4) + ax12 * *(ff + i + 3) + ax13 * *(ff + i + 2) + ax14 * *(ff + i + 1);
	S11 = ax21 * *(ff + i + 3) - *(ff + i + 2) + ax23 * *(ff + i + 1) + ax24 * *(ff + i);
	S12 = ax31 * *(ff + i + 2) + ax32 * *(ff + i + 1) + *(ff + i) + ax34 * *(ff + i - 1);
	S13 = ax41 * *(ff + i + 1) + ax42 * *(ff + i) + ax43 * *(ff + i - 1) + ax44 * *(ff + i - 2);
	// ! 2 阶导数
	S20 = -*(ff + i + 4) + b12 * *(ff + i + 3) + b13 * *(ff + i + 2) + b14 * *(ff + i + 1);
	S21 = *(ff + i + 2) + b22 * *(ff + i + 1) + *(ff + i);
	S22 = *(ff + i + 1) + b22 * *(ff + i) + *(ff + i - 1);
	S23 = b41 * *(ff + i + 1) + b42 * *(ff + i) + b43 * *(ff + i - 1) - *(ff + i - 2);
	// ! 3 阶导数
	S30 = -*(ff + i + 4) + c12 * (*(ff + i + 3) - *(ff + i + 2)) + *(ff + i + 1);
	S31 = -*(ff + i + 3) + c12 * (*(ff + i + 2) - *(ff + i + 1)) + *(ff + i);
	S32 = -*(ff + i + 2) + c12 * (*(ff + i + 1) - *(ff + i)) + *(ff + i - 1);
	S33 = -*(ff + i + 1) + c12 * (*(ff + i) - *(ff + i - 1)) + *(ff + i - 2);

	S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
	S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
	S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
	S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
	// tau = abs(S0 - S1 - S2 + S3);
	tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
	// !-------WENO Z----------------------
	// ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
	// ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
	// ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
	// ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
	// !-------WENO J-S----------------------
	ax0 = CC0 / ((ep + S0) * (ep + S0));
	ax1 = CC1 / ((ep + S1) * (ep + S1));
	ax2 = CC2 / ((ep + S2) * (ep + S2));
	ax3 = CC3 / ((ep + S3) * (ep + S3));
	// !-----------------------------------------------
	am = ax0 + ax1 + ax2 + ax3;
	// !  4阶差分格式的通量
	q0 = e11 * *(ff + i + 4) + e12 * *(ff + i + 3) + e13 * *(ff + i + 2) + e14 * *(ff + i + 1);
	q1 = e21 * *(ff + i + 3) + e22 * *(ff + i + 2) + e23 * *(ff + i + 1) + e24 * *(ff + i);
	q2 = e31 * *(ff + i + 2) + e32 * *(ff + i + 1) + e33 * *(ff + i) + e34 * *(ff + i - 1);
	q3 = e41 * *(ff + i + 1) + e42 * *(ff + i) + e43 * *(ff + i - 1) + e44 * *(ff + i - 2);
	// !  由4个4阶差分格式组合成1个7阶差分格式
	sf1 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;
	yy = sz1 + sf1;

	return yy;
}

mydouble CADSolver::scheme_df_dx_weno7Z(mydouble *fz, mydouble *ff, long ii)
{
	mydouble S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33, ss1, ss2,
		ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, tau, vz[2];
	mydouble CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
			 ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
			 ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
			 ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
			 ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
			 b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
			 b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
			 d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
	mydouble e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
			 e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
			 e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
			 e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-16; //   !! WENO-JS
	mydouble yy, sz1, sf1;
	long i = ii;

	// /*   weno7-Z================================
	// ! 7th order WENO scheme
	// ! 1  阶导数      *(fz + i - 3)
	for (long j = 0; j < 2; j++)
	{

		S10 = ax11 * *(fz + i - 3) + ax12 * *(fz + i - 2) + ax13 * *(fz + i - 1) + ax14 * *(fz + i);
		S11 = ax21 * *(fz + i - 2) - *(fz + i - 1) + ax23 * *(fz + i) + ax24 * *(fz + i + 1);
		S12 = ax31 * *(fz + i - 1) + ax32 * *(fz + i) + *(fz + i + 1) + ax34 * *(fz + i + 2);
		S13 = ax41 * *(fz + i) + ax42 * *(fz + i + 1) + ax43 * *(fz + i + 2) + ax44 * *(fz + i + 3);
		//  ! 2 阶导数
		S20 = -*(fz + i - 3) + b12 * *(fz + i - 2) + b13 * *(fz + i - 1) + b14 * *(fz + i);
		S21 = *(fz + i - 1) + b22 * *(fz + i) + *(fz + i + 1);
		S22 = *(fz + i) + b22 * *(fz + i + 1) + *(fz + i + 2);
		S23 = b41 * *(fz + i) + b42 * *(fz + i + 1) + b43 * *(fz + i + 2) - *(fz + i + 3);
		// ! 3 阶导数
		S30 = -*(fz + i - 3) + c12 * (*(fz + i - 2) - *(fz + i - 1)) + *(fz + i);
		S31 = -*(fz + i - 2) + c12 * (*(fz + i - 1) - *(fz + i)) + *(fz + i + 1);
		S32 = -*(fz + i - 1) + c12 * (*(fz + i) - *(fz + i + 1)) + *(fz + i + 2);
		S33 = -*(fz + i) + c12 * (*(fz + i + 1) - *(fz + i + 2)) + *(fz + i + 3);

		S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
		S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
		S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
		S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
		// tau = abs(S0 - S1 - S2 + S3);
		tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
		// !-------WENO Z----------------------
		ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
		ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
		ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
		ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
		// !-------WENO J-S----------------------
		// ax0 = CC0 / ((ep + S0) * (ep + S0));
		// ax1 = CC1 / ((ep + S1) * (ep + S1));
		// ax2 = CC2 / ((ep + S2) * (ep + S2));
		// ax3 = CC3 / ((ep + S3) * (ep + S3));
		// !-----------------------------------------------
		am = ax0 + ax1 + ax2 + ax3;
		// !  4阶差分格式的通量
		q0 = e11 * *(fz + i - 3) + e12 * *(fz + i - 2) + e13 * *(fz + i - 1) + e14 * *(fz + i);
		q1 = e21 * *(fz + i - 2) + e22 * *(fz + i - 1) + e23 * *(fz + i) + e24 * *(fz + i + 1);
		q2 = e31 * *(fz + i - 1) + e32 * *(fz + i) + e33 * *(fz + i + 1) + e34 * *(fz + i + 2);
		q3 = e41 * *(fz + i) + e42 * *(fz + i + 1) + e43 * *(fz + i + 2) + e44 * *(fz + i + 3);
		// !  由4个4阶差分格式组合成1个7阶差分格式
		// !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
		vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

		i--;
	}
	ss1 = (vz[0] - vz[1]);
	// -------------------------------
	i = ii; ////////////////////////
			// !      7th order WENO scheme
			// ! 1  阶导数   *(ff + i - 2)
	for (long j = 0; j < 2; j++)
	{
		S10 = ax11 * *(ff + i + 4) + ax12 * *(ff + i + 3) + ax13 * *(ff + i + 2) + ax14 * *(ff + i + 1);
		S11 = ax21 * *(ff + i + 3) - *(ff + i + 2) + ax23 * *(ff + i + 1) + ax24 * *(ff + i);
		S12 = ax31 * *(ff + i + 2) + ax32 * *(ff + i + 1) + *(ff + i) + ax34 * *(ff + i - 1);
		S13 = ax41 * *(ff + i + 1) + ax42 * *(ff + i) + ax43 * *(ff + i - 1) + ax44 * *(ff + i - 2);
		// ! 2 阶导数
		S20 = -*(ff + i + 4) + b12 * *(ff + i + 3) + b13 * *(ff + i + 2) + b14 * *(ff + i + 1);
		S21 = *(ff + i + 2) + b22 * *(ff + i + 1) + *(ff + i);
		S22 = *(ff + i + 1) + b22 * *(ff + i) + *(ff + i - 1);
		S23 = b41 * *(ff + i + 1) + b42 * *(ff + i) + b43 * *(ff + i - 1) - *(ff + i - 2);
		// ! 3 阶导数
		S30 = -*(ff + i + 4) + c12 * (*(ff + i + 3) - *(ff + i + 2)) + *(ff + i + 1);
		S31 = -*(ff + i + 3) + c12 * (*(ff + i + 2) - *(ff + i + 1)) + *(ff + i);
		S32 = -*(ff + i + 2) + c12 * (*(ff + i + 1) - *(ff + i)) + *(ff + i - 1);
		S33 = -*(ff + i + 1) + c12 * (*(ff + i) - *(ff + i - 1)) + *(ff + i - 2);

		S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
		S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
		S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
		S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
		// tau = abs(S0 - S1 - S2 + S3);
		tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
		// !-------WENO Z----------------------
		ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
		ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
		ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
		ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
		// !-------WENO J-S----------------------
		// ax0 = CC0 / ((ep + S0) * (ep + S0));
		// ax1 = CC1 / ((ep + S1) * (ep + S1));
		// ax2 = CC2 / ((ep + S2) * (ep + S2));
		// ax3 = CC3 / ((ep + S3) * (ep + S3));
		// !-----------------------------------------------
		am = ax0 + ax1 + ax2 + ax3;
		// !  4阶差分格式的通量
		q0 = e11 * *(ff + i + 4) + e12 * *(ff + i + 3) + e13 * *(ff + i + 2) + e14 * *(ff + i + 1);
		q1 = e21 * *(ff + i + 3) + e22 * *(ff + i + 2) + e23 * *(ff + i + 1) + e24 * *(ff + i);
		q2 = e31 * *(ff + i + 2) + e32 * *(ff + i + 1) + e33 * *(ff + i) + e34 * *(ff + i - 1);
		q3 = e41 * *(ff + i + 1) + e42 * *(ff + i) + e43 * *(ff + i - 1) + e44 * *(ff + i - 2);
		// !  由4个4阶差分格式组合成1个7阶差分格式
		vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

		i--;
	}
	ss2 = (vz[0] - vz[1]);
	yy = ss1 + ss2;

	return yy;
}

mydouble CADSolver::scheme_df_dx_weno5(mydouble *fz, mydouble *ff, long ii)
{

	mydouble y, ss1, ss2;
	long i = ii;
	// weno_5-js  正通量

	mydouble ep = 1.E-5, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	mydouble sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
	for (long j = 0; j < 2; j++)
	{
		sn0 = 13. / 12. * pow((*(fz + i - 2) - 2. * (*(fz + i - 1)) + *(fz + i)), 2) + 0.25 * pow((*(fz + i - 2) - 4. * (*(fz + i - 1)) + 3. * (*(fz + i))), 2);
		sn1 = 13. / 12. * pow((*(fz + i - 1) - 2. * (*(fz + i)) + *(fz + i + 1)), 2) + 0.25 * pow((*(fz + i - 1) - *(fz + i + 1)), 2);
		sn2 = 13. / 12. * pow((*(fz + i) - 2. * (*(fz + i + 1)) + *(fz + i + 2)), 2) + 0.25 * pow((3. * (*(fz + i)) - 4. * (*(fz + i + 1)) + *(fz + i + 2)), 2);
		az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
		W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
		q03 = 2. / 6. * (*(fz + i - 2)) - 7. / 6. * (*(fz + i - 1)) + 11. / 6. * (*(fz + i));
		q13 = -1. / 6. * (*(fz + i - 1)) + 5. / 6. * (*(fz + i)) + 2. / 6. * (*(fz + i + 1));
		q23 = 2. / 6. * (*(fz + i)) + 5. / 6. * (*(fz + i + 1)) - 1. / 6. * (*(fz + i + 2));
		vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
		i--;
	}
	ss1 = (vz[0] - vz[1]);
	// weno_5-js 负通量
	//  double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	//  double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
	i = ii, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
	for (long j = 0; j < 2; j++)
	{
		sn0 = 13. / 12. * pow((*(ff + i + 1) - 2. * (*(ff + i + 2)) + *(ff + i + 3)), 2) + 0.25 * pow((3. * (*(ff + i + 1)) - 4. * (*(ff + i + 2)) + *(ff + i + 3)), 2);
		sn1 = 13. / 12. * pow((*(ff + i) - 2. * (*(ff + i + 1)) + *(ff + i + 2)), 2) + 0.25 * pow((*(ff + i) - *(ff + i + 2)), 2);
		sn2 = 13. / 12. * pow((*(ff + i - 1) - 2. * (*(ff + i)) + *(ff + i + 1)), 2) + 0.25 * pow((*(ff + i - 1) - 4. * (*(ff + i)) + 3. * (*(ff + i + 1))), 2);
		az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
		W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
		q03 = 11. / 6. * (*(ff + i + 1)) - 7. / 6. * (*(ff + i + 2)) + 2. / 6. * (*(ff + i + 3));
		q13 = 2. / 6. * (*(ff + i)) + 5. / 6. * (*(ff + i + 1)) - 1. / 6. * (*(ff + i + 2));
		q23 = -1. / 6. * (*(ff + i - 1)) + 5. / 6. * (*(ff + i)) + 2. / 6. * (*(ff + i + 1));
		vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
		i--;
	}
	ss2 = (vz[0] - vz[1]);
	y = ss1 + ss2;
	return y;
}

mydouble CADSolver::scheme_df_dx_TENO5(mydouble *ffz, mydouble *fff, long ii)
{

	mydouble y, ss1, ss2, vz[2];
	long i = ii;
	mydouble varepsilon = 1.0e-16;
	S_ct = 1.0e-10;

	// mydouble S_ct = 1.0e-7;
	// mydouble S_ct = 1.0e-3; // High_jet case

	mydouble q0, q1, q2, q3, IS0, IS1, IS2, Is3, eps, eta0, eta1, eta2, rm;
	mydouble alpha0, alpha1, alpha2, alpha3;
	mydouble omega0, omega1, omega2, omega3;
	mydouble isumAlpha, absTerm0, absTerm1, d0, d1, d2, fhatp, fhatm, dsum;

	mydouble is1par = 13.0 / 12.0;
	mydouble is2par = 0.25;

	mydouble a[3][3] = {{2.0 / 6.0, -7.0 / 6.0, 11.0 / 6.0},
						{-1.0 / 6.0, 5.0 / 6.0, 2.0 / 6.0},
						{2.0 / 6.0, 5.0 / 6.0, -1.0 / 6.0}};

	mydouble cw[4] = {0.1, 0.6, 0.3};
	// TENO5 正通量
	for (long j = 0; j < 2; j++)
	{
		// polinomia
		q0 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
		q1 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
		q2 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];

		// smoothness index
		IS0 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
		IS1 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
		IS2 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);

		// alpha
		absTerm0 = abs(IS0 - IS2);
		alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
		alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
		alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);

		// omega
		isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2);
		omega0 = alpha0 * isumAlpha;
		omega1 = alpha1 * isumAlpha;
		omega2 = alpha2 * isumAlpha;

		// delta
		d0 = omega0 < S_ct ? 0. : 1.;
		d1 = omega1 < S_ct ? 0. : 1.;
		d2 = omega2 < S_ct ? 0. : 1.;

		dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2]);

		omega0 = d0 * cw[0] * dsum;
		omega1 = d1 * cw[1] * dsum;
		omega2 = d2 * cw[2] * dsum;

		// WENO plus reconstruction
		vz[j] = omega0 * q0 + omega1 * q1 + omega2 * q2;

		i--;
	}
	ss1 = (vz[0] - vz[1]);

	// TENO5 负通量
	i = ii;
	for (long j = 0; j < 2; j++)
	{
		// --- WENO minus part --- !

		// polinomia
		q0 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
		q1 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
		q2 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];

		// smoothness index
		IS0 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
		IS1 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
		IS2 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);

		// alpha
		absTerm0 = abs(IS0 - IS2);
		alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
		alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
		alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);

		// omega
		isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2);
		omega0 = alpha0 * isumAlpha;
		omega1 = alpha1 * isumAlpha;
		omega2 = alpha2 * isumAlpha;

		// delta

		d0 = omega0 < S_ct ? 0. : 1.;
		d1 = omega1 < S_ct ? 0. : 1.;
		d2 = omega2 < S_ct ? 0. : 1.;

		dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2]);

		omega0 = d0 * cw[0] * dsum;
		omega1 = d1 * cw[1] * dsum;
		omega2 = d2 * cw[2] * dsum;

		// WENO minus reconstruction
		vz[j] = omega0 * q0 + omega1 * q1 + omega2 * q2;

		i--;
	}
	ss2 = (vz[0] - vz[1]);
	y = ss1 + ss2;
	return y;
}

mydouble CADSolver::scheme_df_dx_TENO7(mydouble *ffz, mydouble *fff, long ii)
{

	mydouble y, ss1, ss2, vz[2];
	long i = ii;
	mydouble varepsilon = 1.0e-16;

	// mydouble S_ct = 1.0e-7;

	mydouble q0, q1, q2, q3, q4, IS0, IS1, IS2, IS3, IS4, eps, eta0, eta1, eta2, rm;
	mydouble alpha0, alpha1, alpha2, alpha3, alpha4;
	mydouble omega0, omega1, omega2, omega3, omega4;
	mydouble isumAlpha, absTerm0, absTerm1, d0, d1, d2, d3, d4, fhatp, fhatm, dsum, tao5;

	mydouble is1par = 13.0 / 12.0;
	mydouble is2par = 0.25;

	mydouble a[5][4] = {{2.0 / 6.0, -7.0 / 6.0, 11.0 / 6.0, 0},
						{-1.0 / 6.0, 5.0 / 6.0, 2.0 / 6.0, 0},
						{2.0 / 6.0, 5.0 / 6.0, -1.0 / 6.0, 0},
						{3.0 / 12.0, 13.0 / 12.0, -5.0 / 12.0, 1.0 / 12.0},
						{-3.0 / 12.0, 13.0 / 12.0, -23.0 / 12.0, 25.0 / 12.0}};

	// mydouble cw[4] = {0.1, 0.6, 0.3};
	// mydouble cw[4] = {1. / 20., 9. / 20., 6. / 20., 4. / 20.};
	// mydouble cw[5] = {3. / 35., 18. / 35., 9. / 35., 4. / 35., 1. / 35.};
	mydouble cw[5] = {18. / 35., 9. / 35., 3. / 35., 4. / 35., 1. / 35.};

	// TENO5 正通量
	for (long j = 0; j < 2; j++)
	{
		// polinomia
		// q0 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
		// q1 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
		// q2 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];

		q2 = a[0][0] * ffz[i - 2] + a[0][1] * ffz[i - 1] + a[0][2] * ffz[i + 0];
		q0 = a[1][0] * ffz[i - 1] + a[1][1] * ffz[i + 0] + a[1][2] * ffz[i + 1];
		q1 = a[2][0] * ffz[i + 0] + a[2][1] * ffz[i + 1] + a[2][2] * ffz[i + 2];

		q3 = a[3][0] * ffz[i + 0] + a[3][1] * ffz[i + 1] + a[3][2] * ffz[i + 2] + a[3][3] * ffz[i + 3];
		q4 = a[4][0] * ffz[i - 3] + a[4][1] * ffz[i - 2] + a[4][2] * ffz[i - 1] + a[4][3] * ffz[i];

		// smoothness index
		// IS0 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
		// IS1 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
		// IS2 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);

		IS2 = is1par * pow(ffz[i - 2] - 2.0 * ffz[i - 1] + ffz[i + 0], 2) + is2par * pow(ffz[i - 2] - 4.0 * ffz[i - 1] + 3.0 * ffz[i + 0], 2);
		IS0 = is1par * pow(ffz[i - 1] - 2.0 * ffz[i + 0] + ffz[i + 1], 2) + is2par * pow(ffz[i - 1] - ffz[i + 1], 2);
		IS1 = is1par * pow(ffz[i + 0] - 2.0 * ffz[i + 1] + ffz[i + 2], 2) + is2par * pow(3.0 * ffz[i + 0] - 4.0 * ffz[i + 1] + ffz[i + 2], 2);

		IS3 = 1. / 36. * pow(-11. * ffz[i + 0] + 18.0 * ffz[i + 1] - 9. * ffz[i + 2] + 2. * ffz[i + 3], 2) +
			  13. / 12. * pow(2. * ffz[i + 0] - 5.0 * ffz[i + 1] + 4. * ffz[i + 2] - ffz[i + 3], 2) +
			  781. / 720. * pow(-ffz[i + 0] + 3.0 * ffz[i + 1] - 3. * ffz[i + 2] + ffz[i + 3], 2);

		IS4 = 1. / 36. * pow(-2. * ffz[i - 3] + 9.0 * ffz[i - 2] - 18. * ffz[i - 1] + 11. * ffz[i + 0], 2) +
			  13. / 12. * pow(-ffz[i - 3] + 4.0 * ffz[i - 2] - 5. * ffz[i - 1] + 2. * ffz[i + 0], 2) +
			  781. / 720. * pow(-ffz[i - 3] + 3.0 * ffz[i - 2] - 3. * ffz[i - 1] + ffz[i + 0], 2);
		// alpha

		tao5 = (32154783380.0 * pow(ffz[i + 0], 2) + 18133963560.0 * pow(ffz[i + 1], 2) + 2927992563.0 * pow(ffz[i + 2], 2) -
				969999969.0 * ffz[i + 2] * ffz[i + 3] + 84070496.0 * pow(ffz[i + 3], 2) - 12546315963.0 * ffz[i + 2] * ffz[i - 1] +
				1902531828.0 * ffz[i + 3] * ffz[i - 1] + 18133963560.0 * pow(ffz[i - 1], 2) + 4550242446.0 * ffz[i + 2] * ffz[i - 2] -
				676871859.0 * ffz[i + 3] * ffz[i - 2] - 14296379553.0 * ffz[i - 1] * ffz[i - 2] + 2927992563.0 * pow(ffz[i - 2], 2) -
				676871859.0 * ffz[i + 2] * ffz[i - 3] + 99022657.0 * ffz[i + 3] * ffz[i - 3] + 2283428883.0 * ffz[i - 1] * ffz[i - 3] -
				969999969.0 * ffz[i - 2] * ffz[i - 3] + 84070496.0 * pow(ffz[i - 3], 2) +
				3.0 * ffz[i + 1] * (-4765459851.0 * ffz[i + 2] + 761142961.0 * ffz[i + 3] + 11273559435.0 * ffz[i - 1] - 4182105321.0 * ffz[i - 2] + 634177276.0 * ffz[i - 3]) -
				4.0 * ffz[i + 0] * (11857967655.0 * ffz[i + 1] - 4520834943.0 * ffz[i + 2] + 701563133.0 * ffz[i + 3] + 11857967655.0 * ffz[i - 1] - 4520834943.0 * ffz[i - 2] + 701563133.0 * ffz[i - 3])) /
			   59875200.0;

		// alpha
		// absTerm0 = abs(IS0 - IS2);
		// = β3 − 1/6(β0 + β2 + 4β1),
		// absTerm0 = abs(IS3 - (IS1 + IS0 + 4. * IS2) / 6.);
		absTerm0 = abs(tao5 - (IS1 + IS2 + 4. * IS0) / 6.);

		alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
		alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
		alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);
		alpha3 = 1.0 + pow(absTerm0 / (IS3 + varepsilon), 6);
		alpha4 = 1.0 + pow(absTerm0 / (IS4 + varepsilon), 6);

		// omega
		isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3 + alpha4);
		omega0 = alpha0 * isumAlpha;
		omega1 = alpha1 * isumAlpha;
		omega2 = alpha2 * isumAlpha;
		omega3 = alpha3 * isumAlpha;
		omega4 = alpha4 * isumAlpha;

		// delta
		d0 = omega0 < S_ct ? 0. : 1.;
		d1 = omega1 < S_ct ? 0. : 1.;
		d2 = omega2 < S_ct ? 0. : 1.;
		d3 = omega3 < S_ct ? 0. : 1.;
		d4 = omega4 < S_ct ? 0. : 1.;

		dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3] + d4 * cw[4]);

		omega0 = d0 * cw[0] * dsum;
		omega1 = d1 * cw[1] * dsum;
		omega2 = d2 * cw[2] * dsum;
		omega3 = d3 * cw[3] * dsum;
		omega4 = d4 * cw[4] * dsum;

		// WENO plus reconstruction
		vz[j] = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3 + omega4 * q4;

		i--;
	}
	ss1 = (vz[0] - vz[1]);

	// TENO5 负通量
	i = ii;
	for (long j = 0; j < 2; j++)
	{
		// --- WENO minus part --- !

		// polinomia
		// q0 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
		// q1 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
		// q2 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];

		q2 = a[0][0] * fff[i + 3] + a[0][1] * fff[i + 2] + a[0][2] * fff[i + 1];
		q0 = a[1][0] * fff[i + 2] + a[1][1] * fff[i + 1] + a[1][2] * fff[i + 0];
		q1 = a[2][0] * fff[i + 1] + a[2][1] * fff[i + 0] + a[2][2] * fff[i - 1];
		q3 = a[3][0] * fff[i + 1] + a[3][1] * fff[i + 0] + a[3][2] * fff[i - 1] + a[3][3] * fff[i - 2];
		q4 = a[4][0] * fff[i + 4] + a[4][1] * fff[i + 3] + a[4][2] * fff[i + 2] + a[4][3] * fff[i + 1];

		// smoothness index
		// IS0 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
		// IS1 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
		// IS2 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);

		IS2 = is1par * pow(fff[i + 3] - 2.0 * fff[i + 2] + fff[i + 1], 2) + is2par * pow(fff[i + 3] - 4.0 * fff[i + 2] + 3.0 * fff[i + 1], 2);
		IS0 = is1par * pow(fff[i + 2] - 2.0 * fff[i + 1] + fff[i + 0], 2) + is2par * pow(fff[i + 2] - fff[i + 0], 2);
		IS1 = is1par * pow(fff[i + 1] - 2.0 * fff[i + 0] + fff[i - 1], 2) + is2par * pow(3.0 * fff[i + 1] - 4.0 * fff[i + 0] + fff[i - 1], 2);

		IS3 = 1. / 36. * pow(-11. * fff[i + 1] + 18.0 * fff[i + 0] - 9. * fff[i - 1] + 2. * fff[i - 2], 2) +
			  13. / 12. * pow(2. * fff[i + 1] - 5.0 * fff[i + 0] + 4. * fff[i - 1] - fff[i - 2], 2) +
			  781. / 720. * pow(-fff[i + 1] + 3.0 * fff[i + 0] - 3. * fff[i - 1] + fff[i - 2], 2);
		IS4 = 1. / 36. * pow(-2. * fff[i + 4] + 9.0 * fff[i + 3] - 18. * fff[i + 2] + 11. * fff[i + 1], 2) +
			  13. / 12. * pow(-fff[i + 4] + 4.0 * fff[i + 3] - 5. * fff[i + 2] + 2. * fff[i + 1], 2) +
			  781. / 720. * pow(-fff[i + 4] + 3.0 * fff[i + 3] - 3. * fff[i + 2] + fff[i + 1], 2);
		// alpha
		tao5 = (32154783380.0 * pow(fff[i + 1], 2) + 18133963560.0 * pow(fff[i + 0], 2) + 2927992563.0 * pow(fff[i - 1], 2) -
				969999969.0 * fff[i - 1] * fff[i - 2] + 84070496.0 * pow(fff[i - 2], 2) - 12546315963.0 * fff[i - 1] * fff[i + 2] +
				1902531828.0 * fff[i - 2] * fff[i + 2] + 18133963560.0 * pow(fff[i + 2], 2) + 4550242446.0 * fff[i - 1] * fff[i + 3] -
				676871859.0 * fff[i - 2] * fff[i + 3] - 14296379553.0 * fff[i + 2] * fff[i + 3] + 2927992563.0 * pow(fff[i + 3], 2) -
				676871859.0 * fff[i - 1] * fff[i + 4] + 99022657.0 * fff[i - 2] * fff[i + 4] + 2283428883.0 * fff[i + 2] * fff[i + 4] -
				969999969.0 * fff[i + 3] * fff[i + 4] + 84070496.0 * pow(fff[i + 4], 2) +
				3.0 * fff[i + 0] * (-4765459851.0 * fff[i - 1] + 761142961.0 * fff[i - 2] + 11273559435.0 * fff[i + 2] - 4182105321.0 * fff[i + 3] + 634177276.0 * fff[i + 4]) -
				4.0 * fff[i + 1] * (11857967655.0 * fff[i + 0] - 4520834943.0 * fff[i - 1] + 701563133.0 * fff[i - 2] + 11857967655.0 * fff[i + 2] - 4520834943.0 * fff[i + 3] + 701563133.0 * fff[i + 4])) /
			   59875200.0;

		// absTerm0 = abs(IS3 - (IS1 + IS0 + 4. * IS2) / 6.);
		// absTerm1 = abs(IS4 - (IS1 + IS0 + 4. * IS2) / 6.);

		absTerm0 = abs(tao5 - (IS1 + IS2 + 4. * IS0) / 6.);
		// absTerm0 = abs(IS0 - IS2);
		alpha0 = 1.0 + pow(absTerm0 / (IS0 + varepsilon), 6);
		alpha1 = 1.0 + pow(absTerm0 / (IS1 + varepsilon), 6);
		alpha2 = 1.0 + pow(absTerm0 / (IS2 + varepsilon), 6);
		alpha3 = 1.0 + pow(absTerm0 / (IS3 + varepsilon), 6);
		alpha4 = 1.0 + pow(absTerm0 / (IS4 + varepsilon), 6);

		// omega
		isumAlpha = 1.0 / (alpha0 + alpha1 + alpha2 + alpha3 + alpha4);
		omega0 = alpha0 * isumAlpha;
		omega1 = alpha1 * isumAlpha;
		omega2 = alpha2 * isumAlpha;
		omega3 = alpha3 * isumAlpha;
		omega4 = alpha4 * isumAlpha;

		// delta

		d0 = omega0 < S_ct ? 0. : 1.;
		d1 = omega1 < S_ct ? 0. : 1.;
		d2 = omega2 < S_ct ? 0. : 1.;
		d3 = omega3 < S_ct ? 0. : 1.;
		d4 = omega4 < S_ct ? 0. : 1.;

		dsum = 1.0 / (d0 * cw[0] + d1 * cw[1] + d2 * cw[2] + d3 * cw[3] + d4 * cw[4]);

		omega0 = d0 * cw[0] * dsum;
		omega1 = d1 * cw[1] * dsum;
		omega2 = d2 * cw[2] * dsum;
		omega3 = d3 * cw[3] * dsum;
		omega4 = d4 * cw[4] * dsum;

		// WENO minus reconstruction
		vz[j] = omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3 + omega4 * q4;

		i--;
	}
	ss2 = (vz[0] - vz[1]);
	y = ss1 + ss2;
	return y;
}

//___________________Hy__________________

void CADSolver::shock_sensor_Hy(void)
{

	mydouble max_hy, f_hy;
	long nf = 3;

	max_hy = 0.;

#pragma omp parallel for shared(max_hy)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
#pragma omp critical
			{
				flagh_hy[iCell][jCell][4] = 0.;
				if (max_hy < flagh_hy[iCell][jCell][nf])
					max_hy = flagh_hy[iCell][jCell][nf];
			}
		}
	}
	f_hy = 0.7 * max_hy;

#pragma omp parallel for
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			D_out[iCell][jCell][0] = flagh_hy[iCell][jCell][3];
			D_out[iCell][jCell][1] = 0.0; // if (flagh_hy[iCell][jCell][nf] > f_hy)
			// if (flagh_hy[iCell][jCell][nf] > f_hy || iCell < nb + 5 || iCell > nCell_L - 6)
			if (flagh_hy[iCell][jCell][nf] > f_hy)
			{
				D_out[iCell][jCell][1] = 1.0;
				for (int k1 = -1; k1 < 2; k1++)
					for (int k2 = -1; k2 < 2; k2++)
						flagh_hy[iCell + k1][jCell + k2][4] = 1.0;
			}
			if (D_duvpc[iCell][jCell][0] < 1E-5 || D_duvpc[iCell][jCell][3] < 1E-5)
			{
				D_out[iCell][jCell][1] = 1.0;
				for (int k1 = -3; k1 < 4; k1++)
					for (int k2 = -3; k2 < 4; k2++)
						flagh_hy[iCell + k1][jCell + k2][4] = 1.0;
			}
		}
	}
}

void CADSolver::Write_hy(long m)
{

	string Title = "hy.plt";
	long nc = 0;
	ofstream ofsy;
	if (m < 2)
	{
		ofsy.open(Title, ios::out);
		ofsy << " VARIABLES=x,t " << endl;
		for (long i = nb + nc; i < nCell_L - nc; i++) // 间断  int i = nb; i < n1 - nb; i++
		{
			if (flagh_hy[i][nb][4] > 0.5)
				ofsy.precision(10), ofsy << D_Coord[i][nb][0] << "  " << TimeNow << endl;
		}
		ofsy.close();
	}
	else
	{
		ofsy.open(Title, ios::app);
		for (long i = nb + nc; i < nCell_L - nc; i++) // 间断  int i = nb; i < n1 - nb; i++
		{
			if (flagh_hy[i][nb][4] > 0.5)
				ofsy.precision(10), ofsy << D_Coord[i][nb][0] << "  " << TimeNow << endl;
		}
		ofsy.close();
	}
}

//___________________ Viscous_Splitting __________________
void CADSolver::Viscous_mu(void)
{
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble Tn;
		for (long iCell = 0; iCell < nCell_2L; iCell++)
		{

			Tn = D_duvpc[iCell][jCell][5];
			if (Case_name == 24 || Case_name == 26 || Case_name == 27 || Case_name == 83) //
				D_mu[iCell][jCell] = 1. / Re;
			else
				D_mu[iCell][jCell] = pow(Tn / Ts, 1.5) * (Ts + Tc) / (Tn + Tc) / Re;
			D_ka[iCell][jCell] = C_p * D_mu[iCell][jCell] / Pr;

			// if (Case_name == 82)
			// {
			// 	D_mu[iCell][jCell] = 0.00001;
			// 	D_ka[iCell][jCell] = 0.00001;
			// }
		}
	}
}

void CADSolver::Viscous_Splitting(void)
{
	// #pragma omp barrier
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];
	Viscous_mu();

#pragma omp parallel for
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble ux[8], uy[8], vx[8], vy[8], Tx[8], Ty[8];
		mydouble sux, suy, svx, svy, sTx, sTy, tau11, tau12, tau22, mu, ka;
		long i = Viscous_cd / 2, j = 0;
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			j = 0;
			for (long ww = -i; ww < i + 1; ww++)
			{
				ux[j] = D_duvpc[iCell + ww][jCell][1];
				vx[j] = D_duvpc[iCell + ww][jCell][2];
				Tx[j] = D_duvpc[iCell + ww][jCell][5];

				uy[j] = D_duvpc[iCell][jCell + ww][1];
				vy[j] = D_duvpc[iCell][jCell + ww][2];
				Ty[j] = D_duvpc[iCell][jCell + ww][5];
				j = j + 1;
			}

			sux = scheme_ALW_df_dx(ux, Viscous_cd, i) / dx;
			svx = scheme_ALW_df_dx(vx, Viscous_cd, i) / dx;
			sTx = scheme_ALW_df_dx(Tx, Viscous_cd, i) / dx;

			suy = scheme_ALW_df_dx(uy, Viscous_cd, i) / dy;
			svy = scheme_ALW_df_dx(vy, Viscous_cd, i) / dy;
			sTy = scheme_ALW_df_dx(Ty, Viscous_cd, i) / dy;

			mu = D_mu[iCell][jCell];
			ka = D_ka[iCell][jCell];

			tau12 = mu * (suy + svx);
			tau11 = 2. * mu * sux - 2. / 3. * mu * (sux + svy);
			tau22 = 2. * mu * svy - 2. / 3. * mu * (sux + svy);

			D_gx[iCell][jCell][1] = tau11;
			D_gx[iCell][jCell][2] = tau12;
			D_gx[iCell][jCell][3] = ka * sTx + tau11 * ux[i] + tau12 * vx[i];

			D_gy[iCell][jCell][1] = tau12;
			D_gy[iCell][jCell][2] = tau22;
			D_gy[iCell][jCell][3] = ka * sTy + tau12 * ux[i] + tau22 * vx[i];
		}
	}
	// #pragma omp barrier
	Viscous_Schemes();
}

void CADSolver::Viscous_Schemes(void)
{

#pragma omp parallel for
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble fx[8], fy[8];
		mydouble ss1 = 0., ss2 = 0., sss;
		long i = Viscous_cd / 2, j = 0;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (long nn = 1; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < i + 1; ww++)
				{
					fx[j] = D_gx[iCell + ww][jCell][nn];
					fy[j] = D_gy[iCell][jCell + ww][nn];
					j = j + 1;
				}

				ss1 = scheme_ALW_df_dx(fx, Viscous_cd, i) / dx;
				ss2 = scheme_ALW_df_dx(fy, Viscous_cd, i) / dy;

				D_du[iCell][jCell][nn] -= (ss1 + ss2);
			}
		}
	}
}

//___________________ Viscous_Splitting_LW2 __________________

void CADSolver::Viscous_Splitting_LW(void)
{
	// #pragma omp barrier
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];
	Boundary_u_t();

	// #pragma omp barrier

#pragma omp parallel for // num_threads(NumberThreads)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble ux[8], uy[8], vx[8], vy[8], Tx[8], Ty[8], mu, ka;
		mydouble sux, suy, svx, svy, sTx, sTy, tau11t, tau12t, tau22t, tau11_0, tau12_0, tau22_0, ut_0, vt_0, u_0, v_0;
		long i = LW_Viscous_cd / 2, j = 0;
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			j = 0;
			for (long ww = -i; ww < i + 1; ww++)
			{
				ux[j] = D_ut[iCell + ww][jCell];
				vx[j] = D_vt[iCell + ww][jCell];
				Tx[j] = D_Tt[iCell + ww][jCell];

				uy[j] = D_ut[iCell][jCell + ww];
				vy[j] = D_vt[iCell][jCell + ww];
				Ty[j] = D_Tt[iCell][jCell + ww];
				j = j + 1;
			}

			sux = scheme_ALW_df_dx(ux, LW_Viscous_cd, i) / dx;
			svx = scheme_ALW_df_dx(vx, LW_Viscous_cd, i) / dx;
			sTx = scheme_ALW_df_dx(Tx, LW_Viscous_cd, i) / dx;

			suy = scheme_ALW_df_dx(uy, LW_Viscous_cd, i) / dy;
			svy = scheme_ALW_df_dx(vy, LW_Viscous_cd, i) / dy;
			sTy = scheme_ALW_df_dx(Ty, LW_Viscous_cd, i) / dy;

			mu = D_mu[iCell][jCell];
			ka = D_ka[iCell][jCell];

			tau12t = mu * (suy + svx);
			tau11t = 2. * mu * sux - 2. / 3. * mu * (sux + svy);
			tau22t = 2. * mu * svy - 2. / 3. * mu * (sux + svy);

			tau11_0 = D_gx[iCell][jCell][1];
			tau12_0 = D_gx[iCell][jCell][2];
			tau22_0 = D_gy[iCell][jCell][2];
			ut_0 = D_ut[iCell][jCell];
			vt_0 = D_vt[iCell][jCell];

			u_0 = D_duvpc[iCell][jCell][1];
			v_0 = D_duvpc[iCell][jCell][2];

			D_gx_t[iCell][jCell][1] = tau11t;
			D_gx_t[iCell][jCell][2] = tau12t;
			D_gx_t[iCell][jCell][3] = ka * sTx + ut_0 * tau11_0 + tau11t * u_0 + vt_0 * tau12_0 + tau12t * v_0;

			D_gy_t[iCell][jCell][1] = tau12t;
			D_gy_t[iCell][jCell][2] = tau22t;
			D_gy_t[iCell][jCell][3] = ka * sTy + ut_0 * tau12_0 + tau12t * u_0 + vt_0 * tau22_0 + tau22t * v_0;
		}
	}
	// #pragma omp barrier
	Viscous_Schemes_LW();
}

void CADSolver::Viscous_Schemes_LW(void)
{

#pragma omp parallel for // num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble fx[8], fy[8];
		mydouble ss1 = 0., ss2 = 0., sss;
		long i = LW_Viscous_cd / 2, j = 0;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (long nn = 1; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < i + 1; ww++)
				{
					fx[j] = D_gx_t[iCell + ww][jCell][nn];
					fy[j] = D_gy_t[iCell][jCell + ww][nn];
					j = j + 1;
				}

				ss1 = scheme_ALW_df_dx(fx, LW_Viscous_cd, i) / dx;
				ss2 = scheme_ALW_df_dx(fy, LW_Viscous_cd, i) / dy;

				D_ddu[iCell][jCell][nn] -= (ss1 + ss2);
			}
		}
	}
}

void CADSolver::Viscous_Splitting_LW_cd(void)
{
	// #pragma omp barrier
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];
	Boundary_u_t();

	// #pragma omp barrier

#pragma omp parallel for // num_threads(NumberThreads)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble ux[8], uy[8], vx[8], vy[8], Tx[8], Ty[8], mu, ka;
		mydouble sux, suy, svx, svy, sTx, sTy, tau11t, tau12t, tau22t, tau11_0, tau12_0, tau22_0, ut_0, vt_0, u_0, v_0;
		long i = LW_Viscous_cd / 2, j = 0;
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			j = 0;
			for (long ww = -i; ww < i + 1; ww++)
			{
				ux[j] = D_ut[iCell + ww][jCell];
				vx[j] = D_vt[iCell + ww][jCell];
				Tx[j] = D_Tt[iCell + ww][jCell];

				uy[j] = D_ut[iCell][jCell + ww];
				vy[j] = D_vt[iCell][jCell + ww];
				Ty[j] = D_Tt[iCell][jCell + ww];
				j = j + 1;
			}

			sux = scheme_ALW_df_dx(ux, LW_Viscous_cd, i) / dx;
			svx = scheme_ALW_df_dx(vx, LW_Viscous_cd, i) / dx;
			sTx = scheme_ALW_df_dx(Tx, LW_Viscous_cd, i) / dx;

			suy = scheme_ALW_df_dx(uy, LW_Viscous_cd, i) / dy;
			svy = scheme_ALW_df_dx(vy, LW_Viscous_cd, i) / dy;
			sTy = scheme_ALW_df_dx(Ty, LW_Viscous_cd, i) / dy;

			mu = D_mu[iCell][jCell];
			ka = D_ka[iCell][jCell];

			tau12t = mu * (suy + svx);
			tau11t = 2. * mu * sux - 2. / 3. * mu * (sux + svy);
			tau22t = 2. * mu * svy - 2. / 3. * mu * (sux + svy);

			tau11_0 = D_gx[iCell][jCell][1];
			tau12_0 = D_gx[iCell][jCell][2];
			tau22_0 = D_gy[iCell][jCell][2];
			ut_0 = D_ut[iCell][jCell];
			vt_0 = D_vt[iCell][jCell];

			u_0 = D_duvpc[iCell][jCell][1];
			v_0 = D_duvpc[iCell][jCell][2];

			D_gx_t[iCell][jCell][0] = 0.;
			D_gx_t[iCell][jCell][1] = tau11t;
			D_gx_t[iCell][jCell][2] = tau12t;
			D_gx_t[iCell][jCell][3] = ka * sTx + ut_0 * tau11_0 + tau11t * u_0 + vt_0 * tau12_0 + tau12t * v_0;

			D_gy_t[iCell][jCell][0] = 0.;
			D_gy_t[iCell][jCell][1] = tau12t;
			D_gy_t[iCell][jCell][2] = tau22t;
			D_gy_t[iCell][jCell][3] = ka * sTy + ut_0 * tau12_0 + tau12t * u_0 + vt_0 * tau22_0 + tau22t * v_0;
		}
	}
	// #pragma omp barrier
	Viscous_Schemes_LW_cd();
}

void CADSolver::Viscous_Schemes_LW_cd(void)
{

#pragma omp parallel for // num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble fx[8], fy[8];
		mydouble ss1 = 0., ss2 = 0., sss;
		long i = LW_Viscous_cd / 2, j = 0;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (long nn = 0; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < i + 1; ww++)
				{
					fx[j] = D_utx0[iCell + ww][jCell][nn] - D_gx_t[iCell + ww][jCell][nn];
					fy[j] = D_uty0[iCell][jCell + ww][nn] - D_gy_t[iCell][jCell + ww][nn];
					j = j + 1;
				}

				ss1 = scheme_ALW_df_dx(fx, LW_Viscous_cd, i) / dx;
				ss2 = scheme_ALW_df_dx(fy, LW_Viscous_cd, i) / dy;

				D_ddu[iCell][jCell][nn] = (ss1 + ss2);
			}
		}
	}
}

//___________________ Viscous_Splitting_LW3 __________________

void CADSolver::Viscous_Splitting_LW3(void)
{
	// #pragma omp barrier
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];
	Boundary_u_tt();

#pragma omp parallel for // num_threads(NumberThreads)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble ux[8], uy[8], vx[8], vy[8], Tx[8], Ty[8], mu, ka;
		mydouble sux, suy, svx, svy, sTx, sTy, tau11_0, tau12_0, tau22_0, u_t, v_t, u0, v0;
		mydouble u_tt, v_tt, tau11_t, tau12_t, tau22_t, tau11_tt, tau12_tt, tau22_tt;
		long i = LW_Viscous_cd / 2, j = 0;
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			j = 0;
			for (long ww = -i; ww < i + 1; ww++)
			{
				ux[j] = D_utt[iCell + ww][jCell];
				vx[j] = D_vtt[iCell + ww][jCell];
				Tx[j] = D_Ttt[iCell + ww][jCell];

				uy[j] = D_utt[iCell][jCell + ww];
				vy[j] = D_vtt[iCell][jCell + ww];
				Ty[j] = D_Ttt[iCell][jCell + ww];
				j = j + 1;
			}

			sux = scheme_ALW_df_dx(ux, LW_Viscous_cd, i) / dx;
			svx = scheme_ALW_df_dx(vx, LW_Viscous_cd, i) / dx;
			sTx = scheme_ALW_df_dx(Tx, LW_Viscous_cd, i) / dx;

			suy = scheme_ALW_df_dx(uy, LW_Viscous_cd, i) / dy;
			svy = scheme_ALW_df_dx(vy, LW_Viscous_cd, i) / dy;
			sTy = scheme_ALW_df_dx(Ty, LW_Viscous_cd, i) / dy;

			mu = D_mu[iCell][jCell];
			ka = D_ka[iCell][jCell];

			tau12_tt = mu * (suy + svx);
			tau11_tt = 2. * mu * sux - 2. / 3. * mu * (sux + svy);
			tau22_tt = 2. * mu * svy - 2. / 3. * mu * (sux + svy);

			tau11_0 = D_gx[iCell][jCell][1];
			tau12_0 = D_gx[iCell][jCell][2];
			tau22_0 = D_gy[iCell][jCell][2];

			tau11_t = D_gx_t[iCell][jCell][1];
			tau12_t = D_gx_t[iCell][jCell][2];
			tau22_t = D_gy_t[iCell][jCell][2];

			u_t = D_ut[iCell][jCell];
			v_t = D_vt[iCell][jCell];
			u_tt = D_utt[iCell][jCell];
			v_tt = D_vtt[iCell][jCell];

			u0 = D_duvpc[iCell][jCell][1];
			v0 = D_duvpc[iCell][jCell][2];

			D_gx[iCell][jCell][1] = tau11_tt;
			D_gx[iCell][jCell][2] = tau12_tt;
			D_gx[iCell][jCell][3] = ka * sTx + u_tt * tau11_0 + 2 * tau11_t * u_t + u0 * tau11_tt +
									v_tt * tau12_0 + 2 * tau12_t * v_t + v0 * tau12_tt;

			D_gy[iCell][jCell][1] = tau12_tt;
			D_gy[iCell][jCell][2] = tau22_tt;
			D_gy[iCell][jCell][3] = ka * sTy + u_tt * tau12_0 + 2 * tau12_t * u_t + u0 * tau12_tt +
									v_tt * tau22_0 + 2 * tau22_t * v_t + v0 * tau22_tt;
		}
	}
	// #pragma omp barrier
	Viscous_Schemes_LW3();
}

void CADSolver::Viscous_Schemes_LW3(void)
{

#pragma omp parallel for // num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble fx[8], fy[8];
		mydouble ss1 = 0., ss2 = 0., sss;
		long i = LW_Viscous_cd / 2, j = 0;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (long nn = 1; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < i + 1; ww++)
				{
					fx[j] = D_gx[iCell + ww][jCell][nn];
					fy[j] = D_gy[iCell][jCell + ww][nn];
					j = j + 1;
				}

				ss1 = scheme_ALW_df_dx(fx, LW_Viscous_cd, i) / dx;
				ss2 = scheme_ALW_df_dx(fy, LW_Viscous_cd, i) / dy;

				D_dddu[iCell][jCell][nn] -= (ss1 + ss2);
			}
		}
	}
}

void CADSolver::Viscous_Splitting_LW3_cd(void)
{
	// #pragma omp barrier
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];
	Boundary_u_tt();

#pragma omp parallel for // num_threads(NumberThreads)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble ux[8], uy[8], vx[8], vy[8], Tx[8], Ty[8], mu, ka;
		mydouble sux, suy, svx, svy, sTx, sTy, tau11_0, tau12_0, tau22_0, u_t, v_t, u0, v0;
		mydouble u_tt, v_tt, tau11_t, tau12_t, tau22_t, tau11_tt, tau12_tt, tau22_tt;
		long i = LW_Viscous_cd / 2, j = 0;
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			j = 0;
			for (long ww = -i; ww < i + 1; ww++)
			{
				ux[j] = D_utt[iCell + ww][jCell];
				vx[j] = D_vtt[iCell + ww][jCell];
				Tx[j] = D_Ttt[iCell + ww][jCell];

				uy[j] = D_utt[iCell][jCell + ww];
				vy[j] = D_vtt[iCell][jCell + ww];
				Ty[j] = D_Ttt[iCell][jCell + ww];
				j = j + 1;
			}

			sux = scheme_ALW_df_dx(ux, LW_Viscous_cd, i) / dx;
			svx = scheme_ALW_df_dx(vx, LW_Viscous_cd, i) / dx;
			sTx = scheme_ALW_df_dx(Tx, LW_Viscous_cd, i) / dx;

			suy = scheme_ALW_df_dx(uy, LW_Viscous_cd, i) / dy;
			svy = scheme_ALW_df_dx(vy, LW_Viscous_cd, i) / dy;
			sTy = scheme_ALW_df_dx(Ty, LW_Viscous_cd, i) / dy;

			mu = D_mu[iCell][jCell];
			ka = D_ka[iCell][jCell];

			tau12_tt = mu * (suy + svx);
			tau11_tt = 2. * mu * sux - 2. / 3. * mu * (sux + svy);
			tau22_tt = 2. * mu * svy - 2. / 3. * mu * (sux + svy);

			tau11_0 = D_gx[iCell][jCell][1];
			tau12_0 = D_gx[iCell][jCell][2];
			tau22_0 = D_gy[iCell][jCell][2];

			tau11_t = D_gx_t[iCell][jCell][1];
			tau12_t = D_gx_t[iCell][jCell][2];
			tau22_t = D_gy_t[iCell][jCell][2];

			u_t = D_ut[iCell][jCell];
			v_t = D_vt[iCell][jCell];
			u_tt = D_utt[iCell][jCell];
			v_tt = D_vtt[iCell][jCell];

			u0 = D_duvpc[iCell][jCell][1];
			v0 = D_duvpc[iCell][jCell][2];

			D_gx[iCell][jCell][0] = 0.;
			D_gx[iCell][jCell][1] = tau11_tt;
			D_gx[iCell][jCell][2] = tau12_tt;
			D_gx[iCell][jCell][3] = ka * sTx + u_tt * tau11_0 + 2 * tau11_t * u_t + u0 * tau11_tt +
									v_tt * tau12_0 + 2 * tau12_t * v_t + v0 * tau12_tt;

			D_gy[iCell][jCell][0] = 0.;
			D_gy[iCell][jCell][1] = tau12_tt;
			D_gy[iCell][jCell][2] = tau22_tt;
			D_gy[iCell][jCell][3] = ka * sTy + u_tt * tau12_0 + 2 * tau12_t * u_t + u0 * tau12_tt +
									v_tt * tau22_0 + 2 * tau22_t * v_t + v0 * tau22_tt;
		}
	}
	// #pragma omp barrier
	Viscous_Schemes_LW3_cd();
}

void CADSolver::Viscous_Schemes_LW3_cd(void)
{

#pragma omp parallel for // num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble fx[8], fy[8];
		mydouble ss1 = 0., ss2 = 0., sss;
		long i = LW_Viscous_cd / 2, j = 0;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (long nn = 0; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < i + 1; ww++)
				{
					fx[j] = D_utx0[iCell + ww][jCell][nn] - D_gx[iCell + ww][jCell][nn];
					fy[j] = D_uty0[iCell][jCell + ww][nn] - D_gy[iCell][jCell + ww][nn];
					j = j + 1;
				}

				ss1 = scheme_ALW_df_dx(fx, LW_Viscous_cd, i) / dx;
				ss2 = scheme_ALW_df_dx(fy, LW_Viscous_cd, i) / dy;

				D_dddu[iCell][jCell][nn] = (ss1 + ss2);
			}
		}
	}
}
