// #pragma once
#include "xBurgers.hpp"
using namespace std;

//_________________L-W-Approximate_______________

void CADSolver::comput_du_ALW(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble s_d, s_u, s_p, s_v, ffx[5], ffy[5], fxt[8], fyt[8], local[5], fkt, fxdt, fydt;
		mydouble DT1;
		DT1 = DT * 0.00001;
		// mydouble du1[4], un0[4];
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			for (long kt = 1; kt < 5; kt++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					local[iVar] = D_un[iCell][jCell][iVar] - DT1 * (kt * 1.0 - 2.0) * D_du[iCell][jCell][iVar];
				}

				s_d = local[0];
				s_u = local[1] / s_d;
				s_v = local[2] / s_d;
				s_p = (Gamma - 1.0) * (local[3] - 0.5 * s_d * (s_u * s_u + s_v * s_v));

				ffx[0] = local[1], ffy[0] = local[2];
				ffx[1] = local[1] * s_u + s_p, ffy[1] = local[2] * s_u;
				ffx[2] = local[1] * s_v, ffy[2] = local[2] * s_v + s_p;
				ffx[3] = (local[3] + s_p) * s_u, ffy[3] = (local[3] + s_p) * s_v;
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					D_fxat[iCell][jCell][iVar][kt] = ffx[iVar];
					D_fyat[iCell][jCell][iVar][kt] = ffy[iVar];
				}
			}

			for (short iVar = 0; iVar < nVar; iVar++)
			{
				for (long kt = 1; kt < 5; kt++)
				{
					// fxt[kt] = Element[iCell][jCell]->GetSolution_u(kt, 1, iVar);
					// fyt[kt] = Element[iCell][jCell]->GetSolution_u(kt, 2, iVar);
					fxt[kt] = D_fxat[iCell][jCell][iVar][kt];
					fyt[kt] = D_fyat[iCell][jCell][iVar][kt];
				}

				fxt[0] = (-2. * fxt[1] - 3. * fxt[2] + 6. * fxt[3] - fxt[4]) / 6. / DT1;
				fyt[0] = (-2. * fyt[1] - 3. * fyt[2] + 6. * fyt[3] - fyt[4]) / 6. / DT1;

				// fxt[0] = (-fxt[2] + fxt[3]) / DT;
				// fyt[0] = (-fyt[2] + fyt[3]) / DT;

				// D_fxat[iCell][jCell][iVar][0] = fxt[0];
				// D_fyat[iCell][jCell][iVar][0] = fyt[0];

				D_utx0[iCell][jCell][iVar] = fxt[0];
				D_uty0[iCell][jCell][iVar] = fyt[0];

				// D_utx0[iCell][jCell][iVar] = fxt[0];
				// Element[iCell][jCell]->SetSolution_u(0, 1, iVar, fxt[0]);
				// Element[iCell][jCell]->SetSolution_u(0, 2, iVar, fyt[0]);
			}
		}
	}

	// Boundary_dutxy_LWA_periodic();
	// Boundary_dutxy_all();
}

void CADSolver::scheme_du_dx(void)
{
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];

	// #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble fx[8], fy[8];
		mydouble ss1 = 0., ss2 = 0., sss;
		long i = 3, j = 0;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < 4; ww++)
				{
					fx[j] = D_fxat[iCell + ww][jCell][nn][0];
					fy[j] = D_fyat[iCell][jCell + ww][nn][0];
					j = j + 1;
				}

				ss1 = scheme_ALW_df_dx(fx, SpaceScheme_cd, i) / dx;
				ss2 = scheme_ALW_df_dx(fy, SpaceScheme_cd, i) / dy;
				sss = ss1 + ss2;
				D_ddu[iCell][jCell][nn] = sss;
				// Element[iCell][jCell]->SetSolution_ddu(0, nn, sss);
			}
		}
	}
}

mydouble CADSolver::scheme_ALW_df_dx(mydouble *fz, long scheme_name, long i)
{
	mydouble ss1 = 0.;

	if (scheme_name == 50)
	{ // /*  //weno_5-js 负通量  ture https://doi.org/10.1007/s42967-019-0001-3
		double ep = 1.E-5, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
		double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
		for (int j = 0; j < 2; j++)
		{
			sn0 = 13. / 12. * pow((*(fz + i + 1) - 2. * (*(fz + i + 2)) + *(fz + i + 3)), 2) + 0.25 * pow((3. * (*(fz + i + 1)) - 4. * (*(fz + i + 2)) + *(fz + i + 3)), 2);
			sn1 = 13. / 12. * pow((*(fz + i) - 2. * (*(fz + i + 1)) + *(fz + i + 2)), 2) + 0.25 * pow((*(fz + i) - *(fz + i + 2)), 2);
			sn2 = 13. / 12. * pow((*(fz + i - 1) - 2. * (*(fz + i)) + *(fz + i + 1)), 2) + 0.25 * pow((*(fz + i - 1) - 4. * (*(fz + i)) + 3. * (*(fz + i + 1))), 2);
			az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
			W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);

			q03 = 11. / 6. * (*(fz + i + 1)) - 7. / 6. * (*(fz + i + 2)) + 2. / 6. * (*(fz + i + 3));

			q13 = 2. / 6. * (*(fz + i)) + 5. / 6. * (*(fz + i + 1)) - 1. / 6. * (*(fz + i + 2));

			q23 = -1. / 6. * (*(fz + i - 1)) + 5. / 6. * (*(fz + i)) + 2. / 6. * (*(fz + i + 1));
			vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
			i--;
		}
		ss1 = (vz[0] - vz[1]);
		// */
	}

	if (scheme_name == 51)
	{
		double ep = 1.E-5, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
		double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
		for (int j = 0; j < 2; j++)
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
		ss1 = (vz[0] - vz[1]) / dx;
	}

	switch (scheme_name)
	{

	case 1:
		ss1 = (*(fz + i + 1) - *(fz + i));
		break;

	case 11:
		ss1 = (*(fz + i) - *(fz + i - 1));
		break;

	case 2:
		ss1 = (*(fz + i + 1) - *(fz + i - 1)) / 2.;
		break;

	case 5:
		// 5阶 迎风
		ss1 = (-*(fz + i - 3) * 2. + *(fz + i - 2) * 15. - *(fz + i - 1) * 60. + *(fz + i) * 20. + *(fz + i + 1) * 30. - *(fz + i + 2) * 3.) / 60;
		break;

	case 50:
		// weno_5-js 负通量
		break;

	case 51:
		// weno_5-js 正通量
		break;

	case 4:
		// ss1 = (*(fz + i - 2) - *(fz + i - 1) * 8. + *(fz + i + 1) * 8. - *(fz + i + 2)) / 12.;
		ss1 = (fz[i - 2] - fz[i + 2] + (fz[i + 1] - fz[i - 1]) * 8.) / 12.;
		break;

	case 6:
		ss1 = (-fz[i - 3] + fz[i - 2] * 9. - fz[i - 1] * 45. +
			   fz[i + 1] * 45. - fz[i + 2] * 9. + fz[i + 3]) /
			  60.; // 正通量

		break;

	default:
		cout << "scheme_ALW is not avaliable for this Solver." << endl;
		break;
	}

	return ss1;
}
// _________________L-W-Old_____________________________ 

void CADSolver::comput_du_LW_Old(void)
{
#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble q1 = 0., q2 = 0., q3 = 0., q4 = 0.;
		mydouble Fq[4][4] = {0}, Gq[4][4] = {0}, qt[4] = {0}, ft[4] = {0}, gt[4] = {0};
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			q1 = D_un[iCell][jCell][0];
			q2 = D_un[iCell][jCell][1];
			q3 = D_un[iCell][jCell][2];
			q4 = D_un[iCell][jCell][3];

			Fq[0][0] = 0.;
			Fq[0][1] = 1.;
			Fq[0][2] = 0.;
			Fq[0][3] = 0.;

			Fq[1][0] = (q2 * q2 * (-3. + Gamma) + q3 * q3 * (Gamma_1)) / (2. * q1 * q1);
			Fq[1][1] = -(q2 * (-3. + Gamma) / q1);
			Fq[1][2] = (q3 - q3 * Gamma) / q1;
			Fq[1][3] = Gamma_1;

			Fq[2][0] = -(q2 * q3 / (q1 * q1));
			Fq[2][1] = q3 / q1;
			Fq[2][2] = q2 / q1;
			Fq[2][3] = 0.;

			Fq[3][0] = (q2 * (q2 * q2 * (Gamma_1) + q3 * q3 * (Gamma_1)-q1 * q4 * Gamma)) / (q1 * q1 * q1);
			Fq[3][1] = -((3. * q2 * q2 * (Gamma_1) + q3 * q3 * (Gamma_1)-2. * q1 * q4 * Gamma) / (2. * q1 * q1));
			Fq[3][2] = -(q2 * q3 * (Gamma_1) / (q1 * q1));
			Fq[3][3] = (q2 * Gamma) / q1;

			Gq[0][0] = 0.;
			Gq[0][1] = 0.;
			Gq[0][2] = 1.;
			Gq[0][3] = 0.;

			Gq[1][0] = -((q2 * q3) / (q1 * q1));
			Gq[1][1] = (q3) / q1;
			Gq[1][2] = (q2) / q1;
			Gq[1][3] = 0.;

			Gq[2][0] = (q3 * q3 * (-3. + Gamma) + q2 * q2 * (Gamma_1)) / (2. * q1 * q1);
			Gq[2][1] = (q2 - q2 * Gamma) / q1;
			Gq[2][2] = -((q3 * (-3. + Gamma)) / q1);
			Gq[2][3] = Gamma_1;

			Gq[3][0] = (q3 * (q2 * q2 * (Gamma_1) + q3 * q3 * (Gamma_1)-q1 * q4 * Gamma)) / (q1 * q1 * q1);
			Gq[3][1] = -((q2 * q3 * (Gamma_1)) / (q1 * q1));
			Gq[3][2] = -((q2 * q2 * (Gamma_1) + 3. * q3 * q3 * (Gamma_1)-2. * q1 * q4 * Gamma) / (2. * q1 * q1));
			Gq[3][3] = (q3 * Gamma) / q1;

			qt[0] = -D_du[iCell][jCell][0];
			qt[1] = -D_du[iCell][jCell][1];
			qt[2] = -D_du[iCell][jCell][2];
			qt[3] = -D_du[iCell][jCell][3];

			for (long i = 0; i < 4; i++)
			{
				ft[i] = 0., gt[i] = 0.;
				for (long k = 0; k < 4; k++)
				{
					ft[i] += Fq[i][k] * qt[k];
					gt[i] += Gq[i][k] * qt[k];
				}
			}

			for (short iVar = 0; iVar < 4; iVar++)
			{
				D_utx0[iCell][jCell][iVar] = ft[iVar];
				D_uty0[iCell][jCell][iVar] = gt[iVar];
			}
			// __________out_____________
			// D_out[iCell][jCell][0] = uut;
			// D_out[iCell][jCell][1] = ppt;
		}
	}

	// Boundary_dutxy_all();
}

// _________________L-W-2-NEW_____________________________

void CADSolver::comput_du_LW_New(void)
{
	// #pragma omp barrier

	Boundary_du(); // Boundary  D_ddu

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb - 4; jCell < nCellYo_L + 5; jCell++)
	{
		mydouble q1 = 0., q2 = 0., q3 = 0., q4 = 0., q1t, q2t, q3t, q4t, dd, uu, vv, pp, ddt, uut, vvt, ppt, ss1;
		mydouble ft[4], gt[4];
		for (long iCell = nb - 4; iCell < nCell_L + 5; iCell++)
		{

			q1 = D_un[iCell][jCell][0];
			q2 = D_un[iCell][jCell][1];
			q3 = D_un[iCell][jCell][2];
			q4 = D_un[iCell][jCell][3];

			q1t = -D_du[iCell][jCell][0];
			q2t = -D_du[iCell][jCell][1];
			q3t = -D_du[iCell][jCell][2];
			q4t = -D_du[iCell][jCell][3];

			ddt = q1t;

			dd = D_duvpc[iCell][jCell][0];
			uu = D_duvpc[iCell][jCell][1];
			vv = D_duvpc[iCell][jCell][2];
			pp = D_duvpc[iCell][jCell][3];

			uut = (q2t - ddt * uu) / dd;
			vvt = (q3t - ddt * vv) / dd;
			ss1 = uu * uut + vv * vvt;
			ppt = (q4t - 0.5 * ddt * (uu * uu + vv * vv) - dd * ss1) * Gamma_1;

			if (Viscous_flag == 1)
			{
				// D_Tt[iCell][jCell] = ((q4t * dd - q4 * ddt) / (dd * dd) - ss1) / C_v;
				D_Tt[iCell][jCell] = (ppt * dd - pp * ddt) * Gamma / (dd * dd);
			}
			D_ut[iCell][jCell] = uut;
			D_vt[iCell][jCell] = vvt;
			D_pt[iCell][jCell] = ppt;

			D_utx0[iCell][jCell][0] = q2t;
			D_utx0[iCell][jCell][1] = q2 * uut + q2t * uu + ppt;
			D_utx0[iCell][jCell][2] = q2 * vvt + q2t * vv;
			D_utx0[iCell][jCell][3] = q4 * uut + q4t * uu + pp * uut + ppt * uu;

			D_uty0[iCell][jCell][0] = q3t;
			D_uty0[iCell][jCell][1] = D_utx0[iCell][jCell][2];
			D_uty0[iCell][jCell][2] = q3 * vvt + q3t * vv + ppt;
			D_uty0[iCell][jCell][3] = q4 * vvt + q4t * vv + pp * vvt + ppt * vv;

			// ft[0] = q2t;
			// ft[1] = q2 * uut + q2t * uu + ppt;
			// ft[2] = q2 * vvt + q2t * vv;
			// ft[3] = q4 * uut + q4t * uu + pp * uut + ppt * uu;

			// gt[0] = q3t;
			// gt[1] = ft[2];
			// gt[2] = q3 * vvt + q3t * vv + ppt;
			// gt[3] = q4 * vvt + q4t * vv + pp * vvt + ppt * vv;

			// for (short iVar = 0; iVar < 4; iVar++)
			// {
			// 	D_utx0[iCell][jCell][iVar] = ft[iVar];
			// 	D_uty0[iCell][jCell][iVar] = gt[iVar];
			// }
		}
	}
	// #pragma omp barrier
	// Boundary_dutxy_all();
}

void CADSolver::scheme_du_dx_New(void)
{
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];
	// #pragma omp barrier

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble fx[8], fy[8];
		mydouble ss1 = 0., ss2 = 0., sss;
		long i = 3, j = 0;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (long nn = 0; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < 4; ww++)
				{
					fx[j] = D_utx0[iCell + ww][jCell][nn];
					fy[j] = D_uty0[iCell][jCell + ww][nn];
					j = j + 1;
				}

				ss1 = scheme_ALW_df_dx(fx, SpaceScheme_cd, i) / dx;
				if (Case_1D_or_2D == 2)
					ss2 = scheme_ALW_df_dx(fy, SpaceScheme_cd, i) / dy;
				D_ddu[iCell][jCell][nn] = ss1 + ss2;
				if (positivity_flag_LW == 1)
				{
					// ss1 = scheme_ALW_df_dx(fx, 1, i) / dx;
					// ss2 = scheme_ALW_df_dx(fy, 1, i) / dy;
					// ss1 = (fx[i + 1] - fx[i]) / dx;
					if (Case_1D_or_2D == 2)
						ss2 = (fy[i + 1] - fy[i]) / dy;
					ss1 = 0.5 * (fx[i + 1] - fx[i - 1]) / dx;
					// ss1 = (fx[i] - fx[i - 1]) / dx;
					// ss2 = (fy[i] - fy[i - 1]) / dy;
					D_ddu_up1[iCell][jCell][nn] = ss1 + ss2;
				}
			}
		}
	}
	// #pragma omp barrier
	if (positivity_flag_LW == 1)
		comput_ddu_1up_positivity();
}

// _________________L-W-3-NEW_____________________________

void CADSolver::comput_ddu_LW_New(void)
{
	// #pragma omp barrier

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble q1 = 0., q2 = 0., q3 = 0., q4 = 0., q1t, q2t, q3t, q4t, ss1;
		mydouble q1tt, q2tt, q3tt, q4tt, duu_tt, dvv_tt, duv_tt, ee, ee_t, ee_tt;
		mydouble d0, u0, v0, p0, u_t, v_t, p_t, u_tt, v_tt, p_tt;
		mydouble ft[4], gt[4];
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			q1 = D_un[iCell][jCell][0];
			q2 = D_un[iCell][jCell][1];
			q3 = D_un[iCell][jCell][2];
			q4 = D_un[iCell][jCell][3];

			q1t = -D_du[iCell][jCell][0];
			q2t = -D_du[iCell][jCell][1];
			q3t = -D_du[iCell][jCell][2];
			q4t = -D_du[iCell][jCell][3];

			q1tt = -D_ddu[iCell][jCell][0];
			q2tt = -D_ddu[iCell][jCell][1];
			q3tt = -D_ddu[iCell][jCell][2];
			q4tt = -D_ddu[iCell][jCell][3];

			// D_ddu[iCell][jCell][nn] = ss1 + ss2;
			d0 = D_duvpc[iCell][jCell][0];
			u0 = D_duvpc[iCell][jCell][1];
			v0 = D_duvpc[iCell][jCell][2];
			p0 = D_duvpc[iCell][jCell][3];

			u_t = D_ut[iCell][jCell];
			v_t = D_vt[iCell][jCell];
			p_t = D_pt[iCell][jCell];

			u_tt = (q2tt - q1tt * u0 - 2. * q1t * u_t) / d0;
			v_tt = (q3tt - q1tt * v0 - 2. * q1t * v_t) / d0;
			duu_tt = q2 * u_tt + 2. * q2t * u_t + q2tt * u0; // (d*u*u)_tt
			dvv_tt = q3 * v_tt + 2. * q3t * v_t + q3tt * v0; // (d*v*v)_tt
			duv_tt = q2 * v_tt + 2. * q2t * v_t + q2tt * v0; // (d*u*v)_tt
			p_tt = (q4tt - 0.5 * (duu_tt + dvv_tt)) * Gamma_1;

			ee = q4 + p0, ee_t = q4t + p_t, ee_tt = q4tt + p_tt;

			if (Viscous_flag == 1)
			{
				// D_Ttt[iCell][jCell] = (p_tt*d0-p0*q1tt)    /d0/d0/C_v - (u0*u0+v0*v0)*q1tt/d0)    ;
				D_Ttt[iCell][jCell] = Gamma * ((p_tt * d0 - p0 * q1tt) * d0 - 2. * (p_t * d0 - p0 * q1t) * q1t) / (d0 * d0 * d0);
			}

			D_utt[iCell][jCell] = u_tt;
			D_vtt[iCell][jCell] = v_tt;

			D_utx0[iCell][jCell][0] = q2tt;
			D_utx0[iCell][jCell][1] = duu_tt + p_tt;
			D_utx0[iCell][jCell][2] = duv_tt;
			D_utx0[iCell][jCell][3] = ee_tt * u0 + 2. * ee_t * u_t + ee * u_tt;

			D_uty0[iCell][jCell][0] = q3tt;
			D_uty0[iCell][jCell][1] = duv_tt;
			D_uty0[iCell][jCell][2] = dvv_tt + p_tt;
			D_uty0[iCell][jCell][3] = ee_tt * v0 + 2. * ee_t * v_t + ee * v_tt;
		}
	}
	// #pragma omp barrier
	// Boundary_dutxy_all();
}

void CADSolver::scheme_ddu_dx_New(void)
{
	// mydouble *fx, *fy;
	// fx = new mydouble[10], fy = new mydouble[10];
	// #pragma omp barrier

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble fx[8], fy[8];
		mydouble ss1 = 0., ss2 = 0., sss;
		long i = 3, j = 0;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (long nn = 0; nn < 4; nn++)
			{
				j = 0;
				for (long ww = -i; ww < 4; ww++)
				{
					fx[j] = D_utx0[iCell + ww][jCell][nn];
					fy[j] = D_uty0[iCell][jCell + ww][nn];
					j = j + 1;
				}

				ss1 = scheme_ALW_df_dx(fx, SpaceScheme_cd, i) / dx;
				if (Case_1D_or_2D == 2)
					ss2 = scheme_ALW_df_dx(fy, SpaceScheme_cd, i) / dy;
				D_dddu[iCell][jCell][nn] = ss1 + ss2;
				if (positivity_flag_LW == 1)
				{
					// ss1 = scheme_ALW_df_dx(fx, 1, i) / dx;
					// ss2 = scheme_ALW_df_dx(fy, 1, i) / dy;
					// ss1 = (fx[i + 1] - fx[i]) / dx;
					if (Case_1D_or_2D == 2)
						ss2 = (fy[i + 1] - fy[i]) / dy;
					// ss1 = (fx[i] - fx[i - 1]) / dx;
					ss1 = 0.5 * (fx[i + 1] - fx[i - 1]) / dx;
					// ss2 = (fy[i] - fy[i - 1]) / dy;
					D_dddu_up1[iCell][jCell][nn] = ss1 + ss2;
				}
			}
		}
	}
	// #pragma omp barrier
	if (positivity_flag_LW == 1)
		comput_dddu_1up_positivity();
}

// _________________1up_保正性_positivity-preserving____________________________

void CADSolver::comput_du_1up_positivity0(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble q1 = 0., q2 = 0., q3 = 0., q4 = 0., q1t, q2t, q3t, q4t, dd, uu, vv, pp, ddt, uut, vvt, ppt, d_n1, p_n1, du0, du1;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			q1t = -D_du[iCell][jCell][0];
			q2t = -D_du[iCell][jCell][1];
			q3t = -D_du[iCell][jCell][2];
			q4t = -D_du[iCell][jCell][3];

			ddt = q1t;

			dd = D_duvpc[iCell][jCell][0];
			uu = D_duvpc[iCell][jCell][1];
			vv = D_duvpc[iCell][jCell][2];
			pp = D_duvpc[iCell][jCell][3];

			uut = (q2t - ddt * uu) / dd;
			vvt = (q3t - ddt * vv) / dd;
			ppt = (q4t - 0.5 * ddt * (uu * uu + vv * vv) - dd * (uu * uut + vv * vvt)) * Gamma_1;

			d_n1 = dd + DT * ddt;
			p_n1 = pp + DT * ppt;
			if (d_n1 < 1E-15 || p_n1 < 1E-15)
			{
				for (short iVar = 0; iVar < 1; iVar++)
				{
					for (short i = 0; i < 1; i++)
					{
						du0 = D_du[iCell + i][jCell][iVar];
						du1 = D_du_up1[iCell + i][jCell][iVar];
						D_du[iCell + i][jCell][iVar] = 0.2 * du0 + 0.8 * du1;
					}
				}
				for (short iVar = 1; iVar < nVar; iVar++)
				{
					for (short i = 0; i < 1; i++)
					{
						du0 = D_du[iCell + i][jCell][iVar];
						du1 = D_du_up1[iCell + i][jCell][iVar];
						D_du[iCell + i][jCell][iVar] = 0.2 * du0 + 0.8 * du1;
					}
				}
			}
		}
	}
}

void CADSolver::comput_du_1up_positivity(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb - 3; jCell < nCellYo_L + 4; jCell++)
	{
		mydouble dd, uu, vv, pp, du0, du1;
		mydouble uuu[4];
		for (long iCell = nb - 3; iCell < nCell_L + 4; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{
				uuu[iVar] = D_un[iCell][jCell][iVar] - DT * D_du[iCell][jCell][iVar];
			}

			dd = uuu[0], uu = uuu[1] / dd;
			vv = uuu[2] / dd;
			pp = (Gamma_1) * (uuu[3] - 0.5 * dd * (uu * uu + vv * vv));

			if (dd < 1.E-5 || pp < 1.E-5)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du0 = D_du[iCell][jCell][iVar];
					du1 = D_du_up1[iCell][jCell][iVar];
					D_du[iCell][jCell][iVar] = 0.15 * du0 + 0.85 * du1;  //Blast_wave
					// D_du[iCell][jCell][iVar] = 0.00 * du0 + 1. * du1;
				}
			}
		}
	}
}

void CADSolver::comput_ddu_1up_positivity(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble dd, uu, vv, pp, du0, du1;
		mydouble uuu[4];
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{
				uuu[iVar] = D_un[iCell][jCell][iVar] - DT * D_du[iCell][jCell][iVar] - 0.5 * DT * DT * D_ddu[iCell][jCell][iVar];
			}

			dd = uuu[0], uu = uuu[1] / dd;
			vv = uuu[2] / dd;
			pp = (Gamma_1) * (uuu[3] - 0.5 * dd * (uu * uu + vv * vv));

			if (dd < 1.E-6 || pp < 1.E-6)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					// du0 = D_ddu[iCell + i][jCell][iVar];
					// du1 = D_ddu_up1[iCell + i][jCell][iVar];
					D_ddu[iCell][jCell][iVar] = D_ddu_up1[iCell][jCell][iVar];
				}
			}
		}
	}
}

void CADSolver::comput_dddu_1up_positivity(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble dd, uu, vv, pp, du0, du1;
		mydouble uuu[4];
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{
				uuu[iVar] = D_un[iCell][jCell][iVar] - DT * D_du[iCell][jCell][iVar] - 0.5 * DT * DT * D_ddu[iCell][jCell][iVar] - DT * DT * DT * D_dddu[iCell][jCell][iVar] / 6.0;
			}

			dd = uuu[0], uu = uuu[1] / dd;
			vv = uuu[2] / dd;
			pp = (Gamma_1) * (uuu[3] - 0.5 * dd * (uu * uu + vv * vv));

			if (dd < 1.E-6 || pp < 1.E-6)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					// du0 = D_ddu[iCell + i][jCell][iVar];
					// du1 = D_ddu_up1[iCell + i][jCell][iVar];
					D_dddu[iCell][jCell][iVar] = D_dddu_up1[iCell][jCell][iVar];
				}
			}
		}
	}
}
