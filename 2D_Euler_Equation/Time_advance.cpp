// #pragma once
#include "xBurgers.hpp"
using namespace std;

void CADSolver::TimeIntegration_RK11(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, un0, un1, du0, un2;

	// #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble local;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (short iVar = 0; iVar < nVar; iVar++)
			{

				du = D_du[iCell][jCell][iVar];
				un0 = D_un[iCell][jCell][iVar];

				local = un0 - DT * du;

				D_un[iCell][jCell][iVar] = local;
			}
		}
	}
}

void CADSolver::TimeIntegration(unsigned short iRK_Step) // RK33
{

	mydouble local;
	mydouble du, un0, un1, du0, un2;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un[iCell][jCell][iVar];

					local = un0 - DT * du;

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = un0;
					D_du0[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];
					un1 = D_un[iCell][jCell][iVar];

					local = 0.75 * un0 + 0.25 * (un1 - DT * du);

					D_un1[iCell][jCell][iVar] = un1;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];
					un1 = D_un1[iCell][jCell][iVar];
					un2 = D_un[iCell][jCell][iVar];
					local = 1.0 / 3.0 * un0 + 2.0 / 3.0 * (un2 - DT * du);
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_LW3_ALW()
{

	// // #pragma omp barrier

	// #pragma omp parallel for private(a21, v2, v1, vv1, vv2, aa21)
	// mydouble du, ddu, dddu, uu, uu0, local;
	// mydouble DT2 = DT * DT * 0.5, DT3 = DT * DT * DT / 6.;
	// mydouble a21, v2, v1, vv1, vv2, aa21;
	// a21 = 0.594223212099088, v2 = 0.306027487008159;
	// a21 = 2. / 3., v2 = 3. / 8., aa21 = a21 * a21 * 0.5;
	// v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;
#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble du, ddu, dddu, uu, uu0, local;
		mydouble DT2 = DT * DT * 0.5, DT3 = DT * DT * DT / 6.;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (short iVar = 0; iVar < nVar; iVar++)
			{

				du = D_du[iCell][jCell][iVar];
				ddu = D_ddu[iCell][jCell][iVar];
				dddu = D_dddu[iCell][jCell][iVar];
				uu0 = D_un[iCell][jCell][iVar];

				local = uu0 - DT * du - DT2 * ddu - DT3 * dddu;
				// local = uu0 - DT * du - DT2 * ddu;
				D_un[iCell][jCell][iVar] = local;
			}
		}
	}

	// // #pragma omp barrier
}

void CADSolver::TimeIntegration_TDMS23_ALW(unsigned short iRK_Step)
{

	mydouble du, ddu, uu0, Tu1, TL_n1, TLt_n1, local;
	mydouble h_b1, h_b2, h_v1, h_v2, h_vv1, h_vv2;
	// mydouble a21, v2, v1, vv1, vv2, b0, b1, b2, v0, vv0;
	// v0 = 0.8, b1 = -0.2;
	// v1 = 1 - b1 - v0;
	// vv0 = (7 - b1 - 3 * v0) / 6.;
	// vv1 = (2 - 2 * b1 - 3 * v0) / 6.;
	// b0 = -1 - b1;
	// vv2 = 0.0, v2 = 0., b2 = 0.;

	h_b1 = 0.4091 * pow(DT_w1, -1.861); // h_b1 = 0.;  //Tylor
	h_b2 = 1.0 - h_b1;
	h_vv1 = 0.0;
	h_v1 = 1. / 3. * (h_b1 + 1. / pow(DT_w1, 3));
	h_v2 = 1. - 1. / (3.0 * DT_w1 * DT_w1) + (2.0 * h_b1 * DT_w1) / 3.0;
	h_vv2 = (2. + 3.0 * DT_w1 - h_b1 * pow(DT_w1, 3)) / (6.0 * DT_w1);

	// #pragma omp parallel for private(du, ddu, uu0, Tu1, TL_n1, TLt_n1, local, h_b1, h_b2, h_v1, h_v2, h_vv1, h_vv2)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{

				uu0 = D_un[iCell][jCell][iVar];
				du = D_du[iCell][jCell][iVar];
				ddu = D_ddu[iCell][jCell][iVar];

				Tu1 = D_un0[iCell][jCell][iVar];
				TL_n1 = D_du0[iCell][jCell][iVar];
				TLt_n1 = D_ddu0[iCell][jCell][iVar];

				local = h_b1 * Tu1 + h_b2 * uu0 -
						DT * (h_v1 * DT_w1 * TL_n1 + h_v2 * du) -
						DT * DT * (h_vv1 * DT_w1 * DT_w1 * TLt_n1 + h_vv2 * ddu);

				D_un[iCell][jCell][iVar] = local;
				D_un0[iCell][jCell][iVar] = uu0;
				D_du0[iCell][jCell][iVar] = du;
				D_ddu0[iCell][jCell][iVar] = ddu;
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMS34_ALW(unsigned short iRK_Step) // Tylor SSP
{

	mydouble local;

	mydouble a21, v2, v1, vv1, vv2, b0, b1, b2, b3, v0, vv0, uu, vv, w1, w2, v3, vv3, w11, w22, ww1, ww2;

	w2 = DT_w1, w1 = DT_w2;
	if (w1 > 1.3)
		w1 = 1.3;
	if (w2 > 1.3)
		w2 = 1.3;
	b1 = -0.283 + 0.088 * w1 * w1 + 6.246 * exp(-1.227 * w1 - 1.354 * w2);
	b2 = 0.205 + 0.3156 * pow(w1, -1.508) * pow(w2, 1.065) - 3.1 * exp(-2.123 * w1);

	v1 = 22.23 + 8.837 * w1 * w1 - 14.05 * exp(0.7255 * w1 - 0.1072 * w2) - 9.275 * sin(0.5058 * w2);

	b3 = 1. - b1 - b2;
	vv2 = 0.;
	uu = w1 + w2, vv = (b1 + b2) * w2;
	vv1 = (3. + 12. * pow(uu, 2) * v1 * pow(w1, 2) - b1 * pow(uu, 3) * (3. * w1 - w2) + 4. * w2 + b2 * pow(w2, 4)) / (12. * uu * pow(w1, 2) * (3. * w1 + w2));

	v2 = (1. + w2 + w1 * (1. + pow(uu, 3) * v1 - 6. * pow(uu, 2) * vv1 * w1 + b2 * pow(w2, 3))) / ((3. * w1 - w2) * pow(w2, 3));

	vv3 = (3. + (b1 - 12. * vv1) * pow(w1, 4) + 2. * (b1 - 6. * vv1) * pow(w1, 3) * w2 + 2. * w1 * (2. + 3. * w2 - vv * pow(w2, 2)) + w2 * (8. + 6. * w2 - vv * pow(w2, 2))) / (12. * uu * w2);

	v3 = (-1. - 2. * uu - (b1 - 2. * v1) * pow(w1, 4) + 2. * (3. + 2. * vv) * w1 * pow(w2, 2) + (2. + vv) * pow(w2, 3) + pow(w1, 3) * (-4. * b1 * w2 + 6. * v1 * w2)) / (2. * pow(w2, 2) * (3. * w1 + w2));
	w2 = DT_w1, w1 = DT_w2;
	// #pragma omp parallel for private(vv1, vv2, vv3, b1, b2, b3, v2, v1, v3, w1, w2, local)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		mydouble du, ddu, du0, ddu0, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{
				uu0 = D_un[iCell][jCell][iVar];
				du = D_du[iCell][jCell][iVar];
				ddu = D_ddu[iCell][jCell][iVar];

				Tu1 = D_un1[iCell][jCell][iVar];
				TL_n1 = D_dun1[iCell][jCell][iVar];
				TLt_n1 = D_ddun1[iCell][jCell][iVar];

				Tu2 = D_un2[iCell][jCell][iVar];
				TL_n2 = D_dun2[iCell][jCell][iVar];
				TLt_n2 = D_ddun2[iCell][jCell][iVar];

				local = b1 * Tu2 + b2 * Tu1 + b3 * uu0 -
						DT * (v1 * w1 * TL_n2 + v2 * w2 * TL_n1 + v3 * du) -
						DT * DT * (vv1 * w1 * w1 * TLt_n2 + vv2 * w2 * w2 * TLt_n1 + vv3 * ddu);

				D_un[iCell][jCell][iVar] = local;

				D_un2[iCell][jCell][iVar] = D_un1[iCell][jCell][iVar];
				D_dun2[iCell][jCell][iVar] = D_dun1[iCell][jCell][iVar];
				D_ddun2[iCell][jCell][iVar] = D_ddun1[iCell][jCell][iVar];

				D_un1[iCell][jCell][iVar] = uu0;
				D_dun1[iCell][jCell][iVar] = du;
				D_ddun1[iCell][jCell][iVar] = ddu;
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMS34_ALW_1(unsigned short iRK_Step) // SD SSP
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, v2, v1, vv1, vv2, b0, b1, b2, v0, vv0;
	mydouble h_b1, h_b2, h_v1, h_v2, h_vv1, h_vv2, h_b3, h_v3, h_vv3, h_w1, h_w2;

	h_w2 = DT_w1, h_w1 = DT_w2;
	// if (h_w1 > 1.3)
	// 	h_w1 = 1.3;
	// if (h_w2 > 1.3)
	// 	h_w2 = 1.3;
	h_b1 = 22.27 - 22.26 * pow(h_w1, 0.00209) * exp(-0.0006788 * h_w1) * pow(h_w2, -0.000214);
	h_b2 = 0.06399 + 0.2848 * pow(h_w1, -0.2672) * pow(h_w2, -1.599);

	h_b3 = 1. - h_b1 - h_b2, h_v1 = 0., h_vv2 = 0.;

	h_vv1 = (3. + 4. * h_w2 + h_b2 * pow(h_w2, 4) - h_b1 * (3. * h_w1 - h_w2) * pow(h_w1 + h_w2, 3)) /
			(12. * h_w1 * h_w1 * (h_w1 + h_w2) * (3. * h_w1 + h_w2));
	h_v3 = (-2 - 3 * h_w2 + pow(h_w2, 3) + 6 * h_vv1 * h_w1 * h_w1 * (h_w1 + h_w2) * (4 * h_w1 + h_w2) + h_b1 * h_w1 * h_w1 * (h_w1 + h_w2) * (2 * h_w1 + 3 * h_w2)) / pow(h_w2, 3);
	h_v2 = (-1 - h_w2 + h_b1 * h_w1 * pow(h_w1 + h_w2, 3) + 6 * h_vv1 * h_w1 * h_w1 * (h_w1 + h_w2) * (2 * h_w1 + h_w2)) / pow(h_w2, 4);
	h_vv3 = (-h_vv1) * h_w1 * h_w1 + h_v2 * h_w2 * h_w2 + 0.5 * (1. - h_b1 * (h_w1 + h_w2) * (h_w1 + h_w2) - h_b2 * h_w2 * h_w2);

	h_w2 = DT_w1, h_w1 = DT_w2;

	// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, v2, v1, vv1, vv2, b0, b1, b2, v0, vv0, local)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{

			for (short iVar = 0; iVar < nVar; iVar++)
			{
				uu0 = D_un[iCell][jCell][iVar];
				du = D_du[iCell][jCell][iVar];
				ddu = D_ddu[iCell][jCell][iVar];

				Tu1 = D_un1[iCell][jCell][iVar];
				TL_n1 = D_dun1[iCell][jCell][iVar];
				TLt_n1 = D_ddun1[iCell][jCell][iVar];

				Tu2 = D_un2[iCell][jCell][iVar];
				TL_n2 = D_dun2[iCell][jCell][iVar];
				TLt_n2 = D_ddun2[iCell][jCell][iVar];

				local = h_b1 * Tu2 + h_b2 * Tu1 + h_b3 * uu0 -
						DT * (h_v1 * h_w1 * TL_n2 + h_v2 * h_w2 * TL_n1 + h_v3 * du) -
						DT * DT * (h_vv1 * h_w1 * h_w1 * TLt_n2 + h_vv2 * h_w2 * h_w2 * TLt_n1 + h_vv3 * ddu);

				D_un[iCell][jCell][iVar] = local;

				D_un2[iCell][jCell][iVar] = D_un1[iCell][jCell][iVar];
				D_dun2[iCell][jCell][iVar] = D_dun1[iCell][jCell][iVar];
				D_ddun2[iCell][jCell][iVar] = D_ddun1[iCell][jCell][iVar];

				D_un1[iCell][jCell][iVar] = uu0;
				D_dun1[iCell][jCell][iVar] = du;
				D_ddun1[iCell][jCell][iVar] = ddu;
			}
		}
	}
}

/*

void CADSolver::TimeIntegration_TDMSRK223_ALW0(unsigned short iRK_Step)
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1;
	DT_w1 = 1.0;
	w1 = DT_w1;
	a22 = 0.5416 * pow(w1, -0.04681);
	b1 = 0.1397 * exp(0.8707 * w1) - 0.2405 * pow(w1, 1.319);
	v1 = 0.6519 * exp(-2.359 * w1) + 0.01983 * pow(w1, 1.807);
	v3 = 0.3203 * exp(0.4813 * w1) * pow(w1, -0.2711);

	d31 = 0., vv1 = 0., aa21 = 0., a21 = 0.;
	b2 = 1. - b1, d32 = 1. - d31;
	v2 = 1. - v3 + b1 * w1 - v1 * w1;
	aa22 = a22 * a22 * 0.5;
	vv3 = (1. - 3. * a22 * a22 * v3 + (b1 - 3. * v1) * w1 * w1 * w1) / (6. * a22);
	vv2 = -((1. + 3. * a22 * (-1. + a22 * v3) + 3. * a22 * (b1 - 2. * v1) * w1 * w1 + (b1 - 3. * v1) * w1 * w1 * w1) / (6. * a22));

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = d31 * Tu1 + d32 * uu0 - DT * (a21 * w1 * TL_n1 + a22 * du) -
							DT * DT * (aa21 * w1 * w1 * TLt_n1 + aa22 * ddu);

					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = b1 * Tu1 + b2 * uu0 - DT * (v1 * w1 * TL_n1 + v2 * du0 + v3 * du) -
							DT * DT * (vv1 * w1 * w1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu);

					D_un[iCell][jCell][iVar] = local;

					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
				}
			}
		}
	}
}
*/
void CADSolver::TimeIntegration_TDMSRK223_ALW(unsigned short iRK_Step) // TDTSRK23
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1, w2;
	mydouble TL_ns1, TLt_ns1;
	// DT_w1 = 1.0;
	// w1 = DT_w1;
	// a22 = 0.5416 * pow(w1, -0.04681);
	// b1 = 0.1397 * exp(0.8707 * w1) - 0.2405 * pow(w1, 1.319);
	// v1 = 0.6519 * exp(-2.359 * w1) + 0.01983 * pow(w1, 1.807);
	// v3 = 0.3203 * exp(0.4813 * w1) * pow(w1, -0.2711);

	// d31 = 0., vv1 = 0., aa21 = 0., a21 = 0.;
	// b2 = 1. - b1, d32 = 1. - d31;
	// v2 = 1. - v3 + b1 * w1 - v1 * w1;
	// aa22 = a22 * a22 * 0.5;
	// vv3 = (1. - 3. * a22 * a22 * v3 + (b1 - 3. * v1) * w1 * w1 * w1) / (6. * a22);
	// vv2 = -((1. + 3. * a22 * (-1. + a22 * v3) + 3. * a22 * (b1 - 2. * v1) * w1 * w1 + (b1 - 3. * v1) * w1 * w1 * w1) / (6. * a22));

	a21 = 0.5321899654552226, v2 = 0.3890659867744687;
	aa21 = a21 * a21 * 0.5;
	ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.1;

	v1 = 1. - v2 - w1 - w2;
	vv1 = -((1. - 3. * w1 - 3. * w2 + 6. * ww1 + 3. * a21 * (-1. - 2. * w1 + a21 * (v2 + w2) + 2. * ww1) + 6. * ww2) / (6. * a21));
	vv2 = (1. - 3. * w1 - 3. * w2 + 6. * ww1 + 6. * ww2 - 3. * a21 * (-2. * w2 + a21 * (v2 + w2) + 2. * ww2)) / (6. * a21);

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					// Tu1 = D_un1[iCell][jCell][iVar];
					// TL_n1 = D_dun1[iCell][jCell][iVar];
					// TLt_n1 = D_ddun1[iCell][jCell][iVar];

					// local = d31 * Tu1 + d32 * uu0 - DT * (a21 * w1 * TL_n1 + a22 * du) -
					// 		DT * DT * (aa21 * w1 * w1 * TLt_n1 + aa22 * ddu);
					local = uu0 - DT * (a21 * du) - DT * DT * (aa21 * ddu);
					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					TL_ns1 = D_duns1[iCell][jCell][iVar];
					TLt_ns1 = D_dduns1[iCell][jCell][iVar];

					// local = b1 * Tu1 + b2 * uu0 - DT * (v1 * w1 * TL_n1 + v2 * du0 + v3 * du) -
					// 		DT * DT * (vv1 * w1 * w1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu);

					local = uu0 - DT * (v1 * du0 + v2 * du + w1 * TL_n1 + w2 * TL_ns1) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu + ww1 * TLt_n1 + ww2 * TLt_ns1);

					D_un[iCell][jCell][iVar] = local;

					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
					D_duns1[iCell][jCell][iVar] = du;
					D_dduns1[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
}

/*

void CADSolver::TimeIntegration_TDMSRK224_ALW(unsigned short iRK_Step)
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1, a21_d31, w1_2, w1_3, w1_4;
	DT_w1 = 1.0;
	w1 = DT_w1;
	a21 = 0.07982 * pow(w1, -2.688) + 0.01096;
	b1 = 0.06514 * pow(w1, -1.418) + 0.0125;
	d31 = 0.06017 * pow(w1, -1.835) + 0.02575;
	v3 = 0.3985 * exp(-0.2794 * w1) * pow(w1, 0.2483);
	vv3 = 0.1651 * exp(0.1884 * w1) * pow(w1, -0.2552);

	b2 = 1. - b1, d32 = 1. - d31, aa21 = 0.0;
	a21_d31 = pow(3. * a21 - d31, 1. / 3.) * w1;
	w1_2 = w1 * w1, w1_3 = w1_2 * w1, w1_4 = w1_3 * w1;

	a22 = w1 * (d31 - a21) + a21_d31;
	aa22 = 0.5 * ((2. * a21 - d31) * w1 * w1 + a21_d31 * a21_d31);
	v1 = (1. + 2. * w1 + b1 * w1_4 - 2. * a21_d31 * (6. * vv3 * w1 + a21_d31 * (6. * vv3 + 3. * v3 * w1 + 2. * v3 * a21_d31))) / (2. * w1_4);
	v2 = 1. - v3 + (b1 - v1) * w1;
	vv1 = v1 * 0.5 - (1. + b1 * w1_3 - 6. * vv3 * a21_d31 - 3. * v3 * a21_d31 * a21_d31) / (6. * w1_3);
	vv2 = vv1 * w1_2 + (2. + (3. - 6. * vv3) * w1 - b1 * w1_3 - 6. * a21_d31 * (2. * vv3 + v3 * w1 + v3 * a21_d31)) / (6. * w1);

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = d31 * Tu1 + d32 * uu0 - DT * (a21 * w1 * TL_n1 + a22 * du) -
							DT * DT * (aa21 * w1 * w1 * TLt_n1 + aa22 * ddu);

					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = b1 * Tu1 + b2 * uu0 - DT * (v1 * w1 * TL_n1 + v2 * du0 + v3 * du) -
							DT * DT * (vv1 * w1 * w1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu);

					D_un[iCell][jCell][iVar] = local;

					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
				}
			}
		}
	}
}
*/
void CADSolver::TimeIntegration_TDMSRK224_ALW(unsigned short iRK_Step) // TDTSRK24
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1, w2;
	mydouble TL_ns1, TLt_ns1;
	// DT_w1 = 1.0;
	// w1 = DT_w1;
	// a22 = 0.5416 * pow(w1, -0.04681);
	// b1 = 0.1397 * exp(0.8707 * w1) - 0.2405 * pow(w1, 1.319);
	// v1 = 0.6519 * exp(-2.359 * w1) + 0.01983 * pow(w1, 1.807);
	// v3 = 0.3203 * exp(0.4813 * w1) * pow(w1, -0.2711);

	w2 = -0.1, ww2 = 0.0, a21 = 0.4680145029983404, v1 = 1.0;
	v2 = -w2, w1 = 1. - v1;
	vv1 = -((1. - 14. * a21 + 2 * v1 + 6. * a21 * v1 - 2 * w2 - 6. * a21 * w2 + 12. * a21 * ww2) / (12. * a21));
	vv2 = -((-1. - 2 * v1 + 2 * w2 - 12. * a21 * a21 * w2 - 12. * a21 * ww2 + 12. * a21 * a21 * ww2) / (12. * a21 * (1. + a21)));
	ww1 = -((-5. - 4. * a21 + 4. * v1 + 6. * a21 * v1 - 4. * w2 + 6. * a21 * w2 + 12. * ww2 - 12. * a21 * ww2) / (12. * (1. + a21)));

	aa21 = a21 * a21 * 0.5;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					// Tu1 = D_un1[iCell][jCell][iVar];
					// TL_n1 = D_dun1[iCell][jCell][iVar];
					// TLt_n1 = D_ddun1[iCell][jCell][iVar];

					// local = d31 * Tu1 + d32 * uu0 - DT * (a21 * w1 * TL_n1 + a22 * du) -
					// 		DT * DT * (aa21 * w1 * w1 * TLt_n1 + aa22 * ddu);
					local = uu0 - DT * (a21 * du) - DT * DT * (aa21 * ddu);
					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					TL_ns1 = D_duns1[iCell][jCell][iVar];
					TLt_ns1 = D_dduns1[iCell][jCell][iVar];

					// local = b1 * Tu1 + b2 * uu0 - DT * (v1 * w1 * TL_n1 + v2 * du0 + v3 * du) -
					// 		DT * DT * (vv1 * w1 * w1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu);

					local = uu0 - DT * (v1 * du0 + v2 * du + w1 * TL_n1 + w2 * TL_ns1) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu + ww1 * TLt_n1 + ww2 * TLt_ns1);

					D_un[iCell][jCell][iVar] = local;

					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
					D_duns1[iCell][jCell][iVar] = du;
					D_dduns1[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
}

/*

void CADSolver::TimeIntegration_TDMSRK235_ALW(unsigned short iRK_Step)
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1, a21_d31, w1_2, w1_3, w1_4;

	mydouble a32, a31, aa31, aa32, d33, b3, a33, v4, aa33, vv4;

	mydouble arr[20] = {0.014767605406112, 0.167600981747149, 0.817631412846739, 0.030575274114486, 0.133045703046640,
						0.836379022838873, 0.000013886700706, 0.199211349061870, 0.569735915417622, 0.027032987291601, 0.105732896401668,
						0.780704117613404, 0.280726249968939, 0.010423447853488, 0.000000000000000, 0.238971864574849, 0.005532322014941,
						0.031144942360066, 0.126836083983140, 0.208085846602162}; //,(0.841322457462487)

	d31 = arr[0], d32 = arr[1], d33 = arr[2], b1 = arr[3], b2 = arr[4], b3 = arr[5];
	a31 = arr[6], a32 = arr[7], a33 = arr[8];
	v1 = arr[9], v2 = arr[10], v3 = arr[11], v4 = arr[12];
	aa31 = arr[13], aa32 = arr[14], aa33 = arr[15];
	vv1 = arr[16], vv2 = arr[17], vv3 = arr[18], vv4 = arr[19];

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					Tu2 = D_un2[iCell][jCell][iVar];
					TL_n2 = D_dun2[iCell][jCell][iVar];
					TLt_n2 = D_ddun2[iCell][jCell][iVar];

					local = d31 * Tu2 + d32 * Tu1 + d33 * uu0 - DT * (a31 * TL_n2 + a32 * TL_n1 + a33 * du) -
							DT * DT * (aa31 * TLt_n2 + aa32 * TLt_n1 + aa33 * ddu);

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					Tu2 = D_un2[iCell][jCell][iVar];
					TL_n2 = D_dun2[iCell][jCell][iVar];
					TLt_n2 = D_ddun2[iCell][jCell][iVar];

					local = b1 * Tu2 + b2 * Tu1 + b3 * uu0 - DT * (v1 * TL_n2 + +v2 * TL_n1 + v3 * du0 + v4 * du) -
							DT * DT * (vv1 * TLt_n2 + vv2 * TLt_n1 + vv3 * ddu0 + vv4 * ddu);

					D_un[iCell][jCell][iVar] = local;

					D_un2[iCell][jCell][iVar] = D_un1[iCell][jCell][iVar];
					D_dun2[iCell][jCell][iVar] = D_dun1[iCell][jCell][iVar];
					D_ddun2[iCell][jCell][iVar] = D_ddun1[iCell][jCell][iVar];

					// D_un1[iCell][jCell][iVar] = uu0;
					// D_dun1[iCell][jCell][iVar] = du0;
					// D_ddun1[iCell][jCell][iVar] = ddu0;

					D_un1[iCell][jCell][iVar] = D_un0[iCell][jCell][iVar];
					D_dun1[iCell][jCell][iVar] = D_du0[iCell][jCell][iVar];
					D_ddun1[iCell][jCell][iVar] = D_ddu0[iCell][jCell][iVar];
				}
			}
		}
	}
}

*/

void CADSolver::TimeIntegration_TDMSRK235_ALW(unsigned short iRK_Step) // TDTSRK25
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, ww1, ww2, w11, w22, v1, v2, v3, vv1, vv2, vv3, d31, d32, w1, w2;
	mydouble TL_ns1, TLt_ns1;

	a21 = 0.7650141887498161; // 1. / 210. * (93. + sqrt(2139.));
	vv1 = (-31. + 6. * a21 * (-31. + 85. * a21)) / (360. * pow(a21, 2));
	vv2 = 31. / (360. * pow(a21, 2));
	v1 = -(1. / 2.) + 31. / (30. * a21), v2 = 0.;
	w1 = 3. / 2. - 31. / (30. * a21), w2 = 0.;
	ww1 = (31. + 6. * a21 * (-31. + 35. * a21)) / (360. * pow(a21, 2));
	ww2 = -(31. / (360. * pow(a21, 2)));

	aa21 = a21 * a21 * 0.5;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					// Tu1 = D_un1[iCell][jCell][iVar];
					// TL_n1 = D_dun1[iCell][jCell][iVar];
					// TLt_n1 = D_ddun1[iCell][jCell][iVar];

					// local = d31 * Tu1 + d32 * uu0 - DT * (a21 * w1 * TL_n1 + a22 * du) -
					// 		DT * DT * (aa21 * w1 * w1 * TLt_n1 + aa22 * ddu);
					local = uu0 - DT * (a21 * du) - DT * DT * (aa21 * ddu);
					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					TL_ns1 = D_duns1[iCell][jCell][iVar];
					TLt_ns1 = D_dduns1[iCell][jCell][iVar];

					// local = b1 * Tu1 + b2 * uu0 - DT * (v1 * w1 * TL_n1 + v2 * du0 + v3 * du) -
					// 		DT * DT * (vv1 * w1 * w1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu);

					local = uu0 - DT * (v1 * du0 + v2 * du + w1 * TL_n1 + w2 * TL_ns1) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu + ww1 * TLt_n1 + ww2 * TLt_ns1);

					D_un[iCell][jCell][iVar] = local;

					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
					D_duns1[iCell][jCell][iVar] = du;
					D_dduns1[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TDMSRK326_ALW(unsigned short iRK_Step)
{

	mydouble local, b1, b2, aa21;
	mydouble du, ddu, du0, ddu0, uu, uu0, Tu1 = 0., Tu2 = 0., TL_n1 = 0., TL_n2 = 0., TLt_n1 = 0., TLt_n2 = 0.;
	mydouble a21, aa22, a22, v1, v2, v3, vv1, vv2, vv3, d31, d32, a21_d31, du1, ddu1;
	mydouble b3, d41, d42, a41, a42, a43, aa41, aa42, aa43, v4, vv4, a32, a31, aa31, aa32;

	mydouble arr[24] = {3.4192786601275000E-01, 9.2331503301445700E-02, 6.5807213398725000E-01, 9.0766849669855400E-01,
						1.0426387973076500E-01, 8.9573612026923500E-01, 2.5790165593231500E-01, 4.9998434064656400E-01,
						6.6076932775015300E-02, 7.2129498199869600E-01, 0.0000000000000000E+00, 7.5009584555929800E-02,
						6.9357361470979600E-01, 4.4055647477858000E-33, 3.3568068046503900E-01, 5.9967923477784500E-02,
						1.1348038264019000E-01, 1.2601479587349600E-02, 1.2645394691931600E-01, 1.2239634140739500E-01,
						1.3819314026467700E-02, 1.3801717989153900E-01, 5.8941699256348000E-02, 7.8787813242481900E-02}; // 9.2651413378796700E-01
	d31 = arr[0], d41 = arr[1], d32 = arr[2], d42 = arr[3], b1 = arr[4], b2 = arr[5];
	a31 = arr[6], a32 = arr[7], a41 = arr[8], a42 = arr[9], a43 = arr[10];
	v1 = arr[11], v2 = arr[12], v3 = arr[13], v4 = arr[14];
	aa31 = arr[15], aa32 = arr[16], aa41 = arr[17], aa42 = arr[18], aa43 = arr[19];
	vv1 = arr[20], vv2 = arr[21], vv3 = arr[22], vv4 = arr[23];

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = d31 * Tu1 + d32 * uu0 - DT * (a31 * TL_n1 + a32 * du) -
							DT * DT * (aa31 * TLt_n1 + aa32 * ddu);

					D_un[iCell][jCell][iVar] = local;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = d41 * Tu1 + d42 * uu0 - DT * (a41 * TL_n1 + a42 * du0 + a43 * du) -
							DT * DT * (aa41 * TLt_n1 + aa42 * ddu0 + aa43 * ddu);

					D_un[iCell][jCell][iVar] = local;

					D_du1[iCell][jCell][iVar] = du;
					D_ddu1[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					du1 = D_du1[iCell][jCell][iVar];
					ddu1 = D_ddu1[iCell][jCell][iVar];

					Tu1 = D_un1[iCell][jCell][iVar];
					TL_n1 = D_dun1[iCell][jCell][iVar];
					TLt_n1 = D_ddun1[iCell][jCell][iVar];

					local = b1 * Tu1 + b2 * uu0 - DT * (v1 * TL_n1 + v2 * du0 + v3 * du1 + v4 * du) -
							DT * DT * (vv1 * TLt_n1 + vv2 * ddu0 + vv3 * ddu1 + vv4 * ddu);

					D_un[iCell][jCell][iVar] = local;

					D_un1[iCell][jCell][iVar] = uu0;
					D_dun1[iCell][jCell][iVar] = du0;
					D_ddun1[iCell][jCell][iVar] = ddu0;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TD23_ALW(unsigned short iRK_Step)
{

	// // #pragma omp barrier

	// #pragma omp parallel for private(a21, v2, v1, vv1, vv2, aa21)
	mydouble du, ddu, du0, ddu0, uu, uu0, local;
	mydouble a21, v2, v1, vv1, vv2, aa21;
	a21 = 0.594223212099088, v2 = 0.306027487008159; // SD
	// a21 = 2. / 3., v2 = 3. / 8.;   //Taylor
	aa21 = a21 * a21 * 0.5;
	v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;
	// #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (short iVar = 0; iVar < nVar; iVar++)
			{
				if (iRK_Step == 0)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					uu = uu0 - a21 * DT * du - aa21 * DT * DT * ddu;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
				else
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];
					uu0 = D_un0[iCell][jCell][iVar];

					local = uu0 - DT * (v1 * du0 + v2 * du) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu);
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}

	// // #pragma omp barrier
}

void CADSolver::TimeIntegration_TDTDRK_initial(unsigned short iRK_Step, double a_int)
{

	// // #pragma omp barrier

	// #pragma omp parallel for private(a21, v2, v1, vv1, vv2, aa21)
	mydouble du, ddu, du0, ddu0, uu, uu0, local;
	mydouble a21, v2, v1, vv1, vv2, aa21;
	a21 = a_int;
	aa21 = a21 * a21 * 0.5;
	// #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (short iVar = 0; iVar < nVar; iVar++)
			{
				if (iRK_Step == 0)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];

					uu = uu0 - a21 * DT * du - aa21 * DT * DT * ddu;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
				else
				{
					uu0 = D_un0[iCell][jCell][iVar];

					D_un[iCell][jCell][iVar] = uu0;

					D_dun1[iCell][jCell][iVar] = D_du0[iCell][jCell][iVar];
					D_ddun1[iCell][jCell][iVar] = D_ddu0[iCell][jCell][iVar];

					D_duns1[iCell][jCell][iVar] = D_du[iCell][jCell][iVar];
					D_dduns1[iCell][jCell][iVar] = D_ddu[iCell][jCell][iVar];
				}
			}
		}
	}

	// // #pragma omp barrier
}

void CADSolver::TimeIntegration_TD24_ALW(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2;

	a21 = 0.5, v1 = 1., v2 = 0., vv1 = 1. / 6., vv2 = 1. / 3.; // C0 = 1. / 3., C1 = 2. / 3.;
	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];
					// ddu = 0.;
					uu = uu0 - a21 * DT * du - a21 * a21 * 0.5 * DT * DT * ddu;

					// uu = uu0 - DT * du - 0.5 * DT * DT * ddu;
					// uu = uu0- DT * du;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					local = uu0 - DT * (v1 * du0 + v2 * du) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TD34_ALW(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, du1, ddu1, du2, ddu2;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2, a31, a32, aa31, aa32, v3, vv3, aa21;

	a21 = 0.532230958301739;
	a31 = 0.457385979790557;
	a32 = 0.207902718086617;
	aa31 = 0.0553261314403881;
	v3 = 0.374215233358609;

	aa21 = pow(a21, 2) / 2.;
	aa32 = (pow(a31, 2) - 2. * a21 * a32 + 2. * a31 * a32 + pow(a32, 2) - 2. * aa31) / 2.;
	uu = a31 + a32;
	v2 = -(((3. * pow(a21, 2) * a32 + 6. * a21 * aa31 - 3. * a21 * pow(uu, 2) + pow(uu, 3)) * v3) / pow(a21, 3));
	v1 = 1. - v2 - v3;
	vv1 = -(-1. + 2. * uu - 2. * a21 * (-1. + 3. * uu + a21 * (a21 - 3. * uu) * v2) + 2. * (3. * a21 - uu) * pow(uu, 2) * v3) / (12. * a21 * uu);
	vv2 = -(-1. + 2. * uu + 4. * pow(a21, 3) * v2 - 6. * pow(a21, 2) * uu * v2 - 2. * pow(uu, 3) * v3) / (12. * a21 * (a21 - uu));
	vv3 = -(-1. + 2. * a21 - 2. * pow(a21, 3) * v2 - 6. * a21 * pow(uu, 2) * v3 + 4. * pow(uu, 3) * v3) / (12. * uu * (-a21 + uu));

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];
					// ddu = 0.;
					uu = uu0 - a21 * DT * du - DT * DT * aa21 * ddu;

					// uu = uu0 - DT * du - 0.5 * DT * DT * ddu;
					// uu = uu0- DT * du;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					// local = uu0 - DT * (v1 * du0 + v2 * du) -
					// 		DT * DT * (vv1 * ddu0 + vv2 * ddu);
					local = uu0 - DT * (a31 * du0 + a32 * du) - DT * DT * (aa31 * ddu0 + aa32 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
					D_du1[iCell][jCell][iVar] = du;
					D_ddu1[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					du1 = D_du1[iCell][jCell][iVar];
					ddu1 = D_ddu1[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					local = uu0 - DT * (v1 * du0 + v2 * du1 + v3 * du) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu1 + vv3 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TD35_ALW(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, du1, ddu1, du2, ddu2;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2, a31, a32, aa31, aa32, v3, vv3, aa21, aa12;

	// a21 = 0.532230958301739;
	// a31 = 0.457385979790557;
	// a32 = 0.207902718086617;
	// aa31 = 0.0553261314403881;
	// v3 = 0.374215233358609;

	// aa21 = pow(a21, 2) / 2.;
	// aa32 = (pow(a31, 2) - 2. * a21 * a32 + 2. * a31 * a32 + pow(a32, 2) - 2. * aa31) / 2.;
	// uu = a31 + a32;
	// v2 = -(((3. * pow(a21, 2) * a32 + 6. * a21 * aa31 - 3. * a21 * pow(uu, 2) + pow(uu, 3)) * v3) / pow(a21, 3));
	// v1 = 1. - v2 - v3;
	// vv1 = -(-1. + 2. * uu - 2. * a21 * (-1. + 3. * uu + a21 * (a21 - 3. * uu) * v2) + 2. * (3. * a21 - uu) * pow(uu, 2) * v3) / (12. * a21 * uu);
	// vv2 = -(-1. + 2. * uu + 4. * pow(a21, 3) * v2 - 6. * pow(a21, 2) * uu * v2 - 2. * pow(uu, 3) * v3) / (12. * a21 * (a21 - uu));
	// vv3 = -(-1. + 2. * a21 - 2. * pow(a21, 3) * v2 - 6. * a21 * pow(uu, 2) * v3 + 4. * pow(uu, 3) * v3) / (12. * uu * (-a21 + uu));

	// a21 = 0.1; //
	a21 = 0.752;
	aa21 = a21 * a21 / 2.;
	aa12 = -1. + 2. * a21;
	vv1 = (1. + 2. * a21 * (-4. + 5. * a21)) / (12. * a21 * (-3. + 5. * a21));
	vv2 = 1. / (12. * a21 * (3. + 10. * (-1. + a21) * a21));
	vv3 = (25. * pow(aa12, 3)) / (12. * (-3. + 5. * a21) * (3. + 10. * (-1 + a21) * a21));
	v1 = 1., v2 = 0., v3 = 0.;
	a31 = (3. - 5. * a21) / (5. - 10. * a21), a32 = 0.;
	aa31 = ((-3. + 5. * a21) * (-3. + 10. * a21) * (1. + 5. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));
	aa32 = ((-3. + 5. * a21) * (3. + 10. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];
					// ddu = 0.;
					uu = uu0 - a21 * DT * du - DT * DT * aa21 * ddu;

					// uu = uu0 - DT * du - 0.5 * DT * DT * ddu;
					// uu = uu0- DT * du;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					// local = uu0 - DT * (v1 * du0 + v2 * du) -
					// 		DT * DT * (vv1 * ddu0 + vv2 * ddu);
					local = uu0 - DT * (a31 * du0 + a32 * du) - DT * DT * (aa31 * ddu0 + aa32 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
					D_du1[iCell][jCell][iVar] = du;
					D_ddu1[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					du1 = D_du1[iCell][jCell][iVar];
					ddu1 = D_ddu1[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					local = uu0 - DT * (v1 * du0 + v2 * du1 + v3 * du) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu1 + vv3 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_TD46_ALW(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, du1, ddu1, du2, ddu2, du3, ddu3;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2, a31, a32, aa31, aa32, v3, vv3, aa21, aa12;
	mydouble a41, aa42, aa43, vv4, a42, a43, aa41;

	a21 = 0.235740260201908, aa41 = 0.0337154180289742;
	a31 = 0.766058541677957, aa42 = 0.0725363672383309;
	a41 = 0.504428678818826, aa43 = 0.0209723607400122;
	v1 = 1., vv1 = 0.0740479031154873;
	aa21 = 0.0277867351400587, vv2 = 0.253930692818331;
	aa31 = 0.00713461664453855, vv3 = 0.0765679011859131;
	aa32 = 0.286288227994331, vv4 = 0.0954535028801300;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];
					uu0 = D_un[iCell][jCell][iVar];
					// ddu = 0.;
					uu = uu0 - a21 * DT * du - DT * DT * aa21 * ddu;

					// uu = uu0 - DT * du - 0.5 * DT * DT * ddu;
					// uu = uu0- DT * du;

					D_un[iCell][jCell][iVar] = uu;
					D_un0[iCell][jCell][iVar] = uu0;
					D_du0[iCell][jCell][iVar] = du;
					D_ddu0[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					// local = uu0 - DT * (v1 * du0 + v2 * du) -
					// 		DT * DT * (vv1 * ddu0 + vv2 * ddu);
					local = uu0 - DT * (a31 * du0) - DT * DT * (aa31 * ddu0 + aa32 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
					D_du1[iCell][jCell][iVar] = du;
					D_ddu1[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					du1 = D_du1[iCell][jCell][iVar];
					ddu1 = D_ddu1[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					// local = uu0 - DT * (v1 * du0 + v2 * du) -
					// 		DT * DT * (vv1 * ddu0 + vv2 * ddu);
					local = uu0 - DT * (a41 * du0) - DT * DT * (aa41 * ddu0 + aa42 * ddu1 + aa43 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
					D_dun2[iCell][jCell][iVar] = du;
					D_ddun2[iCell][jCell][iVar] = ddu;
				}
			}
		}
	}
	else if (iRK_Step == 3)
	{
		// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					ddu = D_ddu[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					ddu0 = D_ddu0[iCell][jCell][iVar];

					du1 = D_du1[iCell][jCell][iVar];
					ddu1 = D_ddu1[iCell][jCell][iVar];

					du2 = D_dun2[iCell][jCell][iVar];
					ddu2 = D_ddun2[iCell][jCell][iVar];

					uu0 = D_un0[iCell][jCell][iVar];

					local = uu0 - DT * (v1 * du0) -
							DT * DT * (vv1 * ddu0 + vv2 * ddu1 + vv3 * ddu2 + vv4 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_RK22(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0, un0, un1;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2;
	a21 = 0.594223212099088, v2 = 0.306027487008159;
	ww1 = 0.0, ww2 = 0.0, w11 = 0.0, w22 = 0.0;
	v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

	if (iRK_Step == 0)
	{
		// // #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un[iCell][jCell][iVar];

					local = un0 - DT * du;

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = un0;
					D_du0[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// // #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{

			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];
					un1 = D_un[iCell][jCell][iVar];

					// local = Element[iCell][jCell]->GetSolution_Old(0, iVar) -
					// 		DT * (v1 * du0 + v2 * du) -
					// 		DT * DT * (vv1 * ddu0 + vv2 * ddu);
					// local = 0.75 * Element[iCell][jCell]->GetSolution_Old(0, iVar) +
					// 		1.0 / 4.0 * (Element[iCell][jCell]->GetSolution(0, 0, iVar) - DT * Element[iCell][jCell]->GetSolution_du(0, iVar));

					local = 0.5 * (un0 + un1) - 0.5 * DT * du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_RK54(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, du0, un0, un1, un2, un3, un4, du1;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un[iCell][jCell][iVar];

					local = un0 - 0.391752226571890 * DT * du;

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = un0;
					D_du0[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					un0 = D_un0[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un1 = D_un[iCell][jCell][iVar];
					du = D_du[iCell][jCell][iVar];

					local = 0.444370493651235 * un0 + 0.555629506348765 * un1 - 0.368410593050371 * DT * du;

					D_un1[iCell][jCell][iVar] = un1;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];
					un1 = D_un1[iCell][jCell][iVar];
					un2 = D_un[iCell][jCell][iVar];

					local = 0.620101851488403 * un0 + 0.379898148511597 * un2 - 0.251891774271694 * DT * du;

					D_un2[iCell][jCell][iVar] = un2;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 3)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];
					un3 = D_un[iCell][jCell][iVar];

					local = 0.178079954393132 * un0 + 0.821920045606868 * un3 - 0.544974750228521 * DT * du;
					D_un3[iCell][jCell][iVar] = un3;
					D_du1[iCell][jCell][iVar] = du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 4)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];

					un4 = D_un[iCell][jCell][iVar];
					un3 = D_un3[iCell][jCell][iVar];
					un2 = D_un2[iCell][jCell][iVar];

					du1 = D_du1[iCell][jCell][iVar];

					local = 0.517231671970585 * un2 + 0.096059710526147 * un3 - 0.063692468666290 * DT * du1 +
							0.386708617503269 * un4 - 0.226007483236906 * DT * du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_RK65(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, du0, un0, un1, un2, un3, un4, du1, du2, du3, du4;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un[iCell][jCell][iVar];

					local = un0 - 1. / 4. * DT * du;

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = un0;
					D_du0[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					un0 = D_un0[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					du = D_du[iCell][jCell][iVar];

					local = un0 - 1. / 8. * DT * (du + du0);

					D_du1[iCell][jCell][iVar] = du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];

					local = un0 - 0.5 * DT * du;

					D_du2[iCell][jCell][iVar] = du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 3)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					du1 = D_du1[iCell][jCell][iVar];
					du2 = D_du2[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];

					local = un0 - DT * (3. / 16. * du0 - 3. / 8. * du1 + 3. / 8. * du2 + 9. / 16. * du);

					// u_n4[nb + j] = u_n0[nb + j] + 3. / 16. * dt * L_n0[nb + j] - 3. / 8. * dt * L_n1[nb + j] +
					// 			   3. / 8. * dt * L_n2[nb + j] + 9. / 16. * dt * L_n3[nb + j];

					D_du3[iCell][jCell][iVar] = du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 4)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					du1 = D_du1[iCell][jCell][iVar];
					du2 = D_du2[iCell][jCell][iVar];
					du3 = D_du3[iCell][jCell][iVar];

					local = un0 - DT * (-3. / 7. * du0 + 8. / 7. * du1 + 6. / 7. * du2 - 12. / 7. * du3 + 8. / 7. * du);

					// u_n5[nb + j] = u_n0[nb + j] - 3. / 7. * dt * L_n0[nb + j] + 8. / 7. * dt * L_n1[nb + j] + 6. / 7. * dt * L_n2[nb + j] -
					//                12. / 7. * dt * L_n3[nb + j] + 8. / 7. * dt * L_n4[nb + j];

					D_un[iCell][jCell][iVar] = local;
					D_du4[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 5)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					du1 = D_du1[iCell][jCell][iVar];
					du2 = D_du2[iCell][jCell][iVar];
					du3 = D_du3[iCell][jCell][iVar];
					du4 = D_du4[iCell][jCell][iVar];

					// u_nn[nb + j] = u_n0[nb + j] + 7. / 90. * dt * L_n0[nb + j] + 0. * dt * L_n1[nb + j] +
					//                16. / 45. * dt * L_n2[nb + j] + 2. / 15. * dt * L_n3[nb + j] +
					//                16. / 45. * dt * L_n4[nb + j] + 7. / 90. * dt * df[nb + j];
					// 将上式换为 local
					local = un0 - DT * (7. / 90. * du0 + 0. * du1 + 16. / 45. * du2 + 2. / 15. * du3 + 16. / 45. * du4 + 7. / 90. * du);

					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::TimeIntegration_RK76(unsigned short iRK_Step)
{

	mydouble local;
	mydouble du, du0, un0, un1, un2, un3, un4, du1, du2, du3, du4, du5;

	if (iRK_Step == 0)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un[iCell][jCell][iVar];

					local = un0 - 1. / 6. * DT * du;

					D_un[iCell][jCell][iVar] = local;

					D_un0[iCell][jCell][iVar] = un0;
					D_du0[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 1)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{

				for (short iVar = 0; iVar < nVar; iVar++)
				{
					un0 = D_un0[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					du = D_du[iCell][jCell][iVar];

					local = un0 - DT * (12. / 169. * du0 + 27. / 169. * du);

					D_du1[iCell][jCell][iVar] = du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 2)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					du1 = D_du1[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];

					local = un0 - DT * (107952. / 571787. * du0 - 406107. / 571787. * du1 + 566826. / 571787. * du);

					D_du2[iCell][jCell][iVar] = du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 3)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{
					du = D_du[iCell][jCell][iVar];
					du0 = D_du0[iCell][jCell][iVar];
					du1 = D_du1[iCell][jCell][iVar];
					du2 = D_du2[iCell][jCell][iVar];

					un0 = D_un0[iCell][jCell][iVar];

					local = un0 - DT * (1371923. / 6669000. * du0 - 819. / 3800. * du1 + 19411847. / 88236000. * du2 + 561846173. / 1147068000. * du);

					D_du3[iCell][jCell][iVar] = du;
					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
	else if (iRK_Step == 4)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					du1 = D_du1[iCell][jCell][iVar];
					du2 = D_du2[iCell][jCell][iVar];
					du3 = D_du3[iCell][jCell][iVar];

					local = un0 - DT * (-1563412. / 5835375. * du0 + 468. / 475. * du1 - 14488201. / 168199875. * du2 -
										1711096709. / 6846562125. * du3 + 648832. / 1549583. * du);

					D_un[iCell][jCell][iVar] = local;
					D_du4[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 5)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					du1 = D_du1[iCell][jCell][iVar];
					du2 = D_du2[iCell][jCell][iVar];
					du3 = D_du3[iCell][jCell][iVar];
					du4 = D_du4[iCell][jCell][iVar];

					local = un0 - DT * (120115. / 277641. * du0 - 117. / 113. * du1 + 219237109. / 296102601. * du2 +
										29855183083. / 44628054003. * du3 - 3009600. / 9215941. * du4 + 297825. / 572797. * du);

					D_un[iCell][jCell][iVar] = local;
					D_un3[iCell][jCell][iVar] = du;
				}
			}
		}
	}
	else if (iRK_Step == 6)
	{
		// #pragma omp parallel for num_threads(NumberThreads)
		for (long jCell = nb; jCell < nCellYo_L; jCell++)
		{
			mydouble local;
			for (long iCell = nb; iCell < nCell_L; iCell++)
			{
				for (short iVar = 0; iVar < nVar; iVar++)
				{

					du = D_du[iCell][jCell][iVar];
					un0 = D_un0[iCell][jCell][iVar];

					du0 = D_du0[iCell][jCell][iVar];
					du1 = D_du1[iCell][jCell][iVar];
					du2 = D_du2[iCell][jCell][iVar];
					du3 = D_du3[iCell][jCell][iVar];
					du4 = D_du4[iCell][jCell][iVar];
					du5 = D_un3[iCell][jCell][iVar];

					local = un0 - DT * (137. / 1872. * du0 + 371293. / 1164612. * du2 + 3939040643. / 23169727152. * du3 +
										19000. / 104859. * du4 + 45125. / 243312. * du5 + 113. / 1584. * du);

					D_un[iCell][jCell][iVar] = local;
				}
			}
		}
	}
}

void CADSolver::store_du_duu()
{

	mydouble local;
	mydouble du, ddu, du0, ddu0, uu, uu0;
	mydouble a21, v2, ww1, ww2, w11, w22, v1, vv1, vv2;

	// a21 = 0.5, v1 = 1., v2 = 0., vv1 = 1. / 6., vv2 = 1. / 3.; // C0 = 1. / 3., C1 = 2. / 3.;

	// #pragma omp parallel for private(du, ddu, du0, ddu0, uu, uu0, a21, v2, v1, vv1, vv2, local)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			for (short iVar = 0; iVar < nVar; iVar++)
			{

				D_un2[iCell][jCell][iVar] = D_un1[iCell][jCell][iVar];
				D_dun2[iCell][jCell][iVar] = D_dun1[iCell][jCell][iVar];
				D_ddun2[iCell][jCell][iVar] = D_ddun1[iCell][jCell][iVar];

				D_un1[iCell][jCell][iVar] = D_un0[iCell][jCell][iVar];
				D_dun1[iCell][jCell][iVar] = D_du0[iCell][jCell][iVar];
				D_ddun1[iCell][jCell][iVar] = D_ddu0[iCell][jCell][iVar];

				D_duns1[iCell][jCell][iVar] = D_du[iCell][jCell][iVar];
				D_dduns1[iCell][jCell][iVar] = D_ddu[iCell][jCell][iVar];
			}
		}
	}
}
// ___________________________
void CADSolver::Boundary_uu(void)
{
	// #pragma omp barrier

	switch (BoundaryScheme)
	{
	case 0:
		Boundary_uu_periodic();
		break;
	case 1:
		Boundary_uu_Expand();
		break;

	case 11:
		Boundary_uu_Expand_Shock_vortex();
		break;

	case 2:
		Boundary_uu_reflection_x();
		break;
	case 24:
		Boundary_uu_Viscous_shock();
		break;

	case 22:
		Boundary_uu_DoubleMach();
		break;

	case 23:
		Boundary_uu_High_jet();
		break;
	case 26:
		Boundary_uu_Couette_flow();
		break;
	case 27:
		Boundary_uu_Laminar();
		break;

	default:
		cout << "Boundary_uu is not avaliable for this Solver." << endl;
		break;
	}

	// #pragma omp barrier
}

void CADSolver::Boundary_uu_periodic(void)
{
	mydouble uuu;
	// #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nCell + iCell][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell - nCell][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}

	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++) // y方向
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][jCell + nCellYo][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][jCell - nCellYo][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_uu_Expand(void)
{
	mydouble uuu;
	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nb][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nCell_L - 1][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}

	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nb][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_uu_Expand_Shock_vortex(void)
{
	mydouble uuu;
	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nb][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nCell_L - 1][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}

	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nb][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_uu_reflection_x(void)
{
	mydouble uuu;
	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_un[nb * 2 - iCell - 1][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
			uuu = D_un[nb * 2 - iCell - 1][jCell][1];
			D_un[iCell][jCell][1] = -uuu;
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_un[nCell_L * 2 - 1 - iCell][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
			uuu = D_un[nCell_L * 2 - 1 - iCell][jCell][1];
			D_un[iCell][jCell][1] = -uuu;
		}
	}

	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nb][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_uu_Viscous_shock(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu, density, velocity_x, velocity_y, pressure;
		long nb1;
		for (long iCell = 0; iCell < nb; iCell++) //_________lift__________
		{
			nb1 = nb * 2 - iCell - 1;
			density = D_un[nb1][jCell][0], velocity_x = D_un[nb1][jCell][1] / density;
			velocity_y = D_un[nb1][jCell][2] / density;
			pressure = (Gamma - 1.0) * (D_un[nb1][jCell][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));

			D_un[iCell][jCell][0] = density;
			D_un[iCell][jCell][1] = 0.0;
			D_un[iCell][jCell][2] = 0.0;
			D_un[iCell][jCell][3] = pressure / Gamma_1;
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++) //_________right__________
		{

			nb1 = nCell_L * 2 - 1 - iCell;
			density = D_un[nb1][jCell][0], velocity_x = D_un[nb1][jCell][1] / density;
			velocity_y = D_un[nb1][jCell][2] / density;
			pressure = (Gamma - 1.0) * (D_un[nb1][jCell][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));
			D_un[iCell][jCell][0] = density;
			D_un[iCell][jCell][1] = 0.0;
			D_un[iCell][jCell][2] = 0.0;
			D_un[iCell][jCell][3] = pressure / Gamma_1;
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu, density, velocity_x, velocity_y, pressure;
		long nb1;
		for (long jCell = 0; jCell < nb; jCell++) //_________bottom__________
		{
			nb1 = nb * 2 - jCell - 1;
			density = D_un[iCell][nb1][0], velocity_x = D_un[iCell][nb1][1] / density;
			velocity_y = D_un[iCell][nb1][2] / density;
			pressure = (Gamma - 1.0) * (D_un[iCell][nb1][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));

			// D_un[iCell][jCell][0] = density;
			// D_un[iCell][jCell][1] = 0.0;
			// D_un[iCell][jCell][2] = 0.0;
			// D_un[iCell][jCell][3] = pressure / Gamma_1;
			D_un[iCell][jCell][0] = density;
			D_un[iCell][jCell][1] = -density * velocity_x;
			D_un[iCell][jCell][2] = 0.0;
			D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * velocity_y * velocity_y;
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++) //_________top__________
		{

			nb1 = nCellYo_L * 2 - 1 - jCell;
			density = D_un[iCell][nb1][0], velocity_x = D_un[iCell][nb1][1] / density;
			velocity_y = D_un[iCell][nb1][2] / density;
			pressure = (Gamma - 1.0) * (D_un[iCell][nb1][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));

			D_un[iCell][jCell][0] = density;
			D_un[iCell][jCell][1] = D_un[iCell][nb1][1];
			D_un[iCell][jCell][2] = 0.0;
			D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * (velocity_y * velocity_y);
		}
	}
}

void CADSolver::Boundary_uu_Couette_flow(void)
{
	// #pragma omp parallel for num_threads(NumberThreads)
	// 	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	// 	{
	// 		mydouble uuu;
	// 		for (long iCell = 0; iCell < nb; iCell++)
	// 			for (long nn = 0; nn < 4; nn++)
	// 			{

	// 				uuu = D_un[nCell + iCell][jCell][nn];
	// 				D_un[iCell][jCell][nn] = uuu;
	// 			}

	// 		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
	// 			for (long nn = 0; nn < 4; nn++)
	// 			{

	// 				uuu = D_un[iCell - nCell][jCell][nn];
	// 				D_un[iCell][jCell][nn] = uuu;
	// 			}
	// 	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	{
		mydouble pp, dd, uu, vv, yp;
		mydouble H, u1, mu, ka, T0, T1, Ta;
		for (long iCell = 0; iCell < nb; iCell++) // 入口
		{
			yp = D_Coord[iCell][jCell][1];
			H = 2., u1 = 0.1;
			uu = yp * u1 / H, vv = 0.;
			T0 = 0.8, T1 = 0.85;
			Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);
			pp = 1. / Gamma, dd = Gamma * pp / Ta;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
		for (long iCell = nCell_L; iCell < nCell_2L; iCell++) // 出口
		{
			yp = D_Coord[iCell][jCell][1];
			H = 2., u1 = 0.1;
			uu = yp * u1 / H, vv = 0.;
			T0 = 0.8, T1 = 0.85;
			Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);
			pp = 1. / Gamma, dd = Gamma * pp / Ta;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble pp, dd, uu, vv, yp;
		mydouble H, u1, mu, ka, T0, T1, Ta;
		for (long jCell = 0; jCell < nb + 1; jCell++) //_________bottom__________
		{

			H = 2., u1 = 0.1;
			uu = 0., vv = 0.;
			T0 = 0.8, T1 = 0.85;
			// Ta = 0.8;
			// T0 + yp *(T1 - T0) / H + Pr *Gamma_1 *u1 *u1 *yp / (2. * H) * (1 - yp / H);

			pp = 1. / Gamma, dd = Gamma * pp / T0;

			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++) //_________top__________
		{
			yp = D_Coord[iCell][nCellYo_L][1];
			H = 2., u1 = 0.1;
			uu = yp * u1 / H, vv = 0.;
			T0 = 0.8, T1 = 0.85;
			Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);
			pp = 1. / Gamma, dd = Gamma * pp / Ta;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		for (long jCell = nb; jCell < nb + 4; jCell++) //_________1 - 3__________
		{
			yp = D_Coord[iCell][jCell][1];
			H = 2., u1 = 0.1;
			uu = yp * u1 / H, vv = 0.;
			T0 = 0.8, T1 = 0.85;
			Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);
			pp = 1. / Gamma, dd = Gamma * pp / Ta;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
		for (long jCell = nCellYo_L - 3; jCell < nCellYo_L; jCell++) //_________L-3  L __________
		{
			yp = D_Coord[iCell][jCell][1];
			H = 2., u1 = 0.1;
			uu = yp * u1 / H, vv = 0.;
			T0 = 0.8, T1 = 0.85;
			Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);
			pp = 1. / Gamma, dd = Gamma * pp / Ta;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
	}
}

void CADSolver::Boundary_uu_DoubleMach(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble alfa, dd, uu, vv, pp;
		alfa = 60. / 180. * PI_Number;
		dd = 8.0, uu = 8.25 * sin(alfa), vv = -8.25 * cos(alfa), pp = 116.5;

		for (long iCell = 0; iCell < nb; iCell++) //_________lift__inlet________
		{

			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++) //_________right_out_________
		{

			for (long nn = 0; nn < 4; nn++)
			{
				D_un[iCell][jCell][nn] = D_un[nCell_L - 1][jCell][nn];
			}
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble alfa, dd, uu, vv, pp, x0, x1, xp, yp, density, velocity_x, velocity_y, pressure;
		long nb1;

		alfa = 60. / 180. * PI_Number;
		dd = 8.0, uu = 8.25 * sin(alfa), vv = -8.25 * cos(alfa), pp = 116.5;
		x0 = 1. / 6.;

		for (long jCell = 0; jCell < nb; jCell++) //_________bottom__________
		{
			xp = D_Coord[iCell][jCell][0];
			if (xp < x0)
			{
				D_un[iCell][jCell][0] = dd;
				D_un[iCell][jCell][1] = dd * uu;
				D_un[iCell][jCell][2] = dd * vv;
				D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
			}
			else
			{

				nb1 = nb * 2 - jCell - 1;
				density = D_un[iCell][nb1][0], velocity_x = D_un[iCell][nb1][1] / density;
				velocity_y = D_un[iCell][nb1][2] / density;

				D_un[iCell][jCell][0] = density;
				D_un[iCell][jCell][1] = density * velocity_x;
				D_un[iCell][jCell][2] = 0.0;
				D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * velocity_y * velocity_y;
			}
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++) //_________top__________
		{
			yp = D_Coord[iCell][jCell][1];
			xp = D_Coord[iCell][jCell][0];
			x1 = x0 + yp / tan(alfa) + 10. * TimeNow / sin(alfa);

			if (xp < x1)
			{
				dd = 8.0, uu = 8.25 * sin(alfa), vv = -8.25 * cos(alfa), pp = 116.5;
			}
			else
			{
				dd = 1.4, uu = 0.0, vv = 0.0, pp = 1.0;
			}

			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}
	}
}

void CADSolver::Boundary_uu_High_jet(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu, yp, uu, vv, pp, dd;
		for (long iCell = 0; iCell < nb; iCell++)
		{
			yp = D_Coord[iCell][jCell][1];
			if (abs(yp) <= 0.050001)
			{
				uu = 30.0, vv = 0., pp = 0.4127, dd = 5.;
			}
			else
			{
				uu = 0.0, vv = 0., pp = 0.4127, dd = 5.;
			}

			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[nCell_L - 1][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu;
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nb][nn];
				D_un[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_uu_Laminar0(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu, yp, uu, vv, pp, dd;
		for (long iCell = 0; iCell < nb + 6; iCell++) //----------------Inlet boundary -----
		{

			uu = Ma, vv = 0., pp = 1., dd = 1.;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_un[nCell_L - 1][jCell][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu, density, velocity_x, velocity_y, pressure, xx;
		mydouble yp, uu, vv, pp, dd;
		long nb1;
		for (long jCell = 0; jCell < nb; jCell++) //_________bottom__________
		{
			xx = D_Coord[iCell][jCell][0];
			nb1 = nb * 2 - jCell - 1;
			density = D_un[iCell][nb1][0], velocity_x = D_un[iCell][nb1][1] / density;
			velocity_y = D_un[iCell][nb1][2] / density;
			pressure = (Gamma - 1.0) * (D_un[iCell][nb1][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));

			if (xx >= (nx_end - nx_begin) * 0.5)
			{
				// D_un[iCell][jCell][0] = density;
				// D_un[iCell][jCell][1] = 0.0;
				// D_un[iCell][jCell][2] = 0.0;
				// D_un[iCell][jCell][3] = pressure / Gamma_1;

				D_un[iCell][jCell][0] = density;
				D_un[iCell][jCell][1] = -density * velocity_x;
				D_un[iCell][jCell][2] = 0.0;
				D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * velocity_y * velocity_y;
			}
			else
			{
				D_un[iCell][jCell][0] = density;
				D_un[iCell][jCell][1] = density * velocity_x;
				D_un[iCell][jCell][2] = 0.0;
				D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * velocity_y * velocity_y;
			}
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++) //_________top__________
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
		// {
		// 	uu = 0.15, vv = 0., pp = 1., dd = 1.;
		// 	D_un[iCell][jCell][0] = dd;
		// 	D_un[iCell][jCell][1] = dd * uu;
		// 	D_un[iCell][jCell][2] = dd * vv;
		// 	D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		// }
	}
}

void CADSolver::Boundary_uu_Laminar(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu, yp, uu, vv, pp, dd, u, v, p, T, c, ir, d;
		mydouble d1, d2, d3, u1, u2, u3, v1, v2, v3, p1, p2, p3, T1, T2, c1, c2, ir2;
		mydouble deltav[5], deltaw[5];
		mydouble LL[5][5], RR[5][5];
		mydouble d1k, d2k, d3k, d1kh, d2kh, d3kh;
		mydouble rhoh, rhohc, rhohci;

		long nb1;
		for (long iCell = 0; iCell < nb + 6; iCell++) //----------------Inlet boundary -----
		{

			uu = Ma, vv = 0., pp = 1., dd = 1.;
			D_un[iCell][jCell][0] = dd;
			D_un[iCell][jCell][1] = dd * uu;
			D_un[iCell][jCell][2] = dd * vv;
			D_un[iCell][jCell][3] = pp / Gamma_1 + 0.5 * dd * (uu * uu + vv * vv);
		}

		// for (long iCell = nCell_L; iCell < nCell_2L; iCell++) //----------------Outlet boundary -----
		{
			nb1 = nCell_L - 1;
			d2 = D_un[nb1][jCell][0], u2 = D_un[nb1][jCell][1] / d2;
			v2 = D_un[nb1][jCell][2] / d2;
			p2 = (Gamma - 1.0) * (D_un[nb1][jCell][3] - 0.5 * d2 * (u2 * u2 + v2 * v2));
			T2 = p2 * Gamma / d2;
			c2 = sqrt(Gamma * p2 / d2);
			// nb1 = nCell_L - 2;
			// d1 = D_un[nb1][jCell][0], u1 = D_un[nb1][jCell][1] / d1;
			// v1 = D_un[nb1][jCell][2] / d1;
			// p1 = (Gamma - 1.0) * (D_un[nb1][jCell][3] - 0.5 * d1 * (u1 * u1 + v1 * v1));
			// T1 = p1 * Gamma / d1;
			// c1 = sqrt(Gamma * p1 / d1);

			d = d2;
			u = u2;
			v = v2;

			p = 2. - p2;

			for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			{
				D_un[iCell][jCell][0] = d;
				D_un[iCell][jCell][1] = d * u;
				D_un[iCell][jCell][2] = d * v;
				D_un[iCell][jCell][3] = p / (Gamma_1) + 0.5 * d * (u * u + v * v);
			}
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu, density, velocity_x, velocity_y, pressure, xx;
		mydouble yp, uu, vv, pp, dd;
		long nb1;
		for (long jCell = 0; jCell < nb; jCell++) //_________bottom__________
		{
			xx = D_Coord[iCell][jCell][0];
			nb1 = nb * 2 - jCell - 1;
			density = D_un[iCell][nb1][0], velocity_x = D_un[iCell][nb1][1] / density;
			velocity_y = D_un[iCell][nb1][2] / density;
			pressure = (Gamma - 1.0) * (D_un[iCell][nb1][3] - 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y));

			// if (xx >= (nx_end - nx_begin) * 0.5)
			if (xx >= 0.0)
			{
				// D_un[iCell][jCell][0] = density;
				// D_un[iCell][jCell][1] = 0.0;
				// D_un[iCell][jCell][2] = 0.0;
				// D_un[iCell][jCell][3] = pressure / Gamma_1;

				D_un[iCell][jCell][0] = density;
				D_un[iCell][jCell][1] = -density * velocity_x;
				D_un[iCell][jCell][2] = 0.0;
				D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * velocity_y * velocity_y;
			}
			else
			{
				D_un[iCell][jCell][0] = density;
				D_un[iCell][jCell][1] = density * velocity_x;
				D_un[iCell][jCell][2] = 0.0;
				D_un[iCell][jCell][3] = D_un[iCell][nb1][3] - 0.5 * density * velocity_y * velocity_y;
			}
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++) //_________top__________
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_un[iCell][nCellYo_L - 1][nn];
				D_un[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_dutxy_all(void)
{
	switch (BoundaryScheme)
	{
	case 0:
		Boundary_dutxy_periodic();
		break;
	case 1:
		Boundary_dutxy_expand();
		break;

	case 2:
		Boundary_dutxy_reflection_x();
		break;

	default:
		cout << "BoundaryScheme is not avaliable for this Solver." << endl;
		break;
	}
}

void CADSolver::Boundary_dutxy_periodic(void)
{
	mydouble uuu;
	// #pragma omp parallel for num_threads(NumberThreads) // x方向
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nCell + iCell][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[iCell - nCell][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}
	}

	// #pragma omp parallel for num_threads(NumberThreads) // y方向
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][jCell + nCellYo][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][jCell - nCellYo][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_dutxy_expand(void)
{
	mydouble uuu;
	// #pragma omp parallel for num_threads(NumberThreads) // x方向
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nb][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nCell_L - 1][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}
	}

	// #pragma omp parallel for num_threads(NumberThreads) // y方向
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][nb][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][nCellYo_L - 1][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_dutxy_reflection_x(void)
{
	mydouble uuu;
	// #pragma omp parallel for num_threads(NumberThreads) // x方向
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nb * 2 - iCell - 1][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}
			uuu = D_utx0[nb * 2 - iCell - 1][jCell][1];
			D_utx0[iCell][jCell][1] = uuu;
		}
		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_utx0[nCell_L * 2 - 1 - iCell][jCell][nn];
				D_utx0[iCell][jCell][nn] = uuu;
			}
			uuu = D_utx0[nCell_L * 2 - 1 - iCell][jCell][1];
			D_utx0[iCell][jCell][1] = uuu;
		}
	}

	// #pragma omp parallel for num_threads(NumberThreads) // y方向
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][nb][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_uty0[iCell][nCellYo_L - 1][nn];
				D_uty0[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_dutxy_LWA_periodic(void)
{
	mydouble uuu;

	// #pragma omp parallel for num_threads(NumberThreads) // x方向
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_fxat[nCell + iCell][jCell][nn][0];
				D_fxat[iCell][jCell][nn][0] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_fxat[iCell - nCell][jCell][nn][0];
				D_fxat[iCell][jCell][nn][0] = uuu;
			}
	}

	// #pragma omp parallel for num_threads(NumberThreads) // y方向
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_fyat[iCell][jCell + nCellYo][nn][0];
				D_fyat[iCell][jCell][nn][0] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_fyat[iCell][jCell - nCellYo][nn][0];
				D_fyat[iCell][jCell][nn][0] = uuu;
			}
	}
}

// --------------------Boundary_du----------------------
void CADSolver::Boundary_du(void)
{
	// #pragma omp barrier

	switch (BoundaryScheme)
	{
	case 0:
		Boundary_du_periodic();
		break;
	case 1:
		Boundary_du_Expand();
		break;
	case 2:
		Boundary_du_reflection_x();
		break;
		// case 24:
		// 	Boundary_uu_Viscous_shock();
		// 	break;

		// case 22:
		// 	Boundary_uu_DoubleMach();
		// 	break;

		// case 23:
		// 	Boundary_uu_High_jet();
		// 	break;
		// case 26:
		// 	Boundary_uu_Couette_flow();
		// 	break;
		// case 27:
		// 	Boundary_uu_Laminar();
		// 	break;

	default:
		cout << "Boundary_uu is not avaliable for this Solver." << endl;
		break;
	}

	// #pragma omp barrier
}
void CADSolver::Boundary_du_periodic(void)
{

#pragma omp parallel for							  // num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	{
		mydouble uuu;
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_du[nCell + iCell][jCell][nn];
				D_du[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_du[iCell - nCell][jCell][nn];
				D_du[iCell][jCell][nn] = uuu;
			}
	}

#pragma omp parallel for							// num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++) // y方向
	{
		mydouble uuu;
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_du[iCell][jCell + nCellYo][nn];
				D_du[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_du[iCell][jCell - nCellYo][nn];
				D_du[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_du_Expand(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu;
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_du[nb][jCell][nn];
				D_du[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_du[nCell_L - 1][jCell][nn];
				D_du[iCell][jCell][nn] = uuu;
			}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu;
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_du[iCell][nb][nn];
				D_du[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_du[iCell][nCellYo_L - 1][nn];
				D_du[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_du_reflection_x(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu;
		for (long iCell = 0; iCell < nb; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_du[nb * 2 - iCell - 1][jCell][nn];
				D_du[iCell][jCell][nn] = uuu;
			}
			uuu = D_du[nb * 2 - iCell - 1][jCell][1];
			D_du[iCell][jCell][1] = -uuu;
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_ddu[nCell_L * 2 - 1 - iCell][jCell][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}
			uuu = D_du[nCell_L * 2 - 1 - iCell][jCell][1];
			D_du[iCell][jCell][1] = -uuu;
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu;
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_du[iCell][nb][nn];
				D_du[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_du[iCell][nCellYo_L - 1][nn];
				D_du[iCell][jCell][nn] = uuu;
			}
	}
}

// --------------------Boundary_ddu----------------------

void CADSolver::Boundary_ddu(void)
{
	// #pragma omp barrier

	switch (BoundaryScheme)
	{
	case 0:
		Boundary_ddu_periodic();
		break;
	case 1:
		Boundary_ddu_Expand();
		break;
	case 2:
		Boundary_ddu_reflection_x();
		break;
		// case 24:
		// 	Boundary_uu_Viscous_shock();
		// 	break;

		// case 22:
		// 	Boundary_uu_DoubleMach();
		// 	break;

		// case 23:
		// 	Boundary_uu_High_jet();
		// 	break;
		// case 26:
		// 	Boundary_uu_Couette_flow();
		// 	break;
		// case 27:
		// 	Boundary_uu_Laminar();
		// 	break;

	default:
		cout << "Boundary_uu is not avaliable for this Solver." << endl;
		break;
	}

	// #pragma omp barrier
}

void CADSolver::Boundary_ddu_periodic(void)
{

#pragma omp parallel for							  // num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	{
		mydouble uuu;
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_ddu[nCell + iCell][jCell][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_ddu[iCell - nCell][jCell][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}
	}

#pragma omp parallel for							// num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++) // y方向
	{
		mydouble uuu;
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_ddu[iCell][jCell + nCellYo][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_ddu[iCell][jCell - nCellYo][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_ddu_Expand(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu;
		for (long iCell = 0; iCell < nb; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_ddu[nb][jCell][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_ddu[nCell_L - 1][jCell][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu;
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_ddu[iCell][nb][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_ddu[iCell][nCellYo_L - 1][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}
	}
}

void CADSolver::Boundary_ddu_reflection_x(void)
{

#pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		mydouble uuu;
		for (long iCell = 0; iCell < nb; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_ddu[nb * 2 - iCell - 1][jCell][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}
			// uuu = D_un[nb * 2 - iCell - 1][jCell][1];
			// D_un[iCell][jCell][1] = -uuu;
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			for (long nn = 0; nn < 4; nn++)
			{
				uuu = D_ddu[nCell_L * 2 - 1 - iCell][jCell][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}
			// uuu = D_un[nCell_L * 2 - 1 - iCell][jCell][1];
			// D_un[iCell][jCell][1] = -uuu;
		}
	}

#pragma omp parallel for num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{
		mydouble uuu;
		for (long jCell = 0; jCell < nb; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_ddu[iCell][nb][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
			for (long nn = 0; nn < 4; nn++)
			{

				uuu = D_ddu[iCell][nCellYo_L - 1][nn];
				D_ddu[iCell][jCell][nn] = uuu;
			}
	}
}

// --------------------Boundary_u_t-v_t-T_t--------------------

void CADSolver::Boundary_u_t(void)
{
	// #pragma omp barrier

	switch (BoundaryScheme)
	{
	case 0:
		Boundary_u_t_periodic();
		break;
		// case 1:
		// 	Boundary_uu_Expand();
		// 	break;
		// case 2:
		// 	Boundary_uu_reflection_x();
		// 	break;
		// case 24:
		// 	Boundary_uu_Viscous_shock();
		// 	break;

		// case 22:
		// 	Boundary_uu_DoubleMach();
		// 	break;

		// case 23:
		// 	Boundary_uu_High_jet();
		// 	break;
		// case 26:
		// 	Boundary_uu_Couette_flow();
		// 	break;
		// case 27:
		// 	Boundary_uu_Laminar();
		// 	break;

	default:
		cout << "Boundary_uu is not avaliable for this Solver." << endl;
		break;
	}

	// #pragma omp barrier
}

void CADSolver::Boundary_u_t_periodic(void)
{

#pragma omp parallel for							  // num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	{
		for (long iCell = 0; iCell < nb; iCell++)
		{
			D_ut[iCell][jCell] = D_ut[nCell + iCell][jCell];
			D_vt[iCell][jCell] = D_vt[nCell + iCell][jCell];
			D_Tt[iCell][jCell] = D_Tt[nCell + iCell][jCell];
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			D_ut[iCell][jCell] = D_ut[iCell - nCell][jCell];
			D_vt[iCell][jCell] = D_vt[iCell - nCell][jCell];
			D_Tt[iCell][jCell] = D_Tt[iCell - nCell][jCell];
		}
	}

#pragma omp parallel for							// num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++) // y方向
	{
		for (long jCell = 0; jCell < nb; jCell++)
		{
			D_ut[iCell][jCell] = D_ut[iCell][jCell + nCellYo];
			D_vt[iCell][jCell] = D_vt[iCell][jCell + nCellYo];
			D_Tt[iCell][jCell] = D_Tt[iCell][jCell + nCellYo];
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
		{
			D_ut[iCell][jCell] = D_ut[iCell][jCell - nCellYo];
			D_vt[iCell][jCell] = D_vt[iCell][jCell - nCellYo];
			D_Tt[iCell][jCell] = D_Tt[iCell][jCell - nCellYo];
		}
	}
}

void CADSolver::Boundary_u_tt(void)
{
	// #pragma omp barrier

	switch (BoundaryScheme)
	{
	case 0:
		Boundary_u_tt_periodic();
		break;
		// case 1:
		// 	Boundary_uu_Expand();
		// 	break;
		// case 2:
		// 	Boundary_uu_reflection_x();
		// 	break;
		// case 24:
		// 	Boundary_uu_Viscous_shock();
		// 	break;

		// case 22:
		// 	Boundary_uu_DoubleMach();
		// 	break;

		// case 23:
		// 	Boundary_uu_High_jet();
		// 	break;
		// case 26:
		// 	Boundary_uu_Couette_flow();
		// 	break;
		// case 27:
		// 	Boundary_uu_Laminar();
		// 	break;

	default:
		cout << "Boundary_uu is not avaliable for this Solver." << endl;
		break;
	}

	// #pragma omp barrier
}

void CADSolver::Boundary_u_tt_periodic(void)
{

#pragma omp parallel for							  // num_threads(NumberThreads)
	for (long jCell = 0; jCell < nCellYo_2L; jCell++) // X方向
	{
		for (long iCell = 0; iCell < nb; iCell++)
		{
			D_utt[iCell][jCell] = D_utt[nCell + iCell][jCell];
			D_vtt[iCell][jCell] = D_vtt[nCell + iCell][jCell];
			D_Ttt[iCell][jCell] = D_Ttt[nCell + iCell][jCell];
		}

		for (long iCell = nCell_L; iCell < nCell_2L; iCell++)
		{
			D_utt[iCell][jCell] = D_utt[iCell - nCell][jCell];
			D_vtt[iCell][jCell] = D_vtt[iCell - nCell][jCell];
			D_Ttt[iCell][jCell] = D_Ttt[iCell - nCell][jCell];
		}
	}

#pragma omp parallel for							// num_threads(NumberThreads)
	for (long iCell = 0; iCell < nCell_2L; iCell++) // y方向
	{
		for (long jCell = 0; jCell < nb; jCell++)
		{
			D_utt[iCell][jCell] = D_utt[iCell][jCell + nCellYo];
			D_vtt[iCell][jCell] = D_vtt[iCell][jCell + nCellYo];
			D_Ttt[iCell][jCell] = D_Ttt[iCell][jCell + nCellYo];
		}

		for (long jCell = nCellYo_L; jCell < nCellYo_2L; jCell++)
		{
			D_utt[iCell][jCell] = D_utt[iCell][jCell - nCellYo];
			D_vtt[iCell][jCell] = D_vtt[iCell][jCell - nCellYo];
			D_Ttt[iCell][jCell] = D_Ttt[iCell][jCell - nCellYo];
		}
	}
}
