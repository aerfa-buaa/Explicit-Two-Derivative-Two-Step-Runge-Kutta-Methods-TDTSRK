// #pragma once
#include "xBurgers.hpp"
using namespace std;

void CADSolver::Initial_date(void)
{

	mydouble xp, yp;
	mydouble unp, un;
	mydouble pressure, density, velocity_x, velocity_y;

	cout << "--------------- Initial ----------------" << endl;

	// // #pragma omp parallel for num_threads(NumberThreads)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		unsigned short iVar;
		mydouble xp, yp;
		mydouble unp, un, uv;
		mydouble pressure, density, velocity_x, velocity_y, cc;
		mydouble r, deltaT, ka, d1, u1, p1, T1, v1, TT;
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			xp = D_Coord[iCell][jCell][0];
			yp = D_Coord[iCell][jCell][1];

			// ka = 5.; // ka = 5.;
			// d1 = 1.00, u1 = 1.0, p1 = 1.00, T1 = p1 * Gamma / d1;
			// v1 = 1.0;
			// r = (pow(xp - 0.0, 2) + pow(yp - 0.0, 2)) * 1.0; // * 1.0;
			// deltaT = -(Gamma - 1.0) * ka * ka * exp(1.0 - r) / Gamma / PI_Number / PI_Number / 8.0;
			//
			// TT = (T1 + deltaT);
			// density = pow((1.0 + deltaT), (1.0 / (Gamma - 1.0)));
			// velocity_x = u1 - 0.50 * ka / PI_Number * exp(0.50 * (1 - r)) * (yp - 0.0);
			// velocity_y = v1 + 0.50 * ka / PI_Number * exp(0.50 * (1 - r)) * (xp - 0.0);
			// pressure = pow((1.0 + deltaT), (Gamma / (Gamma - 1.0)));
			//
			// density = 2.0 + 0.2 * sin(PI_Number * (xp + yp));
			// velocity_x = 0.7;
			// velocity_y = 0.3;
			// pressure = 1.0;
			//
			// density = 2.0 + 0.5 * sin(2. * PI_Number * (yp)/ 2.  );
			// velocity_x = 0.0;
			// velocity_y = 1.0;
			// pressure = 1.0;

			switch (Case_name)
			{
			case 1:
				case_sod(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;

			case 11:
				case_sod_2D(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;
			case 7:
				case_Shock_vortex_2D(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;

			case 2:
				case_Blast_wave(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;
			case 22:
				case_DoubleMach(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;
			case 23:
				case_High_jet(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;

			case 24:
				case_Viscous_shock(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;

			case 25:
				case_Kelvin_Helmholtz(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;

			case 26:
				case_Couette_flow(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;

			case 27:
				case_Laminar(xp, yp);
				density = S_d, velocity_x = S_u, velocity_y = S_v, pressure = S_p;
				break;

			case 80:

				density = 2. + sin(PI_Number * xp);
				velocity_x = 1.0;
				velocity_y = 0.0;
				pressure = 1.0;

				break;

			case 81:

				density = 2. + sin(PI_Number * yp);
				velocity_x = 0.0;
				velocity_y = 1.0;
				pressure = 1.0;

				break;

			case 82:
				density = 2.0 + 0.5 * sin(PI_Number * (xp + yp));
				velocity_x = 0.5;
				velocity_y = 0.5;
				pressure = 1.0;
				break;

			case 83:

				ka = 5.; // ka = 5.;
				d1 = 1.00, u1 = 1.0, p1 = 1.00, T1 = p1 * Gamma / d1;
				v1 = 1.0;
				r = (pow(xp - 0.0, 2) + pow(yp - 0.0, 2)) * 1.0; // * 1.0;
				deltaT = -(Gamma - 1.0) * ka * ka * exp(1.0 - r) / Gamma / PI_Number / PI_Number / 8.0;

				TT = (T1 + deltaT);
				density = pow((1.0 + deltaT), (1.0 / (Gamma - 1.0)));
				velocity_x = u1 - 0.50 * ka / PI_Number * exp(0.50 * (1. - r)) * (yp - 0.0);
				velocity_y = v1 + 0.50 * ka / PI_Number * exp(0.50 * (1. - r)) * (xp - 0.0);
				pressure = pow((1.0 + deltaT), (Gamma / (Gamma - 1.0)));
				break;

			default:
				cout << "Initial condition is not avaliable for this Solver." << endl;
				break;
			}

			iVar = 0;
			D_un[iCell][jCell][iVar] = density;

			iVar = 1;
			un = density * velocity_x;
			D_un[iCell][jCell][iVar] = un;

			iVar = 2;
			uv = density * velocity_y;
			D_un[iCell][jCell][iVar] = uv;

			iVar = 3;
			unp = pressure / Gamma_1 + 1.0 / 2.0 * density * (velocity_x * velocity_x + velocity_y * velocity_y);
			D_un[iCell][jCell][iVar] = unp;

			cc = sqrt(Gamma * pressure / density);

			D_duvpc[iCell][jCell][0] = density;
			D_duvpc[iCell][jCell][1] = velocity_x;
			D_duvpc[iCell][jCell][2] = velocity_y;
			D_duvpc[iCell][jCell][3] = pressure;
			D_duvpc[iCell][jCell][4] = cc;

			D_Exact[iCell][jCell][0] = density;
			D_Exact[iCell][jCell][1] = velocity_x;
			D_Exact[iCell][jCell][2] = velocity_y;
			D_Exact[iCell][jCell][3] = pressure;
			D_Exact[iCell][jCell][4] = cc;

			// Element[iCell][jCell]->SetSolution_Old(0, 0, 0);
		}
	}
}

void CADSolver::case_sod(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble d1, u1, p1, d2, u2, p2;

	// cout << "--------------- Sod ----------------" << endl;
	StopTime = 0.2;

	d1 = 1.0, u1 = 0.00, p1 = 1.0;
	d2 = 0.125, u2 = 0.00, p2 = 0.1;

	if (xp < (nx_end - nx_begin) / 2.0)
	{
		uu = u1, pp = p1, dd = d1, vv = 0.;
	}
	else
	{
		uu = u2, pp = p2, dd = d2, vv = 0.;
	}

	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
}

void CADSolver::case_sod_2D(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble d1, u1, p1, d2, u2, p2;

	// cout << "--------------- Sod ----------------" << endl;
	StopTime = 0.2;

	if (xp >= (nx_end - nx_begin) / 2.0 && yp >= (ny_end - ny_begin) / 2.0)
	{
		dd = 1.0, uu = 0.75, vv = -0.5, pp = 1.0;
	}
	if (xp < (nx_end - nx_begin) / 2.0 && yp >= (ny_end - ny_begin) / 2.0)
	{
		dd = 2.0, uu = 0.75, vv = 0.5, pp = 1.0;
	}
	if (xp < (nx_end - nx_begin) / 2.0 && yp < (ny_end - ny_begin) / 2.0)
	{
		dd = 1.0, uu = -0.75, vv = 0.5, pp = 1.0;
	}
	if (xp >= (nx_end - nx_begin) / 2.0 && yp < (ny_end - ny_begin) / 2.0)
	{
		dd = 3.0, uu = -0.75, vv = -0.5, pp = 1.0;
	}

	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
}

void CADSolver::case_Shock_vortex_2D(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble d1, u1, p1, d2, u2, p2, rc, ka, mu, M1, y, T1, v1, M2, T2, v2;

	Ma = 1.1;
	rc = 0.05;
	ka = 0.3;
	mu = 0.204;
	d1 = 1.0, u1 = sqrt(Gamma) * Ma, p1 = 1.0;
	T1 = p1 * Gamma / d1, v1 = 0.;

	M1 = Ma;
	y = Gamma;
	d2 = ((y + 1.) * M1 * M1) / ((y - 1.) * M1 * M1 + 2.) * d1;
	u2 = u1 / (((y + 1.) * M1 * M1) / ((y - 1.) * M1 * M1 + 2.));
	p2 = (2. * y) / (y + 1.) * (M1 * M1 - 1.) * p1 + p1;
	T2 = p2 * Gamma / d2, v2 = 0.;

	// cout << "--------------- case_Shock_vortex_2D ----------------" << endl;

	if (xp < 0.5)
	{
		mydouble r = sqrt((xp - 0.25) * (xp - 0.25) + (yp - 0.5) * (yp - 0.5));
		mydouble ta = r / rc;
		mydouble deltaT = -(Gamma - 1.0) * ka * ka * exp(2. * mu * (1. - ta * ta)) / Gamma / mu / 4.0;

		T1 = (T1 + deltaT);
		d1 = pow((1. + deltaT), (1.0 / (Gamma - 1.0)));
		p1 = T1 * d1 / Gamma;
		u1 = u1 + ka / rc * exp(mu * (1. - ta * ta)) * (yp - 0.5);
		v1 = v1 - ka / rc * exp(mu * (1. - ta * ta)) * (xp - 0.25);

		uu = u1, pp = p1, dd = d1, vv = v1;
	}
	else
	{
		uu = u2, pp = p2, dd = d2, vv = 0.;
	}

	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
}

void CADSolver::case_Viscous_shock(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble d1, u1, p1, d2, u2, p2;

	// cout << "--------------- case_Viscous_shock ----------------" << endl;

	d1 = 120.0, u1 = 0.00, p1 = d1 / Gamma;
	d2 = 1.2, u2 = 0.00, p2 = d2 / Gamma;

	if (xp < (nx_end - nx_begin) / 2.0)
	{
		uu = u1, pp = p1, dd = d1, vv = 0.;
	}
	else
	{
		uu = u2, pp = p2, dd = d2, vv = 0.;
	}

	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
}

void CADSolver::case_Kelvin_Helmholtz(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble d1, u1, p1, d2, u2, p2;

	// cout << "--------------- case_Kelvin_Helmholtz ----------------" << endl;

	if (abs(yp) < 0.250001)
	{
		uu = -0.5, vv = 0.01 * sin(2. * PI_Number * xp), pp = 2.5, dd = 2.;
	}
	else
	{
		uu = 0.5, vv = 0.01 * sin(2. * PI_Number * xp), pp = 2.5, dd = 1.;
	}

	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
}

void CADSolver::case_Couette_flow(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble H, u1, mu, ka, T0, T1, Ta;

	// cout << "--------------- case_Kelvin_Helmholtz ----------------" << endl;
	H = 2., u1 = 0.1;
	uu = yp * u1 / H, vv = 0.;
	T0 = 0.8, T1 = 0.85;
	Ta = T0 + yp * (T1 - T0) / H + Pr * Gamma_1 * u1 * u1 * yp / (2. * H) * (1 - yp / H);

	pp = 1. / Gamma, dd = Gamma * pp / Ta;
	// D_ka[iCell][jCell] = C_p * D_mu[iCell][jCell] / Pr;
	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
	// cout << "yp = " << yp << " uu = " << uu << endl;
}

void CADSolver::case_Laminar(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble H, u1, mu, ka, T0, T1, Ta;

	// cout << "--------------- case_Kelvin_Helmholtz ----------------" << endl;
	uu = Ma, vv = 0.;
	pp = 1., dd = 1.;
	// D_ka[iCell][jCell] = C_p * D_mu[iCell][jCell] / Pr;
	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
	// cout << "yp = " << yp << " uu = " << uu << endl;
}
void CADSolver::case_High_jet(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble d1, u1, p1, d2, u2, p2;

	uu = 0.0, vv = 0., pp = 0.4127, dd = 0.5;

	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
}

void CADSolver::case_DoubleMach(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble alfa, x0, x1, VV2, T2, d2, u2, v2, r, Ma;

	// cout << "--------------- case_DoubleMach ----------------" << endl;

	alfa = 60. / 180. * PI_Number;

	x0 = 1. / 6.;
	x1 = x0 + yp / tan(alfa);
	if (xp < x1)
	{
		dd = 8.0, uu = 8.25 * sin(alfa), vv = -8.25 * cos(alfa), pp = 116.5;
	}
	else
	{
		dd = 1.4, uu = 0.0, vv = 0.0, pp = 1.0;
	}

	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
}

void CADSolver::case_Blast_wave(mydouble xp, mydouble yp)
{

	mydouble pp, dd, uu, vv;
	mydouble d1, u1, p1, d2, u2, p2, a1;

	// cout << "--------------- Blast_wave ----------------" << endl;
	StopTime = 0.38;

	if (xp <= 1.0)
	{
		d1 = 1., u1 = 0., p1 = 1000.;
		uu = u1, pp = p1, dd = d1, vv = 0.;
	}
	else if (xp <= 9.0)
	{
		d1 = 1., u1 = 0., p1 = 0.01;
		uu = u1, pp = p1, dd = d1, vv = 0.;
	}
	else
	{
		d1 = 1., u1 = 0., p1 = 100.;
		uu = u1, pp = p1, dd = d1, vv = 0.;
	}

	S_d = dd, S_u = uu, S_v = vv, S_p = pp;
}

void CADSolver::MesherAndAllocation_xy(void)
{

	mydouble xp, yp, nxx, nyy;

	nxx = double(nCell);
	nyy = double(nCellYo);
	dx = (nx_end - nx_begin) / nxx;
	dy = (ny_end - ny_begin) / nyy;

	// // #pragma omp parallel for
	for (long iCell = 0; iCell < nCell_2L; iCell++)
	{

		for (long jCell = 0; jCell < nCellYo_2L; jCell++)
		{
			{
				xp = nx_begin + (double(iCell) - double(nb)) * dx;
				yp = ny_begin + (double(jCell) - double(nb)) * dy;

				D_Coord[iCell][jCell][0] = xp;
				D_Coord[iCell][jCell][1] = yp;

				// Element[iCell][jCell]->SetCoordinate(0, xp); // SetCoordinate_xy(0,0,0, xp);
				// Element[iCell][jCell]->SetCoordinate(1, yp); // SetCoordinate_xy(0,0,0, yp);
				// Element[iCell][jCell]->SetCoordinate(0, 1, 0, dx); // SetCoordinate_dxdy(0, dx);
				// Element[iCell][jCell]->SetCoordinate(1, 1, 0, dy);
			}
		}
	}
}

void CADSolver::Mesher_out(void)
{

	ofstream gridflow;
	stringstream meshflow;
	string meshname = "int_Mesh_";
	string buffer = ".dat";
	string mesh;
	meshflow << meshname << nCell << "_" << nCellYo;
	mesh = meshflow.str() + buffer;
	meshflow.clear();

	gridflow.open(mesh.c_str(), ios_base::out);
	gridflow.precision(15);

	gridflow << "TITLE = \"Visualization of the mesh\"" << endl;
	gridflow << "VARIABLES = ";
	char varmesh[] = "\"x\",\"y\"";
	gridflow << varmesh;
	gridflow << endl;
	gridflow << "ZONE  I = " << nCell_2L << " J =  " << nCellYo_2L << " F = POINT" << endl;

	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nCell_2L; iCell++)
		{
			{
				// xp = nx_begin + (iCell - nb) * dx;
				// yp = ny_begin + (jCell - nb) * dy;

				gridflow << scientific << D_Coord[iCell][jCell][0] << "\t";
				gridflow << scientific << D_Coord[iCell][jCell][1] << "\t";
				gridflow << endl;
			}
		}
	}

	gridflow.close();

	// cout << TotalVol << endl;
}

CADSolver::~CADSolver(void)
{
}

CADSolver::CADSolver(void)
{

	/*
	Aspeed          = 1.0      ;
	ReExpo          = 2        ;
	Alpha           = 4.5      ;
	CFL	            = 0.75     ;
	nMax            = 10000    ;
	nCell           = 100      ;
	ResLimitMag     = 5        ;
	ResLimitMin     = -7       ; */

	PI_Number = 4.0 * atan(1.0);
	Gamma = 1.4;
}

void CADSolver::comput_err1(void) // uu_to_cc
{

	mydouble errl = 0, un = 0., err_max = 0., ea = 0.;
	// #pragma omp parallel for private(pressure, density, velocity_x, velocity_y, cc, lmax, jCell, iCell)
	// #pragma omp parallel for reduction(+ : errl, un)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			// ea = fabs(Element[iCell][jCell]->GetSolution(0, 0, 0) - Element[iCell][jCell]->GetSolution_Exact(0));
			ea = fabs(D_un[iCell][jCell][0] - D_Exact[iCell][jCell][0]);
			errl += ea;
			un += 1.0;
			// if (err_max < ea)
			// 	err_max = ea;
		}
	}
	// #pragma omp parallel for shared(err_max)
	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			// #pragma omp critical
			{
				// ea = fabs(Element[iCell][jCell]->GetSolution(0, 0, 0) - Element[iCell][jCell]->GetSolution_Exact(0));
				ea = fabs(D_un[iCell][jCell][0] - D_Exact[iCell][jCell][0]);
				if (err_max < ea)
					err_max = ea;
			}
		}
	}

	errl = errl / un;
	cout << "ele_all = " << un << endl;
	cout << "nx = " << nCell << " ny = " << nCellYo << "   err1 = " << scientific << errl << "   err_inf = " << scientific << err_max << endl;
}

void CADSolver::out_date(void)
{

	ofstream gridflow;
	stringstream meshflow;
	string meshname = "final_date_";
	string buffer = ".dat";
	string mesh;
	meshflow << meshname << nCell << "_" << nCellYo;
	mesh = meshflow.str() + buffer;
	meshflow.clear();

	gridflow.open(mesh.c_str(), ios_base::out);
	gridflow.precision(15);

	gridflow << "TITLE = \"Visualization of the mesh\"" << endl;
	gridflow << "VARIABLES = ";
	char varmesh[] = " \"x\",\"y\" ,\"D\" ,\"u\" ,\"v\" ,\"p\" ,\"D_ext\" ";
	gridflow << varmesh;
	gridflow << endl;
	gridflow << "ZONE  I = " << nCell << " J =  " << nCellYo << " F = POINT" << endl;

	for (long jCell = nb; jCell < nCellYo_L; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			{

				gridflow << scientific << D_Coord[iCell][jCell][0] << "\t";
				gridflow << scientific << D_Coord[iCell][jCell][1] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][0] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][1] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][2] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][3] << "\t";
				// gridflow << scientific << D_duvpc[iCell][jCell][5] << "\t";
				gridflow << scientific << D_Exact[iCell][jCell][0] << "\t";
				// gridflow << scientific << D_out[iCell][jCell][0] << "\t";

				gridflow << endl;
			}
		}
	}

	gridflow.close();

	// cout << TotalVol << endl;

	out_date_1D();
	out_date_expand();
}

void CADSolver::out_date_1D(void)
{

	ofstream gridflow;
	stringstream meshflow;
	string meshname = "final_1D_date_";
	string buffer = ".dat";
	string mesh;
	meshflow << meshname << nCell << "_" << nCellYo;
	mesh = meshflow.str() + buffer;
	meshflow.clear();

	gridflow.open(mesh.c_str(), ios_base::out);
	gridflow.precision(15);

	gridflow << "TITLE = \"Visualization of the mesh\"" << endl;
	gridflow << "VARIABLES = ";
	char varmesh[] = " \"x\" ,\"D\" ,\"u\" ,\"p\",\"out0\",\"out1\" ,\"D_ext\" ";
	gridflow << varmesh;
	gridflow << endl;
	// gridflow << "ZONE  I = " << nCell << " F = POINT" << endl;

	for (long jCell = nb; jCell < nb + 1; jCell++)
	{
		for (long iCell = nb; iCell < nCell_L; iCell++)
		{
			{

				gridflow << scientific << D_Coord[iCell][jCell][0] << "\t";
				// gridflow << scientific << D_Coord[iCell][jCell][1] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][0] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][1] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][3] << "\t";
				gridflow << scientific << D_out[iCell][jCell][0] << "\t";
				gridflow << scientific << D_out[iCell][jCell][1] << "\t";
				gridflow << scientific << D_Exact[iCell][jCell][0] << "\t";

				gridflow << endl;
			}
		}
	}

	gridflow.close();
}

void CADSolver::out_date_expand(void)
{

	ofstream gridflow;
	stringstream meshflow;
	string meshname = "final_date_Expand_";
	string buffer = ".dat";
	string mesh;
	meshflow << meshname << nCell << "_" << nCellYo;
	mesh = meshflow.str() + buffer;
	meshflow.clear();

	gridflow.open(mesh.c_str(), ios_base::out);
	gridflow.precision(15);

	gridflow << "TITLE = \"Visualization of the mesh\"" << endl;
	gridflow << "VARIABLES = ";
	char varmesh[] = " \"x\",\"y\" ,\"D\" ,\"u\" ,\"v\" ,\"D_ext\",\"out0\" ";
	gridflow << varmesh;
	gridflow << endl;
	gridflow << "ZONE  I = " << nCell_2L << " J =  " << nCellYo_2L << " F = POINT" << endl;

	for (long jCell = 0; jCell < nCellYo_2L; jCell++)
	{
		for (long iCell = 0; iCell < nCell_2L; iCell++)
		{
			{

				gridflow << scientific << D_Coord[iCell][jCell][0] << "\t";
				gridflow << scientific << D_Coord[iCell][jCell][1] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][0] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][1] << "\t";
				gridflow << scientific << D_duvpc[iCell][jCell][2] << "\t";
				gridflow << scientific << D_Exact[iCell][jCell][0] << "\t";
				gridflow << scientific << D_out[iCell][jCell][0] << "\t";
				// gridflow << scientific << D_gx_t[iCell][jCell][1] << "\t";
				// gridflow << scientific << D_gx[iCell][jCell][1]  << "\t";

				gridflow << endl;
			}
		}
	}

	gridflow.close();

	// cout << TotalVol << endl;

	out_date_1D();
}
