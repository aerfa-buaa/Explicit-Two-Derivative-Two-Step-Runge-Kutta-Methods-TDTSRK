#include <math.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>
#include <iterator>
#include <valarray>
using namespace std;
// example:line 3-4-L-W   增加CFL数输出误差   Line Excample-1  精简版本只有1up,2cneter  U_t+U_x=0
class cfda
{

private:
    int w;

public:
    double PI = 3.141592653589793;
    // double PI = 3.1415926535897932385;

    int n;                                         // 网格结点数（有限差分）
    int nb;                                        // 边界网格结点数（有限差分）
    double t, tt, CFL;                             // 计算总时 间t,实际计算时间tt;
    int nt, kt, ij, ijj, onet_out = 2;             // 计算时间步长;
    double nx_begin, nx_end;                       // 网格长度;
    int n1;                                        // 总的网格结点数（有限差分）
    double *Ltt_n0, *x;                            // 声明变长数组
    double *L_n0, *L_n1, *L_n2, *L_n3, *L_n4;      // f=L
    double *Lt_n0, *Lt_n1, *Lt_n2, *Lt_n3, *Lt_n4; // L_t

    double *TL_n0, *TL_n1, *TL_n2, *TL_n3, *TL_n4;      // Tf=L
    double *TLt_n0, *TLt_n1, *TLt_n2, *TLt_n3, *TLt_n4; // TL_t

    double *u_n0, *u_n1, *u_n2, *u_n3, *u_n4, *u_n5, *u_n6, *u_n7, *u_n8, *u_n9; // 推进步中间时刻
    double *f_n0, *f_n1, *f_n2, *f_n3, *f_n4, *f_n5, *f_n6, *f_n7, *f_n8, *f_n9; // 推进步中间时刻
    double *u_nn, *f_nn, *Lttx_n0;                                               // 下一时刻（n+1)
    double *Tu1, *Tu2;
    // double  u_excat[513 + 6 * 2];
    double *df, *u_excat; // du[513 + 6 * 2],df[513 + 6 * 2],
    ////-----------------------------------------
    double dx;
    double dt;
    double *dd; // ss
    double A1, A2, A3, A4, A5, B0, B1, B2, B3, C0, C1, C2, C3, D0, D1, D2, D3, maxuu, TVmax, k, b;
    double a21, a32, a31, aa12, aa21, aa31, aa32, w1, w2, v1, v2, v3, vv1, vv2, vv3, ww1, ww2;
    double put_obj[2][1300];   //
    double put_one_n[2][1300]; // onet_out = 10;
    ////-----------------------------------------
    void intc()
    {
        dx = (nx_end - nx_begin) / (1.0 * (n - 1));

        // for (int i = 0; i < n1; i++)
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     u_n0[i] = intc_fun(x[i]);
        //     df[i] = 1;
        // }

        for (int i = 0; i < n1; i++) // 间断
        {
            x[i] = nx_begin + (i - nb) * dx;
            u_n0[i] = 0.0;
            if (x[i] >= -0.5 && x[i] <= 0.0)
            {
                u_n0[i] = 1.0;
            }
            // df[i] = 1;
        }
        // Write_obj(0);
    }

    double intc_fun(double xxf)
    {
        double yyf;
        // yyf = (sin(2 * PI * xxf) / 2.0 + 0.5) * 1.0; //  / PIbergure
        // yyf = sin(PI * xxf) / 2.0 + 0.25;//bergure
        yyf = sin(PI * xxf); // line
        return yyf;
    }

    // //-----------------------------------------
    void border()
    {
        borderfun(u_n0), borderfun(u_n1), borderfun(u_n2), borderfun(u_n3), borderfun(u_n4), borderfun(u_n5);
        borderfun(u_n6), borderfun(u_n7), borderfun(u_n8), borderfun(u_n9), borderfun(u_nn);
        borderfun(f_n0), borderfun(f_n1), borderfun(f_n2), borderfun(f_n3), borderfun(f_n4), borderfun(f_n5);
        borderfun(f_n6), borderfun(f_n7), borderfun(f_n8), borderfun(f_n9), borderfun(f_nn);
        borderfun(L_n0), borderfun(L_n1), borderfun(L_n2), borderfun(L_n3), borderfun(L_n4);
        borderfun(Lt_n0), borderfun(Lt_n1), borderfun(Lt_n2), borderfun(Lt_n3), borderfun(Lt_n4);
        borderfun(Ltt_n0), borderfun(Lttx_n0);
    }

    void borderfun(double *ffbc)
    {

        for (int i = 0; i < nb; i++)
        {

            *(ffbc + i) = *(ffbc + n1 - 2 * nb + i - 1);  // 周期边条
            *(ffbc + n1 - nb + i) = *(ffbc + nb + i + 1); // 周期边条

            // *(ffbc + i) = 0.;
            // *(ffbc + n1 - nb + i) = 0.;

            // f_n0[i] = f_n0[nb];                   //恒定边条
            // f_n0[n1 - 1 - i] = f_n0[n1 - 1 - nb]; //恒定边条

            // *(ffbc + i) = *(ffbc + nb);                   //恒定边条
            // *(ffbc + n1 - 1 - i) = *(ffbc + n1 - 1 - nb); //恒定边条

            // *(ffbc + i) = 0.0;                   //恒定边条
            // *(ffbc + n1 - 1 - i) = 0.0; //恒定边条
            // ffbc++;
        }
    }

    // //-----------------------------------------
    void carr(double *p1, double *p2, int i)
    {
        while (i-- > 0)
        {
            *p1++ = *p2++;
        }
    }
    // //-----------------------------------------
    void compt_2_stage_int()
    {
        int i, j;
        a21 = 0.516075574501012; // 0.525  0.680; //0.765;   // 1. / 210. * (93. + sqrt(2139.));
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);

            Tu1[j] = f_n0[j];
            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_SGLM_int()
    {
        int i, j;
        a21 = 0.85; // 0.532  0.468; //0.765;   // 1. / 210. * (93. + sqrt(2139.));
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            // Tu1[j] = u_n0[j] + 0.270774947096390 * dt * L_n0[j] + 0.104435936229161 * dt * dt * Lt_n0[j];
            // Tu1[j] = u_n0[j] + 0.527015480634761 * dt * L_n0[j] + 0.027140536548443 * dt * dt * Lt_n0[j];
            // Tu1[j] = u_n0[j] + 0.032910530894256 * dt * L_n0[j] + 0.000541551521871 * dt * dt * Lt_n0[j];
            Tu1[j] = u_n0[j] + 0.093858080521929 * dt * L_n0[j]; //- 0.008889300639912 * dt * dt * Lt_n0[j];
        }
    }

    void compt23_LLttnn_line_SSP()
    {

        int i, j;

        a21 = 0.594223212099088, v2 = 0.306027487008159;
        ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            // Lt_n0[j] = -1. * dfdxcent(L_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt24_LLttnn_line_SSP()
    {
        int i, j;

        A1 = 0.5, B0 = 1., B1 = 0., C0 = 1. / 3., C1 = 2. / 3.;
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A1 * A1 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j] + B1 * L_n1[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt34_LLttnn_line_SSP()
    {
        int i, j;
        // A1 = 1. / 4., A2 = 1. / 2., B0 = 1. / 3.; //稳定性足够  B0 = 1. / 5.;  B0 = 1. / 3.;
        // A1 = 1. / 4., A2 = 2. / 3., B0 = 1. / 3.;
        A1 = 0.5, A2 = 0.8, B0 = 0.6;
        B1 = -((pow(A2, 3) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3)));
        B2 = -((pow(A1, 3) * (-1. + B0)) / (pow(A1, 3) - pow(A2, 3)));
        C0 = -((-1. + 2. * A1 + 2. * A2 - 6. * A1 * A2) / (6. * A1 * A2)) + (A2 * (-pow(A1, 3) + A1 * A2 * A2) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        C1 = -((-1. + 2. * A2) / (6. * A1 * (A1 - A2))) + (A1 * pow(A2, 3) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        C2 = -((1. - 2. * A1) / (6. * (A1 - A2) * A2)) - (pow(A1, 3) * A2 * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A1 * A1 / 2.0 * dt * dt * Lt_n0[j];
            u_n2[j] = u_n0[j] + A2 * dt * L_n0[j] + A2 * A2 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]), f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            L_n2[j] = -dfdx(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            Lt_n2[j] = dfdxx(f_n2, j);

            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j] + B1 * L_n1[j] + B2 * L_n2[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j] + C2 * Lt_n2[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
            // Tu2[j] = Tu1[j];
            // Tu1[j] = u_nn[j];
        }
    }

    void compt35_LLttnn_line_SSP()
    {
        int i, j;

        a21 = 0.7504;
        aa12 = -1. + 2. * a21;
        vv1 = (1. + 2. * a21 * (-4. + 5. * a21)) / (12. * a21 * (-3. + 5. * a21));
        vv2 = 1. / (12. * a21 * (3. + 10. * (-1. + a21) * a21));
        vv3 = (25. * pow(aa12, 3)) / (12. * (-3. + 5. * a21) * (3. + 10. * (-1 + a21) * a21));
        v1 = 1., v2 = 0., v3 = 0.;
        a31 = (3. - 5. * a21) / (5. - 10. * a21), a32 = 0.;
        aa31 = ((-3. + 5. * a21) * (-3. + 10. * a21) * (1. + 5. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));
        aa32 = ((-3. + 5. * a21) * (3. + 10. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -1. * dfdx(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n2[j] = dfdxx(f_n2, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + v3 * L_n2[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j]);

            // TL_n1a[j] = L_n1[j], TL_n1[j] = L_n0[j];
            // TLt_n1a[j] = Lt_n1[j], TLt_n1[j] = Lt_n0[j];
        }
    }

    // -------------------------------------
    void compt_TDTS23_LLttnn_line()
    {
        int i, j;
        // a21 = 0.532, v2 = 0.389; // a21 = 0.365, v2 = 0.57;   a21 = 0.389, v2 = 0.532;
        // ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.1;

        // 0.473394419966701	0.765968261590665	-0.150000000000000	0.00699024210121857	-0.0893626815573664	-0.150000000000000
        // 0.526068099413695	0.389919021445930	-5.01236625573457e-27	4.29145094976513e-28	2.32184837825380e-17	0.0840128791403755

        a21 = 0.525, v2 = 0.448;
        a21 = 0.52489965455223, v2 = 0.447659867744687;
        ww1 = 0.00, ww2 = 0.00, w1 = 0.04, w2 = 0.00;
        // 0.516075574501012	0.504049314006236	0	0	0.0798751114927521	-0.100000000000000	-1.19756489027303	-1.24373977374248
        // a21 = 0.516075574501012, v2 = 0.504049314006236;
        // ww1 = 0.0, ww2 = 0.0, w1 = 0.1798751114927521, w2 = 0.1;

        v1 = 1. - v2 - w1 - w2;
        vv1 = -((1 - 3 * w1 - 3 * w2 + 6 * ww1 + 3 * a21 * (-1 - 2 * w1 + a21 * (v2 + w2) + 2 * ww1) + 6 * ww2) / (6. * a21));
        vv2 = (1 - 3 * w1 - 3 * w2 + 6 * ww1 + 6 * ww2 - 3 * a21 * (-2 * w2 + a21 * (v2 + w2) + 2 * ww2)) / (6. * a21);

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDTS24_LLttnn_line()
    {
        int i, j;

        // w2 = -0.1, ww2 = 0.0, a21 = 0.468, v1 = 1.0;    //C=0.712
        // w2 = -0.14941, ww2 = 0.14873, a21 = 0.67, v1 = 0.87;
        // w2 = -0.15, ww2 = 0.15, a21 = 0.698, v1 = 0.85; // C=0.884

        // 0.680145029983404 0.907813449803004 - 0.150000000000000 0.150000000000000 - 0.908679710215592

        // 0.590356938016122	0.698462083485141	-0.200000000000000	0.200000000000000	-1.04688189288802
        // 0.68014503	0.90781345	-0.15	0.15
        // 0.627237328	0.78967575	-0.15	0.2	-0.985327181

        a21 = 0.680, v1 = 0.908;
        a21 = 0.680145029983405, v1 = 0.90781344980300;
        w2 = -0.15, ww2 = 0.15;

        v2 = -w2, w1 = 1. - v1;
        vv1 = -((1. - 14. * a21 + 2 * v1 + 6. * a21 * v1 - 2 * w2 - 6. * a21 * w2 + 12. * a21 * ww2) / (12. * a21));
        vv2 = -((-1. - 2 * v1 + 2 * w2 - 12. * a21 * a21 * w2 - 12. * a21 * ww2 + 12. * a21 * a21 * ww2) / (12. * a21 * (1. + a21)));
        ww1 = -((-5. - 4. * a21 + 4. * v1 + 6. * a21 * v1 - 4. * w2 + 6. * a21 * w2 + 12. * ww2 - 12. * a21 * ww2) / (12. * (1. + a21)));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDTS25_LLttnn_line()
    {
        int i, j;
        a21 = 0.765; // a21 = 0.971;
        a21 = 0.76543209876543;
        vv1 = (-31. + 6. * a21 * (-31. + 85. * a21)) / (360. * pow(a21, 2));
        vv2 = 31. / (360. * pow(a21, 2)), v1 = -(1. / 2.) + 31. / (30. * a21);
        v2 = 0., w1 = 3. / 2. - 31. / (30. * a21), w2 = 0.,
        ww1 = (31. + 6. * a21 * (-31. + 35. * a21)) / (360. * pow(a21, 2));
        ww2 = -(31. / (360. * pow(a21, 2)));
        b = 1.0;
        // a21 = 0.85, b = 1.1;
        // vv1 = (-32 + b + 6 * a21 * (-32 + b + 5 * a21 * (16 + b))) / (360 * pow(a21, 2));
        // vv2 = -((-32 + b) / (360 * pow(a21, 2)));
        // v1 = -((-32 + b + 15 * a21 * b) / (30 * a21));
        // v2 = 0;
        // w1 = (-32 - 15 * a21 * (-4 + b) + b) / (30 * a21);
        // w2 = 0;
        // ww1 = (32 - b + 6 * a21 * (-32 - 5 * a21 * (-8 + b) + b)) / (360 * pow(a21, 2));
        // ww2 = (-32 + b) / (360 * pow(a21, 2));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);

            u_nn[j] = b * u_n0[j] + (1 - b) * Tu1[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);
            Tu1[j] = u_n0[j];
            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDTS35_LLttnn_line()
    {
        int i, j;

        a21 = 0.381, vv3 = 0.21, a31 = 0.61;
        v2 = 0., v3 = 0., a32 = 0.;
        vv1 = (-31. + 2. * a21 * (-4 + 85 * a21) - 60. * (a21 - a31) * (1. + a31) * (1. + 2. * a31 + a21 * (2 + 6 * a31)) * vv3) / (60. * a21 * (1. + 2. * a21));
        vv2 = (31 - 60. * a31 * (1. + a31) * (1. + 2. * a31) * vv3) / (60. * a21 * (1. + a21) * (1. + 2. * a21));
        v1 = (13 - 5 * a21 + 60. * (a21 - a31) * a31 * (1. + a31) * vv3) / (5 + 10. * a21);
        ww1 = (-27 + 2. * a21 * (6 + 35 * a21) - 60. * (a21 - a31) * a31 * (3 + 4. * a31 + a21 * (4 + 6 * a31)) * vv3) / (60. * (1. + a21) * (1. + 2. * a21));
        w1 = (-8 + 15 * a21 - 60. * (a21 - a31) * a31 * (1. + a31) * vv3) / (5 + 10. * a21);
        aa31 = (-31 * pow(a21, 2) + 60. * a31 * (pow(a21, 2) + 3. * a21 * (1. + 2. * a21 * (2 + a21)) * a31 - (1. + 3. * a21) * pow(a31, 2)) * vv3) / (360. * a21 * (1. + a21) * (1. + 2. * a21) * vv3);
        aa32 = (31 * pow(a21, 2) - 60. * (a21 - a31) * a31 * (a21 + a31 + 3. * a21 * a31) * vv3) / (360. * a21 * (1. + a21) * (1. + 2. * a21) * vv3);

        // a21 = 0.81, a31 = 0.31;
        // a32 = -0.1, v2 = 0., v3 = 0., aa12 = a31 + a32, aa21 = a21 * a21 / 2.0;
        // vv1 = (-70. * pow(aa12, 2) - 3. * (-9 + 4. * a31 + 4. * a32) + 4. * a21 * (-3 + 15. * pow(a31, 2) + 2. * a31 * (7 + 15. * a32) + a32 * (14 + 15. * a32)) + 10. * pow(a21, 2) * (-7 + 30. * pow(a31, 2) + 6. * a32 * (1 + 5. * a32) + a31 * (6 + 60. * a32))) / (60. * a21 * (aa12) * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32)));
        // vv2 = -((-27. + 70. * pow(a31, 2) + 4. * a31 * (3. + 35. * a32) + 2. * a32 * (6. + 35. * a32)) / (60. * a21 * (a21 - a31 - a32) * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32))));
        // vv3 = (-27. + 2. * a21 * (6 + 35. * a21)) / (60. * (a21 - a31 - a32) * (aa12) * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32)));
        // v1 = (12. + 25. * a31 + 25. * a32 + 5. * a21 * (5. + 4. * a31 + 4. * a32)) / (5. * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32)));
        // w1 = (3. - 5. * a31 - 5. * a32 + 5. * a21 * (-1. + 2. * a31 + 2. * a32)) / (5. * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32)));
        // aa31 = (-70. * pow(a21, 4) * a32 + 9. * pow(aa12, 3) - a21 * pow(aa12, 2) * (27. + 4. * a31 + 4. * a32) + 2. * pow(a21, 3) * (-6. * a32 + 35. * pow(aa12, 2)) + pow(a21, 2) * (16. * pow(a31, 2) + 2. * a32 * (9. + 8. * a32) + a31 * (-9. + 32. * a32))) / (2. * a21 * (-27. + 2. * a21 * (6. + 35. * a21)));
        // aa32 = (-12. * pow(a21, 3) * a32 - 70. * pow(a21, 4) * a32 - 9. * pow(aa12, 3) + 4. * a21 * pow(aa12, 3) + pow(a21, 2) * (-4. * pow(a31, 2) + a31 * (9. - 8. * a32) - 4. * (-9 + a32) * a32)) / (2. * a21 * (-27. + 2. * a21 * (6. + 35. * a21)));

        // w1 = 0.05, a21 = 0.32;
        // ww1 = 0., v2 = 0., v3 = 0., a32 = 0., aa12 = (1. - 2. * a21 + 4. * w1 + 6. * a21 * w1), aa21 = a21 * a21 / 2.0;
        // vv1 = (1. - 8. * a21 + 10. * pow(a21, 2) - 88. * w1 + 64 * a21 * w1 + 300. * pow(a21, 2) * w1 + 10. * pow(w1, 2) + 60. * a21 * pow(w1, 2) + 60. * pow(a21, 2) * pow(w1, 2)) / (12. * a21 * (-3. + 5. * a21 + 15. * w1 + 20. * a21 * w1));
        // vv2 = (-1. + 88. * w1 - 10. * pow(w1, 2)) / (12. * a21 * (-3. + 10. * a21 - 10. * pow(a21, 2) + 15. * w1 + 40. * a21 * w1 + 30. * pow(a21, 2) * w1));
        // vv3 = (25. * pow(aa12, 3)) / (12. * (9. - 45. * a21 + 80. * pow(a21, 2) - 50. * pow(a21, 3) - 90. * w1 + 45. * a21 * w1 + 160. * pow(a21, 2) * w1 - 50. * pow(a21, 3) * w1 + 225. * pow(w1, 2) + 900. * a21 * pow(w1, 2) + 1250. * pow(a21, 2) * pow(w1, 2) + 600. * pow(a21, 3) * pow(w1, 2)));
        // v1 = 1. - w1;
        // a31 = (3. - 5. * a21 - 15. * w1 - 20. * a21 * w1) / (5. * aa12);
        // aa31 = (-9. + 90. * a21 - 320. * pow(a21, 2) + 475. * pow(a21, 3) - 250. * pow(a21, 4) + 135. * w1 - 540. * a21 * w1 + 960. * pow(a21, 2) * w1 + 100. * pow(a21, 3) * w1 - 1250. * pow(a21, 4) * w1 - 675. * pow(w1, 2) - 1350. * a21 * pow(w1, 2) - 1800. * pow(a21, 2) * pow(w1, 2) + 50. * pow(a21, 3) * pow(w1, 2) + 2000. * pow(a21, 4) * pow(w1, 2) + 1125. * pow(w1, 3) + 9000. * a21 * pow(w1, 3) + 25000. * pow(a21, 2) * pow(w1, 3) + 29000. * pow(a21, 3) * pow(w1, 3) + 12000. * pow(a21, 4) * pow(w1, 3)) / (250. * a21 * pow(aa12, 3));
        // aa32 = -(((-1. + 5. * w1) * (9 - 45. * a21 + 80. * pow(a21, 2) - 50. * pow(a21, 3) - 90. * w1 + 45. * a21 * w1 + 160. * pow(a21, 2) * w1 - 50. * pow(a21, 3) * w1 + 225. * pow(w1, 2) + 900. * a21 * pow(w1, 2) + 1250. * pow(a21, 2) * pow(w1, 2) + 600. * pow(a21, 3) * pow(w1, 2))) / (250. * a21 * pow(aa12, 3)));

        // a21 = 0.7520, aa12 = (-1. + 2. * a21);
        // v2 = 0., v3 = 0., w1 = 0., a32 = 0., aa21 = a21 * a21 / 2.0;
        // vv1 = (1. + 2. * a21 * (-4. + 5. * a21)) / (12. * a21 * (-3. + 5. * a21));
        // vv2 = 1. / (12. * a21 * (3. + 10. * (-1. + a21) * a21));
        // vv3 = (25. * pow(aa12, 3)) / (12. * (-3. + 5. * a21) * (3. + 10. * (-1. + a21) * a21));
        // v1 = 1.;
        // a31 = (3. - 5. * a21) / (5. - 10. * a21);
        // aa31 = ((-3. + 5. * a21) * (-3. + 10. * a21) * (1. + 5. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));
        // aa32 = ((-3. + 5. * a21) * (3. + 10. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -dfdx(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n2[j] = dfdxx(f_n2, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + v3 * L_n2[j] + w1 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j] + ww1 * TLt_n1[j]);

            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    //------------------------------------------------
    void compt_SGLM3_line()
    {
        int i, j;
        double a21, a31, a32, aa21, aa31, aa32, b11, b12, b13, b21, b22, b23, bb11, bb12, bb13, bb21, bb22, bb23, u11, u12, u21, u22, u31, u32, v;

        a21 = 0.367036327407744, a31 = 0.305289583777252, a32 = 0.594276924770852;
        aa21 = 0.087361998687504, aa31 = 0.072665036743429, aa32 = 0;
        b11 = 0.378319518763646, b12 = 0.539413445719846, b13 = 0;
        b21 = 0.224801836087141, b22 = 0.437599416854594, b23 = 0.526106658638147;
        bb11 = 0.065956620927602, bb12 = 0, bb13 = 0.156810544725088;
        bb21 = 0.120798754439142, bb22 = 0, bb23 = 0;
        u11 = 0.883562994058718, u12 = 0.116437005941282, u21 = 0.756325839981176;
        u22 = 0.243674160018824, u31 = 0.629088685903634, u32 = 0.370911314096366;
        v = 0.303820705714043;
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            u_n1[j] = u11 * u_n0[j] + u12 * Tu1[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u21 * u_n0[j] + u22 * Tu1[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n3[j] = u31 * u_n0[j] + u32 * Tu1[j] + dt * (a31 * L_n1[j] + a32 * L_n2[j]) + dt * dt * (aa31 * Lt_n1[j] + aa32 * Lt_n2[j]);
            f_n3[j] = f_u(u_n3[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n3[j] = -dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j]);

            Tu1[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j]);
        }
    }

    void compt_SGLM4_line()
    {
        int i, j;
        double a21, a31, a32, a41, a42, a43, aa21, aa31, aa32, aa41, aa42, aa43, b11, b12, b13, b14, b21, b22, b23, b24;
        double bb11, bb12, bb13, bb14, bb21, bb22, bb23, bb24, u11, u12, u21, u22, u31, u32, u41, u42, v;

        a21 = 0.264016422097677, a31 = 0.208728827122131, a32 = 0.280067373408790;
        a41 = 0.156642291124758, a42 = 0.210178898836894, a43 = 0.292626140134642;
        aa21 = 0.069927142101271, aa31 = 0.055054745190716, aa32 = 0.060229350201948;
        aa41 = 0.041316293215779, aa42 = 0.045199618752556, aa43 = 0.049611867625336;
        b11 = 0.149050523837038, b12 = 0.203464159352573, b13 = 0.239948323340150;
        b14 = 0.081851422369239, b21 = 0.148523308458399, b22 = 0.199285041091082;
        b23 = 0.277458929910524, b24 = 0.576062630073757, bb11 = 0.050710732728729;
        bb12 = 0.042125352335649, bb13 = 0.006683865350595, bb14 = 0;
        bb21 = 0.042847486964744, bb22 = 0.042856861132357, bb23 = 0.047040417154307;
        bb24 = 0, u11 = 0.734518464168064, u12 = 0.265481535831936;
        u21 = 0.597488310192903, u22 = 0.402511689807097, u31 = 0.471457186431044;
        u32 = 0.528542813568956, u41 = 0.353808982055841, u42 = 0.646191017944159;
        v = 0.617981032945615;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            u_n1[j] = u11 * u_n0[j] + u12 * Tu1[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u21 * u_n0[j] + u22 * Tu1[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n3[j] = u31 * u_n0[j] + u32 * Tu1[j] + dt * (a31 * L_n1[j] + a32 * L_n2[j]) + dt * dt * (aa31 * Lt_n1[j] + aa32 * Lt_n2[j]);
            f_n3[j] = f_u(u_n3[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n3[j] = -dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n4[j] = u41 * u_n0[j] + u42 * Tu1[j] + dt * (a41 * L_n1[j] + a42 * L_n2[j] + a43 * L_n3[j]) +
                      dt * dt * (aa41 * Lt_n1[j] + aa42 * Lt_n2[j] + aa43 * Lt_n3[j]);

            f_n4[j] = f_u(u_n4[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n4[j] = -dfdx(f_n4, j);
            Lt_n4[j] = dfdxx(f_n4, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j] + b14 * L_n4[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j] + bb14 * Lt_n4[j]);

            Tu1[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j] + b24 * L_n4[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j] + bb24 * Lt_n4[j]);
        }
    }

    void compt_SGLM5_line()
    {
        int i, j;
        double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, aa21, aa31, aa32, aa41, aa42, aa43, aa51, aa52, aa53, aa54,
            b11, b12, b13, b14, b15, b21, b22, b23, b24, b25;
        double bb11, bb12, bb13, bb14, bb15, bb21, bb22, bb23, bb24, bb25, u11, u12, u21, u22, u31, u32, u41, u42, u51, u52, v;

        a21 = 0.061516399923018, a31 = 0.043600570185111, a32 = 0.444527820608227;
        a41 = 0.165875098133232, a42 = 0.383629364090464, a43 = 0.136281919092130;
        a51 = 0.242661894606699, a52 = 0.229098986759035, a53 = 0.081385974328652;
        a54 = 0.424475494404930, aa21 = 0, aa31 = 0.073569636769899;
        aa32 = 0.021058600222758, aa41 = 0.041456427578800, aa42 = 0.011866503266506;
        aa43 = 0.093912282949104, aa51 = 0.075476164472320, aa52 = 0.007086537500003;
        aa53 = 0.056083321251691, aa54 = 0, b11 = 0.187235442217099;
        b12 = 0.322862243971992, b13 = 0.248777699245909, b14 = 0.220697796096535;
        b15 = 0.002570885968929, b21 = 0.188276056945535, b22 = 0.322742421729971;
        b23 = 0.252236989947234, b24 = 0.239243854132961, b25 = 0.012555275639019;
        bb11 = 0.060195951745607, bb12 = 0.011759266769746, bb13 = 0.052645476899545;
        bb14 = 0.063174241824732, bb15 = 0.003372973897714, bb21 = 0.059952119791420;
        bb22 = 0.011801850994192, bb23 = 0.051926886958876, bb24 = 0.054650157061697;
        bb25 = 0.021527446840733, u11 = 0, u12 = 1;
        u21 = 0.913453478890676, u22 = 0.086546521109324, u31 = 0.647422355128817;
        u32 = 0.352577644871183, u41 = 0.535921231366958, u42 = 0.464078768633042;
        u51 = 0.320045915619410, u52 = 0.679954084380590, v = 0.542559843744499;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            u_n1[j] = u11 * u_n0[j] + u12 * Tu1[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u21 * u_n0[j] + u22 * Tu1[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n3[j] = u31 * u_n0[j] + u32 * Tu1[j] + dt * (a31 * L_n1[j] + a32 * L_n2[j]) + dt * dt * (aa31 * Lt_n1[j] + aa32 * Lt_n2[j]);
            f_n3[j] = f_u(u_n3[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n3[j] = -dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n4[j] = u41 * u_n0[j] + u42 * Tu1[j] + dt * (a41 * L_n1[j] + a42 * L_n2[j] + a43 * L_n3[j]) +
                      dt * dt * (aa41 * Lt_n1[j] + aa42 * Lt_n2[j] + aa43 * Lt_n3[j]);

            f_n4[j] = f_u(u_n4[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n4[j] = -dfdx(f_n4, j);
            Lt_n4[j] = dfdxx(f_n4, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n5[j] = u51 * u_n0[j] + u52 * Tu1[j] + dt * (a51 * L_n1[j] + a52 * L_n2[j] + a53 * L_n3[j] + a54 * L_n4[j]) +
                      dt * dt * (aa51 * Lt_n1[j] + aa52 * Lt_n2[j] + aa53 * Lt_n3[j] + aa54 * Lt_n4[j]);

            f_n5[j] = f_u(u_n5[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n5, j);
            Lt_n0[j] = dfdxx(f_n5, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j] + b14 * L_n4[j] + b15 * L_n0[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j] + bb14 * Lt_n4[j] + bb15 * Lt_n0[j]);

            Tu1[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j] + b24 * L_n4[j] + b25 * L_n0[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j] + bb24 * Lt_n4[j] + bb25 * Lt_n0[j]);
        }
    }

    void compt_SGLM6_line()
    {
        int i, j;
        double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, aa21, aa31, aa32, aa41, aa42, aa43, aa51, aa52, aa53, aa54,
            a61, a62, a63, a64, a65, aa61, aa62, aa63, aa64, aa65, b11, b12, b13, b14, b15, b16, b21, b22, b23, b24, b25, b26;
        double bb11, bb12, bb13, bb14, bb15, bb16, bb21, bb22, bb23, bb24, bb25, bb26, u11, u12, u21, u22, u31, u32, u41, u42, u51, u52, u61, u62, v;

        a21 = 0.195519007044922, a31 = 0.182698834641781, a32 = 0.224355689307196;
        a41 = 0.200798240963995, a42 = 0.202087232337342, a43 = 0.214237056626004;
        a51 = 0.180742417657866, a52 = 0.206337715995125, a53 = 0.311032694918027;
        a54 = 0, a61 = 0.179400758954363, a62 = 0.231484391689903;
        a63 = 0.258646300437170, a64 = 0, a65 = 0.209570365940087;
        aa21 = 0.028978386926393, aa31 = 0.027078275413061, aa32 = 0.027503486292418;
        aa41 = 0.024390617200713, aa42 = 0.024773623711643, aa43 = 0.026908793160004;
        aa51 = 0.028375394881790, aa52 = 0.025294685242906, aa53 = 0.009443683013244;
        aa54 = 0.026792613259658, aa61 = 0.027844686173205, aa62 = 0.020895996142354;
        aa63 = 0.007469018252919, aa64 = 0.021190304375858, aa65 = 0.019434763451396;

        b11 = 0.162799979155626, b12 = 0.210064073093223, b13 = 0.234712565126676;
        b14 = 0., b15 = 0.190177853235071, b16 = 0.043897073679799,

        b21 = 0.186288536307596, b22 = 0.201414279006715, b23 = 0.299976348757000;
        b24 = 0., b25 = 0.372215734846917, b26 = 0.015729879093722;

        bb11 = 0.028681836912939, bb12 = 0.018962393226422, bb13 = 0.006777875539520;
        bb14 = 0.019229467761441, bb15 = 0.017636375137000, bb16 = 0.060702028317952,

        bb21 = 0.027537753865251, bb22 = 0.024389959300741, bb23 = 0.009092524569629,
        bb24 = 0.025796343863547, bb25 = 0.000782358307294, bb26 = 0.0;

        u11 = 0.382332042901846, u12 = 0.617667957098154, u21 = 0.336433540760703;
        u22 = 0.663566459239297, u31 = 0.314373608788154, u32 = 0.685626391211846;
        u41 = 0.325150848279669, u42 = 0.674849151720331, u51 = 0.321246382913682;
        u52 = 0.678753617086318, u61 = 0.337338926920657, u62 = 0.662661073079343;
        v = 0.676780216227852;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            u_n1[j] = u11 * u_n0[j] + u12 * Tu1[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u21 * u_n0[j] + u22 * Tu1[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n3[j] = u31 * u_n0[j] + u32 * Tu1[j] + dt * (a31 * L_n1[j] + a32 * L_n2[j]) + dt * dt * (aa31 * Lt_n1[j] + aa32 * Lt_n2[j]);
            f_n3[j] = f_u(u_n3[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n3[j] = -dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n4[j] = u41 * u_n0[j] + u42 * Tu1[j] + dt * (a41 * L_n1[j] + a42 * L_n2[j] + a43 * L_n3[j]) +
                      dt * dt * (aa41 * Lt_n1[j] + aa42 * Lt_n2[j] + aa43 * Lt_n3[j]);

            f_n4[j] = f_u(u_n4[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n4[j] = -dfdx(f_n4, j);
            Lt_n4[j] = dfdxx(f_n4, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n5[j] = u51 * u_n0[j] + u52 * Tu1[j] + dt * (a51 * L_n1[j] + a52 * L_n2[j] + a53 * L_n3[j] + a54 * L_n4[j]) +
                      dt * dt * (aa51 * Lt_n1[j] + aa52 * Lt_n2[j] + aa53 * Lt_n3[j] + aa54 * Lt_n4[j]);

            f_n5[j] = f_u(u_n5[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n5, j);
            Lt_n0[j] = dfdxx(f_n5, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n6[j] = u61 * u_n0[j] + u62 * Tu1[j] + dt * (a61 * L_n1[j] + a62 * L_n2[j] + a63 * L_n3[j] + a64 * L_n4[j] + a65 * L_n0[j]) +
                      dt * dt * (aa61 * Lt_n1[j] + aa62 * Lt_n2[j] + aa63 * Lt_n3[j] + aa64 * Lt_n4[j] + aa65 * Lt_n0[j]);
            f_n6[j] = f_u(u_n6[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            u_n7[j] = -dfdx(f_n6, j);
            f_n7[j] = dfdxx(f_n6, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j] + b14 * L_n4[j] + b15 * L_n0[j] + b16 * u_n7[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j] + bb14 * Lt_n4[j] + bb15 * Lt_n0[j] + bb16 * f_n7[j]);

            Tu1[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j] + b24 * L_n4[j] + b25 * L_n0[j] + b26 * u_n7[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j] + bb24 * Lt_n4[j] + bb25 * Lt_n0[j] + bb26 * f_n7[j]);
        }
    }

    void compt_SGLM3_line_int()
    {
        int i, j;
        double a21, a31, a32, aa21, aa31, aa32, b11, b12, b13, b21, b22, b23, bb11, bb12, bb13, bb21, bb22, bb23, u11, u12, u21, u22, u31, u32, v;

        a21 = 0.367036327407744, a31 = 0.305289583777252, a32 = 0.594276924770852;
        aa21 = 0.087361998687504, aa31 = 0.072665036743429, aa32 = 0;
        b11 = 0.378319518763646, b12 = 0.539413445719846, b13 = 0;
        b21 = 0.224801836087141, b22 = 0.437599416854594, b23 = 0.526106658638147;
        bb11 = 0.065956620927602, bb12 = 0, bb13 = 0.156810544725088;
        bb21 = 0.120798754439142, bb22 = 0, bb23 = 0;
        u11 = 0.883562994058718, u12 = 0.116437005941282, u21 = 0.756325839981176;
        u22 = 0.243674160018824, u31 = 0.629088685903634, u32 = 0.370911314096366;
        v = 0.303820705714043;
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            u_n1[j] = u_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n3[j] = u_n0[j] + dt * (a31 * L_n1[j] + a32 * L_n2[j]) + dt * dt * (aa31 * Lt_n1[j] + aa32 * Lt_n2[j]);
            f_n3[j] = f_u(u_n3[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n3[j] = -dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = u_n0[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j]);

            Tu1[j] = u_n0[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j]);
        }
    }

    void compt_SGLM5_line_int()
    {
        int i, j;
        double a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, aa21, aa31, aa32, aa41, aa42, aa43, aa51, aa52, aa53, aa54,
            b11, b12, b13, b14, b15, b21, b22, b23, b24, b25;
        double bb11, bb12, bb13, bb14, bb15, bb21, bb22, bb23, bb24, bb25, u11, u12, u21, u22, u31, u32, u41, u42, u51, u52, v;

        a21 = 0.061516399923018, a31 = 0.043600570185111, a32 = 0.444527820608227;
        a41 = 0.165875098133232, a42 = 0.383629364090464, a43 = 0.136281919092130;
        a51 = 0.242661894606699, a52 = 0.229098986759035, a53 = 0.081385974328652;
        a54 = 0.424475494404930, aa21 = 0, aa31 = 0.073569636769899;
        aa32 = 0.021058600222758, aa41 = 0.041456427578800, aa42 = 0.011866503266506;
        aa43 = 0.093912282949104, aa51 = 0.075476164472320, aa52 = 0.007086537500003;
        aa53 = 0.056083321251691, aa54 = 0, b11 = 0.187235442217099;
        b12 = 0.322862243971992, b13 = 0.248777699245909, b14 = 0.220697796096535;
        b15 = 0.002570885968929, b21 = 0.188276056945535, b22 = 0.322742421729971;
        b23 = 0.252236989947234, b24 = 0.239243854132961, b25 = 0.012555275639019;
        bb11 = 0.060195951745607, bb12 = 0.011759266769746, bb13 = 0.052645476899545;
        bb14 = 0.063174241824732, bb15 = 0.003372973897714, bb21 = 0.059952119791420;
        bb22 = 0.011801850994192, bb23 = 0.051926886958876, bb24 = 0.054650157061697;
        bb25 = 0.021527446840733, u11 = 0, u12 = 1;
        u21 = 0.913453478890676, u22 = 0.086546521109324, u31 = 0.647422355128817;
        u32 = 0.352577644871183, u41 = 0.535921231366958, u42 = 0.464078768633042;
        u51 = 0.320045915619410, u52 = 0.679954084380590, v = 0.542559843744499;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            u_n1[j] = u_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n3[j] = u_n0[j] + dt * (a31 * L_n1[j] + a32 * L_n2[j]) + dt * dt * (aa31 * Lt_n1[j] + aa32 * Lt_n2[j]);
            f_n3[j] = f_u(u_n3[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n3[j] = -dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n4[j] = u_n0[j] + dt * (a41 * L_n1[j] + a42 * L_n2[j] + a43 * L_n3[j]) +
                      dt * dt * (aa41 * Lt_n1[j] + aa42 * Lt_n2[j] + aa43 * Lt_n3[j]);

            f_n4[j] = f_u(u_n4[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n4[j] = -dfdx(f_n4, j);
            Lt_n4[j] = dfdxx(f_n4, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_n5[j] = u_n0[j] + dt * (a51 * L_n1[j] + a52 * L_n2[j] + a53 * L_n3[j] + a54 * L_n4[j]) +
                      dt * dt * (aa51 * Lt_n1[j] + aa52 * Lt_n2[j] + aa53 * Lt_n3[j] + aa54 * Lt_n4[j]);

            f_n5[j] = f_u(u_n5[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n5, j);
            Lt_n0[j] = dfdxx(f_n5, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = u_n0[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j] + b14 * L_n4[j] + b15 * L_n0[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j] + bb14 * Lt_n4[j] + bb15 * Lt_n0[j]);

            Tu1[j] = u_n0[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j] + b24 * L_n4[j] + b25 * L_n0[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j] + bb24 * Lt_n4[j] + bb25 * Lt_n0[j]);
        }
    }

    //------------------------------------------------
    void f_eq_u(double *ff, double *uu)
    {
        for (int i = 0; i < n1; i++)
        {

            *ff++ = f_u(*uu++);
        }
    }
    double f_u(double xx)
    {
        double y;
        // y = xx * xx / 2.0;//bergure
        y = xx; // line
        return y;
    }
    void Store_obj(double *ss, double *ff)
    {
        for (int ik = 0; ik < n1; ik++) //
        {
            *ss++ = *ff++;
        }
    }

    void Write_obj(int j)
    {
        if (j < -0.9)
        {
            // string title = "result_t_burgers_excat"; //+ to_string(t);
            string title = to_string(n);
            string Title = title + ".plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=x,u " << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (int j = 0; j < kt; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_obj[j][i];
                }

                ofs << endl;
            }
            ofs.close();
        }

        if (j < 0.5 && j > -0.8)
        {
            // string title = "result_t_burgers_excat"; //+ to_string(t);
            string title = to_string(n);
            string Title = "excat.txt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << "Title=" << title << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                ofs << "  ";
                ofs.precision(18);
                ofs << *(u_excat + i);
                ofs << endl;
            }
            ofs.close();
        }

        if (j > 1.5)
        {
            string title = "result_one_t_burgers_every-" + to_string(j);
            string Title = title + ".plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << "Title=" << title << endl;
            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (j = 0; j < onet_out; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_one_n[j][i];
                }

                ofs << endl;
            }
            ofs.close();
        }
    }
    ////-----------------------------------------
    void RK_compt_t1()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j = 0, i;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n0, nb + j); //-uu_t0[nb + j]; -dfdx
            df[nb + j] = -dfdx(f_n0, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_nn[nb + j] = f_u(u_nn[nb + j]);
            j++;
        }
    }

    void RK3_compt_t1()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j = 0, i;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n0, nb + j); //-uu_t0[nb + j]; -dfdx
            df[nb + j] = -dfdx(f_n0, nb + j);
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
    }
    void RK3_compt_t2()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t2=uu_t0*0.75 +uu_t1*0.25 + 0.25*dt*df
        int j = 0, i;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n1, nb + j); //-uu_t1[nb + j];
            df[nb + j] = -dfdx(f_n1, nb + j);
            u_n2[nb + j] = u_n0[nb + j] * 0.75 + u_n1[nb + j] * 0.25 + 0.25 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
    }
    void RK3_compt_tt()
    {
        // u=1.0/3.0*uu_t0 + 2.0/3.0 *uu_t2+ 2.0/3.0*dt*df
        int j = 0, i;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n2, nb + j);
            df[nb + j] = -dfdx(f_n2, nb + j);
            u_nn[nb + j] = u_n0[nb + j] * 1.0 / 3.0 + u_n2[nb + j] * 2.0 / 3.0 + 2.0 / 3.0 * dt * df[nb + j];
            j++;
        }
    }

    void SSPRK54_compt()
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.39175222700392 * dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.44437049406734 * u_n0[nb + j] + 0.55562950593266 * u_n1[nb + j] + 0.36841059262959 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.62010185138540 * u_n0[nb + j] + 0.37989814861460 * u_n2[nb + j] + 0.25189177424738 * dt * df[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            Ltt_n0[nb + j] = df[nb + j];
            u_n4[nb + j] = 0.17807995410773 * u_n0[nb + j] + 0.82192004589227 * u_n3[nb + j] + 0.54497475021237 * dt * df[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n4, nb + j);
            u_nn[nb + j] = 0.00683325884039 * u_n0[nb + j] + 0.51723167208978 * u_n2[nb + j] + 0.12759831133288 * u_n3[nb + j] +
                           0.34833675773694 * u_n4[nb + j] + 0.08460416338212 * dt * Ltt_n0[nb + j] + 0.22600748319395 * dt * df[nb + j];
            j++;
        }
    }

    void SSPRK54_compt_Shu()
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.391752226571890 * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.444370493651235 * u_n0[nb + j] + 0.555629506348765 * u_n1[nb + j] + 0.368410593050371 * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.620101851488403 * u_n0[nb + j] + 0.379898148511597 * u_n2[nb + j] + 0.251891774271694 * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = 0.178079954393132 * u_n0[nb + j] + 0.821920045606868 * u_n3[nb + j] + 0.544974750228521 * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n4, nb + j);
            u_nn[nb + j] = 0.517231671970585 * u_n2[nb + j] + 0.096059710526147 * u_n3[nb + j] + 0.063692468666290 * dt * L_n3[nb + j] +
                           0.386708617503269 * u_n4[nb + j] + 0.226007483236906 * dt * df[nb + j];

            j++;
        }
    }

    void RK65LS_compt() // On fifth order Runge-Kutta methods
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 6. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] - 0.21627570527696895373 * dt * L_n0[nb + j] + 0.54960903861030228706 * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 0.08482881411262012706 * dt * L_n0[nb + j] + 0.04162653285051884260 * dt * L_n1[nb + j] +
                           0.37354465303686103035 * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] - 0.08651098424575942561 * dt * L_n0[nb + j] + 0.37955562705964599292 * dt * L_n1[nb + j] +
                           0.01753570971622337002 * dt * L_n2[nb + j] + 0.35608631413655672933 * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = -dfdx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] - 0.12499755969423778621 * dt * L_n0[nb + j] + 0.72695084642093284094 * dt * L_n1[nb + j] - 0.38363171852137430626 * dt * L_n2[nb + j] +
                           0.29492374551818501854 * dt * L_n3[nb + j] + 0.32008801960982756632 * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n5, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 0.07892564703041163884 * dt * L_n0[nb + j] + 0.15537176484794180580 * dt * L_n1[nb + j] +
                           0.08925647030411638840 * dt * L_n2[nb + j] + 0.51074352969588361160 * dt * L_n3[nb + j] -
                           0.30537176484794180580 * dt * L_n4[nb + j] + 0.47107435296958836116 * dt * df[nb + j];
            j++;
        }
    }

    void RK65LS_compt2() // Numerical Methods for Ordinary Differential Equations
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 4. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] + 1. / 8. * dt * L_n0[nb + j] + 1. / 8. * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 0.0 * dt * L_n0[nb + j] + 0.0 * dt * L_n1[nb + j] +
                           1. / 2. * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] + 3. / 16. * dt * L_n0[nb + j] - 3. / 8. * dt * L_n1[nb + j] +
                           3. / 8. * dt * L_n2[nb + j] + 9. / 16. * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = -dfdx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] - 3. / 7. * dt * L_n0[nb + j] + 8. / 7. * dt * L_n1[nb + j] + 6. / 7. * dt * L_n2[nb + j] -
                           12. / 7. * dt * L_n3[nb + j] + 8. / 7. * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n5, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 7. / 90. * dt * L_n0[nb + j] + 0. * dt * L_n1[nb + j] +
                           16. / 45. * dt * L_n2[nb + j] + 2. / 15. * dt * L_n3[nb + j] +
                           16. / 45. * dt * L_n4[nb + j] + 7. / 90. * dt * df[nb + j];
            j++;
        }
    }

    //---------------------------------------------------------

    int Perimeter()
    {
        cout << " yes ";
        return 2 * t;
    }

    double dfdx(double *g, int i)
    {
        double y;
        // 1up
        y = (*(g + i) - *(g + i - 1)) / dx;

        return y;
    }

    double dfdxx(double *g, int i)
    {
        double y;
        // 2ord
        y = (*(g + i + 1) - *(g + i) * 2.0 + *(g + i - 1)) / dx / dx;
        // y = (*(g + i ) - *(g + i - 1) * 2.0 + *(g + i-2)) / dx / dx;
        return y;
    }

    void compt_CFL()
    {
        double *umax;
        // umax = max_element(u_n0 + 1, u_n0 + n1 - 1);
        // cout << *umax;
        // dt = CFL * dx / (*umax);
        // dt=0.02/1.0;
        dt = CFL * dx;
    }

    void maxTV()
    {
        double TV = 0.0;
        for (int i = 0; i < n1 - 1; i++) // 间断
        {
            TV = abs(u_nn[i + 1] - u_nn[i]) + TV;
        }
        TV = abs(TV - 2.0);
        // cout << " \n TV: " << TV;
        if (tt > 0.0)
        {
            if (TV > TVmax)
            {
                TVmax = TV;
            }
        }
    }

    void TVabs(int m)
    {
        double TV = 0.0;
        TV = TVmax;

        if (TV < 2.0E-16)
        {
            TV = 2.0E-16;
        }
        cout << "\n  CFL: " << CFL << ",  TV: " << TV << "   \n";
        string Title = "TV_CFL.plt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=CFL,TV " << endl;
            ofs.precision(4), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            ofs.precision(4), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
    }
};
////-----------------------------------------
void comput_tt(cfda &);
int comput_main(int, double, double, int);
int main()
{
    // int n, oder_n = 2, begin_cell = 1600;
    int n, oder_n = 50, begin_cell = 3200;
    double t = 100.0, CFL = 1.1; // 计算时间总长，CFL数  CFL = 0.66 1.45
    n = begin_cell + 1;
    for (int m = 1; m < oder_n; m++)
    {
        comput_main(n, t, CFL, m);
        CFL = CFL + 0.005;
    }
    // system("pause");
    return 2;
}
int comput_main(int n, double t, double CFL, int m)
{
    int nb = 6; // 边界网格结点数（有限差分）
    int n1 = n + 6 * 2;
    double main_df[n + 6 * 2], u_excat[n + 6 * 2];
    double Ltt_n0[n + 6 * 2], x[n + 6 * 2];                                                          // 声明变长数组
    double L_n0[n + 6 * 2], L_n1[n + 6 * 2], L_n2[n + 6 * 2], L_n3[n + 6 * 2], L_n4[n + 6 * 2];      // f=L
    double Lt_n0[n + 6 * 2], Lt_n1[n + 6 * 2], Lt_n2[n + 6 * 2], Lt_n3[n + 6 * 2], Lt_n4[n + 6 * 2]; // L_t
    double u_n0[n + 6 * 2], u_n1[n + 6 * 2], u_n2[n + 6 * 2], u_n3[n + 6 * 2], u_n4[n + 6 * 2];      // 推进步中间时刻
    double u_n5[n + 6 * 2], u_n6[n + 6 * 2], u_n7[n + 6 * 2], u_n8[n + 6 * 2], u_n9[n + 6 * 2];      // 推进步中间时刻
    double f_n0[n + 6 * 2], f_n1[n + 6 * 2], f_n2[n + 6 * 2], f_n3[n + 6 * 2], f_n4[n + 6 * 2];      // 推进步中间时刻
    double f_n5[n + 6 * 2], f_n6[n + 6 * 2], f_n7[n + 6 * 2], f_n8[n + 6 * 2], f_n9[n + 6 * 2];      // 推进步中间时刻
    double u_nn[n + 6 * 2], f_nn[n + 6 * 2], Lttx_n0[n + 6 * 2];                                     // 下一时刻（n+1)

    double TL_n0[n + 6 * 2], TL_n1[n + 6 * 2], TL_n2[n + 6 * 2], TL_n3[n + 6 * 2], TL_n4[n + 6 * 2];
    double TLt_n0[n + 6 * 2], TLt_n1[n + 6 * 2], TLt_n2[n + 6 * 2], TLt_n3[n + 6 * 2], TLt_n4[n + 6 * 2];
    double Tu1[n + 6 * 2], Tu2[n + 6 * 2];
    class cfda cff;
    cff.kt = 1;
    cff.A1 = 1 / 3.0;
    cff.A2 = 2 / 3.0;
    cff.t = t;           // 1.0 / cff.PI; //计算总时间
    cff.nx_begin = -0.7; // 网格右端点
    cff.nx_end = 1.3;    // 网格左端点
    cff.CFL = CFL, cff.TVmax = 0.0;

    cff.n = n, cff.nb = nb, cff.n1 = n1;
    cff.Ltt_n0 = Ltt_n0, cff.x = x;
    cff.u_n0 = u_n0, cff.u_n1 = u_n1, cff.u_n2 = u_n2, cff.u_n3 = u_n3, cff.u_n4 = u_n4;                     // 推进步中间时刻
    cff.u_n5 = u_n5, cff.u_n6 = u_n6, cff.u_n7 = u_n7, cff.u_n8 = u_n8, cff.u_n9 = u_n9;                     // 推进步中间时刻
    cff.f_n0 = f_n0, cff.f_n1 = f_n1, cff.f_n2 = f_n2, cff.f_n3 = f_n3, cff.f_n4 = f_n4;                     // 推进步中间时刻
    cff.f_n5 = f_n5, cff.f_n6 = f_n6, cff.f_n7 = f_n7, cff.f_n8 = f_n8, cff.f_n9 = f_n9;                     // 推进步中间时刻
                                                                                                             // 声明变长数组
    cff.u_nn = u_nn, cff.f_nn = f_nn, cff.Lttx_n0 = Lttx_n0;                                                 // 下一时刻（n+1)
    cff.L_n0 = L_n0, cff.L_n1 = L_n1, cff.L_n2 = L_n2, cff.L_n3 = L_n3, cff.L_n4 = L_n4;                     // f=L
    cff.Lt_n0 = Lt_n0, cff.Lt_n1 = Lt_n1, cff.Lt_n2 = Lt_n2, cff.Lt_n3 = Lt_n3, cff.Lt_n4 = Lt_n4;           // L_t
                                                                                                             // 存储一时刻（n-1)
    cff.TL_n0 = TL_n0, cff.TL_n1 = TL_n1, cff.TL_n2 = TL_n2, cff.TL_n3 = TL_n3, cff.TL_n4 = TL_n4;           // f=TL
    cff.TLt_n0 = TLt_n0, cff.TLt_n1 = TLt_n1, cff.TLt_n2 = TLt_n2, cff.TLt_n3 = TLt_n3, cff.TLt_n4 = TLt_n4; // TL_t

    cff.Tu1 = Tu1, cff.Tu2 = Tu2;
    cff.df = main_df, cff.u_excat = u_excat, cff.intc();
    for (cff.ij = 0; cff.ij < cff.kt; cff.ij++)
    {
        // cff.nt = 1000 * pow(2, cff.ij); //计算时间步数

        comput_tt(cff);

        cff.Store_obj(cff.put_obj[cff.ij], cff.u_nn); // 存储数据
        // cff.Write_obj(cff.ij + 1);                    //存储数据
    }
    cff.Write_obj(-1);
    cff.TVabs(m);
    return 2;
}
////-----------------------------------------
void comput_tt(cfda &cff1)
{
    int n_all, i3 = 0, tt_flag = 0;
    cff1.tt = 0, n_all = cff1.n1, cff1.border();
    cff1.f_eq_u(cff1.f_n0, cff1.u_n0), cff1.TVmax = 0.0;
    for (int i1 = 0; i1 < 500; i1++) // 500 200 50
    {
        cff1.compt_CFL();
        cff1.tt = cff1.tt + cff1.dt;
        if (i1 < 2)
        {
            if (i1 < 1)
                cff1.compt_2_stage_int(); //   cff1.compt_SGLM_int(); //  cff1.compt_2_stage_int();

            cff1.dt = cff1.dt * 0.5;

            cff1.compt23_LLttnn_line_SSP();
            // cff1.compt24_LLttnn_line_SSP(); // cff1.compt34_LLttnn_line_SSP();
            // cff1.compt35_LLttnn_line_SSP();

            // cff1.compt_SGLM5_line_int();
            // cff1.compt_SGLM3_line_int();

            // cff1.RK3_compt_t1(); // 3-RK
            // cff1.RK3_compt_t2(); // 3-RK
            // cff1.RK3_compt_tt(); // 3-RK

            // cff1.SSPRK54_compt_Shu();// cff1.SSPRK54_compt(); //
            // cff1.RK65LS_compt2();// cff1.RK65LS_compt();
        }
        else
        {
            cff1.compt_TDTS23_LLttnn_line();
            // cff1.compt_TDTS24_LLttnn_line();
            // cff1.compt_TDTS25_LLttnn_line();

            // cff1.compt_SGLM3_line();
            // cff1.compt_SGLM4_line();
            // cff1.compt_SGLM5_line();
            // cff1.compt_SGLM6_line();
            // cff1.compt_TDTS35_LLttnn_line();
        }

        //___________________________________________________________________
        if (cff1.tt > cff1.t)
        {
            cff1.f_eq_u(cff1.f_nn, cff1.u_nn);
            cout.precision(18); // 精度为18，正常为6
            cout << " \n comput time: " << cff1.tt - cff1.dt << " comput time n:   " << i1;
            break;
        }
        cff1.carr(cff1.u_n0, cff1.u_nn, n_all);
        cff1.f_eq_u(cff1.f_n0, cff1.u_nn);

        // cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1;

        cff1.border();
        cff1.maxTV();

        if (i3 * cff1.t / cff1.onet_out < cff1.tt + 0.00000001)
        {
            // cff1.Store_obj(cff1.put_one_n[i3], cff1.u_nn); //存储数据
            i3++;
        }
    }
}
