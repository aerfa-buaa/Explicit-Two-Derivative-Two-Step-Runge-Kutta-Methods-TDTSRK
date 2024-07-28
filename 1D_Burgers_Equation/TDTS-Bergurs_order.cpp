#include <math.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>
#include <iterator>
#include <valarray>
using namespace std;
// example:bergurs TDTS   增加CFL数输出误差   EXP2
class cfda
{

private:
    int w;

public:
    double PI = 3.141592653589793;
    // double PI = 3.1415926535897932385;

    int n;                             // 网格结点数（有限差分）
    int nb;                            // 边界网格结点数（有限差分）
    double t, tt, CFL;                 // 计算总时 间t,实际计算时间tt;
    int nt, kt, ij, ijj, onet_out = 2; // 计算时间步长;
    double nx_begin, nx_end;           // 网格长度;
    int n1, ii1;                       // 总的网格结点数（有限差分）
    double *Ltt_n0, *x;                // 声明变长数组

    double *L_n0, *L_n1, *L_n2, *L_n3, *L_n4;      // f=L
    double *Lt_n0, *Lt_n1, *Lt_n2, *Lt_n3, *Lt_n4; // L_t

    double *TL_n0, *TL_n1, *TL_n2, *TL_n3, *TL_n4;      // Tf=L
    double *TLt_n0, *TLt_n1, *TLt_n2, *TLt_n3, *TLt_n4; // TL_t

    double *u_n0, *u_n1, *u_n2, *u_n3, *u_n4, *u_n5, *u_n6, *u_n7, *u_n8, *u_n9; // 推进步中间时刻
    double *f_n0, *f_n1, *f_n2, *f_n3, *f_n4, *f_n5, *f_n6, *f_n7, *f_n8, *f_n9; // 推进步中间时刻

    double *Ltx_n0, *Ltx_n1, *Ltx_n2; // L_tx
    double *u_nn, *f_nn, *Lttx_n0;    // 下一时刻（n+1)

    double *Tu1, *Tu2, *fz, *ff;
    // double  u_excat[513 + 6 * 2];
    double *df, *u_excat; // du[513 + 6 * 2],df[513 + 6 * 2],
    ////-----------------------------------------
    double dx;
    double dt;
    double *dd; // ss
    double A1, A2, A3, A4, A5, B0, B1, B2, B3, C0, C1, C2, C3, D0, D1, D2, D3, maxuu, TVmax, k;
    double a21, a32, a31, aa12, aa21, aa31, aa32, w1, w2, v1, v2, v3, vv1, vv2, vv3, ww1, ww2;
    double put_obj[2][13000];   //
    double put_one_n[2][13000]; // onet_out = 10;
    ////-----------------------------------------
    void intc()
    {
        dx = (nx_end - nx_begin) / (1.0 * (n - 1));

        for (int i = 0; i < n1; i++)
        {
            x[i] = nx_begin + (i - nb) * dx;
            u_n0[i] = intc_fun(x[i]);
            u_excat[i] = excat_bergurs(x[i]);
        }

        // for (int i = 0; i < n1; i++) //间断
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     u_n0[i] = 0.0;
        //     if (x[i] >= -0.5 && x[i] <= 0.5)
        //     {
        //         u_n0[i] = 1.0;
        //     }
        // }

        //  Write_obj(0);
    }
    double intc_fun(double xxf)
    {
        double yyf;
        // yyf = (sin(2 * PI * xxf) / 2.0 + 0.5) * 1.0; //  / PI bergure
        //  yyf = 0.1 * sin(PI * xxf) + 0.1; //bergure
        yyf = 0.5 * sin(PI * xxf) + 0.5; // bergure
        return yyf;
    }

    double excat_bergurs(double xxf)
    {
        double sigma = 1E-16, vak, val = 0.3; //
        int SMax = 2000, iMax = 0;
        while (iMax < SMax)
        {
            vak = intc_fun(xxf - val * t);
            if (fabs(val - vak) < sigma)
            {
                // cout << iMax << endl;
                break;
            }
            val = vak;
            iMax++;
        }
        return val;
    }

    // //-----------------------------------------
    void border()
    {
        borderfun(u_n0), borderfun(u_n1), borderfun(u_n2), borderfun(u_n3);
        borderfun(u_n4), borderfun(u_n5), borderfun(u_n6), borderfun(u_nn);
        borderfun(f_n0), borderfun(f_n1), borderfun(f_n2), borderfun(f_n3);
        borderfun(f_n4), borderfun(f_n5), borderfun(f_n6), borderfun(f_nn);

        borderfun(L_n0), borderfun(L_n1), borderfun(L_n2), borderfun(L_n3);
        borderfun(Ltx_n0), borderfun(Ltx_n1), borderfun(Ltx_n2);
    }

    void borderfun(double *ffbc)
    {
        for (int i = 0; i < nb; i++)
        {
            *(ffbc + i) = *(ffbc + n1 - 2 * nb + i - 1);  // 周期边条
            *(ffbc + n1 - nb + i) = *(ffbc + nb + i + 1); // 周期边条

            // *(ffbc + i) = 0;
            //*(ffbc + n1 - nb + i) = 0;
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
    void RK33_compt()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -df_dx(f_n0, nb + j);
            // df[nb + j] = dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -df_dx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] * 0.75 + 0.25 * u_n1[nb + j] + 0.25 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -df_dx(f_n2, nb + j);
            u_nn[nb + j] = u_n0[nb + j] * 1.0 / 3.0 + u_n2[nb + j] * 2.0 / 3.0 + 2.0 / 3.0 * dt * df[nb + j];
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
            L_n0[nb + j] = -df_dx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.391752226571890 * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -df_dx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.444370493651235 * u_n0[nb + j] + 0.555629506348765 * u_n1[nb + j] + 0.368410593050371 * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -df_dx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.620101851488403 * u_n0[nb + j] + 0.379898148511597 * u_n2[nb + j] + 0.251891774271694 * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -df_dx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = 0.178079954393132 * u_n0[nb + j] + 0.821920045606868 * u_n3[nb + j] + 0.544974750228521 * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -df_dx(f_n4, nb + j);
            u_nn[nb + j] = 0.517231671970585 * u_n2[nb + j] + 0.096059710526147 * u_n3[nb + j] + 0.063692468666290 * dt * L_n3[nb + j] +
                           0.386708617503269 * u_n4[nb + j] + 0.226007483236906 * dt * df[nb + j];

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
            L_n0[nb + j] = -df_dx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 4. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -df_dx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] + 1. / 8. * dt * L_n0[nb + j] + 1. / 8. * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -df_dx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 0.0 * dt * L_n0[nb + j] + 0.0 * dt * L_n1[nb + j] +
                           1. / 2. * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -df_dx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] + 3. / 16. * dt * L_n0[nb + j] - 3. / 8. * dt * L_n1[nb + j] +
                           3. / 8. * dt * L_n2[nb + j] + 9. / 16. * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = -df_dx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] - 3. / 7. * dt * L_n0[nb + j] + 8. / 7. * dt * L_n1[nb + j] + 6. / 7. * dt * L_n2[nb + j] -
                           12. / 7. * dt * L_n3[nb + j] + 8. / 7. * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -df_dx(f_n5, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 7. / 90. * dt * L_n0[nb + j] + 0. * dt * L_n1[nb + j] +
                           16. / 45. * dt * L_n2[nb + j] + 2. / 15. * dt * L_n3[nb + j] +
                           16. / 45. * dt * L_n4[nb + j] + 7. / 90. * dt * df[nb + j];
            j++;
        }
    }

    //------------------------------------------
    void compt_2_stage_int()
    {
        int i, j;
        a21 = 0.525; // 0.525  0.680; //0.765;   // 1. / 210. * (93. + sqrt(2139.));
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * df_dx(f_n0, j);
            Ltx_n0[j] = -L_n0[j] * u_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = df_dx(Ltx_n0, j);
            Lt_n0[j] = dfdxcent(Ltx_n0, j);

            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * df_dx(f_n1, j);
            Ltx_n1[j] = -L_n1[j] * u_n1[j];
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            // Lt_n1[j] = df_dx(Ltx_n1, j);
            Lt_n1[j] = dfdxcent(Ltx_n1, j);

            TL_n1[j] = TL_n3[j], TL_n0[j] = TL_n2[j];
            TLt_n1[j] = TLt_n3[j], TLt_n0[j] = TLt_n2[j];

            TL_n3[j] = L_n1[j], TL_n2[j] = L_n0[j];
            TLt_n3[j] = Lt_n1[j], TLt_n2[j] = Lt_n0[j];
        }
    }

    void compt23_LLtt_SSP()
    {
        int i, j;
        a21 = 0.594223212099088, v2 = 0.306027487008159;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

        i = n1 - nb * 2;
        border();
        borderfun(f_n0);
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * df_dx(f_n0, j);
            Ltx_n0[j] = -L_n0[j] * u_n0[j];
        }

        i = n1 - nb * 2;
        border();
        borderfun(Ltx_n0);
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = df_dx(Ltx_n0, j);
            Lt_n0[j] = dfdxcent(Ltx_n0, j);

            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 * 0.5 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        borderfun(f_n1);
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * df_dx(f_n1, j);
            Ltx_n1[j] = -L_n1[j] * u_n1[j];
        }

        i = n1 - nb * 2;
        border();
        borderfun(Ltx_n1);
        while (i-- > 0)
        {
            j = i + nb;
            // Lt_n1[j] = df_dx(Ltx_n1, j);
            Lt_n1[j] = dfdxcent(Ltx_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j]);

            // u_nn[j] = u_n0[j] + dt * L_n0[j] + 0.5 * dt * dt * Lt_n0[j];
        }
    }

    void compt24_LLtt_SSP()
    {
        int i, j;
        A1 = 0.5, B0 = 1., B1 = 0., C0 = 1. / 3., C1 = 2. / 3.;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * df_dx(f_n0, j);
            Ltx_n0[j] = -L_n0[j] * u_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = df_dx(Ltx_n0, j);
            Lt_n0[j] = dfdxcent(Ltx_n0, j);

            u_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A1 * A1 * 0.5 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * df_dx(f_n1, j);
            Ltx_n1[j] = -L_n1[j] * u_n1[j];
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            // Lt_n1[j] = df_dx(Ltx_n1, j);
            Lt_n1[j] = dfdxcent(Ltx_n1, j);

            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j]);
            // u_nn[j] = u_n0[j] + dt * (L_n0[j]);
            // u_nn[j] = u_n0[j] + dt * L_n0[j] + 0.5 * dt * dt * Lt_n0[j];
        }
    }

    void compt35_LLtt_SSP()
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
            L_n0[j] = -1. * df_dx(f_n0, j);
            Ltx_n0[j] = -L_n0[j] * u_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = df_dx(Ltx_n0, j);
            Lt_n0[j] = dfdxcent(Ltx_n0, j);

            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n1[j] = -1. * df_dx(f_n1, j);
            Ltx_n1[j] = -L_n1[j] * u_n1[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = df_dx(Ltx_n0, j);
            Lt_n1[j] = dfdxcent(Ltx_n1, j);

            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -1. * df_dx(f_n2, j);
            Ltx_n2[j] = -L_n2[j] * u_n2[j];
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            // Lt_n1[j] = df_dx(Ltx_n1, j);
            Lt_n2[j] = dfdxcent(Ltx_n2, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + v3 * L_n2[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j]);

            // u_nn[j] = u_n0[j] + dt * (L_n0[j]);
        }
    }

    //------------------------------------------------

    void compt_TDTS23_LLttnn()
    {
        int i, j;
        // a21 = 0.532, v2 = 0.389;
        // ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.1;
        a21 = 0.525, v2 = 0.448;
        ww1 = 0.0, ww2 = 0.0, w1 = 0.04, w2 = 0.0;

        v1 = 1. - v2 - w1 - w2;
        vv1 = -((1 - 3 * w1 - 3 * w2 + 6 * ww1 + 3 * a21 * (-1 - 2 * w1 + a21 * (v2 + w2) + 2 * ww1) + 6 * ww2) / (6. * a21));
        vv2 = (1 - 3 * w1 - 3 * w2 + 6 * ww1 + 6 * ww2 - 3 * a21 * (-2 * w2 + a21 * (v2 + w2) + 2 * ww2)) / (6. * a21);

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * df_dx(f_n0, j);
            Ltx_n0[j] = -L_n0[j] * u_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = df_dx(Ltx_n0, j);
            Lt_n0[j] = dfdxcent(Ltx_n0, j);

            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * df_dx(f_n1, j);
            Ltx_n1[j] = -L_n1[j] * u_n1[j];
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            // Lt_n1[j] = df_dx(Ltx_n1, j);
            Lt_n1[j] = dfdxcent(Ltx_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            // u_nn[j] = u_n0[j] + dt * (L_n0[j]);
            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }
    void compt_TDTS24_LLttnn()
    {
        int i, j;

        // w2 = -0.1, ww2 = 0.0, a21 = 0.468, v1 = 1.0;
        // w2 = -0.15, ww2 = 0.15, a21 = 0.698, v1 = 0.85; // C=0.884
        a21 = 0.680, v1 = 0.908;
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
            L_n0[j] = -1. * df_dx(f_n0, j);
            Ltx_n0[j] = -L_n0[j] * u_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = df_dx(Ltx_n0, j);
            Lt_n0[j] = dfdxcent(Ltx_n0, j);

            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * df_dx(f_n1, j);
            Ltx_n1[j] = -L_n1[j] * u_n1[j];
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            // Lt_n1[j] = df_dx(Ltx_n1, j);
            Lt_n1[j] = dfdxcent(Ltx_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            // u_nn[j] = u_n0[j] + dt * (L_n0[j]);
            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }
    void compt_TDTS25_LLttnn()
    {
        int i, j;

        a21 = 0.765; // 1. / 210. * (93. + sqrt(2139.));
        vv1 = (-31. + 6. * a21 * (-31. + 85. * a21)) / (360. * pow(a21, 2));
        vv2 = 31. / (360. * pow(a21, 2));
        v1 = -(1. / 2.) + 31. / (30. * a21), v2 = 0.;
        w1 = 3. / 2. - 31. / (30. * a21), w2 = 0.;
        ww1 = (31. + 6. * a21 * (-31. + 35. * a21)) / (360. * pow(a21, 2));
        ww2 = -(31. / (360. * pow(a21, 2)));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * df_dx(f_n0, j);
            Ltx_n0[j] = -L_n0[j] * u_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            // Lt_n0[j] = df_dx(Ltx_n0, j);
            Lt_n0[j] = dfdxcent(Ltx_n0, j);

            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * df_dx(f_n1, j);
            Ltx_n1[j] = -L_n1[j] * u_n1[j];
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            // Lt_n1[j] = df_dx(Ltx_n1, j);
            Lt_n1[j] = dfdxcent(Ltx_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            // u_nn[j] = u_n0[j] + dt * (L_n0[j]);
            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
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
        y = xx * xx / 2.0; // bergure
        // y = xx; //line
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
    int Perimeter()
    {
        cout << " yes ";
        return 2 * t;
    }

    void Splitting(double *uuu)
    {
        double lmax = 0.0, lmax0;
        for (int i = 0; i < n1; i++) // comput L(n0)
        {
            lmax0 = abs(uuu[i]); // 当地最大特征值
            if (lmax0 > lmax)
            {
                lmax = lmax0;
            }
        }
        //  lmax = lmax *1.1;
        for (int i = 0; i < n1; i++) // comput L(n0)
        {
            fz[i] = 0.5 * (0.5 * uuu[i] * uuu[i] + lmax * uuu[i]); // 正通量
            ff[i] = 0.5 * (0.5 * uuu[i] * uuu[i] - lmax * uuu[i]); // 负通量
        }

        borderfun(ff);
        borderfun(fz);
    }

    double df_dx(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i) - *(g + i - 1)) / dx;

        // 8up
        // y = (-*(g + i - 3) * 5. + *(g + i - 2) * 60. - *(g + i - 1) * 420. - *(g + i) * 378 + *(g + i + 1) * 1050. - *(g + i + 2) * 420. + *(g + i + 3) * 140. - *(g + i + 4) * 30. + *(g + i + 5) * 3.) / 840. / dx;

        y = -(-*(g + i + 3) * 5. + *(g + i + 2) * 60. - *(g + i + 1) * 420. - *(g + i) * 378 + *(g + i - 1) * 1050. - *(g + i - 2) * 420. + *(g + i - 3) * 140. - *(g + i - 4) * 30. + *(g + i - 5) * 3.) / 840. / dx;

        /*  //weno_5-js  ture
        double ep = 1.E-15, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i - 2) - 2. * (*(g + i - 1)) + *(g + i)), 2) + 0.25 * pow((*(g + i - 2) - 4. * (*(g + i - 1)) + 3. * (*(g + i))), 2);
            sn1 = 13. / 12. * pow((*(g + i - 1) - 2. * (*(g + i)) + *(g + i + 1)), 2) + 0.25 * pow((*(g + i - 1) - *(g + i + 1)), 2);
            sn2 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i + 1)) + *(g + i + 2)), 2) + 0.25 * pow((3. * (*(g + i)) - 4. * (*(g + i + 1)) + *(g + i + 2)), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 2. / 6. * (*(g + i - 2)) - 7. / 6. * (*(g + i - 1)) + 11. / 6. * (*(g + i));
            q13 = -1. / 6. * (*(g + i - 1)) + 5. / 6. * (*(g + i)) + 2. / 6. * (*(g + i + 1));
            q23 = 2. / 6. * (*(g + i)) + 5. / 6. * (*(g + i + 1)) - 1. / 6. * (*(g + i + 2));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
        */

        /*  //weno_9
        double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4, vz[2];
        double uz1, uz2, uz3, uz4, uz5;

        double beta0, beta1, beta2, beta3, beta4, omt0, omt1, omt2, omt3, omt4, omts;

        const double g0 = 1.0 / 126.0, g1 = 10. / 63., g2 = 10. / 21., g3 = 20. / 63., g4 = 5. / 126., eps = 1e-6;

        for (int j = 0; j < 2; j++)
        {
            uim4 = *(g + i - 4), uim3 = *(g + i - 3), uim2 = *(g + i - 2), uim1 = *(g + i - 1);
            ui = *(g + i), uip1 = *(g + i + 1), uip2 = *(g + i + 2), uip3 = *(g + i + 3), uip4 = *(g + i + 4);

            beta0 = uim4 * (22658. * uim4 - 208501. * uim3 + 364863. * uim2 - 288007. * uim1 + 86329. * ui) +
                    uim3 * (482963. * uim3 - 1704396. * uim2 + 1358458. * uim1 - 411487. * ui) +
                    uim2 * (1521393. * uim2 - 2462076. * uim1 + 758823. * ui) +
                    uim1 * (1020563. * uim1 - 649501. * ui) + 107918. * pow(ui, 2);
            beta1 = uim3 * (6908. * uim3 - 60871. * uim2 + 99213. * uim1 - 70237. * ui + 18079. * uip1) +
                    uim2 * (138563. * uim2 - 464976. * uim1 + 337018. * ui - 88297. * uip1) +
                    uim1 * (406293. * uim1 - 611976 * ui + 165153. * uip1) +
                    ui * (242723. * ui - 140251. * uip1) + 22658. * pow(uip1, 2);
            beta2 = uim2 * (6908. * uim2 - 51001. * uim1 + 67923. * ui - 38947. * uip1 + 8209. * uip2) +
                    uim1 * (104963. * uim1 - 299076. * ui + 179098. * uip1 - 38947. * uip2) +
                    ui * (231153. * ui - 299076. * uip1 + 67923. * uip2) +
                    uip1 * (104963. * uip1 - 51001. * uip2) + 6908. * pow(uip2, 2);
            beta3 = uim1 * (22658. * uim1 - 140251. * ui + 165153. * uip1 - 88297. * uip2 + 18079. * uip3) +
                    ui * (242723. * ui - 611976. * uip1 + 337018. * uip2 - 70237. * uip3) +
                    uip1 * (406293. * uip1 - 464976. * uip2 + 99213. * uip3) +
                    uip2 * (138563. * uip2 - 60871. * uip3) + 6908. * pow(uip3, 2);
            beta4 = ui * (107918. * ui - 649501. * uip1 + 758823. * uip2 - 411487. * uip3 + 86329. * uip4) +
                    uip1 * (1020563. * uip1 - 2462076. * uip2 + 1358458. * uip3 - 288007 * uip4) +
                    uip2 * (1521393. * uip2 - 1704396. * uip3 + 364863. * uip4) +
                    uip3 * (482963. * uip3 - 208501. * uip4) + 22658. * pow(uip4, 2);

            uz1 = (1. / 5.) * uim4 + (-21. / 20.) * uim3 + (137. / 60.) * uim2 + (-163. / 60.) * uim1 + (137. / 60.) * ui;
            uz2 = (-1. / 20.) * uim3 + (17. / 60.) * uim2 + (-43. / 60.) * uim1 + (77. / 60.) * ui + (1. / 5.) * uip1;
            uz3 = (1. / 30.) * uim2 + (-13. / 60.) * uim1 + (47. / 60.) * ui + (9. / 20.) * uip1 + (-1. / 20.) * uip2;
            uz4 = (-1. / 20.) * uim1 + (9. / 20.) * ui + (47. / 60.) * uip1 + (-13. / 60.) * uip2 + (1. / 30.) * uip3;
            uz5 = (1. / 5.) * ui + (77. / 60.) * uip1 + (-43. / 60.) * uip2 + (17. / 60.) * uip3 + (-1. / 20.) * uip4;

            omt0 = g0 * pow(eps + beta0, -2), omt1 = g1 * pow(eps + beta1, -2), omt2 = g2 * pow(eps + beta2, -2), omt3 = g3 * pow(eps + beta3, -2);
            omt4 = g4 * pow(eps + beta4, -2), omts = omt0 + omt1 + omt2 + omt3 + omt4;
            vz[j] = (omt0 * uz1 + omt1 * uz2 + omt2 * uz3 + omt3 * uz4 + omt4 * uz5) / omts;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
          */

        return y;
    }

    double df_dx1(double *fz, double *ff, int ii)
    {
        double y, ss1, ss2;
        int i = ii;
        // 1up
        //  ss1 = (*(fz + i) - *(fz + i - 1)) / dx; //正通量
        //  ss2 = (*(ff + i + 1) - *(ff + i)) / dx; //负通量

        // /*  //weno_5-js  正通量
        double ep = 1.E-15, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
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
        // weno_5-js 负通量
        //  double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        //  double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        i = ii, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        for (int j = 0; j < 2; j++)
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
        ss2 = (vz[0] - vz[1]) / dx;
        // */

        y = ss1 + ss2;
        return y;
    }

    double dfdxcent(double *g, int i)
    {
        double y;
        // 1 up
        // y = (*(g + i) - *(g + i - 1)) / dx;

        // y = (*(g + i + 1) - *(g + i)) / dx;

        // 2 center
        // y = (*(g + i + 1) - *(g + i - 1)) / dx * 0.5;

        // 4 center
        // y = (*(g + i - 2) - *(g + i - 1) * 8. + *(g + i + 1) * 8. - *(g + i + 2)) / dx / 12.;

        // 8up_center
        y = (*(g + i - 4) * 3. - *(g + i - 3) * 32. + *(g + i - 2) * 168. - *(g + i - 1) * 672. - *(g + i) * 0.0 + *(g + i + 1) * 672. - *(g + i + 2) * 168. + *(g + i + 3) * 32. - *(g + i + 4) * 3.) / 840. / dx;

        // 8up
        // y = (-*(g + i - 3) * 5. + *(g + i - 2) * 60. - *(g + i - 1) * 420. - *(g + i) * 378 + *(g + i + 1) * 1050. - *(g + i + 2) * 420. + *(g + i + 3) * 140. - *(g + i + 4) * 30. + *(g + i + 5) * 3.) / 840. / dx;

        // y = -(-*(g + i + 3) * 5. + *(g + i + 2) * 60. - *(g + i + 1) * 420. - *(g + i) * 378 + *(g + i - 1) * 1050. - *(g + i - 2) * 420. + *(g + i - 3) * 140. - *(g + i - 4) * 30. + *(g + i - 5) * 3.) / 840. / dx;

        /*  //weno_5-js 负通量  ture https://doi.org/10.1007/s42967-019-0001-3
        double ep = 1.E-10, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i + 1) - 2. * (*(g + i + 2)) + *(g + i + 3)), 2) + 0.25 * pow((3. * (*(g + i + 1)) - 4. * (*(g + i + 2)) + *(g + i + 3)), 2);
            sn1 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i + 1)) + *(g + i + 2)), 2) + 0.25 * pow((*(g + i) - *(g + i + 2)), 2);
            sn2 = 13. / 12. * pow((*(g + i - 1) - 2. * (*(g + i)) + *(g + i + 1)), 2) + 0.25 * pow((*(g + i - 1) - 4. * (*(g + i)) + 3. * (*(g + i + 1))), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 11. / 6. * (*(g + i + 1)) - 7. / 6. * (*(g + i + 2)) + 2. / 6. * (*(g + i + 3));
            q13 = 2. / 6. * (*(g + i)) + 5. / 6. * (*(g + i + 1)) - 1. / 6. * (*(g + i + 2));
            q23 = -1. / 6. * (*(g + i - 1)) + 5. / 6. * (*(g + i)) + 2. / 6. * (*(g + i + 1));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
        */

        /*  //weno_5-js  true
        double ep = 1.E-15, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i - 2) - 2. * (*(g + i - 1)) + *(g + i)), 2) + 0.25 * pow((*(g + i - 2) - 4. * (*(g + i - 1)) + 3. * (*(g + i))), 2);
            sn1 = 13. / 12. * pow((*(g + i - 1) - 2. * (*(g + i)) + *(g + i + 1)), 2) + 0.25 * pow((*(g + i - 1) - *(g + i + 1)), 2);
            sn2 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i + 1)) + *(g + i + 2)), 2) + 0.25 * pow((3. * (*(g + i)) - 4. * (*(g + i + 1)) + *(g + i + 2)), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 2. / 6. * (*(g + i - 2)) - 7. / 6. * (*(g + i - 1)) + 11. / 6. * (*(g + i));
            q13 = -1. / 6. * (*(g + i - 1)) + 5. / 6. * (*(g + i)) + 2. / 6. * (*(g + i + 1));
            q23 = 2. / 6. * (*(g + i)) + 5. / 6. * (*(g + i + 1)) - 1. / 6. * (*(g + i + 2));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
        */

        /*  //weno_9
        double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4, vz[2];
        double uz1, uz2, uz3, uz4, uz5;

        double beta0, beta1, beta2, beta3, beta4, omt0, omt1, omt2, omt3, omt4, omts;

        const double g0 = 1.0 / 126.0, g1 = 10. / 63., g2 = 10. / 21., g3 = 20. / 63., g4 = 5. / 126., eps = 1e-6;

        for (int j = 0; j < 2; j++)
        {
            uim4 = *(g + i - 4), uim3 = *(g + i - 3), uim2 = *(g + i - 2), uim1 = *(g + i - 1);
            ui = *(g + i), uip1 = *(g + i + 1), uip2 = *(g + i + 2), uip3 = *(g + i + 3), uip4 = *(g + i + 4);

            beta0 = uim4 * (22658. * uim4 - 208501. * uim3 + 364863. * uim2 - 288007. * uim1 + 86329. * ui) +
                    uim3 * (482963. * uim3 - 1704396. * uim2 + 1358458. * uim1 - 411487. * ui) +
                    uim2 * (1521393. * uim2 - 2462076. * uim1 + 758823. * ui) +
                    uim1 * (1020563. * uim1 - 649501. * ui) + 107918. * pow(ui, 2);
            beta1 = uim3 * (6908. * uim3 - 60871. * uim2 + 99213. * uim1 - 70237. * ui + 18079. * uip1) +
                    uim2 * (138563. * uim2 - 464976. * uim1 + 337018. * ui - 88297. * uip1) +
                    uim1 * (406293. * uim1 - 611976 * ui + 165153. * uip1) +
                    ui * (242723. * ui - 140251. * uip1) + 22658. * pow(uip1, 2);
            beta2 = uim2 * (6908. * uim2 - 51001. * uim1 + 67923. * ui - 38947. * uip1 + 8209. * uip2) +
                    uim1 * (104963. * uim1 - 299076. * ui + 179098. * uip1 - 38947. * uip2) +
                    ui * (231153. * ui - 299076. * uip1 + 67923. * uip2) +
                    uip1 * (104963. * uip1 - 51001. * uip2) + 6908. * pow(uip2, 2);
            beta3 = uim1 * (22658. * uim1 - 140251. * ui + 165153. * uip1 - 88297. * uip2 + 18079. * uip3) +
                    ui * (242723. * ui - 611976. * uip1 + 337018. * uip2 - 70237. * uip3) +
                    uip1 * (406293. * uip1 - 464976. * uip2 + 99213. * uip3) +
                    uip2 * (138563. * uip2 - 60871. * uip3) + 6908. * pow(uip3, 2);
            beta4 = ui * (107918. * ui - 649501. * uip1 + 758823. * uip2 - 411487. * uip3 + 86329. * uip4) +
                    uip1 * (1020563. * uip1 - 2462076. * uip2 + 1358458. * uip3 - 288007 * uip4) +
                    uip2 * (1521393. * uip2 - 1704396. * uip3 + 364863. * uip4) +
                    uip3 * (482963. * uip3 - 208501. * uip4) + 22658. * pow(uip4, 2);

            uz1 = (1. / 5.) * uim4 + (-21. / 20.) * uim3 + (137. / 60.) * uim2 + (-163. / 60.) * uim1 + (137. / 60.) * ui;
            uz2 = (-1. / 20.) * uim3 + (17. / 60.) * uim2 + (-43. / 60.) * uim1 + (77. / 60.) * ui + (1. / 5.) * uip1;
            uz3 = (1. / 30.) * uim2 + (-13. / 60.) * uim1 + (47. / 60.) * ui + (9. / 20.) * uip1 + (-1. / 20.) * uip2;
            uz4 = (-1. / 20.) * uim1 + (9. / 20.) * ui + (47. / 60.) * uip1 + (-13. / 60.) * uip2 + (1. / 30.) * uip3;
            uz5 = (1. / 5.) * ui + (77. / 60.) * uip1 + (-43. / 60.) * uip2 + (17. / 60.) * uip3 + (-1. / 20.) * uip4;

            omt0 = g0 * pow(eps + beta0, -2), omt1 = g1 * pow(eps + beta1, -2), omt2 = g2 * pow(eps + beta2, -2), omt3 = g3 * pow(eps + beta3, -2);
            omt4 = g4 * pow(eps + beta4, -2), omts = omt0 + omt1 + omt2 + omt3 + omt4;
            vz[j] = (omt0 * uz1 + omt1 * uz2 + omt2 * uz3 + omt3 * uz4 + omt4 * uz5) / omts;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
          */

        return y;
    }

    double dfdx(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i + 1) - *(g + i)) / dx;
        // 1up
        y = -(*(g + i) - *(g + i - 1)) / dx;
        return y;
    }

    double dfdxx(double *g, int i)
    {
        double y;
        // 2ord
        y = (*(g + i + 1) - *(g + i) * 2.0 + *(g + i - 1)) / dx / dx;
        return y;
    }

    void compt_CFL()
    {
        double *umax;
        umax = max_element(u_n0 + 1, u_n0 + n1 - 1);
        // cout << *umax;
        //  dt = CFL * dx / (*umax);
        // dt=0.02/1.0;
        dt = CFL * dx;
        // dt = CFL * dx / 0.5;
    }

    void maxTV()
    {
        double TV = 0.0;
        for (int i = 0; i < n1 - 1; i++) // 间断int i = 0; i < n1 - 1; i++   int i = nb; i < n1 - nb; i++
        {
            TV = abs(u_nn[i + 1] - u_nn[i]) + TV;
        }
        TV = abs(TV - 2.0);
        // cout << " \n TV: " << TV;
        // if (tt > 0.0)
        if (ii1 > -10)
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
        TVmax = 0.0;

        for (int i = nb; i < n1 - nb; i++) // 间断  int i = nb; i < n1 - nb; i++
        {
            TV = TV + abs(u_nn[i] - *(u_excat + i));
        }
        TV = TV / (n1 - nb - nb);

        // for (int i = nb; i < n1 - nb; i++) // 间断  int i = nb; i < n1 - nb; i++
        // {
        //     TV = abs(u_nn[i] - *(u_excat + i));
        //     if (TV > TVmax)
        //     {
        //         TVmax = TV;
        //     }
        // }
        // TV = TVmax;

        cout << "\n  Grid: " << n << ",  TV: " << TV;
        string Title = "TV_CFL.plt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=CFL,TV " << endl;
            // ofs.precision(4), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            // ofs.precision(4), ofs << CFL, ofs << "  ";
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
    int n, oder_n = 6, begin_cell = 80;
    double t = 0.2, CFL = 0.8; // 计算时间总长，CFL数  CFL = 0.66 1.45
    // n = begin_cell + 1;
    for (int m = 1; m < oder_n; m++)
    {
        n = pow(2, m - 1) * begin_cell + 1;
        comput_main(n, t, CFL, m);
        // CFL = CFL + 0.01;
    }

    // cout << " \n All comput time: ";
    // printf("%d ms", clock()); //输出运行所费时间，单位毫秒ms
    // system("pause");

    return 2;
}

int comput_main(int n, double t, double CFL, int m)
{
    //
    int nb = 6; // 边界网格结点数（有限差分）
    int n1 = n + 6 * 2;
    double main_df[n + 6 * 2], u_excat[n + 6 * 2];
    double Ltt_n0[n + 6 * 2], x[n + 6 * 2]; // 声明变长数组

    double L_n0[n + 6 * 2], L_n1[n + 6 * 2], L_n2[n + 6 * 2], L_n3[n + 6 * 2], L_n4[n + 6 * 2];      // f=L
    double Lt_n0[n + 6 * 2], Lt_n1[n + 6 * 2], Lt_n2[n + 6 * 2], Lt_n3[n + 6 * 2], Lt_n4[n + 6 * 2]; // L_t
    double u_n0[n + 6 * 2], u_n1[n + 6 * 2], u_n2[n + 6 * 2], u_n3[n + 6 * 2], u_n4[n + 6 * 2];      // 推进步中间时刻
    double u_n5[n + 6 * 2], u_n6[n + 6 * 2], u_n7[n + 6 * 2], u_n8[n + 6 * 2], u_n9[n + 6 * 2];      // 推进步中间时刻
    double f_n0[n + 6 * 2], f_n1[n + 6 * 2], f_n2[n + 6 * 2], f_n3[n + 6 * 2], f_n4[n + 6 * 2];      // 推进步中间时刻
    double f_n5[n + 6 * 2], f_n6[n + 6 * 2], f_n7[n + 6 * 2], f_n8[n + 6 * 2], f_n9[n + 6 * 2];      // 推进步中间时刻
    double u_nn[n + 6 * 2], f_nn[n + 6 * 2], Lttx_n0[n + 6 * 2];                                     // 下一时刻（n+1)

    double Ltx_n0[n + 6 * 2], Ltx_n1[n + 6 * 2], Ltx_n2[n + 6 * 2]; // L_tx
    double TL_n0[n + 6 * 2], TL_n1[n + 6 * 2], TL_n2[n + 6 * 2], TL_n3[n + 6 * 2], TL_n4[n + 6 * 2];
    double TLt_n0[n + 6 * 2], TLt_n1[n + 6 * 2], TLt_n2[n + 6 * 2], TLt_n3[n + 6 * 2], TLt_n4[n + 6 * 2];
    double Tu1[n + 6 * 2], Tu2[n + 6 * 2], fz[n + 6 * 2], ff[n + 6 * 2];
    class cfda cff;
    cff.kt = 1;
    cff.A1 = 1. / 3.0;
    cff.A2 = 2. / 3.0;
    cff.t = t;           // 1.0 / cff.PI; //计算总时间
    cff.nx_begin = -1.0; // 网格右端点
    cff.nx_end = 1.0;    // 网格左端点
    cff.CFL = CFL, cff.TVmax = 0.0;

    cff.n = n, cff.nb = nb, cff.n1 = n1;
    cff.Ltt_n0 = Ltt_n0, cff.x = x; // 声明变长数组

    cff.u_n0 = u_n0, cff.u_n1 = u_n1, cff.u_n2 = u_n2, cff.u_n3 = u_n3, cff.u_n4 = u_n4;           // 推进步中间时刻
    cff.u_n5 = u_n5, cff.u_n6 = u_n6, cff.u_n7 = u_n7, cff.u_n8 = u_n8, cff.u_n9 = u_n9;           // 推进步中间时刻
    cff.f_n0 = f_n0, cff.f_n1 = f_n1, cff.f_n2 = f_n2, cff.f_n3 = f_n3, cff.f_n4 = f_n4;           // 推进步中间时刻
    cff.f_n5 = f_n5, cff.f_n6 = f_n6, cff.f_n7 = f_n7, cff.f_n8 = f_n8, cff.f_n9 = f_n9;           // 推进步中间时刻
                                                                                                   // 声明变长数组
    cff.u_nn = u_nn, cff.f_nn = f_nn, cff.Lttx_n0 = Lttx_n0;                                       // 下一时刻（n+1)
    cff.L_n0 = L_n0, cff.L_n1 = L_n1, cff.L_n2 = L_n2, cff.L_n3 = L_n3, cff.L_n4 = L_n4;           // f=L
    cff.Lt_n0 = Lt_n0, cff.Lt_n1 = Lt_n1, cff.Lt_n2 = Lt_n2, cff.Lt_n3 = Lt_n3, cff.Lt_n4 = Lt_n4; // L_t

    cff.Ltx_n0 = Ltx_n0, cff.Ltx_n1 = Ltx_n1, cff.Ltx_n2 = Ltx_n2; // L_tx

    cff.TL_n0 = TL_n0, cff.TL_n1 = TL_n1, cff.TL_n2 = TL_n2, cff.TL_n3 = TL_n3, cff.TL_n4 = TL_n4;           // f=TL
    cff.TLt_n0 = TLt_n0, cff.TLt_n1 = TLt_n1, cff.TLt_n2 = TLt_n2, cff.TLt_n3 = TLt_n3, cff.TLt_n4 = TLt_n4; // TL_t

    cff.Tu1 = Tu1, cff.Tu2 = Tu2, cff.fz = fz, cff.ff = ff;
    cff.df = main_df, cff.u_excat = u_excat;
    cff.intc();
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
    // double th;
    cff1.tt = 0;
    n_all = cff1.n1;
    cff1.border();
    cff1.f_eq_u(cff1.f_n0, cff1.u_n0);
    cff1.compt_CFL();

    for (int i1 = 0; i1 < 10000000; i1++) // 200 50
    {
        cff1.ii1 = i1;
        cff1.tt = cff1.tt + cff1.dt;
        cff1.compt_CFL();
        if (cff1.tt > cff1.t + 1E-10)
        {
            cff1.f_eq_u(cff1.f_nn, cff1.u_nn);
            cout.precision(18); // 精度为18，正常为6
            cout << " \n comput time: " << cff1.tt - cff1.dt << " comput time n:   " << i1;
            break;
        }
        // cout << " \n comput time: " << cff1.tt << "    comput time n:   " << cff1.dt;

        if (i1 < 2)
        {

            cff1.compt_2_stage_int();

            cff1.dt = cff1.dt * 0.5;

            // cff1.RK33_compt(); // 3-RK
            // cff1.SSPRK54_compt_Shu(); // cff1.RK54LS_compt();
            cff1.RK65LS_compt2(); //   cff1.RK65LS_compt();

            // cff1.compt23_LLtt_SSP();
            // cff1.compt24_LLtt_SSP();
            // cff1.compt35_LLtt_SSP();

            if (i1 < 1)
                cff1.dt = 0.0;
        }
        else
        {
            cff1.compt_TDTS23_LLttnn();
            // cff1.compt_TDTS24_LLttnn();
            // cff1.compt_TDTS25_LLttnn();
        }
        //___________________________________________________________________
        cff1.carr(cff1.u_n0, cff1.u_nn, n_all);
        cff1.f_eq_u(cff1.f_n0, cff1.u_n0);

        // cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1;

        cff1.border();
        // cff1.maxTV();
        if (i3 * cff1.t / cff1.onet_out < cff1.tt + 0.00000001)
        {
            // cff1.Store_obj(cff1.put_one_n[i3], cff1.u_nn); //存储数据

            i3++;
        }
    }
}
