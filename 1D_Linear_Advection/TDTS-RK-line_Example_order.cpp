#include <math.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>
#include <iterator>
#include <valarray>
using namespace std;
// example:line 3-4-L-W   测试精度 算例3 weno, 不计算Tv  简化函数  计算误差值
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
    double A1, A2, A3, AA1, AA2, AA3, B0, B1, B2, B3, bb1, bb2, bb3, bbb1, bbb2, bbb3, C0, C1, C2, C3, D0, D1, D2, D3, maxuu, TVmax, k;
    double a21, a31, a32, aa12, aa21, aa31, aa32, aaa31, aaa32, w1, w2, v1, v2, v3, vv1, vv2, vv3, ww1, ww2;
    double put_obj[2][13000];   //
    double put_one_n[2][13000]; // onet_out = 10;
    ////-----------------------------------------
    void intc()
    {
        dx = (nx_end - nx_begin) / (1.0 * (n - 1.0));
        for (int i = 0; i < n1; i++)
        {
            x[i] = nx_begin + (i - nb) * dx;
            u_n0[i] = intc_fun(x[i]);
            df[i] = 1;
        }

        // for (int i = 0; i < n1; i++) //间断
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     u_n0[i] = 0.1;
        //     if (x[i] >= 0.1 && x[i] <= 0.5)
        //     {
        //         u_n0[i] = 1.0;
        //     }
        //     // df[i] = 1;
        // }
        // Write_obj(0);
    }

    double intc_fun(double xxf)
    {
        double yyf;
        // yyf = (sin(2 * PI * xxf) / 2.0 + 0.5) * 1.0; //  / PI bergure
        // yyf = sin(PI * xxf) / 2.0 + 0.25;//bergure
        yyf = 0.5 * sin(PI * xxf) + 0.5; // line
        // yyf = 0.5 * sin(PI * xxf) + 0.5; //line
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
    void compt_ThD_13_line_SSP()
    {
        int i, j;
        a21 = 0.594223212099088, v2 = 0.306027487008159;
        ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            u_nn[j] = u_n0[j] + dt * L_n0[j] + dt * dt * Lt_n0[j] / 2.0 + dt * dt * dt * TLt_n0[j] / 6.0;
        }
    }

    void compt_ThD_24_line_SSP()
    {
        int i, j;
        a21 = 0.61, A1 = a21, A2 = a21 * a21 / 2.0, A3 = a21 * a21 * a21 / 6.0;
        // B2 = 0.40031, bb2 = 0.090365;
        B2 = 0.40487, bb2 = 0.084628;
        B1 = 1. - B2;
        bb1 = 1. / 2. - a21 * B2 - bb2;
        bbb1 = 1. / 24. * (4. - 1. / a21 - 8. * a21 * a21 * B2 - 12. * a21 * bb2);
        bbb2 = 1. / (24 * a21) - 1. / 6. * a21 * (a21 * B2 + 3. * bb2);

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            f_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A2 * dt * dt * Lt_n0[j] + A3 * dt * dt * dt * TLt_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            TLt_n1[j] = dfdxxx(f_n1, j);
            u_nn[j] = u_n0[j] + dt * (B1 * L_n0[j] + B2 * L_n1[j]) +
                      dt * dt * (bb1 * Lt_n0[j] + bb2 * Lt_n1[j]) +
                      dt * dt * dt * (bbb1 * TLt_n0[j] + bbb2 * TLt_n1[j]);
        }
    }

    void compt_ThD_25_line_SSP()
    {
        int i, j;
        a21 = 0.75, A1 = a21, A2 = a21 * a21 / 2.0, A3 = a21 * a21 * a21 / 6.0;

        B1 = 1., B2 = 0.;
        bb1 = (2. - 5. * a21 + 10. * a21 * a21 * a21) / (20. * a21 * a21 * a21), bb2 = (-2. + 5. * a21) / (20. * a21 * a21 * a21);
        bbb1 = (3. + 10. * (-1. + a21) * a21) / (60. * a21 * a21), bbb2 = (3. - 5. * a21) / (60. * a21 * a21);

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            f_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A2 * dt * dt * Lt_n0[j] + A3 * dt * dt * dt * TLt_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            TLt_n1[j] = dfdxxx(f_n1, j);
            u_nn[j] = u_n0[j] + dt * (B1 * L_n0[j] + B2 * L_n1[j]) +
                      dt * dt * (bb1 * Lt_n0[j] + bb2 * Lt_n1[j]) +
                      dt * dt * dt * (bbb1 * TLt_n0[j] + bbb2 * TLt_n1[j]);
        }
    }

    void compt_ThD_36_line_SSP()
    {
        int i, j;

        // a31 = 0.615, AA1 = a31, aa31 = a31 * a31 / 2.0, aaa31 = a31 * a31 * a31 / 6.0, a32 = 0., aa32 = 0., aaa32 = 0.;
        // a21 = 0.5, A1 = a21, A2 = a21 * a21 / 2.0, A3 = a21 * a21 * a21 / 6.0;
        // B1 = 1., B2 = 0., B3 = 0.;
        // bb1 = 1. / 20. * (6. - 1. / a31 / a31),
        // bb2 = (4. * a31 * a31) / (-5 + 20. * a31 * a31);
        // bb3 = 1. / (20. * a31 * a31 - 80. * a31 * a31 * a31 * a31);
        // bbb1 = (-1. + a31 + 2. * a31 * a31) / (30. * a31 * (1. + 2. * a31));
        // bbb2 = (1. - 2. * a31 * a31) / (15. - 60. * a31 * a31);
        // bbb3 = -(1. / (60. * a31 - 240. * a31 * a31 * a31));

        a21 = 3. / 7. - sqrt(2) / 7., A1 = a21, A2 = a21 * a21 / 2.0, A3 = a21 * a21 * a21 / 6.0;
        a31 = (1. - 2. * a21) / (2. - 5. * a21), a32 = 0.;
        aa31 = a31 * a31 / 2.0, aa32 = 0.;
        aaa31 = (71. + 61. * sqrt(2)) / 14406., aaa32 = pow(-1. + 2. * a21, 3) / (6. * pow(-2. + 5. * a21, 3)) - aaa31;
        B1 = 1., B2 = 0., B3 = 0.;
        bb1 = 0.5, bb2 = 0., bb3 = 0.;
        bbb1 = (1. + 5. * a21 * (-2. + 3. * a21)) / (120. * a21 * (-1. + 2. * a21));
        bbb2 = 1. / (120. * a21 * (1. + a21 * (-4. + 5. * a21)));
        bbb3 = pow(-2. + 5. * a21, 3) / (120. * (-1. + 2. * a21) * (1. + a21 * (-4. + 5. * a21)));
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            f_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A2 * dt * dt * Lt_n0[j] + A3 * dt * dt * dt * TLt_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            TLt_n1[j] = dfdxxx(f_n1, j);
            f_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) +
                      dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]) +
                      dt * dt * dt * (aaa31 * TLt_n0[j] + aaa32 * TLt_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;

            L_n2[j] = dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            TLt_n2[j] = dfdxxx(f_n2, j);

            u_nn[j] = u_n0[j] + dt * (B1 * L_n0[j] + B2 * L_n1[j] + B3 * L_n2[j]) +
                      dt * dt * (bb1 * Lt_n0[j] + bb2 * Lt_n1[j] + bb3 * Lt_n2[j]) +
                      dt * dt * dt * (bbb1 * TLt_n0[j] + bbb2 * TLt_n1[j] + bbb3 * TLt_n2[j]);
        }
    }

    // //-----------------------------------------
    void compt_2_stage_int()
    {
        int i, j;

        A1 = 0.765; // 0.525  0.680 0.765
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
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
            L_n1[j] = dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];

            TL_n1[j] = TL_n3[j], TL_n0[j] = TL_n2[j];
            TLt_n1[j] = TLt_n3[j], TLt_n0[j] = TLt_n2[j];

            TL_n3[j] = L_n1[j], TL_n2[j] = L_n0[j];
            TLt_n3[j] = Lt_n1[j], TLt_n2[j] = Lt_n0[j];
        }
    }

    void compt23_LLttnn_line()
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
            L_n0[j] = dfdx(f_n0, j);
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
            L_n1[j] = dfdx(f_n1, j);
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

    void compt24_LLttnn_line()
    {
        int i, j;

        A1 = 0.5, B0 = 1., B1 = 0., C0 = 1. / 3., C1 = 2. / 3.;
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
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
            L_n1[j] = dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j] + B1 * L_n1[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt34_LLttnn_line()
    {
        int i, j;
        // A1 = 1. / 4., A2 = 1. / 2., B0 = 1. / 3.;
        A1 = 1. / 4., A2 = 2. / 3., B0 = 1. / 3.; // 精度足够 line稳定性足够,  Bergurs稳定性足够
        // A1 = 0.5, A2 = 0.8, B0 = 0.6;
        B1 = -((pow(A2, 3) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3)));
        B2 = -((pow(A1, 3) * (-1. + B0)) / (pow(A1, 3) - pow(A2, 3)));
        C0 = -((-1. + 2. * A1 + 2. * A2 - 6. * A1 * A2) / (6. * A1 * A2)) + (A2 * (-pow(A1, 3) + A1 * A2 * A2) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        C1 = -((-1. + 2. * A2) / (6. * A1 * (A1 - A2))) + (A1 * pow(A2, 3) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        C2 = -((1. - 2. * A1) / (6. * (A1 - A2) * A2)) - (pow(A1, 3) * A2 * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        border();
        i = n1 - nb * 2;
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdx(L_n0, j);
            u_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A1 * A1 / 2.0 * dt * dt * Lt_n0[j];
            u_n2[j] = u_n0[j] + A2 * dt * L_n0[j] + A2 * A2 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            L_n2[j] = dfdx(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdx(L_n1, j), Lt_n2[j] = dfdx(L_n2, j); // Lt_n1[j] = -1. * dfdxcent(L_n1, j),  Lt_n2[j] = -1. * dfdxcent(L_n2, j);
            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j] + B1 * L_n1[j] + B2 * L_n2[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j] + C2 * Lt_n2[j]);

            // TL_n1[j] = TL_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = TLt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt35_LLttnn_line()
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
            L_n0[j] = dfdx(f_n0, j);
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
            L_n1[j] = dfdx(f_n1, j);
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
            L_n2[j] = dfdx(f_n2, j);
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

    //------------------------------------------------

    void compt_TDTS23_LLttnn_line()
    {

        int i, j;

        // a21 = 0.532, v2 = 0.389;
        // ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.1;
        a21 = 0.522, v2 = 0.448;
        ww1 = 0.00, ww2 = 0.00, w1 = 0.04, w2 = 0.00;
        v1 = 1. - v2 - w1 - w2;
        vv1 = -((1 - 3 * w1 - 3 * w2 + 6 * ww1 + 3 * a21 * (-1 - 2 * w1 + a21 * (v2 + w2) + 2 * ww1) + 6 * ww2) / (6. * a21));
        vv2 = (1 - 3 * w1 - 3 * w2 + 6 * ww1 + 6 * ww2 - 3 * a21 * (-2 * w2 + a21 * (v2 + w2) + 2 * ww2)) / (6. * a21);

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
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
            L_n1[j] = dfdx(f_n1, j);
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
            L_n0[j] = dfdx(f_n0, j);
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
            L_n1[j] = dfdx(f_n1, j);
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

    void compt_TDTS25_LLttnn_line()
    {
        int i, j;
        a21 = 0.765; // 1. / 210. * (93. + sqrt(2139.));
        vv1 = (-31. + 6. * a21 * (-31. + 85. * a21)) / (360. * pow(a21, 2));
        vv2 = 31. / (360. * pow(a21, 2));
        v1 = -(1. / 2.) + 31. / (30. * a21), v2 = 0.;
        w1 = 3. / 2. - 31. / (30. * a21), w2 = 0.;
        ww1 = (31. + 6. * a21 * (-31. + 35. * a21)) / (360. * pow(a21, 2));
        ww2 = -(31. / (360. * pow(a21, 2)));

        // a21 = 0.5;
        // w2 = -0.1;
        // ww2 = -0.0, v1 = 1.0;
        // v2 = -w2, w1 = 1. - v1;
        // vv1 = -((1. - 14. * a21 + 2. * v1 + 6. * a21 * v1 - 2. * w2 - 6. * a21 * w2 + 12. * a21 * ww2) / (12. * a21));
        // vv2 = -((-1. - 2. * v1 + 2. * w2 - 12. * a21 * a21 * w2 - 12. * a21 * ww2 + 12. * a21 * a21 * ww2) / (12. * a21 * (1. + a21)));
        // ww1 = -((-5. - 4. * a21 + 4. * v1 + 6. * a21 * v1 - 4. * w2 + 6. * a21 * w2 + 12. * ww2 - 12. * a21 * ww2) / (12. * (1. + a21)));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
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
            L_n1[j] = dfdx(f_n1, j);
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

    void compt_TDTS35_LLttnn_line()
    {
        int i, j;

        // a21 = 0.5;
        // w2 = -0.1;
        // ww2 = -0.0, v1 = 1.0;
        // v2 = -w2, w1 = 1. - v1;
        // vv1 = -((1. - 14. * a21 + 2. * v1 + 6. * a21 * v1 - 2. * w2 - 6. * a21 * w2 + 12. * a21 * ww2) / (12. * a21));
        // vv2 = -((-1. - 2. * v1 + 2. * w2 - 12. * a21 * a21 * w2 - 12. * a21 * ww2 + 12. * a21 * a21 * ww2) / (12. * a21 * (1. + a21)));
        // ww1 = -((-5. - 4. * a21 + 4. * v1 + 6. * a21 * v1 - 4. * w2 + 6. * a21 * w2 + 12. * ww2 - 12. * a21 * ww2) / (12. * (1. + a21)));

        // a21 = 0.381, vv3 = 0.21, a31 = 0.61;
        // v2 = 0., v3 = 0., a32 = 0., aa21 = a21 * a21 / 2.0;
        // vv1 = (-31. + 2. * a21 * (-4 + 85 * a21) - 60. * (a21 - a31) * (1. + a31) * (1. + 2. * a31 + a21 * (2 + 6 * a31)) * vv3) / (60. * a21 * (1. + 2. * a21));
        // vv2 = (31 - 60. * a31 * (1. + a31) * (1. + 2. * a31) * vv3) / (60. * a21 * (1. + a21) * (1. + 2. * a21));
        // v1 = (13 - 5 * a21 + 60. * (a21 - a31) * a31 * (1. + a31) * vv3) / (5 + 10. * a21);
        // ww1 = (-27 + 2. * a21 * (6 + 35 * a21) - 60. * (a21 - a31) * a31 * (3 + 4. * a31 + a21 * (4 + 6 * a31)) * vv3) / (60. * (1. + a21) * (1. + 2. * a21));
        // w1 = (-8 + 15 * a21 - 60. * (a21 - a31) * a31 * (1. + a31) * vv3) / (5 + 10. * a21);
        // aa31 = (-31 * pow(a21, 2) + 60. * a31 * (pow(a21, 2) + 3. * a21 * (1. + 2. * a21 * (2 + a21)) * a31 - (1. + 3. * a21) * pow(a31, 2)) * vv3) / (360. * a21 * (1. + a21) * (1. + 2. * a21) * vv3);
        // aa32 = (31 * pow(a21, 2) - 60. * (a21 - a31) * a31 * (a21 + a31 + 3. * a21 * a31) * vv3) / (360. * a21 * (1. + a21) * (1. + 2. * a21) * vv3);

        // a21 = 0.81, a31 = 0.31;
        // a32 = -0.1, v2 = 0., v3 = 0., aa12 = a31 + a32, aa21 = a21 * a21 / 2.0;
        // vv1 = (-70. * pow(aa12, 2) - 3. * (-9 + 4. * a31 + 4. * a32) + 4. * a21 * (-3 + 15. * pow(a31, 2) + 2. * a31 * (7 + 15. * a32) + a32 * (14 + 15. * a32)) + 10. * pow(a21, 2) * (-7 + 30. * pow(a31, 2) + 6. * a32 * (1 + 5. * a32) + a31 * (6 + 60. * a32))) / (60. * a21 * (aa12) * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32)));
        // vv2 = -((-27. + 70. * pow(a31, 2) + 4. * a31 * (3. + 35. * a32) + 2. * a32 * (6. + 35. * a32)) / (60. * a21 * (a21 - a31 - a32) * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32))));
        // vv3 = (-27. + 2. * a21 * (6 + 35. * a21)) / (60. * (a21 - a31 - a32) * (aa12) * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32)));
        // v1 = (12. + 25. * a31 + 25. * a32 + 5. * a21 * (5. + 4. * a31 + 4. * a32)) / (5. * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32)));
        // w1 = (3. - 5. * a31 - 5. * a32 + 5. * a21 * (-1. + 2. * a31 + 2. * a32)) / (5. * (3. + 4. * a31 + 4. * a32 + a21 * (4. + 6. * a31 + 6. * a32)));
        // aa31 = (-70. * pow(a21, 4) * a32 + 9. * pow(aa12, 3) - a21 * pow(aa12, 2) * (27. + 4. * a31 + 4. * a32) + 2. * pow(a21, 3) * (-6. * a32 + 35. * pow(aa12, 2)) + pow(a21, 2) * (16. * pow(a31, 2) + 2. * a32 * (9. + 8. * a32) + a31 * (-9. + 32. * a32))) / (2. * a21 * (-27. + 2. * a21 * (6. + 35. * a21)));
        // aa32 = (-12. * pow(a21, 3) * a32 - 70. * pow(a21, 4) * a32 - 9. * pow(aa12, 3) + 4. * a21 * pow(aa12, 3) + pow(a21, 2) * (-4. * pow(a31, 2) + a31 * (9. - 8. * a32) - 4. * (-9 + a32) * a32)) / (2. * a21 * (-27. + 2. * a21 * (6. + 35. * a21)));

        w1 = 0.05, a21 = 0.32;
        ww1 = 0., v2 = 0., v3 = 0., a32 = 0., aa12 = (1. - 2. * a21 + 4. * w1 + 6. * a21 * w1), aa21 = a21 * a21 / 2.0;
        vv1 = (1. - 8. * a21 + 10. * pow(a21, 2) - 88. * w1 + 64 * a21 * w1 + 300. * pow(a21, 2) * w1 + 10. * pow(w1, 2) + 60. * a21 * pow(w1, 2) + 60. * pow(a21, 2) * pow(w1, 2)) / (12. * a21 * (-3. + 5. * a21 + 15. * w1 + 20. * a21 * w1));
        vv2 = (-1. + 88. * w1 - 10. * pow(w1, 2)) / (12. * a21 * (-3. + 10. * a21 - 10. * pow(a21, 2) + 15. * w1 + 40. * a21 * w1 + 30. * pow(a21, 2) * w1));
        vv3 = (25. * pow(aa12, 3)) / (12. * (9. - 45. * a21 + 80. * pow(a21, 2) - 50. * pow(a21, 3) - 90. * w1 + 45. * a21 * w1 + 160. * pow(a21, 2) * w1 - 50. * pow(a21, 3) * w1 + 225. * pow(w1, 2) + 900. * a21 * pow(w1, 2) + 1250. * pow(a21, 2) * pow(w1, 2) + 600. * pow(a21, 3) * pow(w1, 2)));
        v1 = 1. - w1;
        a31 = (3. - 5. * a21 - 15. * w1 - 20. * a21 * w1) / (5. * aa12);
        aa31 = (-9. + 90. * a21 - 320. * pow(a21, 2) + 475. * pow(a21, 3) - 250. * pow(a21, 4) + 135. * w1 - 540. * a21 * w1 + 960. * pow(a21, 2) * w1 + 100. * pow(a21, 3) * w1 - 1250. * pow(a21, 4) * w1 - 675. * pow(w1, 2) - 1350. * a21 * pow(w1, 2) - 1800. * pow(a21, 2) * pow(w1, 2) + 50. * pow(a21, 3) * pow(w1, 2) + 2000. * pow(a21, 4) * pow(w1, 2) + 1125. * pow(w1, 3) + 9000. * a21 * pow(w1, 3) + 25000. * pow(a21, 2) * pow(w1, 3) + 29000. * pow(a21, 3) * pow(w1, 3) + 12000. * pow(a21, 4) * pow(w1, 3)) / (250. * a21 * pow(aa12, 3));
        aa32 = -(((-1. + 5. * w1) * (9 - 45. * a21 + 80. * pow(a21, 2) - 50. * pow(a21, 3) - 90. * w1 + 45. * a21 * w1 + 160. * pow(a21, 2) * w1 - 50. * pow(a21, 3) * w1 + 225. * pow(w1, 2) + 900. * a21 * pow(w1, 2) + 1250. * pow(a21, 2) * pow(w1, 2) + 600. * pow(a21, 3) * pow(w1, 2))) / (250. * a21 * pow(aa12, 3)));

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
            L_n0[j] = dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdx(L_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdx(L_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = dfdx(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n2[j] = dfdx(L_n2, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + v3 * L_n2[j] + w1 * TL_n0[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j] + ww1 * TLt_n0[j]);

            // u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n1[j] + w2 * TL_n1a[j]) +
            //           dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n1[j] + ww2 * TLt_n1a[j]);

            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }
    //------------------------------------------------
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
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u21 * u_n0[j] + u22 * Tu1[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = dfdx(f_n2, j);
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
            L_n3[j] = dfdx(f_n3, j);
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
            L_n4[j] = dfdx(f_n4, j);
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
            L_n0[j] = dfdx(f_n5, j);
            Lt_n0[j] = dfdxx(f_n5, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j] + b14 * L_n4[j] + b15 * L_n0[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j] + bb14 * Lt_n4[j] + bb15 * Lt_n0[j]);

            Tu1[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j] + b24 * L_n4[j] + b25 * L_n0[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j] + bb24 * Lt_n4[j] + bb25 * Lt_n0[j]);
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
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = dfdx(f_n2, j);
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
            L_n3[j] = dfdx(f_n3, j);
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
            L_n4[j] = dfdx(f_n4, j);
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
            L_n0[j] = dfdx(f_n5, j);
            Lt_n0[j] = dfdxx(f_n5, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = u_n0[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j] + b14 * L_n4[j] + b15 * L_n0[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j] + bb14 * Lt_n4[j] + bb15 * Lt_n0[j]);

            Tu1[j] = u_n0[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j] + b24 * L_n4[j] + b25 * L_n0[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j] + bb24 * Lt_n4[j] + bb25 * Lt_n0[j]);
        }
    }
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
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u21 * u_n0[j] + u22 * Tu1[j] + a21 * dt * L_n1[j] + aa21 * dt * dt * Lt_n1[j];
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = dfdx(f_n2, j);
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
            L_n3[j] = dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);
            u_nn[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b11 * L_n1[j] + b12 * L_n2[j] + b13 * L_n3[j]) +
                      dt * dt * (bb11 * Lt_n1[j] + bb12 * Lt_n2[j] + bb13 * Lt_n3[j]);

            Tu1[j] = (1 - v) * u_n0[j] + v * Tu1[j] + dt * (b21 * L_n1[j] + b22 * L_n2[j] + b23 * L_n3[j]) +
                     dt * dt * (bb21 * Lt_n1[j] + bb22 * Lt_n2[j] + bb23 * Lt_n3[j]);
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
            L_n0[j] = dfdx(u_n0, j);
            Lt_n0[j] = dfdxx(u_n0, j);
            Ltt_n0[j] = dfdxxx(u_n0, j);
            Tu1[j] = u_n0[j] + 0.270774947096390 * dt * L_n0[j] + 0.104435936229161 * dt * dt * Lt_n0[j] - 0.141473113475966 * dt * dt * dt * Ltt_n0[j];
            // Tu1[j] = u_n0[j] + 0.527015480634761 * dt * L_n0[j] + 0.027140536548443 * dt * dt * Lt_n0[j];
            // Tu1[j] = u_n0[j] + 0.032910530894256 * dt * L_n0[j] + 0.000541551521871 * dt * dt * Lt_n0[j];
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
    void RK22_compt()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n1, nb + j);
            u_nn[nb + j] = u_n0[nb + j] * 0.5 + u_n1[nb + j] * 0.5 + 0.5 * dt * df[nb + j];
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
            // df[nb + j] =dfdx(f_n0, nb + j); //-uu_t0[nb + j]; -dfdx
            df[nb + j] = dfdx(f_n0, nb + j);
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        j = 0, i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] =dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            df[nb + j] = dfdx(f_n1, nb + j);
            u_n2[nb + j] = u_n0[nb + j] * 0.75 + u_n1[nb + j] * 0.25 + 0.25 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        j = 0, i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] =dfdx(f_n2, nb + j);
            df[nb + j] = dfdx(f_n2, nb + j);
            u_nn[nb + j] = u_n0[nb + j] * 1.0 / 3.0 + u_n2[nb + j] * 2.0 / 3.0 + 2.0 / 3.0 * dt * df[nb + j];
            j++;
        }
    }

    void RK33_compt()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] * 0.75 + 0.25 * u_n1[nb + j] + 0.25 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n2, nb + j);
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
            L_n0[nb + j] = dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.391752226571890 * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.444370493651235 * u_n0[nb + j] + 0.555629506348765 * u_n1[nb + j] + 0.368410593050371 * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.620101851488403 * u_n0[nb + j] + 0.379898148511597 * u_n2[nb + j] + 0.251891774271694 * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = 0.178079954393132 * u_n0[nb + j] + 0.821920045606868 * u_n3[nb + j] + 0.544974750228521 * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n4, nb + j);
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
            L_n0[nb + j] = dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 4. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] + 1. / 8. * dt * L_n0[nb + j] + 1. / 8. * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 0.0 * dt * L_n0[nb + j] + 0.0 * dt * L_n1[nb + j] +
                           1. / 2. * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] + 3. / 16. * dt * L_n0[nb + j] - 3. / 8. * dt * L_n1[nb + j] +
                           3. / 8. * dt * L_n2[nb + j] + 9. / 16. * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = dfdx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] - 3. / 7. * dt * L_n0[nb + j] + 8. / 7. * dt * L_n1[nb + j] + 6. / 7. * dt * L_n2[nb + j] -
                           12. / 7. * dt * L_n3[nb + j] + 8. / 7. * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n5, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 7. / 90. * dt * L_n0[nb + j] + 0. * dt * L_n1[nb + j] +
                           16. / 45. * dt * L_n2[nb + j] + 2. / 15. * dt * L_n3[nb + j] +
                           16. / 45. * dt * L_n4[nb + j] + 7. / 90. * dt * df[nb + j];
            j++;
        }
    }

    void RK76LS_compt()
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 6. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] + 12. / 169. * dt * L_n0[nb + j] + 27. / 169. * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 107952. / 571787. * dt * L_n0[nb + j] - 406107. / 571787. * dt * L_n1[nb + j] + 566826. / 571787. * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] + 1371923. / 6669000. * dt * L_n0[nb + j] - 819. / 3800. * dt * L_n1[nb + j] +
                           19411847. / 88236000. * dt * L_n2[nb + j] + 561846173. / 1147068000. * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = dfdx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] - 1563412. / 5835375. * dt * L_n0[nb + j] + 468. / 475. * dt * L_n1[nb + j] - 14488201. / 168199875. * dt * L_n2[nb + j] -
                           1711096709. / 6846562125. * dt * L_n3[nb + j] + 648832. / 1549583. * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            Lt_n1[nb + j] = dfdx(f_n5, nb + j); //-uu_t1[nb + j];
            u_n6[nb + j] = u_n0[nb + j] + 120115. / 277641. * dt * L_n0[nb + j] - 117. / 113. * dt * L_n1[nb + j] + 219237109. / 296102601. * dt * L_n2[nb + j] +
                           29855183083. / 44628054003. * dt * L_n3[nb + j] - 3009600. / 9215941. * dt * L_n4[nb + j] + 297825. / 572797. * dt * Lt_n1[nb + j];
            f_n6[nb + j] = f_u(u_n6[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = dfdx(f_n6, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 137. / 1872. * dt * L_n0[nb + j] + 371293. / 1164612. * dt * L_n2[nb + j] +
                           3939040643. / 23169727152. * dt * L_n3[nb + j] + 19000. / 104859. * dt * L_n4[nb + j] +
                           45125. / 243312. * dt * Lt_n1[nb + j] + 113. / 1584. * dt * df[nb + j];
            j++;
        }
    }

    //---------------------------------------------------------
    void compt_CFL()
    {
        double *umax;
        umax = max_element(u_n0 + 1, u_n0 + n1 - 1);
        // cout << *umax;
        //  dt = CFL * dx / (*umax);
        // dt=0.02/1.0;
        dt = CFL * dx;
    }
    ////-----------------------------------------
    int Perimeter()
    {
        cout << " yes ";
        return 2 * t;
    }

    double df_dx(double *g, int i)
    {
        double y;
        // 1up
        //  y = (*(g + i) - *(g + i - 1)) / dx;
        // 5up
        //  y = (*(g + i - 3) / (-30.) + *(g + i - 2) / 4. - *(g + i - 1) + *(g + i) / 3. + *(g + i + 1) / 2. - *(g + i + 2) / 20.) / dx;
        // 7up
        y = (*(g + i - 4) * 3. - *(g + i - 3) * 28. + *(g + i - 2) * 126. - *(g + i - 1) * 420. + *(g + i) * 105. + *(g + i + 1) * 252. - *(g + i + 2) * 42. + *(g + i + 3) * 4.) / 420. / dx;
        // 13 order
        // y = (*(g + i - 6) * 1.803751803752541E-004 - *(g + i - 5) * 2.597402597403613E-003 + *(g + i - 4) * 1.785714285714940E-002 -
        //      *(g + i - 3) * 7.936507936510601E-002 + *(g + i - 2) * 0.267857142857228 - *(g + i - 1) * 0.857142857144015 +
        //      *(g + i) * 2.101706824838899E-012 +
        //      *(g + i + 1) * 0.85714285714181 - *(g + i + 2) * 0.267857142857150 + *(g + i + 3) * 7.936507936510066E-002 -
        //      *(g + i + 4) * 1.785714285714922E-002 + *(g + i + 5) * 2.597402597403614E-003 - *(g + i + 6) * 1.803751803752537E-004) /
        //     dx;

        //  1.803751803752541E-004 -2.597402597403613E-003  1.785714285714940E-002
        //  -7.936507936510601E-002  0.267857142857228      -0.857142857144015
        //   2.101706824838899E-012
        //   0.857142857141818      -0.267857142857150  7.936507936510066E-002
        //    -1.785714285714922E-002  2.597402597403614E-003 -1.803751803752537E-004

        /*  //weno_5-js    ep = 1.E-6,
        double ep = 1.E-15, C03 = 3. / 10., C13 = 3. / 5., C23 = 1. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i + 1)) + *(g + i + 2)), 2) + 0.25 * pow((3. * (*(g + i)) - 4. * (*(g + i + 1)) + *(g + i + 2)), 2);
            sn1 = 13. / 12. * pow((*(g + i - 1) - 2. * (*(g + i)) + *(g + i + 1)), 2) + 0.25 * pow((*(g + i - 1) - *(g + i + 1)), 2);
            sn2 = 13. / 12. * pow((*(g + i - 2) - 2. * (*(g + i - 1)) + *(g + i)), 2) + 0.25 * pow((*(g + i - 2) - 4. * (*(g + i - 1)) + 3. * (*(g + i))), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 1. / 3. * (*(g + i)) + 5. / 6. * (*(g + i + 1)) - 1. / 6. * (*(g + i + 2));
            q13 = -1. / 6. * (*(g + i - 1)) + 5. / 6. * (*(g + i)) + 1. / 3. * (*(g + i + 1));
            q23 = 1. / 3. * (*(g + i - 2)) - 7. / 6. * (*(g + i - 1)) + 11. / 6. * (*(g + i));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
         */

        /*  //weno_5-js 负通量
        double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i + 2) - 2. * (*(g + i + 1)) + *(g + i)), 2) + 0.25 * pow((3. * (*(g + i)) - 4. * (*(g + i + 1)) + *(g + i + 2)), 2);
            sn1 = 13. / 12. * pow((*(g + i + 1) - 2. * (*(g + i)) + *(g + i - 1)), 2) + 0.25 * pow((*(g + i + 1) - *(g + i - 1)), 2);
            sn2 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i - 1)) + *(g + i - 2)), 2) + 0.25 * pow((*(g + i - 2) - 4. * (*(g + i - 1)) + 3. * (*(g + i))), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 1. / 3. * (*(g + i + 2)) - 7. / 6. * (*(g + i + 1)) + 11. / 6. * (*(g + i));
            q13 = -1. / 6. * (*(g + i + 1)) + 5. / 6. * (*(g + i)) + 1. / 3. * (*(g + i - 1));
            q23 = 1. / 3. * (*(g + i)) + 5. / 6. * (*(g + i - 1)) - 1. / 6. * (*(g + i - 2));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
         */

        /* //weno-z
        double pe = 2.0, epb = 1.E-6, tau, sum_a, aaz[3], wwz[3], qqz[3], b[3];
        double cc0 = 1.0 / 10.0, cc1 = 3.0 / 5.0, cc2 = 3.0 / 10.0, ccz[3], vz[2]; //!理想权重
        ccz[0] = cc0, ccz[1] = cc1, ccz[2] = cc2;
        for (int j = 0; j < 2; j++)
        {
            b[2] = 13. / 12. * pow((*(g + i) - 2. * (*(g + i + 1)) + *(g + i + 2)), 2) + 0.25 * pow((3. * (*(g + i)) - 4. * (*(g + i + 1)) + *(g + i + 2)), 2);
            b[1] = 13. / 12. * pow((*(g + i - 1) - 2. * (*(g + i)) + *(g + i + 1)), 2) + 0.25 * pow((*(g + i - 1) - *(g + i + 1)), 2);
            b[0] = 13. / 12. * pow((*(g + i - 2) - 2. * (*(g + i - 1)) + *(g + i)), 2) + 0.25 * pow((*(g + i - 2) - 4. * (*(g + i - 1)) + 3. * (*(g + i))), 2);
            tau = abs(b[2] - b[0]);
            for (int jj = 0; jj < 3; jj++)
            {
                aaz[jj] = ccz[jj] * (1.0 + pow((tau / (b[jj] + epb)), 2));
            }
            sum_a = aaz[0] + aaz[1] + aaz[2];
            wwz[0] = aaz[0] / sum_a, wwz[1] = aaz[1] / sum_a, wwz[2] = aaz[2] / sum_a; //!正式权重
            qqz[0] = 1.0 / 3.0 * (*(g + i - 2)) - 7.0 / 6.0 * (*(g + i - 1)) + 11.0 / 6.0 * (*(g + i));
            qqz[1] = -1.0 / 6.0 * (*(g + i - 1)) + 5.0 / 6.0 * (*(g + i)) + 1.0 / 3.0 * (*(g + i + 1));
            qqz[2] = 1.0 / 3.0 * (*(g + i)) + 5.0 / 6.0 * (*(g + i + 1)) - 1.0 / 6.0 * (*(g + i + 2));

            vz[j] = wwz[0] * qqz[0] + wwz[1] * qqz[1] + wwz[2] * qqz[2];
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
         */

        /*  //weno_7
        double uim3, uim2, uim1, ui, uip1, uip2, uip3, vz[2];
        double uz1, uz2, uz3, uz4, beta0, beta1, beta2, beta3;
        double omt0, omt1, omt2, omt3, omts;
        double g0 = 1. / 35., g1 = 12. / 35., g2 = 18. / 35., g3 = 4. / 35., eps = 1e-6;

        for (int j = 0; j < 2; j++)
        {
            uim3 = *(g + i - 3), uim2 = *(g + i - 2), uim1 = *(g + i - 1);
            ui = *(g + i), uip1 = *(g + i + 1), uip2 = *(g + i + 2), uip3 = *(g + i + 3);

            beta0 = uim3 * (547. * uim3 - 3882. * uim2 + 4642. * uim1 - 1854. * ui) +
                    uim2 * (7043. * uim2 - 17246. * uim1 + 7042. * ui) +
                    uim1 * (11003. * uim1 - 9402. * ui) + 2107. * pow(ui, 2);
            beta1 = uim2 * (267. * uim2 - 1642. * uim1 + 1602. * ui - 494. * uip1) +
                    uim1 * (2843. * uim1 - 5966. * ui + 1922. * uip1) +
                    ui * (3443. * ui - 2522. * uip1) + 547. * pow(uip1, 2);
            beta2 = uim1 * (547. * uim1 - 2522. * ui + 1922. * uip1 - 494. * uip2) +
                    ui * (3443. * ui - 5966. * uip1 + 1602. * uip2) +
                    uip1 * (2843. * uip1 - 1642. * uip2) + 267. * pow(uip2, 2);
            beta3 = ui * (2107. * ui - 9402. * uip1 + 7042. * uip2 - 1854. * uip3) +
                    uip1 * (11003. * uip1 - 17246. * uip2 + 4642. * uip3) +
                    uip2 * (7043. * uip2 - 3882. * uip3) + 547. * pow(uip3, 2);
            uz1 = (-1. / 4.) * uim3 + (13. / 12.) * uim2 - (23. / 12.) * uim1 + (25. / 12.) * ui;
            uz2 = (1. / 12.) * uim2 - (5. / 12.) * uim1 + (13. / 12.) * ui + (1. / 4.) * uip1;
            uz3 = (-1. / 12.) * uim1 + (7. / 12.) * ui + (7. / 12.) * uip1 - (1. / 12.) * uip2;
            uz4 = (1. / 4.) * ui + (13. / 12.) * uip1 - (5. / 12.) * uip2 + (1. / 12.) * uip3;

            omt0 = g0 * pow(eps + beta0, -2), omt1 = g1 * pow(eps + beta1, -2), omt2 = g2 * pow(eps + beta2, -2);
            omt3 = g3 * pow(eps + beta3, -2), omts = omt0 + omt1 + omt2 + omt3;
            vz[j] = (omt0 * uz1 + omt1 * uz2 + omt2 * uz3 + omt3 * uz4) / omts;
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

    double dfdxcent(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i) - *(g + i - 1)) / dx;
        // y = (  *(g + i - 2) *0.5 - *(g + i - 1)*2.0 + *(g + i) *1.5) / dx;
        // 3up
        //  y = (*(g + i - 3) / (-3.) + *(g + i - 2) *1.5 - *(g + i - 1)/3.0 + *(g + i) *11.0/ 6.0) / dx;
        // 5up
        // y = (*(g + i - 3) / (-30.) + *(g + i - 2) / 4. - *(g + i - 1) + *(g + i) / 3. + *(g + i + 1) / 2. - *(g + i + 2) / 20.) / dx;
        // 7up
        // y = (*(g + i - 4) * 3. - *(g + i - 3) * 28. + *(g + i - 2) * 126. - *(g + i - 1) * 420. + *(g + i) * 105. + *(g + i + 1) * 252. - *(g + i + 2) * 42. + *(g + i + 3) * 4.) / 420. / dx;
        // y = (*(g + i - 5) * (-4.) + *(g + i - 4) * 35. - *(g + i - 3) * 140. + *(g + i - 2) * 350. - *(g + i - 1) * 700. + *(g + i) * 3221. + *(g + i + 1) * 140. - *(g + i + 2) * 10.) / 420. / dx;
        // y = (45. * (*(g + i + 1) - *(g + i - 1)) - 21. * (*(g + i + 2) - *(g + i - 2)) + (*(g + i + 3) - *(g + i - 3))) / 60. / dx;
        // y = -(*(g + i + 4) * 3. - *(g + i + 3) * 28. + *(g + i + 2) * 126. - *(g + i + 1) * 420. + *(g + i) * 105. + *(g + i - 1) * 252. - *(g + i - 2) * 42. + *(g + i - 3) * 4.) / 420. / dx;

        // y = (*(g + i + 1) - *(g + i - 1)) / 2.0 / dx;
        // 3 up
        // y = (*(g + i - 1) / (-3.) + *(g + i) / (-2.) + *(g + i + 1) + *(g + i + 2) / (-6.)) / dx;
        // 4 center
        y = (*(g + i - 2) / (12.) + *(g + i - 1) * (-8. / 12.) + *(g + i + 1) * (8. / 12.) + *(g + i + 2) / (-12.)) / dx;

        /*  //weno_5-js 负通量
        double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i + 2) - 2. * (*(g + i + 1)) + *(g + i)), 2) + 0.25 * pow((3. * (*(g + i)) - 4. * (*(g + i + 1)) + *(g + i + 2)), 2);
            sn1 = 13. / 12. * pow((*(g + i + 1) - 2. * (*(g + i)) + *(g + i - 1)), 2) + 0.25 * pow((*(g + i + 1) - *(g + i - 1)), 2);
            sn2 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i - 1)) + *(g + i - 2)), 2) + 0.25 * pow((*(g + i - 2) - 4. * (*(g + i - 1)) + 3. * (*(g + i))), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 1. / 3. * (*(g + i + 2)) - 7. / 6. * (*(g + i + 1)) + 11. / 6. * (*(g + i));
            q13 = -1. / 6. * (*(g + i + 1)) + 5. / 6. * (*(g + i)) + 1. / 3. * (*(g + i - 1));
            q23 = 1. / 3. * (*(g + i)) + 5. / 6. * (*(g + i - 1)) - 1. / 6. * (*(g + i - 2));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
         */

        /*  //weno_5-js
        double ep = 1.E-6, C03 = 3. / 10., C13 = 3. / 5., C23 = 1. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i + 1)) + *(g + i + 2)), 2) + 0.25 * pow((3. * (*(g + i)) - 4. * (*(g + i + 1)) + *(g + i + 2)), 2);
            sn1 = 13. / 12. * pow((*(g + i - 1) - 2. * (*(g + i)) + *(g + i + 1)), 2) + 0.25 * pow((*(g + i - 1) - *(g + i + 1)), 2);
            sn2 = 13. / 12. * pow((*(g + i - 2) - 2. * (*(g + i - 1)) + *(g + i)), 2) + 0.25 * pow((*(g + i - 2) - 4. * (*(g + i - 1)) + 3. * (*(g + i))), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 1. / 3. * (*(g + i)) + 5. / 6. * (*(g + i + 1)) - 1. / 6. * (*(g + i + 2));
            q13 = -1. / 6. * (*(g + i - 1)) + 5. / 6. * (*(g + i)) + 1. / 3. * (*(g + i + 1));
            q23 = 1. / 3. * (*(g + i - 2)) - 7. / 6. * (*(g + i - 1)) + 11. / 6. * (*(g + i));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
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
        // y = (*(g + i) - *(g + i - 1)) / dx;
        // 5up
        // y = (   *(g + i - 2) * 3. - *(g + i - 1) * 30. - *(g + i) * 20. + *(g + i + 1) * 60. - *(g + i + 2) * 15. + *(g + i + 3) * 2. ) / 60. / dx;

        // 7up
        // y = (-*(g + i - 3) * 4. + *(g + i - 2) * 42. - *(g + i - 1) * 252. - *(g + i) * 105. + *(g + i + 1) * 420. - *(g + i + 2) * 126. + *(g + i + 3) * 28. - *(g + i + 4) * 3.) / 420. / dx;

        // 8up
        y = (-*(g + i - 3) * 5. + *(g + i - 2) * 60. - *(g + i - 1) * 420. - *(g + i) * 378 + *(g + i + 1) * 1050. - *(g + i + 2) * 420. + *(g + i + 3) * 140. - *(g + i + 4) * 30. + *(g + i + 5) * 3.) / 840. / dx;

        // 9up
        //  y = (*(g + i - 4) - *(g + i - 3) * 12. + *(g + i - 2) * 72. - *(g + i - 1) * 336. - *(g + i) * 100.8 + *(g + i + 1) * 504. - *(g + i + 2) * 168. + *(g + i + 3) * 48. - *(g + i + 4) * 9. + *(g + i + 5) * 0.8) / 504. / dx;

        return y;
    }

    double dfdxx(double *g, int i)
    {
        double y;
        // 2ord
        // y = (*(g + i + 1) - *(g + i) * 2.0 + *(g + i - 1)) / dx / dx;
        // 6ord
        // y = (*(g + i + 3) * 2. - *(g + i + 2) * 27.0 + *(g + i + 1) * 270.0 - *(g + i) * 490.0 + *(g + i - 1) * 270.0 - *(g + i - 2) * 27.0 + *(g + i - 3) * 2.) / dx / dx / 180.0;
        // 8ord
        y = (-*(g + i - 4) * 63. + *(g + i - 3) * 896. - *(g + i - 2) * 7056.0 + *(g + i - 1) * 56448.0 - *(g + i) * 100450.0 + *(g + i + 1) * 56448.0 - *(g + i + 2) * 7056.0 + *(g + i + 3) * 896. - *(g + i + 4) * 63.) / dx / dx / 35280.0;

        return y;
    }

    double dfdxxx(double *g, int i)
    {
        double y;
        // 3ord
        // y = (*(g + i + 2) - *(g + i + 1) * 3.0 + *(g + i) * 3.0 - *(g + i - 1)) / dx / dx / dx;
        // 5ord
        // y = (-*(g + i + 3) + *(g + i + 2) * 7.0 - *(g + i + 1) * 14.0 + *(g + i) * 10.0 - *(g + i - 1) - *(g + i - 2)) / dx / dx / dx / 4.0;
        // 6ord
        // y = (*(g + i - 3) - *(g + i - 2) * 8.0 + *(g + i - 1) * 13.0 + *(g + i) * 0.0 - *(g + i + 1) * 13.0 + *(g + i + 2) * 8.0 - *(g + i + 3)) / dx / dx / dx / 8.0;
        // 8ord
        y = (-*(g + i - 4) * 7. + *(g + i - 3) * 72. - *(g + i - 2) * 338.0 + *(g + i - 1) * 488.0 - *(g + i) * 0.0 - *(g + i + 1) * 488.0 + *(g + i + 2) * 338.0 - *(g + i + 3) * 72. + *(g + i + 4) * 7.) / dx / dx / 240.0;

        return y;
    }

    void TVabs(int m)
    {
        double TV = 0.0;
        TVmax = 0.0;

        // for (int i = nb; i < n1 - nb; i++) // 间断  int i = nb; i < n1 - nb; i++
        // {
        //     TV = TV + abs(u_nn[i] - 0.5 * sin(PI * x[i]) - 0.5);
        // }
        // TV = TV / (n1 - nb - nb);

        for (int i = nb; i < n1 - nb; i++) // 无穷  间断  int i = nb; i < n1 - nb; i++
        {
            TV = abs(u_nn[i] - 0.5 * sin(PI * x[i]) - 0.5);
            if (TV > TVmax)
            {
                TVmax = TV;
            }
        }
        TV = TVmax;

        cout << "\n  CFL: " << CFL << ",  TV: " << TV << "   \n";

        string Title = "error.txt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
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
    int n, oder_n = 6, begin_cell = 40;
    double t = 2.0, CFL = 0.5; // 计算时间总长，CFL数  CFL = 0.66 1.45
    // n = begin_cell + 1;
    for (int m = 1; m < oder_n; m++)
    {
        // CFL = CFL + 0.1;
        n = pow(2, m - 1) * begin_cell + 1;
        comput_main(n, t, CFL, m);
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

    double TL_n0[n + 6 * 2], TL_n1[n + 6 * 2], TL_n2[n + 6 * 2], TL_n3[n + 6 * 2], TL_n4[n + 6 * 2];
    double TLt_n0[n + 6 * 2], TLt_n1[n + 6 * 2], TLt_n2[n + 6 * 2], TLt_n3[n + 6 * 2], TLt_n4[n + 6 * 2];
    double Tu1[n + 6 * 2], Tu2[n + 6 * 2];
    class cfda cff;
    cff.kt = 1;
    cff.A1 = 1 / 3.0;
    cff.A2 = 2 / 3.0;
    cff.t = t;          // 1.0 / cff.PI; //计算总时间
    cff.nx_begin = 0.0; // 网格右端点
    cff.nx_end = 2.0;   // 网格左端点
    cff.CFL = CFL, cff.TVmax = 0.0;

    cff.n = n, cff.nb = nb, cff.n1 = n1;
    cff.Ltt_n0 = Ltt_n0, cff.x = x; // 声明变长数组

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
    cff.df = main_df, cff.u_excat = u_excat,
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
    // double th = cff1.dt + 1E-10;
    cff1.tt = 0.0;
    n_all = cff1.n1;
    cff1.border();
    cff1.f_eq_u(cff1.f_n0, cff1.u_n0);
    cff1.compt_CFL();

    for (int i1 = 0; i1 < 20000000; i1++) // 200 50
    {
        cff1.tt = cff1.tt + cff1.dt;
        cff1.compt_CFL();
        // cout << " \n i1 = " << i1 ;
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
            cff1.compt_2_stage_int(); // cff1.compt_SGLM_int();
            cff1.dt = cff1.dt * 0.5;

            // cff1.compt_ThD_13_line_SSP();
            // cff1.compt_ThD_24_line_SSP();
            // cff1.compt_ThD_25_line_SSP();
            // cff1.compt_ThD_36_line_SSP();

            // cff1.compt23_LLttnn_line();
            // cff1.compt24_LLttnn_line();
            // cff1.compt35_LLttnn_line();

            // cff1.RK3_compt_t1(); // 3-RK
            // cff1.RK33_compt(); // 3-RK
            // cff1.SSPRK54_compt_Shu(); //// // cff1.RK54LS_compt();
            cff1.RK65LS_compt2(); // //// cff1.RK65LS_compt();

            // cff1.RK76LS_compt(); // cff1.RK22_compt(); // cff1.compt34_LLttnn_line();

            if (i1 < 1)
                cff1.dt = 0.0;
        }
        else
        {
            // cff1.compt_TDTS23_LLttnn_line();
            // cff1.compt_TDTS24_LLttnn_line();
            cff1.compt_TDTS25_LLttnn_line();

            // cff1.compt_TDTS35_LLttnn_line();

            // cff1.compt_SGLM3_line();
            // cff1.compt_SGLM5_line();
        }
        cff1.carr(cff1.u_n0, cff1.u_nn, n_all);
        cff1.f_eq_u(cff1.f_n0, cff1.u_nn);
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
