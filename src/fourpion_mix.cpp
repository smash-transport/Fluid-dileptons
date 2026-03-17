#include <cmath>
#include <algorithm>
#include <iostream>

#include "dilepton.h"
#include "setup.h"

namespace FluidDileptons {

namespace Rates {

    constexpr double pion_mass = 0.139;  // pion mass in GeV
    constexpr double tadpole_at_Tc = 0.49513467;  // tadpi(T_c = 0.175 GeV)

    // Tabulated pion tadpole in temperature grid [0.050, 0.200] with 0.001 GeV steps
    static const std::vector<double> tadpole_tab = {
        0.052484015, 0.055966439, 0.0595307621, 0.0631723997, 0.0668868373,
        0.0706695309, 0.074516008, 0.0784218749, 0.0823827521, 0.0863944054,
        0.0904526758, 0.0945535052, 0.0986929449, 0.1028671600, 0.1070724190,
        0.1113051630, 0.1155619050, 0.1198393050, 0.1241341670, 0.1284433720,
        0.1327639820, 0.1370931680, 0.1414282320, 0.1457665990, 0.1501058190,
        0.1544435520, 0.1587776230, 0.1631059280, 0.1674264950, 0.1717374690,
        0.1760371100, 0.1803237510, 0.1845958700, 0.1888520320, 0.1930909010,
        0.1973112360, 0.2015118890, 0.2056917890, 0.2098499840, 0.2139855730,
        0.2180977480, 0.2221857740, 0.2262489930, 0.2302868290, 0.2342987420,
        0.2382842890, 0.2422430800, 0.2461747870, 0.2500791290, 0.2539559060,
        0.2578049480, 0.2616261410, 0.2654194210, 0.2691847700, 0.2729222130,
        0.2766318180, 0.2803136990, 0.2839679880, 0.2875948710, 0.2911945660,
        0.2947673140, 0.2983134060, 0.3018331440, 0.3053268640, 0.3087949290,
        0.3122377250, 0.3156556630, 0.3190491750, 0.3224187150, 0.3257647570,
        0.3290877980, 0.3323883370, 0.3356669000, 0.3389240400, 0.3421603080,
        0.3453762740, 0.3485725200, 0.3517496430, 0.3549082490, 0.3580490380,
        0.3611724700, 0.3642792650, 0.3673700680, 0.3704455290, 0.3735063200,
        0.3765530970, 0.3795865380, 0.3826073250, 0.3856161430, 0.3886136840,
        0.3916005660, 0.3945776490, 0.3975455580, 0.4005050040, 0.4034567000,
        0.4064013600, 0.4093397110, 0.4122724720, 0.4152003700, 0.4181241380,
        0.4210445010, 0.4239621970, 0.4268779620, 0.4297925360, 0.4327066590,
        0.4356210740, 0.4385365290, 0.4414537700, 0.4443735480, 0.4472966130,
        0.4502237200, 0.4531556240, 0.4560930830, 0.4590368570, 0.4619877120,
        0.4649464070, 0.4679137090, 0.4708903920, 0.4738772250, 0.4768749820,
        0.4798844410, 0.4829063830, 0.4859415890, 0.4889908450, 0.4920549410,
        0.4951346700, 0.5027521360, 0.5104316360, 0.5181732020, 0.5259768660,
        0.5338426600, 0.5417706160, 0.5497607660, 0.5578131390, 0.5659277680,
        0.5741046810, 0.5823439070, 0.5906454740, 0.5990094140, 0.6074357520,
        0.6159245160, 0.6244757340, 0.6330894320, 0.6417656360, 0.6505043720,
        0.6593056650, 0.6681695420, 0.6770960260, 0.6860851410
    };

    // Interpolation function for pion tadpole diagram
    double tadpole_interpolate(double T) {
        // Linear interpolation in the tabulated data
        constexpr double Tmin = 0.05;
        constexpr double dT = 0.001;

        const int k = int((T - Tmin) / dT);
        if (k < 0)
            return tadpole_tab[0] * T / Tmin; // linear extrapolation to 0
        if (k >= tadpole_tab.size() - 1)
            return tadpole_tab[tadpole_tab.size() - 1]; // constant above Tmax

        const double t0 = Tmin + k*dT;
        const double p0 = tadpole_tab[k], p1 = tadpole_tab[k + 1];

        return p0 + (p1 - p0) * (T - t0) / dT;
    }

    // 3-pion axial-vector correlator (a1 → pi + rho removed)
    double Im_PiA3pi(double s) {
        double img_corr = 0;
        constexpr double s_threshold = 9.0 * pion_mass * pion_mass;

        if (s > s_threshold) {
            constexpr double A = 0.0305;
            constexpr double B = 47.692;
            constexpr double C = 1.1981;
            constexpr double D = 0.38559008;

            const double phase_supp = std::pow(1.0 - s_threshold * s_threshold / (s * s), B);
            const double denom = (s - C) * (s - C) + D;
            img_corr = (M_PI * A / s) * phase_supp / denom;
        }

        return img_corr;
    }

    // Vector correlator (4-pion contribution)
    double Im_PiV4pi(double s) {
        constexpr double s_threshold = 16.0 * pion_mass * pion_mass;

        double im_corr = 0;
        if (s > s_threshold) {
            constexpr double A = 0.22;
            constexpr double B = 0.2;
            constexpr double C = 1.15;
            constexpr double D = 0.11;

            const double sqrts = std::sqrt(s);
            const double phase_supp = std::pow(1.0 - s_threshold / s, 2.0);
            const double num = 1.0 + A / std::log(1.0 + sqrts / B);
            const double denom = 1.0 + std::exp((C - sqrts) / D);

            im_corr = 1 / (8.0 * M_PI) * phase_supp * num / denom;
        }
        return im_corr;
    }

    // 5-pion axial-vector correlator
    double Im_PiA5pi(double s) {
        double img_corr = 0;
        double s_threshold = 25.0 * pion_mass * pion_mass;

        if (s > s_threshold) {
            constexpr double A = 0.22;
            constexpr double B = 0.2;
            constexpr double C = 1.6;
            constexpr double D = 0.14;

            const double phase_supp = std::pow(1.0 - s_threshold / s, B);
            const double sqrts = std::sqrt(s);
            const double num = 1.0 + A / std::log(1.0 + sqrts / B);
            const double denom = 1.0 + std::exp((C - sqrts) / D);

            img_corr = 1.0 / (8.0 * M_PI) * phase_supp * num / denom;
        }
        return img_corr;
    }

    double ImD_multipi(const Parameters& par) {
        const double s = par.m * par.m, T = par.T;
        const double mixing =  tadpole_interpolate(T) / (2.0 * tadpole_at_Tc);
        double correlator = (1.0 - mixing) * Im_PiV4pi(s) + mixing * (0.5 * Im_PiA3pi(s) + Im_PiA5pi(s));
        return s*correlator/3.;
    }

    /*
     * Distribution function dR/dM (without constant factors)
     * not used and replaced by the actual McLerran-Toimela formula, but Stephen Endres used it.
     *   double dRds_multipi(const Parameters& par) {
     *      // WF = 3.0 * Im(PiVmix)
     *      // distm4pi_mix = mass^2 * T * BesselK1(mass/T) * WF
     *      const double s = par.m * par.m;
     *      return 3.0 * s * par.T * Im_PiVmix(s, par.T) * std::cyl_bessel_k(1, par.m / par.T);
     *  }
     */


} // namespace Rates

} // namespace FluidDileptons
