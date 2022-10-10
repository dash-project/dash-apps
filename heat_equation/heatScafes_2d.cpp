#include "ScaFES.hpp"
#include "heatScafes_2d.hpp"

#include "minimon.h"

/** Space dimension of problem. */
const int DIM = 2;

/** Main program for HeatEqnFDM. */
int main(int argc, char *argv[]) {
    MiniMon minimon{};
    minimon.enter();
    ScaFES::Parameters paramsCl(argc, argv);
    ScaFES::GridGlobal<DIM> gg(paramsCl);

    std::vector<std::string> nameDatafield(1);
    nameDatafield[0] = "F";
    std::vector<int> stencilWidth(1);
    stencilWidth[0] = 1;
    std::vector<bool> isKnownDf(1);
    isKnownDf[0] = false;
    std::vector<int> nLayers(1);
    nLayers[0] = 0;
    std::vector<double> defaultValue(1);
    defaultValue[0] = 0.0;
    std::vector<ScaFES::WriteHowOften> writeToFile(1, ScaFES::WriteHowOften::LIKE_GIVEN_AT_CL);
    std::vector<bool> computeError(1);
    computeError[0] = false;

    std::vector<double> geomparamsInit;

    HeatEqnFDM<double, DIM> ppp(paramsCl, gg, false, nameDatafield, stencilWidth,
                                isKnownDf, nLayers, defaultValue, writeToFile,
                                computeError, geomparamsInit);
    ppp.iterateOverTime();
    minimon.leave("total");

    minimon.print(0);

    return 0;
}
