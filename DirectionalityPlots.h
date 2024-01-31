#include "NeutrinoDirectionality.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TMarker.h"

using std::array;

std::string plotDirectory;

struct FinalValues
{
    float phiTrue, phiTrueError, thetaTrue, thetaTrueError;
    array<float, DatasetSize> phi;
    array<float, DatasetSize> phiError;
    array<float, DatasetSize> phiErrorSystematics;
    array<float, DatasetSize> theta;
    array<float, DatasetSize> thetaError;
    array<float, DatasetSize> thetaErrorSystematics;
    array<float, DatasetSize> tilt;
    array<float, DatasetSize> tiltSystematics;
};