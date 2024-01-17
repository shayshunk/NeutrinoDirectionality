#include "GeneralHeader.h"

void CovarianceCalc()
{
    /*
    Sigma x: 0.281277
    Sigma y: 0.281054
    Sigma z: 0.281899
    Now adding systematics.
    Sigma x: 0.384729
    Sigma y: 0.487331
    Sigma z: 0.300112
    Sigma x unbiased: 0.449773
    Sigma y unbiased: 0.453999
    Sigma z unbiased: 0.281899
    Now adding systematics.
    Sigma x unbiased: 0.520764
    Sigma y unbiased: 0.603834
    Sigma z unbiased: 0.300112
    */

    /* 
    Mean x: 6.47142, Mean y: 5.43068, Mean z: -1.5185
    Mean Unbiased x: 8.19875, Mean Unbiased y: 6.95785, Mean Unbiased z: -1.5185
    */

    double p_values[2][3] = { {6.47142, 5.43068, -1.5185}, {8.19875, 6.95785, -1.519}}; // uncorrected then corrected
    double sigma_values[4][3] = { {0.281277, 0.281054, 0.281899}, {0.384729, 0.487331, 0.300112}, {0.449773, 0.453999, 0.281899}, {0.520764, 0.603834, 0.300112}}; // uncorrected and then corrected, stats then stats + sys each
    double angles[2][2] = { {40.0027, -10.1897}, {40.3195, -8.03783}};

    double px = p_values[0][0], py = p_values[0][1], pz = p_values[0][2];
    double sx = sigma_values[0][0], sy = sigma_values[0][1], sz = sigma_values[0][2];

    cout << "Sigma[0][1] is: " << sy << endl;


    double Cov[2][2];

    Cov[0][0] = ((sx * sx) * (py * py) / (pow(px, 4))) + (sy * sy) / (px * px);
    Cov[0][1] = ((py * pz) / (px * pow(px * px + py * py, (3./2.)))) * (sy * sy - sx * sx);
    Cov[1][0] = ((py * pz) / (px * pow(px * px + py * py, (3./2.)))) * (sy * sy - sx * sx); 
    Cov[1][1] = ((px * px * pz * pz * sx * sx + py * py * pz * pz * sy * sy)/ (pow(px * px + py * py, 3))) + (sz * sz) / (px * px + py * py);

    cout << "Without adding the angles:\n";
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            cout << "At row " << i << " and column " << j << ", we have: " << Cov[i][j] << "\n";
        }
    }

    Cov[0][0] = Cov[0][0] * pow(cos(angles[0][0] * pi/180), 4);
    Cov[0][1] = Cov[0][1] * pow(cos(angles[0][0] * pi/180), 2) * pow(cos(angles[0][1] * pi/180), 2);
    Cov[1][0] = Cov[1][0] * pow(cos(angles[0][0] * pi/180), 2) * pow(cos(angles[0][1] * pi/180), 2);
    Cov[1][1] = Cov[1][1] * pow(cos(angles[0][1] * pi/180), 4);

    cout << "After adding the angles:\n";
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            cout << "At row " << i << " and column " << j << ", we have: " << Cov[i][j] << "\n";
        }
    }

    double a = Cov[0][0], b = Cov[0][1], c = Cov[1][0], d = Cov[1][1];

    double lambda1 = ((a + d) + sqrt(pow(a, 2) - 2 * a * d + 4 * b * c + pow(d, 2))) / 2;
    double lambda2 = ((a + d) - sqrt(pow(a, 2) - 2 * a * d + 4 * b * c + pow(d, 2))) / 2;

    cout << "Eigenvalue 1: " << lambda1 << "\n";
    cout << "Eigenvalue 2: " << lambda2 << "\n";

    double phiErr = sqrt(2.291 * lambda1) * 180/pi;
    double thetaErr = sqrt(2.291 * lambda2) * 180/pi;

    cout << "Phi error: " << phiErr << " and theta error: " << thetaErr << "\n";

    float v1 = lambda1 - d;
    float v2 = c;

    float scaler = 1.0 / c;
    v1 = v1 * scaler;
    v2 = v2 * scaler;

    float tilt = atan(v1 / v2) * 180.0 / pi;
    tilt = 360 - tilt;

    cout << "v1: " << v1 << '\n';
    cout << "v2: " << v2 << '\n';
    cout << "Tilt: " << tilt << '\n';

}