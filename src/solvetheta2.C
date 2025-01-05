#include <cmath>
#include <iostream>

double solveTheta2(double L1, double t2, double t3, double theta1, double deltaL, double initialTheta2 = 0.0, double tolerance = 1e-6, int maxIterations = 100) {
    double theta2 = initialTheta2;

    for (int iter = 0; iter < maxIterations; ++iter) {
        // 関数 f(θ2) を計算
        double x_diff = L1 * cos(theta1) + t3 * cos(theta2) - t2 * cos(theta1);
        double y_diff = L1 * sin(theta1) + t3 * sin(theta2) - t2 * sin(theta1);
        double f = x_diff * x_diff + y_diff * y_diff - deltaL * deltaL;

        // f'(θ2) を計算
        double dfdtheta2 = 2 * (x_diff * (-t3 * sin(theta2)) + y_diff * (t3 * cos(theta2)));

        // 更新
        double theta2_next = theta2 - f / dfdtheta2;

        // 収束判定
        if (std::abs(f) < tolerance) {
            return theta2_next;
        }

        theta2 = theta2_next;
    }

    // 収束しなかった場合
    std::cerr << "Newton's method did not converge." << std::endl;
    return theta2;
}
