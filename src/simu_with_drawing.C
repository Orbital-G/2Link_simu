#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <sstream>
#include <Eigen/Dense>

using namespace Eigen;

class RobotArm {
public:
    double l1, l2;         // リンクの長さ
    double springOffset;   // バネ1とリンク1固定箇所の距離
    double theta1, theta2; // 関節角度
    double offsetDistance; // バネリンク接続点のオフセット量

    RobotArm(double l1, double l2, double springOffset, double offsetDistance)
        : l1(l1), l2(l2), springOffset(springOffset), theta1(M_PI / -1.5), theta2(M_PI / -6),
          offsetDistance(offsetDistance) {}

    void computeSpring1Position(Vector3d& start, Vector3d& end) {
        start << 400.0 - springOffset, 300.0, 0.0;
        double x0 = 400, y0 = 300;
        double linkX = x0 + (l1 * 0.5) * cos(theta1) * 100;
        double linkY = y0 - (l1 * 0.5) * sin(theta1) * 100;
        Vector3d linkMid(linkX, linkY, 0.0);

        double offsetX = -offsetDistance * sin(theta1);
        double offsetY = offsetDistance * cos(theta1);
        Vector3d offsetVec(offsetX, offsetY, 0.0);

        end = linkMid + offsetVec;
    }

    void computeSpring2Position(Vector3d& start, Vector3d& end) {
        double x0 = 400, y0 = 300;
        double link1X = x0 + (l1 * 0.25) * cos(theta1) * 100;
        double link1Y = y0 - (l1 * 0.25) * sin(theta1) * 100;
        start << link1X, link1Y, 0.0;

        double x1 = x0 + l1 * cos(theta1) * 100;
        double y1 = y0 - l1 * sin(theta1) * 100;
        double link2X = x1 + (l2 * 0.8) * cos(theta1 + theta2) * 100;
        double link2Y = y1 - (l2 * 0.8) * sin(theta1 + theta2) * 100;
        end << link2X, link2Y, 0.0;
    }

    void computeLinkPositions(Vector3d& joint1, Vector3d& joint2) {
        double x0 = 400, y0 = 300;
        joint1 << x0 + l1 * cos(theta1) * 100,
                  y0 - l1 * sin(theta1) * 100,
                  0.0;
        joint2 << joint1(0) + l2 * cos(theta1 + theta2) * 100,
                  joint1(1) - l2 * sin(theta1 + theta2) * 100,
                  0.0;
    }

    double computeSpringLength(const Vector3d& start, const Vector3d& end) {
        return (end - start).norm(); 
    }

    bool computeForceEllipse(double& a, double& b, double& angle) {
        double J11 = -l1 * sin(theta1) - l2 * sin(theta1 + theta2);
        double J12 = -l2 * sin(theta1 + theta2);
        double J21 = l1 * cos(theta1) + l2 * cos(theta1 + theta2);
        double J22 = l2 * cos(theta1 + theta2);

        // ヤコビ行列の転置を掛けた行列
        Matrix2d JTJ;
        JTJ << J11 * J11 + J21 * J21, J11 * J12 + J21 * J22,
            J11 * J12 + J21 * J22, J12 * J12 + J22 * J22;

        // 固有値分解
        SelfAdjointEigenSolver<Matrix2d> solver(JTJ);
        if (solver.info() != Eigen::Success) {
            return false;
        }

        Vector2d eigenvalues = solver.eigenvalues();
        Matrix2d eigenvectors = solver.eigenvectors();

        // 固有値
        double lambda1 = eigenvalues(1); // 最大固有値
        double lambda2 = eigenvalues(0); // 最小固有値

        // 固有ベクトル（長軸方向）
        Vector2d singularVector = eigenvectors.col(1);

        // リンク2の方向ベクトル
        Vector2d link2Direction(std::cos(theta1 + theta2), std::sin(theta1 + theta2));

        // 特異値ベクトルの符号をリンク2方向に合わせる
        if (singularVector.dot(link2Direction) < 0) {
            singularVector = -singularVector;
        }

        // 相対角度を計算
        double relativeAngle = std::atan2(singularVector(1), singularVector(0)) -
                            std::atan2(link2Direction(1), link2Direction(0));

        // 楕円体の長軸・短軸
        a = 1.0 / std::sqrt(lambda1) * 30; // スケール調整
        b = 1.0 / std::sqrt(lambda2) * 30;

        // 楕円体の角度をリンク2基準に設定
        angle = relativeAngle * 180 / M_PI;

        return true;
    }


    double computeForceInNegativeYDirection() {
        // ヤコビ行列の計算
        double J11 = -l1 * sin(theta1) - l2 * sin(theta1 + theta2);
        double J12 = -l2 * sin(theta1 + theta2);
        double J21 = l1 * cos(theta1) + l2 * cos(theta1 + theta2);
        double J22 = l2 * cos(theta1 + theta2);

        Matrix2d J;
        J << J11, J12,
             J21, J22;

        // J * J^T の計算
        Matrix2d JJ_T = J * J.transpose();

        // 特異値分解
        SelfAdjointEigenSolver<Matrix2d> eigensolver(JJ_T);
        if (eigensolver.info() != Success) return 0.0;

        Vector2d singularValues = eigensolver.eigenvalues().cwiseSqrt();
        Matrix2d singularVectors = eigensolver.eigenvectors();

        // -y方向（0, -1）と特異値ベクトルとのドット積を使った計算
        Vector2d negativeYDirection(0, -1);

        double maxForce = 0.0;
        for (int i = 0; i < 2; ++i) {
            Vector2d singularVector = singularVectors.col(i);
            double projection = std::abs(singularVector.dot(negativeYDirection)) * singularValues(i);
            if (projection > maxForce) {
                maxForce = projection;
            }
        }

        return maxForce;
    }


};

void drawLink(sf::RenderWindow& window, const Vector3d& start, const Vector3d& end, sf::Color color, float thickness) {
    Vector3d diff = end - start;
    float length = diff.head<2>().norm();

    sf::RectangleShape link(sf::Vector2f(length, thickness));
    link.setOrigin(0, thickness / 2);
    link.setPosition(start(0), start(1));
    link.setRotation(std::atan2(diff(1), diff(0)) * 180 / M_PI);
    link.setFillColor(color);

    window.draw(link);
}

int main() {
    sf::RenderWindow window(sf::VideoMode(800, 600), "2-Link Robot Arm with Springs");
    RobotArm arm(1.5, 1.0, 150, 10);
    const double angleStep = 0.0005;

    sf::Font font;
        if (!font.loadFromFile("/Library/Fonts/Arial Unicode.ttf")) {
            std::cerr << "Error: Could not load font." << std::endl;
            return -1;
        }


    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
        }

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) arm.theta1 += angleStep;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) arm.theta1 -= angleStep;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) arm.theta2 += angleStep;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) arm.theta2 -= angleStep;

        window.clear(sf::Color::White);

        Vector3d spring1Start, spring1End;
        arm.computeSpring1Position(spring1Start, spring1End);
        drawLink(window, spring1Start, spring1End, sf::Color::Blue, 5);

        Vector3d spring2Start, spring2End;
        arm.computeSpring2Position(spring2Start, spring2End);
        drawLink(window, spring2Start, spring2End, sf::Color::Magenta, 5);

        Vector3d joint1, joint2;
        arm.computeLinkPositions(joint1, joint2);
        drawLink(window, Vector3d(400.0, 300.0, 0.0), joint1, sf::Color::Black, 10);
        drawLink(window, joint1, joint2, sf::Color::Black, 10);

        double a, b, angle;
        if (arm.computeForceEllipse(a, b, angle)) {
            sf::CircleShape forceEllipse(1.0);
            forceEllipse.setRadius(1.0);
            forceEllipse.setScale(a, b);
            forceEllipse.setOrigin(1.0, 1.0);
            forceEllipse.setPosition(joint2(0), joint2(1));
            forceEllipse.setRotation(angle);
            forceEllipse.setFillColor(sf::Color(0, 255, 100, 100));

            window.draw(forceEllipse);
        }

        // -y方向に最も出せる力を計算
        double maxForceNegativeY = arm.computeForceInNegativeYDirection();

        // テキスト描画
        sf::Text text;
        text.setFont(font);
        text.setCharacterSize(20);
        text.setFillColor(sf::Color::Black);

        std::ostringstream oss;
        oss << "Theta1: " << arm.theta1 * 180 / M_PI << " degrees\n";
        oss << "Theta2: " << arm.theta2 * 180 / M_PI << " degrees\n";

        double link1Length = arm.computeSpringLength(Vector3d(400.0, 300.0, 0.0), joint1);
        double link2Length = arm.computeSpringLength(joint1, joint2);

        oss << "Link1 Length: " << link1Length << "\n";
        oss << "Link2 Length: " << link2Length << "\n";
        oss << "Spring1 Length: " << arm.computeSpringLength(spring1Start, spring1End) << "\n";
        oss << "Spring2 Length: " << arm.computeSpringLength(spring2Start, spring2End) << "\n";

        oss << "Max Force in -Y Direction: " << maxForceNegativeY << "\n"; // -y方向の力を追加

        text.setString(oss.str());
        text.setPosition(10, 10);
        window.draw(text);


        window.display();
    }

    return 0;
}
