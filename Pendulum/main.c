#include <stdio.h>
#include <math.h>
#include "raylib.h"

#define WIN_WIDTH 640
#define WIN_HEIGHT 480

void CalculateAlpha(
	double theta1, double theta2, double l1, double l2,
	double m1, double m2, double* alpha1, double* alpha2
)
{
	double cosThetaDiff = cos(theta1 - theta2);
	*alpha1 = (l2 / l1) * (m2 / (m1 + m2)) * cosThetaDiff;
	*alpha2 = (l1 / l2) * cosThetaDiff;
}

void CalculateF(
	double theta1, double theta2, double omega1, double omega2,
	double l1, double l2, double g, double m1, double m2, double* f1, double* f2)
{
	double sinThetaDiff = sin(theta1 - theta2);
	*f1 = -((l2 / l1) * (m2 / (m1 + m2)) * omega2 * omega2 * sinThetaDiff) - (g / l1) * sin(theta1);
	*f2 = (l1 / l2) * omega1 * omega1 * sinThetaDiff - (g / l2) * sin(theta2);
}

void CalculateG(
	double alpha1, double alpha2, double f1, double f2, double* g1, double* g2
) 
{
	double denominator = 1 - (alpha1 * alpha2);

	// Ensuring the denominator is not zero to avoid division by zero.
	if (fabs(denominator) > 1e-6) {
		*g1 = (f1 - (alpha1 * f2)) / denominator;
		*g2 = ((-alpha2 * f1) + f2) / denominator;
	}
	else {
		// Handle the special case where denominator is zero or very close to zero.
		*g1 = 0;
		*g2 = 0;
	}
}

void ComputeDerivatives(
	double theta1, double theta2, double omega1, double omega2,
	double l1, double l2, double m1, double m2, double g,
	double* g1, double* g2) 
{
	double alpha1, alpha2, f1, f2;

	CalculateAlpha(theta1, theta2, l1, l2, m1, m2, &alpha1, &alpha2);
	CalculateF(theta1, theta2, omega1, omega2, l1, l2, g, m1, m2, &f1, &f2);
	CalculateG(alpha1, alpha2, f1, f2, g1, g2);
}

void RungeKuttaStep(
	double* theta1, double* theta2, double* omega1, double* omega2,
	double dt, double l1, double l2, double g, double m1, double m2) 
{
	double k1_theta1, k1_theta2, k1_omega1, k1_omega2;
	double k2_theta1, k2_theta2, k2_omega1, k2_omega2;
	double k3_theta1, k3_theta2, k3_omega1, k3_omega2;
	double k4_theta1, k4_theta2, k4_omega1, k4_omega2;

	// Compute k1 values
	ComputeDerivatives(
		*theta1, *theta2, *omega1, *omega2,
		l1, l2, g, m1, m2, &k1_omega1, &k1_omega2
	);
	k1_theta1 = *omega1;
	k1_theta2 = *omega2;

	// Compute k2 values
	ComputeDerivatives(
		*theta1 + k1_theta1 * dt / 2, *theta2 + k1_theta2 * dt / 2,
		*omega1 + k1_omega1 * dt / 2, *omega2 + k1_omega2 * dt / 2,
		l1, l2, g, m1, m2, &k2_omega1, &k2_omega2
	);
	k2_theta1 = *omega1 + k1_omega1 * dt / 2;
	k2_theta2 = *omega2 + k1_omega2 * dt / 2;

	// Compute k3 values
	ComputeDerivatives(
		*theta1 + k2_theta1 * dt / 2, *theta2 + k2_theta2 * dt / 2,
		*omega1 + k2_omega1 * dt / 2, *omega2 + k2_omega2 * dt / 2,
		l1, l2, g, m1, m2, &k3_omega1, &k3_omega2
	);
	k3_theta1 = *omega1 + k2_omega1 * dt / 2;
	k3_theta2 = *omega2 + k2_omega2 * dt / 2;

	// Compute k4 values
	ComputeDerivatives(
		*theta1 + k3_theta1 * dt, *theta2 + k3_theta2 * dt,
		*omega1 + k3_omega1 * dt, *omega2 + k3_omega2 * dt,
		l1, l2, g, m1, m2, &k4_omega1, &k4_omega2
	);
	k4_theta1 = *omega1 + k3_omega1 * dt;
	k4_theta2 = *omega2 + k3_omega2 * dt;

	// Update state
	*theta1 += (k1_theta1 + 2 * k2_theta1 + 2 * k3_theta1 + k4_theta1) * dt / 6;
	*theta2 += (k1_theta2 + 2 * k2_theta2 + 2 * k3_theta2 + k4_theta2) * dt / 6;
	*omega1 += (k1_omega1 + 2 * k2_omega1 + 2 * k3_omega1 + k4_omega1) * dt / 6;
	*omega2 += (k1_omega2 + 2 * k2_omega2 + 2 * k3_omega2 + k4_omega2) * dt / 6;
}

int main()
{
	InitWindow(WIN_WIDTH, WIN_HEIGHT, "Double Pendulum Simulation");

	double theta1 = PI / 2;	// Initial angle for pendulum 1
	double theta2 = PI / 2; // Initial angle for pendulum 2
	double omega1 = 0;      // Initial angular velocity for pendulum 1
	double omega2 = 0;      // Initial angular velocity for pendulum 2
	double l1 = 100;        // Length of pendulum 1
	double l2 = 100;        // Length of pendulum 2
	double m1 = 1;          // Mass of pendulum 1
	double m2 = 1;          // Mass of pendulum 2
	double g = 9.80665;     // Acceleration due to gravity
	double dt = 0.1;        // Time step
	Vector2 origin = { WIN_WIDTH / 2, WIN_HEIGHT / 2 };

	SetTargetFPS(60);

	while (!WindowShouldClose())
	{
		RungeKuttaStep(&theta1, &theta2, &omega1, &omega2, dt, l1, l2, g, m1, m2);

		// Calculate screen coordinates for pendulum parts
		Vector2 pendulum1End = { origin.x + l1 * sin(theta1), origin.y + l1 * cos(theta1) };
		Vector2 pendulum2End = { pendulum1End.x + l2 * sin(theta2), pendulum1End.y + l2 * cos(theta2) };

		BeginDrawing();

		DrawLineV(origin, pendulum1End, BLACK);			// Draw first arm
		DrawCircleV(pendulum1End, 10, RED);             // Draw weight of first arm
		DrawLineV(pendulum1End, pendulum2End, BLACK);	// Draw second arm
		DrawCircleV(pendulum2End, 10, BLUE);            // Draw weight of second arm

		ClearBackground(LIGHTGRAY);

		EndDrawing();
	}

	CloseWindow();

	return 0;
}