#pragma once

#include "../include/glm/glm.hpp"
#include <math.h>

typedef glm::vec2 Vector;

#define MAX(a, b) (((a) > (b))? (a): (b))
#define MIN(a, b) (((a) > (b))? (b): (a))

static Vector randomVector(float radius)
{
	float u = rand() / (float)RAND_MAX;
	float v = 1.f - (rand() / (float)RAND_MAX) * 2;

	float incl = acosf(v);
	float azi = 2 * 3.14f * u;

	float rx = radius * sinf(incl) * cos(azi);
	float ry = radius * sinf(incl) * sin(azi);

	return Vector(rx, ry);
}