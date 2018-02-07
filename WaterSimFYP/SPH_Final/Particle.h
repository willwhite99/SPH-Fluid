#pragma once

struct Cell;

#include "navier_impl_helper.h"
#include "Grid.h"
#include <vector>

struct Particle
{
	Vector Position;
	Vector Velocity;
	Vector Force;
	float Rho, Pressure;
	Cell* Cell;
	std::vector<Particle*> Neighbours;
};