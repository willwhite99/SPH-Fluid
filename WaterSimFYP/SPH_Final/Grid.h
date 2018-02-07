#pragma once

#include "Particle.h"

struct Cell
{
	int X;
	int Y;
	bool XCorner;
	bool YCorner;
	std::vector<Particle*> Particles;
};

class Grid
{
public:
	Grid(int width, int height, int cWidth, int cHeight);
	Cell* getCell(int x, int y) const { return &Cells[x][y]; }
	void populate(Particle* particles, int size);
	void clear();

	int getCellX() const { return CellsX; }
	int getCellY() const { return CellsY; }
private:
	Cell** Cells;

	float CellSizeWidth;
	float CellSizeHeight;
	int CellsX;
	int CellsY;
};