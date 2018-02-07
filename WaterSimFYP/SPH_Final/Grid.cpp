#include "Grid.h"
#include <omp.h>

Grid::Grid(int width, int height, int cWidth, int cHeight)
{
	CellSizeWidth = width / (float)cWidth;
	CellSizeHeight = height / (float)cHeight;

	CellsX = cWidth;
	CellsY = cHeight;

	Cells = new Cell*[cWidth];
	for (int i = 0; i < cWidth; i++)
		Cells[i] = new Cell[cHeight];

	// initialise grid cells
	for (int i = 0; i < cWidth; i++)
		for (int j = 0; j < cHeight; j++)
		{
			Cells[i][j].X = i;
			Cells[i][j].Y = j;

			Cells[i][j].XCorner = false;
			Cells[i][j].YCorner = false;

			if (i == 0 || i == cWidth - 1)
				Cells[i][j].XCorner = true;
			if (j == 0 || j == cHeight - 1)
				Cells[i][j].YCorner = true;
		}
}

void Grid::populate(Particle* particles, int size)
{
	// add particles to their cells
	for (int i = 0; i < size; i++)
	{
		int x = (int)floor(particles[i].Position.x / CellSizeWidth);
		int y = (int)floor(particles[i].Position.y / CellSizeHeight);

		x = MIN(x, CellsX - 1);
		y = MIN(y, CellsY - 1);
		x = MAX(x, 0);
		y = MAX(y, 0);

		Cells[x][y].Particles.push_back(&particles[i]);
		// add cell to particle, so we can find its neighbours
		particles[i].Cell = &Cells[x][y];
	}
}

void Grid::clear()
{
#pragma omp parallel for
	for (int i = 0; i < CellsX; i++)
		for (int j = 0; j < CellsY; j++)
			Cells[i][j].Particles.clear();
}