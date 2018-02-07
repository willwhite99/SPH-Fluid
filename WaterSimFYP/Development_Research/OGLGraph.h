/*
	Header file for OGLGraph class
	Copyright (C) 2006 Gabriyel Wong (gabriyel@gmail.com)

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#include <time.h>
#include <List>
#include "GL/freeglut.h"

class OGLGraph
{
public:
	static OGLGraph* Instance();
	~OGLGraph() {}
	void setup( int, int, int, int, int, int, int, int );
	void update( float );
	void draw();
	std::list<float> getData() const { return _data; }
	void clear() { _data.clear(); }
protected:
//	OGLGraph();

private:
	static OGLGraph* _instance;
	OGLGraph() {}
	int _cacheSize;
	int _pWidth, _pHeight;	// in percentages of existing window sizes
	int _width, _height;
	int _channels;
	int _scaleX, _scaleY;
	int _offsetX, _offsetY;
	std::list<float> _data;
	std::list<float> _movingAverage;
	std::list<float>::const_iterator _iterData;
	std::list<float>::const_iterator _iterMA;
};
