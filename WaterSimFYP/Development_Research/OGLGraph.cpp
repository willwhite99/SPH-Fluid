/*
	Implementation of the OGLGraph class functions
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
#include "OGLGraph.h"

OGLGraph* OGLGraph::_instance = 0;

OGLGraph* OGLGraph::Instance()
{
	if (_instance == 0) 
	{
		_instance = new OGLGraph();
    }
    return _instance;
}

void
OGLGraph::setup( int width, int height, int offsetX, int offsetY, 
				 int scaleX, int scaleY, int channels, int cacheSize )
{
	_offsetX = offsetX;
	_offsetY = offsetY;
	_width = width;
	_height = height;
	_scaleX = scaleX;
	_scaleY = scaleY;
	_channels = channels;	// ToDo
	_cacheSize = cacheSize;
}

void
OGLGraph::update( float data )
{
	if ( _data.size() > (unsigned int)(_cacheSize) )
	{ 
		return;
		//_data.pop_front();
		//_movingAverage.pop_front();
		//
		//std::list<float>::const_iterator it;

		//// Computes moving average
		//float subTotal = 0.f;
		//_data.push_back( data );
		//for ( it=_data.begin(); it != _data.end(); it++ )
		//{
		//	subTotal += *it;
		//}
		//float maTemp = (float) (subTotal/(_data.size()));
		//_movingAverage.push_back( maTemp );

	}
	else
	{
        _data.push_back( data );

		// Computes moving average
		std::list<float>::const_iterator it;
		float subTotal = 0.f;
		for ( it=_data.begin(); it != _data.end(); it++ )
		{
			subTotal += *it;
		}
		float maTemp = (float) (subTotal/(_data.size()));
		_movingAverage.push_back( maTemp );
	}
}

void
OGLGraph::draw()
{
	int cntX = 0;
	int cntY = 0;
	// Set up the display
	glMatrixMode(GL_PROJECTION);  //We'll talk about this one more as we go 
	glLoadIdentity(); 
	glOrtho(0, _width, 0, _height, 0, 1.0); 

	glMatrixMode(GL_MODELVIEW); 
	glLoadIdentity();

	// Draw the axes - OpenGL screen coordinates start from bottom-left (0,0)
	glDisable (GL_LINE_STIPPLE);
	glBegin( GL_LINES );
	glColor3ub( 255, 255, 255 );	// Green
	// Draw x-axis
	glVertex3f( _offsetX, _offsetY, 0 );
	glVertex3f( _width-_offsetX, _offsetY, 0 );
	// Draw y-axis
	glVertex3f( _offsetX, _offsetY, 0 );
	glVertex3f( _offsetX, 50, 0 );
	glEnd();

	// Draw the data points
	glBegin( GL_LINE_STRIP );
	glColor3ub( 0, 255, 0 );	// Green
	for ( _iterData=_data.begin(); _iterData != _data.end(); _iterData++ )
	{
        glVertex3f( 10+(cntX*_scaleX), ((*_iterData))*_scaleY, 0 );
		cntX++;
	}
	glEnd();

	// Draw the moving average values
	glEnable (GL_LINE_STIPPLE);
	glLineStipple (3, 0xAAAA);
	glBegin( GL_LINE_STRIP );
	glColor3ub( 255, 0, 0);	// Red
	// Draw moving average graph
	//for (_iterMA=_movingAverage.begin(); _iterMA != _movingAverage.end(); _iterMA++)
	//{
 //       glVertex3f( 10+(cntY*_scaleX), (10+(*_iterMA))*_scaleY, 0);
	//	cntY++;
	//}
	glEnd();
}
