/*
	Demo program for the OpenGL Graph component
	Author: Gabriyel Wong (gabriyel@gmail.com)
*/
#define _USE_MATH_DEFINES
#include <windows.h>   // Standard Header For Most Programs
#include <math.h>
#include "GL/freeglut.h"   // The GL Utility Toolkit (Glut) Header
#include "OGLGraph.h"
#include <vector>

OGLGraph* myGraph;
float last_time = 0.f;
int frame = 0;
std::vector<float> temp;

int nx = 101;
float dx = 2.f/(nx-1);
int nt = 20;
float c = 1.f;

float vis = 0.2f;
float dt = 0.2f*dx*dx/vis;
bool enabled = false;

void computeBurgers(std::vector<float>& tempe)
{

	std::vector<float> orgCopy(tempe);
	for (int i = 1; i < nx - 1; i++)
	{
		// part deriv of position in resp to time + non-linear convection = diffsuion
		// rearranged to get position
		tempe[i] = orgCopy[i] - orgCopy[i] * dt / dx*(orgCopy[i] - orgCopy[i - 1]) + (vis * dt / (dx*dx))*(orgCopy[i + 1] - (2 * orgCopy[i]) + orgCopy[i - 1]);
		//tempe[i] = orgCopy[i] -c * dt / dx*(orgCopy[i] - orgCopy[i - 1]); linear convection
	}
	int last = nx - 1;
	tempe[last] = orgCopy[last] - orgCopy[last] * dt / dx*(orgCopy[last] - orgCopy[last - 1]) + (vis * dt / (dx*dx))*(orgCopy[0] - (2 * orgCopy[last]) + orgCopy[last - 1]);
}

float computeConditionals(float t, int x, float nu)
{
	float phi = exp(-pow((x - 4 * t), 2) / (4 * nu*(t + 1))) + exp(-pow((x - 4 * t - 2 * M_PI), 2) / (4 * nu*(t + 1)));
	float phiprime = -(-8 * t + 2 * x)*exp(-pow((-4 * t + x), 2) / (4 * nu*(t + 1))) / (4 * nu*(t + 1)) - (-8 * t + 2 * x - 12.5663706143592)*
		exp(-pow((-4 * t + x - 6.28318530717959), 2) / (4 * nu*(t + 1))) / (4 * nu*(t + 1));

	return -2 * nu*(phiprime / phi) + 4;
}

void init ( GLvoid )     // Create Some Everyday Functions
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.f);				// Black Background
	//glClearDepth(1.0f);								// Depth Buffer Setup
	myGraph = OGLGraph::Instance();
	myGraph->setup( 500, 100, 10, 10, 3, 20, 1, 235 );

	for (int i = 0; i < nx; i++)
	{
		temp.push_back((i > 40 && i < 50) ? 2 : 1);
		//temp.push_back(computeConditionals(0, i, vis));
	}
}

void display ( void )   // Create The Display Function
{
	frame++;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (float) 800 / (float) 600, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	if (enabled)
		computeBurgers(temp);

	myGraph->clear();
	for (int i = 0; i < temp.size(); i++)
		myGraph->update(temp[i]);
	myGraph->draw();

	glutSwapBuffers ( );
}

void reshape ( int w, int h )   // Create The Reshape Function (the viewport)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (float) w / (float) h, 0.1, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, w, h);
}

void keyboard ( unsigned char key, int x, int y )  // Create Keyboard Function
{
	switch ( key ) 
	{
	case 27:        // When Escape Is Pressed...
		exit ( 0 );   // Exit The Program
		break;        // Ready For Next Case
	case ' ':
		enabled = !enabled;
		break;
	default:        // Now Wrap It Up
		break;
	}
}

void idle(void)
{
  glutPostRedisplay();
}

void main ( int argc, char** argv )   // Create Main Function For Bringing It All Together
{
  glutInit( &argc, argv ); // Erm Just Write It =)
  init();

  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE ); // Display Mode
  glutInitWindowSize( 500, 250 ); // If glutFullScreen wasn't called this is the window size
  glutCreateWindow( "OpenGL Graph Component" ); // Window Title (argv[0] for current directory as title)
  glutDisplayFunc( display );  // Matching Earlier Functions To Their Counterparts
  glutReshapeFunc( reshape );
  glutKeyboardFunc( keyboard );
  glutIdleFunc(idle);
  glutMainLoop( );          // Initialize The Main Loop
}

