#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GL/glew.h"
#include "GL/freeglut.h"

#define GLM_FORCE_RADIANS
#include "include/glm/glm.hpp"
#include "include/glm/gtc/matrix_transform.hpp"
#include "include/glm/gtc/type_ptr.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

GLuint window, vvbo, evbo, vao, program;
int x = 800, y = 600;

#define N 80
#define INITIAL_HEIGHT 40.f

GLfloat u[N][N];
GLfloat v[N][N];
GLfloat scale = 1.f;

GLfloat camera_x = 0.f;
GLfloat camera_z = 50.f;
GLfloat camera_rot = 0.f;

int c = 1;
float dx = 2.f / (N - 1);
float dy = dx;
float vis = 0.01f;
float dt = 0.0009f*dx*dy/vis;
bool enabled = false;
bool wireframe = false;

struct Vertex
{
	glm::vec3 pos;
};

void idle()
{
	glutPostRedisplay();
}

void reshape(int x, int y)
{
	glViewport(0, 0, x, y);
}

void display()
{
	if (enabled)
	{
		//rebuild graph
		GLfloat uCopy[N][N];
		GLfloat vCopy[N][N];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				uCopy[i][j] = u[i][j];
				vCopy[i][j] = v[i][j];
			}

		for (int i = 1; i < N; i++)
			for (int j = 1; j < N; j++)
			{
			// Burgers equation
			u[i][j] += -dt / dx*(uCopy[i][j] * (uCopy[i][j] - uCopy[i - 1][j])) - dt / dy*(vCopy[i][j] * (uCopy[i][j] - uCopy[i][j - 1])) +
				((vis * dt) / (dx*dx))*(uCopy[i + 1][j] - (2 * uCopy[i][j]) + uCopy[i - 1][j]) + ((vis * dt) / (dy*dy))*(uCopy[i][j + 1] - (2 * uCopy[i][j]) + uCopy[i][j - 1]);

			v[i][j] += -dt / dx*(uCopy[i][j] * (vCopy[i][j] - vCopy[i - 1][j])) - dt / dy*(vCopy[i][j] * (vCopy[i][j] - vCopy[i][j - 1])) +
				((vis * dt) / (dx*dx))*(vCopy[i + 1][j] - (2 * vCopy[i][j]) + vCopy[i - 1][j]) + ((vis * dt) / (dy*dy))*(vCopy[i][j + 1] - (2 * vCopy[i][j]) + vCopy[i][j - 1]);


			// ----------- Old equations, non-linear and diffusuion are sepearate here -----------
			// Non-linear convection
			//graph[i][j] = posCopy[i][j] - (posCopy[i][j] * (dt / dx)*(posCopy[i][j] - posCopy[i - 1][j])) - (vCopy[i][j] * (dt / dy)*(posCopy[i][j] - posCopy[i][j - 1]));
			//graphv[i][j] = vCopy[i][j] - (posCopy[i][j] * (dt / dx)*(vCopy[i][j] - vCopy[i - 1][j])) - (vCopy[i][j] * (dt / dy)*(vCopy[i][j] - vCopy[i][j - 1]));

			// Diffusion only
			//graph[i][j] += ((vis * dt) / (dx*dx))*(posCopy[i + 1][j] - (2 * posCopy[i][j]) + posCopy[i - 1][j]) + ((vis * dt) / (dy*dy))*(posCopy[i][j + 1] - (2 * posCopy[i][j]) + posCopy[i][j - 1]);
			}

		for (int z = 0; z < N; z++)
		{
			u[0][z] = 1.f;
			u[(N - 1)][z] = 1.f;
			u[z][0] = 1.f;
			u[z][(N - 1)] = 1.f;

			v[0][z] = 1;
			v[(N - 1)][z] = 1;
			v[z][0] = 1;
			v[z][(N - 1)] = 1;
		}
	}

	glClearColor(0.f, 0.f, 0.55f, 0.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(75.f, x / (float)y, 0.1f, 1000.f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(camera_x, 60.f, camera_z,
		0, 20, 0, 0, 1, 0);

	glScalef(scale, scale, scale);

	// draw surface
	glBegin(GL_TRIANGLES);
	float start = -(float)N / 2;
	for (int n = 0; n < N; n++)
	{
		for (int z = 0; z < N; z++)
		{
			// colours
			float nz = u[n][z] / 255.f;
			float nz1 = u[n][z + 1] / 255.f;
			float n1z = u[n+ 1][z] / 255.f;
			float n1z1 = u[n+1][z+1] / 255.f;

			// triangle drawing
			glColor3f(nz1*4.f, 0.f, 0.15f); glVertex3f(start + (z + 1), u[n][z + 1], start + (n));
			glColor3f(nz*4.f, 0.f, 0.15f);  glVertex3f(start + (z), u[n][z], start + (n));
			glColor3f(n1z*4.f, 0.f, 0.15f); glVertex3f(start + (z), u[n + 1][z], start + (n + 1));

			glColor3f(nz1*4.f, 0.f, 0.15f);  glVertex3f(start + (z + 1), u[n][z + 1], start + (n));
			glColor3f(n1z*4.f, 0.f, 0.15f); glVertex3f(start + (z), u[n + 1][z], start + (n + 1));
			glColor3f(n1z1*4.f, 0.f, 0.15f); glVertex3f(start + (z + 1), u[n + 1][z + 1], start + (n + 1));
		}
	}
	glEnd();
	
	if (wireframe)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glBegin(GL_LINES);
		glColor4f(1.f, 1.f, 1.f, 0.33f);

		for (int n = 1; n < N + 1; n++)
		{
			glVertex3f(start, u[0][n - 1], start + n - 1);
			glVertex3f(start, u[0][n], start + n);

			for (int z = 1; z < N + 1; z++)
			{
				glVertex3f(start + z - 1, u[n][z - 1], start + n);
				glVertex3f(start + z, u[n][z], start + n);

				glVertex3f(start + z, u[n][z], start + n);
				glVertex3f(start + z, u[n - 1][z], start + n - 1);
			}
		}

		glEnd();
	}

	glutSwapBuffers();

	glDisable(GL_BLEND);
}

void special(int key, int x, int y) {
	switch (key)
	{
	case GLUT_KEY_LEFT:
		camera_rot += 0.02f;

		camera_x -= cosf(camera_rot);
		camera_z -= sinf(camera_rot);

		break;

	case GLUT_KEY_RIGHT:
		camera_rot -= 0.02f;

		camera_x += cosf(camera_rot);
		camera_z += sinf(camera_rot);

		break;

	case GLUT_KEY_F1:
		wireframe = !wireframe;
		break;
	}

	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) 
	{
	case ' ':
		enabled = !enabled;
		break;
	case ',':
		scale += 0.1f;
		break;

	case '.':
		scale -= 0.1f;
		break;
	}
	glutPostRedisplay();
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(x, y);
	window = glutCreateWindow("3-D Window");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(special);
	glutReshapeFunc(reshape);
	GLenum glew_status = glewInit();

	// initial conditions
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{

			float in = (float)i / N;
			float jn = (float)j / N;
			u[i][j] = ((in > 0.25f && in <= .5f && jn > 0.25f && jn <= .5f) ? INITIAL_HEIGHT : 1);
			v[i][j] = ((in > 0.25f && in <= .5f && jn > 0.25f && jn <= .5f) ? INITIAL_HEIGHT : 1);
		}
	}

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glutMainLoop();
	return 0;
}