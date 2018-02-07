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

#define N 40
#define INITIAL_HEIGHT 40.f

GLfloat u[N][N];
GLfloat v[N][N];
GLfloat p[N][N];
GLfloat b[N][N];
GLfloat scale = 1.f;

GLfloat camera_x = 0.f;
GLfloat camera_y = 52.f;
GLfloat camera_z = 1.f;
GLfloat camera_rot = 0.f;

float dx = 2.f / (N - 1);
float dy = dx;
float vis = 0.1f;
float rho = 1.f;
float dt = 0.001f;
bool enabled = false;
bool wireframe = false;

void idle()
{
	glutPostRedisplay();
}

void reshape(int x, int y)
{
	glViewport(0, 0, x, y);
}

// build the b parameter of the poisson equation
void buildB()
{
	for (int i = 1; i < N - 1; i++)
		for (int j = 1; j < N - 1; j++)
		{
		b[i][j] = rho * (1 / dt*((u[i + 1][j] - u[i - 1][j]) / (2 * dx) + (v[i][j + 1] - v[i][j - 1]) / (2 * dy)) - ((u[i + 1][j] - u[i - 1][j]) / (2 * dx)) * ((u[i + 1][j] - u[i - 1][j]) / (2 * dx))
			- 2 * ((u[i][j + 1] - u[i][j - 1]) / (2 * dy)) * ((u[i][j + 1] - u[i][j - 1]) / (2 * dy)) - ((v[i][j + 1] - v[i][j - 1]) / (2 * dy)) * ((v[i][j + 1] - v[i][j - 1]) / (2 * dy)));
		}
}

void poisson()
{
	GLfloat pCopy[N][N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
		pCopy[i][j] = p[i][j];
		}
	float dys = dy*dy;
	float dxs = dx*dx;

	for (int i = 1; i < N - 1; i++)
	{
		for (int j = 1; j < N - 1; j++)
		{
			p[i][j] = (((pCopy[i + 1][j] + pCopy[i - 1][j]) * dys + (pCopy[i][j + 1] + pCopy[i][j - 1]) * dxs) / 2 * (dxs + dys)) - (dxs*dys) / (2 * (dxs + dys)) * b[i][j];
		}
	}

	for (int i = 0; i < N; i++)
	{
		p[N - 1][i] = p[N - 2][i];
		p[0][i] = p[1][i];
		p[i][0] = p[i][1];
		p[i][N - 1] = 0;
	}
}

void display()
{
	if (enabled)
	{
		GLfloat uCopy[N][N];
		GLfloat vCopy[N][N];
		for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			uCopy[i][j] = u[i][j];
			vCopy[i][j] = v[i][j];
		}

		buildB();
		poisson();

		float dys = dy*dy;
		float dxs = dx*dx;

		for (int i = 0; i < N - 1; i++)
		{
			for (int j = 0; j < N - 1; j++)
			{
				u[i][j] += 
				   -uCopy[i][j] * dt / dx*(uCopy[i][j] - uCopy[i - 1][j]) -
					vCopy[i][j] * dt / dy*(uCopy[i][j] - uCopy[i][j - 1]) -
					dt / (2 * rho*dx)*(p[i+1][j] - p[i-1][j]) +
					vis*(dt / dxs * (uCopy[i + 1][j] - 2 * uCopy[i][j] + uCopy[i - 1][j]) +
					(dt / dys * (uCopy[i][j + 1] - 2 * uCopy[i][j] + uCopy[i][j - 1])));

				v[i][j] += 
				   -uCopy[i][j] * dt / dx*(vCopy[i][j] - vCopy[i - 1][j]) -
					vCopy[i][j] * dt / dy*(vCopy[i][j] - vCopy[i][j - 1]) -
					dt / (2 * rho*dy)*(p[i][j + 1] - p[i][j - 1]) +
					vis*(dt / dxs * (vCopy[i + 1][j] - 2 * vCopy[i][j] + vCopy[i - 1][j]) +
					(dt / dys * (vCopy[i][j + 1] - 2 * vCopy[i][j] + vCopy[i][j - 1])));
			}
		}

		for (int i = 0; i < N; i++)
		{
			u[0][i] = u[i][0] = 0;
			v[0][i] = 1.f;
			u[0][i] = v[N - 1][i] = v[i][0] = v[i][N - 1] = u[N - 1][0] = 0;
		}
	}

	glClearColor(0.f, 0.f, 0.55f, 0.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(75.f, x / (float)y, 0.1f, 1000.f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(camera_x, camera_y, camera_z,
		0, 0, 0, 0, 1, 0);

	glScalef(scale, scale, scale);

	// draw surface
	glBegin(GL_TRIANGLES);
	float start = -(float)N / 2;
	for (int n = 0; n < N; n++)
	{
		for (int z = 0; z < N; z++)
		{
			// colours
			float nz = b[n][z] / 120.f;
			float nz1 = b[n][z + 1] / 120.f;
			float n1z = b[n + 1][z] / 120.f;
			float n1z1 = b[n + 1][z + 1] / 120.f;

			// triangle drawing
			if (nz1 < 0.f)
				glColor3f(0.f, (-nz1)*4.f, 0.15f);
			else
			glColor3f(nz1*4.f, 0.f, 0.15f); glVertex3f(start + (z + 1), 0.f, start + (n));
			if (nz < 0.f)
				glColor3f(0.f, (-nz)*4.f, 0.15f);
			else
				glColor3f(nz*4.f, 0.f, 0.15f); glVertex3f(start + (z), 0.f, start + (n));
			if (n1z < 0.f)
				glColor3f(0.f, (-n1z)*4.f, 0.15f);
			else
				glColor3f(n1z*4.f, 0.f, 0.15f); glVertex3f(start + (z), 0.f, start + (n + 1));

			if (nz1 < 0.f)
				glColor3f(0.f, (-nz1)*4.f, 0.15f);
			else
				glColor3f(nz1*4.f, 0.f, 0.15f);   glVertex3f(start + (z + 1), 0.f, start + (n));
			if (n1z < 0.f)
				glColor3f(0.f, (-n1z)*4.f, 0.15f);
			else
				glColor3f(n1z*4.f, 0.f, 0.15f); glVertex3f(start + (z), 0.f, start + (n + 1));
			if (n1z1 < 0.f)
				glColor3f(0.f, (-n1z1)*4.f, 0.15f);
			else
				glColor3f(n1z1*4.f, 0.f, 0.15f); glVertex3f(start + (z + 1), 0.f, start + (n + 1));
		}
	}
	glEnd();

	glBegin(GL_LINES);
		glColor4f(1.f, 1.f, 1.f, 0.33f);

		for (int n = 1; n < N + 1; n++)
		{
			for (int z = 1; z < N + 1; z++)
			{
				glVertex3f(start + z, 0.01f, start + n);
				float uf = u[n][z] * 2;
				float vf = v[n][z] * 2;
				glVertex3f(start + z + vf, 0.01f, start + n + uf);
			}
		}

		glEnd();

	glutSwapBuffers();

	glDisable(GL_BLEND);
}

void special(int key, int x, int y) {
	switch (key)
	{
	case GLUT_KEY_LEFT:
		//camera_rot += 0.2f;

		//camera_x -= cosf(camera_rot) / 10.f;
		//camera_z -= sinf(camera_rot) / 10.f;

		break;

	case GLUT_KEY_RIGHT:
		//camera_rot -= 0.2f;

		//camera_x += cosf(camera_rot) / 10.f;
		//camera_z += sinf(camera_rot) / 10.f;

		break;

	case GLUT_KEY_UP:
		camera_y += 0.5f;
		break;
	case GLUT_KEY_DOWN:
		camera_y -= 0.5f;
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
			u[i][j] = 0;
			v[i][j] = 0;
			p[i][j] = 0;
			b[i][j] = 0;
		}
	}

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glutMainLoop();
	return 0;
}