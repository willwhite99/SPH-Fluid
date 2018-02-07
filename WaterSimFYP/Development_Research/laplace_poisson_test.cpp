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

GLfloat p[N][N];
GLfloat b[N][N];
GLfloat scale = 1.f;

GLfloat camera_x = 0.f;
GLfloat camera_y = 25.f;
GLfloat camera_z = 50.f;
GLfloat camera_rot = 0.f;

int c = 1;
float xmin = 0.f;
float xmax = 50.f;
float dx = (xmax-xmin) / (N - 1);
float dy = dx;
float vis = 0.01f;
float dt = 0.0009f*dx*dy / vis;
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

void laplace(GLfloat p[N][N], GLfloat dx, GLfloat dy, GLfloat ltarget)
{
	GLfloat lnorm = 40;
	GLfloat pCopy[N][N];

	while (lnorm > ltarget)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
			pCopy[i][j] = p[i][j];
			}

		for (int i = 1; i < N - 1; i++)
			for (int j = 1; j < N - 1; j++)
			{
			p[i][j] = (dy*dy*(pCopy[i + 1][j] + pCopy[i - 1][j]) + dx*dx*(pCopy[i][j + 1] + pCopy[i][j - 1])) / (2 * (dx*dx + dy*dy));
			}

		p[0][0] = (dy*dy*(pCopy[1][0] + pCopy[N - 1][0]) + dx*dx*(pCopy[0][1] + pCopy[0][N - 1])) / (2 * (dx*dx + dy*dy));
		p[N - 1][N - 1] = (dy*dy*(pCopy[0][N - 1] + pCopy[N - 2][N - 1]) + dx*dx*(pCopy[N - 1][0] + pCopy[N - 1][N - 2])) / (2 * (dx*dx + dy*dy));

		for (int i = 0; i < N; i++)
		{
			p[i][0] = 0;
			float in = ((float)i / N);
			p[i][N - 1] = in * 40;
			p[0][i] = p[1][i];
			p[N - 1][i] = p[N - 1][i];
		}

		float absP = 0.f;
		float absPCopy = 0.f;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
			absP += abs(p[i][j]);
			absPCopy += abs(pCopy[i][j]);
			}

		lnorm = (absP - absPCopy) / absPCopy;
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

	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 0; j < N - 1; j++)
		{
			float dys = dy*dy;
			float dxs = dx*dx;
			p[i][j] = (dys / (2 * (dxs + dys))*(pCopy[i + 1][j] + pCopy[i - 1][j]) + dxs / (2 * (dxs + dys))*(pCopy[i][j + 1] + pCopy[i][j - 1]) - b[i][j] * dxs*dys / (2 * (dxs + dys)));

			p[0][j] = p[N - 1][j] = p[i][0] = p[i][N - 1] = 0;
		}
	}
}

void display()
{
	if (enabled)
	{
		poisson();
	}

	glClearColor(0.f, 0.f, 0.55f, 0.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(75.f, x / (float)y, 0.1f, 1000.f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(camera_x, camera_y, camera_z,
		0, 10 + camera_y - 25.f, 0, 0, 1, 0);

	glScalef(scale, scale, scale);

	// draw surface
	glBegin(GL_TRIANGLES);
	float start = -(float)N / 2;
	for (int n = 0; n < N; n++)
	{
		for (int z = 0; z < N; z++)
		{
			// colours
			float nz = p[n][z] / 255.f;
			float nz1 = p[n][z + 1] / 255.f;
			float n1z = p[n + 1][z] / 255.f;
			float n1z1 = p[n + 1][z + 1] / 255.f;

			// triangle drawing
			glColor3f(nz1*4.f, 0.f, 0.15f); glVertex3f(start + (z + 1), p[n][z + 1], start + (n));
			glColor3f(nz*4.f, 0.f, 0.15f);  glVertex3f(start + (z), p[n][z], start + (n));
			glColor3f(n1z*4.f, 0.f, 0.15f); glVertex3f(start + (z), p[n + 1][z], start + (n + 1));

			glColor3f(nz1*4.f, 0.f, 0.15f);  glVertex3f(start + (z + 1), p[n][z + 1], start + (n));
			glColor3f(n1z*4.f, 0.f, 0.15f); glVertex3f(start + (z), p[n + 1][z], start + (n + 1));
			glColor3f(n1z1*4.f, 0.f, 0.15f); glVertex3f(start + (z + 1), p[n + 1][z + 1], start + (n + 1));
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
			glVertex3f(start, p[0][n - 1], start + n - 1);
			glVertex3f(start, p[0][n], start + n);

			for (int z = 1; z < N + 1; z++)
			{
				glVertex3f(start + z - 1, p[n][z - 1], start + n);
				glVertex3f(start + z, p[n][z], start + n);

				glVertex3f(start + z, p[n][z], start + n);
				glVertex3f(start + z, p[n - 1][z], start + n - 1);
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
			p[i][j] = 0;
	}

	// laplace setup
	//for (int i = 0; i < N; i++)
	//{

	//	float in = ((float)i / N);
	//	p[i][N - 1] = in * 40;
	//}

	// poisson setup
	b[N / 4][N / 4] = 100.f;
	b[3 * N / 4][3 * N / 4] = -100.f;

	// apply laplace
	//laplace(p, dx, dy, 0.0001f);

	// apply poisson
	//for (int i = 0; i < 100; i++)
	//	poisson();
	//for (int i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < N; j++)
	//		p[i][j] = b[i][j];
	//}

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glutMainLoop();
	return 0;
}