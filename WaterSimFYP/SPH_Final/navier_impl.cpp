#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "../GL/glew.h"
#include "../GL/freeglut.h"

#define GLM_FORCE_RADIANS
#include "../include/glm/glm.hpp"
#include "../include/glm/gtx/norm.hpp"
#include "../include/glm/gtc/matrix_transform.hpp"
#include "../include/glm/gtc/type_ptr.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>
#include "navier_impl_helper.h"
#include "Grid.h"
#include "Particle.h"
#include "HighResTimer.h"

// SIMD header (SSE)
#include <xmmintrin.h>

GLuint window;
int x = 800, y = 900;

#define N 10000
#define N_STACK 3333 // N / 3; saves a division
#define T 1 / 60.f

std::vector<Particle> particles;
Grid* grid;
GLfloat scale = 1.f;

// gather stats after specified amount of frames
// comment in/out to remove/add this feature
//
//#define ACCUMALATE_DATA

#ifdef ACCUMALATE_DATA

#define FRAMES_TEST 1800
float stats[FRAMES_TEST];
float density_stats[FRAMES_TEST];
float pressure_build_stats[FRAMES_TEST];
float pressure_stats[FRAMES_TEST];
float viscosity_stats[FRAMES_TEST];
int frame = 0;

#endif

// viscosity
float vis = .12f;
// radius of support
float rad = 5.5f;
float h = rad * 1.15f; // radius of particles
float mass = 1.f;
// constants for ideal gas equation
float k = rad / 1000.f;
float rest_density = 2.f;

float h2 = h*h;
bool enabled = false;

bool mouse_enabled = false;
int mouse_x = 0;
int mouse_y = 0;

void reshape(int x_, int y_)
{
	glViewport(0, 0, x_, y_);
	x = x_;
	y = y_;
}

float density_time = 0.f;
void density()
{
	std::chrono::time_point<HighResClock> pre_update, post_update;
	pre_update = HighResClock::now();

#pragma omp parallel for
	for (int i = 0; i < (int)particles.size(); i++)
	// regarding the int conversion, OpenMP does not allow for
	// for loops to initialise with an unsigned int
	{
		Particle* parts[N_STACK]; // a particle could have a very high number of neighbours
		
		// calculate this particle's cell position
		int x = MAX(particles[i].Cell->X - 1, 0);
		int y = MAX(particles[i].Cell->Y - 1, 0);
		// are we on a corner?
		int xx = x + (particles[i].Cell->XCorner) ? 2 : 3;
		int yy = x + (particles[i].Cell->YCorner) ? 2 : 3;

		int size = 0;
		// populate potential neighbour list
		for (int a = x; a < x + xx; a++)
			for (int b = y; b < y + yy; b++)
				for (unsigned int p = 0; p < grid->getCell(a, b)->Particles.size(); p++)
					parts[size++] = grid->getCell(a, b)->Particles[p];

		for (int j = 0; j < size; j++)
		{
			if (parts[j] == &particles[i])
				continue;

			Vector diff = parts[j]->Position - particles[i].Position;
			float diff_len = diff.x * diff.x + diff.y * diff.y;
			// is this particle within the range of support?
			if (diff_len < h2)
			{
				__m128 m1, m2, m3, m4, m5, m6, mout;

				// q = r / h
				// rho += (1 - q)^2
				m1 = _mm_set_ps1(diff_len);
				m2 = _mm_sqrt_ps(m1);
				m3 = _mm_set_ps1(h);
				m4 = _mm_set_ps1(1);
				m5 = _mm_div_ps(m2, m3);
				m6 = _mm_sub_ps(m4, m5);
				mout = _mm_mul_ps(m6, m6);

				// our spike weighting kernel
				float q2;

				_mm_storeu_ps(&q2, mout);

				particles[i].Rho += q2;

				particles[i].Neighbours.push_back(parts[j]);
			}
		}
	}

	post_update = HighResClock::now();

	std::chrono::microseconds m = std::chrono::duration_cast<std::chrono::microseconds>(post_update - pre_update);
	std::chrono::duration<float> time = m;

	density_time = time.count();
}

float pressure_build = 0.f;
void pressure()
{
	std::chrono::time_point<HighResClock> pre_update, post_update;
	pre_update = HighResClock::now();
	// use an equation of state to replace poisson
	// poisson takes too long to calculate
	// we will use an ideal gas equation
	// P = k * rho
	__m128 mk, rk;
	mk = _mm_set_ps1(k);
	rk = _mm_set_ps1(rest_density);
#pragma omp parallel for
	for (int i = 0; i < (int)particles.size() - 3; i += 4)
	{
		float rhos[4];
		rhos[0] = particles[i].Rho;
		rhos[1] = particles[i+1].Rho;
		rhos[2] = particles[i+2].Rho;
		rhos[3] = particles[i+3].Rho;

		__m128 mrhos = _mm_loadu_ps(rhos);

		__m128 msub, mout;
		msub = _mm_sub_ps(mrhos, rk);
		mout = _mm_mul_ps(mk, msub);

		float press[4];
		_mm_storeu_ps(press, mout);

		for (int j = 0; j < 4; j++)
			particles[i+j].Pressure = press[j];
		// we should also add an additional artifical pressure term
		// here to solve particle clustering
	}

	// calculate remaining particle's pressure
#pragma omp parallel for
	for (int i = N - (N % 4); i < (int)particles.size(); i++)
	{
		particles[i].Pressure = k * (particles[i].Rho - rest_density);
	}

	post_update = HighResClock::now();

	std::chrono::microseconds m = std::chrono::duration_cast<std::chrono::microseconds>(post_update - pre_update);
	std::chrono::duration<float> time = m;

	pressure_build = time.count();
}

float pressure_time = 0.f;
void pressure_force()
{
	std::chrono::time_point<HighResClock> pre_update, post_update;
	pre_update = HighResClock::now();
#pragma omp parallel for
	for (int i = 0; i < (int)particles.size(); i++)
	{		
		// predefine SIMD variables
		__m128 one = _mm_set_ps1(1);
		__m128 mh = _mm_set_ps1(h);
		__m128 mP = _mm_set_ps1(particles[i].Pressure);

		// calculate pressure force for each neighbour
		// however, we can group them into 4's with SIMD
		for (unsigned int j = 0; j < particles[i].Neighbours.size(); j+=4)
		{
			if (j + 4 > particles[i].Neighbours.size())
				break;

			Vector diffs[4];
			diffs[0] = particles[i].Neighbours[j]->Position - particles[i].Position;
			diffs[1] = particles[i].Neighbours[j+1]->Position - particles[i].Position;
			diffs[2] = particles[i].Neighbours[j+2]->Position - particles[i].Position;
			diffs[3] = particles[i].Neighbours[j+3]->Position - particles[i].Position;

			float l1[4], l2[4];
			l1[0] = diffs[0].x;
			l1[1] = diffs[0].y;
			l1[2] = diffs[1].x;
			l1[3] = diffs[1].y;
			l2[0] = diffs[2].x;
			l2[1] = diffs[2].y;
			l2[2] = diffs[3].x;
			l2[3] = diffs[3].y;

			__m128 len1 = _mm_loadu_ps(l1);
			__m128 len2 = _mm_loadu_ps(l2);

			float lengths[4];
			// get the final length for each particle, 
			//end parameter determines the input and outputs for the function
			lengths[0] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(len1, len1, 0x31)));
			lengths[1] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(len1, len1, 0xC1)));
			lengths[2] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(len2, len2, 0x31)));
			lengths[3] = _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(len2, len2, 0xC1)));

			float press[4];
			press[0] = particles[i].Neighbours[j]->Pressure;
			press[1] = particles[i].Neighbours[j+1]->Pressure;
			press[2] = particles[i].Neighbours[j+2]->Pressure;
			press[3] = particles[i].Neighbours[j+3]->Pressure;

			__m128 mLen = _mm_loadu_ps(lengths);
			__m128 mPress = _mm_loadu_ps(press);

			__m128 m3, m33, m6, mout;
			m3 = _mm_div_ps(mLen, mh);
			m33 = _mm_sub_ps(one, m3);
			m6 = _mm_add_ps(mP, mPress);
			mout = _mm_mul_ps(m6, m33);

			// dm = (1 - q) * pressurei + pressurej
			float dm[4];
			_mm_storeu_ps(dm, mout);

			for (int q = 0; q < 4; q++)
			{
				Vector D = (diffs[q] / lengths[q]) * dm[q];
				particles[i].Neighbours[j+q]->Force += D;
				particles[i].Force -= D;
			}
		}
		// calculate remaining particles
		if (particles[i].Neighbours.size() % 4 > 0)
		for (unsigned int j = particles[i].Neighbours.size() - (particles[i].Neighbours.size() % 4); j < particles[i].Neighbours.size(); j++)
		{
			Vector diff = particles[i].Neighbours[j]->Position - particles[i].Position;
			float diff_len = glm::length(diff);

			__m128 m1, m3, m33, m5, m6, mout;
			m1 = _mm_set_ps1(diff_len);
			m3 = _mm_div_ps(m1, mh);
			m33 = _mm_sub_ps(one, m3);
			m5 = _mm_set_ps1(particles[i].Neighbours[j]->Pressure);
			m6 = _mm_add_ps(mP, m5);
			mout = _mm_mul_ps(m6, m33);

			float dm;
			_mm_storeu_ps(&dm, mout);

			Vector D = (diff / diff_len) * dm;
			particles[i].Neighbours[j]->Force += D;
			particles[i].Force -= D;
		}
	}

	post_update = HighResClock::now();

	std::chrono::microseconds m = std::chrono::duration_cast<std::chrono::microseconds>(post_update - pre_update);
	std::chrono::duration<float> time = m;

	pressure_time = time.count();
}

float viscosity_time = 0.f;
void viscosity(float dt)
{
	std::chrono::time_point<HighResClock> pre_update, post_update;
	pre_update = HighResClock::now();
#pragma omp parallel for
	for (int i = 0; i < (int)particles.size(); i++)
		for (unsigned int j = 0; j < particles[i].Neighbours.size(); j++)
		{
			Vector diff = particles[i].Neighbours[j]->Position - particles[i].Position;
			float diff_len = glm::length(diff);

			float q = diff_len / h;
			Vector diff_norm = (diff / diff_len);
			Vector vel_diff = particles[i].Velocity - particles[i].Neighbours[j]->Velocity;
			float u = vel_diff.x * diff_norm.x + vel_diff.y * diff_norm.y;
			if (u > 0)
			{
				__m128 m1, m2, m3, m4, m5, m6, m7, mout;
				float q1 = 1 - q;
				m2 = _mm_set_ps1(vis); // vis should be 2 terms, linear and quadratic
				m3 = _mm_set_ps1(u);
				m4 = _mm_mul_ps(m2, m3);
				m5 = _mm_mul_ps(m3, m3);
				m6 = _mm_mul_ps(m5, m2);
				m7 = _mm_add_ps(m4, m6);
				m1 = _mm_set_ps1(q1);
				mout = _mm_mul_ps(m1, m7);

				float qout;
				_mm_storeu_ps(&qout, mout);
		
				// impulse calculation
				// I = (1 - q) * (u * viscosity_linear + u * u * viscosity_quadratic) * diff_norm / 2
				Vector I = qout * diff_norm * 0.5f;

				particles[i].Velocity -= I;
				particles[i].Neighbours[j]->Velocity += I;
			}
		}
	post_update = HighResClock::now();

	std::chrono::microseconds m = std::chrono::duration_cast<std::chrono::microseconds>(post_update - pre_update);
	std::chrono::duration<float> time = m;

	viscosity_time = time.count();
}

void addParticlesFromMouse()
{
#define PARTICLES_ADD 10

	for (int i = 0; i < PARTICLES_ADD; i++)
	{
		Particle p;
		p.Position = Vector(mouse_x, mouse_y);
		p.Velocity = randomVector(1.f);

		particles.push_back(p);
	}
}

void update_particles(float dt)
{
	grid->clear();

	// integrate position and velocity
	// here
	// ---------------------
#pragma omp parallel for
	for (int i = 0; i < (int)particles.size(); i++)
	{
		Vector oldPos = particles[i].Position;
		Vector oldVel = particles[i].Velocity;

		// dt = 1, constant
		particles[i].Velocity += particles[i].Force;
		particles[i].Position += (oldVel + particles[i].Velocity) * 0.5f;
		
		particles[i].Force = Vector(0, -0.009f); // apply constant gravity

		// basic spring force if particle goes out of bounds
		if (particles[i].Position.x < 0)
			particles[i].Force.x -= (particles[i].Position.x - 0) / 8;
		else if (particles[i].Position.x > x)
			particles[i].Force.x -= (particles[i].Position.x - x) / 8;
		if (particles[i].Position.y < h)
			particles[i].Force.y -= (particles[i].Position.y - h) / 8;
		else if (particles[i].Position.y > y)
			particles[i].Force.y -= (particles[i].Position.y - y) / 8;

		// bit of a hack here, to prevent the simulation from going unstable
		// only happes if particles start out of bounds (spring issue)
		float max_vel = 7.5f;
		float vel_mag  = glm::length2(particles[i].Velocity);
		// If the velocity is greater than the max velocity, then cap it.
		if (vel_mag > max_vel)
			particles[i].Velocity = particles[i].Velocity * .5f;

		particles[i].Rho = 0.f;
		particles[i].Neighbours.clear();
	}
	
	grid->populate(particles.data(), particles.size());
	
	density();

	pressure();

	pressure_force();

	viscosity(dt);

	// if mouse is down add particles
	if (mouse_enabled)
		addParticlesFromMouse();
}

void idle()
{
	if (enabled)
	{
		std::chrono::time_point<HighResClock> pre_update, post_update;
		pre_update = HighResClock::now();
		update_particles(T);
		post_update = HighResClock::now();

		std::chrono::microseconds m = std::chrono::duration_cast<std::chrono::microseconds>(post_update - pre_update);
		std::chrono::duration<float> time = m;

		float milliseconds = time.count() * 1000.f;

#ifdef ACCUMALATE_DATA
		stats[frame] = milliseconds;
		density_stats[frame] = density_time * 1000.f;
		pressure_build_stats[frame] = pressure_build * 1000.f;
		pressure_stats[frame] = pressure_time * 1000.f;
		viscosity_stats[frame] = viscosity_time * 1000.f;
		frame++;
#endif

		printf("Update(%i): %4.3f / D: %4.3f / PB: %4.3f / P: %4.3f / V: %3.4f\n", particles.size(), milliseconds, density_time * 1000.f, pressure_build * 1000.f, pressure_time * 1000.f, viscosity_time * 1000.f);
	}

#ifdef ACCUMALATE_DATA
	if (frame == FRAMES_TEST)
		exit(0); // finish test
#endif

	glutPostRedisplay();
}

void display()
{
	glClearColor(0.6f, 1.f, 1.f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, x, 0, y, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glScalef(scale, scale, scale);

	// draw particles
	glPointSize(h * 1.25f);
	glBegin(GL_POINTS);
	glColor3f(0.06f, 0.21f, 0.909f);

	for (unsigned int i = 0; i < particles.size(); i++)
	{
		glVertex2f(particles[i].Position.x, particles[i].Position.y);
	}

	glEnd();

	glutSwapBuffers();
}

void special(int key, int x, int y) 
{

	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
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
	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}

void finished()
{
#ifdef ACCUMALATE_DATA
	// if we finish, lets process the data
	FILE* f;
	fopen_s(&f, "./stats.txt", "w");

	for (int i = 0; i < FRAMES_TEST; i++)
	{
		float combine = density_stats[i] + pressure_build_stats[i] + pressure_stats[i] + viscosity_stats[i];
		fprintf(f, "%f	%f	%f	%f	%f\n", stats[i],density_stats[i],pressure_build_stats[i],pressure_stats[i],viscosity_stats[i]);
	}

	fclose(f);

#endif
}

void mouse(int button, int state, int _x, int _y)
{
	if (button == GLUT_LEFT_BUTTON)
	{
		mouse_x = _x;
		mouse_y = y - _y;

		mouse_enabled = (state == GLUT_DOWN) ? true : false;
	}
}

void mouse_move(int _x, int _y)
{
	mouse_x = _x;
	mouse_y = y - _y;
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
	glutInitWindowSize(x, y);
	window = glutCreateWindow("Water Simulation Demo");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(special);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutMotionFunc(mouse_move);

	int row_amount = (int)((x - (50 * 2)) / (h * 1.25f));
	int row_num = (int)ceil(N / (float)row_amount);

	grid = new Grid(x, y, x / 11, y / 11);
	Vector centre(x / 2, 0);

	// set initial position
	int n_amount = 0;
	particles.resize(N);
	for (int i = 0; i < row_num; i++)
	{
		for (int j = 0; j < row_amount; j++)
		{
			particles[(i * row_amount) + j].Position = Vector(100 + j * (h * 1.1f), 50 + (i * (h * 1.1f)));
			// slight attraction to the centre removes particles stacking
			particles[(i * row_amount) + j].Velocity += 0.001f * glm::normalize(centre - particles[(i * row_amount) + j].Position);
			if (++n_amount >= N)
				break;
		}
		if (n_amount >= N)
			break;
	}

	grid->populate(particles.data(), N);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	// we want to exit!
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

	atexit(finished);

	glutMainLoop();
	return 0;
}