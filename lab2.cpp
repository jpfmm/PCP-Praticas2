#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <omp.h>
#include <string>
#include <vector>

#define SIZE 1536	// Only powers of 2 to simplify the code
#define MAX_THREADS 48 // maximum number of threads to run
#define REPETITIONS 1 // Number of measurements repetitions
#define TIME_RESOLUTION 1000000	// time measuring resolution (us)


/*
 *	Matrix-matrix arithmetic
 *
 */

using namespace std;

float m1[SIZE][SIZE], m2[SIZE][SIZE], result[SIZE][SIZE];
long long unsigned initial_time;
vector<long long unsigned> sequential_measurements;
vector<long long unsigned> measurements[MAX_THREADS/4 + 1];	// position 0 is for the sequential version
timeval t;

void clearCache (void) {
	double clearcache [30000000];

	for (unsigned i = 0; i < 30000000; ++i)
		clearcache[i] = i;
}

long long unsigned minElement (int thread) {
	long long unsigned best = 0;

	if (thread == - 1) {
		for (auto i : sequential_measurements) {
		if (i < best || !best)
			best = i;
		}
	} else {
		for (auto i : measurements[thread]) {
			if (i < best || !best)
				best = i;
		}
	}

	return best;
}

void writeResults (void) {
	ofstream file ("timings.dat");

	// write the optimised version
	file << 0 << " " << minElement(-1) << endl;

	// write the optimised version
	for (unsigned i = 0; (i * 4) <= MAX_THREADS; ++i)
		if ((i * 4) == 0)
			file << 1 << " " << minElement(i) << endl;
		else
			file << i * 4 << " " << minElement(i) << endl;

	file.close();
}

void start (void) {
	gettimeofday(&t, NULL);
	initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}

long long unsigned stop (int thread) {
	gettimeofday(&t, NULL);
	long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;

	if (thread == -1)
		sequential_measurements.push_back(final_time - initial_time);
	else
		measurements[thread].push_back(final_time - initial_time);

	return final_time - initial_time;
}

void fillMatrices (void) {

	for (unsigned i = 0; i < SIZE; ++i) {
		for (unsigned j = 0; j < SIZE; ++j) {
			m1[i][j] = ((float) rand()) / ((float) RAND_MAX);
			m2[i][j] = ((float) rand()) / ((float) RAND_MAX);
		}
	}
}

__attribute__((noinline))
void compute (float *x, float *y, float *r){
	*r += (*x) * (*y);
}

void matrixMult (void) {
	for (unsigned i = 0; i < SIZE; ++i) {
		for (unsigned j = 0; j < SIZE; ++j) {
			float x = 0;
			for (unsigned k = 0; k < SIZE; ++k) {
				compute(&m1[i][k], &m2[k][j], &x);
			}
			result[i][j] = x;
		}
	}
}

// assumes a square matrix witha size multiple of the block size
void matrixMultLab23 (void) {
	float x;
	unsigned s = 16;	// block size s * s

	#pragma omp parallel for
	for(unsigned jj = 0; jj < SIZE; jj += s) {
		for(unsigned kk = 0; kk < SIZE; kk += s) {
			for(unsigned i = 0; i < SIZE; ++i) {
				for(unsigned j = jj; j < jj + s; ++j) {
					x = 0;
					for(unsigned k = kk; k < kk + s; ++k) {
						x += m1[i][k] * m2[k][j];
					}
					#pragma omp critical
					result[i][j] += x;
				}
			}
		}
	}
}

int main (int argc, char *argv[]) {
	long long unsigned time;

	fillMatrices();

	// run the original code
	for (unsigned i = 0; i < REPETITIONS; ++i) {
		clearCache();

		start();
		matrixMult();
		stop(-1);
	}

	// optimised run! (SET MAX_THREADS to 1 for the HW 2.1 - it's not parallel yet)
	for (unsigned i = 0; (i * 4) <= MAX_THREADS; ++i) {

		omp_set_num_threads(i * 4);

		for (unsigned j = 0; j < REPETITIONS; ++j) {
			clearCache();

			start();
			// insert here your new function call
			matrixMultLab23();
			stop(i);
		}
	}

	writeResults();

	return 1;
}
