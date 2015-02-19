
#pragma once
#include <windows.h>
#include <string>
#include <iostream>
#include <fstream>
#include <strstream>
#include <stdarg.h>
#include <process.h>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <vector>

using namespace std;

#define uint16 unsigned int
#define uint32 unsigned int // as long as Visual C compiler int = long
#define int16 int
#define int32 long

#define M_PI       3.14159265358979323846

enum el_component 
{
	COMP_UNDEF,

	E_X,
	E_Y,
	E_Z,
	E_XY,
	E_XZ,
	E_YZ,
	E_VOL,
	E_1,
	E_2,
	E_3,

	S_X,
	S_Y,
	S_Z,
	S_XY,
	S_XZ,
	S_YZ,
	S_P,
	S_1,
	S_2,
	S_3,

	U_X,
	U_Y,
	U_Z,

	POS_X,
	POS_Y,
	POS_Z,

	COMP_LAST
};

const char* const el_component_labels[]={"UNDEFINED","E_X","E_Y","E_Z","E_XY","E_XZ","E_YZ","E_VOL","E_1","E_2","E_3","S_X","S_Y","S_Z","S_XY","S_XZ","S_YZ","S_P","S_1","S_2","S_3","U_X","U_Y","U_Z","POS_X","POS_Y","POS_Z","COMP_LAST"};

// use it in material matrix creator
#define ANALYSIS_3D 1
#define ANALYSIS_2D_PLANE_STRESS 2
#define ANALYSIS_2D_PLANE_STRAIN 3
#define ANALYSIS_2D_AXISYMMETRIC 4

#define GP_MEAN 100 //среднее значение по элементу

//ключ закрепления != 0  !
#define D_UX 1
#define D_UY 2
#define D_UZ 3

#define F_X 4
#define F_Y 5
#define F_Z 6

const char SYS_VERSION[] = "1.1";
const char SYS_DATA[] = "30.01.12";
const char log_file_name[]="log.txt";
const bool debug_mode = true;

class FE_Storage_Interface;

class Log_opts 
{
public:
	Log_opts() 
	{
		output_lock=CreateMutex(NULL, FALSE, NULL);
		ofstream file(log_file_name,ios::trunc);
		file.close();
		
	};
	HANDLE output_lock;
};

void warning(const char* logline, ...);
void debug(const char* logline, ...);
void debug(string &str);
void log(string &str);
void log(const char* logline, ...);
void error(const char* logline, ...);
void echo (const char* logline, ...);
void echolog(const char* logline, ...);
uint32 tick();

string IntToStr (uint32 dig);
int32 npow(int16 dig, uint16 power);

vector<string> read_tokens(char *input);
void del_spaces (string &str);

bool read_ans_data (const char *filename, FE_Storage_Interface *storage);

class Timer
{
public:
	Timer(bool _start = false) : start_time(0), end_time(0)
	{
		if (_start)
			start();
	}
	void start()
	{
		start_time = clock();
		end_time = start_time;
	}
	double stop()
	{
		end_time = clock();
		return time();
	}
	double time()
	{
		return ((double)end_time - start_time) / CLOCKS_PER_SEC;
	}
private:
	clock_t start_time;
	clock_t end_time;
};

uint16 str2dof (string dof_name);
