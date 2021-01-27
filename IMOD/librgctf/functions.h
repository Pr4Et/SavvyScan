// Stuff from functions.h
inline bool IsEven(int number_to_check)
{
	  if ( number_to_check % 2== 0 ) return true;
	  else return false;
}
std::vector<size_t> rankSort(const std::vector<float>& v_temp);
inline float deg_2_rad(float degrees)
{
  return degrees * PI / 180.;
}

inline int myroundint(double a)
{
	if (a > 0) return int(a + 0.5); else return int(a - 0.5);
}

inline int myroundint(float a)
{
	if (a > 0) return int(a + 0.5);	else return int(a - 0.5);
}

/* My replacement for wxPrint */
typedef void (*CharArgType)(const char *);

void internalSetPrintFunc(CharArgType func);
void wxPrintf(const char *format, ...);

typedef int (*WriteSliceType)(const char *, float *, int , int);

void internalSetWriteSliceFunc(WriteSliceType func);

int ctfNumOMPthreads(int optimalThreads);
int ctfOMPthreadNum();
double ctfWallTime(void);


