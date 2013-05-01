
#ifndef PARALLEL_H
#define PARALLEL_H




#ifdef PARALLEL_THREADS_WIN32
#include <windows.h>
#else

#ifdef PARALLEL_THREADS_PTHREAD
#include <pthread.h>
#else
#define PARALLEL_THREADS_PASSTHROUGH
#endif

#endif


/*
some examples


-------------
// to call a member function in parallel

class someclass {
  void foo() {
    ParallelExecutor(idealNumThreads, makeClassFunctor(this, &someclass::bar));
  }
  void bar(int nt, int id) {
    // do some parallel work
  }
};

-------------
// to call a global function / functor

void pe2(int nt, int id, int bi, const char *s) {
  // do some parallel work
}

int bi=0;
char *s="";
ParallelExecutor(idealNumThreads, pe2, bi, s);



// To add more versions for more parameters to your functions, you have to define the appropriate macro's below,
// and instantiate them at the bottom.  For the class functors to work, you need 2 extra variables since the 
// id and numThreads are not treated as special.  That's why the #defines go up to 7, but the parallel functions
// can only take 5 variables.


*/



#define PE_TEMPLATE_DEF0
#define PE_TEMPLATE_DEF1 <typename T1>
#define PE_TEMPLATE_DEF2 <typename T1, typename T2>
#define PE_TEMPLATE_DEF3 <typename T1, typename T2, typename T3>
#define PE_TEMPLATE_DEF4 <typename T1, typename T2, typename T3, typename T4>
#define PE_TEMPLATE_DEF5 <typename T1, typename T2, typename T3, typename T4, typename T5>
#define PE_TEMPLATE_DEF6 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
#define PE_TEMPLATE_DEF7 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
#define PE_TEMPLATE_DEF8 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
#define PE_TEMPLATE_DEF9 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
#define PE_TEMPLATE_DEF10 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
#define PE_TEMPLATE_DEF11 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
#define PE_TEMPLATE_DEF12 <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>

#define PE_TEMPLATE_ADD_DEF0
#define PE_TEMPLATE_ADD_DEF1 , typename T1
#define PE_TEMPLATE_ADD_DEF2 , typename T1, typename T2
#define PE_TEMPLATE_ADD_DEF3 , typename T1, typename T2, typename T3
#define PE_TEMPLATE_ADD_DEF4 , typename T1, typename T2, typename T3, typename T4
#define PE_TEMPLATE_ADD_DEF5 , typename T1, typename T2, typename T3, typename T4, typename T5
#define PE_TEMPLATE_ADD_DEF6 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6
#define PE_TEMPLATE_ADD_DEF7 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7
#define PE_TEMPLATE_ADD_DEF8 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8
#define PE_TEMPLATE_ADD_DEF9 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9
#define PE_TEMPLATE_ADD_DEF10 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10
#define PE_TEMPLATE_ADD_DEF11 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11
#define PE_TEMPLATE_ADD_DEF12 , typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12

#define PE_TEMPLATE_ADD_CALL0 
#define PE_TEMPLATE_ADD_CALL1 , T1
#define PE_TEMPLATE_ADD_CALL2 , T1, T2
#define PE_TEMPLATE_ADD_CALL3 , T1, T2, T3
#define PE_TEMPLATE_ADD_CALL4 , T1, T2, T3, T4
#define PE_TEMPLATE_ADD_CALL5 , T1, T2, T3, T4, T5
#define PE_TEMPLATE_ADD_CALL6 , T1, T2, T3, T4, T5, T6
#define PE_TEMPLATE_ADD_CALL7 , T1, T2, T3, T4, T5, T6, T7
#define PE_TEMPLATE_ADD_CALL8 , T1, T2, T3, T4, T5, T6, T7, T8
#define PE_TEMPLATE_ADD_CALL9 , T1, T2, T3, T4, T5, T6, T7, T8, T9
#define PE_TEMPLATE_ADD_CALL10 , T1, T2, T3, T4, T5, T6, T7, T8, T9, T10
#define PE_TEMPLATE_ADD_CALL11 , T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11
#define PE_TEMPLATE_ADD_CALL12 , T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12

#define PE_TEMPLATE_CALL0
#define PE_TEMPLATE_CALL1 <T1>
#define PE_TEMPLATE_CALL2 <T1,T2>
#define PE_TEMPLATE_CALL3 <T1,T2,T3>
#define PE_TEMPLATE_CALL4 <T1,T2,T3,T4>
#define PE_TEMPLATE_CALL5 <T1,T2,T3,T4,T5>
#define PE_TEMPLATE_CALL6 <T1,T2,T3,T4,T5,T6>
#define PE_TEMPLATE_CALL7 <T1,T2,T3,T4,T5,T6,T7>
#define PE_TEMPLATE_CALL8 <T1,T2,T3,T4,T5,T6,T7,T8>
#define PE_TEMPLATE_CALL9 <T1,T2,T3,T4,T5,T6,T7,T8,T9>
#define PE_TEMPLATE_CALL10 <T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>
#define PE_TEMPLATE_CALL11 <T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>
#define PE_TEMPLATE_CALL12 <T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12>

#define PE_PARAMETER_ADD_DEF0(T,p)
#define PE_PARAMETER_ADD_DEF1(T,p) ,T##1 p##1
#define PE_PARAMETER_ADD_DEF2(T,p) ,T##1 p##1, T##2 p##2
#define PE_PARAMETER_ADD_DEF3(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3
#define PE_PARAMETER_ADD_DEF4(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4
#define PE_PARAMETER_ADD_DEF5(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5
#define PE_PARAMETER_ADD_DEF6(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6
#define PE_PARAMETER_ADD_DEF7(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7
#define PE_PARAMETER_ADD_DEF8(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8
#define PE_PARAMETER_ADD_DEF9(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9
#define PE_PARAMETER_ADD_DEF10(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10
#define PE_PARAMETER_ADD_DEF11(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10, T##11 p##11
#define PE_PARAMETER_ADD_DEF12(T,p) ,T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10, T##11 p##11, T##12 p##12

#define PE_PARAMETER_CALL0(p)
#define PE_PARAMETER_CALL1(p) p##1
#define PE_PARAMETER_CALL2(p) p##1, p##2
#define PE_PARAMETER_CALL3(p) p##1, p##2, p##3
#define PE_PARAMETER_CALL4(p) p##1, p##2, p##3, p##4
#define PE_PARAMETER_CALL5(p) p##1, p##2, p##3, p##4, p##5
#define PE_PARAMETER_CALL6(p) p##1, p##2, p##3, p##4, p##5, p##6
#define PE_PARAMETER_CALL7(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7
#define PE_PARAMETER_CALL8(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8
#define PE_PARAMETER_CALL9(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9
#define PE_PARAMETER_CALL10(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10
#define PE_PARAMETER_CALL11(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10, p##11
#define PE_PARAMETER_CALL12(p) p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10, p##11, p##12

#define PE_PARAMETER_ADD_CALL0(p)
#define PE_PARAMETER_ADD_CALL1(p) ,p##1
#define PE_PARAMETER_ADD_CALL2(p) ,p##1, p##2
#define PE_PARAMETER_ADD_CALL3(p) ,p##1, p##2, p##3
#define PE_PARAMETER_ADD_CALL4(p) ,p##1, p##2, p##3, p##4
#define PE_PARAMETER_ADD_CALL5(p) ,p##1, p##2, p##3, p##4, p##5
#define PE_PARAMETER_ADD_CALL6(p) ,p##1, p##2, p##3, p##4, p##5, p##6
#define PE_PARAMETER_ADD_CALL7(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7
#define PE_PARAMETER_ADD_CALL8(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8
#define PE_PARAMETER_ADD_CALL9(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9
#define PE_PARAMETER_ADD_CALL10(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10
#define PE_PARAMETER_ADD_CALL11(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10, p##11
#define PE_PARAMETER_ADD_CALL12(p) ,p##1, p##2, p##3, p##4, p##5, p##6, p##7, p##8, p##9, p##10, p##11, p##12

#define PE_PARAMETER_DEF0(T,p)
#define PE_PARAMETER_DEF1(T,p) T##1 p##1
#define PE_PARAMETER_DEF2(T,p) T##1 p##1, T##2 p##2
#define PE_PARAMETER_DEF3(T,p) T##1 p##1, T##2 p##2, T##3 p##3
#define PE_PARAMETER_DEF4(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4
#define PE_PARAMETER_DEF5(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5
#define PE_PARAMETER_DEF6(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6
#define PE_PARAMETER_DEF7(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7
#define PE_PARAMETER_DEF8(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8
#define PE_PARAMETER_DEF9(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9
#define PE_PARAMETER_DEF10(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10
#define PE_PARAMETER_DEF11(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10, T##11 p##11
#define PE_PARAMETER_DEF12(T,p) T##1 p##1, T##2 p##2, T##3 p##3, T##4 p##4, T##5 p##5, T##6 p##6, T##7 p##7, T##8 p##8, T##9 p##9, T##10 p##10, T##11 p##11, T##12 p##12

#define PE_MEMBER_DEF0(T,p)
#define PE_MEMBER_DEF1(T,p) T##1 p##1;
#define PE_MEMBER_DEF2(T,p) T##1 p##1; T##2 p##2;
#define PE_MEMBER_DEF3(T,p) T##1 p##1; T##2 p##2; T##3 p##3;
#define PE_MEMBER_DEF4(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4;
#define PE_MEMBER_DEF5(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5;
#define PE_MEMBER_DEF6(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6;
#define PE_MEMBER_DEF7(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7;
#define PE_MEMBER_DEF8(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8;
#define PE_MEMBER_DEF9(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8; T##9 p##9;
#define PE_MEMBER_DEF10(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8; T##9 p##9; T##10 p##10;
#define PE_MEMBER_DEF11(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8; T##9 p##9; T##10 p##10; T##11 p##11;
#define PE_MEMBER_DEF12(T,p) T##1 p##1; T##2 p##2; T##3 p##3; T##4 p##4; T##5 p##5; T##6 p##6; T##7 p##7; T##8 p##8; T##9 p##9; T##10 p##10; T##11 p##11; T##12 p##12;

#define PE_INITIALIZE0(m,p)
#define PE_INITIALIZE1(m,p) : m##1(p##1)
#define PE_INITIALIZE2(m,p) : m##1(p##1), m##2(p##2)
#define PE_INITIALIZE3(m,p) : m##1(p##1), m##2(p##2), m##3(p##3)
#define PE_INITIALIZE4(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4)
#define PE_INITIALIZE5(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5)
#define PE_INITIALIZE6(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6)
#define PE_INITIALIZE7(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7)
#define PE_INITIALIZE8(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8)
#define PE_INITIALIZE9(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8), m##9(p##9)
#define PE_INITIALIZE10(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8), m##9(p##9), m##10(p##10)
#define PE_INITIALIZE11(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8), m##9(p##9), m##10(p##10), m##11(p##11)
#define PE_INITIALIZE12(m,p) : m##1(p##1), m##2(p##2), m##3(p##3), m##4(p##4), m##5(p##5), m##6(p##6), m##7(p##7), m##8(p##8), m##9(p##9), m##10(p##10), m##11(p##11), m##12(p##12)

#define PE_PARAMETER_CONSTRUCT0(p)
#define PE_PARAMETER_CONSTRUCT1(p) (p##1)
#define PE_PARAMETER_CONSTRUCT2(p) (p##1,p##2)
#define PE_PARAMETER_CONSTRUCT3(p) (p##1,p##2,p##3)
#define PE_PARAMETER_CONSTRUCT4(p) (p##1,p##2,p##3,p##4)
#define PE_PARAMETER_CONSTRUCT5(p) (p##1,p##2,p##3,p##4,p##5)
#define PE_PARAMETER_CONSTRUCT6(p) (p##1,p##2,p##3,p##4,p##5,p##6)
#define PE_PARAMETER_CONSTRUCT7(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7)
#define PE_PARAMETER_CONSTRUCT8(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8)
#define PE_PARAMETER_CONSTRUCT9(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9)
#define PE_PARAMETER_CONSTRUCT10(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9,p##10)
#define PE_PARAMETER_CONSTRUCT11(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9,p##10,p##11)
#define PE_PARAMETER_CONSTRUCT12(p) (p##1,p##2,p##3,p##4,p##5,p##6,p##7,p##8,p##9,p##10,p##11,p##12)


#ifdef PARALLEL_THREADS_PASSTHROUGH
#define PARALLEL_CRITICAL_SECTION int
inline void ParallelInitializeCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
}
inline void ParallelEnterCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
}
inline void ParallelLeaveCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
}
inline void ParallelDeleteCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
}
#endif

#ifdef PARALLEL_THREADS_WIN32
#define PARALLEL_CRITICAL_SECTION CRITICAL_SECTION
inline void ParallelInitializeCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
  InitializeCriticalSection(mutex);
}
inline void ParallelEnterCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
  EnterCriticalSection(mutex);
}
inline void ParallelLeaveCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
  LeaveCriticalSection(mutex);
}
inline void ParallelDeleteCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
  DeleteCriticalSection(mutex);
}
#endif

#ifdef PARALLEL_THREADS_PTHREAD
#define PARALLEL_CRITICAL_SECTION pthread_mutex_t
inline void ParallelInitializeCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
  pthread_mutex_init(mutex, NULL);
}
inline void ParallelEnterCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
  pthread_mutex_lock(mutex);
}
inline void ParallelLeaveCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
  pthread_mutex_unlock(mutex);
}
inline void ParallelDeleteCriticalSection(PARALLEL_CRITICAL_SECTION *mutex) {
}
#endif


// classes to hold a set of parameters - surely some template generalization of this exists somewhere
class ParallelExecutorParameters0 {
  public:
  ParallelExecutorParameters0() { }
};
#define PE_PARAMETER_CLASS(np)                                  \
  template PE_TEMPLATE_DEF##np                                  \
	class ParallelExecutorParameters##np {                        \
    public:                                                     \
		ParallelExecutorParameters##np(PE_PARAMETER_DEF##np(T,&_p))	\
    PE_INITIALIZE##np(p,_p) { }                                 \
		PE_MEMBER_DEF##np(T,&p)                                     \
	};


// class to hold the information about what we're executing to pass to the threadmain
class ParallelExecutorClassBase {
  public:
  int numThreads;
  int id;
  PARALLEL_CRITICAL_SECTION *cs;
  virtual void run() = 0;
};

#define PE_BASE_CLASS(np)                                           \
  template <typename Functor PE_TEMPLATE_ADD_DEF##np>               \
  class ParallelExecutorClassBase##np : public ParallelExecutorClassBase { \
    public:                                                         \
    Functor f;                                                      \
    ParallelExecutorParameters##np PE_TEMPLATE_CALL##np *params;    \
    virtual void run() {                                            \
      f(numThreads, id, cs PE_PARAMETER_ADD_CALL##np(params->p));   \
    }                                                               \
  };




// functions of type void* f(void*) - these are the threadmain's
#ifdef PARALLEL_THREADS_PASSTHROUGH
inline int ParallelExecutorStub(void *arg) {
  ParallelExecutorClassBase *pec = (ParallelExecutorClassBase*)arg;
  pec->run();
  return 0;
}
#endif
#ifdef PARALLEL_THREADS_WIN32
inline DWORD __stdcall ParallelExecutorStub(void *arg) {
  ParallelExecutorClassBase *pec = (ParallelExecutorClassBase*)arg;
  pec->run();
  return 0;
}
#endif
#ifdef PARALLEL_THREADS_PTHREAD
inline void* ParallelExecutorStub(void *arg) {
  ParallelExecutorClassBase *pec = (ParallelExecutorClassBase*)arg;
  pec->run();
  return 0;
}
#endif



// make sure you don't try to recursively call in parallel - I've already made this mistake once!
inline void ParallelExecutorCheck(bool start) {
  /*
  static int count = 0;
  static PARALLEL_CRITICAL_SECTION cs;
  ParallelInitializeCriticalSection(&cs);

  ParallelEnterCriticalSection(&cs);
  if (start) {
    if (count > 0)
      std::cerr<<"starting parallel executor recursively?!?!"<<std::endl;
    count++;
  } else {
    count--;
  }
  ParallelLeaveCriticalSection(&cs);
  ParallelDeleteCriticalSection(&cs);
  */
}


// the real work - create the set of threads, have them run, then join them back together
#ifdef PARALLEL_THREADS_PASSTHROUGH
#define PARALLEL_EXECUTOR(np)                                           \
  template <typename Functor PE_TEMPLATE_ADD_DEF##np>                   \
	void ParallelExecutor(const int reqThreads, Functor f PE_PARAMETER_ADD_DEF##np(T,&p)) { \
		ParallelExecutorCheck(true);                                        \
    ParallelExecutorParameters##np PE_TEMPLATE_CALL##np params PE_PARAMETER_CONSTRUCT##np(p); \
    std::vector< ParallelExecutorClassBase##np<Functor PE_TEMPLATE_ADD_CALL##np> > pec(1); \
    PARALLEL_CRITICAL_SECTION cs;                                       \
    ParallelInitializeCriticalSection(&cs);                             \
    pec[0].numThreads = reqThreads;                                     \
    pec[0].id = 0;                                                      \
    pec[0].cs = &cs;                                                    \
    pec[0].f = f;                                                       \
    pec[0].params = &params;                                            \
    ParallelExecutorStub(&pec[0]);                                      \
    ParallelDeleteCriticalSection(&cs);                                 \
    ParallelExecutorCheck(false);                                       \
  }
#endif

#ifdef PARALLEL_THREADS_WIN32
#define PARALLEL_EXECUTOR(np)                                           \
  template <typename Functor PE_TEMPLATE_ADD_DEF##np>                   \
	void ParallelExecutor(const int reqThreads, Functor f PE_PARAMETER_ADD_DEF##np(T,&p)) { \
		ParallelExecutorCheck(true);                                        \
    ParallelExecutorParameters##np PE_TEMPLATE_CALL##np params PE_PARAMETER_CONSTRUCT##np(p); \
    std::vector< ParallelExecutorClassBase##np<Functor PE_TEMPLATE_ADD_CALL##np> > pec(reqThreads); \
    PARALLEL_CRITICAL_SECTION cs;                                       \
    ParallelInitializeCriticalSection(&cs);                             \
    for (int i=0; i<reqThreads; i++) {                                  \
	    pec[i].numThreads = reqThreads;                                   \
	    pec[i].id = i;                                                    \
      pec[i].cs = &cs;                                                  \
	    pec[i].f = f;                                                     \
	    pec[i].params = &params;                                          \
    }                                                                   \
    std::vector<HANDLE> threads(reqThreads-1);                          \
    for (int i=1; i<reqThreads; i++) {                                  \
      threads[i-1] = CreateThread(NULL, 0, &ParallelExecutorStub, (void*)&pec[i], 0, NULL); \
    }                                                                   \
    ParallelExecutorStub(&pec[0]);                                      \
    for (int i=1; i<reqThreads; i++) {                                  \
      WaitForSingleObject(threads[i-1], INFINITE);                      \
      CloseHandle(threads[i-1]);                                        \
    }                                                                   \
    ParallelDeleteCriticalSection(&cs);                                 \
    ParallelExecutorCheck(false);                                       \
  }
#endif

#ifdef PARALLEL_THREADS_PTHREAD
#define PARALLEL_EXECUTOR(np)                                           \
  template <typename Functor PE_TEMPLATE_ADD_DEF##np>                   \
	void ParallelExecutor(const int reqThreads, Functor f PE_PARAMETER_ADD_DEF##np(T,&p)) { \
		ParallelExecutorCheck(true);                                        \
    ParallelExecutorParameters##np PE_TEMPLATE_CALL##np params PE_PARAMETER_CONSTRUCT##np(p); \
    std::vector< ParallelExecutorClassBase##np<Functor PE_TEMPLATE_ADD_CALL##np> > pec(reqThreads); \
    PARALLEL_CRITICAL_SECTION cs;                                       \
    ParallelInitializeCriticalSection(&cs);                             \
    for (int i=0; i<reqThreads; i++) {                                  \
	    pec[i].numThreads = reqThreads;                                   \
	    pec[i].id = i;                                                    \
      pec[i].cs = &cs;                                                  \
	    pec[i].f = f;                                                     \
	    pec[i].params = &params;                                          \
    }                                                                   \
    std::vector<pthread_t> threads(reqThreads-1);                       \
    for (int i=1; i<reqThreads; i++) {                                  \
      pthread_create(&threads[i-1], NULL, &ParallelExecutorStub, (void*)&pec[i]); \
    }                                                                   \
    ParallelExecutorStub(&pec[0]);                                      \
    for (int i=1; i<reqThreads; i++) {                                  \
	    int *status;                                                      \
      pthread_join(threads[i-1], (void**)&status);                      \
    }                                                                   \
    ParallelDeleteCriticalSection(&cs);                                 \
    ParallelExecutorCheck(false);                                       \
  }
#endif


// instantiate them for up to 5 parameters

PE_PARAMETER_CLASS(1)
PE_PARAMETER_CLASS(2)
PE_PARAMETER_CLASS(3)
PE_PARAMETER_CLASS(4)
PE_PARAMETER_CLASS(5)
PE_PARAMETER_CLASS(6)
PE_PARAMETER_CLASS(7)
PE_PARAMETER_CLASS(8)
PE_PARAMETER_CLASS(9)
PE_PARAMETER_CLASS(10)
PE_PARAMETER_CLASS(11)
PE_PARAMETER_CLASS(12)

PE_BASE_CLASS(0);
PE_BASE_CLASS(1);
PE_BASE_CLASS(2);
PE_BASE_CLASS(3);
PE_BASE_CLASS(4);
PE_BASE_CLASS(5);
PE_BASE_CLASS(6);
PE_BASE_CLASS(7);
PE_BASE_CLASS(8);
PE_BASE_CLASS(9);
PE_BASE_CLASS(10);
PE_BASE_CLASS(11);
PE_BASE_CLASS(12);

PARALLEL_EXECUTOR(0)
PARALLEL_EXECUTOR(1)
PARALLEL_EXECUTOR(2)
PARALLEL_EXECUTOR(3)
PARALLEL_EXECUTOR(4)
PARALLEL_EXECUTOR(5)
PARALLEL_EXECUTOR(6)
PARALLEL_EXECUTOR(7)
PARALLEL_EXECUTOR(8)
PARALLEL_EXECUTOR(9)
PARALLEL_EXECUTOR(10)
PARALLEL_EXECUTOR(11)
PARALLEL_EXECUTOR(12)






// I think the std::mem_fun stuff is similar, but doesn't take any number of parameters.
// to use with the parallel executor functions, you need instantiation for 2 extra
// variables since the thread id and num threads parameters are not treated as special.
template <typename Class, typename Functor>
class ClassFunctor {
 public:
  ClassFunctor() { }
  ClassFunctor(Class *_c, Functor _f) : c(_c), f(_f) { }
  
  void operator()() { (c->*f)(); }


#define PE_CLASS_FUNCTOR(np)                      \
  template PE_TEMPLATE_DEF##np                    \
	void operator() (PE_PARAMETER_DEF##np(T,&p)) {  \
    (c->*f)(PE_PARAMETER_CALL##np(p));            \
  }                                               \

PE_CLASS_FUNCTOR(1)
PE_CLASS_FUNCTOR(2)
PE_CLASS_FUNCTOR(3)
PE_CLASS_FUNCTOR(4)
PE_CLASS_FUNCTOR(5)
PE_CLASS_FUNCTOR(6)
PE_CLASS_FUNCTOR(7)
PE_CLASS_FUNCTOR(8)
PE_CLASS_FUNCTOR(9)
PE_CLASS_FUNCTOR(10)
PE_CLASS_FUNCTOR(11)
PE_CLASS_FUNCTOR(12)


 private:
  Class *c;
  Functor f;
};


// this will autodetect/deduce the types, so it's easy to place inline (create a temp object) without
// having to figure out what the damn syntax is for member function pointers
template <typename Class, typename Functor>
ClassFunctor<Class, Functor> makeClassFunctor(Class *c, Functor f) {
  return ClassFunctor<Class,Functor>(c,f);
}


// wait for all threads to enter this function before any leave
inline void ParallelSyncThreads(int nt, int id, PARALLEL_CRITICAL_SECTION *mutex, volatile int *syncToken) {

  while ((*syncToken) != id) { }
  (*syncToken)++;

  while ((*syncToken) != id+nt) { }
  (*syncToken)++;

  if (id==nt-1)
    (*syncToken) = 0;
}





class ParallelThreadPool {
  public:

  ParallelThreadPool() {
    mStop = false;
    mSyncToken = 0;
    ParallelInitializeCriticalSection(&mMutex);
  }
  ~ParallelThreadPool() {
    ParallelDeleteCriticalSection(&mMutex);
  }

  int GetNumThreads() const { return (int)mThreads.size() + 1; }
  
  void StartThreads(int nt) {
#ifdef PARALLEL_THREADS_PASSTHROUGH
    nt = 1;
#endif

    mThreads.resize(nt-1);
    mRunInfo.resize(nt-1);
    mRun = new volatile int[nt-1];
    mStarters.resize(nt-1);

    for (int i=0; i<nt-1; i++) {
      mRunInfo[i] = NULL;
      mRun[i] = 0;
      mStarters[i].pool = this;
      mStarters[i].id = i+1;
#ifdef PARALLEL_THREADS_PASSTHROUGH
      mThreads[i] = 0;
#endif
#ifdef PARALLEL_THREADS_WIN32
      mThreads[i] = CreateThread(NULL, 0, &SThreadMain, (void*)&mStarters[i], 0, NULL);
#endif
#ifdef PARALLEL_THREADS_PTHREADS
      pthread_create(&mThreads[i], NULL, &SThreadMain, (void*)&mStarters[i]);
#endif
    }
  }


  void StopThreads() {
    mStop = true;
    for (int i=0; i<(int)mThreads.size(); i++) {
#ifdef PARALLEL_THREADS_PASSTHROUGH
      // do nothing
#endif
#ifdef PARALLEL_THREADS_WIN32
      WaitForSingleObject(mThreads[i], INFINITE);
      CloseHandle(mThreads[i]);
#endif
#ifdef PARALLEL_THREADS_PTHREADS
      int *status;
      pthread_join(mThreads[i], (void**)&status);
#endif
    }

    mThreads.clear();
    mRunInfo.clear();
    delete [] mRun;
    mStarters.clear();
  }


#define PARALLEL_THREAD_RUN(np)                                         \
  template <typename Functor PE_TEMPLATE_ADD_DEF##np>                   \
  void ParallelRun(Functor f PE_PARAMETER_ADD_DEF##np(T,&p)) {          \
    ParallelExecutorParameters##np PE_TEMPLATE_CALL##np params PE_PARAMETER_CONSTRUCT##np(p); \
    for (int i=0; i<(int)mThreads.size(); i++) {                        \
      ParallelExecutorClassBase##np<Functor PE_TEMPLATE_ADD_CALL##np> *pec = new ParallelExecutorClassBase##np<Functor PE_TEMPLATE_ADD_CALL##np>(); \
      pec->numThreads = (int)mThreads.size()+1;                         \
      pec->id = i+1;                                                    \
      pec->cs = &mMutex;                                                \
      pec->f = f;                                                       \
      pec->params = &params;                                            \
      mRunInfo[i] = pec;                                                \
      mRun[i] = 1;                                                      \
    }                                                                   \
    ParallelExecutorClassBase##np<Functor PE_TEMPLATE_ADD_CALL##np> pec0; \
    pec0.numThreads = (int)mThreads.size()+1;                           \
    pec0.id = 0;                                                        \
    pec0.cs = &mMutex;                                                  \
    pec0.f = f;                                                         \
    pec0.params = &params;                                              \
    ParallelExecutorStub((void*)&pec0);                                 \
    for (int i=0; i<(int)mThreads.size(); i++) {                        \
      while (1) {                                                       \
        if (mRun[i] == 0) break;                                        \
      }                                                                 \
    }                                                                   \
  }
  PARALLEL_THREAD_RUN(0);
  PARALLEL_THREAD_RUN(1);
  PARALLEL_THREAD_RUN(2);
  PARALLEL_THREAD_RUN(3);
  PARALLEL_THREAD_RUN(4);
  PARALLEL_THREAD_RUN(5);
  PARALLEL_THREAD_RUN(6);
  PARALLEL_THREAD_RUN(7);
  PARALLEL_THREAD_RUN(8);
  PARALLEL_THREAD_RUN(9);
  PARALLEL_THREAD_RUN(10);
  PARALLEL_THREAD_RUN(11);
  PARALLEL_THREAD_RUN(12);


private:

  class ParallelThreadPoolStarter {
    public:
    ParallelThreadPool *pool;
    int id;
  };
#ifdef PARALLEL_THREADS_PASSTHROUGH
  static int SThreadMain(void *_starter) {
#endif
#ifdef PARALLEL_THREADS_WIN32
  static DWORD __stdcall SThreadMain(void *_starter) {
#endif
#ifdef PARALLEL_THREADS_PTHREADS
  static void* SThreadMain(void *_starter) {
#endif
    ParallelThreadPoolStarter *starter = (ParallelThreadPoolStarter*)_starter;
    starter->pool->ThreadMain(starter->id);
    return 0;
  }

  void ThreadMain(int id) {
    while (!mStop) {
      if (mRun[id-1] != 0) {
        ParallelExecutorStub(mRunInfo[id-1]);
        delete mRunInfo[id-1];
        mRunInfo[id-1] = NULL;
        mRun[id-1] = 0;
      }
    }
  }



#ifdef PARALLEL_THREADS_PASSTHROUGH
  std::vector<int> mThreads;
#endif
#ifdef PARALLEL_THREADS_WIN32
  std::vector<HANDLE> mThreads;
#endif
#ifdef PARALLEL_THREADS_PTHREADS
  std::vector<pthread_t> mThreads;
#endif

  PARALLEL_CRITICAL_SECTION mMutex;
  volatile int mSyncToken;
  volatile int mStop;

  volatile int *mRun;
  std::vector<ParallelExecutorClassBase*> mRunInfo;
  std::vector<ParallelThreadPoolStarter> mStarters;
};




extern int idealNumThreads;



#endif

