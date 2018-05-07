#ifndef SHAREFILE_INCLUDED
#define SHAREFILE_INCLUDED
#ifdef  MAIN_FILE
int global;
int nElemglobal;
int rownumglobal;
#else
extern int global;
extern int nElemglobal;
extern int rownumglobal;
#endif
#endif