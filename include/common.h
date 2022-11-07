#ifndef COMMON_H_
#define COMMON_H_

#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <omp.h>
#include <stdexcept>
#include <map>
#include <memory>
#include <algorithm>
#include <functional>
#include <string>
#include <cstring>

#include "boost/container/vector.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/make_shared.hpp"
#include "our_make_shared.hpp"

#ifdef WITH_MPI_
#include <mpi.h>
#endif

#define VECTOR          boost::container::vector
#define SHARED_PTR      boost::shared_ptr
#define MAKE_SHARED     boost::make_shared
#define ALLOCATE_SHARED boost::allocate_shared

#if __cplusplus < 201103L
template<typename T> std::string to_string(const T &val)
{
  std::ostringstream of;
  of << val;
  return of.str();
}
#define TO_STRING     to_string
#else
#define TO_STRING     std::to_string
#endif

#include "types.h"
#include "logs.h"
#include "stdstreams.h"
#include "datadir.h"
#include "utils.h"
#include "rand_num_generators.h"
#include "common_parameters.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef MAX
#define MAX(_X,_Y) ((_X) > (_Y) ? (_X) : (_Y))
#endif

#ifndef MIN
#define MIN(_X,_Y) ((_X) < (_Y) ? (_X) : (_Y))
#endif

#ifndef SQR
#define SQR(_X) ((_X)*(_X))
#endif

#define SIGNRE(_x)  (creal(_x)>0? 1.0 : (creal(_x)<0? -1.0: 0.0))

#define SIGN(_x)    ((_x)>0? 1.0 : ((_x)<0? -1.0: 0.0))

class AssertFailure : public std::runtime_error
{
public:
  AssertFailure(const std::string &msg) : std::runtime_error(msg) {}
  AssertFailure(const char *msg) : std::runtime_error(msg) {}
};

class FileFailure : public std::runtime_error
{
public:
  FileFailure(const std::string &msg) : std::runtime_error(msg) {}
  FileFailure(const char *msg) : std::runtime_error(msg) {}
};

class MPIFailure : public std::runtime_error
{
public:
  MPIFailure(const std::string &msg) : std::runtime_error(msg) {}
  MPIFailure(const char *msg) : std::runtime_error(msg) {}
};

#define THROW_EXCEPTION(exception, msg) \
    throw exception("[" + std::string(__FILE__) + ":" + TO_STRING(__LINE__) + "]: " + msg);

#define THROW_EXCEPTION_VERB(exception, msg)               \
  {                                                        \
    std::string __exception_msg = "[" + std::string(__FILE__) + ":" + TO_STRING(__LINE__) + "]: " + msg; \
    std::string __log_msg = std::string("EXCEPTION THROWN: ") +__exception_msg; \
    std::string __error_log_msg = std::string("EXCEPTION THROWN: ") + msg; \
    pStdLogs.Write(__log_msg.c_str()); \
    WRITE_ERROR(__error_log_msg.c_str());               \
    throw exception(__exception_msg);                       \
  }

#define ASSERT(x)                                       \
  if (!(x))                                             \
  {                                                     \
    THROW_EXCEPTION(AssertFailure, "Assert failed.");   \
  }

#define ASSERTVERB(x, msg)                              \
  if (!(x))                                             \
  {                                                     \
    THROW_EXCEPTION(AssertFailure, msg);                \
  }

#ifdef ENABLE_DEBUG_ASSERT
#define DEBUG_ASSERT(x)  ASSERT(x);
#else
#define DEBUG_ASSERT(x)  { (void) 0; }
#endif

#define ASSERT_FILE(x, file_name)                                               \
  if ((x) == NULL)                                                              \
  {                                                                             \
    THROW_EXCEPTION(FileFailure, "Can't open file " + std::string(file_name));  \
  }

#define SAFE_FOPEN(file, file_name, spec)  \
  {                                        \
    file = fopen(file_name, spec);         \
    ASSERT_FILE(file, file_name);          \
  }

#define SAFE_FSEEK(f, pos, origin)                                      \
  {                                                                     \
    if(fseek((f), (pos), (origin)) != 0)              \
      THROW_EXCEPTION(FileFailure, "Failed to perform seek operation in file"); \
  }

#define SAFE_FWRITE(ptr, size, num, f)                                  \
  {                                                                     \
    if(fwrite((ptr), (size), (num), (f)) != (size_t)(num))              \
      THROW_EXCEPTION(FileFailure, "Failed to write data to the file"); \
  }

#define SAFE_FREAD(ptr, size, num, f)                                    \
  {                                                                      \
    if(fread((ptr), (size), (num), (f)) != (size_t)(num))                \
      THROW_EXCEPTION(FileFailure, "Failed to read data from the file"); \
  }

#define SAFE_FPRINTF(f, p, ...)                                          \
  {                                                                      \
    if(fprintf((f), (p), ##__VA_ARGS__) < 0)                             \
      THROW_EXCEPTION(FileFailure, "Failed to write data to the file");  \
  }

#define SAFE_FSCANF(f, p, ...)                                           \
  {                                                                      \
    if(fscanf((f), (p), ##__VA_ARGS__) < 0)                              \
      THROW_EXCEPTION(FileFailure, "Failed to read data from the file"); \
  }

#ifdef WITH_MPI_
#define SAFE_MPI_CALL(f)                                 \
  {                                                      \
    int mpi_err = (f);                                   \
    if(mpi_err != MPI_SUCCESS)                           \
      THROW_EXCEPTION(MPIFailure, "MPI routine failed"); \
  }
#else
#define SAFE_MPI_CALL(f) { }
#endif

void
common_AppInit(int argc, char **argv, const std::string &app_name);

void
common_AppFin();

#endif /* COMMON_H_ */
