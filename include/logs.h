#ifndef LOGS_H_
#define LOGS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <stdarg.h>

class Logs
{
private:
  FILE *stdout_f;
  FILE *stderr_f;

  void PrintErrorHeader(FILE *f);
  void PrintErrorHeader(const std::string &file, int line, FILE *f);

public:
  Logs();
  ~Logs();

  void SetInfoOut(const std::string &file, const std::string &attr);
  void ResetInfoOut();
  void SetErrOut(const std::string &file, const std::string &attr);
  void ResetErrOut();

  void Write(const char *fmt, ...);
  void WriteError(const char *fmt, ...);
  void WriteError(const std::string &file, int line, const char *fmt, ...);

  void WriteVf(const char *fmt, va_list va);
  void WriteErrorVf(const char *fmt, va_list va);
  void WriteErrorVf(const std::string &file, int line, const char *fmt, va_list va);

  void CloseAll();
};

#endif /* LOGS_H_ */
