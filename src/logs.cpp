#include <stdarg.h>
#include <logs.h>
#include <common.h>
#include <utils.h>
#include <time.h>
#include <unistd.h>

Logs::Logs() : 
  stdout_f(stdout), stderr_f(stderr)
{
}

Logs::~Logs()
{
  CloseAll();
}

void Logs::PrintErrorHeader(FILE *f)
{
  fprintf(f, "ERROR [");
  Utils::PrintDateTime(f);
  fprintf(f, "]");
}

void Logs::PrintErrorHeader(const std::string &file, int line, FILE *f)
{
  PrintErrorHeader(f);
  fprintf(f, " (%s:%d) ", file.c_str(), line);
}

void Logs::SetInfoOut(const std::string &file, const std::string &attr)
{
  FILE *new_f;
  SAFE_FOPEN(new_f, file.c_str(), attr.c_str());

  stdout_f = new_f;
}

void Logs::ResetInfoOut()
{
  fclose(stdout_f);
  stdout_f = stdout;
}

void Logs::SetErrOut(const std::string &file, const std::string &attr)
{
  FILE *new_f;
  SAFE_FOPEN(new_f, file.c_str(), attr.c_str());

  stdout_f = new_f;
}

void Logs::ResetErrOut()
{
  fclose(stderr_f);
  stderr_f = stderr;
}

void Logs::WriteVf(const char *fmt, va_list va)
{
  vfprintf(stdout_f, fmt, va);
  fflush(stdout_f);
}

void Logs::WriteErrorVf(const char *fmt, va_list va)
{
  PrintErrorHeader(stderr_f);
  vfprintf(stderr_f, fmt, va);
  fprintf(stderr_f, "\n");
  fflush(stderr_f);

  PrintErrorHeader(stdout_f);
  vfprintf(stdout_f, fmt, va);
  fprintf(stdout_f, "\n");
  fflush(stdout_f);
}

void Logs::WriteErrorVf(const std::string &file, int line, const char *fmt, va_list va)
{
  PrintErrorHeader(file, line, stderr_f);
  vfprintf(stderr_f, fmt, va);
  fprintf(stderr_f, "\n");
  fflush(stderr_f);

  PrintErrorHeader(file, line, stdout_f);
  vfprintf(stdout_f, fmt, va);
  fprintf(stdout_f, "\n");
  fflush(stdout_f);
}

void Logs::Write(const char *fmt, ...)
{
  va_list va;
  va_start(va, fmt);
  WriteVf(fmt, va);
  va_end(va);
}

void Logs::WriteError(const char *fmt, ...)
{
  va_list va;
  va_start(va, fmt);
  WriteErrorVf(fmt, va);
  va_end(va);
}

void Logs::WriteError(const std::string &file, int line, const char *fmt, ...)
{
  va_list va;
  va_start(va, fmt);
  WriteErrorVf(file, line, fmt, va);
  va_end(va);
}

void Logs::CloseAll()
{
  fclose(stdout_f);
  fclose(stderr_f);
}
