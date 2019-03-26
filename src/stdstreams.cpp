#include <stdarg.h>
#include <stdstreams.h>
#include <common.h>
#include <utils.h>
#include <time.h>
#include <unistd.h>

static bool std_streams_class_created = false;

StdStreams pStdLogs;

StdStreams::StdStreams():
  stdout_f(stdout), stderr_f(stderr)
{
  if (std_streams_class_created)
    THROW_EXCEPTION(std::logic_error, "Second instance of StdStreams class is not allowed");

  std_streams_class_created = true;

  std::ios_base::sync_with_stdio(true);

  stdout_dup = dup(STDOUT_FILENO);
  stderr_dup = dup(STDERR_FILENO);
}

StdStreams::~StdStreams()
{
  CloseAll();
}

void StdStreams::SetStdOut(const std::string &file, const std::string &attr)
{
  FILE *new_stdout = freopen(file.c_str(), attr.c_str(), stdout);
  ASSERT_FILE(new_stdout, file.c_str());

  if (fileno(stdout_f) != STDOUT_FILENO)
    fclose(stdout_f);

  stdout_f = stdout;
}

void StdStreams::ResetStdOut()
{
  if (fileno(stdout_f) != STDOUT_FILENO)
    fclose(stdout_f);

  dup2(stdout_dup, STDOUT_FILENO);
  stdout = fdopen(STDOUT_FILENO, "w");
  stdout_f = stdout;
}

void StdStreams::SetStdErr(const std::string &file, const std::string &attr)
{
  FILE *new_stderr = freopen(file.c_str(), attr.c_str(), stderr);
  ASSERT_FILE(new_stderr, file.c_str());

  if (fileno(stderr_f) != STDERR_FILENO)
    fclose(stderr_f);

  stderr_f = stderr;
}

void StdStreams::ResetStdErr()
{
  if (fileno(stderr_f) != STDERR_FILENO)
    fclose(stderr_f);

  dup2(stderr_dup, STDERR_FILENO);
  stderr = fdopen(STDERR_FILENO, "w");
  stderr_f = stderr;
}

void StdStreams::Write(const char *fmt, ...)
{
  va_list va;
  va_start(va, fmt);
  log.WriteVf(fmt, va);
  va_end(va);
}

void StdStreams::WriteError(const std::string &file, int line, const char *fmt, ...)
{
  va_list va;
  va_start(va, fmt);
  log.WriteErrorVf(file, line, fmt, va);
  va_end(va);
}

void StdStreams::CloseAll()
{
  fclose(stdout_f);
  fclose(stderr_f);
}
