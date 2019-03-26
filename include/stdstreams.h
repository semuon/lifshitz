#ifndef STDLOGS_H_
#define STDLOGS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <logs.h>

class StdStreams
{
private:
  FILE *stdout_f;
  FILE *stderr_f;
  int stdout_dup;
  int stderr_dup;
  Logs log;

public:
  StdStreams();
  ~StdStreams();

  void SetStdOut(const std::string &file, const std::string &attr);
  void ResetStdOut();
  void SetStdErr(const std::string &file, const std::string &attr);
  void ResetStdErr();

  void Write(const char *fmt, ...);
  void WriteError(const std::string &file, int line, const char *fmt, ...);

  void CloseAll();
};

// Global instance of logs
extern StdStreams pStdLogs;

#define WRITE_ERROR(fmt, ...)  pStdLogs.WriteError(__FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define TERMINATE(fmt, ...)    { WRITE_ERROR(fmt, ##__VA_ARGS__); std::terminate(); }

#endif /* LOGS_H_ */
