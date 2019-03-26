#include <utils.h>
#include <string.h>
#include <string>
#include <ctime>
#include <iomanip>
#include <sstream>

#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

using std::string;
using std::cout;
using std::endl;

off_t 
Utils::GetFileSize(FILE *file)
{
  off_t pos = ftell(file);
  fseek(file, 0, SEEK_END);
  off_t size = ftell(file);
  fseek(file, pos, SEEK_SET);
  return size;
}

string
Utils::ReadStringFromFile(FILE *file, const int max_length)
{
  string res;

  if (feof(file))
    return res;

  int count = 0;

  char c = fgetc(file);
  while(c != '\r' && c != '\n')
  {
    if (count >= max_length)
      return res;

    res = res + string(1, c);

    if (feof(file))
      break;

    c = fgetc(file);
  }

  return res;
}

bool
Utils::IsAbsolutePath(const string &str)
{
  if (str.empty())
  {
    return false;
  }
  else
  {
    return (str[0] == '/');
  }
}

string
Utils::ConcatPaths(const string &dir, const string &name)
{
  uint c1 = dir.length();

  while (c1 > 0 && dir[c1 - 1] == '/')
    c1--;

  string res = dir.substr(0, c1) + string("/") + name;

  return res;
}

string
Utils::GetAbsFileName(const string &dir, const string &name)
{
  if (Utils::IsAbsolutePath(name))
    return name;
  else
    return Utils::ConcatPaths(dir, name);
}

string
Utils::GetDateTime()
{
  const uint buf_size = 200;

  char datetime[buf_size + 1];

  time_t t;
  time(&t);

  strftime(datetime, buf_size, "%d-%m-%Y %H:%M:%S", localtime (&t));

  string datestr(datetime);

  return datestr;
}

void
Utils::PrintDateTime(FILE *file)
{
  string str = Utils::GetDateTime();
  fprintf (file, "%s", str.c_str());
}

void
Utils::PrintGreetings()
{
  cout << "--------------------------------------------------" << endl;
  cout << "-                Program started                 -" << endl;
  cout << "--------------------------------------------------" << endl;
  cout << "Start time: " << Utils::GetDateTime();
  cout << endl << endl;
}

void
Utils::PrintEnd()
{
  cout << endl << endl;
  cout << "End time: " << Utils::GetDateTime() << endl;
  cout << "--------------------------------------------------" << endl;
  cout << "-                Program completed               -" << endl;
  cout << "--------------------------------------------------" << endl;
}

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
 * windows and linux. */

int64_t
Utils::GetTimeMs64()
{
#ifdef WIN32
  /* Windows */
  FILETIME ft;
  LARGE_INTEGER li;

  /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
    * to a LARGE_INTEGER structure. */
  GetSystemTimeAsFileTime(&ft);
  li.LowPart = ft.dwLowDateTime;
  li.HighPart = ft.dwHighDateTime;

  int64_t ret = li.QuadPart;
  ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
  ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

  return ret;
#else
  /* Linux */
  struct timeval tv;

  gettimeofday(&tv, NULL);

  int64_t ret = tv.tv_usec;
  /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
  ret /= 1000;

  /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
  ret += (tv.tv_sec * 1000);

  return ret;
#endif
}

int64_t
Utils::GetTimeMks64()
{
  const int64_t c_billion = 1000000000L;
  int64_t ret;

  struct timespec tp;
  ASSERT(clock_gettime(CLOCK_REALTIME, &tp) == 0);
  ret = tp.tv_sec * c_billion + tp.tv_nsec;

  return ret;
}

bool
Utils::CopyFile(const string &src, const string &dst)
{
  if (src.empty() || dst.empty())
    return false;

  const int buf_size = 4096;

  size_t nread, nwritten;

  bool ret = true;
  char *buf = new char[buf_size];

  FILE *fsrc, *fdst;

  SAFE_FOPEN(fsrc, src.c_str(), "rb");
  SAFE_FOPEN(fdst, dst.c_str(), "wb");

  nread = fread(buf, sizeof(char), buf_size, fsrc);
  if (ferror(fsrc))
  {
    nread = 0;
    ret = false;
  }

  while(nread > 0)
  {
    nwritten = fwrite(buf, sizeof(char), nread, fdst);
    if (nwritten != nread)
    {
      ret = false;
      break;
    }

    nread = fread(buf, sizeof(char), buf_size, fsrc);
    if (ferror(fsrc))
    {
      ret = false;
      break;
    }
  }

  fclose(fsrc);
  fclose(fdst);

  delete buf;

  return ret;
}

string
Utils::FileNameIfSet(const string &fname)
{
  if (fname.empty())
    return "(none)";
  else
    return fname;
}
