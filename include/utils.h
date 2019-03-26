#ifndef UTILS_H_
#define UTILS_H_

#include <common.h>

class Utils
{
private:
  Utils() {}

public:
  static off_t GetFileSize(FILE *file);
  static std::string ReadStringFromFile(FILE *file, const int max_length = 65536);
  static bool IsAbsolutePath(const std::string &str);
  static std::string ConcatPaths(const std::string &p1, const std::string &p2);
  static std::string GetAbsFileName(const std::string &name, const std::string &dir);
  static std::string GetDateTime();
  static void PrintDateTime(FILE *file);
  static void PrintGreetings();
  static void PrintEnd();
  static int64_t GetTimeMs64();
  static int64_t GetTimeMks64();
  static bool CopyFile(const std::string &src, const std::string &dst);
  static std::string FileNameIfSet(const std::string &fname);

  template <typename T> static void Chop(VECTOR<T> A, const double prec = 1e-8);
};

template <typename T> void Utils::Chop(VECTOR<T> A, const double prec)
{
  for(uint i = 0; i < A.size(); i++)
  {
    if (abs(A[i]) < prec)
      A[i] = 0;
  }
};

#endif /* UTILS_H_ */
