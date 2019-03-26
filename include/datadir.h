#ifndef DATADIR_H_
#define DATADIR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>

class DataDir
{
private:
  const std::string default_dir;

  bool is_overwrite;
  std::string dir;
  std::string suff;
  std::string pref;

public:
  DataDir();
  ~DataDir();

  bool GetOverwrite() const;
  void SetOverwrite(bool overwrite);

  std::string DefaultDirectory() const;

  std::string GetDirectory() const;
  void SetDirectory(const std::string &directory);

  std::string GetSuffix() const;
  void SetSuffix(const std::string &suffix);

  std::string GetPrefix() const;
  void SetPrefix(const std::string &prefix);

  std::string GetFullPath(const std::string &name) const;
  FILE* OpenFile(const std::string &name, const std::string &attr) const;
};

#endif /* DATADIR_H_ */
