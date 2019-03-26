#include <datadir.h>
#include <common.h>
#include <utils.h>

using std::string;

DataDir::DataDir():
  default_dir("."), is_overwrite(true), dir(default_dir), suff(""), pref("")
{
}

DataDir::~DataDir()
{
}

bool DataDir::GetOverwrite() const
{
  return is_overwrite;
}

void DataDir::SetOverwrite(bool overwrite)
{
  is_overwrite = overwrite;
}

string DataDir::DefaultDirectory() const
{
  return default_dir;
}

string DataDir::GetDirectory() const
{
  return dir;
}

void DataDir::SetDirectory(const string &directory)
{
  if (directory.empty())
  {
    dir = default_dir;
  }
  else
  {
    dir = directory;
  }
}

string DataDir::GetSuffix() const
{
  return suff;
}

void DataDir::SetSuffix(const string &suffix)
{
  suff = suffix;
}

string DataDir::GetPrefix() const
{
  return pref;
}

void DataDir::SetPrefix(const string &prefix)
{
  pref = prefix;
}

string DataDir::GetFullPath(const string &name) const
{
  const string c_dot(".");

  string fname = (pref.empty()) ? name : pref + c_dot + name;
  fname = (suff.empty()) ? fname : fname + c_dot + suff;
  string fpath = Utils::GetAbsFileName(dir, fname);

  return fpath;
}

FILE* DataDir::OpenFile(const string &name, const string &attr) const
{
  ASSERT(!attr.empty());

  FILE *f;

  string fpath = GetFullPath(name);
  SAFE_FOPEN(f, fpath.c_str(), attr.c_str());

  return f;
}
