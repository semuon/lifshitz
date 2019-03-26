#ifndef PROFILER_H_
#define PROFILER_H_

#include <common.h>

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                           Outline                                //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

/*
class Timer
{
public:
  Timer(const std::string &name);
  ~Timer();

  void Start();
  void Stop();

  void PrintStatistics() const;

};
*/

/*
class Profiler
{
public:
  Profiler(const std::string &name);
  ~Profiler();

  void StartSection(const std::string &name);
  void EndSection(const std::string &name);

  void Reset();

  void StartTimer(const std::string &name);
  void StopTimer(const std::string &name);
  void IncrementCounter(const std::string &name);
  template<typename T> void SaveValue(const std::string &name, const T &value);

  void PrintStatistics() const;

};
*/

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                          Internal classes                        //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

/*
class IValue
{
public:
  virtual ~IValue() {}

  virtual void PrintStatistics() const {}

};
*/

/*
template <typename T>
class Value : public IValue
{
public:
  Value(const std::string &name);
  ~Value();

  void AddValue(const T &value);

  void PrintStatistics() const;

};
*/

/*
class ProfilingSection
{
public:
  ProfilingSection(const std::string &name);
  ~ProfilingSection();

  void StartTimer(const std::string &name);
  void StopTimer(const std::string &name);

  template<typename T> void SaveValue(const std::string &name, const T &value);

  void IncrementCounter(const std::string &name);

  void StopEverything();

  void PrintStatistics() const;

};
*/

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                             Timer                                //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

class Timer
{
private:
  const double scaling_to_ms; // = 1000000.0;
  int64_t begin;
  int64_t end;
  int64_t interval;
  int64_t max;
  int64_t min;
  int64_t tot;
  int64_t num;
  bool is_running;
  std::string m_name;

public:
  Timer(const std::string &name) : scaling_to_ms(1000000.0),
    begin(0), end(0), interval(0), max(0),
    min(-1), tot(0), num(0), is_running(false),
    m_name(name)
  {
  }

  ~Timer()
  {
    //if (is_running) WRITE_ERROR("WARNING: Timer \"%s\" was running when destroyed!", m_name.c_str());
  }

  void Start()
  {
    if (!is_running)
    {
      num++;
      begin = Utils::GetTimeMks64();
      is_running = true;
    }
  }

  void Stop()
  {
    if (is_running)
    {
      end = Utils::GetTimeMks64();
      interval = end - begin;
      tot += interval;
      if (max < interval) max = interval;
      if (min > interval || min < 0) min = interval;
      is_running = false;
    }
  }

  void PrintStatistics() const
  {
    std::cout << "      Timer \"" << m_name << "\":" << std::endl;
    if (is_running) std::cout << "        WARNING: Timer is still running!" << std::endl;
    std::cout << "        Times called: " << num << std::endl;
    std::cout << "        Total execution time: " << tot / scaling_to_ms << " ms" << std::endl;
    std::cout << "        Mean execution time: " << (tot / scaling_to_ms) / num << " ms" << std::endl;
    std::cout << "        Maximal execution time: " << (max / scaling_to_ms) << " ms" << std::endl;
    std::cout << "        Minimal execution time: " << (min / scaling_to_ms) << " ms" << std::endl;
  }
};

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                             Values                               //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

class IValue
{
public:
  virtual ~IValue() {}

  virtual void PrintStatistics() const {}
};

template <typename T>
class Value : public IValue
{
private:
  std::string m_name;
  VECTOR<T> vals;

public:
  Value(const std::string &name):
  m_name(name), vals()
  {
  }

  ~Value()
  {
  }

  void AddValue(const T &value)
  {
    vals.push_back(value);
  }

  void PrintStatistics() const
  {
    if (vals.size() > 0)
    {
      std::cout << "      Saved values of \"" << m_name << "\":" << std::endl;
      std::cout << "        ";
      for(uint i = 0; i < vals.size(); i++)
        std::cout << vals[i] << " ";
       std::cout << std::endl;
      std::cout << "        Total number: " << vals.size() << std::endl;
    }
    else
    {
      std::cout << "      No saved values of \"" << m_name << "\"" << std::endl;
    }
  }
};

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                        Profiling section                         //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

class ProfilingSection
{
private:
  typedef SHARED_PTR<IValue> IValuePtr;

  std::string m_name;
  std::map<std::string, Timer> timers;
  std::map<std::string, int64_t> counters;
  std::map<std::string, IValuePtr> values;

  template<typename T> T& GetItem(std::map<std::string, T> &tmap, const std::string &name)
  {
    typename std::map<std::string, T>::iterator item = tmap.find(name);
    if (item == tmap.end())
    {
      std::pair< typename std::map<std::string, T>::iterator, bool > pair = tmap.insert(typename std::map<std::string, T>::value_type (name, T(name) ));
      item = pair.first;
    }

    return item->second;
  }

  template<typename T> IValuePtr& GetItemIValue(std::map<std::string, IValuePtr> &tmap, const std::string &name)
  {
    typename std::map<std::string, IValuePtr>::iterator item = tmap.find(name);
    if (item == tmap.end())
    {
      IValuePtr new_val(new Value<T>(name));
      std::pair< typename std::map<std::string, IValuePtr>::iterator, bool > pair = tmap.insert(typename std::map<std::string, IValuePtr>::value_type (name, new_val));
      item = pair.first;
    }

    return item->second;
  }

public:
  ProfilingSection(const std::string &name) :
  m_name(name), timers(), counters(), values()
  {
  }

  ~ProfilingSection()
  {
  }

  void StartTimer(const std::string &name)
  {
    Timer &timer = GetItem(timers, name);
    timer.Start();
  }

  void StopTimer(const std::string &name)
  {
    Timer &timer = GetItem(timers, name);
    timer.Stop();
  }

  template<typename T> void SaveValue(const std::string &name, const T &value)
  {
    IValuePtr &val = GetItemIValue<T>(values, name);
    Value<T> &valt = *(static_cast<Value<T>*>(val.get()));
    valt.AddValue(value);
  }

  void IncrementCounter(const std::string &name)
  {
    int64_t &ctr = counters[name];
    ctr++;
  }

  void StopEverything()
  {
    for(std::map<std::string, Timer>::iterator it = timers.begin(); it != timers.end(); ++it)
      it->second.Stop();
  }

  void PrintStatistics() const
  {
    std::cout << "    Profiling section \"" << m_name << "\":" << std::endl;
    for(std::map<std::string, Timer>::const_iterator it = timers.begin(); it != timers.end(); ++it)
      it->second.PrintStatistics();

    for(std::map<std::string, IValuePtr>::const_iterator it = values.begin(); it != values.end(); ++it)
      it->second.get()->PrintStatistics();

    if (counters.size() > 0)
    {
      std::cout << "      Counters:" << std::endl;
      for(std::map<std::string, int64_t>::const_iterator it = counters.begin(); it != counters.end(); ++it)
        std::cout << "        " << it->first << ": " << it->second << std::endl;
    }
    else
    {
      std::cout << "      No counters" << std::endl;
    }
  }
};

// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
//                            Profiler                              //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //
// ---------------------------------------------------------------- //

class Profiler
{
private:
  typedef std::map<std::string, ProfilingSection> ProfMap;
  typedef std::map<std::string, ProfMap::iterator> ActiveProfMap;

  std::string m_name;
  ProfMap sections;
  ActiveProfMap active_sections;

public:
  Profiler(const std::string &name) :
  m_name(name), sections(), active_sections()
  {
  }

  ~Profiler()
  {
  }

  void StartSection(const std::string &name)
  {
    ProfMap::iterator item = sections.find(name);
    if (item == sections.end())
    {
      std::pair< ProfMap::iterator, bool > pair = sections.insert( ProfMap::value_type (name, ProfilingSection(name) ));
      item = pair.first;
    }

    ActiveProfMap::iterator active_item = active_sections.find(name);
    if (active_item == active_sections.end())
      active_sections.insert( ActiveProfMap::value_type (name, item ));

    const std::string timer_name = "Section \"" + name + "\" execution time";
    StartTimer(timer_name);
  }

  void EndSection(const std::string &name)
  {
    const std::string timer_name = "Section \"" + name + "\" execution time";
    StopTimer(timer_name);

    ActiveProfMap::iterator active_item = active_sections.find(name);
    if (active_item != active_sections.end())
      active_item->second->second.StopEverything();

    active_sections.erase(name);
  }

  void Reset()
  {
    active_sections.clear();
    sections.clear();
  }

  void StartTimer(const std::string &name)
  {
    for(ActiveProfMap::iterator it = active_sections.begin(); it != active_sections.end(); ++it)
      it->second->second.StartTimer(name);
  }

  void StopTimer(const std::string &name)
  {
    for(ActiveProfMap::iterator it = active_sections.begin(); it != active_sections.end(); ++it)
      it->second->second.StopTimer(name);
  }

  void IncrementCounter(const std::string &name)
  {
    for(ActiveProfMap::iterator it = active_sections.begin(); it != active_sections.end(); ++it)
      it->second->second.IncrementCounter(name);
  }

  template<typename T> void SaveValue(const std::string &name, const T &value)
  {
    for(ActiveProfMap::iterator it = active_sections.begin(); it != active_sections.end(); ++it)
      it->second->second.SaveValue(name, value);
  }

  void PrintStatistics() const
  {
    std::cout << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "==============        PROFILER REPORT BELOW        ==============" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << std::endl;

    std::cout << "  Profiler \"" << m_name << "\":" << std::endl;

    if (sections.size() > 0)
    {
      for(ProfMap::const_iterator it = sections.begin(); it != sections.end(); ++it)
        it->second.PrintStatistics();
    }
    else
    {
      std::cout << "   No data collected" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "===============           END OF REPORT           ===============" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << "=================================================================" << std::endl;
    std::cout << std::endl;
  }

};

extern Profiler pGlobalProfiler;

#endif
