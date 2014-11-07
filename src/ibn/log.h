/*
 * =====================================================================================
 *
 *       Filename:  log.h
 *
 *    Description:  Control log output.
 *		    информации.
 *
 *        Version:  1.0
 *        Created:  02/08/2008 05:13:33 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */
#ifndef IBN_LOG_H
#define IBN_LOG_H

#include <iostream>

namespace ibn
{
  class log 
  {
    int LOG_LEVEL;
    std::ostream * log_stream; /*  output log stream */
    mutable std::ostream * os; /*  main ostream */
    mutable int message_rank; /*  rank of the current messages */
    std::ostream &(*head)(std::ostream&); /* head function  */
    mutable bool headed; /* whether head is printed or not */
    public:
    log (int level=0 /* log level */, std::ostream * =0/*  log stream */, std::ostream &(*h)(std::ostream&)=0 );
    ~log(void);
    void set_level(int); /* set log level to control output */
    void set_log(std::ostream *)  ; /*  set log stream */
    void set_head(std::ostream &(*h)(std::ostream&)); /*  set head function */
    void set_stream(std::ostream *)const; /*  set output stream */
    const log & operator()(int level) const; /*  set level for current message */
    const log &  operator<<(std::ostream & ) const; /*  connect with the output stream */
    const log &  operator<<(std::ios&(*pf)(std::ios &))const; /*  redirect ostream manipulators */
    const log &  operator<<(std::ios_base&(*pf)(std::ios_base &))const; /*  redirect ostream manipulators hex, etc... */
    const log &  operator<<(std::ostream&(*pf)(std::ostream &))const; /*  redirect manipulator such as endl */
    template<class T> const log & operator<<(const T  & message) const; /*  output messages */
    private:
    bool need_output(void) const; /*  check the output level  */
  };

  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  log
   * Description:  Costructor
   *--------------------------------------------------------------------------------------
   */
  inline log::log(int level /*  output level */, std::ostream * ls /*  log stream */, std::ostream &(*h)(std::ostream&) /*  head */) : 
    LOG_LEVEL(level),
    log_stream(ls), 
    head(h)
  {
    headed=false;
    os = 0;
  }

  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  ~log
   * Description:  Destructor. Nothing to destruct.
   *--------------------------------------------------------------------------------------
   */
  inline log::~log(void){}


  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  set_level
   * Description:  Set output log level. The message will be printed if rank
   * of message less or equal the log level.
   *--------------------------------------------------------------------------------------
   */
  inline void log::set_level(int l)
  {
    LOG_LEVEL =l;
  }


  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  set_log
   * Description:  Set log stream. Usualy this log stream should be the file.
   *--------------------------------------------------------------------------------------
   */
  inline void log::set_log(std::ostream * logstr)
  {
    log_stream = logstr;
  }


  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  set_head
   * Description: Set head function to output the head for every log line 
   *--------------------------------------------------------------------------------------
   */
  inline void log::set_head(std::ostream &(*h)(std::ostream&))
  {
    head = h;
  }


  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  set_stream
   * Description:  set current ouput stream
   *--------------------------------------------------------------------------------------
   */
  inline void log::set_stream(std::ostream * o) const
  {
    os = o;
  }


  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  operator()(int)
   * Description:  Init rank of the following messages.
   *--------------------------------------------------------------------------------------
   */
  inline const log & log::operator()(int l) const
  {
    message_rank=l;
    headed=false;
    set_stream(0);
    return *this;
  }

  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  operator<<(ostream&)
   * Description:  Connect the log with the output stream using commont C++ ostream
   * interface
   *--------------------------------------------------------------------------------------
   */
  inline const log & log::operator<<(std::ostream & o) const
  {
    set_stream(&o);/*  save the main ostream */
    return *this;
  }


  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  need_ouput
   * Description:  Check the rank of the messaged comparing with the output level
   * Return true if following messages are needed to pring.
   *--------------------------------------------------------------------------------------
   */
  inline bool log::need_output(void) const
  {
    /* For negative log level print the only messages which 
     * equal local log level */
    if(LOG_LEVEL<=0)
    {
      if(message_rank!=-LOG_LEVEL) 
        return false;
    }
    else 
    {
      /*  For positive log level  print only
       *  local log level less the global log level*/
      if(message_rank>LOG_LEVEL) 
        return false;
    }
    return true;
  }


  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  operator<<(const T &)
   * Description:  Print any messages if needed.
   *--------------------------------------------------------------------------------------
   */
  template<class T> 
  inline  const log & log::operator<<(const T  & message) const
  {
    if(!need_output()) return *this;
    if(head!=0 && !headed) //print head if head method is set and head had not been printed
    {
      if(os) *os << head;
      if(log_stream) *log_stream << head;
      headed = true;
    }
    if(os) *os << message;
    if(log_stream) *log_stream << message;
    return *this;
  }


  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  operator<<(ios_base manipulator)
   * Description:  Apply ios_base stream manipulators such as hex.
   *--------------------------------------------------------------------------------------
   */
  inline const log &  log::operator<<(std::ios_base&(*pf)(std::ios_base &)) const
  {
    if(!need_output()) return *this;
    if(log_stream) pf(*log_stream);
    if(os) pf(*os);
    return *this;
  }

  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  operator<<(ios manipulator)
   * Description:  Apply ios stream manipulators
   *--------------------------------------------------------------------------------------
   */
  inline const log &  log::operator<<(std::ios&(*pf)(std::ios &)) const
  {
    if(!need_output()) return *this;
    if(log_stream) pf(*log_stream);
    if(os) pf(*os);
    return *this;
  }

  /*
   *--------------------------------------------------------------------------------------
   *       Class:  log
   *      Method:  operator<<(ostream manipulator)
   * Description:  Apply ostream  manipulators such as endl
   *--------------------------------------------------------------------------------------
   */
  inline const log &  log::operator<<(std::ostream&(*pf)(std::ostream &)) const
  {
    if(!need_output()) return *this;
    if(log_stream) pf(*log_stream);
    if(os) pf(*os);
    return *this;
  }
}

namespace std
{

  /* 
   * ===  FUNCTION  ======================================================================
   *         Name:  operator<<
   *  Description:  Connect ostream with the log
   * =====================================================================================
   */
  inline const ibn::log & operator<<(std::ostream &os, const ibn::log & l)
  {
    l.set_stream(&os);
    return l;
  }
}
#endif
