#ifndef IBN_REGEXP_H
#define IBN_REGEXP_H

/*
 * =====================================================================================
 * 
 *        Filename:  Regexp.h
 * 
 *     Description:  Моя библиотека для работы с регулярными выражениями.
 * 
 *         Version:  1.0
 *         Created:  10.03.2004 17:53:02 NOVT
 *        Revision:  none
 *        Compiler:  gcc
 * 
 *          Author:  Ivan Nikolaev ()
 *         Company:  BINP
 *           Email:  I.B.Nikolaev@inp.nsk.su
 * 
 * =====================================================================================
 */


#include <sys/types.h>
#include <stdlib.h>

#include <regex.h>
#include <string>
#include <iostream>


namespace ibn
{
  class regex 
  {
    public:
      regex(const char * rg,int cflag=0);   
      ~regex(void)	{ /*regfree(preg);*/ }
      inline bool ismatch(const char * buf, int eflags);
      inline bool ismatch(std::string str);
    private:
      regex_t preg;
      size_t nmatch;
  };


  inline regex::regex(const char * rg, int cflag)
  {
    /*
       cflags may be the bitwise-or of one or more of the following:

       REG_EXTENDED 
       Use POSIX Extended Regular Expression syntax when interpreting regex.  
       If not set, POSIX Basic Regular Expression syntax is used.

       REG_ICASE
       Do not differentiate case.  Subsequent regexec searches using this pattern
       buffer will be case insensitive.

       REG_NOSUB
       Support  for  substring  addressing  of  matches  is not required.  The 
       nmatch and pmatch parameters to regexec are ignored if the pattern
       buffer supplied was compiled with this flag set.

       REG_NEWLINE
       Match-any-character operators don't match a newline.

       A non-matching list ([^...])  not containing a newline does not match a newline.

       Match-beginning-of-line operator (^) matches the empty string immediately after a newline, regardless of  whether  eflags,  the  execution
       flags of regexec, contains REG_NOTBOL.

       Match-end-of-line operator ($) matches the empty string immediately before a newline, regardless of whether eflags contains REG_NOTEOL.
       */
    int rv = regcomp(&preg,rg,cflag);
    if (  rv )
    {
      std::cerr << "Error compiling regular expression " << rg << "\n";
      exit(1);
    }
  }

  inline bool regex::ismatch(const char * buf, int eflag=0)
  {
    /*
       eflags may be the bitwise-or of one or both of REG_NOTBOL and  REG_NOTEOL  which  cause  changes  in
       matching behaviour described below.

       REG_NOTBOL
       The  match-beginning-of-line  operator  always  fails to match (but see the compilation flag REG_NEWLINE above) This flag may be used when
       different portions of a string are passed to regexec and the beginning of the string should not be interpreted as  the  beginning  of  the
       line.

       REG_NOTEOL
       The match-end-of-line operator always fails to match (but see the compilation flag REG_NEWLINE above) 
       */ 
    regmatch_t pmatch;
    if ( !regexec ( &preg, buf, 1, &pmatch,eflag ) ) return true;
    return false;
  }

}
#endif
