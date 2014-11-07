/*
 * =====================================================================================
 *
 *       Filename:  valer.h
 *
 *    Description:  Value and error pair
 *
 *        Version:  1.0
 *        Created:  07/12/2010 11:12:46 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physycs, Novosibirsk, Russia
 *
 * =====================================================================================
 */
#ifndef IBN_VALER_H
#define IBN_VALER_H

#include <cmath>

namespace ibn
{
  template < class T>
  struct valer
  { 
    T value,error; //physical value is the value and the error
    valer(void) {value=0; error=0;} // I dont know should i init valer or not
    valer(const T & v) : value(v), error(0){}
    valer(const T & v, const T & e): value(v), error(e) {}

    struct reference
    {
      T & value, & error;
      inline reference(T & v,T & e) : value(v), error(e) {} 
      inline reference(valer <T> & v): value(v.value), error(v.error) {}
      inline void operator=(valer <T> & pv) { value=pv.value; error = pv.error; }
      inline void operator=(const valer <T> & pv) { value=pv.value; error = pv.error; }

      inline reference & operator+=(const T & t) { value+=t; return *this;}
      inline reference & operator-=(const T & t) { value-=t; return *this;}
      inline reference & operator*=(const T & t) { value*=t; error*=t; return *this;}
      inline reference & operator/=(const T & t) { value/=t; error/=t; return *this;}

      inline reference & operator+=(const valer<T> & v) { *this=(valer<T>(*this)+=v); return *this;}
      inline reference & operator-=(const valer<T> & v) { *this=(valer<T>(*this)-=v); return *this;}
      inline reference & operator*=(const valer<T> & v) { *this=(valer<T>(*this)*=v); return *this;}
      inline reference & operator/=(const valer<T> & v) { *this=(valer<T>(*this)/=v); return *this;}

      inline valer<T> operator + (const T & y) { return valer<T>(*this)+=y; } 
      inline valer<T> operator + (const reference & y) { return valer<T>(*this)+=valer<T>(y); } 
      inline valer<T> operator * (const T & y) { return valer<T>(*this)*=y; } 
      inline valer<T> operator * (const reference & y) { return valer<T>(*this)*=valer<T>(y); } 
      inline valer<T> operator / (const T & y) { return valer<T>(*this)/=y; } 
      inline valer<T> operator / (const reference & y) { return valer<T>(*this)/=valer<T>(y); } 
    };

    struct const_reference
    {
      const T & value, & error;
      inline const_reference( const T & v, const T & e) : value(v), error(e) {} 
      inline const_reference( const valer <T> & pv): value(pv.value), error(pv.error) {}
      inline const_reference( const valer <T>::reference & pv): value(pv.value), error(pv.error) {}
      inline const_reference( const const_reference & pv): value(pv.value), error(pv.error) {}
    };

    inline valer(const valer<T> & v): value(v.value), error(v.error){}
    inline valer(const valer<T>::reference & v) : value(v.value), error(v.error) {}
    inline valer(const valer<T>::const_reference & v): value(v.value), error(v.error) {}

    inline void operator=(const T & v){value=v; error=0;}
    inline void operator=(const valer<T> & v) { value=v.value; error=v.error;}
    inline void operator=(const valer<T>::reference & v){ value=v.value; error=v.error;};
    inline void operator=(const valer<T>::const_reference & v){ value=v.value; error=v.error;};

    inline valer <T>& operator+=(const T & x)
    {
      value+=x;
      return *this;
    }

    inline valer <T>& operator+=(const valer <T> & x)
    {
      error=sqrt(error*error + x.error*x.error);
      value+=x.value;
      return *this;
    }

    
    inline valer <T>& operator-=(const T & x)
    {
      value-=x;
      return *this;
    }

    inline valer <T>& operator-=(const valer <T> & x)
    {
      error=sqrt(error*error + x.error*x.error);
      value-=x.value;
      return *this;
    }

    inline valer <T>& operator*=(const T & x)
    {
      value*=x;
      error*=x;
      return *this;
    }

    inline valer <T> & operator*=(const valer <T> & x)
    {
      error=sqrt(error*error*x.value*x.value + value*value*x.error*x.error);
      value*=x.value;
      return *this;
    }

    inline valer <T>& operator/=(const T & x)
    {
      value/=x;
      error/=x;
      return *this;
    }

    inline valer <T> & operator/=(const valer <T> & x)
    {
      T x2=x.value*x.value;
      error=sqrt(error*error/x2 + value*value/(x2*x2)*x.error*x.error);
      value/=x.value;
      return *this;
    }
    //inline valer<T>& operator+=(const valer<T>::reference & x)
    //{
    //  return *this+=valer<T>(x);
    //}

    //inline valer<T> & operator+=(const valer<T>::const_reference & x)
    //{
    //  return *this+=valer<T>(x);
    //}


    //inline valer <T>& operator-=(const valer<T>::const_reference & x)
    //{
    //  return *this+=valer<T>(x);
    //}

    //inline valer <T> operator+(const valer <T> & x)
    //{
    //  valer <T> result(*this);
    //  result+=x;
    //  return result;
    //}

    //inline valer <T> operator+(const T & x)
    //{
    //  return *this+valer<T>(x);
    //}

    //inline valer <T> operator+(const valer_reference<T> & x)
    //{
    //  return *this+valer<T>(x);
    //}

    //inline valer <T>  operator+(const valer<T>::const_reference & x)
    //{
    //  return *this+valer<T>(x);
    //}

    //inline valer <T> operator-(const valer <T> & x)
    //{
    //  valer <T> result(*this);
    //  result-=x;
    //  return result;
    //}

    //inline valer <T>  operator-(const T & x)
    //{
    //  return *this-valer<T>(x);
    //}


    //inline valer <T>  operator-(const valer<T>::const_reference & x)
    //{
    //  valer <T> result(*this);
    //  result-=x;
    //  return result;
    //  return *this-valer<T>(x);
    //}
    



    //inline valer<T> & operator*=(const valer<T>::const_reference & x)
    //{
    //  return *this*=valer<T>(x);
    //}

    //inline valer <T> operator*(const valer <T> & x)
    //{
    //  valer <T> result(*this);
    //  result*=x;
    //  return result;
    //}

    //inline valer <T>  operator*(const T & x)
    //{
    //  return *this*valer<T>(x);
    //}


    //inline valer <T>  operator*(const valer<T>::const_reference & x)
    //{
    //  valer <T> result(*this);
    //  result*=x;
    //  return result;
    //  //return *this*valer<T>(x);
    //}




    //inline valer<T> & operator/=(const valer<T>::const_reference & x)
    //{
    //  return *this/=valer<T>(x);
    //}

    //inline valer <T> operator/(const valer <T> & x)
    //{
    //  valer <T> result(*this);
    //  result/=x;
    //  return result;
    //}

    //inline valer <T>  operator/(const T & x)
    //{
    //  return *this/valer<T>(x);
    //}


    //inline valer <T>  operator/(const valer<T>::const_reference & x)
    //{
    //  return *this/valer<T>(x);
    //}

    private:
  };


  template <class T> inline valer<T> operator + (const valer<T> & x, const valer<T> & y) { return valer<T>(x)+=y; } 
  template <class T> inline valer<T> operator - (const valer<T> & x, const valer<T> & y) { return valer<T>(x)-=y; }
  template <class T> inline valer<T> operator * (const valer<T> & x, const valer<T> & y) { return valer<T>(x)*=y; } 
  template <class T> inline valer<T> operator / (const valer<T> & x, const valer<T> & y) { return valer<T>(x)/=y; }

  template <class T> inline valer<T> operator + (const valer<T> & x, const T & y) { return valer<T>(x)+=y; } 
  template <class T> inline valer<T> operator - (const valer<T> & x, const T & y) { return valer<T>(x)-=y; }
  template <class T> inline valer<T> operator * (const valer<T> & x, const T & y) { return valer<T>(x)*=y; } 
  template <class T> inline valer<T> operator / (const valer<T> & x, const T & y) { return valer<T>(x)/=y; }

  template <class T> inline valer<T> operator + (const T & x, const valer<T> & y) { return valer<T>(x)+=y; } 
  template <class T> inline valer<T> operator - (const T & x, const valer<T> & y) { return valer<T>(x)-=y; }
  template <class T> inline valer<T> operator * (const T & x, const valer<T> & y) { return valer<T>(x)*=y; } 
  template <class T> inline valer<T> operator / (const T & x, const valer<T> & y) { return valer<T>(x)/=y; }


  // reference operations
  template <class T> inline valer<T> operator + (const valer<T> & x, const typename valer<T>::reference & y) { return valer<T>(x)+=y; } 
  template <class T> inline valer<T> operator - (const valer<T> & x, const typename valer<T>::reference & y) { return valer<T>(x)-=y; }
  template <class T> inline valer<T> operator * (const valer<T> & x, const typename valer<T>::reference & y) { return valer<T>(x)*=y; } 
  template <class T> inline valer<T> operator / (const valer<T> & x, const typename valer<T>::reference & y) { return valer<T>(x)/=y; }

  template <class T> inline valer<T> operator + (const typename valer<T>::reference & x, const valer<T> & y) { return valer<T>(x)+=y; } 
  template <class T> inline valer<T> operator - (const typename valer<T>::reference & x, const valer<T> & y) { return valer<T>(x)-=y; }
  template <class T> inline valer<T> operator * (const typename valer<T>::reference & x, const valer<T> & y) { return valer<T>(x)*=y; } 
  template <class T> inline valer<T> operator / (const typename valer<T>::reference & x, const valer<T> & y) { return valer<T>(x)/=y; }

  //template <class T> inline valer<T> operator + (const typename valer<T>::reference & x, const typename valer<T>::reference & y) { return valer<T>(x)+=y; } 
  //template <class T> inline valer<T> operator - (const typename valer<T>::reference & x, const typename valer<T>::reference & y) { return valer<T>(x)-=y; }
  //template <class T> inline valer<T> operator * (const typename valer<T>::reference & x, const typename valer<T>::reference & y) { return valer<T>(x)*=y; } 
  template <class T> inline valer<T> operator / (typename valer<T>::reference  x, typename valer<T>::reference  y) { return valer<T>(x)/=y; }

  template <class T> inline bool operator == (const valer<T> & x, const valer<T> & y) { return x.value == y.value; } 
  template <class T> inline bool operator <  (const valer<T> & x, const valer<T> & y) { return x.value <  y.value; } 
  template <class T> inline bool operator <= (const valer<T> & x, const valer<T> & y) { return x.value <= y.value; } 
  template <class T> inline bool operator >  (const valer<T> & x, const valer<T> & y) { return x.value >  y.value; } 
  template <class T> inline bool operator >= (const valer<T> & x, const valer<T> & y) { return x.value >= y.value; } 

  template <class T> T sq ( const T & t ) { return t*t; }
  template <class T> T cb ( const T & t ) { return t*t*t; }

  template <class T>
  inline valer <T> pow(const valer <T>& x, const valer<T> &p)
  {
    if(x.value==0.0) return 0;
    if(p.value==0.0) return valer <T> (1.0, std::log(x.value)*p.error);
    T r = std::pow(x.value, p.value);
    return valer<T>( r,  sqrt( sq(x.error*p.value*r/x.value) + sq(p.error*std::log(x.value)*r)));
  }
}
#endif //IBN_VALER_H
