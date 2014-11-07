#ifndef IBN_AVERAGER_H
#define IBN_AVERAGER_H
#include <queue>
#include <stdexcept>
#include <algorithm>
#include <limits>

namespace ibn 	
{
  template <class T,class W=T>
    class averager
    {
      public:
        struct pair
        { 
          T data; 
          W weight; 
          pair(T d, W w) 
          { 
            data=d;
            weight=w; 
          }
          bool operator<(const pair & p) const { return data<p.data; }
          bool operator>(const pair & p) const { return data>p.data; }
          bool operator<=(const pair & p) const { return data<=p.data; }
          bool operator>=(const pair & p) const { return data>=p.data; }
        };
        typedef std::deque <pair> Data_t;
        typedef typename Data_t::iterator iterator;
        typedef typename Data_t::const_iterator const_iterator;
        typedef typename Data_t::reverse_iterator reverse_iterator;
        typedef typename Data_t::const_reverse_iterator const_reverse_iterator;
        Data_t array;

        averager(unsigned p=0) { resize(p); }
        inline void reset(void); 
        inline void resize(unsigned p=0) { max_size=p; reset(); }
        inline void add(const T &, W w=1);

        inline void pop_back(void); //delete last element
        inline void pop_front(void); //delete first element

        inline void pop(const T &, W w=1); //delete some element

        inline averager & operator+=(T data) { add(data); return *this;}

        inline T average(void)	{ return xsum/N; }
        inline T average2(void) { return x2sum/N; }
        inline T rms(void)	
        { 
          double r=average2()-sq(average());
          if(r<=0) r=-r;
          return sqrt(r); 
        }
        inline T sigma(void)	{ return rms()*sqrt(static_cast<T>(N)/(N-1)); }
        inline T sigma_average(void) { return rms()*sqrt(static_cast<T>(1)/(N-1)); }

        inline T waverage(void) { return xwsum/wsum; }
        inline T waverage2(void) { return x2wsum/wsum; }
        inline T chi2(void) { return x2wsum - sq(xwsum)/wsum; }
        inline unsigned ndf(void) { return N-1; }

        inline T wsigma(void)	{ return sqrt(chi2()/wsum); }
        inline T wsigma_average(void);
        inline T scale(void)	{ return chi2()/(N-1); }
        inline T min(void) const { return max_size==0 ? _min : min_element(array.begin(),array.end())->data;}
        inline T max(void) const { return max_size==0 ? _max : max_element(array.begin(),array.end())->data;}
        inline T front(void) const { return max_size==0 ? _front : array.front().data; }
        inline T back(void) const { return max_size==0 ? _back : array.back().data ;}

        inline const_iterator begin(void) const 
        { 
          if(max_size==0) throw std::runtime_error("averager array is empty for max_element=0");
          return  array.begin(); 
        }

        inline const_iterator end(void) const 
        { 
          if(max_size==0) throw std::runtime_error("averager array is empty for max_element=0");
          return  array.end(); 
        }

        inline const_reverse_iterator rbegin(void) const 
        { 
          if(max_size==0) throw std::runtime_error("averager array is empty for max_element=0");
          return  array.rbegin(); 
        }

        inline const_reverse_iterator rend(void) const 
        { 
          if(max_size==0) throw std::runtime_error("averager array is empty for max_element=0");
          return  array.rend(); 
        }


        inline size_t size(void) { return N; }
        inline W get_wsum(void) { return wsum; }

        inline void recalculate(void) 
        {
          xsum=0;
          x2sum=0;
          wsum=0;
          xwsum=0;
          x2wsum=0;
          for(iterator i =array.begin(); i!=array.end(); i++) 
          {
            xsum+=i->data;
            x2sum+=i->data*i->data;
            xwsum+=i->data*i->weight;
            x2wsum+=i->data*i->data*i->weight;
            wsum+=i->weight;
          }
        }
      private:
        size_t max_size;    //max allowed size
        size_t N;	    //current size
        W wsum; //weight accumulator
        W xwsum,x2wsum; //x,x2 accumulator
        T xsum,x2sum;
        T _min,_max; // minimum and maximum elementy in raw
        T _front,_back; //first and last element


        inline T sq( T x ) { return x*x; }
      public:
        inline Data_t & get_array(void) { return array; }
    };

  template <class T, class W>  
    inline void averager<T,W>::reset(void)
    { 
      wsum=0; 
      xwsum=0;
      x2wsum=0;
      N=0; 
      xsum=0; 
      x2sum=0; 
      //_min=numeric_limits<T>::max();
      //_max=-numeric_limits<T>::max();
      //_front=numeric_limits<T>::quiet_NaN();
      //_back=numeric_limits<T>::quiet_NaN();
      array.erase(array.begin(),array.end());
    }

  template <class T, class W>  
    inline void averager<T,W>::add( const T & data, W  w)	{
      if(w<=0) return; //Это ещё вопрос, нужно ли сохранять такие данные с нулевым весом.
      if( max_size !=0 )  {
        if( array.size() >= max_size) {
          //Если массив заполнен, то удалим самые первые данные
          //xwsum  -= (array.front().data*array.front().weight);
          //x2wsum -= (array.front().data*array.front().data*array.front().weight);
          //wsum  -=  array.front().weight;
          //
          //xsum -= array.front().data;
          //x2sum -= sq(array.front().data);
          //N--; 
          pop(array.back().data, array.back().weight);
          array.pop_front(); //само удаление
        }
        array.push_back(pair(data,w)); //запомним новые данные
      }
      //добавим новые данные в статистику
      xwsum+=(data*w);
      x2wsum+=(data*data*w);
      wsum+=w;
      //cout << xwsum << " " << x2sum << " " << wsum << endl;

      xsum+=data;
      x2sum+=(data*data);
      //fill max min
      if(data > _max) _max=data;
      if(data < _min) _min=data;
      //fill front back
      _back=data;
      if(N==0) _front=data;
      N++;
    }

  template <class T, class W>  
    inline void averager<T,W>::pop_back(void)	{
      if( max_size !=0 )  {
        if(N==0) return;
        pop(array.back().data, array.back().weight);
        array.pop_back(); //само удаление
      }
    }

  template <class T, class W>  
    inline void averager<T,W>::pop_front(void)	{
      if( max_size !=0 )  {
        //Remove first data
        if(N==0) return;
        pop(array.front().data, array.front().weight);
        array.pop_front(); 
      }
    }

  template <class T, class W>  
    inline void averager<T,W>::pop(const T & data, W weight)	{
      if(N==0) return;
      xwsum  -= (data*weight);
      x2wsum -= (data*data*weight);
      wsum  -=  weight;
      xsum -= data;
      x2sum -= sq(data);
      N--;
    }

  template <class T, class W> 
    inline T averager<T,W>::wsigma_average(void)	{ 
      if (N<=0)	{ 
      }
      if(N==1 || scale() <=1.0 )  return 1./sqrt(wsum);
      return sqrt(chi2()/wsum/(N-1)); 
    }

}

#endif
