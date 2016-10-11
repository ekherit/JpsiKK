/*
 * =====================================================================================
 *
 *       Filename:  hashlist.h
 *
 *    Description:  This is hashlists of usefull channels
 *
 *        Version:  1.0
 *        Created:  07.10.2016 16:43:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once
#include <string>
#include <vector>
#include <list>

#include <TCut.h>

struct channel_view_t
{
   long n, count;
   std::string fstate;
   std::string title;
   std::string hash_string;
};


template < class Container   >
inline TCut make_cut( const Container  & ChanViewList)
{
  TCut hash_cut = "";
  for(const auto &  v : ChanViewList )
  {
      std::list<std::string> hash_string_list;
      boost::split(hash_string_list, v.hash_string, boost::is_any_of(", "));
      for(auto & s : hash_string_list)
      {
        hash_cut = hash_cut || TCut(("hash==0x"+s).c_str());
      }
  }
  return hash_cut;
}; 


const std::list<channel_view_t> ChanKK =
{
  { 75,    3064,  "π+π-K+K-     ",  "Ψ(2S) → π+π-(J/Ψ → K+K-) ", "3031ea55,cdffabb3,ecdbe915"}
};

const std::list<channel_view_t> ChanUU =
{
  { 74,     947,  "μ+μ-π+π-     ",  "Ψ(2S) → π+π-(J/Ψ → μ+μ-) ", "1397e7e7,8ac60398,aca004d7,cf7de4a7,daa0af44"},
};

const std::list<channel_view_t> ChanNoJPsi =
{
  {  1   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc1(1P) → π-(-K*(892)0 → π+K-)K+)                          ",  "183789a,d0c6c2a9"},
  {  4   ,    1  , "ɣɣπ+π+π-K- ",    "Ψ(2S) → π+(K*(892)0 → (π0 → ɣɣ)(K0 → (-K0s → π+π-)))K-                ",  "202c4645,112c7a5f"},
  {  6   ,    1  , "ɣπ+π-π-K+  ",    "Ψ(2S) → ɣ(-χc2(1P) → π-(-K0 → (-K0s → π+π-))K+)                       ",  "31c0aeed,4ba80f38"},
  {  7   ,    1  , "π+π+π-K-   ",    "Ψ(2S) → K-(K2∗(1430)+ → π+(K0 → (-K0s → π+π-)))                       ",  "349d88c1,227b430"},
  {  8   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc0(1P) → π-(-K*(892)0 → π+K-)K+)                          ",  "4f9ddb0a,9ed86139"},
  { 11   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc1(1P) → (K*(892)0 → π-K+)(-K*(892)0 → π+K-))             ",  "6df8df1e,59b99e9a"},
  { 12   ,    1  , "ɣɣπ+π-K+K- ",    "Ψ(2S) → K-(K1(1270)+ → (ω(782) → (π0 → ɣɣ)π+π-)K+)                    ",  "7ecfea62,901f2b09"},
  { 13   ,    1  , "ɣɣπ+π+π-K- ",    "Ψ(2S) → K-(K1(1270)+ → (π0 → ɣɣ)(K0∗(1430)+ → π+(K0 → (-K0s → π+π-))))",  "82e597b2,552e1a87"},
  { 14   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc2(1P) → (K*(892)0 → π-K+)(-K*(892)0 → π+K-))             ",  "859ea6c1,b1dfe745"},
  { 16   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc0(1P) → (K*(892)0 → π-K+)(-K*(892)0 → π+K-))             ",  "9db952f9,a9f8137d"},
  { 18   ,    1  , "ɣɣπ+π-π-K+ ",    "Ψ(2S) → K+(K1(1270)- → (π0 → ɣɣ)(K*(892)- → π-(-K0 → (-K0s → π+π-)))) ",  "a7ee33a3,7025be96"},
  { 19   ,    1  , "ɣπ+π+π-π-  ",    "Ψ(2S) → ɣ(-χc0(1P) → π+π-(f2(1270) → π+π-))                           ",  "abf8ac4c,8460d425"},
  { 20   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc0(1P) → π+(K*(892)0 → π-K+)K-)                           ",  "b1401950,6005a363"},
  { 21   ,    1  , "ɣπ+π+π-K-  ",    "Ψ(2S) → ɣ(-χc2(1P) → π+(K0 → (-K0s → π+π-))K-)                        ",  "b359e7a3,c9314676"},
  { 23   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc2(1P) → π+π-K+K-)                                        ",  "c0aab2b1,37ae5c8d"},
  { 24   ,    1  , "ɣπ+π-π-K+  ",    "Ψ(2S) → ɣ(-χc1(1P) → K+(K*(892)- → π-(-K0 → (-K0s → π+π-))))          ",  "c9fcaf60,15a795e9"},
  { 26   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc0(1P) → π+π-K+K-)                                        ",  "ccf8f4be,3bfc1a82"},
  { 27   ,    1  , "ɣɣπ+π+π-K- ",    "Ψ(2S) → K-(K1(1270)+ → (π0 → ɣɣ)(K*(892)+ → π+(K0 → (-K0s → π+π-))))  ",  "ccfb3422,1b30b917"},
  { 28   ,    1  , "ɣπ+π+π-π-  ",    "Ψ(2S) → ɣ(-χc2(1P) → (ρ(770)0 → π+π-)π+π-)                            ",  "cfa021dd,e03859b4"},
  { 29   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc1(1P) → π+π-K+K-)                                        ",  "d6d8305c,21dcde60"},
  { 30   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc2(1P) → π-(-K*(892)0 → π+K-)K+)                          ",  "df2e94dd,e6b2eee"},
  { 31   ,    1  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(-χc0(1P) → K+(K1(1270)- → π-(-K0*(1430)0 → π+K-)))          ",  "e4b16026,1b59b91"},
  { 32   ,    1  , "ɣπ+π+π-π-  ",    "Ψ(2S) → ɣπ+π+π-π-                                                     ",  "e6f0a031"},
  { 33   ,    1  , "ɣπ+π-π-K+  ",    "Ψ(2S) → ɣπ-(-K0s → π+π-)K+                                            ",  "e7aec6e2,861f5eda"},
  { 36   ,    1  , "ɣɣπ+π+π-π- ",    "Ψ(2S) → (π0 → ɣɣ)π+π+π-π-                                             ",  "f46414d5"},
  { 39   ,    2  , "ɣπ+π-K+K-  ",    "Ψ(2S) → ɣ(K*(892)0 → π-K+)(-K*(892)0 → π+K-)                          ",  "35c7a381"},
  { 40   ,    2  , "ɣπ+π+π-π-  ",    "Ψ(2S) → ɣ(-χc0(1P) → (ρ(770)0 → π+π-)π+π-)                            ",  "5f136e0a,708b1663"},
  { 41   ,    2  , "ɣπ+π-K+K-  ",    "Ψ(2S) → K+(K1(1270)- → (ρ(770)0 → ɣπ+π-)K-)                           ",  "6f6bf2a3,2eb71455"},
  { 43   ,    2  , "ɣπ+π+π-π-  ",    "Ψ(2S) → ɣ(-χc2(1P) → π+π+π-π-)                                        ",  "a8b93aa8,5fbdd494"},
  { 44   ,    2  , "π+π-K+K-   ",    "Ψ(2S) → π+(K2*(1430)0 → π-K+)K-                                       ",  "bcdd2654,4200e40e"},
  { 48   ,    3  , "ɣπ+π-π-K+  ",    "Ψ(2S) → ɣ(-χc1(1P) → π-(K0 → (-K0s → π+π-))K+)                        ",  "2ea2390e,54ca98db"},
  { 55   ,    4  , "ɣπ+π+π-π-  ",    "Ψ(2S) → ɣ(-χc1(1P) → (ρ(770)0 → π+π-)π+π-)                            ",  "110dcd9a,3e95b5f3"},
  { 57   ,    4  , "π+π-K+K-   ",    "Ψ(2S) → (ρ(770)0 → π+π-)K+K-                                          ",  "a5f18375"},
  { 59   ,    5  , "ɣπ+π+π-π-  ",    "Ψ(2S) → ɣ(-χc0(1P) → π+π+π-π-)                                        ",  "a4eb7ca7,53ef929b"},
  { 61   ,    6  , "ɣɣπ+π+π-π- ",    "Ψ(2S) → π+(b1(1235)- → π-(ω(782) → (π0 → ɣɣ)π+π-))                    ",  "7a564e5a,fac211ae"},
  { 63   ,    8  , "ɣπ+π+π-π-  ",    "Ψ(2S) → ɣ(-χc0(1P) → (ρ(770)0 → π+π-)(ρ(770)0 → π+π-))                ",  "27940cb5,13d54d31"},
  { 64   ,   15  , "π+π-K+K-   ",    "Ψ(2S) → (-K*(892)0 → π+K-)(K2*(1430)0 → π-K+)                         ",  "32d836e2,e39d8cd1"},
  { 65   ,   19  , "π+π-K+K-   ",    "Ψ(2S) → π+π-K+K-                                                      ",  "3f291a5a"},
  { 66   ,   21  , "π+π-K+K-   ",    "Ψ(2S) → K-(K1(1270)+ → (ρ(770)0 → π+π-)K+)                            ",  "552c11f1,d7b558bf"},
  { 67   ,   26  , "π+π-K+K-   ",    "Ψ(2S) → (K*(892)0 → π-K+)(-K*(892)0 → π+K-)                           ",  "e01958bf"},
  { 69   ,   57  , "π+π-K+K-   ",    "Ψ(2S) → K-(K1(1270)+ → π+π-K+)                                        ",  "329ac69e,f922962d"},
  { 70   ,   93  , "π+π-K+K-   ",    "Ψ(2S) → π-(-K*(892)0 → π+K-)K+                                        ",  "e8092c85,16d4eedf,e8092c85"},
  { 71   ,  116  , "π+π-K+K-   ",    "Ψ(2S) → K-(K1(1270)+ → π+(K*(892)0 → π-K+))                           ",  "163b9509,603b3e9a,c77e2f3a"},
  { 73   ,  255  , "π+π-K+K-   ",    "Ψ(2S) → K-(K1(1270)+ → π+(K0*(1430)0 → π-K+))                         ",  "4854d053,7c5b6e7b,99116a60"}
};

const std::list<channel_view_t> ChanJPsi =
{
  {  2,       1,  "ɣɣɣɣπ+π+π-K- ",  "Ψ(2S) → π+π-(J/Ψ → K-(K*(892)+ → π+(K0 → (-K0s → (π0 → ɣɣ)(π0 → ɣɣ)))))", "152f1c20,d714a2bd"},
  {  3,       1,  "ɣɣπ+π+π-π-   ",  "Ψ(2S) → π+π-(J/Ψ → (π0 → ɣɣ)(ρ(770)0 → π+π-))                          ", "167c1bc8"},
  {  5,       1,  "π+π+π-π-π-K+ ",  "Ψ(2S) → π+π-(J/Ψ → π-(-K0 → (-K0s → π+π-))K+)                          ", "270c3648,6925d619"},
  {  9,       1,  "-K0lπ+π-π-K+ ",  "Ψ(2S) → π+π-(J/Ψ → π-(-K0 → -K0l)K+)                                   ", "5721f99e,6cceb77a"},
  { 10,       1,  "μ+μ-ɣɣ       ",  "Ψ(2S) → ɣ(-χc2(1P) → ɣ(J/Ψ → μ+μ-))                                    ", "5b016575,ac058b49"},
  { 15,       1,  "ɣπ+π-K+K-    ",  "Ψ(2S) → π+π-(J/Ψ → ɣ(f0(1710) → K+K-))                                 ", "8c471f87"},
  { 17,       1,  "π+π+π+π-π-K- ",  "Ψ(2S) → π+π-(J/Ψ → π+(K0 → (-K0s → π+π-))K-)                           ", "a5957f06,ebbc9f57"},
  { 22,       1,  "μ+μ-ɣπ+π-    ",  "Ψ(2S) → (η → ɣπ+π-)(J/Ψ → μ+μ-)                                        ", "b6e90d8a"},
  { 25,       1,  "e+e-ɣɣɣπ+π-  ",  "Ψ(2S) → (η → e+e-ɣ)(J/Ψ → π-(ρ(770)+ → (π0 → ɣɣ)π+))                   ", "cc5b2bea,4ec262a4"},
  { 34,       1,  "π+π-n-n      ",  "Ψ(2S) → π+π-(J/Ψ → n-n)                                                ", "e990b7b9"},
  { 35,       1,  "e+e-ɣπ+π+π-π-",  "Ψ(2S) → π+π-(J/Ψ → π-(ρ(770)+ → (π0 → e+e-ɣ)π+))                       ", "ecab451c,ad77a3ea"},
  { 37,       1,  "π+π+π+π-π-π- ",  "Ψ(2S) → π+π-(J/Ψ → π+π+π-π-)                                           ", "f5bb9d48"},
  { 38,       2,  "-K0lπ+π-π-K+ ",  "Ψ(2S) → π+π-(J/Ψ → (-K0 → -K0l)(K*(892)0 → π-K+))                      ", "2555a9e3,7181a332"},
  { 42,       2,  "ɣɣK+K-       ",  "Ψ(2S) → ɣ(-χc2(1P) → ɣ(J/Ψ → K+K-))                                    ", "78a768c7,8fa386fb"},
  { 45,       2,  "π+π+π+π-π-K- ",  "Ψ(2S) → π+π-(J/Ψ → K-(K*(892)+ → π+(K0 → (-K0s → π+π-))))              ", "f3b62be6,c50c1717"},
  { 46,       2,  "-K0lπ+π+π-K- ",  "Ψ(2S) → π+π-(J/Ψ → K-(K*(892)+ → π+(K0 → -K0l)))                       ", "f3b81e60,e36bbd16"},
  { 47,       3,  "ɣπ+π+π-π-    ",  "Ψ(2S) → π+π-(J/Ψ → ɣ(f2(1270) → π+π-))                                 ", "19f87eb1"},
  { 49,       3,  "-K0lπ+π-π-K+ ",  "Ψ(2S) → π+π-(J/Ψ → K+(K*(892)- → π-(-K0 → -K0l)))                      ", "57600321,47b3a057"},
  { 50,       3,  "ɣɣπ+π-K+K-   ",  "Ψ(2S) → π+π-(J/Ψ → K-(K*(892)+ → (π0 → ɣɣ)K+))                         ", "5d977b97,df0e32d9"},
  { 51,       3,  "π+π+π-π-     ",  "Ψ(2S) → π+π-(J/Ψ → π+π-)                                               ", "6d47f0a5,9089b143"},
  { 52,       3,  "e+e-π+π-     ",  "Ψ(2S) → π+π-(J/Ψ → e+e-)                                               ", "acaee554"},
  { 53,       3,  "ɣπ+π+π-π-    ",  "Ψ(2S) → π+π-(J/Ψ → ɣ(f0(1710) → π+π-))                                 ", "d1310577"},
  { 54,       3,  "ɣπ+π+π-π-    ",  "Ψ(2S) → π+π-(J/Ψ → ɣ(f2(2150) → π+π-))                                 ", "e4ebc641"},
  { 56,       4,  "ɣɣπ+π-K+K-   ",  "Ψ(2S) → π+π-(J/Ψ → (π0 → ɣɣ)K+K-)                                      ", "84fce654"},
  { 58,       4,  "-K0lπ+π+π-K- ",  "Ψ(2S) → π+π-(J/Ψ → π+(K0 → -K0l)K-)                                    ", "d4c719cb,ef28572f"},
  { 60,       6,  "π+π+π-π-π-K+ ",  "Ψ(2S) → π+π-(J/Ψ → K+(K*(892)- → π-(-K0 → (-K0s → π+π-))))             ", "79d29da3,4f68a152"},
  { 62,       6,  "ɣɣπ+π+π-π-   ",  "Ψ(2S) → π+π-(J/Ψ → (π0 → ɣɣ)π+π-)                                      ", "d98afca4"},
  { 68,      34,  "ɣπ+π+π-π-    ",  "Ψ(2S) → π+π-(J/Ψ → ɣ(f4(2050) → π+π-))                                 ", "1ef1d66d"},
  { 72,     228,  "ɣɣπ+π+π-π-   ",  "Ψ(2S) → π+π-(J/Ψ → π+(ρ(770)- → (π0 → ɣɣ)π-))                          ", "cfe7a549,4d7eec07,59476175,8a4d7453,4d7eec07,cfe7a549,dbde283b"},
};


inline TCut make_cut(std::string hash_list_str)
{
  TCut cut = "";
  std::list<std::string> hash_list;
  boost::split(hash_list,hash_list_str,boost::is_any_of(", "));
  std::cout << " Initial hash_list: ";
  for ( auto & x : hash_list)
  {
    std::cout << x << ",";
  }
  std::cout << std::endl;
  for( auto & x : hash_list)
  {
    std::string s=x;
    if(x == "jpsi")  { cut = cut || make_cut (ChanJPsi);   continue; };
    if(x == "nojpsi") { cut = cut || make_cut (ChanNoJPsi); continue; };
    if(x == "KK") { cut = cut || make_cut (ChanKK); continue; };
    if(x == "uu" || x == "UU") { cut = cut || make_cut (ChanUU); continue; };
    cut = cut || TCut(("hash==0x"+s).c_str());
  }
  return cut;
};

  std::string jpsi = "152f1c20,d714a2bd,167c1bc8,270c3648,6925d619,5721f99e,6cceb77a,8c471f87,a5957f06,ebbc9f57,e990b7b9,ecab451c,ad77a3ea,f5bb9d48,2555a9e3,7181a332,f3b62be6,c50c1717,f3b81e60,e36bbd16,19f87eb1,57600321,47b3a057,5d977b97,df0e32d9,6d47f0a5,9089b143,acaee554,d1310577,e4ebc641,84fce654,d4c719cb,ef28572f,79d29da3,4f68a152,d98afca4,1ef1d66d,cfe7a549,4d7eec07,59476175,8a4d7453,4d7eec07,cfe7a549,dbde283b,1397e7e7,8ac60398,aca004d7,cf7de4a7,daa0af44,3031ea55,cdffabb3,ecdbe915";

  
  std::string bg_from_jpsi = "cfe7a549,4d7eec07,59476175,8a4d7453,4d7eec07,cfe7a549,1ef1d66d,79d29da3,4f68a152,d98afca4,d1310577,84fce654,5d977b97,df0e32d9,57600321,47b3a057,e4ebc641,acaee554,6d47f0a5,9089b143,2555a9e3,7181a332,f5bb9d48,e990b7b9,78a768c7,8fa386fb,5721f99e,6cceb77a,270c3648,6925d619,19f87eb1,167c1bc8";

  std::string nojpsi_bg = "183789a,d0c6c2a9,35c7a381,5f136e0a,708b1663,6df8df1e,59b99e9a,6f6bf2a3,2eb71455,859ea6c1,b1dfe745,a7ee33a3,7025be96,a8b93aa8,5fbdd494,bcdd2654,4200e40e,ccf8f4be,3bfc1a82,cfa021dd,e03859b4,d6d8305c,21dcde60,e6f0a031,e7aec6e2,861f5eda,f46414d5,110dcd9a,3e95b5f3,a5f18375,2ea2390e,54ca98db,7a564e5a,fac211ae,27940cb5,13d54d31,32d836e2,e39d8cd1,552c11f1,d7b558bf,3f291a5a,e01958bf,329ac69e,f922962d,16d4eedf,e8092c85,163b9509,c77e2f3a,4854d053,7c5b6e7b,99116a60";
    
  std::string uu_str = "1397e7e7,8ac60398,aca004d7,cf7de4a7,daa0af44";
  std::string KK_str = "3031ea55,cdffabb3,ecdbe915";
  std::string pipi_str = "67c1bc8,19f87eb1,270c3648,6925d619,5f136e0a,708b1663,a8b93aa8,5fbdd494,cfa021dd,e03859b4,e6f0a031,f46414d5,f5bb9d48,110dcd9a,3e95b5f3,6d47f0a5,9089b143,e4ebc641,7a564e5a,fac211ae,d1310577,d98afca4,27940cb5,13d54d31,79d29da3,4f68a152,1ef1d66d,cfe7a549,4d7eec07,59476175,8a4d7453,4d7eec07,cfe7a549";
  std::string pipiKK_str = "183789a,d0c6c2a9,35c7a381,6df8df1e,59b99e9a,6f6bf2a3,2eb71455,859ea6c1,b1dfe745,bcdd2654,4200e40e,ccf8f4be,3bfc1a82,d6d8305c,21dcde60,a5f18375,5d977b97,df0e32d9,84fce654,32d836e2,e39d8cd1,552c11f1,d7b558bf,3f291a5a,e01958bf,329ac69e,f922962d,16d4eedf,e8092c85,163b9509,c77e2f3a,4854d053,7c5b6e7b,99116a60";

  std::string jpsirho="cfe7a549,4d7eec07,59476175,8a4d7453,4d7eec07,cfe7a549,dbde283b";
