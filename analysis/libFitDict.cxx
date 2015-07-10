// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME libFitDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "libFit.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_ModifiedCrystalBall(void *p = 0);
   static void *newArray_ModifiedCrystalBall(Long_t size, void *p);
   static void delete_ModifiedCrystalBall(void *p);
   static void deleteArray_ModifiedCrystalBall(void *p);
   static void destruct_ModifiedCrystalBall(void *p);
   static void streamer_ModifiedCrystalBall(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ModifiedCrystalBall*)
   {
      ::ModifiedCrystalBall *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ModifiedCrystalBall >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ModifiedCrystalBall", ::ModifiedCrystalBall::Class_Version(), "ModifiedCrystalBall.h", 20,
                  typeid(::ModifiedCrystalBall), DefineBehavior(ptr, ptr),
                  &::ModifiedCrystalBall::Dictionary, isa_proxy, 16,
                  sizeof(::ModifiedCrystalBall) );
      instance.SetNew(&new_ModifiedCrystalBall);
      instance.SetNewArray(&newArray_ModifiedCrystalBall);
      instance.SetDelete(&delete_ModifiedCrystalBall);
      instance.SetDeleteArray(&deleteArray_ModifiedCrystalBall);
      instance.SetDestructor(&destruct_ModifiedCrystalBall);
      instance.SetStreamerFunc(&streamer_ModifiedCrystalBall);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ModifiedCrystalBall*)
   {
      return GenerateInitInstanceLocal((::ModifiedCrystalBall*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::ModifiedCrystalBall*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ModifiedCrystalBall::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ModifiedCrystalBall::Class_Name()
{
   return "ModifiedCrystalBall";
}

//______________________________________________________________________________
const char *ModifiedCrystalBall::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModifiedCrystalBall*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ModifiedCrystalBall::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModifiedCrystalBall*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ModifiedCrystalBall::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModifiedCrystalBall*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ModifiedCrystalBall::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModifiedCrystalBall*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ModifiedCrystalBall::Streamer(TBuffer &R__b)
{
   // Stream an object of class ModifiedCrystalBall.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      fX.Streamer(R__b);
      fSigma.Streamer(R__b);
      {
         vector<RooRealProxy> &R__stl =  fStaple;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            RooRealProxy R__t;
            R__t.Streamer(R__b);
            R__stl.push_back(R__t);
         }
      }
      {
         vector<RooRealProxy> &R__stl =  fN;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            RooRealProxy R__t;
            R__t.Streamer(R__b);
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, ModifiedCrystalBall::IsA());
   } else {
      R__c = R__b.WriteVersion(ModifiedCrystalBall::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      fX.Streamer(R__b);
      fSigma.Streamer(R__b);
      {
         vector<RooRealProxy> &R__stl =  fStaple;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<RooRealProxy>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            ((RooRealProxy&)(*R__k)).Streamer(R__b);
            }
         }
      }
      {
         vector<RooRealProxy> &R__stl =  fN;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<RooRealProxy>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            ((RooRealProxy&)(*R__k)).Streamer(R__b);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ModifiedCrystalBall(void *p) {
      return  p ? new(p) ::ModifiedCrystalBall : new ::ModifiedCrystalBall;
   }
   static void *newArray_ModifiedCrystalBall(Long_t nElements, void *p) {
      return p ? new(p) ::ModifiedCrystalBall[nElements] : new ::ModifiedCrystalBall[nElements];
   }
   // Wrapper around operator delete
   static void delete_ModifiedCrystalBall(void *p) {
      delete ((::ModifiedCrystalBall*)p);
   }
   static void deleteArray_ModifiedCrystalBall(void *p) {
      delete [] ((::ModifiedCrystalBall*)p);
   }
   static void destruct_ModifiedCrystalBall(void *p) {
      typedef ::ModifiedCrystalBall current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_ModifiedCrystalBall(TBuffer &buf, void *obj) {
      ((::ModifiedCrystalBall*)obj)->::ModifiedCrystalBall::Streamer(buf);
   }
} // end of namespace ROOT for class ::ModifiedCrystalBall

namespace ROOT {
   static TClass *vectorlERooRealProxygR_Dictionary();
   static void vectorlERooRealProxygR_TClassManip(TClass*);
   static void *new_vectorlERooRealProxygR(void *p = 0);
   static void *newArray_vectorlERooRealProxygR(Long_t size, void *p);
   static void delete_vectorlERooRealProxygR(void *p);
   static void deleteArray_vectorlERooRealProxygR(void *p);
   static void destruct_vectorlERooRealProxygR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RooRealProxy>*)
   {
      vector<RooRealProxy> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RooRealProxy>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RooRealProxy>", -2, "vector", 214,
                  typeid(vector<RooRealProxy>), DefineBehavior(ptr, ptr),
                  &vectorlERooRealProxygR_Dictionary, isa_proxy, 0,
                  sizeof(vector<RooRealProxy>) );
      instance.SetNew(&new_vectorlERooRealProxygR);
      instance.SetNewArray(&newArray_vectorlERooRealProxygR);
      instance.SetDelete(&delete_vectorlERooRealProxygR);
      instance.SetDeleteArray(&deleteArray_vectorlERooRealProxygR);
      instance.SetDestructor(&destruct_vectorlERooRealProxygR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RooRealProxy> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<RooRealProxy>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERooRealProxygR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RooRealProxy>*)0x0)->GetClass();
      vectorlERooRealProxygR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERooRealProxygR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERooRealProxygR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<RooRealProxy> : new vector<RooRealProxy>;
   }
   static void *newArray_vectorlERooRealProxygR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<RooRealProxy>[nElements] : new vector<RooRealProxy>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERooRealProxygR(void *p) {
      delete ((vector<RooRealProxy>*)p);
   }
   static void deleteArray_vectorlERooRealProxygR(void *p) {
      delete [] ((vector<RooRealProxy>*)p);
   }
   static void destruct_vectorlERooRealProxygR(void *p) {
      typedef vector<RooRealProxy> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RooRealProxy>

namespace {
  void TriggerDictionaryInitialization_libFitDict_Impl() {
    static const char* headers[] = {
"libFit.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/root-6.04/include/root",
"/home/nikolaev/work/BES/JpsiKK/analysis/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Your description goes here...)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$libFit.h")))  ModifiedCrystalBall;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "libFit.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"ModifiedCrystalBall", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libFitDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libFitDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libFitDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libFitDict() {
  TriggerDictionaryInitialization_libFitDict_Impl();
}
