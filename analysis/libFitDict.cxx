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
   static TClass *RooMcbPdf_Dictionary();
   static void RooMcbPdf_TClassManip(TClass*);
   static void *new_RooMcbPdf(void *p = 0);
   static void *newArray_RooMcbPdf(Long_t size, void *p);
   static void delete_RooMcbPdf(void *p);
   static void deleteArray_RooMcbPdf(void *p);
   static void destruct_RooMcbPdf(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooMcbPdf*)
   {
      ::RooMcbPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RooMcbPdf));
      static ::ROOT::TGenericClassInfo 
         instance("RooMcbPdf", "RooMcbPdf.h", 19,
                  typeid(::RooMcbPdf), DefineBehavior(ptr, ptr),
                  &RooMcbPdf_Dictionary, isa_proxy, 0,
                  sizeof(::RooMcbPdf) );
      instance.SetNew(&new_RooMcbPdf);
      instance.SetNewArray(&newArray_RooMcbPdf);
      instance.SetDelete(&delete_RooMcbPdf);
      instance.SetDeleteArray(&deleteArray_RooMcbPdf);
      instance.SetDestructor(&destruct_RooMcbPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooMcbPdf*)
   {
      return GenerateInitInstanceLocal((::RooMcbPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RooMcbPdf*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RooMcbPdf_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RooMcbPdf*)0x0)->GetClass();
      RooMcbPdf_TClassManip(theClass);
   return theClass;
   }

   static void RooMcbPdf_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RootEvent_Dictionary();
   static void RootEvent_TClassManip(TClass*);
   static void *new_RootEvent(void *p = 0);
   static void *newArray_RootEvent(Long_t size, void *p);
   static void delete_RootEvent(void *p);
   static void deleteArray_RootEvent(void *p);
   static void destruct_RootEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RootEvent*)
   {
      ::RootEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RootEvent));
      static ::ROOT::TGenericClassInfo 
         instance("RootEvent", "RootEvent.h", 17,
                  typeid(::RootEvent), DefineBehavior(ptr, ptr),
                  &RootEvent_Dictionary, isa_proxy, 0,
                  sizeof(::RootEvent) );
      instance.SetNew(&new_RootEvent);
      instance.SetNewArray(&newArray_RootEvent);
      instance.SetDelete(&delete_RootEvent);
      instance.SetDeleteArray(&deleteArray_RootEvent);
      instance.SetDestructor(&destruct_RootEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RootEvent*)
   {
      return GenerateInitInstanceLocal((::RootEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RootEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RootEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RootEvent*)0x0)->GetClass();
      RootEvent_TClassManip(theClass);
   return theClass;
   }

   static void RootEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *RootMC_Dictionary();
   static void RootMC_TClassManip(TClass*);
   static void *new_RootMC(void *p = 0);
   static void *newArray_RootMC(Long_t size, void *p);
   static void delete_RootMC(void *p);
   static void deleteArray_RootMC(void *p);
   static void destruct_RootMC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RootMC*)
   {
      ::RootMC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::RootMC));
      static ::ROOT::TGenericClassInfo 
         instance("RootMC", "RootMC.h", 17,
                  typeid(::RootMC), DefineBehavior(ptr, ptr),
                  &RootMC_Dictionary, isa_proxy, 0,
                  sizeof(::RootMC) );
      instance.SetNew(&new_RootMC);
      instance.SetNewArray(&newArray_RootMC);
      instance.SetDelete(&delete_RootMC);
      instance.SetDeleteArray(&deleteArray_RootMC);
      instance.SetDestructor(&destruct_RootMC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RootMC*)
   {
      return GenerateInitInstanceLocal((::RootMC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RootMC*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *RootMC_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::RootMC*)0x0)->GetClass();
      RootMC_TClassManip(theClass);
   return theClass;
   }

   static void RootMC_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_Selector(void *p = 0);
   static void *newArray_Selector(Long_t size, void *p);
   static void delete_Selector(void *p);
   static void deleteArray_Selector(void *p);
   static void destruct_Selector(void *p);
   static void streamer_Selector(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Selector*)
   {
      ::Selector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Selector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Selector", ::Selector::Class_Version(), "Selector.h", 19,
                  typeid(::Selector), DefineBehavior(ptr, ptr),
                  &::Selector::Dictionary, isa_proxy, 16,
                  sizeof(::Selector) );
      instance.SetNew(&new_Selector);
      instance.SetNewArray(&newArray_Selector);
      instance.SetDelete(&delete_Selector);
      instance.SetDeleteArray(&deleteArray_Selector);
      instance.SetDestructor(&destruct_Selector);
      instance.SetStreamerFunc(&streamer_Selector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Selector*)
   {
      return GenerateInitInstanceLocal((::Selector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::Selector*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Selector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Selector::Class_Name()
{
   return "Selector";
}

//______________________________________________________________________________
const char *Selector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Selector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Selector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Selector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Selector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Selector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Selector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Selector*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RooMcbPdf(void *p) {
      return  p ? new(p) ::RooMcbPdf : new ::RooMcbPdf;
   }
   static void *newArray_RooMcbPdf(Long_t nElements, void *p) {
      return p ? new(p) ::RooMcbPdf[nElements] : new ::RooMcbPdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooMcbPdf(void *p) {
      delete ((::RooMcbPdf*)p);
   }
   static void deleteArray_RooMcbPdf(void *p) {
      delete [] ((::RooMcbPdf*)p);
   }
   static void destruct_RooMcbPdf(void *p) {
      typedef ::RooMcbPdf current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RooMcbPdf

namespace ROOT {
   // Wrappers around operator new
   static void *new_RootEvent(void *p) {
      return  p ? new(p) ::RootEvent : new ::RootEvent;
   }
   static void *newArray_RootEvent(Long_t nElements, void *p) {
      return p ? new(p) ::RootEvent[nElements] : new ::RootEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_RootEvent(void *p) {
      delete ((::RootEvent*)p);
   }
   static void deleteArray_RootEvent(void *p) {
      delete [] ((::RootEvent*)p);
   }
   static void destruct_RootEvent(void *p) {
      typedef ::RootEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RootEvent

namespace ROOT {
   // Wrappers around operator new
   static void *new_RootMC(void *p) {
      return  p ? new(p) ::RootMC : new ::RootMC;
   }
   static void *newArray_RootMC(Long_t nElements, void *p) {
      return p ? new(p) ::RootMC[nElements] : new ::RootMC[nElements];
   }
   // Wrapper around operator delete
   static void delete_RootMC(void *p) {
      delete ((::RootMC*)p);
   }
   static void deleteArray_RootMC(void *p) {
      delete [] ((::RootMC*)p);
   }
   static void destruct_RootMC(void *p) {
      typedef ::RootMC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RootMC

//______________________________________________________________________________
void Selector::Streamer(TBuffer &R__b)
{
   // Stream an object of class Selector.

   TSelector::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Selector(void *p) {
      return  p ? new(p) ::Selector : new ::Selector;
   }
   static void *newArray_Selector(Long_t nElements, void *p) {
      return p ? new(p) ::Selector[nElements] : new ::Selector[nElements];
   }
   // Wrapper around operator delete
   static void delete_Selector(void *p) {
      delete ((::Selector*)p);
   }
   static void deleteArray_Selector(void *p) {
      delete [] ((::Selector*)p);
   }
   static void destruct_Selector(void *p) {
      typedef ::Selector current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_Selector(TBuffer &buf, void *obj) {
      ((::Selector*)obj)->::Selector::Streamer(buf);
   }
} // end of namespace ROOT for class ::Selector

namespace {
  void TriggerDictionaryInitialization_libFitDict_Impl() {
    static const char* headers[] = {
"libFit.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/root-6.04.10/include/root",
"/home/nikolaev/work/BES/JpsiKK/analysis/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$libFit.h")))  RooMcbPdf;
class __attribute__((annotate("$clingAutoload$libFit.h")))  RootEvent;
class __attribute__((annotate("$clingAutoload$libFit.h")))  RootMC;
class __attribute__((annotate("$clingAutoload$libFit.h")))  Selector;
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
"RooMcbPdf", payloadCode, "@",
"RootEvent", payloadCode, "@",
"RootMC", payloadCode, "@",
"Selector", payloadCode, "@",
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
