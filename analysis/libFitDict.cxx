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
