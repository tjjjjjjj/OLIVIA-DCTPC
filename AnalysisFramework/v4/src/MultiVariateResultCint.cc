//
// File generated by rootcint at Mon Nov 23 19:12:55 2015

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME MultiVariateResultCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "MultiVariateResultCint.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void MultiVariatecLcLMultiVariateResult_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_MultiVariatecLcLMultiVariateResult(void *p = 0);
   static void *newArray_MultiVariatecLcLMultiVariateResult(Long_t size, void *p);
   static void delete_MultiVariatecLcLMultiVariateResult(void *p);
   static void deleteArray_MultiVariatecLcLMultiVariateResult(void *p);
   static void destruct_MultiVariatecLcLMultiVariateResult(void *p);
   static void streamer_MultiVariatecLcLMultiVariateResult(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MultiVariate::MultiVariateResult*)
   {
      ::MultiVariate::MultiVariateResult *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MultiVariate::MultiVariateResult >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MultiVariate::MultiVariateResult", ::MultiVariate::MultiVariateResult::Class_Version(), "./MultiVariate.hh", 83,
                  typeid(::MultiVariate::MultiVariateResult), DefineBehavior(ptr, ptr),
                  &::MultiVariate::MultiVariateResult::Dictionary, isa_proxy, 0,
                  sizeof(::MultiVariate::MultiVariateResult) );
      instance.SetNew(&new_MultiVariatecLcLMultiVariateResult);
      instance.SetNewArray(&newArray_MultiVariatecLcLMultiVariateResult);
      instance.SetDelete(&delete_MultiVariatecLcLMultiVariateResult);
      instance.SetDeleteArray(&deleteArray_MultiVariatecLcLMultiVariateResult);
      instance.SetDestructor(&destruct_MultiVariatecLcLMultiVariateResult);
      instance.SetStreamerFunc(&streamer_MultiVariatecLcLMultiVariateResult);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MultiVariate::MultiVariateResult*)
   {
      return GenerateInitInstanceLocal((::MultiVariate::MultiVariateResult*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MultiVariate::MultiVariateResult*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

      namespace MultiVariate {
//______________________________________________________________________________
TClass *MultiVariateResult::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *MultiVariateResult::Class_Name()
{
   return "MultiVariate::MultiVariateResult";
}

//______________________________________________________________________________
const char *MultiVariateResult::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MultiVariate::MultiVariateResult*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MultiVariateResult::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MultiVariate::MultiVariateResult*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void MultiVariateResult::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MultiVariate::MultiVariateResult*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *MultiVariateResult::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MultiVariate::MultiVariateResult*)0x0)->GetClass();
   return fgIsA;
}

} // namespace MultiVariate
      namespace MultiVariate {
//______________________________________________________________________________
void MultiVariateResult::Streamer(TBuffer &R__b)
{
   // Stream an object of class MultiVariate::MultiVariateResult.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::MultiVariate::MultiVariateResult thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> optimized;
      R__b >> reader;
      R__b >> E;
      R__b >> range;
      R__b >> skewness;
      R__b >> cluster_rms;
      R__b >> cluster_mean;
      R__b >> maxpixel;
      R__b >> maxpixel_over_E;
      R__b >> cluster_rms_over_E;
      R__b >> tmoment3;
      R__b >> tmoment4;
      R__b >> npixel;
      R__b >> npixel_red;
      R__b >> neighbors;
      R__b >> logE;
      R__b >> logrange;
      R__b >> logcluster_rms;
      R__b >> logmaxpixel;
      R__b >> lognpixel;
      R__b >> lognpixel_red;
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << optimized;
      R__b << reader;
      R__b << E;
      R__b << range;
      R__b << skewness;
      R__b << cluster_rms;
      R__b << cluster_mean;
      R__b << maxpixel;
      R__b << maxpixel_over_E;
      R__b << cluster_rms_over_E;
      R__b << tmoment3;
      R__b << tmoment4;
      R__b << npixel;
      R__b << npixel_red;
      R__b << neighbors;
      R__b << logE;
      R__b << logrange;
      R__b << logcluster_rms;
      R__b << logmaxpixel;
      R__b << lognpixel;
      R__b << lognpixel_red;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} // namespace MultiVariate
//______________________________________________________________________________
      namespace MultiVariate {
void MultiVariateResult::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class MultiVariate::MultiVariateResult.
      TClass *R__cl = ::MultiVariate::MultiVariateResult::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "optimized", &optimized);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*reader", &reader);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*current", &current);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "E", &E);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "range", &range);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "skewness", &skewness);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "cluster_rms", &cluster_rms);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "cluster_mean", &cluster_mean);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "maxpixel", &maxpixel);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "maxpixel_over_E", &maxpixel_over_E);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "cluster_rms_over_E", &cluster_rms_over_E);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "tmoment3", &tmoment3);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "tmoment4", &tmoment4);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "npixel", &npixel);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "npixel_red", &npixel_red);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "neighbors", &neighbors);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "logE", &logE);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "logrange", &logrange);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "logcluster_rms", &logcluster_rms);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "logmaxpixel", &logmaxpixel);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "lognpixel", &lognpixel);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "lognpixel_red", &lognpixel_red);
      TObject::ShowMembers(R__insp);
}

} // namespace MultiVariate
namespace ROOT {
   // Wrappers around operator new
   static void *new_MultiVariatecLcLMultiVariateResult(void *p) {
      return  p ? new(p) ::MultiVariate::MultiVariateResult : new ::MultiVariate::MultiVariateResult;
   }
   static void *newArray_MultiVariatecLcLMultiVariateResult(Long_t nElements, void *p) {
      return p ? new(p) ::MultiVariate::MultiVariateResult[nElements] : new ::MultiVariate::MultiVariateResult[nElements];
   }
   // Wrapper around operator delete
   static void delete_MultiVariatecLcLMultiVariateResult(void *p) {
      delete ((::MultiVariate::MultiVariateResult*)p);
   }
   static void deleteArray_MultiVariatecLcLMultiVariateResult(void *p) {
      delete [] ((::MultiVariate::MultiVariateResult*)p);
   }
   static void destruct_MultiVariatecLcLMultiVariateResult(void *p) {
      typedef ::MultiVariate::MultiVariateResult current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_MultiVariatecLcLMultiVariateResult(TBuffer &buf, void *obj) {
      ((::MultiVariate::MultiVariateResult*)obj)->::MultiVariate::MultiVariateResult::Streamer(buf);
   }
} // end of namespace ROOT for class ::MultiVariate::MultiVariateResult

/********************************************************
* MultiVariateResultCint.cc
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableMultiVariateResultCint();

extern "C" void G__set_cpp_environmentMultiVariateResultCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("MultiVariate.hh");
  G__cpp_reset_tagtableMultiVariateResultCint();
}
#include <new>
extern "C" int G__cpp_dllrevMultiVariateResultCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* MultiVariate::MultiVariateResult */
static int G__MultiVariateResultCint_863_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MultiVariate::MultiVariateResult* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new MultiVariate::MultiVariateResult[n];
     } else {
       p = new((void*) gvp) MultiVariate::MultiVariateResult[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new MultiVariate::MultiVariateResult;
     } else {
       p = new((void*) gvp) MultiVariate::MultiVariateResult;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 3:
      G__letint(result7, 103, (long) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->setResult((const char*) G__int(libp->para[0]), (const char*) G__int(libp->para[1])
, (const char*) G__int(libp->para[2])));
      break;
   case 2:
      G__letint(result7, 103, (long) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->setResult((const char*) G__int(libp->para[0]), (const char*) G__int(libp->para[1])));
      break;
   case 1:
      G__letint(result7, 103, (long) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->setResult((const char*) G__int(libp->para[0])));
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 4:
      G__letdouble(result7, 100, (double) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->getClassifier((DmtpcSkimEvent*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (const char*) G__int(libp->para[3])));
      break;
   case 3:
      G__letdouble(result7, 100, (double) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->getClassifier((DmtpcSkimEvent*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2])));
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 5:
      G__letint(result7, 103, (long) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->getCutsClassifier((DmtpcSkimEvent*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (double) G__double(libp->para[3])
, (const char*) G__int(libp->para[4])));
      break;
   case 4:
      G__letint(result7, 103, (long) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->getCutsClassifier((DmtpcSkimEvent*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (double) G__double(libp->para[3])));
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 5:
      G__letdouble(result7, 100, (double) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->getSignalProbability((DmtpcSkimEvent*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (double) G__double(libp->para[3])
, (const char*) G__int(libp->para[4])));
      break;
   case 4:
      G__letdouble(result7, 100, (double) ((MultiVariate::MultiVariateResult*) G__getstructoffset())->getSignalProbability((DmtpcSkimEvent*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (double) G__double(libp->para[3])));
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) MultiVariate::MultiVariateResult::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MultiVariate::MultiVariateResult::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) MultiVariate::MultiVariateResult::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      MultiVariate::MultiVariateResult::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MultiVariate::MultiVariateResult*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MultiVariate::MultiVariateResult::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MultiVariate::MultiVariateResult::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MultiVariate::MultiVariateResult::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MultiVariateResultCint_863_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MultiVariate::MultiVariateResult::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__MultiVariateResultCint_863_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   MultiVariate::MultiVariateResult* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new MultiVariate::MultiVariateResult(*(MultiVariate::MultiVariateResult*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef MultiVariate::MultiVariateResult G__TMultiVariatecLcLMultiVariateResult;
static int G__MultiVariateResultCint_863_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (MultiVariate::MultiVariateResult*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((MultiVariate::MultiVariateResult*) (soff+(sizeof(MultiVariate::MultiVariateResult)*i)))->~G__TMultiVariatecLcLMultiVariateResult();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (MultiVariate::MultiVariateResult*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((MultiVariate::MultiVariateResult*) (soff))->~G__TMultiVariatecLcLMultiVariateResult();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__MultiVariateResultCint_863_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MultiVariate::MultiVariateResult* dest = (MultiVariate::MultiVariateResult*) G__getstructoffset();
   *dest = *(MultiVariate::MultiVariateResult*) libp->para[0].ref;
   const MultiVariate::MultiVariateResult& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* MultiVariate::MultiVariateResult */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncMultiVariateResultCint {
 public:
  G__Sizep2memfuncMultiVariateResultCint(): p(&G__Sizep2memfuncMultiVariateResultCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncMultiVariateResultCint::*p)();
};

size_t G__get_sizep2memfuncMultiVariateResultCint()
{
  G__Sizep2memfuncMultiVariateResultCint a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceMultiVariateResultCint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult))) {
     MultiVariate::MultiVariateResult *G__Lderived;
     G__Lderived=(MultiVariate::MultiVariateResult*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult),G__get_linked_tagnum(&G__MultiVariateResultCintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableMultiVariateResultCint() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* MultiVariate::MultiVariateResult */
static void G__setup_memvarMultiVariatecLcLMultiVariateResult(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult));
   { MultiVariate::MultiVariateResult *p; p=(MultiVariate::MultiVariateResult*)0x1000; if (p) { }
   G__memvar_setup((void*)0,103,0,0,-1,-1,-1,4,"optimized=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__MultiVariateResultCintLN_TMVAcLcLReader),-1,-1,4,"reader=",0,(char*)NULL);
   G__memvar_setup((void*)0,67,0,0,-1,-1,-1,4,"current=",0,"!");
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"E=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"range=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"skewness=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"cluster_rms=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"cluster_mean=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"maxpixel=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"maxpixel_over_E=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"cluster_rms_over_E=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"tmoment3=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"tmoment4=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"npixel=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"npixel_red=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"neighbors=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"logE=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"logrange=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"logcluster_rms=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"logmaxpixel=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"lognpixel=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"lognpixel_red=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__MultiVariateResultCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarMultiVariateResultCint() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncMultiVariatecLcLMultiVariateResult(void) {
   /* MultiVariate::MultiVariateResult */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult));
   G__memfunc_setup("MultiVariateResult",1878,G__MultiVariateResultCint_863_0_1, 105, G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setResult",971,G__MultiVariateResultCint_863_0_2, 103, -1, -1, 0, 3, 1, 1, 0, 
"C - - 10 - job_name C - - 10 '\"Fisher\"' method_name "
"C - - 10 '\"weights/\"' path_to_weights_dir", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getClassifier",1349,G__MultiVariateResultCint_863_0_3, 100, -1, -1, 0, 4, 1, 1, 0, 
"U 'DmtpcSkimEvent' - 10 - ev i - - 0 - cam "
"i - - 0 - track C - - 10 'NULL' method_name", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getCutsClassifier",1764,G__MultiVariateResultCint_863_0_4, 103, -1, -1, 0, 5, 1, 1, 0, 
"U 'DmtpcSkimEvent' - 10 - ev i - - 0 - cam "
"i - - 0 - track d - - 0 - signal_efficency "
"C - - 10 '\"Cuts\"' method_name", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getSignalProbability",2079,G__MultiVariateResultCint_863_0_5, 100, -1, -1, 0, 5, 1, 1, 0, 
"U 'DmtpcSkimEvent' - 10 - ev i - - 0 - cam "
"i - - 0 - track d - - 0 - signal_fraction "
"C - - 10 'NULL' method_name", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setrdvars",990,(G__InterfaceMethod) NULL, 121, -1, -1, 0, 3, 1, 4, 0, 
"U 'DmtpcSkimEvent' - 10 - ev i - - 0 - cam "
"i - - 0 - track", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__MultiVariateResultCint_863_0_7, 85, G__get_linked_tagnum(&G__MultiVariateResultCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&MultiVariate::MultiVariateResult::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__MultiVariateResultCint_863_0_8, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MultiVariate::MultiVariateResult::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__MultiVariateResultCint_863_0_9, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&MultiVariate::MultiVariateResult::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__MultiVariateResultCint_863_0_10, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&MultiVariate::MultiVariateResult::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__MultiVariateResultCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__MultiVariateResultCint_863_0_14, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__MultiVariateResultCint_863_0_15, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MultiVariate::MultiVariateResult::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__MultiVariateResultCint_863_0_16, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MultiVariate::MultiVariateResult::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__MultiVariateResultCint_863_0_17, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MultiVariate::MultiVariateResult::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__MultiVariateResultCint_863_0_18, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MultiVariate::MultiVariateResult::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("MultiVariateResult", 1878, G__MultiVariateResultCint_863_0_19, (int) ('i'), G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult), -1, 0, 1, 1, 1, 0, "u 'MultiVariate::MultiVariateResult' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~MultiVariateResult", 2004, G__MultiVariateResultCint_863_0_20, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__MultiVariateResultCint_863_0_21, (int) ('u'), G__get_linked_tagnum(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult), -1, 1, 1, 1, 1, 0, "u 'MultiVariate::MultiVariateResult' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncMultiVariateResultCint() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {
}

static void G__cpp_setup_global3() {
}

static void G__cpp_setup_global4() {
}

static void G__cpp_setup_global5() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalMultiVariateResultCint() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
  G__cpp_setup_global3();
  G__cpp_setup_global4();
  G__cpp_setup_global5();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {
}

static void G__cpp_setup_func14() {
}

static void G__cpp_setup_func15() {
}

static void G__cpp_setup_func16() {
}

static void G__cpp_setup_func17() {
}

static void G__cpp_setup_func18() {
}

static void G__cpp_setup_func19() {
}

static void G__cpp_setup_func20() {
}

static void G__cpp_setup_func21() {
}

static void G__cpp_setup_func22() {
}

static void G__cpp_setup_func23() {
}

static void G__cpp_setup_func24() {
}

static void G__cpp_setup_func25() {
}

static void G__cpp_setup_func26() {
}

static void G__cpp_setup_func27() {
}

static void G__cpp_setup_func28() {
}

static void G__cpp_setup_func29() {
}

static void G__cpp_setup_func30() {
}

static void G__cpp_setup_func31() {
}

static void G__cpp_setup_func32() {
}

static void G__cpp_setup_func33() {
}

static void G__cpp_setup_func34() {
}

static void G__cpp_setup_func35() {
}

static void G__cpp_setup_func36() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcMultiVariateResultCint() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
  G__cpp_setup_func14();
  G__cpp_setup_func15();
  G__cpp_setup_func16();
  G__cpp_setup_func17();
  G__cpp_setup_func18();
  G__cpp_setup_func19();
  G__cpp_setup_func20();
  G__cpp_setup_func21();
  G__cpp_setup_func22();
  G__cpp_setup_func23();
  G__cpp_setup_func24();
  G__cpp_setup_func25();
  G__cpp_setup_func26();
  G__cpp_setup_func27();
  G__cpp_setup_func28();
  G__cpp_setup_func29();
  G__cpp_setup_func30();
  G__cpp_setup_func31();
  G__cpp_setup_func32();
  G__cpp_setup_func33();
  G__cpp_setup_func34();
  G__cpp_setup_func35();
  G__cpp_setup_func36();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__MultiVariateResultCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__MultiVariateResultCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__MultiVariateResultCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__MultiVariateResultCintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__MultiVariateResultCintLN_TMVA = { "TMVA" , 110 , -1 };
G__linked_taginfo G__MultiVariateResultCintLN_TMVAcLcLReader = { "TMVA::Reader" , 99 , -1 };
G__linked_taginfo G__MultiVariateResultCintLN_DmtpcSkimEvent = { "DmtpcSkimEvent" , 99 , -1 };
G__linked_taginfo G__MultiVariateResultCintLN_MultiVariate = { "MultiVariate" , 110 , -1 };
G__linked_taginfo G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult = { "MultiVariate::MultiVariateResult" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableMultiVariateResultCint() {
  G__MultiVariateResultCintLN_TClass.tagnum = -1 ;
  G__MultiVariateResultCintLN_TBuffer.tagnum = -1 ;
  G__MultiVariateResultCintLN_TMemberInspector.tagnum = -1 ;
  G__MultiVariateResultCintLN_TObject.tagnum = -1 ;
  G__MultiVariateResultCintLN_TMVA.tagnum = -1 ;
  G__MultiVariateResultCintLN_TMVAcLcLReader.tagnum = -1 ;
  G__MultiVariateResultCintLN_DmtpcSkimEvent.tagnum = -1 ;
  G__MultiVariateResultCintLN_MultiVariate.tagnum = -1 ;
  G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableMultiVariateResultCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_TObject);
   G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_TMVA);
   G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_TMVAcLcLReader);
   G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_DmtpcSkimEvent);
   G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_MultiVariate);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__MultiVariateResultCintLN_MultiVariatecLcLMultiVariateResult),sizeof(MultiVariate::MultiVariateResult),-1,29952,(char*)NULL,G__setup_memvarMultiVariatecLcLMultiVariateResult,G__setup_memfuncMultiVariatecLcLMultiVariateResult);
}
extern "C" void G__cpp_setupMultiVariateResultCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupMultiVariateResultCint()");
  G__set_cpp_environmentMultiVariateResultCint();
  G__cpp_setup_tagtableMultiVariateResultCint();

  G__cpp_setup_inheritanceMultiVariateResultCint();

  G__cpp_setup_typetableMultiVariateResultCint();

  G__cpp_setup_memvarMultiVariateResultCint();

  G__cpp_setup_memfuncMultiVariateResultCint();
  G__cpp_setup_globalMultiVariateResultCint();
  G__cpp_setup_funcMultiVariateResultCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncMultiVariateResultCint();
  return;
}
class G__cpp_setup_initMultiVariateResultCint {
  public:
    G__cpp_setup_initMultiVariateResultCint() { G__add_setup_func("MultiVariateResultCint",(G__incsetup)(&G__cpp_setupMultiVariateResultCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initMultiVariateResultCint() { G__remove_setup_func("MultiVariateResultCint"); }
};
G__cpp_setup_initMultiVariateResultCint G__cpp_setup_initializerMultiVariateResultCint;

