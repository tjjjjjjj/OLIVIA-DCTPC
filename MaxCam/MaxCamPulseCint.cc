//
// File generated by rootcint at Tue Aug 16 18:39:27 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME MaxCamPulseCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "MaxCamPulseCint.h"

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
   void MaxCamPulse_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_MaxCamPulse(void *p = 0);
   static void *newArray_MaxCamPulse(Long_t size, void *p);
   static void delete_MaxCamPulse(void *p);
   static void deleteArray_MaxCamPulse(void *p);
   static void destruct_MaxCamPulse(void *p);
   static void streamer_MaxCamPulse(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MaxCamPulse*)
   {
      ::MaxCamPulse *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MaxCamPulse >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MaxCamPulse", ::MaxCamPulse::Class_Version(), "./MaxCamPulse.hh", 9,
                  typeid(::MaxCamPulse), DefineBehavior(ptr, ptr),
                  &::MaxCamPulse::Dictionary, isa_proxy, 0,
                  sizeof(::MaxCamPulse) );
      instance.SetNew(&new_MaxCamPulse);
      instance.SetNewArray(&newArray_MaxCamPulse);
      instance.SetDelete(&delete_MaxCamPulse);
      instance.SetDeleteArray(&deleteArray_MaxCamPulse);
      instance.SetDestructor(&destruct_MaxCamPulse);
      instance.SetStreamerFunc(&streamer_MaxCamPulse);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MaxCamPulse*)
   {
      return GenerateInitInstanceLocal((::MaxCamPulse*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MaxCamPulse*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *MaxCamPulse::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *MaxCamPulse::Class_Name()
{
   return "MaxCamPulse";
}

//______________________________________________________________________________
const char *MaxCamPulse::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MaxCamPulse*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MaxCamPulse::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MaxCamPulse*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void MaxCamPulse::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MaxCamPulse*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *MaxCamPulse::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MaxCamPulse*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void MaxCamPulse::Streamer(TBuffer &R__b)
{
   // Stream an object of class MaxCamPulse.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> init;
      R__b >> _nbin;
      R__b >> _pulseHeight;
      R__b >> _pulseHeightTime;
      R__b >> _pulseStartTime;
      R__b >> _pulseStartBin;
      R__b >> _pulseEndTime;
      R__b >> _pulseEndBin;
      R__b >> _pulseIntegral;
      R__b.CheckByteCount(R__s, R__c, MaxCamPulse::IsA());
   } else {
      R__c = R__b.WriteVersion(MaxCamPulse::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << init;
      R__b << _nbin;
      R__b << _pulseHeight;
      R__b << _pulseHeightTime;
      R__b << _pulseStartTime;
      R__b << _pulseStartBin;
      R__b << _pulseEndTime;
      R__b << _pulseEndBin;
      R__b << _pulseIntegral;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void MaxCamPulse::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class MaxCamPulse.
      TClass *R__cl = ::MaxCamPulse::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "init", &init);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_nbin", &_nbin);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_pulseHeight", &_pulseHeight);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_pulseHeightTime", &_pulseHeightTime);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_pulseStartTime", &_pulseStartTime);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_pulseStartBin", &_pulseStartBin);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_pulseEndTime", &_pulseEndTime);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_pulseEndBin", &_pulseEndBin);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_pulseIntegral", &_pulseIntegral);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MaxCamPulse(void *p) {
      return  p ? new(p) ::MaxCamPulse : new ::MaxCamPulse;
   }
   static void *newArray_MaxCamPulse(Long_t nElements, void *p) {
      return p ? new(p) ::MaxCamPulse[nElements] : new ::MaxCamPulse[nElements];
   }
   // Wrapper around operator delete
   static void delete_MaxCamPulse(void *p) {
      delete ((::MaxCamPulse*)p);
   }
   static void deleteArray_MaxCamPulse(void *p) {
      delete [] ((::MaxCamPulse*)p);
   }
   static void destruct_MaxCamPulse(void *p) {
      typedef ::MaxCamPulse current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_MaxCamPulse(TBuffer &buf, void *obj) {
      ((::MaxCamPulse*)obj)->::MaxCamPulse::Streamer(buf);
   }
} // end of namespace ROOT for class ::MaxCamPulse

/********************************************************
* MaxCamPulseCint.cc
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

extern "C" void G__cpp_reset_tagtableMaxCamPulseCint();

extern "C" void G__set_cpp_environmentMaxCamPulseCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("MaxCamPulse.hh");
  G__cpp_reset_tagtableMaxCamPulseCint();
}
#include <new>
extern "C" int G__cpp_dllrevMaxCamPulseCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* MaxCamPulse */
static int G__MaxCamPulseCint_207_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MaxCamPulse* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new MaxCamPulse[n];
     } else {
       p = new((void*) gvp) MaxCamPulse[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new MaxCamPulse;
     } else {
       p = new((void*) gvp) MaxCamPulse;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MaxCamPulse* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new MaxCamPulse((Int_t) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) MaxCamPulse((Int_t) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((MaxCamPulse*) G__getstructoffset())->getBin());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamPulse*) G__getstructoffset())->setPulseHeight((Double_t) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((MaxCamPulse*) G__getstructoffset())->getPulseHeight());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamPulse*) G__getstructoffset())->setPulseHeightTime((Double_t) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((MaxCamPulse*) G__getstructoffset())->getPulseHeightTime());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamPulse*) G__getstructoffset())->setPulseStartTime((Double_t) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((MaxCamPulse*) G__getstructoffset())->getPulseStartTime());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamPulse*) G__getstructoffset())->setPulseStartBin((Double_t) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((MaxCamPulse*) G__getstructoffset())->getPulseStartBin());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamPulse*) G__getstructoffset())->setPulseEndTime((Double_t) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((MaxCamPulse*) G__getstructoffset())->getPulseEndTime());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamPulse*) G__getstructoffset())->setPulseEndBin((Double_t) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((MaxCamPulse*) G__getstructoffset())->getPulseEndBin());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamPulse*) G__getstructoffset())->setPulseIntegral((Double_t) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((MaxCamPulse*) G__getstructoffset())->getPulseIntegral());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 1:
      ((MaxCamPulse*) G__getstructoffset())->print(*(ostream*) libp->para[0].ref);
      G__setnull(result7);
      break;
   case 0:
      ((MaxCamPulse*) G__getstructoffset())->print();
      G__setnull(result7);
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) MaxCamPulse::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamPulse::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) MaxCamPulse::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      MaxCamPulse::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_26(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamPulse*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_27(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamPulse::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_28(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamPulse::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_29(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamPulse::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamPulseCint_207_0_30(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamPulse::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__MaxCamPulseCint_207_0_31(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   MaxCamPulse* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new MaxCamPulse(*(MaxCamPulse*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef MaxCamPulse G__TMaxCamPulse;
static int G__MaxCamPulseCint_207_0_32(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (MaxCamPulse*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((MaxCamPulse*) (soff+(sizeof(MaxCamPulse)*i)))->~G__TMaxCamPulse();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (MaxCamPulse*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((MaxCamPulse*) (soff))->~G__TMaxCamPulse();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__MaxCamPulseCint_207_0_33(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MaxCamPulse* dest = (MaxCamPulse*) G__getstructoffset();
   *dest = *(MaxCamPulse*) libp->para[0].ref;
   const MaxCamPulse& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* MaxCamPulse */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncMaxCamPulseCint {
 public:
  G__Sizep2memfuncMaxCamPulseCint(): p(&G__Sizep2memfuncMaxCamPulseCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncMaxCamPulseCint::*p)();
};

size_t G__get_sizep2memfuncMaxCamPulseCint()
{
  G__Sizep2memfuncMaxCamPulseCint a;
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
extern "C" void G__cpp_setup_inheritanceMaxCamPulseCint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse))) {
     MaxCamPulse *G__Lderived;
     G__Lderived=(MaxCamPulse*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse),G__get_linked_tagnum(&G__MaxCamPulseCintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableMaxCamPulseCint() {

   /* Setting up typedef entry */
   G__search_typename2("Int_t",105,-1,0,-1);
   G__setnewtype(-1,"Signed integer 4 bytes (int)",0);
   G__search_typename2("Double_t",100,-1,0,-1);
   G__setnewtype(-1,"Double 8 bytes",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__MaxCamPulseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__MaxCamPulseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamPulseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__MaxCamPulseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamPulseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__MaxCamPulseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__MaxCamPulseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamPulseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__MaxCamPulseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamPulseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Float_t>",117,G__get_linked_tagnum(&G__MaxCamPulseCintLN_TVectorTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Double_t>",117,G__get_linked_tagnum(&G__MaxCamPulseCintLN_TVectorTlEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* MaxCamPulse */
static void G__setup_memvarMaxCamPulse(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse));
   { MaxCamPulse *p; p=(MaxCamPulse*)0x1000; if (p) { }
   G__memvar_setup((void*)0,103,0,0,-1,-1,-1,4,"init=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,G__defined_typename("Int_t"),-1,4,"_nbin=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"_pulseHeight=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"_pulseHeightTime=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"_pulseStartTime=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"_pulseStartBin=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"_pulseEndTime=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"_pulseEndBin=",0,(char*)NULL);
   G__memvar_setup((void*)0,100,0,0,-1,G__defined_typename("Double_t"),-1,4,"_pulseIntegral=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__MaxCamPulseCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarMaxCamPulseCint() {
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
static void G__setup_memfuncMaxCamPulse(void) {
   /* MaxCamPulse */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse));
   G__memfunc_setup("MaxCamPulse",1088,G__MaxCamPulseCint_207_0_1, 105, G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("MaxCamPulse",1088,G__MaxCamPulseCint_207_0_2, 105, G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse), -1, 0, 1, 1, 1, 0, "i - 'Int_t' 0 - nbin", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getBin",601,G__MaxCamPulseCint_207_0_3, 105, -1, G__defined_typename("Int_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPulseHeight",1454,G__MaxCamPulseCint_207_0_4, 121, -1, -1, 0, 1, 1, 1, 0, "d - 'Double_t' 0 - value", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPulseHeight",1442,G__MaxCamPulseCint_207_0_5, 100, -1, G__defined_typename("Double_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPulseHeightTime",1853,G__MaxCamPulseCint_207_0_6, 121, -1, -1, 0, 1, 1, 1, 0, "d - 'Double_t' 0 - value", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPulseHeightTime",1841,G__MaxCamPulseCint_207_0_7, 100, -1, G__defined_typename("Double_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPulseStartTime",1778,G__MaxCamPulseCint_207_0_8, 121, -1, -1, 0, 1, 1, 1, 0, "d - 'Double_t' 0 - value", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPulseStartTime",1766,G__MaxCamPulseCint_207_0_9, 100, -1, G__defined_typename("Double_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPulseStartBin",1660,G__MaxCamPulseCint_207_0_10, 121, -1, -1, 0, 1, 1, 1, 0, "d - 'Double_t' 0 - value", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPulseStartBin",1648,G__MaxCamPulseCint_207_0_11, 100, -1, G__defined_typename("Double_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPulseEndTime",1531,G__MaxCamPulseCint_207_0_12, 121, -1, -1, 0, 1, 1, 1, 0, "d - 'Double_t' 0 - value", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPulseEndTime",1519,G__MaxCamPulseCint_207_0_13, 100, -1, G__defined_typename("Double_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPulseEndBin",1413,G__MaxCamPulseCint_207_0_14, 121, -1, -1, 0, 1, 1, 1, 0, "d - 'Double_t' 0 - value", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPulseEndBin",1401,G__MaxCamPulseCint_207_0_15, 100, -1, G__defined_typename("Double_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPulseIntegral",1675,G__MaxCamPulseCint_207_0_16, 121, -1, -1, 0, 1, 1, 1, 0, "d - 'Double_t' 0 - value", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPulseIntegral",1663,G__MaxCamPulseCint_207_0_17, 100, -1, G__defined_typename("Double_t"), 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("print",557,G__MaxCamPulseCint_207_0_18, 121, -1, -1, 0, 1, 1, 1, 0, "u 'basic_ostream<char,char_traits<char> >' 'ostream' 1 'std::cout' out", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__MaxCamPulseCint_207_0_19, 85, G__get_linked_tagnum(&G__MaxCamPulseCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&MaxCamPulse::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__MaxCamPulseCint_207_0_20, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamPulse::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__MaxCamPulseCint_207_0_21, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&MaxCamPulse::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__MaxCamPulseCint_207_0_22, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&MaxCamPulse::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__MaxCamPulseCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__MaxCamPulseCint_207_0_26, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__MaxCamPulseCint_207_0_27, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamPulse::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__MaxCamPulseCint_207_0_28, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MaxCamPulse::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__MaxCamPulseCint_207_0_29, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamPulse::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__MaxCamPulseCint_207_0_30, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MaxCamPulse::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("MaxCamPulse", 1088, G__MaxCamPulseCint_207_0_31, (int) ('i'), G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse), -1, 0, 1, 1, 1, 0, "u 'MaxCamPulse' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~MaxCamPulse", 1214, G__MaxCamPulseCint_207_0_32, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__MaxCamPulseCint_207_0_33, (int) ('u'), G__get_linked_tagnum(&G__MaxCamPulseCintLN_MaxCamPulse), -1, 1, 1, 1, 1, 0, "u 'MaxCamPulse' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncMaxCamPulseCint() {
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

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalMaxCamPulseCint() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
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

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcMaxCamPulseCint() {
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
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__MaxCamPulseCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_basic_ostreamlEcharcOchar_traitslEchargRsPgR = { "basic_ostream<char,char_traits<char> >" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_TVectorTlEfloatgR = { "TVectorT<float>" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_TVectorTlEdoublegR = { "TVectorT<double>" , 99 , -1 };
G__linked_taginfo G__MaxCamPulseCintLN_MaxCamPulse = { "MaxCamPulse" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableMaxCamPulseCint() {
  G__MaxCamPulseCintLN_TClass.tagnum = -1 ;
  G__MaxCamPulseCintLN_TBuffer.tagnum = -1 ;
  G__MaxCamPulseCintLN_TMemberInspector.tagnum = -1 ;
  G__MaxCamPulseCintLN_TObject.tagnum = -1 ;
  G__MaxCamPulseCintLN_basic_ostreamlEcharcOchar_traitslEchargRsPgR.tagnum = -1 ;
  G__MaxCamPulseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__MaxCamPulseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__MaxCamPulseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__MaxCamPulseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__MaxCamPulseCintLN_TVectorTlEfloatgR.tagnum = -1 ;
  G__MaxCamPulseCintLN_TVectorTlEdoublegR.tagnum = -1 ;
  G__MaxCamPulseCintLN_MaxCamPulse.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableMaxCamPulseCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_TObject);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_basic_ostreamlEcharcOchar_traitslEchargRsPgR);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_TVectorTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_TVectorTlEdoublegR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__MaxCamPulseCintLN_MaxCamPulse),sizeof(MaxCamPulse),-1,62720,(char*)NULL,G__setup_memvarMaxCamPulse,G__setup_memfuncMaxCamPulse);
}
extern "C" void G__cpp_setupMaxCamPulseCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupMaxCamPulseCint()");
  G__set_cpp_environmentMaxCamPulseCint();
  G__cpp_setup_tagtableMaxCamPulseCint();

  G__cpp_setup_inheritanceMaxCamPulseCint();

  G__cpp_setup_typetableMaxCamPulseCint();

  G__cpp_setup_memvarMaxCamPulseCint();

  G__cpp_setup_memfuncMaxCamPulseCint();
  G__cpp_setup_globalMaxCamPulseCint();
  G__cpp_setup_funcMaxCamPulseCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncMaxCamPulseCint();
  return;
}
class G__cpp_setup_initMaxCamPulseCint {
  public:
    G__cpp_setup_initMaxCamPulseCint() { G__add_setup_func("MaxCamPulseCint",(G__incsetup)(&G__cpp_setupMaxCamPulseCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initMaxCamPulseCint() { G__remove_setup_func("MaxCamPulseCint"); }
};
G__cpp_setup_initMaxCamPulseCint G__cpp_setup_initializerMaxCamPulseCint;

