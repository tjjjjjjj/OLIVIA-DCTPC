//
// File generated by rootcint at Fri Jul 15 15:53:49 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME McDarkDigiCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "McDarkDigiCint.h"

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
   void McDarkDigi_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_McDarkDigi(void *p = 0);
   static void *newArray_McDarkDigi(Long_t size, void *p);
   static void delete_McDarkDigi(void *p);
   static void deleteArray_McDarkDigi(void *p);
   static void destruct_McDarkDigi(void *p);
   static void streamer_McDarkDigi(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::McDarkDigi*)
   {
      ::McDarkDigi *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::McDarkDigi >(0);
      static ::ROOT::TGenericClassInfo 
         instance("McDarkDigi", ::McDarkDigi::Class_Version(), "./McDarkDigi.hh", 8,
                  typeid(::McDarkDigi), DefineBehavior(ptr, ptr),
                  &::McDarkDigi::Dictionary, isa_proxy, 0,
                  sizeof(::McDarkDigi) );
      instance.SetNew(&new_McDarkDigi);
      instance.SetNewArray(&newArray_McDarkDigi);
      instance.SetDelete(&delete_McDarkDigi);
      instance.SetDeleteArray(&deleteArray_McDarkDigi);
      instance.SetDestructor(&destruct_McDarkDigi);
      instance.SetStreamerFunc(&streamer_McDarkDigi);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::McDarkDigi*)
   {
      return GenerateInitInstanceLocal((::McDarkDigi*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::McDarkDigi*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *McDarkDigi::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *McDarkDigi::Class_Name()
{
   return "McDarkDigi";
}

//______________________________________________________________________________
const char *McDarkDigi::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::McDarkDigi*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int McDarkDigi::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::McDarkDigi*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void McDarkDigi::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::McDarkDigi*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *McDarkDigi::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::McDarkDigi*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void McDarkDigi::Streamer(TBuffer &R__b)
{
   // Stream an object of class McDarkDigi.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> _moduleID;
      R__b >> _channelID;
      R__b >> _weight;
      R__b >> _trackIndex;
      R__b.CheckByteCount(R__s, R__c, McDarkDigi::IsA());
   } else {
      R__c = R__b.WriteVersion(McDarkDigi::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << _moduleID;
      R__b << _channelID;
      R__b << _weight;
      R__b << _trackIndex;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void McDarkDigi::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class McDarkDigi.
      TClass *R__cl = ::McDarkDigi::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_moduleID", &_moduleID);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_channelID", &_channelID);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_weight", &_weight);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_trackIndex", &_trackIndex);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_McDarkDigi(void *p) {
      return  p ? new(p) ::McDarkDigi : new ::McDarkDigi;
   }
   static void *newArray_McDarkDigi(Long_t nElements, void *p) {
      return p ? new(p) ::McDarkDigi[nElements] : new ::McDarkDigi[nElements];
   }
   // Wrapper around operator delete
   static void delete_McDarkDigi(void *p) {
      delete ((::McDarkDigi*)p);
   }
   static void deleteArray_McDarkDigi(void *p) {
      delete [] ((::McDarkDigi*)p);
   }
   static void destruct_McDarkDigi(void *p) {
      typedef ::McDarkDigi current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_McDarkDigi(TBuffer &buf, void *obj) {
      ((::McDarkDigi*)obj)->::McDarkDigi::Streamer(buf);
   }
} // end of namespace ROOT for class ::McDarkDigi

/********************************************************
* McDarkDigiCint.cc
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

extern "C" void G__cpp_reset_tagtableMcDarkDigiCint();

extern "C" void G__set_cpp_environmentMcDarkDigiCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("McDarkDigi.hh");
  G__cpp_reset_tagtableMcDarkDigiCint();
}
#include <new>
extern "C" int G__cpp_dllrevMcDarkDigiCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* McDarkDigi */
static int G__McDarkDigiCint_214_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   McDarkDigi* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new McDarkDigi[n];
     } else {
       p = new((void*) gvp) McDarkDigi[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new McDarkDigi;
     } else {
       p = new((void*) gvp) McDarkDigi;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   McDarkDigi* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new McDarkDigi(*(McDarkDigi*) libp->para[0].ref);
   } else {
     p = new((void*) gvp) McDarkDigi(*(McDarkDigi*) libp->para[0].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((McDarkDigi*) G__getstructoffset())->setModuleID((int) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((McDarkDigi*) G__getstructoffset())->getModuleID());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((McDarkDigi*) G__getstructoffset())->setChannelID((int) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((McDarkDigi*) G__getstructoffset())->getChannelID());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((McDarkDigi*) G__getstructoffset())->setWeight((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((McDarkDigi*) G__getstructoffset())->getWeight());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((McDarkDigi*) G__getstructoffset())->setTrackID((int) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((McDarkDigi*) G__getstructoffset())->getTrackID());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) McDarkDigi::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) McDarkDigi::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) McDarkDigi::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      McDarkDigi::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((McDarkDigi*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) McDarkDigi::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) McDarkDigi::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) McDarkDigi::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__McDarkDigiCint_214_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) McDarkDigi::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef McDarkDigi G__TMcDarkDigi;
static int G__McDarkDigiCint_214_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (McDarkDigi*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((McDarkDigi*) (soff+(sizeof(McDarkDigi)*i)))->~G__TMcDarkDigi();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (McDarkDigi*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((McDarkDigi*) (soff))->~G__TMcDarkDigi();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__McDarkDigiCint_214_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   McDarkDigi* dest = (McDarkDigi*) G__getstructoffset();
   *dest = *(McDarkDigi*) libp->para[0].ref;
   const McDarkDigi& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* McDarkDigi */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncMcDarkDigiCint {
 public:
  G__Sizep2memfuncMcDarkDigiCint(): p(&G__Sizep2memfuncMcDarkDigiCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncMcDarkDigiCint::*p)();
};

size_t G__get_sizep2memfuncMcDarkDigiCint()
{
  G__Sizep2memfuncMcDarkDigiCint a;
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
extern "C" void G__cpp_setup_inheritanceMcDarkDigiCint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi))) {
     McDarkDigi *G__Lderived;
     G__Lderived=(McDarkDigi*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi),G__get_linked_tagnum(&G__McDarkDigiCintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableMcDarkDigiCint() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__McDarkDigiCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__McDarkDigiCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__McDarkDigiCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__McDarkDigiCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTBase<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTBaselEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTBase<Double_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTBaselEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TVectorTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Double_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TVectorTlEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixT<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTRow_const<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTRow_constlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTColumn_const<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTColumn_constlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTDiag_const<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTDiag_constlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTFlat_const<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTFlat_constlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTSub_const<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTSub_constlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTSparseRow_const<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTSparseRow_constlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTSparseDiag_const<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTSparseDiag_constlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTRow<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTRowlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTColumn<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTColumnlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTDiag<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTDiaglEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTFlat<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTFlatlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTSub<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTSublEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTSparseRow<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTSparseRowlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTSparseDiag<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TMatrixTSparseDiaglEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TElementActionT<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TElementActionTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TElementPosActionT<Float_t>",117,G__get_linked_tagnum(&G__McDarkDigiCintLN_TElementPosActionTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* McDarkDigi */
static void G__setup_memvarMcDarkDigi(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi));
   { McDarkDigi *p; p=(McDarkDigi*)0x1000; if (p) { }
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"_moduleID=",0,"module number, e.g. camera number or PMT number for this TPC");
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"_channelID=",0,"pixel-bin in CCD, time-bin in charge/PMT readout");
   G__memvar_setup((void*)0,100,0,0,-1,-1,-1,4,"_weight=",0,"ADC value");
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"_trackIndex=",0,"index of track causing this digi");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__McDarkDigiCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarMcDarkDigiCint() {
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
static void G__setup_memfuncMcDarkDigi(void) {
   /* McDarkDigi */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi));
   G__memfunc_setup("McDarkDigi",943,G__McDarkDigiCint_214_0_1, 105, G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("McDarkDigi",943,G__McDarkDigiCint_214_0_2, 105, G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi), -1, 0, 1, 1, 1, 0, "u 'McDarkDigi' - 11 - other", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setModuleID",1087,G__McDarkDigiCint_214_0_3, 121, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - m", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getModuleID",1075,G__McDarkDigiCint_214_0_4, 105, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setChannelID",1170,G__McDarkDigiCint_214_0_5, 121, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - c", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getChannelID",1158,G__McDarkDigiCint_214_0_6, 105, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setWeight",948,G__McDarkDigiCint_214_0_7, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - w", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getWeight",936,G__McDarkDigiCint_214_0_8, 100, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setTrackID",974,G__McDarkDigiCint_214_0_9, 121, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - index", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getTrackID",962,G__McDarkDigiCint_214_0_10, 105, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__McDarkDigiCint_214_0_11, 85, G__get_linked_tagnum(&G__McDarkDigiCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&McDarkDigi::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__McDarkDigiCint_214_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&McDarkDigi::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__McDarkDigiCint_214_0_13, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&McDarkDigi::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__McDarkDigiCint_214_0_14, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&McDarkDigi::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__McDarkDigiCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__McDarkDigiCint_214_0_18, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__McDarkDigiCint_214_0_19, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&McDarkDigi::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__McDarkDigiCint_214_0_20, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&McDarkDigi::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__McDarkDigiCint_214_0_21, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&McDarkDigi::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__McDarkDigiCint_214_0_22, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&McDarkDigi::DeclFileLine) ), 0);
   // automatic destructor
   G__memfunc_setup("~McDarkDigi", 1069, G__McDarkDigiCint_214_0_23, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__McDarkDigiCint_214_0_24, (int) ('u'), G__get_linked_tagnum(&G__McDarkDigiCintLN_McDarkDigi), -1, 1, 1, 1, 1, 0, "u 'McDarkDigi' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncMcDarkDigiCint() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalMcDarkDigiCint() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
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

extern "C" void G__cpp_setup_funcMcDarkDigiCint() {
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
G__linked_taginfo G__McDarkDigiCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTBaselEfloatgR = { "TMatrixTBase<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTBaselEdoublegR = { "TMatrixTBase<double>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TVectorTlEfloatgR = { "TVectorT<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TVectorTlEdoublegR = { "TVectorT<double>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TElementActionTlEfloatgR = { "TElementActionT<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TElementPosActionTlEfloatgR = { "TElementPosActionT<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTlEfloatgR = { "TMatrixT<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTRow_constlEfloatgR = { "TMatrixTRow_const<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTRowlEfloatgR = { "TMatrixTRow<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTDiag_constlEfloatgR = { "TMatrixTDiag_const<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTColumn_constlEfloatgR = { "TMatrixTColumn_const<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTFlat_constlEfloatgR = { "TMatrixTFlat_const<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTSub_constlEfloatgR = { "TMatrixTSub_const<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTSparseRow_constlEfloatgR = { "TMatrixTSparseRow_const<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTSparseDiag_constlEfloatgR = { "TMatrixTSparseDiag_const<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTColumnlEfloatgR = { "TMatrixTColumn<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTDiaglEfloatgR = { "TMatrixTDiag<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTFlatlEfloatgR = { "TMatrixTFlat<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTSublEfloatgR = { "TMatrixTSub<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTSparseRowlEfloatgR = { "TMatrixTSparseRow<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_TMatrixTSparseDiaglEfloatgR = { "TMatrixTSparseDiag<float>" , 99 , -1 };
G__linked_taginfo G__McDarkDigiCintLN_McDarkDigi = { "McDarkDigi" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableMcDarkDigiCint() {
  G__McDarkDigiCintLN_TClass.tagnum = -1 ;
  G__McDarkDigiCintLN_TBuffer.tagnum = -1 ;
  G__McDarkDigiCintLN_TMemberInspector.tagnum = -1 ;
  G__McDarkDigiCintLN_TObject.tagnum = -1 ;
  G__McDarkDigiCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__McDarkDigiCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__McDarkDigiCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__McDarkDigiCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTBaselEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTBaselEdoublegR.tagnum = -1 ;
  G__McDarkDigiCintLN_TVectorTlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TVectorTlEdoublegR.tagnum = -1 ;
  G__McDarkDigiCintLN_TElementActionTlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TElementPosActionTlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTRow_constlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTRowlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTDiag_constlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTColumn_constlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTFlat_constlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTSub_constlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTSparseRow_constlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTSparseDiag_constlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTColumnlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTDiaglEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTFlatlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTSublEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTSparseRowlEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_TMatrixTSparseDiaglEfloatgR.tagnum = -1 ;
  G__McDarkDigiCintLN_McDarkDigi.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableMcDarkDigiCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TObject);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTBaselEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTBaselEdoublegR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TVectorTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TVectorTlEdoublegR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TElementActionTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TElementPosActionTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTRow_constlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTRowlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTDiag_constlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTColumn_constlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTFlat_constlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTSub_constlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTSparseRow_constlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTSparseDiag_constlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTColumnlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTDiaglEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTFlatlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTSublEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTSparseRowlEfloatgR);
   G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_TMatrixTSparseDiaglEfloatgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__McDarkDigiCintLN_McDarkDigi),sizeof(McDarkDigi),-1,30464,(char*)NULL,G__setup_memvarMcDarkDigi,G__setup_memfuncMcDarkDigi);
}
extern "C" void G__cpp_setupMcDarkDigiCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupMcDarkDigiCint()");
  G__set_cpp_environmentMcDarkDigiCint();
  G__cpp_setup_tagtableMcDarkDigiCint();

  G__cpp_setup_inheritanceMcDarkDigiCint();

  G__cpp_setup_typetableMcDarkDigiCint();

  G__cpp_setup_memvarMcDarkDigiCint();

  G__cpp_setup_memfuncMcDarkDigiCint();
  G__cpp_setup_globalMcDarkDigiCint();
  G__cpp_setup_funcMcDarkDigiCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncMcDarkDigiCint();
  return;
}
class G__cpp_setup_initMcDarkDigiCint {
  public:
    G__cpp_setup_initMcDarkDigiCint() { G__add_setup_func("McDarkDigiCint",(G__incsetup)(&G__cpp_setupMcDarkDigiCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initMcDarkDigiCint() { G__remove_setup_func("McDarkDigiCint"); }
};
G__cpp_setup_initMcDarkDigiCint G__cpp_setup_initializerMcDarkDigiCint;

