//
// File generated by rootcint at Mon Jun  6 11:01:32 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME DmtpcLensCorrectionCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "DmtpcLensCorrectionCint.h"

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
   void DmtpcLensCorrection_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_DmtpcLensCorrection(void *p = 0);
   static void *newArray_DmtpcLensCorrection(Long_t size, void *p);
   static void delete_DmtpcLensCorrection(void *p);
   static void deleteArray_DmtpcLensCorrection(void *p);
   static void destruct_DmtpcLensCorrection(void *p);
   static void streamer_DmtpcLensCorrection(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DmtpcLensCorrection*)
   {
      ::DmtpcLensCorrection *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DmtpcLensCorrection >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DmtpcLensCorrection", ::DmtpcLensCorrection::Class_Version(), "./DmtpcLensCorrection.hh", 28,
                  typeid(::DmtpcLensCorrection), DefineBehavior(ptr, ptr),
                  &::DmtpcLensCorrection::Dictionary, isa_proxy, 0,
                  sizeof(::DmtpcLensCorrection) );
      instance.SetNew(&new_DmtpcLensCorrection);
      instance.SetNewArray(&newArray_DmtpcLensCorrection);
      instance.SetDelete(&delete_DmtpcLensCorrection);
      instance.SetDeleteArray(&deleteArray_DmtpcLensCorrection);
      instance.SetDestructor(&destruct_DmtpcLensCorrection);
      instance.SetStreamerFunc(&streamer_DmtpcLensCorrection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DmtpcLensCorrection*)
   {
      return GenerateInitInstanceLocal((::DmtpcLensCorrection*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::DmtpcLensCorrection*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *DmtpcLensCorrection::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *DmtpcLensCorrection::Class_Name()
{
   return "DmtpcLensCorrection";
}

//______________________________________________________________________________
const char *DmtpcLensCorrection::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DmtpcLensCorrection*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DmtpcLensCorrection::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DmtpcLensCorrection*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void DmtpcLensCorrection::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DmtpcLensCorrection*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *DmtpcLensCorrection::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DmtpcLensCorrection*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void DmtpcLensCorrection::Streamer(TBuffer &R__b)
{
   // Stream an object of class DmtpcLensCorrection.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> _order;
      delete [] _coeffs;
      _coeffs = new double[_order+1];
      R__b.ReadFastArray(_coeffs,_order+1);
      R__b.CheckByteCount(R__s, R__c, DmtpcLensCorrection::IsA());
   } else {
      R__c = R__b.WriteVersion(DmtpcLensCorrection::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << _order;
      R__b.WriteFastArray(_coeffs,_order+1);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void DmtpcLensCorrection::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class DmtpcLensCorrection.
      TClass *R__cl = ::DmtpcLensCorrection::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_order", &_order);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_coeffs", &_coeffs);
      TNamed::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DmtpcLensCorrection(void *p) {
      return  p ? new(p) ::DmtpcLensCorrection : new ::DmtpcLensCorrection;
   }
   static void *newArray_DmtpcLensCorrection(Long_t nElements, void *p) {
      return p ? new(p) ::DmtpcLensCorrection[nElements] : new ::DmtpcLensCorrection[nElements];
   }
   // Wrapper around operator delete
   static void delete_DmtpcLensCorrection(void *p) {
      delete ((::DmtpcLensCorrection*)p);
   }
   static void deleteArray_DmtpcLensCorrection(void *p) {
      delete [] ((::DmtpcLensCorrection*)p);
   }
   static void destruct_DmtpcLensCorrection(void *p) {
      typedef ::DmtpcLensCorrection current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_DmtpcLensCorrection(TBuffer &buf, void *obj) {
      ((::DmtpcLensCorrection*)obj)->::DmtpcLensCorrection::Streamer(buf);
   }
} // end of namespace ROOT for class ::DmtpcLensCorrection

/********************************************************
* DmtpcLensCorrectionCint.cc
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

extern "C" void G__cpp_reset_tagtableDmtpcLensCorrectionCint();

extern "C" void G__set_cpp_environmentDmtpcLensCorrectionCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("DmtpcLensCorrection.hh");
  G__cpp_reset_tagtableDmtpcLensCorrectionCint();
}
#include <new>
extern "C" int G__cpp_dllrevDmtpcLensCorrectionCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* DmtpcLensCorrection */
static int G__DmtpcLensCorrectionCint_178_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmtpcLensCorrection* p = NULL;
   char* gvp = (char*) G__getgvp();
   switch (libp->paran) {
   case 3:
     //m: 3
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new DmtpcLensCorrection(
(const char*) G__int(libp->para[0]), (unsigned int) G__int(libp->para[1])
, (double*) G__int(libp->para[2]));
     } else {
       p = new((void*) gvp) DmtpcLensCorrection(
(const char*) G__int(libp->para[0]), (unsigned int) G__int(libp->para[1])
, (double*) G__int(libp->para[2]));
     }
     break;
   case 2:
     //m: 2
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new DmtpcLensCorrection((const char*) G__int(libp->para[0]), (unsigned int) G__int(libp->para[1]));
     } else {
       p = new((void*) gvp) DmtpcLensCorrection((const char*) G__int(libp->para[0]), (unsigned int) G__int(libp->para[1]));
     }
     break;
   case 1:
     //m: 1
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new DmtpcLensCorrection((const char*) G__int(libp->para[0]));
     } else {
       p = new((void*) gvp) DmtpcLensCorrection((const char*) G__int(libp->para[0]));
     }
     break;
   case 0:
     int n = G__getaryconstruct();
     if (n) {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new DmtpcLensCorrection[n];
       } else {
         p = new((void*) gvp) DmtpcLensCorrection[n];
       }
     } else {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new DmtpcLensCorrection;
       } else {
         p = new((void*) gvp) DmtpcLensCorrection;
       }
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcLensCorrection*) G__getstructoffset())->setParameter((unsigned int) G__int(libp->para[0]), (double) G__double(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcLensCorrection*) G__getstructoffset())->setParameters((double*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 5:
      G__letint(result7, 85, (long) ((const DmtpcLensCorrection*) G__getstructoffset())->correctDistortion((TH2*) G__int(libp->para[0]), (TH2*) G__int(libp->para[1])
, (const char*) G__int(libp->para[2]), (double*) G__int(libp->para[3])
, (bool) G__int(libp->para[4])));
      break;
   case 4:
      G__letint(result7, 85, (long) ((const DmtpcLensCorrection*) G__getstructoffset())->correctDistortion((TH2*) G__int(libp->para[0]), (TH2*) G__int(libp->para[1])
, (const char*) G__int(libp->para[2]), (double*) G__int(libp->para[3])));
      break;
   case 3:
      G__letint(result7, 85, (long) ((const DmtpcLensCorrection*) G__getstructoffset())->correctDistortion((TH2*) G__int(libp->para[0]), (TH2*) G__int(libp->para[1])
, (const char*) G__int(libp->para[2])));
      break;
   case 2:
      G__letint(result7, 85, (long) ((const DmtpcLensCorrection*) G__getstructoffset())->correctDistortion((TH2*) G__int(libp->para[0]), (TH2*) G__int(libp->para[1])));
      break;
   case 1:
      G__letint(result7, 85, (long) ((const DmtpcLensCorrection*) G__getstructoffset())->correctDistortion((TH2*) G__int(libp->para[0])));
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((const DmtpcLensCorrection*) G__getstructoffset())->distortCoords((double*) G__int(libp->para[0]), (double*) G__int(libp->para[1])
, (double*) G__int(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((const DmtpcLensCorrection*) G__getstructoffset())->unDistortCoords((double*) G__int(libp->para[0]), (double*) G__int(libp->para[1])
, (double*) G__int(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) DmtpcLensCorrection::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmtpcLensCorrection::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) DmtpcLensCorrection::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      DmtpcLensCorrection::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcLensCorrection*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmtpcLensCorrection::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) DmtpcLensCorrection::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmtpcLensCorrection::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcLensCorrectionCint_178_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) DmtpcLensCorrection::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__DmtpcLensCorrectionCint_178_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   DmtpcLensCorrection* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new DmtpcLensCorrection(*(DmtpcLensCorrection*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef DmtpcLensCorrection G__TDmtpcLensCorrection;
static int G__DmtpcLensCorrectionCint_178_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (DmtpcLensCorrection*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((DmtpcLensCorrection*) (soff+(sizeof(DmtpcLensCorrection)*i)))->~G__TDmtpcLensCorrection();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (DmtpcLensCorrection*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((DmtpcLensCorrection*) (soff))->~G__TDmtpcLensCorrection();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__DmtpcLensCorrectionCint_178_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmtpcLensCorrection* dest = (DmtpcLensCorrection*) G__getstructoffset();
   *dest = *(DmtpcLensCorrection*) libp->para[0].ref;
   const DmtpcLensCorrection& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* DmtpcLensCorrection */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncDmtpcLensCorrectionCint {
 public:
  G__Sizep2memfuncDmtpcLensCorrectionCint(): p(&G__Sizep2memfuncDmtpcLensCorrectionCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncDmtpcLensCorrectionCint::*p)();
};

size_t G__get_sizep2memfuncDmtpcLensCorrectionCint()
{
  G__Sizep2memfuncDmtpcLensCorrectionCint a;
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
extern "C" void G__cpp_setup_inheritanceDmtpcLensCorrectionCint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection))) {
     DmtpcLensCorrection *G__Lderived;
     G__Lderived=(DmtpcLensCorrection*)0x1000;
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection),G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection),G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableDmtpcLensCorrectionCint() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* DmtpcLensCorrection */
static void G__setup_memvarDmtpcLensCorrection(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection));
   { DmtpcLensCorrection *p; p=(DmtpcLensCorrection*)0x1000; if (p) { }
   G__memvar_setup((void*)0,104,0,0,-1,-1,-1,4,"_order=",0,(char*)NULL);
   G__memvar_setup((void*)0,68,0,0,-1,-1,-1,4,"_coeffs=",0,"[_order+1]; ");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarDmtpcLensCorrectionCint() {
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
static void G__setup_memfuncDmtpcLensCorrection(void) {
   /* DmtpcLensCorrection */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection));
   G__memfunc_setup("DmtpcLensCorrection",1954,G__DmtpcLensCorrectionCint_178_0_1, 105, G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection), -1, 0, 3, 1, 1, 0, 
"C - - 10 '\"lenscorr\"' name h - - 0 '2' - "
"D - - 0 '0' polyn", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setParameter",1261,G__DmtpcLensCorrectionCint_178_0_2, 121, -1, -1, 0, 2, 1, 1, 0, 
"h - - 0 - - d - - 0 - val", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setParameters",1376,G__DmtpcLensCorrectionCint_178_0_3, 121, -1, -1, 0, 1, 1, 1, 0, "D - - 10 - vals", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("correctDistortion",1825,G__DmtpcLensCorrectionCint_178_0_4, 85, G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_TH2), -1, 0, 5, 1, 1, 8, 
"U 'TH2' - 10 - in U 'TH2' - 0 '0' out "
"C - - 10 '\"bicubic\"' interpolation D - - 10 '0' camera_center "
"g - - 0 'false' distort", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("distortCoords",1395,G__DmtpcLensCorrectionCint_178_0_5, 105, -1, -1, 0, 3, 1, 1, 8, 
"D - - 10 - center D - - 10 - in "
"D - - 0 - out", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("unDistortCoords",1590,G__DmtpcLensCorrectionCint_178_0_6, 105, -1, -1, 0, 3, 1, 1, 8, 
"D - - 10 - center D - - 10 - in "
"D - - 0 - out", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__DmtpcLensCorrectionCint_178_0_7, 85, G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&DmtpcLensCorrection::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__DmtpcLensCorrectionCint_178_0_8, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmtpcLensCorrection::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__DmtpcLensCorrectionCint_178_0_9, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&DmtpcLensCorrection::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__DmtpcLensCorrectionCint_178_0_10, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&DmtpcLensCorrection::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__DmtpcLensCorrectionCint_178_0_14, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__DmtpcLensCorrectionCint_178_0_15, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmtpcLensCorrection::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__DmtpcLensCorrectionCint_178_0_16, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&DmtpcLensCorrection::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__DmtpcLensCorrectionCint_178_0_17, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmtpcLensCorrection::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__DmtpcLensCorrectionCint_178_0_18, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&DmtpcLensCorrection::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("DmtpcLensCorrection", 1954, G__DmtpcLensCorrectionCint_178_0_19, (int) ('i'), G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection), -1, 0, 1, 1, 1, 0, "u 'DmtpcLensCorrection' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~DmtpcLensCorrection", 2080, G__DmtpcLensCorrectionCint_178_0_20, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__DmtpcLensCorrectionCint_178_0_21, (int) ('u'), G__get_linked_tagnum(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection), -1, 1, 1, 1, 1, 0, "u 'DmtpcLensCorrection' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncDmtpcLensCorrectionCint() {
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
extern "C" void G__cpp_setup_globalDmtpcLensCorrectionCint() {
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

extern "C" void G__cpp_setup_funcDmtpcLensCorrectionCint() {
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
G__linked_taginfo G__DmtpcLensCorrectionCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_TNamed = { "TNamed" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_TH2 = { "TH2" , 99 , -1 };
G__linked_taginfo G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection = { "DmtpcLensCorrection" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableDmtpcLensCorrectionCint() {
  G__DmtpcLensCorrectionCintLN_TClass.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_TBuffer.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_TMemberInspector.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_TObject.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_TNamed.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_TH2.tagnum = -1 ;
  G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableDmtpcLensCorrectionCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_TObject);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_TNamed);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_TH2);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__DmtpcLensCorrectionCintLN_DmtpcLensCorrection),sizeof(DmtpcLensCorrection),-1,62720,(char*)NULL,G__setup_memvarDmtpcLensCorrection,G__setup_memfuncDmtpcLensCorrection);
}
extern "C" void G__cpp_setupDmtpcLensCorrectionCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupDmtpcLensCorrectionCint()");
  G__set_cpp_environmentDmtpcLensCorrectionCint();
  G__cpp_setup_tagtableDmtpcLensCorrectionCint();

  G__cpp_setup_inheritanceDmtpcLensCorrectionCint();

  G__cpp_setup_typetableDmtpcLensCorrectionCint();

  G__cpp_setup_memvarDmtpcLensCorrectionCint();

  G__cpp_setup_memfuncDmtpcLensCorrectionCint();
  G__cpp_setup_globalDmtpcLensCorrectionCint();
  G__cpp_setup_funcDmtpcLensCorrectionCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncDmtpcLensCorrectionCint();
  return;
}
class G__cpp_setup_initDmtpcLensCorrectionCint {
  public:
    G__cpp_setup_initDmtpcLensCorrectionCint() { G__add_setup_func("DmtpcLensCorrectionCint",(G__incsetup)(&G__cpp_setupDmtpcLensCorrectionCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initDmtpcLensCorrectionCint() { G__remove_setup_func("DmtpcLensCorrectionCint"); }
};
G__cpp_setup_initDmtpcLensCorrectionCint G__cpp_setup_initializerDmtpcLensCorrectionCint;
