//
// File generated by rootcint at Tue Aug 16 18:40:35 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME RobinHoodTriangleMakerCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "RobinHoodTriangleMakerCint.h"

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
   void RobinHoodTriangleMaker_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_RobinHoodTriangleMaker(void *p = 0);
   static void *newArray_RobinHoodTriangleMaker(Long_t size, void *p);
   static void delete_RobinHoodTriangleMaker(void *p);
   static void deleteArray_RobinHoodTriangleMaker(void *p);
   static void destruct_RobinHoodTriangleMaker(void *p);
   static void streamer_RobinHoodTriangleMaker(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RobinHoodTriangleMaker*)
   {
      ::RobinHoodTriangleMaker *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RobinHoodTriangleMaker >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RobinHoodTriangleMaker", ::RobinHoodTriangleMaker::Class_Version(), "./RobinHoodTriangleMaker.hh", 11,
                  typeid(::RobinHoodTriangleMaker), DefineBehavior(ptr, ptr),
                  &::RobinHoodTriangleMaker::Dictionary, isa_proxy, 0,
                  sizeof(::RobinHoodTriangleMaker) );
      instance.SetNew(&new_RobinHoodTriangleMaker);
      instance.SetNewArray(&newArray_RobinHoodTriangleMaker);
      instance.SetDelete(&delete_RobinHoodTriangleMaker);
      instance.SetDeleteArray(&deleteArray_RobinHoodTriangleMaker);
      instance.SetDestructor(&destruct_RobinHoodTriangleMaker);
      instance.SetStreamerFunc(&streamer_RobinHoodTriangleMaker);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RobinHoodTriangleMaker*)
   {
      return GenerateInitInstanceLocal((::RobinHoodTriangleMaker*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RobinHoodTriangleMaker*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *RobinHoodTriangleMaker::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *RobinHoodTriangleMaker::Class_Name()
{
   return "RobinHoodTriangleMaker";
}

//______________________________________________________________________________
const char *RobinHoodTriangleMaker::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RobinHoodTriangleMaker*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RobinHoodTriangleMaker::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RobinHoodTriangleMaker*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void RobinHoodTriangleMaker::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RobinHoodTriangleMaker*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *RobinHoodTriangleMaker::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RobinHoodTriangleMaker*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void RobinHoodTriangleMaker::Streamer(TBuffer &R__b)
{
   // Stream an object of class RobinHoodTriangleMaker.

   TObject::Streamer(R__b);
}

//______________________________________________________________________________
void RobinHoodTriangleMaker::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class RobinHoodTriangleMaker.
      TClass *R__cl = ::RobinHoodTriangleMaker::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_list", &_list);
      R__insp.InspectMember(_list, "_list.");
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RobinHoodTriangleMaker(void *p) {
      return  p ? new(p) ::RobinHoodTriangleMaker : new ::RobinHoodTriangleMaker;
   }
   static void *newArray_RobinHoodTriangleMaker(Long_t nElements, void *p) {
      return p ? new(p) ::RobinHoodTriangleMaker[nElements] : new ::RobinHoodTriangleMaker[nElements];
   }
   // Wrapper around operator delete
   static void delete_RobinHoodTriangleMaker(void *p) {
      delete ((::RobinHoodTriangleMaker*)p);
   }
   static void deleteArray_RobinHoodTriangleMaker(void *p) {
      delete [] ((::RobinHoodTriangleMaker*)p);
   }
   static void destruct_RobinHoodTriangleMaker(void *p) {
      typedef ::RobinHoodTriangleMaker current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RobinHoodTriangleMaker(TBuffer &buf, void *obj) {
      ((::RobinHoodTriangleMaker*)obj)->::RobinHoodTriangleMaker::Streamer(buf);
   }
} // end of namespace ROOT for class ::RobinHoodTriangleMaker

/********************************************************
* RobinHoodTriangleMakerCint.cc
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

extern "C" void G__cpp_reset_tagtableRobinHoodTriangleMakerCint();

extern "C" void G__set_cpp_environmentRobinHoodTriangleMakerCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("RobinHoodTriangleMaker.hh");
  G__cpp_reset_tagtableRobinHoodTriangleMakerCint();
}
#include <new>
extern "C" int G__cpp_dllrevRobinHoodTriangleMakerCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* RobinHoodTriangleMaker */
static int G__RobinHoodTriangleMakerCint_188_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   RobinHoodTriangleMaker* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new RobinHoodTriangleMaker[n];
     } else {
       p = new((void*) gvp) RobinHoodTriangleMaker[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new RobinHoodTriangleMaker;
     } else {
       p = new((void*) gvp) RobinHoodTriangleMaker;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 7:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->cylinder(
(TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]), (int) G__int(libp->para[5])
, *((TString*) G__int(libp->para[6])));
      G__setnull(result7);
      break;
   case 6:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->cylinder((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]), (int) G__int(libp->para[5]));
      G__setnull(result7);
      break;
   case 5:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->cylinder((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]));
      G__setnull(result7);
      break;
   case 4:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->cylinder((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3]));
      G__setnull(result7);
      break;
   case 3:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->cylinder((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]));
      G__setnull(result7);
      break;
   case 2:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->cylinder((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1]));
      G__setnull(result7);
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 5:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->sphere((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (int) G__int(libp->para[4]));
      G__setnull(result7);
      break;
   case 4:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->sphere((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3]));
      G__setnull(result7);
      break;
   case 3:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->sphere((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]));
      G__setnull(result7);
      break;
   case 2:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->sphere((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1]));
      G__setnull(result7);
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 5:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->disk((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]), (int) G__int(libp->para[3])
, *((TString*) G__int(libp->para[4])));
      G__setnull(result7);
      break;
   case 4:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->disk((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]), (int) G__int(libp->para[3]));
      G__setnull(result7);
      break;
   case 3:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->disk((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (double) G__double(libp->para[2]));
      G__setnull(result7);
      break;
   case 2:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->disk((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1]));
      G__setnull(result7);
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 5:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->rectangle((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (TVector3*) G__int(libp->para[2]), (int) G__int(libp->para[3])
, *((TString*) G__int(libp->para[4])));
      G__setnull(result7);
      break;
   case 4:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->rectangle((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (TVector3*) G__int(libp->para[2]), (int) G__int(libp->para[3]));
      G__setnull(result7);
      break;
   case 3:
      ((RobinHoodTriangleMaker*) G__getstructoffset())->rectangle((TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (TVector3*) G__int(libp->para[2]));
      G__setnull(result7);
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((RobinHoodTriangleMaker*) G__getstructoffset())->dumpTriangles((const char*) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((RobinHoodTriangleMaker*) G__getstructoffset())->addConductorToRhc((const char*) G__int(libp->para[0]), (double) G__double(libp->para[1])
, (const char*) G__int(libp->para[2]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) RobinHoodTriangleMaker::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) RobinHoodTriangleMaker::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) RobinHoodTriangleMaker::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      RobinHoodTriangleMaker::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((RobinHoodTriangleMaker*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) RobinHoodTriangleMaker::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) RobinHoodTriangleMaker::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) RobinHoodTriangleMaker::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleMakerCint_188_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) RobinHoodTriangleMaker::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef RobinHoodTriangleMaker G__TRobinHoodTriangleMaker;
static int G__RobinHoodTriangleMakerCint_188_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (RobinHoodTriangleMaker*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((RobinHoodTriangleMaker*) (soff+(sizeof(RobinHoodTriangleMaker)*i)))->~G__TRobinHoodTriangleMaker();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (RobinHoodTriangleMaker*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((RobinHoodTriangleMaker*) (soff))->~G__TRobinHoodTriangleMaker();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* RobinHoodTriangleMaker */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncRobinHoodTriangleMakerCint {
 public:
  G__Sizep2memfuncRobinHoodTriangleMakerCint(): p(&G__Sizep2memfuncRobinHoodTriangleMakerCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncRobinHoodTriangleMakerCint::*p)();
};

size_t G__get_sizep2memfuncRobinHoodTriangleMakerCint()
{
  G__Sizep2memfuncRobinHoodTriangleMakerCint a;
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
extern "C" void G__cpp_setup_inheritanceRobinHoodTriangleMakerCint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker))) {
     RobinHoodTriangleMaker *G__Lderived;
     G__Lderived=(RobinHoodTriangleMaker*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker),G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableRobinHoodTriangleMakerCint() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* RobinHoodTriangleMaker */
static void G__setup_memvarRobinHoodTriangleMaker(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker));
   { RobinHoodTriangleMaker *p; p=(RobinHoodTriangleMaker*)0x1000; if (p) { }
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_TList),-1,-1,4,"_list=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarRobinHoodTriangleMakerCint() {
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
static void G__setup_memfuncRobinHoodTriangleMaker(void) {
   /* RobinHoodTriangleMaker */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker));
   G__memfunc_setup("RobinHoodTriangleMaker",2218,G__RobinHoodTriangleMakerCint_188_0_1, 105, G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("cylinder",858,G__RobinHoodTriangleMakerCint_188_0_2, 121, -1, -1, 0, 7, 1, 1, 0, 
"U 'TVector3' - 0 - center U 'TVector3' - 0 - norm "
"d - - 0 '-1.0' ID d - - 0 '-1.0' OD "
"d - - 0 '-1.0' thickness i - - 0 '0' nTriangles "
"u 'TString' - 0 '\"\"' opt", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("sphere",647,G__RobinHoodTriangleMakerCint_188_0_3, 121, -1, -1, 0, 5, 1, 1, 0, 
"U 'TVector3' - 0 - center U 'TVector3' - 0 - norm "
"d - - 0 '-1.0' radius d - - 0 '-1.0' height "
"i - - 0 '0' nTriangles", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("disk",427,G__RobinHoodTriangleMakerCint_188_0_4, 121, -1, -1, 0, 5, 1, 1, 0, 
"U 'TVector3' - 0 - center U 'TVector3' - 0 - norm "
"d - - 0 '0' Diameter i - - 0 '0' nTriangles "
"u 'TString' - 0 '\"\"' opt", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("rectangle",949,G__RobinHoodTriangleMakerCint_188_0_5, 121, -1, -1, 0, 5, 1, 1, 0, 
"U 'TVector3' - 0 - center U 'TVector3' - 0 - side1 "
"U 'TVector3' - 0 - side2 i - - 0 '0' nTriangles "
"u 'TString' - 0 '\"\"' opt", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("dumpTriangles",1375,G__RobinHoodTriangleMakerCint_188_0_6, 121, -1, -1, 0, 1, 1, 1, 0, "C - - 10 - geoFile", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("addConductorToRhc",1722,G__RobinHoodTriangleMakerCint_188_0_7, 121, -1, -1, 0, 3, 1, 1, 0, 
"C - - 10 - rhcFile d - - 0 - volts "
"C - - 10 - geoFile", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__RobinHoodTriangleMakerCint_188_0_8, 85, G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&RobinHoodTriangleMaker::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__RobinHoodTriangleMakerCint_188_0_9, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&RobinHoodTriangleMaker::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__RobinHoodTriangleMakerCint_188_0_10, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&RobinHoodTriangleMaker::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__RobinHoodTriangleMakerCint_188_0_11, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&RobinHoodTriangleMaker::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__RobinHoodTriangleMakerCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__RobinHoodTriangleMakerCint_188_0_15, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__RobinHoodTriangleMakerCint_188_0_16, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&RobinHoodTriangleMaker::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__RobinHoodTriangleMakerCint_188_0_17, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&RobinHoodTriangleMaker::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__RobinHoodTriangleMakerCint_188_0_18, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&RobinHoodTriangleMaker::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__RobinHoodTriangleMakerCint_188_0_19, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&RobinHoodTriangleMaker::DeclFileLine) ), 0);
   // automatic destructor
   G__memfunc_setup("~RobinHoodTriangleMaker", 2344, G__RobinHoodTriangleMakerCint_188_0_20, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncRobinHoodTriangleMakerCint() {
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
extern "C" void G__cpp_setup_globalRobinHoodTriangleMakerCint() {
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

extern "C" void G__cpp_setup_funcRobinHoodTriangleMakerCint() {
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
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_TString = { "TString" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_TList = { "TList" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_TVector3 = { "TVector3" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker = { "RobinHoodTriangleMaker" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableRobinHoodTriangleMakerCint() {
  G__RobinHoodTriangleMakerCintLN_TClass.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_TBuffer.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_TMemberInspector.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_TObject.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_TString.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_TList.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_TVector3.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableRobinHoodTriangleMakerCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_TObject);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_TString);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_TList);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_TVector3);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__RobinHoodTriangleMakerCintLN_RobinHoodTriangleMaker),sizeof(RobinHoodTriangleMaker),-1,29952,(char*)NULL,G__setup_memvarRobinHoodTriangleMaker,G__setup_memfuncRobinHoodTriangleMaker);
}
extern "C" void G__cpp_setupRobinHoodTriangleMakerCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupRobinHoodTriangleMakerCint()");
  G__set_cpp_environmentRobinHoodTriangleMakerCint();
  G__cpp_setup_tagtableRobinHoodTriangleMakerCint();

  G__cpp_setup_inheritanceRobinHoodTriangleMakerCint();

  G__cpp_setup_typetableRobinHoodTriangleMakerCint();

  G__cpp_setup_memvarRobinHoodTriangleMakerCint();

  G__cpp_setup_memfuncRobinHoodTriangleMakerCint();
  G__cpp_setup_globalRobinHoodTriangleMakerCint();
  G__cpp_setup_funcRobinHoodTriangleMakerCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncRobinHoodTriangleMakerCint();
  return;
}
class G__cpp_setup_initRobinHoodTriangleMakerCint {
  public:
    G__cpp_setup_initRobinHoodTriangleMakerCint() { G__add_setup_func("RobinHoodTriangleMakerCint",(G__incsetup)(&G__cpp_setupRobinHoodTriangleMakerCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initRobinHoodTriangleMakerCint() { G__remove_setup_func("RobinHoodTriangleMakerCint"); }
};
G__cpp_setup_initRobinHoodTriangleMakerCint G__cpp_setup_initializerRobinHoodTriangleMakerCint;

