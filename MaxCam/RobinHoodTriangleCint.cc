//
// File generated by rootcint at Mon Jun  6 11:03:36 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME RobinHoodTriangleCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "RobinHoodTriangleCint.h"

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
   void RobinHoodTriangle_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_RobinHoodTriangle(void *p);
   static void deleteArray_RobinHoodTriangle(void *p);
   static void destruct_RobinHoodTriangle(void *p);
   static void streamer_RobinHoodTriangle(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RobinHoodTriangle*)
   {
      ::RobinHoodTriangle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RobinHoodTriangle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RobinHoodTriangle", ::RobinHoodTriangle::Class_Version(), "./RobinHoodTriangle.hh", 9,
                  typeid(::RobinHoodTriangle), DefineBehavior(ptr, ptr),
                  &::RobinHoodTriangle::Dictionary, isa_proxy, 0,
                  sizeof(::RobinHoodTriangle) );
      instance.SetDelete(&delete_RobinHoodTriangle);
      instance.SetDeleteArray(&deleteArray_RobinHoodTriangle);
      instance.SetDestructor(&destruct_RobinHoodTriangle);
      instance.SetStreamerFunc(&streamer_RobinHoodTriangle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RobinHoodTriangle*)
   {
      return GenerateInitInstanceLocal((::RobinHoodTriangle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RobinHoodTriangle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *RobinHoodTriangle::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *RobinHoodTriangle::Class_Name()
{
   return "RobinHoodTriangle";
}

//______________________________________________________________________________
const char *RobinHoodTriangle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RobinHoodTriangle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RobinHoodTriangle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RobinHoodTriangle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void RobinHoodTriangle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RobinHoodTriangle*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *RobinHoodTriangle::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RobinHoodTriangle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void RobinHoodTriangle::Streamer(TBuffer &R__b)
{
   // Stream an object of class RobinHoodTriangle.

   TObject::Streamer(R__b);
}

//______________________________________________________________________________
void RobinHoodTriangle::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class RobinHoodTriangle.
      TClass *R__cl = ::RobinHoodTriangle::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_v[3]", &_v);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_charge", &_charge);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_center", &_center);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "_area", &_area);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RobinHoodTriangle(void *p) {
      delete ((::RobinHoodTriangle*)p);
   }
   static void deleteArray_RobinHoodTriangle(void *p) {
      delete [] ((::RobinHoodTriangle*)p);
   }
   static void destruct_RobinHoodTriangle(void *p) {
      typedef ::RobinHoodTriangle current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RobinHoodTriangle(TBuffer &buf, void *obj) {
      ((::RobinHoodTriangle*)obj)->::RobinHoodTriangle::Streamer(buf);
   }
} // end of namespace ROOT for class ::RobinHoodTriangle

/********************************************************
* RobinHoodTriangleCint.cc
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

extern "C" void G__cpp_reset_tagtableRobinHoodTriangleCint();

extern "C" void G__set_cpp_environmentRobinHoodTriangleCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("RobinHoodTriangle.hh");
  G__cpp_reset_tagtableRobinHoodTriangleCint();
}
#include <new>
extern "C" int G__cpp_dllrevRobinHoodTriangleCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* RobinHoodTriangle */
static int G__RobinHoodTriangleCint_163_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   RobinHoodTriangle* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 4
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new RobinHoodTriangle(
(TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (TVector3*) G__int(libp->para[2]), (double) G__double(libp->para[3]));
   } else {
     p = new((void*) gvp) RobinHoodTriangle(
(TVector3*) G__int(libp->para[0]), (TVector3*) G__int(libp->para[1])
, (TVector3*) G__int(libp->para[2]), (double) G__double(libp->para[3]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   RobinHoodTriangle* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new RobinHoodTriangle(*(RobinHoodTriangle*) libp->para[0].ref);
   } else {
     p = new((void*) gvp) RobinHoodTriangle(*(RobinHoodTriangle*) libp->para[0].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((RobinHoodTriangle*) G__getstructoffset())->getCharge());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((RobinHoodTriangle*) G__getstructoffset())->setCharge((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((RobinHoodTriangle*) G__getstructoffset())->getArea());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((RobinHoodTriangle*) G__getstructoffset())->getVertex((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) RobinHoodTriangle::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) RobinHoodTriangle::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) RobinHoodTriangle::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      RobinHoodTriangle::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((RobinHoodTriangle*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) RobinHoodTriangle::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) RobinHoodTriangle::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) RobinHoodTriangle::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__RobinHoodTriangleCint_163_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) RobinHoodTriangle::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef RobinHoodTriangle G__TRobinHoodTriangle;
static int G__RobinHoodTriangleCint_163_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (RobinHoodTriangle*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((RobinHoodTriangle*) (soff+(sizeof(RobinHoodTriangle)*i)))->~G__TRobinHoodTriangle();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (RobinHoodTriangle*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((RobinHoodTriangle*) (soff))->~G__TRobinHoodTriangle();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__RobinHoodTriangleCint_163_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   RobinHoodTriangle* dest = (RobinHoodTriangle*) G__getstructoffset();
   *dest = *(RobinHoodTriangle*) libp->para[0].ref;
   const RobinHoodTriangle& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* RobinHoodTriangle */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncRobinHoodTriangleCint {
 public:
  G__Sizep2memfuncRobinHoodTriangleCint(): p(&G__Sizep2memfuncRobinHoodTriangleCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncRobinHoodTriangleCint::*p)();
};

size_t G__get_sizep2memfuncRobinHoodTriangleCint()
{
  G__Sizep2memfuncRobinHoodTriangleCint a;
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
extern "C" void G__cpp_setup_inheritanceRobinHoodTriangleCint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle))) {
     RobinHoodTriangle *G__Lderived;
     G__Lderived=(RobinHoodTriangle*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle),G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableRobinHoodTriangleCint() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* RobinHoodTriangle */
static void G__setup_memvarRobinHoodTriangle(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle));
   { RobinHoodTriangle *p; p=(RobinHoodTriangle*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_TVector3),-1,-1,4,"_v[3]=",0,"corners of triangle");
   G__memvar_setup((void*)0,100,0,0,-1,-1,-1,4,"_charge=",0,"charge deposited on triangle");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_TVector3),-1,-1,4,"_center=",0,"center of triangle");
   G__memvar_setup((void*)0,100,0,0,-1,-1,-1,4,"_area=",0,"area of triangle");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarRobinHoodTriangleCint() {
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
static void G__setup_memfuncRobinHoodTriangle(void) {
   /* RobinHoodTriangle */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle));
   G__memfunc_setup("RobinHoodTriangle",1722,G__RobinHoodTriangleCint_163_0_1, 105, G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle), -1, 0, 4, 1, 1, 0, 
"U 'TVector3' - 0 - a U 'TVector3' - 0 - b "
"U 'TVector3' - 0 - c d - - 0 - charge", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("RobinHoodTriangle",1722,G__RobinHoodTriangleCint_163_0_2, 105, G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle), -1, 0, 1, 1, 1, 0, "u 'RobinHoodTriangle' - 1 - t", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getCharge",906,G__RobinHoodTriangleCint_163_0_3, 100, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setCharge",918,G__RobinHoodTriangleCint_163_0_4, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - charge", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getArea",697,G__RobinHoodTriangleCint_163_0_5, 100, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getVertex",958,G__RobinHoodTriangleCint_163_0_6, 85, G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_TVector3), -1, 0, 1, 1, 1, 0, "i - - 0 - i", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__RobinHoodTriangleCint_163_0_7, 85, G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&RobinHoodTriangle::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__RobinHoodTriangleCint_163_0_8, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&RobinHoodTriangle::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__RobinHoodTriangleCint_163_0_9, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&RobinHoodTriangle::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__RobinHoodTriangleCint_163_0_10, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&RobinHoodTriangle::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__RobinHoodTriangleCint_163_0_14, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__RobinHoodTriangleCint_163_0_15, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&RobinHoodTriangle::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__RobinHoodTriangleCint_163_0_16, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&RobinHoodTriangle::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__RobinHoodTriangleCint_163_0_17, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&RobinHoodTriangle::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__RobinHoodTriangleCint_163_0_18, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&RobinHoodTriangle::DeclFileLine) ), 0);
   // automatic destructor
   G__memfunc_setup("~RobinHoodTriangle", 1848, G__RobinHoodTriangleCint_163_0_19, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__RobinHoodTriangleCint_163_0_20, (int) ('u'), G__get_linked_tagnum(&G__RobinHoodTriangleCintLN_RobinHoodTriangle), -1, 1, 1, 1, 1, 0, "u 'RobinHoodTriangle' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncRobinHoodTriangleCint() {
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
extern "C" void G__cpp_setup_globalRobinHoodTriangleCint() {
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

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcRobinHoodTriangleCint() {
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
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__RobinHoodTriangleCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_TVector3 = { "TVector3" , 99 , -1 };
G__linked_taginfo G__RobinHoodTriangleCintLN_RobinHoodTriangle = { "RobinHoodTriangle" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableRobinHoodTriangleCint() {
  G__RobinHoodTriangleCintLN_TClass.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_TBuffer.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_TMemberInspector.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_TObject.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_TVector3.tagnum = -1 ;
  G__RobinHoodTriangleCintLN_RobinHoodTriangle.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableRobinHoodTriangleCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_TObject);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_TVector3);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__RobinHoodTriangleCintLN_RobinHoodTriangle),sizeof(RobinHoodTriangle),-1,62976,(char*)NULL,G__setup_memvarRobinHoodTriangle,G__setup_memfuncRobinHoodTriangle);
}
extern "C" void G__cpp_setupRobinHoodTriangleCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupRobinHoodTriangleCint()");
  G__set_cpp_environmentRobinHoodTriangleCint();
  G__cpp_setup_tagtableRobinHoodTriangleCint();

  G__cpp_setup_inheritanceRobinHoodTriangleCint();

  G__cpp_setup_typetableRobinHoodTriangleCint();

  G__cpp_setup_memvarRobinHoodTriangleCint();

  G__cpp_setup_memfuncRobinHoodTriangleCint();
  G__cpp_setup_globalRobinHoodTriangleCint();
  G__cpp_setup_funcRobinHoodTriangleCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncRobinHoodTriangleCint();
  return;
}
class G__cpp_setup_initRobinHoodTriangleCint {
  public:
    G__cpp_setup_initRobinHoodTriangleCint() { G__add_setup_func("RobinHoodTriangleCint",(G__incsetup)(&G__cpp_setupRobinHoodTriangleCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initRobinHoodTriangleCint() { G__remove_setup_func("RobinHoodTriangleCint"); }
};
G__cpp_setup_initRobinHoodTriangleCint G__cpp_setup_initializerRobinHoodTriangleCint;
