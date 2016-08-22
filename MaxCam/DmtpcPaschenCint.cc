//
// File generated by rootcint at Tue Aug 16 18:38:50 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME DmtpcPaschenCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "DmtpcPaschenCint.h"

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
   void DmtpcPaschen_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_DmtpcPaschen(void *p = 0);
   static void *newArray_DmtpcPaschen(Long_t size, void *p);
   static void delete_DmtpcPaschen(void *p);
   static void deleteArray_DmtpcPaschen(void *p);
   static void destruct_DmtpcPaschen(void *p);
   static void streamer_DmtpcPaschen(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DmtpcPaschen*)
   {
      ::DmtpcPaschen *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::DmtpcPaschen >(0);
      static ::ROOT::TGenericClassInfo 
         instance("DmtpcPaschen", ::DmtpcPaschen::Class_Version(), "./DmtpcPaschen.hh", 12,
                  typeid(::DmtpcPaschen), DefineBehavior(ptr, ptr),
                  &::DmtpcPaschen::Dictionary, isa_proxy, 0,
                  sizeof(::DmtpcPaschen) );
      instance.SetNew(&new_DmtpcPaschen);
      instance.SetNewArray(&newArray_DmtpcPaschen);
      instance.SetDelete(&delete_DmtpcPaschen);
      instance.SetDeleteArray(&deleteArray_DmtpcPaschen);
      instance.SetDestructor(&destruct_DmtpcPaschen);
      instance.SetStreamerFunc(&streamer_DmtpcPaschen);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DmtpcPaschen*)
   {
      return GenerateInitInstanceLocal((::DmtpcPaschen*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::DmtpcPaschen*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *DmtpcPaschen::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *DmtpcPaschen::Class_Name()
{
   return "DmtpcPaschen";
}

//______________________________________________________________________________
const char *DmtpcPaschen::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DmtpcPaschen*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int DmtpcPaschen::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::DmtpcPaschen*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void DmtpcPaschen::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DmtpcPaschen*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *DmtpcPaschen::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::DmtpcPaschen*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void DmtpcPaschen::Streamer(TBuffer &R__b)
{
   // Stream an object of class DmtpcPaschen.

   ::Error("DmtpcPaschen::Streamer", "version id <=0 in ClassDef, dummy Streamer() called"); if (R__b.IsReading()) { }
}

//______________________________________________________________________________
void DmtpcPaschen::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class DmtpcPaschen.
      TClass *R__cl = ::DmtpcPaschen::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_pas", &_pas);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_fname", &_fname);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_delim", &_delim);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_fit", &_fit);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_DmtpcPaschen(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::DmtpcPaschen : new ::DmtpcPaschen;
   }
   static void *newArray_DmtpcPaschen(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::DmtpcPaschen[nElements] : new ::DmtpcPaschen[nElements];
   }
   // Wrapper around operator delete
   static void delete_DmtpcPaschen(void *p) {
      delete ((::DmtpcPaschen*)p);
   }
   static void deleteArray_DmtpcPaschen(void *p) {
      delete [] ((::DmtpcPaschen*)p);
   }
   static void destruct_DmtpcPaschen(void *p) {
      typedef ::DmtpcPaschen current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_DmtpcPaschen(TBuffer &buf, void *obj) {
      ((::DmtpcPaschen*)obj)->::DmtpcPaschen::Streamer(buf);
   }
} // end of namespace ROOT for class ::DmtpcPaschen

/********************************************************
* DmtpcPaschenCint.cc
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

extern "C" void G__cpp_reset_tagtableDmtpcPaschenCint();

extern "C" void G__set_cpp_environmentDmtpcPaschenCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("DmtpcPaschen.hh");
  G__cpp_reset_tagtableDmtpcPaschenCint();
}
#include <new>
extern "C" int G__cpp_dllrevDmtpcPaschenCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* DmtpcPaschen */
static int G__DmtpcPaschenCint_207_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmtpcPaschen* p = NULL;
   char* gvp = (char*) G__getgvp();
   switch (libp->paran) {
   case 3:
     //m: 3
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new DmtpcPaschen(
(const char*) G__int(libp->para[0]), *((TString*) G__int(libp->para[1]))
, *((TString*) G__int(libp->para[2])));
     } else {
       p = new((void*) gvp) DmtpcPaschen(
(const char*) G__int(libp->para[0]), *((TString*) G__int(libp->para[1]))
, *((TString*) G__int(libp->para[2])));
     }
     break;
   case 2:
     //m: 2
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new DmtpcPaschen((const char*) G__int(libp->para[0]), *((TString*) G__int(libp->para[1])));
     } else {
       p = new((void*) gvp) DmtpcPaschen((const char*) G__int(libp->para[0]), *((TString*) G__int(libp->para[1])));
     }
     break;
   case 1:
     //m: 1
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new DmtpcPaschen((const char*) G__int(libp->para[0]));
     } else {
       p = new((void*) gvp) DmtpcPaschen((const char*) G__int(libp->para[0]));
     }
     break;
   case 0:
     int n = G__getaryconstruct();
     if (n) {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new DmtpcPaschen[n];
       } else {
         p = new((void*) gvp) DmtpcPaschen[n];
       }
     } else {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new DmtpcPaschen;
       } else {
         p = new((void*) gvp) DmtpcPaschen;
       }
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_DmtpcPaschen));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmtpcPaschen* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new DmtpcPaschen(*(DmtpcPaschen*) libp->para[0].ref);
   } else {
     p = new((void*) gvp) DmtpcPaschen(*(DmtpcPaschen*) libp->para[0].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_DmtpcPaschen));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         TString* pobj;
         TString xobj = ((DmtpcPaschen*) G__getstructoffset())->GetName();
         pobj = new TString(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((DmtpcPaschen*) G__getstructoffset())->getPaschenCurve());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcPaschen*) G__getstructoffset())->readPaschenData();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcPaschen*) G__getstructoffset())->makeFitFunction();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((DmtpcPaschen*) G__getstructoffset())->getFitFunction());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcPaschen*) G__getstructoffset())->testme();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) DmtpcPaschen::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmtpcPaschen::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) DmtpcPaschen::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      DmtpcPaschen::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((const DmtpcPaschen*) G__getstructoffset())->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcPaschen*) G__getstructoffset())->ShowMembers(*(TMemberInspector*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcPaschen*) G__getstructoffset())->Streamer(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((DmtpcPaschen*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmtpcPaschen::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) DmtpcPaschen::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) DmtpcPaschen::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcPaschenCint_207_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) DmtpcPaschen::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef DmtpcPaschen G__TDmtpcPaschen;
static int G__DmtpcPaschenCint_207_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 0
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (DmtpcPaschen*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((DmtpcPaschen*) (soff+(sizeof(DmtpcPaschen)*i)))->~G__TDmtpcPaschen();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (DmtpcPaschen*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((DmtpcPaschen*) (soff))->~G__TDmtpcPaschen();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__DmtpcPaschenCint_207_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   DmtpcPaschen* dest = (DmtpcPaschen*) G__getstructoffset();
   *dest = *(DmtpcPaschen*) libp->para[0].ref;
   const DmtpcPaschen& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* DmtpcPaschen */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncDmtpcPaschenCint {
 public:
  G__Sizep2memfuncDmtpcPaschenCint(): p(&G__Sizep2memfuncDmtpcPaschenCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncDmtpcPaschenCint::*p)();
};

size_t G__get_sizep2memfuncDmtpcPaschenCint()
{
  G__Sizep2memfuncDmtpcPaschenCint a;
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
extern "C" void G__cpp_setup_inheritanceDmtpcPaschenCint() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableDmtpcPaschenCint() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* DmtpcPaschen */
static void G__setup_memvarDmtpcPaschen(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__DmtpcPaschenCintLN_DmtpcPaschen));
   { DmtpcPaschen *p; p=(DmtpcPaschen*)0x1000; if (p) { }
   G__memvar_setup((void*)0,108,0,0,-1,-1,-1,4,"G__virtualinfo=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TGraph),-1,-1,4,"_pas=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TString),-1,-1,4,"_fname=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TString),-1,-1,4,"_delim=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TF1),-1,-1,4,"_fit=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarDmtpcPaschenCint() {
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
static void G__setup_memfuncDmtpcPaschen(void) {
   /* DmtpcPaschen */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__DmtpcPaschenCintLN_DmtpcPaschen));
   G__memfunc_setup("DmtpcPaschen",1210,G__DmtpcPaschenCint_207_0_1, 105, G__get_linked_tagnum(&G__DmtpcPaschenCintLN_DmtpcPaschen), -1, 0, 3, 1, 1, 0, 
"C - - 10 '\"Paschen_CF4.dat\"' fileName u 'TString' - 0 '\"\"' delim "
"u 'TString' - 0 '\"\"' opt", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DmtpcPaschen",1210,G__DmtpcPaschenCint_207_0_2, 105, G__get_linked_tagnum(&G__DmtpcPaschenCintLN_DmtpcPaschen), -1, 0, 1, 1, 1, 0, "u 'DmtpcPaschen' - 11 - other", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetName",673,G__DmtpcPaschenCint_207_0_3, 117, G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TString), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("getPaschenCurve",1543,G__DmtpcPaschenCint_207_0_4, 85, G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TGraph), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("readPaschenData",1496,G__DmtpcPaschenCint_207_0_5, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("makeFitFunction",1543,G__DmtpcPaschenCint_207_0_6, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getFitFunction",1449,G__DmtpcPaschenCint_207_0_7, 85, G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TF1), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("testme",658,G__DmtpcPaschenCint_207_0_8, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__DmtpcPaschenCint_207_0_9, 85, G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&DmtpcPaschen::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__DmtpcPaschenCint_207_0_10, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmtpcPaschen::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__DmtpcPaschenCint_207_0_11, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&DmtpcPaschen::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__DmtpcPaschenCint_207_0_12, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&DmtpcPaschen::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,G__DmtpcPaschenCint_207_0_13, 85, G__get_linked_tagnum(&G__DmtpcPaschenCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,G__DmtpcPaschenCint_207_0_14, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,G__DmtpcPaschenCint_207_0_15, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__DmtpcPaschenCint_207_0_16, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__DmtpcPaschenCint_207_0_17, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmtpcPaschen::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__DmtpcPaschenCint_207_0_18, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&DmtpcPaschen::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__DmtpcPaschenCint_207_0_19, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&DmtpcPaschen::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__DmtpcPaschenCint_207_0_20, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&DmtpcPaschen::DeclFileLine) ), 0);
   // automatic destructor
   G__memfunc_setup("~DmtpcPaschen", 1336, G__DmtpcPaschenCint_207_0_21, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__DmtpcPaschenCint_207_0_22, (int) ('u'), G__get_linked_tagnum(&G__DmtpcPaschenCintLN_DmtpcPaschen), -1, 1, 1, 1, 1, 0, "u 'DmtpcPaschen' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncDmtpcPaschenCint() {
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
extern "C" void G__cpp_setup_globalDmtpcPaschenCint() {
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

extern "C" void G__cpp_setup_funcDmtpcPaschenCint() {
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
G__linked_taginfo G__DmtpcPaschenCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_TString = { "TString" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_TF1 = { "TF1" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_TGraph = { "TGraph" , 99 , -1 };
G__linked_taginfo G__DmtpcPaschenCintLN_DmtpcPaschen = { "DmtpcPaschen" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableDmtpcPaschenCint() {
  G__DmtpcPaschenCintLN_TClass.tagnum = -1 ;
  G__DmtpcPaschenCintLN_TBuffer.tagnum = -1 ;
  G__DmtpcPaschenCintLN_TMemberInspector.tagnum = -1 ;
  G__DmtpcPaschenCintLN_TString.tagnum = -1 ;
  G__DmtpcPaschenCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcPaschenCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcPaschenCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__DmtpcPaschenCintLN_TF1.tagnum = -1 ;
  G__DmtpcPaschenCintLN_TGraph.tagnum = -1 ;
  G__DmtpcPaschenCintLN_DmtpcPaschen.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableDmtpcPaschenCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_TString);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_TF1);
   G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_TGraph);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__DmtpcPaschenCintLN_DmtpcPaschen),sizeof(DmtpcPaschen),-1,1792,(char*)NULL,G__setup_memvarDmtpcPaschen,G__setup_memfuncDmtpcPaschen);
}
extern "C" void G__cpp_setupDmtpcPaschenCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupDmtpcPaschenCint()");
  G__set_cpp_environmentDmtpcPaschenCint();
  G__cpp_setup_tagtableDmtpcPaschenCint();

  G__cpp_setup_inheritanceDmtpcPaschenCint();

  G__cpp_setup_typetableDmtpcPaschenCint();

  G__cpp_setup_memvarDmtpcPaschenCint();

  G__cpp_setup_memfuncDmtpcPaschenCint();
  G__cpp_setup_globalDmtpcPaschenCint();
  G__cpp_setup_funcDmtpcPaschenCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncDmtpcPaschenCint();
  return;
}
class G__cpp_setup_initDmtpcPaschenCint {
  public:
    G__cpp_setup_initDmtpcPaschenCint() { G__add_setup_func("DmtpcPaschenCint",(G__incsetup)(&G__cpp_setupDmtpcPaschenCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initDmtpcPaschenCint() { G__remove_setup_func("DmtpcPaschenCint"); }
};
G__cpp_setup_initDmtpcPaschenCint G__cpp_setup_initializerDmtpcPaschenCint;

