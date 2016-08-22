//
// File generated by rootcint at Tue Aug 16 18:40:19 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME DmtpcTheoryCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "DmtpcTheoryCint.h"

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

namespace DmtpcTheory {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static void DmtpcTheory_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("DmtpcTheory", 0 /*version*/, "./DmtpcTheory.hh", 10,
                     ::ROOT::DefineBehavior((void*)0,(void*)0),
                     &DmtpcTheory_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static void DmtpcTheory_Dictionary() {
         GenerateInitInstance()->GetClass();
      }

   }
}

/********************************************************
* DmtpcTheoryCint.cc
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

extern "C" void G__cpp_reset_tagtableDmtpcTheoryCint();

extern "C" void G__set_cpp_environmentDmtpcTheoryCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("DmtpcTheory.hh");
  G__cpp_reset_tagtableDmtpcTheoryCint();
}
#include <new>
extern "C" int G__cpp_dllrevDmtpcTheoryCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* DmtpcTheory */
static int G__DmtpcTheoryCint_209_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::ConvertXsectToRate((const Double_t) G__double(libp->para[0]), (const Double_t) G__double(libp->para[1])
, (const Double_t) G__double(libp->para[2]), (const Double_t) G__double(libp->para[3])
, (const Double_t) G__double(libp->para[4])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcTheoryCint_209_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::dRdERbetweenVearthandInfinity2DOverR0((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcTheoryCint_209_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::dRdERbetweenVearthandInfinity2D((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcTheoryCint_209_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::dRdERbetweenVearthandInfinity1DOverR0((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcTheoryCint_209_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::dRdERbetweenVearthandInfinity1D((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcTheoryCint_209_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::dRdERbetweenVearthandVesc1DOverR0((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcTheoryCint_209_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::dRdERbetweenVearthandVesc1D((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcTheoryCint_209_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::dRdThOverR0((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcTheoryCint_209_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) DmtpcTheory::dRdTh((Double_t*) G__int(libp->para[0]), (Double_t*) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* DmtpcTheory */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncDmtpcTheoryCint {
 public:
  G__Sizep2memfuncDmtpcTheoryCint(): p(&G__Sizep2memfuncDmtpcTheoryCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncDmtpcTheoryCint::*p)();
};

size_t G__get_sizep2memfuncDmtpcTheoryCint()
{
  G__Sizep2memfuncDmtpcTheoryCint a;
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
extern "C" void G__cpp_setup_inheritanceDmtpcTheoryCint() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableDmtpcTheoryCint() {

   /* Setting up typedef entry */
   G__search_typename2("Double_t",100,-1,0,-1);
   G__setnewtype(-1,"Double 8 bytes",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__DmtpcTheoryCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* DmtpcTheory */
static void G__setup_memvarDmtpcTheory(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__DmtpcTheoryCintLN_DmtpcTheory));
   {
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarDmtpcTheoryCint() {
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
static void G__setup_memfuncDmtpcTheory(void) {
   /* DmtpcTheory */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__DmtpcTheoryCintLN_DmtpcTheory));
   G__memfunc_setup("ConvertXsectToRate",1847,G__DmtpcTheoryCint_209_0_1, 100, -1, -1, 0, 5, 1, 1, 0, 
"d - 'Double_t' 10 - xsect d - 'Double_t' 10 - v_0 "
"d - 'Double_t' 10 - A d - 'Double_t' 10 - M_dark "
"d - 'Double_t' 10 - rho_D", (char*)NULL, (void*) G__func2void( (double (*)(const Double_t, const Double_t, const Double_t, const Double_t, const Double_t))(&DmtpcTheory::ConvertXsectToRate) ), 0);
   G__memfunc_setup("dRdERbetweenVearthandInfinity2DOverR0",3606,G__DmtpcTheoryCint_209_0_2, 100, -1, -1, 0, 2, 1, 1, 0, 
"D - 'Double_t' 0 - x D - 'Double_t' 0 - par", (char*)NULL, (void*) G__func2void( (double (*)(Double_t*, Double_t*))(&DmtpcTheory::dRdERbetweenVearthandInfinity2DOverR0) ), 0);
   G__memfunc_setup("dRdERbetweenVearthandInfinity2D",3064,G__DmtpcTheoryCint_209_0_3, 100, -1, -1, 0, 2, 1, 1, 0, 
"D - 'Double_t' 0 - x D - 'Double_t' 0 - par", (char*)NULL, (void*) G__func2void( (double (*)(Double_t*, Double_t*))(&DmtpcTheory::dRdERbetweenVearthandInfinity2D) ), 0);
   G__memfunc_setup("dRdERbetweenVearthandInfinity1DOverR0",3605,G__DmtpcTheoryCint_209_0_4, 100, -1, -1, 0, 2, 1, 1, 0, 
"D - 'Double_t' 0 - x D - 'Double_t' 0 - par", (char*)NULL, (void*) G__func2void( (double (*)(Double_t*, Double_t*))(&DmtpcTheory::dRdERbetweenVearthandInfinity1DOverR0) ), 0);
   G__memfunc_setup("dRdERbetweenVearthandInfinity1D",3063,G__DmtpcTheoryCint_209_0_5, 100, -1, -1, 0, 2, 1, 1, 0, 
"D - 'Double_t' 0 - x D - 'Double_t' 0 - par", (char*)NULL, (void*) G__func2void( (double (*)(Double_t*, Double_t*))(&DmtpcTheory::dRdERbetweenVearthandInfinity1D) ), 0);
   G__memfunc_setup("dRdERbetweenVearthandVesc1DOverR0",3164,G__DmtpcTheoryCint_209_0_6, 100, -1, -1, 0, 2, 1, 1, 0, 
"D - 'Double_t' 0 - x D - 'Double_t' 0 - par", (char*)NULL, (void*) G__func2void( (double (*)(Double_t*, Double_t*))(&DmtpcTheory::dRdERbetweenVearthandVesc1DOverR0) ), 0);
   G__memfunc_setup("dRdERbetweenVearthandVesc1D",2622,G__DmtpcTheoryCint_209_0_7, 100, -1, -1, 0, 2, 1, 1, 0, 
"D - 'Double_t' 0 - x D - 'Double_t' 0 - par", (char*)NULL, (void*) G__func2void( (double (*)(Double_t*, Double_t*))(&DmtpcTheory::dRdERbetweenVearthandVesc1D) ), 0);
   G__memfunc_setup("dRdThOverR0",1012,G__DmtpcTheoryCint_209_0_8, 100, -1, -1, 0, 2, 1, 1, 0, 
"D - 'Double_t' 0 - x D - 'Double_t' 0 - par", (char*)NULL, (void*) G__func2void( (double (*)(Double_t*, Double_t*))(&DmtpcTheory::dRdThOverR0) ), 0);
   G__memfunc_setup("dRdTh",470,G__DmtpcTheoryCint_209_0_9, 100, -1, -1, 0, 2, 1, 1, 0, 
"D - 'Double_t' 0 - x D - 'Double_t' 0 - par", (char*)NULL, (void*) G__func2void( (double (*)(Double_t*, Double_t*))(&DmtpcTheory::dRdTh) ), 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncDmtpcTheoryCint() {
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
extern "C" void G__cpp_setup_globalDmtpcTheoryCint() {
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

extern "C" void G__cpp_setup_funcDmtpcTheoryCint() {
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
G__linked_taginfo G__DmtpcTheoryCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcTheoryCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcTheoryCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__DmtpcTheoryCintLN_DmtpcTheory = { "DmtpcTheory" , 110 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableDmtpcTheoryCint() {
  G__DmtpcTheoryCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcTheoryCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcTheoryCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__DmtpcTheoryCintLN_DmtpcTheory.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableDmtpcTheoryCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__DmtpcTheoryCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcTheoryCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcTheoryCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcTheoryCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__DmtpcTheoryCintLN_DmtpcTheory),0,-1,0,(char*)NULL,G__setup_memvarDmtpcTheory,G__setup_memfuncDmtpcTheory);
}
extern "C" void G__cpp_setupDmtpcTheoryCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupDmtpcTheoryCint()");
  G__set_cpp_environmentDmtpcTheoryCint();
  G__cpp_setup_tagtableDmtpcTheoryCint();

  G__cpp_setup_inheritanceDmtpcTheoryCint();

  G__cpp_setup_typetableDmtpcTheoryCint();

  G__cpp_setup_memvarDmtpcTheoryCint();

  G__cpp_setup_memfuncDmtpcTheoryCint();
  G__cpp_setup_globalDmtpcTheoryCint();
  G__cpp_setup_funcDmtpcTheoryCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncDmtpcTheoryCint();
  return;
}
class G__cpp_setup_initDmtpcTheoryCint {
  public:
    G__cpp_setup_initDmtpcTheoryCint() { G__add_setup_func("DmtpcTheoryCint",(G__incsetup)(&G__cpp_setupDmtpcTheoryCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initDmtpcTheoryCint() { G__remove_setup_func("DmtpcTheoryCint"); }
};
G__cpp_setup_initDmtpcTheoryCint G__cpp_setup_initializerDmtpcTheoryCint;

