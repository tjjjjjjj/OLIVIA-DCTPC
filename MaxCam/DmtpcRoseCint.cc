//
// File generated by rootcint at Fri Jul 15 15:53:54 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME DmtpcRoseCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "DmtpcRoseCint.h"

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

namespace DmtpcRose {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static void DmtpcRose_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("DmtpcRose", 0 /*version*/, "./DmtpcRose.hh", 15,
                     ::ROOT::DefineBehavior((void*)0,(void*)0),
                     &DmtpcRose_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static void DmtpcRose_Dictionary() {
         GenerateInitInstance()->GetClass();
      }

   }
}

/********************************************************
* DmtpcRoseCint.cc
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

extern "C" void G__cpp_reset_tagtableDmtpcRoseCint();

extern "C" void G__set_cpp_environmentDmtpcRoseCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("DmtpcRose.hh");
  G__cpp_reset_tagtableDmtpcRoseCint();
}
#include <new>
extern "C" int G__cpp_dllrevDmtpcRoseCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* DmtpcRose */
static int G__DmtpcRoseCint_210_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 2:
      DmtpcRose::DrawRose((TObjArray*) G__int(libp->para[0]), *((TString*) G__int(libp->para[1])));
      G__setnull(result7);
      break;
   case 1:
      DmtpcRose::DrawRose((TObjArray*) G__int(libp->para[0]));
      G__setnull(result7);
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcRoseCint_210_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) DmtpcRose::makeRoseSlices((TH1*) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* DmtpcRose */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncDmtpcRoseCint {
 public:
  G__Sizep2memfuncDmtpcRoseCint(): p(&G__Sizep2memfuncDmtpcRoseCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncDmtpcRoseCint::*p)();
};

size_t G__get_sizep2memfuncDmtpcRoseCint()
{
  G__Sizep2memfuncDmtpcRoseCint a;
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
extern "C" void G__cpp_setup_inheritanceDmtpcRoseCint() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableDmtpcRoseCint() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcRoseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcRoseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcRoseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcRoseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__DmtpcRoseCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* DmtpcRose */
static void G__setup_memvarDmtpcRose(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__DmtpcRoseCintLN_DmtpcRose));
   {
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarDmtpcRoseCint() {
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
static void G__setup_memfuncDmtpcRose(void) {
   /* DmtpcRose */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__DmtpcRoseCintLN_DmtpcRose));
   G__memfunc_setup("DrawRose",807,G__DmtpcRoseCint_210_0_1, 121, -1, -1, 0, 2, 1, 1, 0, 
"U 'TObjArray' - 0 - slices u 'TString' - 0 '\"\"' opt", (char*)NULL, (void*) G__func2void( (void (*)(TObjArray*, TString))(&DmtpcRose::DrawRose) ), 0);
   G__memfunc_setup("makeRoseSlices",1434,G__DmtpcRoseCint_210_0_2, 85, G__get_linked_tagnum(&G__DmtpcRoseCintLN_TObjArray), -1, 0, 1, 1, 1, 0, "U 'TH1' - 0 - h", (char*)NULL, (void*) G__func2void( (TObjArray* (*)(TH1*))(&DmtpcRose::makeRoseSlices) ), 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncDmtpcRoseCint() {
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
extern "C" void G__cpp_setup_globalDmtpcRoseCint() {
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

extern "C" void G__cpp_setup_funcDmtpcRoseCint() {
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
G__linked_taginfo G__DmtpcRoseCintLN_TString = { "TString" , 99 , -1 };
G__linked_taginfo G__DmtpcRoseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__DmtpcRoseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcRoseCintLN_TObjArray = { "TObjArray" , 99 , -1 };
G__linked_taginfo G__DmtpcRoseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__DmtpcRoseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcRoseCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__DmtpcRoseCintLN_TH1 = { "TH1" , 99 , -1 };
G__linked_taginfo G__DmtpcRoseCintLN_DmtpcRose = { "DmtpcRose" , 110 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableDmtpcRoseCint() {
  G__DmtpcRoseCintLN_TString.tagnum = -1 ;
  G__DmtpcRoseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__DmtpcRoseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcRoseCintLN_TObjArray.tagnum = -1 ;
  G__DmtpcRoseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__DmtpcRoseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcRoseCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__DmtpcRoseCintLN_TH1.tagnum = -1 ;
  G__DmtpcRoseCintLN_DmtpcRose.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableDmtpcRoseCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_TString);
   G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_TObjArray);
   G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_TH1);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__DmtpcRoseCintLN_DmtpcRose),0,-1,0,(char*)NULL,G__setup_memvarDmtpcRose,G__setup_memfuncDmtpcRose);
}
extern "C" void G__cpp_setupDmtpcRoseCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupDmtpcRoseCint()");
  G__set_cpp_environmentDmtpcRoseCint();
  G__cpp_setup_tagtableDmtpcRoseCint();

  G__cpp_setup_inheritanceDmtpcRoseCint();

  G__cpp_setup_typetableDmtpcRoseCint();

  G__cpp_setup_memvarDmtpcRoseCint();

  G__cpp_setup_memfuncDmtpcRoseCint();
  G__cpp_setup_globalDmtpcRoseCint();
  G__cpp_setup_funcDmtpcRoseCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncDmtpcRoseCint();
  return;
}
class G__cpp_setup_initDmtpcRoseCint {
  public:
    G__cpp_setup_initDmtpcRoseCint() { G__add_setup_func("DmtpcRoseCint",(G__incsetup)(&G__cpp_setupDmtpcRoseCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initDmtpcRoseCint() { G__remove_setup_func("DmtpcRoseCint"); }
};
G__cpp_setup_initDmtpcRoseCint G__cpp_setup_initializerDmtpcRoseCint;

