//
// File generated by rootcint at Tue Aug 16 18:40:23 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME DmtpcStringToolsCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "DmtpcStringToolsCint.h"

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

namespace DmtpcStringTools {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static void DmtpcStringTools_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("DmtpcStringTools", 0 /*version*/, "./DmtpcStringTools.hh", 6,
                     ::ROOT::DefineBehavior((void*)0,(void*)0),
                     &DmtpcStringTools_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static void DmtpcStringTools_Dictionary() {
         GenerateInitInstance()->GetClass();
      }

   }
}

/********************************************************
* DmtpcStringToolsCint.cc
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

extern "C" void G__cpp_reset_tagtableDmtpcStringToolsCint();

extern "C" void G__set_cpp_environmentDmtpcStringToolsCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("DmtpcStringTools.hh");
  G__cpp_reset_tagtableDmtpcStringToolsCint();
}
#include <new>
extern "C" int G__cpp_dllrevDmtpcStringToolsCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* DmtpcStringTools */
static int G__DmtpcStringToolsCint_162_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      DmtpcStringTools::ltrim(*(string*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcStringToolsCint_162_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      DmtpcStringTools::rtrim(*(string*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcStringToolsCint_162_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      DmtpcStringTools::trim(*(string*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DmtpcStringToolsCint_162_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         string* pobj;
         string xobj = DmtpcStringTools::basename(*((string*) G__int(libp->para[0])));
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* DmtpcStringTools */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncDmtpcStringToolsCint {
 public:
  G__Sizep2memfuncDmtpcStringToolsCint(): p(&G__Sizep2memfuncDmtpcStringToolsCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncDmtpcStringToolsCint::*p)();
};

size_t G__get_sizep2memfuncDmtpcStringToolsCint()
{
  G__Sizep2memfuncDmtpcStringToolsCint a;
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
extern "C" void G__cpp_setup_inheritanceDmtpcStringToolsCint() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableDmtpcStringToolsCint() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* DmtpcStringTools */
static void G__setup_memvarDmtpcStringTools(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_DmtpcStringTools));
   {
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarDmtpcStringToolsCint() {
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
static void G__setup_memfuncDmtpcStringTools(void) {
   /* DmtpcStringTools */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_DmtpcStringTools));
   G__memfunc_setup("ltrim",552,G__DmtpcStringToolsCint_162_0_1, 121, -1, -1, 0, 1, 1, 1, 0, "u 'string' - 1 - str", (char*)NULL, (void*) G__func2void( (void (*)(string&))(&DmtpcStringTools::ltrim) ), 0);
   G__memfunc_setup("rtrim",558,G__DmtpcStringToolsCint_162_0_2, 121, -1, -1, 0, 1, 1, 1, 0, "u 'string' - 1 - str", (char*)NULL, (void*) G__func2void( (void (*)(string&))(&DmtpcStringTools::rtrim) ), 0);
   G__memfunc_setup("trim",444,G__DmtpcStringToolsCint_162_0_3, 121, -1, -1, 0, 1, 1, 1, 0, "u 'string' - 1 - str", (char*)NULL, (void*) G__func2void( (void (*)(string&))(&DmtpcStringTools::trim) ), 0);
   G__memfunc_setup("basename",828,G__DmtpcStringToolsCint_162_0_4, 117, G__get_linked_tagnum(&G__DmtpcStringToolsCintLN_string), -1, 0, 1, 1, 1, 0, "u 'string' - 0 - str", (char*)NULL, (void*) G__func2void( (string (*)(string))(&DmtpcStringTools::basename) ), 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncDmtpcStringToolsCint() {
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
extern "C" void G__cpp_setup_globalDmtpcStringToolsCint() {
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

extern "C" void G__cpp_setup_funcDmtpcStringToolsCint() {
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
G__linked_taginfo G__DmtpcStringToolsCintLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__DmtpcStringToolsCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcStringToolsCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DmtpcStringToolsCintLN_DmtpcStringTools = { "DmtpcStringTools" , 110 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableDmtpcStringToolsCint() {
  G__DmtpcStringToolsCintLN_string.tagnum = -1 ;
  G__DmtpcStringToolsCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcStringToolsCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DmtpcStringToolsCintLN_DmtpcStringTools.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableDmtpcStringToolsCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__DmtpcStringToolsCintLN_string);
   G__get_linked_tagnum_fwd(&G__DmtpcStringToolsCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DmtpcStringToolsCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__DmtpcStringToolsCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__DmtpcStringToolsCintLN_DmtpcStringTools),0,-1,0,(char*)NULL,G__setup_memvarDmtpcStringTools,G__setup_memfuncDmtpcStringTools);
}
extern "C" void G__cpp_setupDmtpcStringToolsCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupDmtpcStringToolsCint()");
  G__set_cpp_environmentDmtpcStringToolsCint();
  G__cpp_setup_tagtableDmtpcStringToolsCint();

  G__cpp_setup_inheritanceDmtpcStringToolsCint();

  G__cpp_setup_typetableDmtpcStringToolsCint();

  G__cpp_setup_memvarDmtpcStringToolsCint();

  G__cpp_setup_memfuncDmtpcStringToolsCint();
  G__cpp_setup_globalDmtpcStringToolsCint();
  G__cpp_setup_funcDmtpcStringToolsCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncDmtpcStringToolsCint();
  return;
}
class G__cpp_setup_initDmtpcStringToolsCint {
  public:
    G__cpp_setup_initDmtpcStringToolsCint() { G__add_setup_func("DmtpcStringToolsCint",(G__incsetup)(&G__cpp_setupDmtpcStringToolsCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initDmtpcStringToolsCint() { G__remove_setup_func("DmtpcStringToolsCint"); }
};
G__cpp_setup_initDmtpcStringToolsCint G__cpp_setup_initializerDmtpcStringToolsCint;

