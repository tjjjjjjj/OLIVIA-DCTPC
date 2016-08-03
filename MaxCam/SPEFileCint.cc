//
// File generated by rootcint at Tue Aug  2 13:55:55 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME SPEFileCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "SPEFileCint.h"

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
   void SPEFile_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void SPEFile_Dictionary();
   static void delete_SPEFile(void *p);
   static void deleteArray_SPEFile(void *p);
   static void destruct_SPEFile(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SPEFile*)
   {
      ::SPEFile *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SPEFile),0);
      static ::ROOT::TGenericClassInfo 
         instance("SPEFile", "./SPEFile.hh", 9,
                  typeid(::SPEFile), DefineBehavior(ptr, ptr),
                  0, &SPEFile_Dictionary, isa_proxy, 0,
                  sizeof(::SPEFile) );
      instance.SetDelete(&delete_SPEFile);
      instance.SetDeleteArray(&deleteArray_SPEFile);
      instance.SetDestructor(&destruct_SPEFile);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SPEFile*)
   {
      return GenerateInitInstanceLocal((::SPEFile*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::SPEFile*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void SPEFile_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::SPEFile*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SPEFile(void *p) {
      delete ((::SPEFile*)p);
   }
   static void deleteArray_SPEFile(void *p) {
      delete [] ((::SPEFile*)p);
   }
   static void destruct_SPEFile(void *p) {
      typedef ::SPEFile current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SPEFile

/********************************************************
* SPEFileCint.cc
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

extern "C" void G__cpp_reset_tagtableSPEFileCint();

extern "C" void G__set_cpp_environmentSPEFileCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("SPEFile.hh");
  G__cpp_reset_tagtableSPEFileCint();
}
#include <new>
extern "C" int G__cpp_dllrevSPEFileCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* SPEFile */
static int G__SPEFileCint_216_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   SPEFile* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new SPEFile((const char*) G__int(libp->para[0]));
   } else {
     p = new((void*) gvp) SPEFile((const char*) G__int(libp->para[0]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__SPEFileCintLN_SPEFile));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__SPEFileCint_216_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((SPEFile*) G__getstructoffset())->Getxpix());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__SPEFileCint_216_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((SPEFile*) G__getstructoffset())->Getypix());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__SPEFileCint_216_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((SPEFile*) G__getstructoffset())->Getdata_type());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__SPEFileCint_216_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 108, (long) ((SPEFile*) G__getstructoffset())->Getnimages());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__SPEFileCint_216_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((SPEFile*) G__getstructoffset())->GetImagetoROOT((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef SPEFile G__TSPEFile;
static int G__SPEFileCint_216_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (SPEFile*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((SPEFile*) (soff+(sizeof(SPEFile)*i)))->~G__TSPEFile();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (SPEFile*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((SPEFile*) (soff))->~G__TSPEFile();
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

/* SPEFile */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncSPEFileCint {
 public:
  G__Sizep2memfuncSPEFileCint(): p(&G__Sizep2memfuncSPEFileCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncSPEFileCint::*p)();
};

size_t G__get_sizep2memfuncSPEFileCint()
{
  G__Sizep2memfuncSPEFileCint a;
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
extern "C" void G__cpp_setup_inheritanceSPEFileCint() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableSPEFileCint() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__SPEFileCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__SPEFileCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__SPEFileCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__SPEFileCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__SPEFileCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__SPEFileCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__SPEFileCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__SPEFileCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__SPEFileCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__SPEFileCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Float_t>",117,G__get_linked_tagnum(&G__SPEFileCintLN_TVectorTlEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TVectorT<Double_t>",117,G__get_linked_tagnum(&G__SPEFileCintLN_TVectorTlEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTBase<Float_t>",117,G__get_linked_tagnum(&G__SPEFileCintLN_TMatrixTBaselEfloatgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("TMatrixTBase<Double_t>",117,G__get_linked_tagnum(&G__SPEFileCintLN_TMatrixTBaselEdoublegR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* SPEFile */
static void G__setup_memvarSPEFile(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__SPEFileCintLN_SPEFile));
   { SPEFile *p; p=(SPEFile*)0x1000; if (p) { }
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__SPEFileCintLN_SPEFilecLcLheader),-1,-1,4,"myheader=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"nbytesperim=",0,(char*)NULL);
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__SPEFileCintLN_basic_fstreamlEcharcOchar_traitslEchargRsPgR),G__defined_typename("fstream"),-1,4,"myfile=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarSPEFileCint() {
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
static void G__setup_memfuncSPEFile(void) {
   /* SPEFile */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__SPEFileCintLN_SPEFile));
   G__memfunc_setup("SPEFile",616,G__SPEFileCint_216_0_1, 105, G__get_linked_tagnum(&G__SPEFileCintLN_SPEFile), -1, 0, 1, 1, 1, 0, "C - - 10 - filename", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Getxpix",745,G__SPEFileCint_216_0_2, 105, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Getypix",746,G__SPEFileCint_216_0_3, 105, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Getdata_type",1243,G__SPEFileCint_216_0_4, 105, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Getnimages",1028,G__SPEFileCint_216_0_5, 108, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetImagetoROOT",1322,G__SPEFileCint_216_0_6, 85, G__get_linked_tagnum(&G__SPEFileCintLN_TH2F), -1, 0, 1, 1, 1, 0, "i - - 0 - a", (char*)NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~SPEFile", 742, G__SPEFileCint_216_0_7, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncSPEFileCint() {
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
extern "C" void G__cpp_setup_globalSPEFileCint() {
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

extern "C" void G__cpp_setup_funcSPEFileCint() {
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
G__linked_taginfo G__SPEFileCintLN_basic_fstreamlEcharcOchar_traitslEchargRsPgR = { "basic_fstream<char,char_traits<char> >" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_TVectorTlEfloatgR = { "TVectorT<float>" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_TVectorTlEdoublegR = { "TVectorT<double>" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_TMatrixTBaselEfloatgR = { "TMatrixTBase<float>" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_TMatrixTBaselEdoublegR = { "TMatrixTBase<double>" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_TH2F = { "TH2F" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_SPEFile = { "SPEFile" , 99 , -1 };
G__linked_taginfo G__SPEFileCintLN_SPEFilecLcLheader = { "SPEFile::header" , 115 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableSPEFileCint() {
  G__SPEFileCintLN_basic_fstreamlEcharcOchar_traitslEchargRsPgR.tagnum = -1 ;
  G__SPEFileCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__SPEFileCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__SPEFileCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__SPEFileCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__SPEFileCintLN_TVectorTlEfloatgR.tagnum = -1 ;
  G__SPEFileCintLN_TVectorTlEdoublegR.tagnum = -1 ;
  G__SPEFileCintLN_TMatrixTBaselEfloatgR.tagnum = -1 ;
  G__SPEFileCintLN_TMatrixTBaselEdoublegR.tagnum = -1 ;
  G__SPEFileCintLN_TH2F.tagnum = -1 ;
  G__SPEFileCintLN_SPEFile.tagnum = -1 ;
  G__SPEFileCintLN_SPEFilecLcLheader.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableSPEFileCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_basic_fstreamlEcharcOchar_traitslEchargRsPgR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_TVectorTlEfloatgR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_TVectorTlEdoublegR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_TMatrixTBaselEfloatgR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_TMatrixTBaselEdoublegR);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_TH2F);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__SPEFileCintLN_SPEFile),sizeof(SPEFile),-1,32768,(char*)NULL,G__setup_memvarSPEFile,G__setup_memfuncSPEFile);
   G__get_linked_tagnum_fwd(&G__SPEFileCintLN_SPEFilecLcLheader);
}
extern "C" void G__cpp_setupSPEFileCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupSPEFileCint()");
  G__set_cpp_environmentSPEFileCint();
  G__cpp_setup_tagtableSPEFileCint();

  G__cpp_setup_inheritanceSPEFileCint();

  G__cpp_setup_typetableSPEFileCint();

  G__cpp_setup_memvarSPEFileCint();

  G__cpp_setup_memfuncSPEFileCint();
  G__cpp_setup_globalSPEFileCint();
  G__cpp_setup_funcSPEFileCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncSPEFileCint();
  return;
}
class G__cpp_setup_initSPEFileCint {
  public:
    G__cpp_setup_initSPEFileCint() { G__add_setup_func("SPEFileCint",(G__incsetup)(&G__cpp_setupSPEFileCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initSPEFileCint() { G__remove_setup_func("SPEFileCint"); }
};
G__cpp_setup_initSPEFileCint G__cpp_setup_initializerSPEFileCint;

