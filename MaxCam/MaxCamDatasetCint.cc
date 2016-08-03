//
// File generated by rootcint at Tue Aug  2 13:54:16 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME MaxCamDatasetCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "MaxCamDatasetCint.h"

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
   void MaxCamDataset_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_MaxCamDataset(void *p = 0);
   static void *newArray_MaxCamDataset(Long_t size, void *p);
   static void delete_MaxCamDataset(void *p);
   static void deleteArray_MaxCamDataset(void *p);
   static void destruct_MaxCamDataset(void *p);
   static void streamer_MaxCamDataset(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MaxCamDataset*)
   {
      ::MaxCamDataset *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MaxCamDataset >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MaxCamDataset", ::MaxCamDataset::Class_Version(), "./MaxCamDataset.hh", 16,
                  typeid(::MaxCamDataset), DefineBehavior(ptr, ptr),
                  &::MaxCamDataset::Dictionary, isa_proxy, 0,
                  sizeof(::MaxCamDataset) );
      instance.SetNew(&new_MaxCamDataset);
      instance.SetNewArray(&newArray_MaxCamDataset);
      instance.SetDelete(&delete_MaxCamDataset);
      instance.SetDeleteArray(&deleteArray_MaxCamDataset);
      instance.SetDestructor(&destruct_MaxCamDataset);
      instance.SetStreamerFunc(&streamer_MaxCamDataset);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MaxCamDataset*)
   {
      return GenerateInitInstanceLocal((::MaxCamDataset*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MaxCamDataset*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *MaxCamDataset::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *MaxCamDataset::Class_Name()
{
   return "MaxCamDataset";
}

//______________________________________________________________________________
const char *MaxCamDataset::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MaxCamDataset*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MaxCamDataset::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MaxCamDataset*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void MaxCamDataset::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MaxCamDataset*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *MaxCamDataset::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MaxCamDataset*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void MaxCamDataset::Streamer(TBuffer &R__b)
{
   // Stream an object of class MaxCamDataset.

   ::Error("MaxCamDataset::Streamer", "version id <=0 in ClassDef, dummy Streamer() called"); if (R__b.IsReading()) { }
}

//______________________________________________________________________________
void MaxCamDataset::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class MaxCamDataset.
      TClass *R__cl = ::MaxCamDataset::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_imageTree", &_imageTree);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_file", &_file);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_img_histo", &_img_histo);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_meshHV", &_meshHV);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_wireHV", &_wireHV);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_pressure", &_pressure);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_config", &_config);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "*_dtime", &_dtime);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MaxCamDataset(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::MaxCamDataset : new ::MaxCamDataset;
   }
   static void *newArray_MaxCamDataset(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::MaxCamDataset[nElements] : new ::MaxCamDataset[nElements];
   }
   // Wrapper around operator delete
   static void delete_MaxCamDataset(void *p) {
      delete ((::MaxCamDataset*)p);
   }
   static void deleteArray_MaxCamDataset(void *p) {
      delete [] ((::MaxCamDataset*)p);
   }
   static void destruct_MaxCamDataset(void *p) {
      typedef ::MaxCamDataset current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_MaxCamDataset(TBuffer &buf, void *obj) {
      ((::MaxCamDataset*)obj)->::MaxCamDataset::Streamer(buf);
   }
} // end of namespace ROOT for class ::MaxCamDataset

/********************************************************
* MaxCamDatasetCint.cc
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

extern "C" void G__cpp_reset_tagtableMaxCamDatasetCint();

extern "C" void G__set_cpp_environmentMaxCamDatasetCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("MaxCamDataset.hh");
  G__cpp_reset_tagtableMaxCamDatasetCint();
}
#include <new>
extern "C" int G__cpp_dllrevMaxCamDatasetCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* MaxCamDataset */
static int G__MaxCamDatasetCint_240_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MaxCamDataset* p = NULL;
   char* gvp = (char*) G__getgvp();
   switch (libp->paran) {
   case 2:
     //m: 2
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new MaxCamDataset((const char*) G__int(libp->para[0]), *((TString*) G__int(libp->para[1])));
     } else {
       p = new((void*) gvp) MaxCamDataset((const char*) G__int(libp->para[0]), *((TString*) G__int(libp->para[1])));
     }
     break;
   case 1:
     //m: 1
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new MaxCamDataset((const char*) G__int(libp->para[0]));
     } else {
       p = new((void*) gvp) MaxCamDataset((const char*) G__int(libp->para[0]));
     }
     break;
   case 0:
     int n = G__getaryconstruct();
     if (n) {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new MaxCamDataset[n];
       } else {
         p = new((void*) gvp) MaxCamDataset[n];
       }
     } else {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new MaxCamDataset;
       } else {
         p = new((void*) gvp) MaxCamDataset;
       }
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamDataset));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MaxCamDataset* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new MaxCamDataset(*(MaxCamDataset*) libp->para[0].ref);
   } else {
     p = new((void*) gvp) MaxCamDataset(*(MaxCamDataset*) libp->para[0].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamDataset));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         TString* pobj;
         TString xobj = ((MaxCamDataset*) G__getstructoffset())->GetName();
         pobj = new TString(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamDataset*) G__getstructoffset())->Write();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((MaxCamDataset*) G__getstructoffset())->wire());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((MaxCamDataset*) G__getstructoffset())->mesh());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((MaxCamDataset*) G__getstructoffset())->pressure());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((MaxCamDataset*) G__getstructoffset())->tree());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((MaxCamDataset*) G__getstructoffset())->chain());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((MaxCamDataset*) G__getstructoffset())->ccdConfig());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((MaxCamDataset*) G__getstructoffset())->timeStamp());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((MaxCamDataset*) G__getstructoffset())->ccdImage());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamDataset*) G__getstructoffset())->fill();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) MaxCamDataset::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamDataset::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) MaxCamDataset::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      MaxCamDataset::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((const MaxCamDataset*) G__getstructoffset())->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamDataset*) G__getstructoffset())->ShowMembers(*(TMemberInspector*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamDataset*) G__getstructoffset())->Streamer(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamDataset*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamDataset::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamDataset::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamDataset::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamDatasetCint_240_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamDataset::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef MaxCamDataset G__TMaxCamDataset;
static int G__MaxCamDatasetCint_240_0_26(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (MaxCamDataset*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((MaxCamDataset*) (soff+(sizeof(MaxCamDataset)*i)))->~G__TMaxCamDataset();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (MaxCamDataset*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((MaxCamDataset*) (soff))->~G__TMaxCamDataset();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__MaxCamDatasetCint_240_0_27(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MaxCamDataset* dest = (MaxCamDataset*) G__getstructoffset();
   *dest = *(MaxCamDataset*) libp->para[0].ref;
   const MaxCamDataset& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* MaxCamDataset */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncMaxCamDatasetCint {
 public:
  G__Sizep2memfuncMaxCamDatasetCint(): p(&G__Sizep2memfuncMaxCamDatasetCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncMaxCamDatasetCint::*p)();
};

size_t G__get_sizep2memfuncMaxCamDatasetCint()
{
  G__Sizep2memfuncMaxCamDatasetCint a;
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
extern "C" void G__cpp_setup_inheritanceMaxCamDatasetCint() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableMaxCamDatasetCint() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<std::bidirectional_iterator_tag,TObject*,std::ptrdiff_t,const TObject**,const TObject*&>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("iterator<bidirectional_iterator_tag,TObject*,long,const TObject**>",117,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* MaxCamDataset */
static void G__setup_memvarMaxCamDataset(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamDataset));
   { MaxCamDataset *p; p=(MaxCamDataset*)0x1000; if (p) { }
   G__memvar_setup((void*)0,108,0,0,-1,-1,-1,4,"G__virtualinfo=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->_imageTree)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TTree),-1,-1,1,"_imageTree=",0,"tree with events");
   G__memvar_setup((void*)((long)(&p->_file)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TFile),-1,-1,1,"_file=",0,"event file");
   G__memvar_setup((void*)((long)(&p->_img_histo)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TH2F),-1,-1,1,"_img_histo=",0,"image histogram");
   G__memvar_setup((void*)((long)(&p->_meshHV)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamChannel),-1,-1,1,"_meshHV=",0,"dift voltage");
   G__memvar_setup((void*)((long)(&p->_wireHV)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamChannel),-1,-1,1,"_wireHV=",0,"anode voltage");
   G__memvar_setup((void*)((long)(&p->_pressure)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamChannel),-1,-1,1,"_pressure=",0,"pressure");
   G__memvar_setup((void*)((long)(&p->_config)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamConfig),-1,-1,1,"_config=",0,"ccd configuration");
   G__memvar_setup((void*)((long)(&p->_dtime)-(long)(p)),85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TDatime),-1,-1,1,"_dtime=",0,"time-stamp");
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarMaxCamDatasetCint() {
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
static void G__setup_memfuncMaxCamDataset(void) {
   /* MaxCamDataset */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamDataset));
   G__memfunc_setup("MaxCamDataset",1277,G__MaxCamDatasetCint_240_0_1, 105, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamDataset), -1, 0, 2, 1, 1, 0, 
"C - - 10 '\"dataset.root\"' fileName u 'TString' - 0 '\"\"' foption", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("MaxCamDataset",1277,G__MaxCamDatasetCint_240_0_2, 105, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamDataset), -1, 0, 1, 1, 1, 0, "u 'MaxCamDataset' - 11 - other", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetName",673,G__MaxCamDatasetCint_240_0_3, 117, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TString), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Write",523,G__MaxCamDatasetCint_240_0_4, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("wire",439,G__MaxCamDatasetCint_240_0_5, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamChannel), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("mesh",429,G__MaxCamDatasetCint_240_0_6, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamChannel), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("pressure",889,G__MaxCamDatasetCint_240_0_7, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamChannel), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("tree",432,G__MaxCamDatasetCint_240_0_8, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TTree), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("chain",515,G__MaxCamDatasetCint_240_0_9, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TChain), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("ccdConfig",896,G__MaxCamDatasetCint_240_0_10, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamConfig), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("timeStamp",948,G__MaxCamDatasetCint_240_0_11, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TDatime), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("ccdImage",781,G__MaxCamDatasetCint_240_0_12, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TH2F), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("fill",423,G__MaxCamDatasetCint_240_0_13, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__MaxCamDatasetCint_240_0_14, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&MaxCamDataset::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__MaxCamDatasetCint_240_0_15, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamDataset::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__MaxCamDatasetCint_240_0_16, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&MaxCamDataset::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__MaxCamDatasetCint_240_0_17, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&MaxCamDataset::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,G__MaxCamDatasetCint_240_0_18, 85, G__get_linked_tagnum(&G__MaxCamDatasetCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,G__MaxCamDatasetCint_240_0_19, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,G__MaxCamDatasetCint_240_0_20, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__MaxCamDatasetCint_240_0_21, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__MaxCamDatasetCint_240_0_22, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamDataset::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__MaxCamDatasetCint_240_0_23, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MaxCamDataset::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__MaxCamDatasetCint_240_0_24, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamDataset::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__MaxCamDatasetCint_240_0_25, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MaxCamDataset::DeclFileLine) ), 0);
   // automatic destructor
   G__memfunc_setup("~MaxCamDataset", 1403, G__MaxCamDatasetCint_240_0_26, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__MaxCamDatasetCint_240_0_27, (int) ('u'), G__get_linked_tagnum(&G__MaxCamDatasetCintLN_MaxCamDataset), -1, 1, 1, 1, 1, 0, "u 'MaxCamDataset' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncMaxCamDatasetCint() {
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
extern "C" void G__cpp_setup_globalMaxCamDatasetCint() {
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

extern "C" void G__cpp_setup_funcMaxCamDatasetCint() {
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
G__linked_taginfo G__MaxCamDatasetCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_TString = { "TString" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR = { "iterator<bidirectional_iterator_tag,TObject*,long,const TObject**,const TObject*&>" , 115 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_TDatime = { "TDatime" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_TFile = { "TFile" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_MaxCamConfig = { "MaxCamConfig" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_MaxCamChannel = { "MaxCamChannel" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_TChain = { "TChain" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_TTree = { "TTree" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_TH2F = { "TH2F" , 99 , -1 };
G__linked_taginfo G__MaxCamDatasetCintLN_MaxCamDataset = { "MaxCamDataset" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableMaxCamDatasetCint() {
  G__MaxCamDatasetCintLN_TClass.tagnum = -1 ;
  G__MaxCamDatasetCintLN_TBuffer.tagnum = -1 ;
  G__MaxCamDatasetCintLN_TMemberInspector.tagnum = -1 ;
  G__MaxCamDatasetCintLN_TString.tagnum = -1 ;
  G__MaxCamDatasetCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__MaxCamDatasetCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__MaxCamDatasetCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR.tagnum = -1 ;
  G__MaxCamDatasetCintLN_TDatime.tagnum = -1 ;
  G__MaxCamDatasetCintLN_TFile.tagnum = -1 ;
  G__MaxCamDatasetCintLN_MaxCamConfig.tagnum = -1 ;
  G__MaxCamDatasetCintLN_MaxCamChannel.tagnum = -1 ;
  G__MaxCamDatasetCintLN_TChain.tagnum = -1 ;
  G__MaxCamDatasetCintLN_TTree.tagnum = -1 ;
  G__MaxCamDatasetCintLN_TH2F.tagnum = -1 ;
  G__MaxCamDatasetCintLN_MaxCamDataset.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableMaxCamDatasetCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TString);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TDatime);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TFile);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_MaxCamConfig);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_MaxCamChannel);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TChain);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TTree);
   G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_TH2F);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__MaxCamDatasetCintLN_MaxCamDataset),sizeof(MaxCamDataset),-1,1792,(char*)NULL,G__setup_memvarMaxCamDataset,G__setup_memfuncMaxCamDataset);
}
extern "C" void G__cpp_setupMaxCamDatasetCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupMaxCamDatasetCint()");
  G__set_cpp_environmentMaxCamDatasetCint();
  G__cpp_setup_tagtableMaxCamDatasetCint();

  G__cpp_setup_inheritanceMaxCamDatasetCint();

  G__cpp_setup_typetableMaxCamDatasetCint();

  G__cpp_setup_memvarMaxCamDatasetCint();

  G__cpp_setup_memfuncMaxCamDatasetCint();
  G__cpp_setup_globalMaxCamDatasetCint();
  G__cpp_setup_funcMaxCamDatasetCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncMaxCamDatasetCint();
  return;
}
class G__cpp_setup_initMaxCamDatasetCint {
  public:
    G__cpp_setup_initMaxCamDatasetCint() { G__add_setup_func("MaxCamDatasetCint",(G__incsetup)(&G__cpp_setupMaxCamDatasetCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initMaxCamDatasetCint() { G__remove_setup_func("MaxCamDatasetCint"); }
};
G__cpp_setup_initMaxCamDatasetCint G__cpp_setup_initializerMaxCamDatasetCint;

