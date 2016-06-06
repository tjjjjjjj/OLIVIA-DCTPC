//
// File generated by rootcint at Mon Jun  6 11:02:42 2016

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME MaxCamSerialCint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "MaxCamSerialCint.h"

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
   void MaxCamSerial_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_MaxCamSerial(void *p = 0);
   static void *newArray_MaxCamSerial(Long_t size, void *p);
   static void delete_MaxCamSerial(void *p);
   static void deleteArray_MaxCamSerial(void *p);
   static void destruct_MaxCamSerial(void *p);
   static void streamer_MaxCamSerial(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MaxCamSerial*)
   {
      ::MaxCamSerial *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MaxCamSerial >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MaxCamSerial", ::MaxCamSerial::Class_Version(), "./MaxCamSerial.hh", 26,
                  typeid(::MaxCamSerial), DefineBehavior(ptr, ptr),
                  &::MaxCamSerial::Dictionary, isa_proxy, 0,
                  sizeof(::MaxCamSerial) );
      instance.SetNew(&new_MaxCamSerial);
      instance.SetNewArray(&newArray_MaxCamSerial);
      instance.SetDelete(&delete_MaxCamSerial);
      instance.SetDeleteArray(&deleteArray_MaxCamSerial);
      instance.SetDestructor(&destruct_MaxCamSerial);
      instance.SetStreamerFunc(&streamer_MaxCamSerial);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MaxCamSerial*)
   {
      return GenerateInitInstanceLocal((::MaxCamSerial*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MaxCamSerial*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *MaxCamSerial::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *MaxCamSerial::Class_Name()
{
   return "MaxCamSerial";
}

//______________________________________________________________________________
const char *MaxCamSerial::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MaxCamSerial*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MaxCamSerial::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MaxCamSerial*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void MaxCamSerial::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MaxCamSerial*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *MaxCamSerial::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MaxCamSerial*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void MaxCamSerial::Streamer(TBuffer &R__b)
{
   // Stream an object of class MaxCamSerial.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, MaxCamSerial::IsA());
   } else {
      R__c = R__b.WriteVersion(MaxCamSerial::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void MaxCamSerial::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class MaxCamSerial.
      TClass *R__cl = ::MaxCamSerial::IsA();
      if (R__cl || R__insp.IsA()) { }
      TNamed::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MaxCamSerial(void *p) {
      return  p ? new(p) ::MaxCamSerial : new ::MaxCamSerial;
   }
   static void *newArray_MaxCamSerial(Long_t nElements, void *p) {
      return p ? new(p) ::MaxCamSerial[nElements] : new ::MaxCamSerial[nElements];
   }
   // Wrapper around operator delete
   static void delete_MaxCamSerial(void *p) {
      delete ((::MaxCamSerial*)p);
   }
   static void deleteArray_MaxCamSerial(void *p) {
      delete [] ((::MaxCamSerial*)p);
   }
   static void destruct_MaxCamSerial(void *p) {
      typedef ::MaxCamSerial current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_MaxCamSerial(TBuffer &buf, void *obj) {
      ((::MaxCamSerial*)obj)->::MaxCamSerial::Streamer(buf);
   }
} // end of namespace ROOT for class ::MaxCamSerial

/********************************************************
* MaxCamSerialCint.cc
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

extern "C" void G__cpp_reset_tagtableMaxCamSerialCint();

extern "C" void G__set_cpp_environmentMaxCamSerialCint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("MaxCamSerial.hh");
  G__cpp_reset_tagtableMaxCamSerialCint();
}
#include <new>
extern "C" int G__cpp_dllrevMaxCamSerialCint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* MaxCamSerial */
static int G__MaxCamSerialCint_178_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MaxCamSerial* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new MaxCamSerial[n];
     } else {
       p = new((void*) gvp) MaxCamSerial[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new MaxCamSerial;
     } else {
       p = new((void*) gvp) MaxCamSerial;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::Open(*((string*) G__int(libp->para[0]))));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 7:
      G__letint(result7, 105, (long) MaxCamSerial::Init(
(SerialHandle) G__int(libp->para[0]), (long) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (char) G__int(libp->para[3])
, (int) G__int(libp->para[4]), (int) G__int(libp->para[5])
, (int) G__int(libp->para[6])));
      break;
   case 6:
      G__letint(result7, 105, (long) MaxCamSerial::Init((SerialHandle) G__int(libp->para[0]), (long) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (char) G__int(libp->para[3])
, (int) G__int(libp->para[4]), (int) G__int(libp->para[5])));
      break;
   case 5:
      G__letint(result7, 105, (long) MaxCamSerial::Init((SerialHandle) G__int(libp->para[0]), (long) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (char) G__int(libp->para[3])
, (int) G__int(libp->para[4])));
      break;
   case 4:
      G__letint(result7, 105, (long) MaxCamSerial::Init((SerialHandle) G__int(libp->para[0]), (long) G__int(libp->para[1])
, (int) G__int(libp->para[2]), (char) G__int(libp->para[3])));
      break;
   case 3:
      G__letint(result7, 105, (long) MaxCamSerial::Init((SerialHandle) G__int(libp->para[0]), (long) G__int(libp->para[1])
, (int) G__int(libp->para[2])));
      break;
   case 2:
      G__letint(result7, 105, (long) MaxCamSerial::Init((SerialHandle) G__int(libp->para[0]), (long) G__int(libp->para[1])));
      break;
   case 1:
      G__letint(result7, 105, (long) MaxCamSerial::Init((SerialHandle) G__int(libp->para[0])));
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      MaxCamSerial::Abort((SerialHandle) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::Send((SerialHandle) G__int(libp->para[0]), (char*) G__int(libp->para[1])
, (float) G__double(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::SerialWrite((SerialHandle) G__int(libp->para[0]), (char*) G__int(libp->para[1])
, (int) G__int(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   switch (libp->paran) {
   case 3:
      {
         string* pobj;
         string xobj = MaxCamSerial::ReadPort((SerialHandle) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (bool) G__int(libp->para[2]));
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
      break;
   case 2:
      {
         string* pobj;
         string xobj = MaxCamSerial::ReadPort((SerialHandle) G__int(libp->para[0]), (int) G__int(libp->para[1]));
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
      break;
   }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamSerial::ReadPort((SerialHandle) G__int(libp->para[0]), (char*) G__int(libp->para[1])
, (int) G__int(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::SerialRead((SerialHandle) G__int(libp->para[0]), (char*) G__int(libp->para[1])
, (int) G__int(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::Close((SerialHandle) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::SetBaudRate((termios*) G__int(libp->para[0]), (long) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::SetDataBits((termios*) G__int(libp->para[0]), (int) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::SetParity((termios*) G__int(libp->para[0]), (char) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::SetStopBits((termios*) G__int(libp->para[0]), (int) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::SetFlowCtl((termios*) G__int(libp->para[0]), (int) G__int(libp->para[1])
, (int) G__int(libp->para[2])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      MaxCamSerial::Wait((float) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      MaxCamSerial::Flush((SerialHandle) G__int(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) MaxCamSerial::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamSerial::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) MaxCamSerial::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      MaxCamSerial::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MaxCamSerial*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_26(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamSerial::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_27(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_28(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MaxCamSerial::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MaxCamSerialCint_178_0_29(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MaxCamSerial::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__MaxCamSerialCint_178_0_30(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   MaxCamSerial* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new MaxCamSerial(*(MaxCamSerial*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef MaxCamSerial G__TMaxCamSerial;
static int G__MaxCamSerialCint_178_0_31(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
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
       delete[] (MaxCamSerial*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((MaxCamSerial*) (soff+(sizeof(MaxCamSerial)*i)))->~G__TMaxCamSerial();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (MaxCamSerial*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((MaxCamSerial*) (soff))->~G__TMaxCamSerial();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__MaxCamSerialCint_178_0_32(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MaxCamSerial* dest = (MaxCamSerial*) G__getstructoffset();
   *dest = *(MaxCamSerial*) libp->para[0].ref;
   const MaxCamSerial& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* MaxCamSerial */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncMaxCamSerialCint {
 public:
  G__Sizep2memfuncMaxCamSerialCint(): p(&G__Sizep2memfuncMaxCamSerialCint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncMaxCamSerialCint::*p)();
};

size_t G__get_sizep2memfuncMaxCamSerialCint()
{
  G__Sizep2memfuncMaxCamSerialCint a;
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
extern "C" void G__cpp_setup_inheritanceMaxCamSerialCint() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial))) {
     MaxCamSerial *G__Lderived;
     G__Lderived=(MaxCamSerial*)0x1000;
     {
       TNamed *G__Lpbase=(TNamed*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial),G__get_linked_tagnum(&G__MaxCamSerialCintLN_TNamed),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial),G__get_linked_tagnum(&G__MaxCamSerialCintLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableMaxCamSerialCint() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__MaxCamSerialCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__MaxCamSerialCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamSerialCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__MaxCamSerialCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamSerialCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__MaxCamSerialCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__MaxCamSerialCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamSerialCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__MaxCamSerialCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MaxCamSerialCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("SerialHandle",105,-1,0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* MaxCamSerial */
static void G__setup_memvarMaxCamSerial(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial));
   { MaxCamSerial *p; p=(MaxCamSerial*)0x1000; if (p) { }
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__MaxCamSerialCintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarMaxCamSerialCint() {
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
static void G__setup_memfuncMaxCamSerial(void) {
   /* MaxCamSerial */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial));
   G__memfunc_setup("MaxCamSerial",1175,G__MaxCamSerialCint_178_0_1, 105, G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Open",402,G__MaxCamSerialCint_178_0_2, 105, -1, G__defined_typename("SerialHandle"), 0, 1, 3, 1, 0, "u 'string' - 0 - device", (char*)NULL, (void*) G__func2void( (SerialHandle (*)(string))(&MaxCamSerial::Open) ), 0);
   G__memfunc_setup("Init",404,G__MaxCamSerialCint_178_0_3, 105, -1, -1, 0, 7, 3, 1, 0, 
"i - 'SerialHandle' 0 - rs232 l - - 0 '19200' baudRate "
"i - - 0 '8' dataBits c - - 0 ''n'' parity "
"i - - 0 '1' stopBits i - - 0 '0' xonxoff "
"i - - 0 '0' rtscts", (char*)NULL, (void*) G__func2void( (int (*)(SerialHandle, long, int, char, int, int, int))(&MaxCamSerial::Init) ), 0);
   G__memfunc_setup("Abort",504,G__MaxCamSerialCint_178_0_4, 121, -1, -1, 0, 1, 3, 1, 0, "i - 'SerialHandle' 0 - rs232", (char*)NULL, (void*) G__func2void( (void (*)(SerialHandle))(&MaxCamSerial::Abort) ), 0);
   G__memfunc_setup("Send",394,G__MaxCamSerialCint_178_0_5, 105, -1, -1, 0, 3, 3, 1, 0, 
"i - 'SerialHandle' 0 - rs232 C - - 0 - sendstring "
"f - - 0 - pause", (char*)NULL, (void*) G__func2void( (int (*)(SerialHandle, char*, float))(&MaxCamSerial::Send) ), 0);
   G__memfunc_setup("SerialWrite",1131,G__MaxCamSerialCint_178_0_6, 105, -1, -1, 0, 3, 3, 1, 0, 
"i - 'SerialHandle' 0 - rs232 C - - 0 - sendstring "
"i - - 0 - nbytes", (char*)NULL, (void*) G__func2void( (int (*)(SerialHandle, char*, int))(&MaxCamSerial::SerialWrite) ), 0);
   G__memfunc_setup("ReadPort",801,G__MaxCamSerialCint_178_0_7, 117, G__get_linked_tagnum(&G__MaxCamSerialCintLN_string), -1, 0, 3, 3, 1, 0, 
"i - 'SerialHandle' 0 - rs232 i - - 0 - readlength "
"g - - 0 'true' trim", (char*)NULL, (void*) G__func2void( (string (*)(SerialHandle, int, bool))(&MaxCamSerial::ReadPort) ), 0);
   G__memfunc_setup("ReadPort",801,G__MaxCamSerialCint_178_0_8, 67, -1, -1, 0, 3, 3, 1, 0, 
"i - 'SerialHandle' 0 - rs232 C - - 0 - readbuf "
"i - - 0 - readlength", (char*)NULL, (void*) G__func2void( (char* (*)(SerialHandle, char*, int))(&MaxCamSerial::ReadPort) ), 0);
   G__memfunc_setup("SerialRead",988,G__MaxCamSerialCint_178_0_9, 105, -1, -1, 0, 3, 3, 1, 0, 
"i - 'SerialHandle' 0 - rs232 C - - 0 - readbuf "
"i - - 0 - readlength", (char*)NULL, (void*) G__func2void( (int (*)(SerialHandle, char*, int))(&MaxCamSerial::SerialRead) ), 0);
   G__memfunc_setup("Close",502,G__MaxCamSerialCint_178_0_10, 105, -1, -1, 0, 1, 3, 1, 0, "i - 'SerialHandle' 0 - rs232", (char*)NULL, (void*) G__func2void( (int (*)(SerialHandle))(&MaxCamSerial::Close) ), 0);
   G__memfunc_setup("SetBaudRate",1076,G__MaxCamSerialCint_178_0_11, 105, -1, -1, 0, 2, 3, 1, 0, 
"U 'termios' - 0 - rs232_attr l - - 0 - baudrate", (char*)NULL, (void*) G__func2void( (int (*)(termios*, long))(&MaxCamSerial::SetBaudRate) ), 0);
   G__memfunc_setup("SetDataBits",1080,G__MaxCamSerialCint_178_0_12, 105, -1, -1, 0, 2, 3, 1, 0, 
"U 'termios' - 0 - rs232_attr i - - 0 - databits", (char*)NULL, (void*) G__func2void( (int (*)(termios*, int))(&MaxCamSerial::SetDataBits) ), 0);
   G__memfunc_setup("SetParity",933,G__MaxCamSerialCint_178_0_13, 105, -1, -1, 0, 2, 3, 1, 0, 
"U 'termios' - 0 - rs232_attr c - - 0 - parity", (char*)NULL, (void*) G__func2void( (int (*)(termios*, char))(&MaxCamSerial::SetParity) ), 0);
   G__memfunc_setup("SetStopBits",1124,G__MaxCamSerialCint_178_0_14, 105, -1, -1, 0, 2, 3, 1, 0, 
"U 'termios' - 0 - rs232_attr i - - 0 - stopbits", (char*)NULL, (void*) G__func2void( (int (*)(termios*, int))(&MaxCamSerial::SetStopBits) ), 0);
   G__memfunc_setup("SetFlowCtl",999,G__MaxCamSerialCint_178_0_15, 105, -1, -1, 0, 3, 3, 1, 0, 
"U 'termios' - 0 - rs232_attr i - - 0 - xonxoff "
"i - - 0 - rtscts", (char*)NULL, (void*) G__func2void( (int (*)(termios*, int, int))(&MaxCamSerial::SetFlowCtl) ), 0);
   G__memfunc_setup("Wait",405,G__MaxCamSerialCint_178_0_16, 121, -1, -1, 0, 1, 3, 1, 0, "f - - 0 - waittime", (char*)NULL, (void*) G__func2void( (void (*)(float))(&MaxCamSerial::Wait) ), 0);
   G__memfunc_setup("Flush",514,G__MaxCamSerialCint_178_0_17, 121, -1, -1, 0, 1, 3, 1, 0, "i - 'SerialHandle' 0 - rs232", (char*)NULL, (void*) G__func2void( (void (*)(SerialHandle))(&MaxCamSerial::Flush) ), 0);
   G__memfunc_setup("Class",502,G__MaxCamSerialCint_178_0_18, 85, G__get_linked_tagnum(&G__MaxCamSerialCintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&MaxCamSerial::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__MaxCamSerialCint_178_0_19, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamSerial::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__MaxCamSerialCint_178_0_20, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&MaxCamSerial::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__MaxCamSerialCint_178_0_21, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&MaxCamSerial::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__MaxCamSerialCintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - insp", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__MaxCamSerialCint_178_0_25, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__MaxCamSerialCint_178_0_26, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamSerial::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__MaxCamSerialCint_178_0_27, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MaxCamSerial::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__MaxCamSerialCint_178_0_28, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MaxCamSerial::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__MaxCamSerialCint_178_0_29, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MaxCamSerial::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("MaxCamSerial", 1175, G__MaxCamSerialCint_178_0_30, (int) ('i'), G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial), -1, 0, 1, 1, 1, 0, "u 'MaxCamSerial' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~MaxCamSerial", 1301, G__MaxCamSerialCint_178_0_31, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__MaxCamSerialCint_178_0_32, (int) ('u'), G__get_linked_tagnum(&G__MaxCamSerialCintLN_MaxCamSerial), -1, 1, 1, 1, 1, 0, "u 'MaxCamSerial' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncMaxCamSerialCint() {
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
extern "C" void G__cpp_setup_globalMaxCamSerialCint() {
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

extern "C" void G__cpp_setup_funcMaxCamSerialCint() {
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
G__linked_taginfo G__MaxCamSerialCintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_TNamed = { "TNamed" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_termios = { "termios" , 115 , -1 };
G__linked_taginfo G__MaxCamSerialCintLN_MaxCamSerial = { "MaxCamSerial" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableMaxCamSerialCint() {
  G__MaxCamSerialCintLN_TClass.tagnum = -1 ;
  G__MaxCamSerialCintLN_TBuffer.tagnum = -1 ;
  G__MaxCamSerialCintLN_TMemberInspector.tagnum = -1 ;
  G__MaxCamSerialCintLN_TObject.tagnum = -1 ;
  G__MaxCamSerialCintLN_TNamed.tagnum = -1 ;
  G__MaxCamSerialCintLN_string.tagnum = -1 ;
  G__MaxCamSerialCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__MaxCamSerialCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__MaxCamSerialCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__MaxCamSerialCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__MaxCamSerialCintLN_termios.tagnum = -1 ;
  G__MaxCamSerialCintLN_MaxCamSerial.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableMaxCamSerialCint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_TClass);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_TObject);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_TNamed);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_string);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_termios);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__MaxCamSerialCintLN_MaxCamSerial),sizeof(MaxCamSerial),-1,61696,(char*)NULL,G__setup_memvarMaxCamSerial,G__setup_memfuncMaxCamSerial);
}
extern "C" void G__cpp_setupMaxCamSerialCint(void) {
  G__check_setup_version(30051515,"G__cpp_setupMaxCamSerialCint()");
  G__set_cpp_environmentMaxCamSerialCint();
  G__cpp_setup_tagtableMaxCamSerialCint();

  G__cpp_setup_inheritanceMaxCamSerialCint();

  G__cpp_setup_typetableMaxCamSerialCint();

  G__cpp_setup_memvarMaxCamSerialCint();

  G__cpp_setup_memfuncMaxCamSerialCint();
  G__cpp_setup_globalMaxCamSerialCint();
  G__cpp_setup_funcMaxCamSerialCint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncMaxCamSerialCint();
  return;
}
class G__cpp_setup_initMaxCamSerialCint {
  public:
    G__cpp_setup_initMaxCamSerialCint() { G__add_setup_func("MaxCamSerialCint",(G__incsetup)(&G__cpp_setupMaxCamSerialCint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initMaxCamSerialCint() { G__remove_setup_func("MaxCamSerialCint"); }
};
G__cpp_setup_initMaxCamSerialCint G__cpp_setup_initializerMaxCamSerialCint;

