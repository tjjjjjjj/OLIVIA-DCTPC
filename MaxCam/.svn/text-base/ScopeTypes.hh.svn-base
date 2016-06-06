#ifndef SCOPE_TYPES_HH
#define SCOPE_TYPES_HH


//#ifdef __CINT__
//typedef int __signed__;
//#endif

enum ScopeClockSource   { SCOPE_CLOCK_INTERNAL, SCOPE_CLOCK_EXTERNAL };
enum ScopeClockEdge     { SCOPE_CLOCK_EDGE_RISING, SCOPE_CLOCK_EDGE_FALLING };
enum ScopeInputCoupling { SCOPE_COUPLING_AC=1, SCOPE_COUPLING_DC };
enum ScopeTriggerSlope  { SCOPE_TRIGGER_SLOPE_UNKNOWN,
			  SCOPE_TRIGGER_SLOPE_POSITIVE, 
			  SCOPE_TRIGGER_SLOPE_NEGATIVE };
enum ScopeTriggerSource { SCOPE_TRIGGER_SOURCE_DISABLE,
			  SCOPE_TRIGGER_SOURCE_EXTERNAL,
			  SCOPE_TRIGGER_SOURCE_CHAN_A, 
			  SCOPE_TRIGGER_SOURCE_CHAN_B};
enum ScopeTriggerEngine { SCOPE_TRIGGER_ENGINE_J,
			  SCOPE_TRIGGER_ENGINE_K};

enum ScopeTriggerEngOp { SCOPE_TRIGGER_ENG_OP_J,
			 SCOPE_TRIGGER_ENG_OP_K,
			 SCOPE_TRIGGER_ENG_OP_J_OR_K,
			 SCOPE_TRIGGER_ENG_OP_J_AND_K,
			 SCOPE_TRIGGER_ENG_OP_J_XOR_K,
			 SCOPE_TRIGGER_ENG_OP_J_AND_NOT_K,
			 SCOPE_TRIGGER_ENG_OP_NOT_J_AND_K};

// Alazar defines U32 and U8 types which are equivalent to
// unsigned int and unsigned char
typedef unsigned int AU32;
typedef unsigned char AU8;

//typedef signed __signed__;
//typedef __signed__ signed;
//#ifndef __signed__
//#define __signed__ = signed;
//#endif

#endif
