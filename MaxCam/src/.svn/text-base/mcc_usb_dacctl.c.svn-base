#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>
#include <sys/types.h>
#include <stdint.h>
#include <asm/types.h>
#include <libhid/pmd.h>
#include <libhid/usb-3100.h>

static uint16_t volts_to_digital(double v, bool bipolar)
{
   if (bipolar) 
    return  (uint16_t)  ((v/10. + 1.) * ((1 << 15)-1)); 

   return  (uint16_t)  ((v/10.) * ((1 << 16) -1)); 
}

int main(int nargs, char ** argv) 
{            
  HIDInterface *hid = 0; 
  hid_return ret; 
  int nInterfaces; 
  uint8_t chan; 
  double val; 
  uint16_t write; 


  if (nargs < 2)
  {
    printf("Usage: dac_control chan [vout = 0]\n"); 
    return 1; 
  }

  chan = atoi(argv[1]); 
  val = nargs < 3 ?  0 :  atof(argv[2]); 

  if (fabs(val) > 10)
  {
    printf("Max voltage: +/-10, aborting.");  
    return 2; 
  }

  write = volts_to_digital(val, val < 0);
  
  ret = hid_init(); 
  if (ret != HID_RET_SUCCESS) {
    fprintf(stderr, "hid_init failed with return code %d\n", ret);
    return -1;
  }

  if (PMD_Find_Interface(&hid, 0, USB3101_PID))
  {
    fprintf(stderr, "No interface found.\n");
    return -1;
  }

  usbAOutConfig_USB31XX(hid,chan, val < 0 ? 1 : 0); 
  usbAOut_USB31XX(hid, chan, write, 0); 

	hid_delete_HIDInterface(&hid);
	ret = hid_cleanup();

	if (ret != HID_RET_SUCCESS) {
	  fprintf(stderr, "hid_cleanup failed with return code %d\n", ret);
	  return 1;
	}
	return 0;


}
