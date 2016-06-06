private static uint CalcCRC(byte[] cmnd) {
  /* cmnd is a byte array containing the command ASCII string ë°≠ cmnd[]=-YÅ°Sinv2.000Å° */ 
  /* An unsigned 32 bit integer is return to the calling program */
  /* only the lower 16 bits contain the crc */
  int i,j;	/* interating indexes for the for loops */
  uint crc;	/* crc variable that will be returned */
  crc=0xffff;	/* initialize crc to hex value 0xffff */
  /* this for loop starts with ASCCII Ü°SÜ¢ and loops through to the last ASCII Ü°0Ü¢ */	
  for (i=0; i<cmnd.Length; i++) {
    /* the ASCII value is times by 0x0100 first then XORED to the current crc value */
    crc=crc^((uint)(cmnd[i]*0x0100));
    for(j=0; j<8; j++) {
      /* the crc is hashed 8 times with this for loop */
      /* if the 15th bit is set (tested by ANDING with hex 0x8000 and testing for 0x8000 result) then crc is shifted left one bit (same as times 2) XORED with hex 0x1021 and ANDED to hex 0xffff to limit the crc to lower 16 bits. If the 15th bit is not set then the crc is shifted left one bit and ANDED with hex 0xffff to limit the crc to lower 16 bits. */
      if((crc&0x8000)==0x8000)
	crc=((crc<<1)^0x1021)&0xffff;
      else
	crc=(crc<<1)&0xffff;
    }
  }/* There are some crc values that are not allowed, 0x00 and 0x0d */
  /* These are byte values so the high byte and the low byte of the crc must be checked and incremented if the bytes are either 0x00 0r 0x0d. */
  if((crc&0xff00)==0x0d00) crc +=0x0100;
  if((crc&0x00ff)==0x000d) crc +=0x0001;
  if((crc&0xff00)==0x0000) crc +=0x0100;
  if((crc&0x00ff)==0x0000) crc +=0x0001;
  return crc;		
}
