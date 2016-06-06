/* DM_UPS.c  
 *
 * Program to poll status of APC UPS.
 * 
 * compile with gcc DM_UPS.c -o DM_UPS
 *
 * Returns:
 *   0 if UPS is on line power,
 *   1 if UPS is on battery power,
 *   2 if UPS is discharged 
 *   3 if UPS is off  
 *   4 if UPS could not be reached
 *   5 if gethostbyname fails
 *   6 for malformed input
 *   7 if a socket could not be created or fcntl fails to set blocking mode 
 *   8 if authentication fails
 *   9 if timeout while connected
 *   10 if timeout but probably because someone else is telnetted in 
 */

//Define returns

#define RETURN_LINE 0 
#define RETURN_BATTERY 1 
#define RETURN_DISCHARGED 2 
#define RETURN_OFF 3
#define RETURN_NOREACH 4
#define RETURN_HOST 5
#define RETURN_BADINPUT 6
#define RETURN_SOCKET 7
#define RETURN_AUTH 8
#define RETURN_TIMEOUT 9
#define RETURN_BUSY 10

// toggle debug output
#define DEBUG 0

// Default port (telnet) 
#define DEFAULT_PORT 23

//secs before we give up trying to connect to UPS
#define CONNECT_TO 5 

//secs before we stop trying to read more from socket
#define TIMEOUT 1 

//usecs
#define SLEEP_AMT 50000 

// Default credentials
static char * DEFAULT_USER = "dmatter"; 
static char * DEFAULT_PASS = "seedark"; 

#include <stdio.h>
#include <netdb.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/fcntl.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>

/* The action to perform once connected */
typedef enum
{
  STATUS,
  START,
  STOP
} ups_action_t; 


/* Use ring buffer of size 3n to scan for text coming from a socket. Once buffer is longer than 2n, chop off the first n 
 * and move buffer over to left. returns 0 if found, 1 if timed out, -1 if bad string. 
 *
 * Will not consume any additional bytes. 
 */ 

int wait_for_string(int sock_fd, char * str) 
{
  uint32_t waited,pos; 
  int str_n = strlen(str); 
  char * buf = (char*) calloc(3*str_n,1); 

  if (str_n == 0) return -1; 

  waited = 0; 
  pos = 0; 
  while (waited < TIMEOUT*1e6)
  {
    if (DEBUG) printf("Sleeping...\n"); 
    waited+=SLEEP_AMT; 
    usleep(SLEEP_AMT); 
    if (pos >= 2*str_n)
    {
      pos -= str_n; 
      memcpy(buf+str_n,buf,pos); 
      memset(buf+2*str_n,0,str_n); 
    }

    while(read(sock_fd,buf+pos,1))  
    {
      if (DEBUG)printf("%c",buf[pos]); 
      pos++; 
      if (strstr(buf,str)!=NULL)
      {
        if (DEBUG) printf("\n"); 
        free(buf);
        return 0; 
      }

      if (pos == 2*str_n)
      {
        if (DEBUG) printf("\n\nTrimming buffer:\nBefore: \t\t %s\n",buf); 
        pos -= str_n; 
        memmove(buf,buf+str_n,pos); 
        memset(buf+2*str_n,0,str_n); 
        if (DEBUG) printf("\nAfter: \t\t %s\n",buf); 
      }
    }
  }

  if (DEBUG) printf("\n\n%s\ntiming out \n",buf); 
  free(buf); 
  return 1; 
}

void usage()
{
  printf("DM_UPS [-start|-stop|(-status)] host [port = %d] [username = %s] [pw = %s]\n",
          DEFAULT_PORT,DEFAULT_USER,DEFAULT_PASS); 
}

/* Print only first line of a buffer*/ 
void printline(FILE * f, char * str)
{
  int i,n; 
  n = strlen(str); 
  for (i = 0; i < n; i++)
  {
    if (str[i] == '\n' || str[i] == '\r')  
    {
      str[i] = '\0'; 
      break; 
    }
  }
  fprintf(f,"%s\n",str); 
}

/* Non blocking connect. Uses select to enforce timeout on nonblocking socket.
 * 
 * Will quit if there's an fcntl error, return negative if timeout or connection error 
 * and returns 0 on success. 
 *
 * timeout set via CONNECT_TO #defien
 */
int connect_nonblock(int sock_fd, struct sockaddr * addr, size_t size)
{
  long arg; 
  int res,valopt; 
  struct timeval t; 
  fd_set set; 
  socklen_t len; 

  //Change socket to nonblocking

  if (arg = fcntl(sock_fd, F_GETFL) <0)
  {
    fprintf(stderr,"fcntl error: %s\n",strerror(errno)); 
    exit(RETURN_SOCKET);
  }

  arg |= O_NONBLOCK; 

  if (fcntl(sock_fd, F_SETFL,arg) <0)
  {
    fprintf(stderr,"fcntl error: %s\n",strerror(errno)); 
    exit(RETURN_SOCKET); 
  }


  //Do connection 
  res = connect(sock_fd, addr, size); 

  //Apply timeout or fail
  if (res < 0)
  {
    if (errno == EINPROGRESS)
    {
      if (DEBUG) printf("einprogress\n"); 
      //Prepare select call
      t.tv_sec = CONNECT_TO; 
      t.tv_usec = 0; 

      FD_ZERO(&set); 
      FD_SET(sock_fd,&set); 

      res = select(sock_fd+1, NULL, &set, NULL, &t); 

      //Check select call result... if failed or timeout then fail 
      if (res <= 0)
      {
        return -1;  
      }

      len = sizeof(int);  
      //if getsockopt fails or valopt is 0 then we have a bad connection and should fail
      if (getsockopt(sock_fd,SOL_SOCKET,SO_ERROR, (void*)(&valopt), &len) || valopt) 
      {
        return -1; 
      }
    }
    //Any other errors we should fail 
    else return -1; 
  }

  // reblock

  if (arg =fcntl(sock_fd, F_GETFL) <0)
  {
    fprintf(stderr,"fcntl error: %s\n",strerror(errno)); 
    exit(RETURN_SOCKET);
  }

  arg &= ~O_NONBLOCK; 

  if (fcntl(sock_fd, F_SETFL,arg) <0)
  {
    fprintf(stderr,"fcntl error: %s\n",strerror(errno)); 
    exit(RETURN_SOCKET); 
  }

  return 0; 
}

int main (int nargs, char ** args)
{
  int sockfd, port, waited, rd,n;   
  char buf[64]; 
  char * hostn; 
  struct hostent * host; 
  struct sockaddr_in addr; 
  ups_action_t action = STATUS; 


  setbuf(stdout,0); 
  memset(buf,0,sizeof(buf));

  if (nargs < 2) 
  {
    usage(); 
    return RETURN_BADINPUT; 
  }

  n = 1;  // keeps track of the argument index

  if (args[n][0]=='-')
  {

    if (strcmp(args[n],"-start") == 0)  
    {
      action = START; 
    }
    else if (strcmp(args[n],"-stop") == 0)  
    {
      action = STOP; 
    }
    else if (strcmp(args[n],"-status") == 0)  
    {
      action = STATUS; 
    }
    else
    {
      fprintf(stderr,"Invalid argument: '%s'; assuming action is status\n",args[n]); 
    }

    n++; 
  }

  hostn = args[n++]; 
  host = (struct hostent *) gethostbyname(hostn); 
  if (host == NULL) 
  {
    fprintf(stderr,"Could not gethostbyname on %s\n",hostn); 
    return RETURN_HOST; 
  }


  port = DEFAULT_PORT; 
  if(nargs > n)
  { 
    port = atoi(args[n]); 
    if (port == 0)
    {
      fprintf(stderr,"Could not parse port string '%s' using default port %d\n",args[n], DEFAULT_PORT); 
      port = DEFAULT_PORT; 
    }
    n++; 
  } 
 
  // Open a socket
  sockfd = socket(AF_INET,SOCK_STREAM,0); 
  if (sockfd < 0)
  {
    fprintf(stderr,"Could not open socket\n"); 
    return RETURN_SOCKET; 
  }

  if (DEBUG) printf("Socket open\n");

  //Set up the addr struct 
  memset(&addr,0,sizeof(addr)); 
  addr.sin_family = AF_INET; 
  addr.sin_port = htons(port); //Use network byte order
  memcpy(&addr.sin_addr.s_addr, host->h_addr, host->h_length); 

  //Attempt to connect to the UPS
  if (connect_nonblock(sockfd,(struct sockaddr *) &addr,sizeof(addr)) < 0)
  {
    fprintf(stderr,"Could not connect to host %s\n",hostn); 
    return RETURN_NOREACH; 
  }

  char * user = DEFAULT_USER;
  char * pw = DEFAULT_PASS;
  if (nargs > n)
    user = args[n++]; 
  if (nargs > n)
    pw = args[n++]; 

  //Try to log in
  
  if (DEBUG) printf("Logging in\n");
  if (wait_for_string(sockfd,"User Name :"))
  {
    fprintf(stderr,"Timeout waiting for username: This probably means that someone else is logged in\n"); 
    return RETURN_BUSY; 
  }

  write(sockfd,user,strlen(user)); 
  write(sockfd,"\r\n",2); 

  if (DEBUG) printf("Sent username: %s\n",user);

  if (wait_for_string(sockfd,"Password  :"))
  {
    fprintf(stderr,"Timeout waiting for password\n"); 
    return RETURN_TIMEOUT; 
  }

  write(sockfd,pw,strlen(pw)); 
  write(sockfd,"\r\n",2); 
  if (DEBUG) printf("Sent password: %s\n",pw);

  if (wait_for_string(sockfd,"apc>"))
  {
    fprintf(stderr,"Auth failure or timeout waiting for prompt\n"); 
    return RETURN_AUTH; 
  }

  if (action == STOP || action == START)
  {
    if (action == STOP)
    {
      strcpy(buf,"ups -c Off");  
    }
    else if (action == START)
    {
      strcpy(buf,"ups -c On");  
    }

    write(sockfd,buf,strlen(buf)); 
    write(sockfd,"\r\n",2); 
  
    if (wait_for_string(sockfd,"E"))
    {
      fprintf(stderr,"Action '%s' timeout",buf); 
      return RETURN_TIMEOUT; 
    }
  }

  // Try to bring up the status
  strcpy(buf,"ups -st");  
  write(sockfd,buf,strlen(buf)); 
  write(sockfd,"\r\n",2); 

  if (wait_for_string(sockfd,"State: "))
  {
    fprintf(stderr,"Timeout or UPS not responding to ups -st\n"); 
    return RETURN_TIMEOUT; 
  }

  memset(buf,0,sizeof(buf)); 
  rd = 0; 
  waited = 0; 
  //Read the state
  while (waited < TIMEOUT * 1e6)
  {
    rd += read(sockfd,buf+rd,sizeof(buf)-rd); 

    if (strstr(buf,"On Line") != NULL)
    {
      printline(stdout,buf); 
      return RETURN_LINE;
    }
    else if (strstr(buf,"On Battery") != NULL)
    {
      printline(stdout,buf); 
      return RETURN_BATTERY; 
    }
    else if (strstr(buf,"Discharged") != NULL)
    {
      printline(stdout,buf); 
      return RETURN_DISCHARGED; 
    }
    else if (strstr(buf,"Off") != NULL)
    {
      printline(stdout,buf); 
      return RETURN_OFF; 
    }

    waited+= SLEEP_AMT;
    usleep(SLEEP_AMT); 
  }

  fprintf(stderr,"Timeout or Bad Status: %s\n",buf); 
  return RETURN_TIMEOUT; 
}

    

