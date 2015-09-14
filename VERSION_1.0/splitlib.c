/* splitlib.c */

#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

#define ERROR -1

#define DECK "#DECK"

#define RECLEN 80

extern char *optarg;
extern int optind, opterr;

char *directive, delims[]=" ,\n";

char card[RECLEN+1];
char poundcard[RECLEN+1];

FILE *lib, *outfile;

main(argc,argv)
int argc;
char *argv[];
{
     int fd;

     switch ( fd = parse(argc,argv) )  {
     case ERROR :
           show_usage("splitlib");
           exit(1);
           break;
     case 0 :
           lib = stdin;
           break;
     default :
           lib = fdopen(fd,"r");
           break;
     }

     /* main loop - read library file until EOF is reached */
     fgets(card, RECLEN, lib);
     directive = strtok(card,delims);
     while ( feof(lib) == 0 )
          wrtfile();
     fclose(lib);
}


int parse(argc,argv)
int argc;
char *argv[];
{
     int c, fd=ERROR;
     char libname[96], outdir[96];

     /* turn off error messages from getopt() */
     opterr = 0;

     /* set default directory equal to current working directory */
     strcpy(outdir,".");

     /* parse command line options */
     while ( (c=getopt(argc, argv, "hd:")) != EOF )  {
     switch (c)  {
     case 'h' :
         return(ERROR);
         break;
     case 'd' :
         if ( strlen(optarg) == 0 )  {
           fprintf(stderr,"\n No directory given. \n\n");
           return(ERROR);
         }
         strcpy(outdir,optarg);
         break;
     default:
         return(ERROR);
         break;
     }
     }

     /* see what's left on the command line, if no file name */
     /* is given use standard input */
     if ( optind >= argc ) 
       fd = 0;
     else {
       strcpy(libname,argv[optind]);
       /* check to see if library is accessible */
       if ( access(libname, 0) != 0 )  {
         fprintf(stderr,"\n %s - Bad file name or path.\n\n", libname);
         return(ERROR);
       }
  
       /* check to see if library can be read */
       if ( access(libname, 4) != 0 )  {
         fprintf(stderr,"\n %s - Permission denied.\n\n", libname);
         return(ERROR);
       }
  
       /* open library file for input */
       if ( (fd=open(libname, O_RDONLY)) == NULL ) {
         fprintf(stderr,"\n Unable to open %s.\n\n", libname);
         return(ERROR);
       }
    }

    /* If directory given change to it */
    if ( strcmp(outdir,".") != 0 )
      if ( chdir(outdir) == ERROR )  {
        fprintf(stderr,"\n %s - Permission denied.\n\n", outdir);
        return(ERROR);
      }

    return(fd);
}

show_usage(name)
char *name;
{
     printf("\tUSAGE:  %s [-h] [-d outdir] [file]\n\n",name);
     printf("\t        -h ; print this summary\n");
     printf("\t        -d outdir ; output directory, default=current directory\n");
     printf("\t        file ; if no file name is given standard input is assumed\n");
}

wrtfile()
{
     char name[95];

     /* build file name */
     strcpy(name,strtok(NULL,delims));

     /* open file for output */
     if ( (outfile=fopen(name,"w")) == NULL ) {
       fprintf(stderr,"\n Unable to open %s.\n\n", name);
       fclose(lib);
       exit(1);
     }

     /* loop until EOF is reached or another # card has been found */
     while ( fgets(card,RECLEN,lib) != NULL ) {
           /* check for a '#' in column 1 */
           if ( card[0] != '#' ) {
             fputs(card, outfile);
             continue;
           }
           strcpy(poundcard,card);
           directive = strtok(poundcard,delims);
           /* check for #DECK */
           if ( (strcmp(directive,DECK) == 0) )
             break;
           /* pass all other # directives through */
           fputs(card, outfile);
     }
     fclose(outfile);
     return;
}
