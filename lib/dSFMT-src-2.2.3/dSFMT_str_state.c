/** 
 * @file dSFMT_state.c 
 * @brief function to save and load the state of a dSFMT
 *
 * @author Andrea C G Mennucci (Scuola Normale Superiore)
 *
 * Copyright (C) 2010 Andrea C G Mennucci
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 *
 * Changes:
 *  James Spencer - define equivalent function to strdup for portability
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "dSFMT-params.h"
#include "dSFMT.h"

/**
 * This function returns a pointer to a new string which is a duplicate of string s.
 * See strdup(3) manpage.
 */
char *local_strdup(const char *s) {
    char *p = malloc(strlen(s) + 1);
    if (p) {
        strcpy(p, s);
    }
    return p;
}

/**
 * This function returns a string that represents the state of the dSFMT.
 * The string is allocated and should be free-d after use.
 *
 * @param dsfmt dsfmt state vector.
 * @param prefix a prefix to start all lines; if it is NULL, it is set to "dsfmt_"
 */
char *dsfmt_state_to_str(dsfmt_t *dsfmt, char *prefix)
{
  int i,l,p=0;
  char *str;
  if(prefix==NULL)
    prefix="dsfmt_";

  l=strlen(dsfmt_get_idstring())+(strlen(prefix)+110)*(DSFMT_N+3)+20;
  str=malloc(l);

  p=sprintf(str,"%sid=%s\n%sidx=%d\n",prefix,dsfmt_get_idstring(),prefix,dsfmt->idx);
  for(i=0;i<=DSFMT_N;i++) {
    w128_t *t=&(dsfmt->status[i]);
    //p+=sprintf(str+p,"%sstate[%03d]=%" PRIu64 " %" PRIu64 "\n", prefix, i, t->u[0], t->u[1]);
    p+=sprintf(str+p,"%sstate[%03d]=%la %la\n", prefix, i, t->d[0], t->d[1]);
    assert(p<l);
  }
  return str;
}


/**
 * This function reads a NULL terminated list of string that represents the state of the dSFMT, and fills the state.
 * It returns NULL if OK, or a string explaining the error, if any. (The error should be freed after use)
 *
 * @param dsfmt dsfmt state vector to be filled.
 * @param strlist the NULL terminated list of strings encoding the state
 * @param prefix the prefix that starts all lines; if it is NULL, it is set to "dsfmt_"
 */
char *dsfmt_strlist_to_state(dsfmt_t *dsfmt, char **strlist, char *prefix)
{
  int i,j,r,num_lines=0;
  char *str;
  if(prefix==NULL)
    prefix="dsfmt_";
  const int pl=strlen(prefix);
  while(strlist[num_lines]!=NULL && num_lines < (DSFMT_N+6))
    num_lines++;
  if(num_lines<=0)
    return local_strdup("dsfmt_strlist_to_state: the list of strings is empty\n");
  if(num_lines!= (DSFMT_N+3))
    return local_strdup("dsfmt_strlist_to_state: the list of strings has wrong length\n");

  str=strlist[0];
  if(strcmp(str+pl+3,dsfmt_get_idstring())) {
    char *s=malloc(strlen(str)+strlen(dsfmt_get_idstring())+110);
    sprintf(s,"dsfmt_strlist_to_state: the id differ, the string reports '%s' while the code requires '%s'\n",
	    str+pl+3,dsfmt_get_idstring());
    return s;
  }

  str=strlist[1];
  r=sscanf(str+pl,"idx=%d", &(dsfmt->idx));
  if(r!=1) {
    char *s=malloc(pl+strlen(str)+110);
    sprintf(s,"dsfmt_strlist_to_state: the second line is malformed, could not parse idx from '%s'\n", str);
    return s;
  }

  for(i=0;i<=DSFMT_N;i++) {
    w128_t *t=&(dsfmt->status[i]);
    str=strlist[i+2];
    //r=sscanf(str+pl,"state[%d]=%" PRIu64 "%" PRIu64, &j, &(t->u[0]), &(t->u[1]));
    r=sscanf(str+pl,"state[%d]=%la %la", &j, &(t->d[0]), &(t->d[1]));
    if(r!=3 || i!=j) {
      char *s=malloc(pl+strlen(str)+110);
      sprintf(s,"dsfmt_strlist_to_state: the line %d is malformed, could not parse from '%s'\n", 
	      i,str);
      return s;
    }
  }
  return NULL;
}


/**
 * This function reads a string that represents the state of the dSFMT, and fills the state.
 * It returns NULL if OK, or a string explaining the error, if any. (It should be freed after use)
 *
 * @param dsfmt dsfmt state vector.
 * @param origstr the string encoding the state
 * @param prefix the prefix that starts all lines; if it is NULL, it is set to "dsfmt_"
 */
char *dsfmt_str_to_state(dsfmt_t *dsfmt, char *origstr, char *prefix)
{
  int i;
  if(prefix==NULL)
    prefix="dsfmt_";
  const int pl=strlen(prefix), sl=strlen(origstr);
  if(sl<=0)
    return local_strdup("dsfmt_str_to_state: the string is empty\n");
  char *str=local_strdup(origstr);
  char **strlist=calloc(4+DSFMT_N,sizeof(char *));
  char *p=strchr(str,'\n');
  char *op=str;
  if(p==NULL) {
      free(str);
      free(strlist);
      return local_strdup("dsfmt_str_to_state: the string is garbage\n");
  }
  *p=0; p++;

  for(i=0;(i<=(2+DSFMT_N)) && (op<(str+sl)) && p; i++) {
    strlist[i]=op;
    if(strncmp(op,prefix,pl)) {
      char *s=malloc(pl+(p-op)+110);
      sprintf(s,"dsfmt_str_to_state: the prefix differs at line %d, the string is '%s' while the prefix should be '%s'\n",
	      i,op,prefix);
      free(str);free(strlist);
      return s;
    }
    op=p;
    if(p<str+sl) {
      p=strchr(p,'\n');
      if(p==NULL)  {
	free(str);
	free(strlist);
	return local_strdup("dsfmt_str_to_state: the string is truncated\n");
      }
      *p=0; p++;
    } else break;
  }
  if(p<str+sl) printf(" warning, garbage at the end of the string: %s\n",p);
  char *err=dsfmt_strlist_to_state(dsfmt, strlist, prefix);
  free(str);
  free(strlist);
  return err;
}

/**
 * This function reads a file descriptor that represents the state of the dSFMT, and fills the state.
 * It returns NULL if OK, or a string explaining the error, if any. (It should be freed after use)
 *
 * @param dsfmt dsfmt state vector.
 * @param origfile the file encoding the state
 * @param prefix the prefix that starts all lines; if it is NULL, it is set to "dsfmt_"
 */

static inline int string_is_comment(char *s)
{
  while(*s==' ' || *s == '\t') s++;
  return *s==0 || *s=='\n' || *s == '#';
}
static inline char *fgets_noncomment(char *s,int l,FILE *f)
{
  char *S=fgets(s,l,f);
  while(S && string_is_comment(s))
    S=fgets(s,l,f);
  return S;
}

char *dsfmt_file_to_state(dsfmt_t *dsfmt, FILE *origfile, char *prefix)
{
  int i;
  if(prefix==NULL)
    prefix="dsfmt_";
  const int pl=strlen(prefix);
  if(feof(origfile))
    return local_strdup("dsfmt_file_to_state: the file is empty\n");
  const int sl=pl+108;
  char **strlist=calloc(4+DSFMT_N,sizeof(char *));
  for(i=0;(i<=(2+DSFMT_N)); i++) { 
    strlist[i]=calloc(sl+4,1);
    char *S=fgets_noncomment(strlist[i],sl,origfile);
    if(!S) {
      for(;i>=0;i--) free(strlist[i]);
      free(strlist);
      return local_strdup("dsfmt_file_to_state: the file is too short\n");
    }
    if(strncmp(strlist[i],prefix,pl)) {
      char *s=malloc(pl+sl+110);
      sprintf(s,"dsfmt_str_to_state: the prefix differs at line %d, the string is '%s' while the prefix should be '%s'\n",
	      i,strlist[i],prefix);
      for(;i>=0;i--) free(strlist[i]);
      free(strlist);
      return s;
    }
    { char *p=strchr(strlist[i],'\n');
      if(p==NULL) {
      char *s=malloc(pl+sl+110);
      sprintf(s,"dsfmt_str_to_state: the line %d in the file is too long: '%s'\n", i,strlist[i]);
      for(;i>=0;i--) free(strlist[i]);
      free(strlist);
      return s;
      }
      *p=0; //delete newline
    }
  }
  char *err=dsfmt_strlist_to_state(dsfmt, strlist, prefix);
  for(i=0;(i<=(3+DSFMT_N)); i++) 
    free(strlist[i]);
  free(strlist);
  return err;
}

