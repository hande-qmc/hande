#include <stdlib.h>
#include "lua.h"

typedef struct
{
  int length;
  int space;
  char *container;
} charbuf;

// Writer to use during lua_dump.
static int buf_writer(lua_State *L, const void* p, size_t sz, void* ud)
{
  charbuf *dat;
  const char *buf;
  int i;

  dat = ud;
  buf = p;

  if ( sz + dat->length > dat->space ) {
    // Increase the size of the buffer, if needed.
    dat->space = ((dat->space*2) > (sz + dat->length))
               ? (dat->space*2) : (sz + dat->length);
    dat->container = realloc(dat->container, dat->space);
    if (!dat->container) return -10;
  }

  // Append the data to write into the buffer.
  for (i=0; i<sz; i++) {
    dat->container[dat->length + i] = buf[i];
  }
  dat->length = dat->length + sz;
  return 0;
}


// Wrapper around lua_dump to write into a memory buffer.
// Return Fortran friendly arguments.
const char* dump_lua_toBuf(lua_State *L, int *length, int *ierr)
{
  charbuf dat;
  char *buf;
  int i;
  int errcode;
  size_t sz;

  dat.length = 0;
  dat.space = 1024;
  dat.container = malloc(dat.space);

  errcode = lua_dump(L, buf_writer, &dat);

  (*ierr) = errcode;
  (*length) = dat.length;
  sz = dat.length;
  buf = malloc(dat.length);
  for (i=0; i<dat.length; i++) {
    buf[i] = dat.container[i];
  }
  free(dat.container);
  return buf;
}
