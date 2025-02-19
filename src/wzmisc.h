#ifndef _WZMISC_H
#define _WZMISC_H

#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

/* utility function for freeing character array */
static inline void free_char_array(char **char_array, int n) {
  int i;
  for (i=0; i<n; ++i)
    free(char_array[i]);
  free(char_array);
}

/******************************************
 * same as strdup, but free/reuse original
 * Caution! src and dest must not overlap 
 * since realloc will release src and make 
 * its content undefined
 *******************************************/
static inline char *strcpy_realloc(char *dest, char *src) {
  dest = realloc(dest, strlen(src) + 1);
  strcpy(dest, src);
  return dest;
}

/* Exit with message */
/* static inline void wzfatal(const char *msg, ...) { */
/*   va_list args; */
/*   va_start(args, msg); */
/*   vfprintf(stderr, msg, args); */
/*   va_end (args); */
/*   fflush(stderr); */
/*   exit(EXIT_FAILURE); */
/* } */

static inline void wzfatal(const char *msg, ...) {
  REprintf("Error: %s\n", msg);
  error("An error occurred.");
}

static inline void wzfread(void *ptr, size_t size, size_t count, FILE *stream) {
  if (fread(ptr, size, count, stream) != count) {
    wzfatal("Reading error");
  }
}

/* convert string to uppercase */
static inline void wzstrupr(char *s) {
  for (;*s;++s) *s = toupper((unsigned char)*s);
}

static inline int is_number(char *s) {
  int i;
  for (i=0;s[i];++i) {
    if (!isdigit(s[i]) && s[i]!='.') {
      return 0;
    }
  }
  return 1;
}

static inline void ensure_number(char *s) {
  int i;
  for (i=0;s[i];++i) {
    if (!isdigit(s[i]) && s[i]!='.') {
      REprintf("[%s:%d] Trying to convert nondigit string to number: %s\n", __func__, __LINE__, s);
      error("Abort.");
    }
  }
}

/* allow negative */
static inline void ensure_number2(char *s) {
  int i;
  for (i=0;s[i];++i) {
    if (i==0 && s[0] == '-') continue;
    if (!isdigit(s[i]) && s[i]!='.') {
      REprintf("[%s:%d] Trying to convert nondigit string to number: %s\n", __func__, __LINE__, s);
      error("Abort.");
    }
  }
}

static inline int strcount_char(char *s, char c) {
  int i, n=0;
  for (i=0; s[i]; ++i)
    if (s[i] == c)
      ++n;
  return n;
}

static inline const char* get_basename(const char* filepath) {
  const char* last_slash = strrchr(filepath, '/'); // Find the last occurrence of '/'
  if (last_slash == NULL) {
    last_slash = strrchr(filepath, '\\'); // If on Windows, find the last occurrence of '\'
  }
  if (last_slash == NULL) {
    return filepath; // If no slash found, return the original filepath as it is the basename.
  }
  return last_slash + 1; // Increment the pointer to point to the character after the slash (basename).
}

/***************
 * min and max *
 ***************/
#define max(a,b)                \
  ({ __typeof__ (a) _a = (a);   \
    __typeof__ (b) _b = (b);    \
    _a > _b ? _a : _b; })

#define min(a,b)                \
  ({ __typeof__ (a) _a = (a);   \
    __typeof__ (b) _b = (b);    \
    _a > _b ? _b : _a; })


#endif /* _WZMISC_H */
