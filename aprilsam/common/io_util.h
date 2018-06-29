/* $COPYRIGHT_UM
$LICENSE_LGPL
*/

#ifndef _IOUTILS_H
#define _IOUTILS_H

#include <stdio.h>
#include <inttypes.h>

int ioutils_write_file(const char *path, const void *buf, long bufsz);

// read a file, allocating a buffer which user must free. Returns 0 on
// success. On failure, will always free any allocated
// memory. NULL-terminates because it's awesome. if maxalloc < 0, any file size is allowed.
//
// NB: ptr_to_buf should be of type char**, uint8_t** or similar. We
// will write the address of the allocated buffer to it.
int ioutils_read_file(const char *path, void *ptr_to_buf, long *_bufsz, long maxalloc);

/**
 * get_unique_filename:
 *
 * Returns an unused filename that can be used for writing without fear
 * of overwriting another file.  The returned string must be freed with
 * free() when no longer needed.  The new file is not opened.  The
 * filename is constructed as path/[date-]basename.XX.extension where
 * XX is an integer chosen to make it unique.  date- is only included
 * if @time_prefix is 1.
 *
 * @path: The path to the directory where the file should be stored.
 *    Pass NULL to use the current working directory or if the path
 *    is already included in @basename.
 * @basename: The desired name of the file without any trailing extension.
 *    This may include a leading path, but only if @time_prefix is zero.
 * @time_prefix: Set to 1 if the basename should be prefixed with the
 *    current date in the format YYYY-MM-DD-.
 * @extension: Filename extension to be used.  Use NULL if no extension
 *    is desired.
 */
char *
get_unique_filename (const char * path, const char * basename,
        uint8_t time_prefix, const char * extension);

int fdprintf(int fd, const char *fmt, ...);

// writes the number of bytes written (always len), unless an error
// occurs (in which a negative count is returned).
int write_fully(int fd, const void *b, int len);
int read_fully(int fd, void *b, int len);

/* returns -1 on error, 0 on timeout, nchar on success. */
int read_timeout(int fd, void *buf, int maxlen, int msTimeout);
int read_fully_timeout(int fd, void *bufin, int len, int msTimeout);
int read_line_timeout(int fd, void *buf, int maxlen, int msTimeout);
int read_line_timeout_ex(int fd, void *buf_in, int maxlen, int msTimeout, int *timed_out);

// returns 0 if a line ending in one of the newlinechars was found. *buflen will be set with the total number of characters in buf.
// returns -1 on timeout
// returns -2 on io error
// returns -3 on buffer too short.
// if maxlen is negative, any length string is allowed.
int read_line_timeout_flags(int fd, char **_buf, int *_buflen, int maxlen, int timeout_ms,
                            const char *newlinechars, const char *stripchars);

int read_available(int fd);
void read_flush(int fd);

int fwrite32(FILE *f, int32_t v);
int fwrite64(FILE *f, int64_t v64);
int fread32(FILE *f, int32_t *v32);
int fread64(FILE *f, int64_t *v64);

#endif
